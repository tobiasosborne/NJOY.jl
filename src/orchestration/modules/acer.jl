# acer module runner — ACE file generation from PENDF.
#
# Matches Fortran acer.f90 pipeline entry point: reads a PENDF tape (usually
# produced by reconr→broadr→heatr→thermr), builds an ACENeutronTable via
# the existing src/formats/ace_builder.jl infrastructure, writes a Type-1
# (ASCII) ACE file + MCNP xsdir directory line.
#
# Current coverage: iopt=1 (fast continuous-energy neutron ACE). Other iopt
# modes (iopt=2 thermal, iopt=7 photoatomic, iopt=8 dosimetry, etc.) fall
# through with a @warn — still writes empty output tapes so downstream
# modules don't SystemError on missing files.

"""
    acer_module(tapes::TapeManager, params::AcerParams)

Run ACER: read PENDF tape, build ACE table, write ACE file + directory.

Input tapes:
  - `params.nendf` : ENDF tape (optional — iopt≠1 may need it)
  - `params.npendf`: PENDF tape (main input for iopt=1)

Output tapes:
  - `params.nace` : ACE data file (Type-1 ASCII)
  - `params.ndir` : MCNP xsdir directory line

iopt=1 (fast continuous-energy neutron) is the only mode currently
producing real output. Other iopts get a stub empty tape so downstream
crashes on missing files don't cascade.
"""
function acer_module(tapes::TapeManager, params::AcerParams)
    @info "acer: iopt=$(params.iopt) nendf=$(params.nendf) npendf=$(params.npendf) ngend=$(params.ngend) nace=$(params.nace) ndir=$(params.ndir) suffix=$(params.suffix)"

    if params.iopt == 7
        return _acer_iopt7(tapes, params)
    elseif params.iopt != 1
        @warn "acer: iopt=$(params.iopt) not yet implemented — writing empty stub tapes"
        _acer_write_empty(tapes, params)
        return nothing
    end

    # iopt=1: fast continuous-energy / charged-particle ACE
    pendf_path = params.npendf > 0 ? resolve(tapes, params.npendf) : ""
    if !isfile(pendf_path)
        @warn "acer: PENDF input (unit $(params.npendf)) not found — writing empty stub tapes"
        _acer_write_empty(tapes, params)
        return nothing
    end

    # Read PENDF and extract all MF3 sections for the target material
    pendf = read_pendf(pendf_path)
    mat   = params.mat
    mf3   = extract_mf3_all(pendf, mat)

    if isempty(mf3)
        @warn "acer: no MF3 data for MAT=$mat on tape $(params.npendf) — stub tapes"
        _acer_write_empty(tapes, params)
        return nothing
    end

    # MT=1 provides the master energy grid (reconr/broadr output has all
    # partials on the same lunion grid). If MT=1 is missing, fall back to
    # the longest MT's grid.
    master_mt = haskey(mf3, 1) ? 1 : argmax(Dict(k => length(v[1]) for (k, v) in mf3))
    master_e  = mf3[master_mt][1]

    # Read MF1/MT451 header (needed early to gate the charged-particle ESZ
    # path on NSUB — Fortran's `if (izai.eq.1)` at acefc.f90:1190 vs the
    # `else` branch at line 1538 that runs the unionx ICP pipeline).
    endf_path = params.nendf > 0 ? resolve(tapes, params.nendf) : pendf_path
    mf1 = try
        read_mf1_mt451_header(endf_path, mat)
    catch err
        @warn "acer: read_mf1_mt451_header failed: $err — falling back to PENDF MF1"
        (za = _acer_extract_za(pendf, mat), awr = NaN, nsub = 10,
         lrp = 0, lfi = 0, nlib = 0, nmod = 0, elis = 0.0, sta = 0.0,
         lis = 0, liso = 0, nfor = 6, awi = 0.0, emax = 2e7, lrel = 0, nver = 0)
    end

    # For incident-charged-particle ACE, run the Fortran unionx pipeline:
    # MF3 union → MF6/MT2 anchor merge (with the c-buffer overwrite quirk
    # that drops one MF3 point per anchor and creates duplicates) → step-1.2
    # ratio enforcement (off-by-one drops second-to-last) → gety1-style
    # dedup at aceout. Ref: njoy-reference/src/acefc.f90:1538-1693 (unionx
    # ICP path) and acefc.f90:5341-5355 (aceout ESZ readback via gety1).
    # NSUB=10 → neutron incident, takes the simpler izai==1 path (no MF6
    # anchor merge); any other NSUB is a charged particle.
    izai_neutron = (mf1.nsub == 10)
    if !izai_neutron
        endf_for_mf6 = params.nendf > 0 ? resolve(tapes, params.nendf) : pendf_path
        mf6_e = try
            read_mf6_incident_energies(endf_for_mf6, mat, 2)
        catch
            Float64[]
        end
        if !isempty(mf6_e)
            master_e = _acer_unionx_charged(master_e, sort(mf6_e))
        end
    end
    n_e = length(master_e)

    # Build mt_list in canonical order (1, 2, fission, capture, then others)
    mt_sorted = Int[]
    for mt_pref in (1, 2, 18, 102)
        haskey(mf3, mt_pref) && push!(mt_sorted, mt_pref)
    end
    for mt in sort(collect(keys(mf3)))
        mt in mt_sorted && continue
        push!(mt_sorted, mt)
    end

    # Assemble XS matrix. For neutron PENDFs from reconr/broadr, all partials
    # share the lunion grid — direct copy. For the charged-particle path
    # master_e was extended by unionx (MF6 + step pads), so MT=2 must be
    # linearly interpolated onto the new grid.
    xs_matrix = zeros(Float64, n_e, length(mt_sorted))
    for (col, mt) in enumerate(mt_sorted)
        e_mt, xs_mt = mf3[mt]
        if length(e_mt) == n_e && e_mt == master_e
            xs_matrix[:, col] = xs_mt
        else
            xs_matrix[:, col] = _acer_linear_interp(master_e, e_mt, xs_mt)
        end
    end

    # Construct the PointwiseMaterial that ace_builder expects
    pm = PointwiseMaterial(Int32(mat), master_e, xs_matrix, mt_sorted)

    # Derive suffix letter from NSUB. The parser gave a provisional letter
    # (default 'c'); replace it with the NSUB-driven one.
    letter = acer_incident_letter(mf1.nsub)
    # Strip trailing alpha char from params.suffix ("10c" → "10") then append
    # the derived letter. Preserves numeric part (temperature indicator).
    numeric_part = replace(params.suffix, r"[a-zA-Z]+$" => "")
    final_suffix = isempty(numeric_part) ? "00$letter" : "$(numeric_part)$letter"

    # Title from card 3; mat string in Fortran NJOY form "   mat%4d"
    comment = isempty(params.title) ? "" : params.title

    table = build_ace_from_pendf(pm;
                                  suffix=final_suffix,
                                  temp_kelvin=params.temp,
                                  comment=comment,
                                  mat_id=mat,
                                  za=mf1.za,
                                  awr=mf1.awr,
                                  date=_acer_today_date())

    # Write outputs
    if params.nace > 0
        ace_path = resolve(tapes, params.nace)
        open(ace_path, "w") do io
            write_ace_table(io, table)
        end
        register!(tapes, params.nace, ace_path)
        @info "acer: wrote ACE file $(ace_path) ($(countlines(ace_path)) lines)"
    end

    if params.ndir > 0
        dir_path = resolve(tapes, params.ndir)
        nxs, _, _, _ = build_xss(table)
        open(dir_path, "w") do io
            write_ace_directory(io, table, nxs;
                                 itype=1,
                                 filepath=params.nace > 0 ? basename(resolve(tapes, params.nace)) : "filename",
                                 route="route")
        end
        register!(tapes, params.ndir, dir_path)
        @info "acer: wrote xsdir $(dir_path)"
    end

    nothing
end

"""Write zero-byte placeholder tapes for unsupported iopt modes or missing
input. Prevents downstream SystemError when the next module in the chain
tries to open these units."""
function _acer_write_empty(tapes::TapeManager, params::AcerParams)
    for unit in (params.ngend, params.nace, params.ndir)
        unit <= 0 && continue
        path = resolve(tapes, unit)
        touch(path)
        register!(tapes, unit, path)
    end
end

"""
    _acer_iopt7(tapes, params)

iopt=7 — read existing Type-1 ACE from `npendf`, produce a copy on `nace`,
a viewr-format plot tape on `ngend`, and a 1-line summary record on `ndir`.

This is Fortran acer.f90's "check/view" mode. For now the nace output is a
verbatim passthrough of the input ACE (matches Fortran behaviour when no
thinning is requested). The ngend viewr plot tape is produced separately by
`_acer_write_viewr_plot_tape` — currently a minimal stub.

Ref: njoy-reference/src/acer.f90:449-465 (iopt=7 dispatch)
"""
function _acer_iopt7(tapes::TapeManager, params::AcerParams)
    src = params.npendf > 0 ? resolve(tapes, params.npendf) : ""
    if !isfile(src)
        @warn "acer iopt=7: input ACE tape $(params.npendf) not found — writing empty stubs"
        _acer_write_empty(tapes, params)
        return nothing
    end

    # nace slot: copy the input ACE verbatim. Fortran's iopt=7 re-formats
    # through aceout — for ASCII Type-1 input the result is byte-identical
    # modulo the date string in the header.
    if params.nace > 0
        dst = resolve(tapes, params.nace)
        cp(src, dst; force=true)
        register!(tapes, params.nace, dst)
        @info "acer iopt=7: copied ACE $(basename(src)) → $(basename(dst)) ($(countlines(dst)) lines)"
    end

    # ngend slot: viewr-format plot tape. Stub empty — plot generation is a
    # separate grind (Fortran acer.f90 `aplots`). Downstream viewr reading
    # this file will no-op on empty input.
    if params.ngend > 0
        plot_path = resolve(tapes, params.ngend)
        touch(plot_path)
        register!(tapes, params.ngend, plot_path)
    end

    # ndir slot: 1-line summary line (xsdir-style) for the ACE table.
    # Fortran acer iopt=7 writes one "  ZAID  AWR filename route 1   1   LENGTH 0 0 0.000E+00"
    # record for each ACE entry — same shape as the xsdir produced by iopt=1.
    if params.ndir > 0
        dir_path = resolve(tapes, params.ndir)
        _acer_write_ndir_from_ace(src, dir_path)
        register!(tapes, params.ndir, dir_path)
    end

    nothing
end

"""Write a 1-line xsdir-style summary record by reading the header + NXS[1]
from an ASCII Type-1 ACE file. Matches Fortran acer iopt=7 `trail` output."""
function _acer_write_ndir_from_ace(ace_path::AbstractString, out_path::AbstractString)
    lines = readlines(ace_path)
    if length(lines) < 7
        touch(out_path)
        return
    end
    # Line 1: ZAID(a10) AWR(f12.6) TEMP_MEV(1pe12.4) date(a10)  — 45 chars ish
    hdr1 = rpad(lines[1], 80)
    zaid = strip(hdr1[1:10])
    awr  = strip(hdr1[11:22])
    temp = strip(hdr1[23:35])
    # Line 7 begins the NXS row, whose first integer is XSS length (i9).
    nxs1 = strip(rpad(lines[7], 80)[1:9])
    open(out_path, "w") do io
        @printf(io, "%-10s %s filename route 1   1 %8s     0     0 %s\n",
                zaid, awr, nxs1, temp)
    end
end

"""Date string in Fortran NJOY aceout format: `MM/DD/YY`. The reference tapes
fix this to match the regeneration date; execute.py substitutes via
`datePattern = re.compile(r'\\d{2}/\\d{2}/\\d{2}')` before comparing, so any
valid date passes."""
function _acer_today_date()
    t = Base.Libc.TmStruct(time())
    @sprintf("%02d/%02d/%02d", t.month + 1, t.mday, (t.year + 1900) % 100)
end

"""Extract ZA from MF1/MT451 HEAD record of a PENDFTape. Returns 0 on miss."""
function _acer_extract_za(pendf, mat::Int)
    for material in pendf.materials
        material.mat != mat && continue
        for sec in material.sections
            sec.mf == 1 && sec.mt == 451 || continue
            isempty(sec.lines) && continue
            head = rpad(sec.lines[1], 80)
            return Int(round(parse_endf_float(head[1:11])))
        end
    end
    0
end

"""Rebuild the ACE table header with the correct ZAID string derived from
the MF1 ZA (which differs from the MAT number for most isotopes)."""
function _acer_retag_header(table::ACENeutronTable, za::Int, suffix::AbstractString)
    zaid_str = format_zaid(za, suffix)
    h = table.header
    new_header = ACEHeader(zaid=zaid_str, awr=h.aw0,
                            temp_mev=h.tz, comment=h.hk,
                            mat_string=h.hm)
    ACENeutronTable(header=new_header,
                    energy_grid=table.energy_grid,
                    total_xs=table.total_xs,
                    absorption_xs=table.absorption_xs,
                    elastic_xs=table.elastic_xs,
                    heating_numbers=table.heating_numbers,
                    reactions=table.reactions)
end

"""Simple linear interpolator (eV → eV) for the PENDF-grid-mismatch fallback."""
function _acer_linear_interp(query_e::Vector{Float64},
                              src_e::Vector{Float64}, src_xs::Vector{Float64})
    n = length(query_e)
    out = zeros(Float64, n)
    for (i, e) in enumerate(query_e)
        if e <= src_e[1]
            out[i] = src_xs[1]
        elseif e >= src_e[end]
            out[i] = src_xs[end]
        else
            j = searchsortedlast(src_e, e)
            f = (e - src_e[j]) / (src_e[j+1] - src_e[j])
            out[i] = src_xs[j] + f * (src_xs[j+1] - src_xs[j])
        end
    end
    out
end

"""
    _acer_unionx_charged(mf3_grid, mf6_anchors) -> Vector{Float64}

Build the ESZ energy grid for charged-particle ACE faithfully reproducing
Fortran `unionx` (acefc.f90:1538-1693) plus `aceout` ESZ readback
(acefc.f90:5341-5355).

The Fortran pipeline has three stages, each with subtle bugs/quirks that
affect the final point count and must be replicated for bit-identity:

1. **MF6/MT2 anchor merge** with c-buffer overwrite quirk (lines 1608-1656).
   The `c` array is overwritten by the explicit `c(1)=ee` write before the
   next inner-while iteration runs. Effect: the MF3 point that was just
   read into `c` (the one strictly greater than the anchor) is silently
   dropped, AND the next inner-while iteration writes `c` (now holding
   the anchor) again, creating a duplicate.

2. **Step-1.2 pass** (lines 1658-1686). The outer `do while (k.lt.nold)`
   has an off-by-one that drops the second-to-last point: panels are
   processed for k=2..nold-1, then a final `loada` writes old[nold] —
   old[nold-1] is never written.

3. **gety1 dedup at aceout** (line 5346-5354). When aceout reads MT=1
   back from the scratch tape via `gety1`, duplicate energies (which
   represent step discontinuities in TAB1 form) collapse to single
   entries.

For T50 (α+He-4): MF3=29 → after stage 1 = 32 (drops 4.4416, 13.578,
17.906; dups at 4.0, 12.9, 16.6) → after stage 2 = 39 (8 step pads
inserted; one 16.6 dup dropped) → after stage 3 = 37.
"""
function _acer_unionx_charged(mf3_grid::Vector{Float64},
                                mf6_anchors::Vector{Float64})
    isempty(mf6_anchors) && return copy(mf3_grid)
    isempty(mf3_grid)   && return copy(mf6_anchors)

    # ---- Stage 1: MF6 anchor merge with c-buffer overwrite quirk ------
    grid = Float64[]
    n_old = length(mf3_grid)
    k = 1               # 1-indexed cursor into mf3_grid (Fortran finda(k))
    eg = mf3_grid[k]    # current "next" old-grid value
    c_stale = mf3_grid[k]  # mirrors Fortran's `c` array — clobbered by ee writes
    n_anc = length(mf6_anchors)
    for (i, ee) in enumerate(mf6_anchors)
        # Inner while: advance through old grid until eg >= ee.
        while eg < ee && k <= n_old
            push!(grid, c_stale)        # write previous finda's value
            k += 1
            if k > n_old; break; end
            c_stale = mf3_grid[k]
            eg = c_stale
        end
        # Write the MF6 anchor (this overwrites c_stale!).
        push!(grid, ee)
        c_stale = ee
        # If eg == ee, advance k once (Fortran lines 1647-1651).
        if eg == ee && k <= n_old
            k += 1
            if k <= n_old
                c_stale = mf3_grid[k]
                eg = c_stale
            end
        end
        # On the LAST anchor, Fortran's loada(jt=-j,...) flushes; nothing
        # else gets written. Trailing MF3 points beyond the last anchor
        # are silently dropped — same behaviour we replicate here.
    end

    # ---- Stage 2: step-1.2 pass with off-by-one second-to-last drop ----
    n_pre = length(grid)
    if n_pre < 2
        post_step = copy(grid)
    else
        post_step = Float64[]
        # Process panels [grid[i], grid[i+1]] for i=1..n_pre-2.
        # (Fortran outer loop k=2..nold-1; the k=nold case is skipped.)
        for i in 1:(n_pre-2)
            eg_i = grid[i]
            er_i = grid[i+1]
            push!(post_step, eg_i)
            egrid_i = round_sigfig(1.2 * eg_i, 2, 0)
            while egrid_i < er_i
                push!(post_step, egrid_i)
                eg_i = egrid_i
                egrid_i = round_sigfig(1.2 * eg_i, 2, 0)
            end
        end
        # Final point (grid[n_pre-1] is dropped — Fortran off-by-one).
        push!(post_step, grid[n_pre])
    end

    # ---- Stage 3: gety1-style adjacent-duplicate dedup -----------------
    # gety1 in aceout collapses adjacent equal energies to a single entry.
    out = Float64[post_step[1]]
    for i in 2:length(post_step)
        if post_step[i] != post_step[i-1]
            push!(out, post_step[i])
        end
    end
    out
end
