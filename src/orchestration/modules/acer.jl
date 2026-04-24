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

    # Assemble XS matrix. All MTs on reconr-output PENDFs share the lunion
    # grid; log length mismatches and pad with interpolation (defensive).
    xs_matrix = zeros(Float64, n_e, length(mt_sorted))
    for (col, mt) in enumerate(mt_sorted)
        e_mt, xs_mt = mf3[mt]
        if length(e_mt) == n_e && e_mt == master_e
            xs_matrix[:, col] = xs_mt
        else
            @warn "acer: MT=$mt grid ($(length(e_mt)) pts) differs from master ($(n_e) pts) — linear interp"
            xs_matrix[:, col] = _acer_linear_interp(master_e, e_mt, xs_mt)
        end
    end

    # Construct the PointwiseMaterial that ace_builder expects
    pm = PointwiseMaterial(Int32(mat), master_e, xs_matrix, mt_sorted)

    # Build the ACE table
    temp_k = params.temp > 0 ? params.temp : 300.0
    zaid_mat = _acer_extract_za(pendf, mat)  # MF1/MT451 HEAD gives ZA for true ZAID
    table = build_ace_from_pendf(pm;
                                  suffix=params.suffix,
                                  temp_kelvin=temp_k,
                                  mat_id=mat)

    # If we extracted a real ZA from MF1, override the header's hz (the
    # builder defaulted to `mat` which for pipeline chains != ZA).
    if zaid_mat > 0 && zaid_mat != mat
        table = _acer_retag_header(table, zaid_mat, params.suffix)
    end

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
