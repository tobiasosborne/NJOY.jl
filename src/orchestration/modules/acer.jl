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

    # NSUB=10 → neutron incident; any other NSUB is a charged particle.
    izai_neutron = (mf1.nsub == 10)
    endf_for_mf6 = params.nendf > 0 ? resolve(tapes, params.nendf) : pendf_path
    if !izai_neutron
        # Charged-particle ACE: run the Fortran unionx pipeline —
        # MF3 union (ALL MF3 sections in tape order) → MF6/MT2 anchor merge
        # (with the c-buffer overwrite quirk that drops one MF3 point per
        # anchor and creates duplicates) → step-1.2 densification (off-by-one
        # drops second-to-last) → gety1-style dedup at aceout. Ref:
        # njoy-reference/src/acefc.f90:1541-1693 (unionx ICP path) and
        # :5341-5355 (aceout ESZ readback).
        #
        # CRITICAL: the base grid is the union of *every* MF3 reaction grid,
        # not just the longest single MT. T54 (proton+H-3) has no MT=1 and
        # three MF3 sections (MT2/50/650) whose union (plus the MF6/MT2
        # anchors) gives the 140-point ESZ grid; using only MT50's grid (the
        # longest) yields 125 points and a structurally wrong tape.
        mf6_e = try
            read_mf6_incident_energies(endf_for_mf6, mat, 2)
        catch
            Float64[]
        end
        # mf3 grids in ascending-MT (= tape) order, as Fortran processes them.
        mf3_grids = [mf3[mt][1] for mt in sort(collect(keys(mf3)))]
        master_e = _acer_unionx_charged(mf3_grids, sort(mf6_e))
    else
        # Neutron incident: union MF6 incident energies for ALL non-total
        # MTs into the master grid. This is the simpler (izai==1) path that
        # Fortran takes at acefc.f90:1190+ — every MF6 tabulation point
        # must land on an ESZ grid point so the LAND/AND pointers can
        # cross-reference. Pre-existing behaviour from before the Phase-5
        # charged-particle restructure; required by T08 (Ni-61 has MF6 for
        # several inelastic MTs whose incident energies aren't all in MT=1's
        # lunion grid).
        for mt in keys(mf3)
            mt == 1 && continue
            mf6_e = try
                read_mf6_incident_energies(endf_for_mf6, mat, mt)
            catch
                Float64[]
            end
            isempty(mf6_e) && continue
            master_e = sort(unique(vcat(master_e, mf6_e)))
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
    # share the lunion grid — direct copy.
    #
    # The charged-particle path is fundamentally different: master_e is the
    # union of MT=2's grid with the MF6 incident-energy anchors (unionx), so
    # it is FINER than each reaction's native MF3 grid, AND the charged-
    # particle MF3 sections retain non-linear ENDF interpolation (INT=5
    # log-log, INT=6 Coulomb penetrability) that reconr never linearised.
    # Fortran acelod (acefc.f90:5494-5505) samples each reaction onto the ESZ
    # grid via gety1 — i.e. terp1 with the panel's own law — then sigfig(s,7,0)
    # BEFORE storing into the SIG block / accumulating the disappear+total
    # columns. Linear interpolation here gave e.g. MT=600 at E=120 eV =
    # 6.08e-83 instead of the INT=6 value 4.07e-107 (T62 SIG-block tail bug).
    xs_matrix = zeros(Float64, n_e, length(mt_sorted))
    mf3_q = Dict{Int, Float64}()   # per-MT QI [eV] for the LQR column (charged path)
    mf3_tab_all = Dict()           # per-MT MF3 TAB1 (charged path; for production scan)
    # per-MT MF3 threshold (first tabulated energy, eV) — used to place the
    # SIG-block ie_start the Fortran way (acelod walks DOWN from nes to the
    # first grid point ≥ threshold, so the threshold zero is INCLUDED). This
    # is what makes T54's MT=650 start at grid index 118 (E=5.383 MeV, xs=0)
    # rather than 119 (first nonzero). Ref: acefc.f90:5594-5610.
    mf3_thresh = Dict{Int, Float64}()
    if !izai_neutron
        mf3_tab = extract_mf3_tab1_all(pendf, mat)
        mf3_tab_all = mf3_tab
        for (col, mt) in enumerate(mt_sorted)
            rec = get(mf3_tab, mt, nothing)
            if rec === nothing
                # Fallback: lin-lin on the bare (e,xs) (should not happen).
                e_mt, xs_mt = mf3[mt]
                xs_matrix[:, col] = _acer_linear_interp(master_e, e_mt, xs_mt)
                continue
            end
            mf3_q[mt] = rec.qi
            mf3_thresh[mt] = rec.tab.x[1]
            # gety1-faithful sample + sigfig(s,7,0). thr6=0 in the ACE path
            # (njoy never sets thr6 in acefc; the endf module global stays 0).
            for (i, e) in enumerate(master_e)
                s = interpolate(rec.tab, e; coulomb_threshold=0.0)
                xs_matrix[i, col] = round_sigfig(s, 7, 0)
            end
            # Force the first stored (threshold) point to zero — Fortran
            # acelod acefc.f90:5627-5628: `test=100; if (j>1 .and. n==0 .and.
            # enext>test) s=0`. The first SIG point of a threshold reaction
            # (threshold > 100 eV, and not the absolute first grid point) is
            # zeroed before it feeds the SIG block, the disappear/total
            # accumulation, AND the acelcp production xs. Without this, T54's
            # MT=50 stored xs[1] at grid[61]=1.02 MeV reads the interpolated
            # 0.023438 instead of the reference 0.0 (and the same value leaks
            # into absorption[61], total[61], and the neutron HPD).
            thr = rec.tab.x[1]
            eps = 1.0e-10
            ie0 = n_e
            while ie0 >= 1 && thr <= (1 + eps) * master_e[ie0]
                ie0 -= 1
            end
            ie0 += 1
            if ie0 > 1 && ie0 <= n_e && thr > 100.0
                xs_matrix[ie0, col] = 0.0
            end
        end
    else
        for (col, mt) in enumerate(mt_sorted)
            e_mt, xs_mt = mf3[mt]
            if length(e_mt) == n_e && e_mt == master_e
                xs_matrix[:, col] = xs_mt
            else
                xs_matrix[:, col] = _acer_linear_interp(master_e, e_mt, xs_mt)
            end
        end
    end

    # Construct the PointwiseMaterial that ace_builder expects
    pm = PointwiseMaterial(Int32(mat), master_e, xs_matrix, mt_sorted)

    # Charged-particle elastic LAW=5 → ACE LAW=14 (LAND/AND block + Coulomb-
    # corrected total/elastic). For neutron incident, Fortran takes a
    # different angular path (acensd at acefc.f90:5874-5877); for charged
    # particles, acecpe at acefc.f90:6492-6671 produces the tabulated
    # (μ, pdf, cdf) tables and replaces elastic_xs / total_xs with the
    # Coulomb-included values via log-log interpolation.
    charged_elastic = nothing
    if !izai_neutron
        ce_endf = params.nendf > 0 ? resolve(tapes, params.nendf) : pendf_path
        ce_hdr, ce_subs = try
            read_mf6_law5(ce_endf, mat, 2)
        catch
            (nothing, MF6Law5Subsection[])
        end
        if !isempty(ce_subs)
            elastic_col = findfirst(==(2), mt_sorted)
            esz_elastic = elastic_col === nothing ? zeros(n_e) :
                            xs_matrix[:, elastic_col]
            # The ACE ESZ total column that acecpe reads is NOT the raw MT=1
            # XS — Fortran builds it by accumulating the *sigfig-7-rounded*
            # reaction cross sections (acefc.f90:5497-5505: MT=2 always, MT≥5
            # added), then sigfig-9 (acefc.f90:6304). For the MT=2-only
            # charged-particle case this is sigfig9(sigfig7(elastic)). The
            # post-loop subtracts the sigfig-7 elastic and adds signow
            # (acefc.f90:6659-6667); if `esz_total` were the raw elastic the
            # residual (raw − sigfig7) corrupts the total at the 9th sigfig
            # (T50 line-30 total 3.55810398 vs 3.5581040). Build the column
            # the Fortran way.
            esz_total = zeros(n_e)
            for (col, mt) in enumerate(mt_sorted)
                (mt == 2 || mt >= 5) || continue
                @views esz_total .+= round_sigfig.(xs_matrix[:, col], 7, 0)
            end
            esz_total = round_sigfig.(esz_total, 9, 0)
            izai = mf1.nsub ÷ 10  # nsub=20040 → izai=2004 (alpha)
            ang, new_total, new_elastic, new_heat = acer_charged_elastic(
                ce_subs, master_e, esz_total, esz_elastic,
                Float64(mf1.awr), Float64(mf1.awi), Int(izai), Int(mf1.za))
            charged_elastic = (angular = ang, total = new_total,
                                elastic = new_elastic, heating = new_heat)
        end
    end

    # Particle-production scan (acelcp). For charged incident, scan MF6 for
    # light-particle productions and build the nprod/iprod/... table + the
    # t201..t207 thresholds. NTYPE=0 (no production found) → no blocks appended
    # and the basic tape is byte-identical to the pre-acelcp output.
    particle_production = nothing
    if !izai_neutron
        izai_p = mf1.nsub ÷ 10
        ce_endf = params.nendf > 0 ? resolve(tapes, params.nendf) : pendf_path
        particle_production = try
            scan_particle_production(ce_endf, mat, Int(izai_p),
                                     Float64(mf1.awi), Float64(mf1.awr),
                                     Int(mf1.za), mf3_tab_all)
        catch err
            @warn "acer: particle-production scan failed: $err — no particle blocks"
            nothing
        end
        # If the scan found no productions, leave it nothing (gate stays closed).
        if particle_production !== nothing && isempty(particle_production.mprod)
            particle_production = nothing
        end
    end

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
                                  date=_acer_today_date(),
                                  charged_elastic=charged_elastic,
                                  incident_charged=!izai_neutron,
                                  mf3_q=mf3_q,
                                  mf3_thresh=mf3_thresh,
                                  particle_production=particle_production)

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

    # Determine the data-class letter from the input ZAID's 10th column so we
    # apply the right fix routine. ht='t' → thrfix (thermal), 'c'/'h'/'o'/'r'/
    # 's'/'a' → acefix (neutron/charged), etc. Ref: acer.f90:504-533.
    src_lines = readlines(src)
    isempty(src_lines) && (touch_all(tapes, params); return nothing)
    hz0   = rpad(src_lines[1], 80)[1:10]   # `read(npend,'(a10)') hz(1:10)`
    ht    = length(strip(hz0)) >= 1 ? hz0[10] : ' '

    # Apply the new ZAID suffix from card-2 `suff` (acer.f90:295 reads it as a
    # real). thrfix/acefix call newsuff only when suff >= 0 (aceth.f90:618);
    # suff < 0 (and our NaN "absent" sentinel) leaves the ZAID unchanged.
    new_hz = hz0
    if !isnan(params.suff) && params.suff >= 0.0
        new_hz = _acer_newsuff(hz0, params.suff)
    end

    # nace slot: rewrite the input ACE with the new ZAID in line 1. Fortran's
    # iopt=7 reads the ACE through thrfix/throut and re-emits it; for ASCII
    # Type-1 input the body is byte-identical, only the header ZAID changes.
    # Ref: aceth.f90:617-618 (newsuff) + 2314-2316 (throut header rewrite).
    if params.nace > 0
        dst = resolve(tapes, params.nace)
        out_lines = copy(src_lines)
        # Header line 1: replace columns 1-10 (the a10 ZAID field) with new_hz,
        # preserving the AWR/TEMP/date tail byte-for-byte.
        tail = length(out_lines[1]) > 10 ? out_lines[1][11:end] : ""
        out_lines[1] = new_hz * tail
        open(dst, "w") do io
            for ln in out_lines
                println(io, ln)
            end
        end
        register!(tapes, params.nace, dst)
        @info "acer iopt=7: rewrote ACE $(basename(src)) → $(basename(dst)) " *
              "($(length(out_lines)) lines, zaid $(strip(hz0))→$(strip(new_hz)))"
    end

    # ngend slot: viewr-format plot tape. Fortran acefix calls aplots only for
    # continuous/charged data classes (ht in {c,h,o,r,s,a}); thermal ('t') and
    # other classes never plot. Ref: acefc.f90:14198-14206 + 14205 (aplots).
    if params.ngend > 0
        plot_path = resolve(tapes, params.ngend)
        if ht in ('c', 'h', 'o', 'r', 's', 'a')
            # hk title = ACE header line 2, cols 1-70 (Fortran hko, a70).
            hk = length(src_lines) >= 2 ? rpad(src_lines[2], 70)[1:70] : ""
            table = read_ace_ascii(src_lines)
            # Render into a buffer first so a not-yet-ported aplots block leaves
            # the on-disk plot tape as the (intentionally deferred) empty stub
            # rather than a truncated partial — and never blocks tape35 below.
            try
                buf = IOBuffer()
                _acer_aplots(buf, table, hk, ht)
                write(plot_path, take!(buf))
                @info "acer iopt=7: wrote viewr plot tape $(basename(plot_path)) " *
                      "(class '$ht', nes=$(Int(table.nxs[NXS_NES])))"
            catch err
                err isa AplotsNotPortedError || rethrow()
                touch(plot_path)
                @warn "acer iopt=7: viewr plot tape $(basename(plot_path)) left as " *
                      "empty stub — $(err.msg)"
            end
        else
            # Non-plottable class: Fortran writes nothing to the plot unit.
            touch(plot_path)
        end
        register!(tapes, params.ngend, plot_path)
    end

    # ndir slot: 1-line xsdir directory record. The format differs by data
    # class: thermal (thrfix→throut, ht='t') uses the i2,' 1 ',i9 form
    # (aceth.f90:2422); neutron/fast (acefix→aceout) uses the i2,i4,1x,i8 form
    # (acefc.f90:13013). Match whichever class the input ACE is.
    if params.ndir > 0
        dir_path = resolve(tapes, params.ndir)
        _acer_write_ndir_from_ace(src_lines, new_hz, ht, dir_path)
        register!(tapes, params.ndir, dir_path)
    end

    nothing
end

"""touch all output tapes (used on the degenerate empty-input early-out)."""
function touch_all(tapes::TapeManager, params::AcerParams)
    for unit in (params.ngend, params.nace, params.ndir)
        unit <= 0 && continue
        p = resolve(tapes, unit)
        touch(p)
        register!(tapes, unit, p)
    end
end

"""
    _acer_newsuff(hz::AbstractString, suff::Float64) -> String

Port of Fortran `newsuff` (acecm.f90:619-684) for the MCNP (mcnpx=0) case:
update the 2-digit cross-section suffix of a 10-char ZAID `hz`, right-shifting
the class identifier into column 10. The new suffix digits are `nint(100*suff)`
written `i2.2` (zero-padded) into the two columns after the last `.`.

For `hz="  hh2o.99t"`, `suff=0.90` → `nint(90.0)=90` → `"  hh2o.90t"`.
"""
function _acer_newsuff(hz::AbstractString, suff::Float64)
    chars = collect(rpad(hz, 13))          # legacy njoy uses a 13-char buffer
    s10   = String(chars[1:13])
    lenhz = length(rstrip(String(chars)))  # len_trim(hz)
    # index(hz,".",.TRUE.) — last '.' (1-based)
    indx = findlast(==('.'), String(chars))
    if indx === nothing || indx < 1 || indx > 7
        # newsuff: "zaid name is nonstandard, no change made" (acecm.f90:662)
        return String(rpad(hz, 10)[1:10])
    end
    # write(hz(indx+1:indx+2),'(i2.2)') nint(100*suff)
    digits = @sprintf("%02d", round(Int, 100 * suff))
    chars[indx+1] = digits[1]
    chars[indx+2] = digits[2]
    # Right-shift so class id lands in column 10 (mcnpx=0 → indx=10).
    target = 10
    if lenhz < target
        idiff = target - lenhz
        for i in target:-1:(target - lenhz + 1)
            chars[i] = chars[i-idiff]
        end
        for i in 1:idiff
            chars[i] = ' '
        end
    end
    String(chars[1:10])
end

"""Write a 1-line xsdir-style directory record from an ASCII Type-1 ACE file's
header lines, using the new (suffix-applied) ZAID `new_hz`. Selects the thermal
vs neutron directory-line format by the data-class letter `ht`.

Thermal (ht='t', aceth.f90:2422):
  '(a10,f12.6,'' filename route'',i2,'' 1 '',i9,2i6,1p,e10.3)'
Neutron/fast (acefc.f90:13013):
  '(a10,f12.6,'' filename route'',i2,i4,1x,i8,2i6,1p,e10.3)'  (i4 = irec1 = 1)

For iopt=7 the ACE file is re-emitted unchanged, so lrec=nern=0 (throut sets
both to 0 for itype=1; aceout likewise for a freshly written ASCII file)."""
function _acer_write_ndir_from_ace(src_lines::Vector{String}, new_hz::AbstractString,
                                    ht::AbstractChar, out_path::AbstractString)
    if length(src_lines) < 7
        touch(out_path)
        return
    end
    # Header line 1: ZAID(a10) AWR(f12.6 / e12.0) TEMP(e12.0) date(a10).
    # throut re-emits aw0 with f12.6 and tz with 1pe10.3 — both come straight
    # from the input ACE header (read with e12.0 each; aceth.f90:580).
    hdr1 = rpad(src_lines[1], 80)
    aw0  = parse_endf_float(hdr1[11:22])
    tz   = parse_endf_float(hdr1[23:34])
    # Line 7 begins the NXS row; its first integer (i9) is len2 = XSS length.
    len2 = parse(Int, strip(rpad(src_lines[7], 80)[1:9]))
    itype = 1
    lrec  = 0
    nern  = 0
    zaid10 = rpad(new_hz, 10)[1:10]
    open(out_path, "w") do io
        if ht == 't'
            # aceth.f90:2422 — i2,' 1 ',i9,2i6,1p,e10.3
            @printf(io, "%s%12.6f filename route%2d 1 %9d%6d%6d%10.3E\n",
                    zaid10, aw0, itype, len2, lrec, nern, tz)
        else
            # acefc.f90:13013 — i2,i4,1x,i8,2i6,1p,e10.3  (i4 = irec1 = 1)
            @printf(io, "%s%12.6f filename route%2d%4d %8d%6d%6d%10.3E\n",
                    zaid10, aw0, itype, 1, len2, lrec, nern, tz)
        end
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
    _acer_unionx_charged(mf3_grids, mf6_anchors) -> Vector{Float64}

Build the ESZ energy grid for charged-particle ACE faithfully reproducing
Fortran `unionx` (acefc.f90:1541-1693) plus `aceout` ESZ readback
(acefc.f90:5341-5355).

`mf3_grids` is the list of per-reaction MF3 energy grids in tape (ascending-
MT) order; `mf6_anchors` is the sorted MF6/MT2 incident-energy list.

The Fortran pipeline has four stages, each with subtle quirks that affect
the final point count and must be replicated for bit-identity:

0. **MF3 union merge** (lines 1542-1606). Each MF3 reaction grid is merged,
   in tape order, into a running union via a two-cursor walk. The terminator
   `if (enext>test .and. k==nold) jt=-jt` fires only when BOTH cursors are
   exhausted, so the union spans the full energy range of all sections.

1. **MF6/MT2 anchor merge** with c-buffer overwrite quirk (lines 1620-1656).
   The `c` array is overwritten by the explicit `c(1)=ee` write before the
   next inner-while iteration runs. Effect: the next inner-while iteration's
   first `loada(j,c)` writes `c` (now holding the previous anchor `ee`)
   again, creating a duplicate that the stage-3 dedup later collapses. On
   the LAST anchor `jt=-jt` terminates the grid — trailing MF3 points above
   the last anchor are dropped (this is what caps T54 at 1.2e7 = the proton
   LAW=5 top energy, even though MT50's MF3 runs to 2e7).

2. **Step-1.2 densification** (lines 1659-1686). The outer `do while
   (k.lt.nold)` has an off-by-one that drops the second-to-last point:
   panels are processed for k=2..nold-1, then a final `loada` writes
   old[nold] — old[nold-1] is never written.

3. **gety1 dedup at aceout** (lines 5346-5354). When aceout reads MT=1 back
   from the scratch tape via `gety1`, adjacent duplicate energies collapse
   to single entries.

For T54 (proton+H-3): MF3 union (MT2∪MT50∪MT650) = 92 pts → +MF6/MT2 anchors
= 107 (capped at 1.2e7) → step-1.2 = 141 → dedup = 140.
For T62 (d+He-3): MF3 union (MT2∪MT600) = 226 → ... → 238.
"""
function _acer_unionx_charged(mf3_grids::Vector{Vector{Float64}},
                                mf6_anchors::Vector{Float64})
    etop = 1.0e10

    # ---- Stage 0: MF3 union merge (acefc.f90:1542-1606) ---------------
    # Merge a new section grid `ergrid` (its successive gety1 energies) into
    # the running union `old`. Two cursors (old via k/eg, new via ip/er),
    # writing the min each step; terminate when both exhausted.
    function merge_one(old::Vector{Float64}, ergrid::Vector{Float64})
        nold = length(old)
        out  = Float64[]
        k = 0; eg = etop
        if nold > 0; k = 1; eg = old[k]; end
        ne = length(ergrid); ip = 1
        er = ip <= ne ? ergrid[ip] : etop
        enext = er; jt = 1
        test = etop - etop / 100
        while jt > 0
            if er > eg
                if enext > test && k == nold; jt = -1; end
                push!(out, eg)
                if k < nold; k += 1; eg = old[k]; else; eg = etop; end
            else
                ip += 1
                enext = ip <= ne ? ergrid[ip] : etop
                if enext > test && k == nold; jt = -1; end
                push!(out, er)
                if er == eg && k < nold; k += 1; eg = old[k]; end
                er = enext
            end
        end
        out
    end

    grid = Float64[]
    for g in mf3_grids
        isempty(g) && continue
        grid = isempty(grid) ? copy(g) : merge_one(grid, g)
    end
    isempty(grid) && return copy(mf6_anchors)

    # ---- Stage 1: MF6/MT2 anchor merge with c-buffer overwrite --------
    if !isempty(mf6_anchors)
        n_old = length(grid)
        merged = Float64[]
        j = 0; k = 1
        c_stale = grid[k]   # mirrors Fortran `c` — clobbered by `c(1)=ee`
        eg = c_stale
        n_anc = length(mf6_anchors)
        for (iee, ee) in enumerate(mf6_anchors)
            while eg < ee
                j += 1; push!(merged, c_stale)       # loada(j,c)
                k += 1
                if k <= n_old; c_stale = grid[k]; eg = c_stale; else; eg = etop; end
            end
            j += 1
            jt = iee == n_anc ? -j : j
            c_stale = ee                              # c(1)=ee  ← clobbers c
            push!(merged, ee)                          # loada(jt,c)
            if eg == ee && jt > 0
                k += 1
                if k <= n_old; c_stale = grid[k]; eg = c_stale; else; eg = etop; end
            end
            # On the last anchor jt<0 → Fortran loop ends; trailing old
            # points (above the last anchor) are dropped.
        end
        grid = merged
    end

    # ---- Stage 2: step-1.2 densification with off-by-one drop ---------
    n_old = length(grid)
    if n_old >= 2
        out = Float64[]
        k = 1; eg = grid[k]
        while k < n_old
            k += 1
            cval = grid[k]
            if k < n_old
                er = cval
                egrid = 0.0
                while egrid < er
                    push!(out, eg)
                    egrid = round_sigfig(1.2 * eg, 2, 0)
                    if egrid < er; eg = egrid; end
                end
                eg = er
            end
        end
        # Final loada(jt=-j,...) writes grid[nold]; grid[nold-1] dropped.
        push!(out, grid[n_old])
        grid = out
    end

    # ---- Stage 3: gety1-style adjacent-duplicate dedup ----------------
    isempty(grid) && return grid
    dedup = Float64[grid[1]]
    for i in 2:length(grid)
        grid[i] != grid[i-1] && push!(dedup, grid[i])
    end
    dedup
end
