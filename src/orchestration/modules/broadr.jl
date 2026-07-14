# broadr module runner -- Doppler broaden pointwise cross sections
#
# Reads PENDF from reconr, broadens primary partials (not total),
# writes broadened PENDF. Matches Fortran broadr.f90 interface.

"""
    broadr_module(tapes::TapeManager, params::BroadrParams)

Run BROADR: read ENDF + PENDF tapes, Doppler-broaden cross sections,
write broadened PENDF tape.

Input tapes: nendf (ENDF), npendf_in (reconr PENDF)
Output tapes: npendf_out (broadened PENDF)
"""
function broadr_module(tapes::TapeManager, params::BroadrParams)
    endf_path = resolve(tapes, params.nendf)
    pendf_in_path = resolve(tapes, params.npendf_in)
    pendf_out_path = resolve(tapes, params.npendf_out)

    # Read input PENDF
    pendf_in = read_pendf(pendf_in_path)
    mf3_data = extract_mf3_all(pendf_in, params.mat)

    # Extract AWR from MF3 HEAD record or MF1
    awr = _extract_awr_from_pendf(pendf_in, params.mat)
    awr <= 0 && error("broadr: could not determine AWR for MAT=$(params.mat)")

    # Read reconr result for mf3_sections (needed by resolve_thnmax)
    reconr_mf3_secs = _read_mf3_section_metadata(endf_path, params.mat)

    # Extract eresh (resolved resonance range upper energy) from MF2/MT151.
    # Fortran caps it at e6pt5 = 6.5e6 eV (broadr.f90:423-424).
    eresh = _extract_eresh_from_pendf(pendf_in, params.mat)
    eresh > 6.5e6 && (eresh = 6.5e6)

    # Read original ENDF MF1/MT451 HEAD for LRP and the ENDF format version
    # (iverf). The PENDF copy has LRP overwritten to 2 by reconr
    # (reconr.f90:5031), so we MUST read the original ENDF here.
    # Ref: broadr.f90:249 (lrp=l1h) and :363-368 (iverf from NX/N1/N2 probe).
    hdr = read_mf1_header_info(endf_path, params.mat)
    lrp = hdr === nothing ? 0 : hdr.lrp
    iverf = hdr === nothing ? 6 : hdr.iverf

    # emin: 1.0 for non-resonance (lrp==0); eresh for resolved (lrp==1).
    # Ref: broadr.f90:356 (emin=1) and :432 (if lrp==1 emin=eresh).
    emin = lrp == 1 ? eresh : 1.0

    # Resolve thnmax
    thnmax = resolve_thnmax(params.thnmax, reconr_mf3_secs, awr; eresh=eresh)
    @info "broadr: MAT=$(params.mat) thnmax=$(thnmax) lrp=$(lrp) iverf=$(iverf) eresh=$(eresh) ntemp=$(length(params.temperatures))"

    # Dosimetry / IRDFF evaluations (e.g. Fe-nat IRDFF-II) carry only
    # partial reactions in MF3 — no MT=1 (total). Fortran broadr handles
    # this by broadening each partial independently; our current pipeline
    # is built around an MT=1 master grid. STUB: pass the input PENDF
    # through unchanged so downstream modules (groupr, errorr) see
    # usable data.
    if !haskey(mf3_data, 1)
        @warn "broadr: MT=1 not present on input PENDF (MAT=$(params.mat)) — " *
              "dosimetry pass-through (broadening skipped)"
        cp(pendf_in_path, pendf_out_path; force=true)
        register!(tapes, params.npendf_out, pendf_out_path)
        return nothing
    end
    energies = mf3_data[1][1]

    # LRF7 channel-partial MT list (Fortran nppmt, broadr.f90:308-330). Empty
    # for non-LRF7 materials. These reactions are broadened even when their MF3
    # threshold exceeds emin (broadr.f90:513-519).
    nppmt = _extract_nppmt_from_endf(endf_path, params.mat)

    # Select the mtr[] partials and project them onto the MT=1 master grid.
    # Ref: broadr.f90:500-524 (membership) + :620 (shared union grid).
    xs_partials, partial_mts = select_broadr_partials(
        mf3_data, energies; emin=emin, lrp=lrp, eresh=eresh,
        iverf=iverf, nppmt=nppmt)
    @info "broadr: selected broadened partials MT=$(partial_mts)"

    # Sequential broadening for each temperature — collect ALL results
    cur_energies = energies
    cur_xs = xs_partials
    cur_total = mf3_data[1][2]
    all_temp_results = NamedTuple[]
    t_old = 0.0  # previous temperature for sequential broadening

    for (it, temp) in enumerate(params.temperatures)
        t_eff = temp - t_old  # Fortran broadr uses T_eff = T_new - T_old

        # T=0K (or already-broadened T_old ≥ T_new): broadening kernel is a
        # delta function — pass through the current grid unchanged. Fortran
        # broadr handles this by short-circuiting the broadn call.
        if t_eff <= 0.0
            @info "broadr: T=$(temp)K ΔT=$(t_eff) ≤ 0 — pass-through (no broadening)"
            push!(all_temp_results, (
                temperature = temp,
                energies = copy(cur_energies),
                total = copy(cur_total),
                partials = copy(cur_xs),
                partial_mts = partial_mts,
            ))
            t_old = temp
            continue
        end

        alpha = awr / (PhysicsConstants.bk * t_eff)
        tol = params.tol
        errmax = 10 * tol
        errint = tol / 20000

        @info "broadr: T=$(temp)K alpha=$(round(alpha, sigdigits=5))"

        # Broaden the mtr[] partials TOGETHER on one shared union grid (b_e).
        # Ref: broadr.f90:620 (bfile3 broadens all nreac reactions at once).
        b_e, b_xs = broadn_grid(cur_energies, cur_xs, alpha, tol, errmax, errint, thnmax)

        # Reconstruct MT=1 (total) on the broadened grid.
        # Ref: broadr.f90:938-980 (label 275, MT=1 / ireac=0 path).
        #   - Thermal range, E ≤ thnmax (f90:944-973): MT1 = Σ of the
        #     individually-broadened partials  scr(k)=Σ tt(1+i)  (f90:974),
        #     summed from the UNROUNDED broadened partials, with the iflag
        #     channel-redundancy dedup (f90:944-973): when MT103 is in the set,
        #     drop its summed-into partials MT600-649; MT104→650-699;
        #     MT105→700-749; MT106→750-799; MT107→800-849 (iverf≥6 ranges,
        #     f90:457-467; iverf<6 ranges f90:469-478). sigfig is applied ONLY
        #     AFTER the sum (f90:980). The Doppler kernel is nonlinear so
        #     sigma1(Σpᵢ) ≠ Σ sigma1(pᵢ); broadening the total directly is wrong
        #     (bead e5n — gave ~5 b vs the correct ~1.93e5 b for B-10).
        #   - Above thnmax (f90:939-940): MT1 is the ORIGINAL (unbroadened) total
        #     read verbatim via gety1 — mirrored here by lin-lin interpolation of
        #     cur_total (the original total carried through the loop).
        #
        # Determine the iflag channel-dedup ranges from the version (broadr.f90
        # :457-478) and which MT103-107 sums are present in partial_mts.
        mpmin, mpmax, mdmin, mdmax, mtmin, mtmax, m3min, m3max, m4min, m4max =
            iverf >= 6 ?
                (600,649, 650,699, 700,749, 750,799, 800,849) :
                (700,718, 720,738, 740,758, 760,768, 780,798)
        has103 = 103 in partial_mts; has104 = 104 in partial_mts
        has105 = 105 in partial_mts; has106 = 106 in partial_mts
        has107 = 107 in partial_mts
        function _iflag(mt::Int)::Bool
            mt >= 201 && return true                         # f90:957
            has103 && mpmin <= mt <= mpmax && return true    # f90:959-961
            has104 && mdmin <= mt <= mdmax && return true    # f90:962-964
            has105 && mtmin <= mt <= mtmax && return true    # f90:965-967
            has106 && m3min <= mt <= m3max && return true    # f90:968-970
            has107 && m4min <= mt <= m4max && return true    # f90:971-973
            return false
        end
        sum_cols = [j for j in axes(b_xs, 2) if !_iflag(partial_mts[j])]

        b_total = Vector{Float64}(undef, length(b_e))
        for i in eachindex(b_e)
            if b_e[i] <= thnmax
                # f90:974 — sum the UNROUNDED broadened partials (dedup applied)
                s = 0.0
                for j in sum_cols
                    s += b_xs[i, j]
                end
                b_total[i] = s
            else
                # f90:940 — verbatim original total (gety1 lin-lin interp)
                idx = searchsortedfirst(cur_energies, b_e[i])
                if idx <= 1
                    b_total[i] = cur_total[1]
                elseif idx > length(cur_energies)
                    b_total[i] = cur_total[end]
                else
                    f = (b_e[i] - cur_energies[idx-1]) /
                        (cur_energies[idx] - cur_energies[idx-1])
                    b_total[i] = cur_total[idx-1] + f * (cur_total[idx] - cur_total[idx-1])
                end
            end
        end

        # sigfig AFTER the sum (f90:980 for MT1; Trap 120 for each partial).
        # MT1 sum above used the UNROUNDED b_xs, so sigfig the partials now.
        for j in axes(b_xs, 2), i in axes(b_xs, 1)
            b_xs[i, j] = round_sigfig(b_xs[i, j], 7, 0)
        end
        for i in eachindex(b_total)
            b_total[i] = round_sigfig(b_total[i], 7, 0)
        end

        cur_energies = b_e
        cur_xs = b_xs
        cur_total = b_total
        t_old = temp

        # Collect this temperature's results
        push!(all_temp_results, (
            temperature = temp,
            energies = copy(cur_energies),
            total = copy(cur_total),
            partials = copy(cur_xs),
            partial_mts = partial_mts,
        ))

        @info "broadr: T=$(temp)K → $(length(b_e)) points"
    end

    # Write multi-temperature PENDF (one MAT block per temperature)
    write_broadr_pendf(pendf_out_path, pendf_in, params.mat,
                       all_temp_results, params.tol)

    @info "broadr: wrote $pendf_out_path"

    return (all_temps=all_temp_results,
            energies=cur_energies, total=cur_total,
            partials=cur_xs, partial_mts=partial_mts)
end

# =========================================================================
# Helpers
# =========================================================================

"""Extract eresh (resolved resonance range EH) from MF2/MT151 in a PENDFTape."""
function _extract_eresh_from_pendf(tape::PENDFTape, mat::Int)
    for material in tape.materials
        material.mat != mat && continue
        for sec in material.sections
            (sec.mf == 2 && sec.mt == 151) || continue
            length(sec.lines) >= 3 || continue
            p = rpad(sec.lines[3], 80)
            return parse_endf_float(p[12:22])
        end
    end
    return Inf
end

"""Extract AWR from the first MF3 HEAD record in a PENDFTape."""
function _extract_awr_from_pendf(tape::PENDFTape, mat::Int)
    for material in tape.materials
        material.mat != mat && continue
        for sec in material.sections
            sec.mf != 3 && continue
            !isempty(sec.lines) || continue
            p = rpad(sec.lines[1], 80)
            return parse_endf_float(p[12:22])
        end
    end
    return 0.0
end

"""Read MF3 section metadata from an ENDF file."""
function _read_mf3_section_metadata(endf_path::String, mat::Int)
    io = open(endf_path, "r")
    try
        return read_mf3_sections(io, mat)
    finally
        close(io)
    end
end

"""
Extract the LRF7 channel-partial MT list (Fortran `nppmt`, broadr.f90:308-330)
from the ENDF MF2/MT151 R-Matrix-Limited (LRU=1, LRF=7) range. Returns the set
of particle-pair MTs, or an empty set for non-LRF7 materials.
"""
function _extract_nppmt_from_endf(endf_path::String, mat::Int)
    isfile(endf_path) || return Set{Int}()
    io = open(endf_path, "r")
    try
        rml = read_rml_data(io)
        rml === nothing && return Set{Int}()
        return Set{Int}(pp.MT for pp in rml.pairs)
    catch e
        @warn "broadr: failed to read LRF7 nppmt list — assuming none" exception=e
        return Set{Int}()
    finally
        close(io)
    end
end

# =========================================================================
# Multi-temperature PENDF writer for broadr
# =========================================================================

"""
    write_broadr_pendf(path, pendf_in, mat, temp_results, tol)

Write a multi-temperature PENDF tape matching Fortran broadr output.
Each temperature produces a complete MAT block (MF1 + MF2 + MF3 + MEND).
"""
function write_broadr_pendf(path::AbstractString, pendf_in::PENDFTape,
                            mat::Int, temp_results::AbstractVector,
                            tol::Float64)
    open(path, "w") do io
        src_mat = nothing
        for m in pendf_in.materials
            m.mat == mat && (src_mat = m; break)
        end
        src_mat === nothing && error("broadr: MAT=$mat not found in input PENDF")

        za = 0.0; awr = 0.0
        for sec in src_mat.sections
            sec.mf == 3 || continue
            isempty(sec.lines) && continue
            p = rpad(sec.lines[1], 80)
            za = parse_endf_float(p[1:11])
            awr = parse_endf_float(p[12:22])
            break
        end

        descriptions = _extract_broadr_descriptions(src_mat.mf1_lines)

        mf2_sections = Tuple{Int, Vector{String}}[]
        for sec in src_mat.sections
            sec.mf == 2 && push!(mf2_sections, (sec.mt, sec.lines))
        end

        mf3_order = Int[]
        mf3_lines = Dict{Int, Vector{String}}()
        for sec in src_mat.sections
            sec.mf == 3 || continue
            push!(mf3_order, sec.mt)
            mf3_lines[sec.mt] = sec.lines
        end
        # BROADR copies every non-MF3 section unchanged (apart from MF2/MT152's
        # temperature update above). Ref: broadr.f90:862-881,989-999.
        passthrough = [sec for sec in src_mat.sections if sec.mf != 2 && sec.mf != 3]

        broadened_mts = Set{Int}()
        if !isempty(temp_results)
            push!(broadened_mts, 1)
            for mt in temp_results[1].partial_mts
                push!(broadened_mts, mt)
            end
        end

        _write_tpid_line(io, pendf_in.tpid, mat)

        for tr in temp_results
            temp = tr.temperature
            b_np = length(tr.energies)

            dir_entries = Tuple{Int,Int,Int,Int}[]
            for (mt2, lines2) in mf2_sections
                push!(dir_entries, (2, mt2, _count_data_lines(lines2), 0))
            end
            for mt3 in mf3_order
                nc = mt3 in broadened_mts ? 3 + cld(b_np, 3) : _count_data_lines(mf3_lines[mt3])
                push!(dir_entries, (3, mt3, nc, 0))
            end
            for sec in passthrough
                push!(dir_entries, (sec.mf, sec.mt, _count_data_lines(sec.lines), 0))
            end
            nxc = length(dir_entries) + 1
            nwd = length(descriptions)

            ns = Ref(1)
            self_nc = 2 + nwd + nxc
            _write_cont_line(io, za, awr, 3, 1, 0, nxc, mat, 1, 451, ns)
            _write_cont_line(io, temp, tol, 0, 0, nwd, nxc, mat, 1, 451, ns)
            for desc in descriptions
                line = rpad(desc, 66)[1:66]
                @printf(io, "%s%4d%2d%3d%5d\n", line, mat, 1, 451, ns[])
                ns[] += 1
            end
            @printf(io, "%22s%11d%11d%11d%11d%4d%2d%3d%5d\n",
                    "", 1, 451, self_nc, 0, mat, 1, 451, ns[])
            ns[] += 1
            for (mf_d, mt_d, nc_d, mod_d) in dir_entries
                @printf(io, "%22s%11d%11d%11d%11d%4d%2d%3d%5d\n",
                        "", mf_d, mt_d, nc_d, mod_d, mat, 1, 451, ns[])
                ns[] += 1
            end
            _write_send(io, mat, 1)
            _write_fend_zero(io, mat)

            for (mt2, lines2) in mf2_sections
                ns[] = 1
                n_copy = _count_data_lines(lines2)
                for li in 1:n_copy
                    p = rpad(lines2[li], 80)
                    if mt2 == 152 && li == 2
                        p = format_endf_float(temp) * p[12:end]
                    end
                    @printf(io, "%s%4d%2d%3d%5d\n", p[1:66], mat, 2, mt2, ns[])
                    ns[] += 1
                end
                _write_send(io, mat, 2)
            end
            _write_fend_zero(io, mat)

            for mt3 in mf3_order
                ns[] = 1
                if mt3 in broadened_mts
                    xs_data = mt3 == 1 ? tr.total :
                        tr.partials[:, findfirst(==(mt3), tr.partial_mts)]
                    src_lines = mf3_lines[mt3]
                    ph = rpad(src_lines[1], 80)
                    l2_head = _broadr_parse_int(ph, 4)
                    pt = rpad(src_lines[2], 80)
                    qm = parse_endf_float(pt[1:11])
                    qi = parse_endf_float(pt[12:22])
                    lr = _broadr_parse_int(pt, 4)

                    _write_cont_line(io, za, awr, 0, l2_head, 0, 0, mat, 3, mt3, ns)
                    _write_cont_line(io, qm, qi, 0, lr, 1, b_np, mat, 3, mt3, ns)
                    # Interpolation record: NBT=NP, INT=2 — as integers matching Fortran tab1io
                    @printf(io, "%11d%11d%44s%4d%2d%3d%5d\n",
                            b_np, 2, "", mat, 3, mt3, ns[])
                    ns[] += 1
                    data = Float64[]
                    sizehint!(data, 2 * b_np)
                    for i in 1:b_np
                        push!(data, tr.energies[i])
                        push!(data, xs_data[i])
                    end
                    _write_data_values(io, data, mat, 3, mt3, ns; pair_data=true)
                else
                    src = mf3_lines[mt3]
                    n_copy = _count_data_lines(src)
                    for li in 1:n_copy
                        p = rpad(src[li], 80)
                        @printf(io, "%s%4d%2d%3d%5d\n", p[1:66], mat, 3, mt3, ns[])
                        ns[] += 1
                    end
                end
                _write_send(io, mat, 3)
            end
            _write_fend_zero(io, mat)

            # Copy MF10/MF12/MF13/etc. verbatim, regenerating only identifiers
            # and sequence numbers. MF14 is absent because RECONR omitted it.
            copied_mf = 0
            for sec in passthrough
                copied_mf != 0 && sec.mf != copied_mf && _write_fend_zero(io, mat)
                copied_mf = sec.mf
                ns[] = 1
                for li in 1:_count_data_lines(sec.lines)
                    p = rpad(sec.lines[li], 80)
                    @printf(io, "%s%4d%2d%3d%5d\n",
                            p[1:66], mat, sec.mf, sec.mt, ns[])
                    ns[] += 1
                end
                _write_send(io, mat, sec.mf)
            end
            copied_mf != 0 && _write_fend_zero(io, mat)
            # MEND (MAT=0): end of material block for this temperature
            @printf(io, "%66s%4d%2d%3d%5d\n", "", 0, 0, 0, 0)
        end

        # TEND (MAT=-1): end of tape
        @printf(io, "%66s%4d%2d%3d%5d\n", "", -1, 0, 0, 0)
    end
end

"""Extract description text lines from MF1/MT451 raw lines."""
function _extract_broadr_descriptions(mf1_lines::Vector{String})
    length(mf1_lines) < 2 && return String[]
    p = rpad(mf1_lines[2], 80)
    nwd = _broadr_parse_int(p, 5)
    nwd <= 0 && return String[]
    descs = String[]
    for i in 1:nwd
        idx = 2 + i
        idx > length(mf1_lines) && break
        push!(descs, rpad(mf1_lines[idx], 66)[1:66])
    end
    return descs
end

"""Count data lines in a stored section, excluding the trailing SEND record."""
function _count_data_lines(lines::Vector{String})
    n = length(lines)
    n == 0 && return 0
    p = rpad(lines[end], 80)
    mt = tryparse(Int, strip(p[73:75]))
    (mt !== nothing && mt == 0) ? n - 1 : n
end

"""Parse integer field (3=L1, 4=L2, 5=N1, 6=N2) from an 80-col ENDF CONT record."""
function _broadr_parse_int(line::AbstractString, field::Int)
    start = 22 + (field - 3) * 11 + 1
    stop = start + 10
    stop > length(line) && return 0
    s = strip(line[start:stop])
    isempty(s) && return 0
    return parse(Int, s)
end
