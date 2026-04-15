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

    # Extract eresh (resolved resonance range upper energy) from MF2/MT151
    eresh = _extract_eresh_from_pendf(pendf_in, params.mat)

    # Resolve thnmax
    thnmax = resolve_thnmax(params.thnmax, reconr_mf3_secs, awr; eresh=eresh)
    @info "broadr: MAT=$(params.mat) thnmax=$(thnmax) ntemp=$(length(params.temperatures))"

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

    # Select partials for convergence testing (Trap 53: NOT total)
    xs_partials, partial_mts = select_broadr_partials(mf3_data)

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

        # Broaden partials via convergence stack
        b_e, b_xs = broadn_grid(cur_energies, cur_xs, alpha, tol, errmax, errint, thnmax)

        # Broaden total separately on the same grid
        b_total = Vector{Float64}(undef, length(b_e))
        for i in eachindex(b_e)
            if b_e[i] <= thnmax
                b_total[i] = sigma1_at(b_e[i], cur_energies, cur_total, alpha)
            else
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

        # Apply round_sigfig(x, 7, 0) to broadened XS (Trap 120)
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
