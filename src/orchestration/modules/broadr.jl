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
    # We re-run reconr's reader to get the section metadata
    reconr_mf3_secs = _read_mf3_section_metadata(endf_path, params.mat)

    # Resolve thnmax
    thnmax = resolve_thnmax(params.thnmax, reconr_mf3_secs, awr)
    @info "broadr: MAT=$(params.mat) thnmax=$(thnmax) ntemp=$(length(params.temperatures))"

    # Get energy grid from MT=1 (total)
    haskey(mf3_data, 1) || error("broadr: MT=1 (total) not found on input PENDF")
    energies = mf3_data[1][1]

    # Select partials for convergence testing (Trap 53: NOT total)
    xs_partials, partial_mts = select_broadr_partials(mf3_data)

    # Sequential broadening for each temperature
    cur_energies = energies
    cur_xs = xs_partials
    cur_total = mf3_data[1][2]

    for (it, temp) in enumerate(params.temperatures)
        alpha = awr / (PhysicsConstants.bk * temp)
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
                # Above thnmax: interpolate from previous grid
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

        @info "broadr: T=$(temp)K → $(length(b_e)) points"
    end

    # Build modified MF3 sections
    modified_mf3 = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
    modified_mf3[1] = (cur_energies, cur_total)
    for (col, mt) in enumerate(partial_mts)
        modified_mf3[mt] = (cur_energies, cur_xs[:, col])
    end

    # Write output PENDF: copy input with broadened sections replaced
    pendf_out = copy_with_modifications(pendf_in, params.mat;
        modified_mf3=modified_mf3,
        temperature=params.temperatures[end])
    write_pendf_tape(pendf_out_path, pendf_out)

    @info "broadr: wrote $pendf_out_path"

    # Return raw broadened data (for downstream modules that need full precision)
    return (energies=cur_energies, total=cur_total, partials=cur_xs, partial_mts=partial_mts)
end

# =========================================================================
# Helpers
# =========================================================================

"""Extract AWR from the first MF3 HEAD record in a PENDFTape."""
function _extract_awr_from_pendf(tape::PENDFTape, mat::Int)
    for material in tape.materials
        material.mat != mat && continue
        for sec in material.sections
            sec.mf != 3 && continue
            !isempty(sec.lines) || continue
            p = rpad(sec.lines[1], 80)
            return parse_endf_float(p[12:22])  # AWR is field 2 of HEAD
        end
    end
    return 0.0
end

"""
Read MF3 section metadata (QI, QM, AWR, etc.) from an ENDF file.
Used by resolve_thnmax to compute thresholds.
"""
function _read_mf3_section_metadata(endf_path::String, mat::Int)
    # Use the existing read_mf3_sections from reconr_types.jl
    io = open(endf_path, "r")
    try
        return read_mf3_sections(io, mat)
    finally
        close(io)
    end
end
