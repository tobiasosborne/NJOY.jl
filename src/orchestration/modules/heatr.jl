# heatr module runner -- Compute KERMA coefficients and damage energy
#
# Reads ENDF supplementary data (MF4, MF5, MF12, MF13), extracts XS from
# input PENDF, calls compute_kerma, adds MT=301/444 to output PENDF.
# Matches Fortran heatr.f90 interface.

"""
    heatr_module(tapes::TapeManager, params::HeatrParams; reconr_result=nothing)

Run HEATR: read ENDF + PENDF tapes, compute heating KERMA and damage energy,
write PENDF with MT=301 (total KERMA) and MT=444 (damage energy) added.

Input tapes: nendf (ENDF), npendf_in (broadened PENDF)
Output tapes: npendf_out (PENDF with KERMA)

The `reconr_result` kwarg provides the reconr NamedTuple (needed for
`_get_legacy_section` to get threshold-adjusted XS). If not provided,
reconr is re-run from the ENDF tape.
"""
function heatr_module(tapes::TapeManager, params::HeatrParams;
                      reconr_result=nothing, broadened_data=nothing)
    endf_path = resolve(tapes, params.nendf)
    pendf_in_path = resolve(tapes, params.npendf_in)
    pendf_out_path = resolve(tapes, params.npendf_out)

    # Read input PENDF
    pendf_in = read_pendf(pendf_in_path)
    mf3_data = extract_mf3_all(pendf_in, params.mat)

    # Extract AWR and ZA from PENDF
    awr = _extract_awr_from_pendf(pendf_in, params.mat)
    za = _extract_za_from_pendf(pendf_in, params.mat)
    Z = extract_Z(za)

    @info "heatr: MAT=$(params.mat) Z=$Z AWR=$awr"

    # Get reconr result for _get_legacy_section (threshold-adjusted XS)
    if reconr_result === nothing
        reconr_result = reconr(endf_path; mat=params.mat, err=0.005)
    end

    # Get broadened energy grid — prefer raw data (full precision) over PENDF
    b_e = if broadened_data !== nothing
        broadened_data.energies
    else
        haskey(mf3_data, 1) || error("heatr: MT=1 not found on input PENDF")
        mf3_data[1][1]
    end

    # Read supplementary ENDF data
    gamma_data = _read_gamma_data(endf_path, params.mat)
    mu_bar = read_mf4_mubar(endf_path, params.mat)
    mf4_leg = read_mf4_legendre(endf_path, params.mat)
    mf13_gamma = _read_mf13_data(endf_path, params.mat)
    mf5_params = read_mf5_evaporation(endf_path, params.mat)

    # Read MF4 Legendre for all discrete inelastic MTs
    mf4_leg_all = Dict{Int, Tuple{Vector{Float64}, Vector{Vector{Float64}}}}()
    mf4_mubar_dict = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
    for mt_inel in 51:91
        leg = read_mf4_legendre(endf_path, params.mat; mt=mt_inel)
        if length(leg[1]) >= 2 && !(length(leg[1]) == 2 && leg[1][1] == 0.0)
            mf4_leg_all[mt_inel] = leg
            mf4_mubar_dict[mt_inel] = (leg[1], Float64[length(c) >= 2 ? c[2] : 0.0 for c in leg[2]])
        end
    end

    # Build PointwiseMaterial on broadened grid with ALL non-redundant MTs
    skip_mts = Set([1, 3, 4, 101, 203, 204, 207, 251, 252, 253])
    heatr_secs = [sec for sec in reconr_result.mf3_sections if !(Int(sec.mt) in skip_mts)]
    heatr_mts = Int32[sec.mt for sec in heatr_secs]
    heatr_qs = Dict{Int,Float64}()
    heatr_qm = Dict{Int,Float64}()
    heatr_lr = Dict{Int,Int32}()
    xs_cols = Matrix{Float64}(undef, length(b_e), length(heatr_secs))

    for (j, sec) in enumerate(heatr_secs)
        mt = Int(sec.mt)
        heatr_qs[mt] = Float64(sec.QI)
        heatr_qm[mt] = Float64(sec.QM)
        heatr_lr[mt] = sec.LR
        if broadened_data !== nothing && mt in broadened_data.partial_mts
            # Use raw broadened XS (full precision, no PENDF round-trip)
            col_idx = findfirst(==(mt), broadened_data.partial_mts)
            xs_cols[:, j] .= broadened_data.partials[:, col_idx]
        elseif mt == 2 && haskey(mf3_data, 2)
            _interp_to_grid!(view(xs_cols, :, j), b_e, mf3_data[2][1], mf3_data[2][2])
        elseif mt == 102 && haskey(mf3_data, 102)
            _interp_to_grid!(view(xs_cols, :, j), b_e, mf3_data[102][1], mf3_data[102][2])
        elseif mt == 18 && haskey(mf3_data, 18)
            _interp_to_grid!(view(xs_cols, :, j), b_e, mf3_data[18][1], mf3_data[18][2])
        else
            # Non-broadened MTs: use reconr-processed data with full Float64 precision
            # This matches the hand-wired pipeline which uses _get_legacy_section
            legacy = _get_legacy_section(reconr_result, mt)
            if legacy !== nothing
                _interp_to_grid!(view(xs_cols, :, j), b_e, legacy[1], legacy[2])
            else
                thresh = sec.tab.x[1]
                for i in eachindex(b_e)
                    xs_cols[i, j] = b_e[i] < thresh ? 0.0 : interpolate(sec.tab, b_e[i])
                end
            end
        end
    end

    pm_broad = PointwiseMaterial(Int32(params.mat), b_e, xs_cols, heatr_mts)

    # Build conbar params from MF5
    conbar_p = (mf5_params !== nothing && mf5_params.theta > 0) ?
               (mf5_params.u, mf5_params.theta) : nothing

    # Call compute_kerma
    kr = compute_kerma(pm_broad;
        awr=awr, Z=Z,
        gamma_data=gamma_data,
        mu_bar_data=mu_bar,
        mf4_legendre_data=mf4_leg,
        conbar_params=conbar_p,
        mf4_mubar_all=isempty(mf4_mubar_dict) ? nothing : mf4_mubar_dict,
        mf4_legendre_all=isempty(mf4_leg_all) ? nothing : mf4_leg_all,
        mf13_gamma=mf13_gamma,
        Q_values=heatr_qs, qm_values=heatr_qm, lr_values=heatr_lr)

    @info "heatr: KERMA computed, peak=$(round(maximum(abs, kr.total_kerma), sigdigits=4))"

    # Add MT=301 and MT=444 to output PENDF
    added_mf3 = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
    !isempty(kr.total_kerma) && (added_mf3[301] = (b_e, kr.total_kerma))
    !isempty(kr.damage_energy) && (added_mf3[444] = (b_e, kr.damage_energy))

    pendf_out = copy_with_modifications(pendf_in, params.mat; added_mf3=added_mf3)
    write_pendf_tape(pendf_out_path, pendf_out)

    @info "heatr: wrote $pendf_out_path"
    nothing
end

# =========================================================================
# Helpers
# =========================================================================

"""Extract ZA from the first MF3 HEAD record."""
function _extract_za_from_pendf(tape::PENDFTape, mat::Int)
    for material in tape.materials
        material.mat != mat && continue
        for sec in material.sections
            sec.mf != 3 && continue
            !isempty(sec.lines) || continue
            return parse_endf_float(rpad(sec.lines[1], 80)[1:11])
        end
    end
    return 0.0
end

"""Read MF12 gamma data as Dict{MT => Vector{(E_gamma, yield)}}."""
function _read_gamma_data(endf_path::String, mat::Int)
    gd = Dict{Int, Vector{Tuple{Float64,Float64}}}()
    gammas = read_mf12_gammas(endf_path, mat; mt=102)
    !isempty(gammas) && (gd[102] = gammas)
    gd
end

"""Read MF13 sections as Vector{(E_gamma, TabulatedFunction)}."""
function _read_mf13_data(endf_path::String, mat::Int)
    io = open(endf_path)
    try
        secs = read_mf13_sections(io, mat)
        return Tuple{Float64, TabulatedFunction}[(Float64(s.QM), s.tab) for s in secs]
    finally
        close(io)
    end
end

"""Interpolate (src_e, src_xs) onto target grid, filling zeros below threshold."""
function _interp_to_grid!(out::AbstractVector, target_e::Vector{Float64},
                          src_e::Vector{Float64}, src_xs::Vector{Float64})
    for i in eachindex(target_e)
        E = target_e[i]
        if E < src_e[1]
            out[i] = 0.0
        elseif E >= src_e[end]
            out[i] = src_xs[end]
        else
            idx = searchsortedfirst(src_e, E)
            if idx <= 1
                out[i] = src_xs[1]
            else
                f = (E - src_e[idx-1]) / (src_e[idx] - src_e[idx-1])
                out[i] = src_xs[idx-1] + f * (src_xs[idx] - src_xs[idx-1])
            end
        end
    end
end
