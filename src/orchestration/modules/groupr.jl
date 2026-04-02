# groupr module runner -- Group-averaged cross sections from PENDF
#
# Matches Fortran groupr.f90: read pointwise data from PENDF,
# compute flux-weighted group averages, write GENDF tape.

"""
    groupr_module(tapes::TapeManager, params::GrouprParams)

Run GROUPR: compute multigroup cross sections from PENDF data.
Write output GENDF tape. Matches Fortran groupr.f90 interface.

Input tapes: nendf (ENDF), npend (PENDF from reconr)
Output tapes: nout (GENDF multigroup output)
"""
function groupr_module(tapes::TapeManager, params::GrouprParams)
    @info "groupr: MAT=$(params.mat) ign=$(params.ign) → tape $(params.nout)"

    endf_path = resolve(tapes, params.nendf)
    pendf_path = resolve(tapes, params.npend)
    nout_path = resolve(tapes, params.nout)

    # Get group structure
    egn = _groupr_group_structure(params)
    ngn = length(egn) - 1
    @info "groupr: $ngn groups, $(length(params.mt_list)) MT requests"

    # Read ZA and AWR
    za, awr = _read_za_awr(endf_path, params.mat)

    # Temperature and sigz
    temp = isempty(params.temperatures) ? 0.0 : params.temperatures[1]
    sigz = isempty(params.sigz) ? 1.0e10 : params.sigz[1]

    # Weight function
    wfn = _groupr_weight_function(params.iwt)

    # Read PENDF data for fission XS (needed for nubar weighting)
    tape = read_pendf(pendf_path)
    mf3 = extract_mf3_all(tape, params.mat)

    # Process each MT request
    mt_results = Vector{NamedTuple{(:mfd,:mt,:name,:records),
                                   Tuple{Int,Int,String,Vector{NTuple{3,Float64}}}}}()

    for (mfd, mtd, name) in params.mt_list
        if mfd == 3
            if mtd == 452 || mtd == 455 || mtd == 456
                # Nubar: read from ENDF MF1, weight by fission XS from PENDF
                nubar_data = _read_nubar(endf_path, params.mat, mtd)
                if nubar_data !== nothing
                    # Get fission XS from PENDF (MT=18 or MT=19)
                    fis_mt = haskey(mf3, 18) ? 18 : (haskey(mf3, 19) ? 19 : 0)
                    records = _groupr_nubar_records(nubar_data, mf3, fis_mt,
                                                    egn, wfn)
                    push!(mt_results, (mfd=mfd, mt=mtd, name=name, records=records))
                end
            else
                # Standard XS from PENDF
                if haskey(mf3, mtd)
                    records = _groupr_xs_records(mf3[mtd], egn, wfn)
                    push!(mt_results, (mfd=mfd, mt=mtd, name=name, records=records))
                end
            end
        end
    end

    # Write GENDF output tape
    open(nout_path, "w") do io
        _write_groupr_tape(io, params.mat, za, awr, egn, mt_results,
                           params.title, temp, sigz)
    end
    register!(tapes, params.nout, nout_path)

    lines = countlines(nout_path)
    @info "groupr: wrote $nout_path ($lines lines, $ngn groups)"
    nothing
end

# =========================================================================
# Group structure
# =========================================================================

function _groupr_group_structure(params::GrouprParams)
    ign = params.ign
    if ign == 2 || ign == 3
        return collect(Float64, LANL_30)
    elseif ign == 4
        return collect(Float64, WIMS_69)
    else
        @warn "groupr: unsupported ign=$ign, using LANL-30"
        return collect(Float64, LANL_30)
    end
end

function _groupr_weight_function(iwt::Int)
    iwt == 2 && return constant_weight
    return inv_e_weight  # iwt=3 and default
end

# =========================================================================
# Nubar reading from ENDF MF1
# =========================================================================

"""Read nubar (MF1/MT452, 455, or 456) from ENDF file as (energies, values)."""
function _read_nubar(endf_path::String, mat::Int, mt::Int)
    open(endf_path) do io
        seekstart(io)
        find_section(io, 1, mt; target_mat=mat) || return nothing
        head = read_cont(io)
        lnu = Int(head.L2)  # LNU: 1=polynomial, 2=tabulated
        if lnu == 1
            # Polynomial representation: ν(E) = c₁ + c₂E + c₃E² + ...
            lst = read_list(io)
            nc = Int(lst.N1)
            coeffs = lst.data[1:nc]
            # Generate dense tabulation for proper integration
            # Use logarithmic spacing to cover full energy range
            energies = Float64[]
            for dec in -5:8
                for sub in [1.0, 2.0, 5.0]
                    e = sub * 10.0^dec
                    push!(energies, e)
                end
            end
            push!(energies, 2e7)
            sort!(unique!(energies))
            values = [sum(coeffs[k] * e^(k-1) for k in 1:nc) for e in energies]
            return (energies=energies, values=values)
        elseif lnu == 2
            # Tabulated
            tab = read_tab1(io)
            return (energies=collect(Float64, tab.x), values=collect(Float64, tab.y))
        end
        nothing
    end
end

# =========================================================================
# Group averaging using existing group_integrate infrastructure
# =========================================================================

"""Compute GENDF records for nubar: (flux, nubar_g, nu_sigf_integral) per group.
Nubar is weighted by fission XS: ν_g = ∫ν·σf·W dE / ∫σf·W dE."""
function _groupr_nubar_records(nubar_data, mf3::Dict, fis_mt::Int,
                                egn::Vector{Float64}, wfn)
    ngn = length(egn) - 1
    nu_e = nubar_data.energies
    nu_v = nubar_data.values

    # Build a merged energy grid from nubar + fission XS
    if fis_mt > 0 && haskey(mf3, fis_mt)
        fis_e, fis_xs = mf3[fis_mt]
    else
        # No fission XS — use flat weighting
        fis_e = nu_e
        fis_xs = ones(Float64, length(nu_e))
    end

    # Merge grids
    merged_e = sort!(unique!(vcat(nu_e, fis_e, egn)))

    # Interpolate nubar and fission XS onto merged grid
    nu_merged = _interp_linlin(merged_e, nu_e, nu_v)
    fis_merged = _interp_linlin(merged_e, fis_e, fis_xs)

    # Build product arrays: W(E), σf(E)·W(E), ν(E)·σf(E)·W(E)
    n = length(merged_e)
    w_vals = [wfn(merged_e[i]) for i in 1:n]
    sigf_w = [fis_merged[i] * w_vals[i] for i in 1:n]
    nu_sigf_w = [nu_merged[i] * fis_merged[i] * w_vals[i] for i in 1:n]

    # Integrate over each group using exact panel integration
    flux = group_integrate(merged_e, w_vals, egn)
    sigf_int = group_integrate(merged_e, sigf_w, egn)
    nu_sigf_int = group_integrate(merged_e, nu_sigf_w, egn)

    # Build records: (flux, nubar_g, sigf_avg)
    # Fortran displa transforms ratio records to (flux, yield, cross_section)
    records = NTuple{3,Float64}[]
    for g in 1:ngn
        nubar_g = sigf_int[g] > 0 ? nu_sigf_int[g] / sigf_int[g] : 0.0
        sigf_avg = flux[g] > 0 ? sigf_int[g] / flux[g] : 0.0
        # Round to 7 sigfigs matching Fortran displa (sigfig(x,7,0))
        nubar_g = round_sigfig(nubar_g, 7)
        sigf_avg = round_sigfig(sigf_avg, 7)
        push!(records, (round_sigfig(flux[g], 7), nubar_g, sigf_avg))
    end
    records
end

"""Compute GENDF records for standard XS: (flux, sigma_g, sigma_flux) per group."""
function _groupr_xs_records(xs_data::Tuple, egn::Vector{Float64}, wfn)
    energies, xs = xs_data
    ngn = length(egn) - 1

    w_vals = [wfn(energies[i]) for i in 1:length(energies)]
    sig_w = [xs[i] * w_vals[i] for i in 1:length(energies)]

    flux = group_integrate(energies, w_vals, egn)
    sig_int = group_integrate(energies, sig_w, egn)

    records = NTuple{3,Float64}[]
    for g in 1:ngn
        sigma_g = flux[g] > 0 ? sig_int[g] / flux[g] : 0.0
        push!(records, (flux[g], sigma_g, sig_int[g]))
    end
    records
end

"""Linear interpolation of (x_data, y_data) at points x_out."""
function _interp_linlin(x_out::Vector{Float64}, x_data::Vector{Float64},
                        y_data::Vector{Float64})
    n = length(x_out)
    y_out = Vector{Float64}(undef, n)
    j = 1
    for i in 1:n
        x = x_out[i]
        if x <= x_data[1]
            y_out[i] = y_data[1]
        elseif x >= x_data[end]
            y_out[i] = y_data[end]
        else
            while j < length(x_data) - 1 && x_data[j+1] < x
                j += 1
            end
            frac = (x - x_data[j]) / (x_data[j+1] - x_data[j])
            y_out[i] = y_data[j] + frac * (y_data[j+1] - y_data[j])
        end
    end
    y_out
end

# =========================================================================
# GENDF output tape writer
# =========================================================================

"""Write GROUPR output tape in GENDF format matching Fortran groupr."""
function _write_groupr_tape(io::IO, mat::Int, za::Float64, awr::Float64,
                            egn::Vector{Float64}, mt_results, title::String,
                            temp::Float64, sigz::Float64)
    ngn = length(egn) - 1

    # TPID record
    @printf(io, "%-66s%4d%2d%3d%5d\n", title, 0, 0, 0, 0)

    seq = 1
    # MF1/MT451 HEAD: ZA, AWR, 0, ntemp=1, NFC=-1, nmod=1
    _write_cont_line(io, za, awr, 0, 1, -1, 1, mat, 1, 451, seq); seq += 1
    # CONT: 0, 0, NGN, 0, NW, 0
    nw_dir = ngn + 4  # 2 (elow,ehigh) + ngn+1 (bounds) + 1 (trailing zero)
    _write_cont_line(io, 0.0, 0.0, ngn, 0, nw_dir, 0, mat, 1, 451, seq); seq += 1

    # Write elow=0, ehigh=max, then group boundaries, then trailing zero
    all_vals = Float64[0.0, 1e10]
    append!(all_vals, egn)
    push!(all_vals, 0.0)

    idx = 1
    while idx <= length(all_vals)
        buf = ""
        for col in 1:6
            idx > length(all_vals) && break
            buf *= format_endf_float(all_vals[idx]); idx += 1
        end
        _write_data_line(io, buf, mat, 1, 451, seq); seq += 1
    end

    # FEND for MF1
    _write_fend_line(io, mat)

    # MF3 sections for each MT
    for res in mt_results
        mt = res.mt
        records = res.records
        seq = 1

        # HEAD: ZA, 0, NL=1, NG2=1, NZ=0, NGN
        _write_cont_line(io, za, 0.0, 1, 1, 0, ngn, mat, 3, mt, seq); seq += 1

        # One record per group: CONT + data
        for (g, rec) in enumerate(records)
            # CONT: 0, 0, NW=3, NG2=1, NW=3, IG
            _write_cont_line(io, 0.0, 0.0, 3, 1, 3, g, mat, 3, mt, seq); seq += 1
            # Data: v1, v2, v3
            buf = format_endf_float(rec[1]) *
                  format_endf_float(rec[2]) *
                  format_endf_float(rec[3])
            _write_data_line(io, buf, mat, 3, mt, seq); seq += 1
        end

        _write_send_line(io, mat, 3)
    end
    _write_fend_line(io, mat)

    # MEND
    _write_fend_line(io, 0)
end
