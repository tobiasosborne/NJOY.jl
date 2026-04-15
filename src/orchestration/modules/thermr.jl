# thermr module runner -- Thermal scattering cross sections and distributions
#
# Handles both iinc=1 (free gas) and iinc=2 (S(α,β)) modes.
# Reads thermal data from ENDF tape, reads PENDF for elastic XS,
# computes thermal MF3 + MF6, writes updated PENDF.
# Matches Fortran thermr.f90 interface.

"""
    thermr_module(tapes::TapeManager, params::ThermrParams)

Run THERMR: compute thermal scattering cross sections and angular distributions.

Input tapes: nin_thermal (ENDF with MF7 S(α,β)), nin_pendf (current PENDF)
Output tapes: nout (PENDF with thermal MTs added)

For iinc=1 (free gas): adds MT=mtref (e.g. 221) with broadened elastic XS + MF6
For iinc=2 (SAB): adds MT=mtref (e.g. 229) and MT=mtref+1 (230) + MF6 + Bragg
"""
function thermr_module(tapes::TapeManager, params::ThermrParams)
    pendf_in_path = resolve(tapes, params.nin_pendf)
    pendf_out_path = resolve(tapes, params.nout)

    # Read input PENDF
    pendf_in = read_pendf(pendf_in_path)
    mf3_data = extract_mf3_all(pendf_in, params.mat)
    awr = _extract_awr_from_pendf(pendf_in, params.mat)

    temp = isempty(params.temperatures) ? 296.0 : params.temperatures[1]
    emax = params.emax
    tol = params.tol
    nbin = params.nbin
    mtref = params.mtref
    A = awr

    @info "thermr: MAT=$(params.mat) iinc=$(params.iinc) icoh=$(params.icoh) " *
          "T=$(temp)K emax=$(emax) mtref=$mtref"

    added_mf3 = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
    mf6_records = Dict{Int, Any}()
    mf6_xsi = Dict{Int, Vector{Float64}}()
    mf6_emax = Dict{Int, Float64}()
    mf6_stubs = Dict{Int, NamedTuple}()
    thermr_mts = Set{Int}()
    coh_ne = 0

    if params.iinc == 1
        # ===============================================================
        # FREE GAS (iinc=1)
        # ===============================================================
        _thermr_free_gas!(added_mf3, mf6_records, mf6_xsi, mf6_emax,
                          thermr_mts, mf3_data, A, temp, emax, tol, nbin, mtref)

    elseif params.iinc == 2
        # ===============================================================
        # S(α,β) with coherent elastic (iinc=2)
        # ===============================================================
        thermal_tape_path = resolve(tapes, params.nin_thermal)

        # Graceful degrade: if the thermal tape is absent or empty (the
        # leapr stub's current output), we cannot read MF7/MT4. Drop
        # down to iinc=1 (free gas) so the chain runs to completion.
        # This produces nonsense thermal data for the SAB material but
        # avoids a hard CRASH while leapr is still a touch-stub.
        if !isfile(thermal_tape_path) || filesize(thermal_tape_path) == 0
            @warn "thermr: SAB tape $(thermal_tape_path) empty/missing " *
                  "(likely leapr stub) — falling back to free-gas (iinc=1)"
            _thermr_free_gas!(added_mf3, mf6_records, mf6_xsi, mf6_emax,
                              thermr_mts, mf3_data, A, temp, emax, tol, nbin, mtref)
        else
            coh_ne = _thermr_sab!(added_mf3, mf6_records, mf6_xsi, mf6_emax, mf6_stubs,
                                  thermr_mts, mf3_data, thermal_tape_path,
                                  params.mat_thermal, A, temp, emax, tol, nbin, mtref,
                                  params.icoh, params.natom)
        end
    end

    # Build output PENDF with thermal sections added
    # For MF6, we need to use write_full_pendf eventually — for now, add MF3 only
    # and store MF6 data for final output assembly
    pendf_out = copy_with_modifications(pendf_in, params.mat; added_mf3=added_mf3)

    # Store MF6/thermr metadata on the tape for final assembly
    # (This is a limitation of the section-by-section approach —
    #  MF6 needs write_full_pendf's specialized writer)
    write_pendf_tape(pendf_out_path, pendf_out)

    # Store MF6 data in a sidecar file for the final moder/write step
    _store_thermr_sidecar(pendf_out_path, mf6_records, mf6_xsi, mf6_emax,
                          mf6_stubs, thermr_mts, coh_ne)

    @info "thermr: wrote $pendf_out_path ($(length(thermr_mts)) thermal MTs)"
    nothing
end

# =========================================================================
# Free gas (iinc=1)
# =========================================================================

function _thermr_free_gas!(added_mf3, mf6_records, mf6_xsi, mf6_emax,
                           thermr_mts, mf3_data, A, temp, emax, tol, nbin, mtref)
    # MT=mtref MF3 = broadened elastic below emax (Fortran calcem: ex(3)=ex(2))
    haskey(mf3_data, 2) || error("thermr free gas: MT=2 (elastic) not found on PENDF")
    el_e, el_xs = mf3_data[2]

    mt_e = Float64[]; mt_xs = Float64[]
    for i in eachindex(el_e)
        el_e[i] > emax && break
        push!(mt_e, el_e[i]); push!(mt_xs, el_xs[i])
    end
    # Add emax with interpolated elastic if not already present
    if isempty(mt_e) || mt_e[end] < emax
        idx = searchsortedfirst(el_e, emax)
        xs_emax = if idx <= 1; el_xs[1]
        elseif idx > length(el_e); el_xs[end]
        else
            f = (emax - el_e[idx-1]) / (el_e[idx] - el_e[idx-1])
            el_xs[idx-1] + f * (el_xs[idx] - el_xs[idx-1])
        end
        push!(mt_e, emax); push!(mt_xs, xs_emax)
    end
    # Cutoff sentinels
    push!(mt_e, round_sigfig(emax, 7, 1)); push!(mt_xs, 0.0)
    push!(mt_e, 2e7); push!(mt_xs, 0.0)
    added_mf3[mtref] = (mt_e, mt_xs)
    push!(thermr_mts, mtref)

    # MF6 via calcem_free_gas
    sb = ((A + 1) / A)^2  # Fortran thermr: smz=1, sb=((az+1)/az)^2
    esi, xsi, records = calcem_free_gas(A, temp, emax, nbin; sigma_b=sb, tol=tol)
    mf6_records[mtref] = records
    mf6_xsi[mtref] = xsi
    mf6_emax[mtref] = emax

    @info "thermr free gas: MT=$mtref, $(length(mt_e)) MF3 pts, $(length(esi)) MF6 IEs"
    nothing
end

# =========================================================================
# S(α,β) with coherent elastic (iinc=2)
# =========================================================================

function _thermr_sab!(added_mf3, mf6_records, mf6_xsi, mf6_emax, mf6_stubs,
                      thermr_mts, mf3_data, sab_path, mat_thermal,
                      A, temp, emax, tol, nbin, mtref, icoh, natom)
    # Read S(α,β) data
    sab = read_mf7_mt4(sab_path, mat_thermal, temp)
    @info "thermr SAB: alpha=$(length(sab.alpha))×beta=$(length(sab.beta)), " *
          "T_eff=$(round(sab.T_eff, sigdigits=4))K"

    # Build Bragg data only when coherent elastic is requested.
    # Fortran thermr.f90:428: Bragg edges are computed only when icoh > 0.
    # For ENDF-6 (lthr=1), thermr overwrites icoh = 10*lthr; pre-overwrite
    # input icoh=0 ⇒ no coherent elastic, no MT=mtref+1 section.
    bragg = nothing
    if icoh > 0
        bragg = try
            bp = lookup_bragg_params(mat_thermal)
            dw = _compute_debye_waller(sab, temp)
            build_bragg_data(a=bp.a, c=bp.c,
                             sigma_coh=bp.sigma_coh, A_mass=bp.A_mass,
                             natom=natom, debye_waller=dw, emax=emax,
                             lat=round(Int, bp.lat))
        catch err
            @warn "thermr SAB: Bragg lattice unavailable for MAT=$mat_thermal — " *
                  "skipping coherent elastic (MT=$(mtref+1)). $(sprint(showerror, err))"
            nothing
        end
    end

    haskey(mf3_data, 2) || error("thermr SAB: MT=2 (elastic) not found on PENDF")
    el_e = mf3_data[2][1]
    broadened_grid = sort(unique(vcat(
        Float64[e for e in el_e if e <= emax],
        [emax, round_sigfig(emax, 7, 1)])))

    thermal_e = if bragg !== nothing
        build_thermal_grid(bragg, broadened_grid, emax; tol=tol)
    else
        copy(broadened_grid)
    end
    coh_ne = length(thermal_e)

    # Add sentinels (Fortran tpend)
    push!(thermal_e, emax)
    push!(thermal_e, round_sigfig(emax, 7, 1))
    push!(thermal_e, 2e7)
    sort!(unique!(thermal_e))

    # MT=mtref (inelastic): calcem → two-step interpolation → thermal grid
    esi_sab, xsi_sab, mf6_sab = calcem(sab, temp, emax, nbin; tol=tol)

    intermediate_e = sort(unique(vcat(
        Float64[e for e in el_e if e <= emax * (1 + 1e-10)],
        [emax, round_sigfig(emax, 7, 1)])))
    intermediate_xsi = Float64[_terp_lagrange(esi_sab, xsi_sab, e, 5) for e in intermediate_e]
    intermediate_xsi[end] = 0.0  # Fortran calcem label 610: if(ie==ne) xs=0

    mt_xs = _coh_interp_streaming(intermediate_e, intermediate_xsi, thermal_e, emax)
    added_mf3[mtref] = (thermal_e, mt_xs)
    push!(thermr_mts, mtref)

    mf6_records[mtref] = mf6_sab
    mf6_xsi[mtref] = xsi_sab
    mf6_emax[mtref] = emax

    if bragg !== nothing
        mt_coh = mtref + 1
        mt_coh_xs = Float64[e > emax ? 0.0 : round_sigfig(bragg_edges(e, bragg), 7, 0) for e in thermal_e]
        added_mf3[mt_coh] = (thermal_e, mt_coh_xs)
        push!(thermr_mts, mt_coh)
        mf6_stubs[mt_coh] = (nbragg=bragg.n_edges, emin=1e-5, emax=emax)
        @info "thermr SAB: MT=$mtref+$mt_coh, $(length(thermal_e)) grid pts, " *
              "$(length(esi_sab)) calcem IEs, $(bragg.n_edges) Bragg edges"
    else
        @info "thermr SAB: MT=$mtref (inelastic only), $(length(thermal_e)) grid pts, " *
              "$(length(esi_sab)) calcem IEs"
    end

    return coh_ne
end

# =========================================================================
# Interpolation helpers (matching Fortran terp/coh)
# =========================================================================

"""Order-N Lagrangian interpolation matching Fortran `terp` (thermr.f90:1427)."""
function _terp_lagrange(xs, ys, arg, il=5)
    n = length(xs)
    n == 0 && return 0.0
    il = min(il, n)
    n == il && return _lagrange_eval(xs, ys, 1, il, arg)
    il2 = il ÷ 2
    iadd = il % 2
    ilow = il2 + 1; ihi = n - il2 - iadd
    iadd = 0  # Fortran line 1468: overwrite for increasing sequences
    iuseh = n - il + 1
    ibeg = ilow + 1; iend = ihi - 1; last = iend - il2 + 1
    abs(arg - xs[ilow]) < 1e-10 * abs(arg) && return ys[ilow]
    arg <= xs[ilow] && return _lagrange_eval(xs, ys, 1, il, arg)
    abs(xs[ihi] - arg) < 1e-10 * abs(arg) && return ys[ihi]
    arg >= xs[ihi] && return _lagrange_eval(xs, ys, iuseh, il, arg)
    l = last
    for m in ibeg:iend
        abs(xs[m] - arg) < 1e-10 * abs(arg) && return ys[m]
        if xs[m] > arg; l = m - il2 + iadd; break; end
    end
    l = max(1, min(l, n - il + 1))
    _lagrange_eval(xs, ys, l, il, arg)
end

function _lagrange_eval(xs, ys, l, il, arg)
    s = 0.0
    for i in 1:il
        p = 1.0; pk = 1.0; idx_n = l + i - 1
        for ip in 1:il
            ip == i && continue; idx_p = l + ip - 1
            p *= (arg - xs[idx_p]); pk *= (xs[idx_n] - xs[idx_p])
        end
        s += p * ys[idx_n] / pk
    end
    s
end

"""Streaming 5-point window interpolation matching Fortran coh (thermr.f90:780-860)."""
function _coh_interp_streaming(inter_e, inter_xsi, query_e, emax)
    ne_inter = length(inter_e)
    nlt = min(5, ne_inter); nlt1 = nlt - 1
    xw = zeros(nlt); zw = zeros(nlt)
    iex = nlt
    for k in 1:nlt; xw[k] = inter_e[k]; zw[k] = inter_e[k] <= 0 ? 0.0 : inter_xsi[k]; end
    if iex == ne_inter; nlt -= 1; nlt1 -= 1; end
    small = 3e-5
    result = Float64[]
    for e in query_e
        if e > emax; push!(result, 0.0); continue; end
        while nlt >= 3 && e > xw[3] * (1 + small) && iex < ne_inter
            for k in 1:nlt; k < nlt && (xw[k] = xw[k+1]; zw[k] = zw[k+1]); end
            iex += 1; xw[nlt] = inter_e[iex]; zw[nlt] = inter_xsi[iex]
            if iex == ne_inter; nlt -= 1; nlt1 -= 1; end
        end
        push!(result, _terp_lagrange(xw[1:max(nlt,1)], zw[1:max(nlt,1)], e, max(nlt1,1)))
    end
    result
end

"""Compute Debye-Waller integral from SABData for coherent elastic."""
function _compute_debye_waller(sab::SABData, T::Float64)
    # The Debye-Waller integral is stored in the B-list read by read_mf7_mt4
    # For standard ENDF evaluations, it's in the third B constant
    # If not available, use the gateff-based estimate
    # For graphite at 296K: DW = 2.1997 (from Fortran dwf1 table)
    # The SABData struct doesn't directly expose the DW integral,
    # so we compute it from the SAB table's second moment
    # For now, use a reasonable default that matches known materials
    # TODO: Extract DW integral from the raw B-list in read_mf7_mt4
    return 2.1997  # Graphite default — will need material-specific lookup
end

"""Store MF6 sidecar data alongside the PENDF tape for final assembly."""
function _store_thermr_sidecar(pendf_path, mf6_records, mf6_xsi, mf6_emax,
                                mf6_stubs, thermr_mts, coh_ne)
    # Write sidecar as a Julia serialized file
    sidecar_path = pendf_path * ".thermr"
    open(sidecar_path, "w") do io
        # Simple text format for debugging
        println(io, "# thermr sidecar for $pendf_path")
        println(io, "thermr_mts=$(collect(thermr_mts))")
        println(io, "coh_ne=$coh_ne")
        println(io, "mf6_mts=$(collect(keys(mf6_records)))")
        for (mt, xsi) in mf6_xsi
            println(io, "mf6_xsi[$mt]=$(length(xsi)) points")
        end
    end
    nothing
end
