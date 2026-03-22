# module_inspector.jl — Deep per-module diagnostic logging
#
# Wraps existing execute_* functions with diagnostic capture: grid stats,
# probe cross sections at 26 standard energies, module-specific extras
# (resonance bounds, broadening parameters, etc.)

using NJOY
using Printf

include("test_executor.jl")

# =========================================================================
# Types
# =========================================================================

"""Diagnostic snapshot from one module execution."""
struct ModuleDiagnostics
    module_name::Symbol
    elapsed_ms::Float64
    n_grid_in::Int
    n_grid_out::Int
    e_range::Tuple{Float64, Float64}
    probe_xs::Dict{Float64, Dict{Int, Float64}}  # energy -> MT -> XS value
    extra::Dict{Symbol, Any}                       # module-specific diagnostics
end

# Standard probe energies covering key physics regions
const STANDARD_PROBES = Float64[
    1e-5, 1e-4, 1e-3, 0.01, 0.0253, 0.05, 0.1, 0.5,
    1.0, 5.0, 10.0, 50.0, 100.0, 500.0,
    1e3, 5e3, 1e4, 5e4, 1e5, 5e5,
    1e6, 2e6, 5e6, 1e7, 1.5e7, 2e7,
]

# =========================================================================
# Probe XS sampling
# =========================================================================

"""Sample cross sections at specified energies for all available MTs."""
function probe_xs_at(result, probe_energies::Vector{Float64}=STANDARD_PROBES)
    probes = Dict{Float64, Dict{Int, Float64}}()

    if result === nothing
        return probes
    end

    # Get energies and MT data
    e, mt_data = _extract_energy_mt_data(result)
    isempty(e) && return probes

    for ep in probe_energies
        (ep < e[1] || ep > e[end]) && continue
        mt_vals = Dict{Int, Float64}()
        for (mt, xs) in mt_data
            mt_vals[mt] = _interp_at_mi(e, xs, ep)
        end
        probes[ep] = mt_vals
    end
    probes
end

"""Extract (energies, Dict(MT => xs_array)) from various result types."""
function _extract_energy_mt_data(result)
    if result isa PointwiseMaterial
        e = result.energies
        mt_data = Dict{Int, Vector{Float64}}()
        for (col, mt) in enumerate(result.mt_list)
            mt_data[mt] = result.cross_sections[:, col]
        end
        return e, mt_data
    elseif hasproperty(result, :energies)
        e = result.energies
        mt_data = Dict{Int, Vector{Float64}}()
        for (mt, sym) in [(1, :total), (2, :elastic), (18, :fission), (102, :capture)]
            hasproperty(result, sym) && (mt_data[mt] = getproperty(result, sym))
        end
        return e, mt_data
    end
    (Float64[], Dict{Int, Vector{Float64}}())
end

# =========================================================================
# Sum rule checking
# =========================================================================

"""
    check_sum_rules(result) -> violations::Vector{Tuple{Float64, Float64}}

Check that total >= elastic + fission + capture at every energy.
Returns (energy, violation_magnitude) pairs where the rule fails.
"""
function check_sum_rules(result)
    violations = Tuple{Float64, Float64}[]
    e, mt_data = _extract_energy_mt_data(result)
    isempty(e) && return violations
    haskey(mt_data, 1) || return violations

    total = mt_data[1]
    partials = zeros(length(e))
    for mt in [2, 18, 102]
        haskey(mt_data, mt) && (partials .+= mt_data[mt])
    end

    for i in eachindex(e)
        deficit = partials[i] - total[i]
        if deficit > max(1e-6, 1e-4 * abs(total[i]))
            push!(violations, (e[i], deficit))
        end
    end
    violations
end

# =========================================================================
# Per-module inspectors
# =========================================================================

"""Run reconr with detailed diagnostics."""
function inspect_reconr(endf_file::String, mc::ModuleCall, mat::Int,
                        err::Float64)
    t0 = time()
    extra = Dict{Symbol, Any}()

    # Run reconr
    r, step = execute_reconr(endf_file, mc; mat_override=mat, err_override=err)
    elapsed_ms = (time() - t0) * 1000

    if r === nothing
        diag = ModuleDiagnostics(:reconr, elapsed_ms, 0, 0, (0.0, 0.0),
                                  Dict{Float64, Dict{Int, Float64}}(), extra)
        return diag, nothing, step
    end

    n_out = length(r.energies)
    e_range = (r.energies[1], r.energies[end])

    # Extract extra diagnostics from the reconr result
    # The reconr() function returns mf2 and mf3_sections in the NamedTuple
    try
        if hasproperty(r, :mf2)
            mf2 = r.mf2
            # Resonance bounds via internal function
            bounds = NJOY._resonance_bounds(mf2)
            extra[:eresl] = bounds[1]
            extra[:eresh] = bounds[2]
            extra[:eresr] = length(bounds) >= 3 ? bounds[3] : bounds[2]

            # Formalism detection
            if !isempty(mf2.isotopes) && !isempty(mf2.isotopes[1].ranges)
                rng = mf2.isotopes[1].ranges[1]
                extra[:formalism] = _formalism_name(Int(rng.LRF))
            end

            # Count resonances (parameters field holds formalism-specific data)
            n_res = 0
            for iso in mf2.isotopes
                for rng in iso.ranges
                    p = rng.parameters
                    if hasproperty(p, :Er)  # SLBW/MLBW
                        n_res += sum(length, p.Er; init=0)
                    elseif hasproperty(p, :resonances)
                        n_res += length(p.resonances)
                    end
                end
            end
            extra[:n_resonances] = n_res
        end

        if hasproperty(r, :mf3_sections)
            extra[:n_mf3_sections] = length(r.mf3_sections)
        end
    catch ex
        extra[:parse_error] = sprint(showerror, ex)
    end

    # Sum rule check
    violations = check_sum_rules(r)
    extra[:sum_rule_violations] = length(violations)
    if !isempty(violations)
        extra[:sum_rule_worst] = violations[1]
    end

    probes = probe_xs_at(r)

    ModuleDiagnostics(:reconr, elapsed_ms, 0, n_out, e_range, probes, extra),
    r, step
end

"""Run broadr with detailed diagnostics."""
function inspect_broadr(pendf::PointwiseMaterial, mc::ModuleCall; awr::Float64=1.0)
    t0 = time()
    extra = Dict{Symbol, Any}()
    n_in = length(pendf.energies)

    p = parse_broadr(mc)
    temp = isempty(p.temperatures) ? 0.0 : p.temperatures[1]
    extra[:temperature] = temp
    extra[:thnmax] = p.thnmax

    if temp > 0 && awr > 0
        k_boltz = 8.617333e-5  # eV/K
        extra[:alpha_param] = awr / (k_boltz * temp)
    end

    r, step = execute_broadr(pendf, mc; awr=awr)
    elapsed_ms = (time() - t0) * 1000

    n_out = r !== nothing ? length(r.energies) : 0
    extra[:n_thinned] = n_in - n_out

    e_range = r !== nothing ? (r.energies[1], r.energies[end]) : (0.0, 0.0)
    probes = probe_xs_at(r)

    ModuleDiagnostics(:broadr, elapsed_ms, n_in, n_out, e_range, probes, extra),
    r, step
end

"""Run heatr with detailed diagnostics."""
function inspect_heatr(pendf::PointwiseMaterial, mc::ModuleCall; awr::Float64=1.0)
    t0 = time()
    extra = Dict{Symbol, Any}()
    n_in = length(pendf.energies)

    r, step = execute_heatr(pendf, mc; awr=awr)
    elapsed_ms = (time() - t0) * 1000

    n_out = r !== nothing ? length(r.energies) : n_in
    e_range = r !== nothing ? (r.energies[1], r.energies[end]) : (pendf.energies[1], pendf.energies[end])
    probes = probe_xs_at(r !== nothing ? r : pendf)

    ModuleDiagnostics(:heatr, elapsed_ms, n_in, n_out, e_range, probes, extra),
    r !== nothing ? r : pendf, step
end

"""Run thermr with detailed diagnostics."""
function inspect_thermr(pendf::PointwiseMaterial, mc::ModuleCall)
    t0 = time()
    extra = Dict{Symbol, Any}()
    n_in = length(pendf.energies)

    p = parse_thermr(mc)
    extra[:temperature] = isempty(p.temperatures) ? 296.0 : p.temperatures[1]
    extra[:emax] = p.emax
    extra[:model] = p.nin_thermal > 0 ? :sab : :free_gas
    extra[:mat_thermal] = p.mat_thermal

    r, step = execute_thermr(pendf, mc)
    elapsed_ms = (time() - t0) * 1000

    n_out = r !== nothing ? length(r.energies) : n_in
    e_range = r !== nothing ? (r.energies[1], r.energies[end]) : (pendf.energies[1], pendf.energies[end])
    probes = probe_xs_at(r !== nothing ? r : pendf)

    ModuleDiagnostics(:thermr, elapsed_ms, n_in, n_out, e_range, probes, extra),
    r !== nothing ? r : pendf, step
end

# =========================================================================
# Formatting
# =========================================================================

"""Pretty-print a ModuleDiagnostics record."""
function format_diagnostics(diag::ModuleDiagnostics; indent::String="  ")
    buf = IOBuffer()
    @printf(buf, "%s%-8s  %d → %d pts  [%.3e, %.3e] eV  (%.0f ms)\n",
            indent, diag.module_name, diag.n_grid_in, diag.n_grid_out,
            diag.e_range[1], diag.e_range[2], diag.elapsed_ms)

    # Module-specific extras
    for (k, v) in sort(collect(diag.extra); by=first)
        @printf(buf, "%s  %s = %s\n", indent, k, _format_val(v))
    end

    # Probe XS table (compact: show a few key energies)
    key_probes = [0.0253, 1.0, 1e3, 1e5, 1e7]
    avail = sort(collect(keys(diag.probe_xs)))
    if !isempty(avail)
        @printf(buf, "%s  Probe XS:\n", indent)
        @printf(buf, "%s    %12s  %12s  %12s  %12s  %12s\n",
                indent, "E (eV)", "MT1", "MT2", "MT18", "MT102")
        for ep in key_probes
            # Find closest available probe
            _, idx = findmin(abs.(avail .- ep))
            actual_e = avail[idx]
            abs(actual_e - ep) / max(ep, 1e-30) > 0.5 && continue

            mt_vals = diag.probe_xs[actual_e]
            @printf(buf, "%s    %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n",
                    indent, actual_e,
                    get(mt_vals, 1, NaN), get(mt_vals, 2, NaN),
                    get(mt_vals, 18, NaN), get(mt_vals, 102, NaN))
        end
    end

    String(take!(buf))
end

_format_val(v::Float64) = @sprintf("%.6g", v)
_format_val(v::Int) = string(v)
_format_val(v::Symbol) = string(v)
_format_val(v::Tuple) = join(map(_format_val, v), ", ")
_format_val(v) = string(v)

function _formalism_name(lrf::Int)
    names = Dict(0 => :none, 1 => :slbw, 2 => :mlbw, 3 => :reich_moore, 7 => :sammy)
    get(names, lrf, :unknown)
end

# =========================================================================
# Interpolation helper
# =========================================================================

function _interp_at_mi(e::Vector{Float64}, xs::Vector{Float64}, ep::Float64)
    n = length(e)
    n == 0 && return 0.0
    ep <= e[1] && return xs[1]
    ep >= e[n] && return xs[n]
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >> 1
        e[mid] <= ep ? (lo = mid) : (hi = mid)
    end
    t = (ep - e[lo]) / max(e[hi] - e[lo], 1e-30)
    xs[lo] + t * (xs[hi] - xs[lo])
end
