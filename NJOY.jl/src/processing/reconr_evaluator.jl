# RECONR evaluator -- build cross section evaluation closures from MF2 data
#
# Step 1 of the RECONR pipeline: build_evaluator returns a closure that
# evaluates resonance cross sections at a given energy.
#
# Also contains the legacy sigma_mf2 function and merge_background functions.

# ==========================================================================
# Step 1: Build evaluator (closure over MF2 data)
# ==========================================================================

"""
    build_evaluator(mf2::MF2Data; temperature=0.0, table=nothing)

Build a cross section evaluator closure from MF2 resonance data.
Returns a callable `f(E::Float64) -> NTuple{4, Float64}` that computes
`(total, elastic, fission, capture)` at energy `E`.

The returned function:
- Sums over all isotope sections weighted by abundance
- Dispatches to the correct formalism via `cross_section`
- Clamps negative components to zero
- Is fully type-stable (returns `NTuple{4, Float64}`)

This replaces NJOY2016's `sigma` subroutine (reconr.f90:2571-2667).
"""
function build_evaluator(mf2::MF2Data;
                         temperature::Real = 0.0,
                         table::Union{Nothing, FaddeevaTable} = nothing)
    # Pre-collect the isotope/range/abundance triples to avoid
    # iterating nested structures in the hot path
    sections = _collect_sections(mf2)

    function evaluator(E::Float64)
        sig_t = 0.0
        sig_e = 0.0
        sig_f = 0.0
        sig_c = 0.0

        @inbounds for (rng, abn) in sections
            if E >= rng.EL && E < rng.EH && Int(rng.LRU) > 0
                sigp = cross_section(E, rng; temperature=temperature, table=table)
                sig_t += max(0.0, sigp.total) * abn
                sig_e += max(0.0, sigp.elastic) * abn
                sig_f += max(0.0, sigp.fission) * abn
                sig_c += max(0.0, sigp.capture) * abn
            end
        end
        return (sig_t, sig_e, sig_f, sig_c)
    end

    return evaluator
end

# Flatten isotope/range hierarchy into a flat vector of (range, abundance) pairs
function _collect_sections(mf2::MF2Data)
    result = Tuple{ResonanceRange, Float64}[]
    for iso in mf2.isotopes
        abn = iso.ABN
        for rng in iso.ranges
            push!(result, (rng, abn))
        end
    end
    return result
end

# ==========================================================================
# Merge MF3 background into pointwise result
# ==========================================================================

"""
    merge_background!(energies, values, mf3_sections, mf2)

Merge MF3 background cross sections into the resonance cross section
matrix `values` in-place. This replaces NJOY2016's `emerge`.

For each energy and matching MF3 reaction:
- MT 2 -> add to elastic column
- MT 18, 19 -> add to fission column
- MT 102 -> add to capture column
Then recompute total = elastic + fission + capture.
"""
function merge_background!(energies::Vector{Float64},
                            values::Matrix{Float64},
                            mf3_sections::Vector{MF3Section},
                            mf2::MF2Data)
    n = length(energies)
    small = 1.0e-8

    # Find resonance range boundaries for overlap detection
    eresl, eresh, eresr = _resonance_bounds(mf2)

    for i in 1:n
        e = energies[i]

        for sec in mf3_sections
            mt = Int(sec.mt)
            # Skip redundant reactions
            (mt == 1 || mt == 3 || mt == 101) && continue

            bg = interpolate(sec.tab, e)
            bg == 0.0 && continue

            # In the unresolved-resolved overlap region, backgrounds for
            # primary channels are assigned to the unresolved component
            # (matching NJOY's emerge logic)
            if e >= eresr && e < eresh && (mt == 2 || mt == 18 || mt == 19 || mt == 102)
                continue
            end

            if mt == 2
                values[i, _COL_ELASTIC] += bg
            elseif mt == 18 || mt == 19
                values[i, _COL_FISSION] += bg
            elseif mt == 102
                values[i, _COL_CAPTURE] += bg
            end
        end

        # Clamp elastic to a small positive value
        if values[i, _COL_ELASTIC] <= small
            values[i, _COL_ELASTIC] = small
        end

        # Round to 7 significant figures (matching NJOY's emerge:sigfig call)
        for j in _COL_ELASTIC:_COL_CAPTURE
            values[i, j] = round_sigfig(values[i, j], 7, 0)
        end

        # Recompute total
        values[i, _COL_TOTAL] = values[i, _COL_ELASTIC] +
                                 values[i, _COL_FISSION] +
                                 values[i, _COL_CAPTURE]
    end
end

# Extract resonance range boundaries from MF2 data
function _resonance_bounds(mf2::MF2Data)
    eresl = Inf
    eresh = 0.0
    eresr = 0.0
    for iso in mf2.isotopes
        for rng in iso.ranges
            Int(rng.LRU) == 0 && continue
            eresl = min(eresl, rng.EL)
            eresh = max(eresh, rng.EH)
            if Int(rng.LRU) == 1
                eresr = max(eresr, rng.EH)
            end
        end
    end
    eresr = clamp(eresr, eresl, eresh)
    return (eresl, eresh, eresr)
end

# ==========================================================================
# Legacy interface functions
# ==========================================================================

"""
    sigma_mf2(E, mf2) -> CrossSections

Evaluate resonance cross sections at energy E by summing over all isotopes.
Legacy interface returning CrossSections struct.
"""
function sigma_mf2(E::Real, mf2::MF2Data)
    sig = CrossSections()
    E_f = Float64(E)
    for iso in mf2.isotopes
        abn = iso.ABN
        for rng in iso.ranges
            if E_f >= rng.EL && E_f < rng.EH && Int(rng.LRU) > 0
                try
                    sigp = cross_section(E_f, rng)
                    total = max(0.0, sigp.total)
                    elastic = max(0.0, sigp.elastic)
                    fission = max(0.0, sigp.fission)
                    capture = max(0.0, sigp.capture)
                    sig = sig + abn * CrossSections(total, elastic, fission, capture)
                catch e
                    @warn "sigma_mf2: cross_section failed at E=$E_f for range [$(rng.EL), $(rng.EH)]" exception=(e, catch_backtrace())
                end
            end
        end
    end
    return sig
end

"""
    merge_background_legacy(energies, res_xs, mf3_sections) -> Vector{CrossSections}

Legacy merge returning Vector{CrossSections}.
"""
function merge_background_legacy(energies::Vector{Float64},
                                  res_xs::Vector{CrossSections},
                                  mf3_sections::Vector{MF3Section})
    n = length(energies)
    result = Vector{CrossSections}(undef, n)
    for i in 1:n
        e = energies[i]
        elastic = res_xs[i].elastic
        fission = res_xs[i].fission
        capture = res_xs[i].capture
        for sec in mf3_sections
            mt = Int(sec.mt)
            (mt == 1 || mt == 3 || mt == 101) && continue
            bg = interpolate(sec.tab, e)
            bg == 0.0 && continue
            if mt == 2
                elastic += bg
            elseif mt == 18 || mt == 19
                fission += bg
            elseif mt == 102
                capture += bg
            end
        end
        if elastic <= 1.0e-8
            elastic = 1.0e-8
        end
        total = elastic + fission + capture
        result[i] = CrossSections(total, elastic, fission, capture)
    end
    return result
end
