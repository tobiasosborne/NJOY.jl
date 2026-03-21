# GROUPR -- Group-averaged cross sections via exact piecewise integration
#
# PROPOSAL B DESIGN: Composable functions, exact integration, pluggable weights.
# The key insight: PENDF data is piecewise-linear, so group integrals are sums
# of analytical panel integrals. We reuse panel_integral from endf/interpolation.jl
# for all weighted integrals, giving exact results without quadrature.
#
# Correspondence to NJOY2016 groupr.f90:
#   gengpn         -> group_structures.jl (built-in structures)
#   sigma*flux     -> group_integrate (exact panel sums)
#   getflx+getwtf  -> weight functions are user-supplied callables
#   panel integ    -> panel_integral (from interpolation.jl, NJOY's gral)
#
# Separation of concerns:
#   1. group_integrate: raw integration of tabulated data over group bounds
#   2. group_average:   flux-weighted averaging (sigma_g = int(s*w)/int(w))
#   3. group_average_shielded: Bondarenko self-shielding on top

# ==========================================================================
# Result type
# ==========================================================================

"""
    MultiGroupXS

Result of group averaging: cross sections collapsed to a group structure.
"""
struct MultiGroupXS
    group_bounds::Vector{Float64}
    mt_list::Vector{Int}
    xs::Matrix{Float64}        # (n_groups, n_reactions)
    flux::Vector{Float64}      # group fluxes: integral of weight over each group
end

# Weight function aliases for backward compatibility.
# Canonical implementations are in weight_functions.jl.
const weight_flat = constant_weight
const weight_inv_e = inv_e_weight
weight_maxwell_fission(E::Real; kT::Real=0.0253) =
    maxwell_inv_e_fission(E; Eb=kT, Tb=kT)

# ==========================================================================
# Core: exact integration over group bounds
# ==========================================================================

"""
    group_integrate(energies, values, group_bounds; law=LinLin) -> Vector{Float64}

Integrate a piecewise-defined function (given by `energies` and `values` arrays)
over each group interval. Uses exact analytical panel integrals.

The data is treated as piecewise with the given interpolation `law` (default
linear-linear, appropriate for PENDF data). Each group integral is the sum of
panel_integral contributions from panels that overlap the group interval.

Returns a vector of length `length(group_bounds)-1`.
"""
function group_integrate(energies::AbstractVector{<:Real},
                         values::AbstractVector{<:Real},
                         bounds::Union{AbstractVector{<:Real}, NTuple};
                         law::InterpolationLaw=LinLin)
    ne = length(energies)
    @assert ne == length(values) && ne >= 2 "need at least 2 data points"
    @assert issorted(energies) "energies must be sorted"

    ng = length(bounds) - 1
    result = zeros(Float64, ng)

    # Walk through groups and data panels simultaneously
    # Panel index: data between energies[ip] and energies[ip+1]
    ip = 1  # current panel start index in data

    for g in 1:ng
        elo = Float64(bounds[g])
        ehi = Float64(bounds[g+1])

        # Advance panel pointer to first panel overlapping this group
        while ip < ne - 1 && energies[ip+1] <= elo
            ip += 1
        end

        jp = ip  # don't modify ip -- next group may reuse it
        s = 0.0
        while jp < ne && energies[jp] < ehi
            # Panel [energies[jp], energies[jp+1]] intersects [elo, ehi]
            xl = Float64(energies[jp])
            xh = Float64(energies[jp+1])
            yl = Float64(values[jp])
            yh = Float64(values[jp+1])

            # Clip to group bounds
            x1 = max(xl, elo)
            x2 = min(xh, ehi)

            if x2 > x1
                s += panel_integral(xl, yl, xh, yh, x1, x2, law)
            end
            jp += 1
        end
        result[g] = s
    end
    result
end

# ==========================================================================
# Weighted group averaging
# ==========================================================================

"""
    _build_weighted_tab(energies, values, weight_fn)

Compute product values[i]*weight_fn(energies[i]) at each data point.
Returns the product array. This allows exact integration of the product
since both factors are evaluated at the same knots.
"""
function _build_weighted_tab(energies::AbstractVector{<:Real},
                             values::AbstractVector{<:Real},
                             weight_fn)
    n = length(energies)
    prod = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        prod[i] = Float64(values[i]) * Float64(weight_fn(Float64(energies[i])))
    end
    prod
end

"""
    _integrate_weight(energies, weight_fn, bounds; law=LinLin) -> Vector{Float64}

Compute the integral of the weight function alone over each group.
Evaluates weight_fn at the data energy points to build a tabulation,
then integrates analytically.
"""
function _integrate_weight(energies::AbstractVector{<:Real},
                           weight_fn, bounds;
                           law::InterpolationLaw=LinLin)
    n = length(energies)
    wvals = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        wvals[i] = Float64(weight_fn(Float64(energies[i])))
    end
    group_integrate(energies, wvals, bounds; law=law)
end

"""
    group_average(energies, xs_cols, mt_list, group_bounds;
                  weight_fn=weight_inv_e, law=LinLin) -> MultiGroupXS

Compute flux-weighted group-averaged cross sections:

    σ_g = ∫_g σ(E)·w(E) dE  /  ∫_g w(E) dE

Arguments:
- `energies`: PENDF energy grid [eV], sorted ascending
- `xs_cols`: matrix of cross sections (n_energies, n_reactions) or single vector
- `mt_list`: reaction MT numbers corresponding to columns of xs_cols
- `group_bounds`: energy group boundaries (ascending), tuple or vector
- `weight_fn`: callable E -> weight value (default: 1/E)
- `law`: interpolation law for panel integrals (default: LinLin for PENDF)
"""
function group_average(energies::AbstractVector{<:Real},
                       xs_cols::AbstractMatrix{<:Real},
                       mt_list::AbstractVector{<:Integer},
                       bounds::Union{AbstractVector{<:Real}, NTuple};
                       weight_fn=weight_inv_e,
                       law::InterpolationLaw=LinLin)
    ne = length(energies)
    nr = length(mt_list)
    @assert size(xs_cols) == (ne, nr) "xs_cols shape must be (n_energies, n_reactions)"

    ng = length(bounds) - 1
    gb = collect(Float64, bounds)

    # Flux integral: ∫ w(E) dE over each group
    flux = _integrate_weight(energies, weight_fn, bounds; law=law)

    # Average each reaction
    xs_avg = Matrix{Float64}(undef, ng, nr)
    for r in 1:nr
        sigma_col = @view xs_cols[:, r]
        weighted = _build_weighted_tab(energies, sigma_col, weight_fn)
        num = group_integrate(energies, weighted, bounds; law=law)
        for g in 1:ng
            xs_avg[g, r] = flux[g] > 0.0 ? num[g] / flux[g] : 0.0
        end
    end

    MultiGroupXS(gb, collect(Int, mt_list), xs_avg, flux)
end

# Single-vector convenience
function group_average(energies::AbstractVector{<:Real},
                       xs::AbstractVector{<:Real},
                       mt::Integer,
                       bounds::Union{AbstractVector{<:Real}, NTuple};
                       weight_fn=weight_inv_e,
                       law::InterpolationLaw=LinLin)
    mat = reshape(collect(Float64, xs), :, 1)
    group_average(energies, mat, [mt], bounds; weight_fn=weight_fn, law=law)
end

# ==========================================================================
# Self-shielded group averaging (Bondarenko method)
# ==========================================================================

"""
    group_average_shielded(energies, total_xs, reaction_xs, mt_list,
                           group_bounds, sigma0;
                           weight_fn=weight_inv_e, law=LinLin) -> MultiGroupXS

Bondarenko self-shielded group averaging. The effective weight becomes:

    w_eff(E) = w(E) · σ₀ / (σ_t(E) + σ₀)

where σ_t is the total cross section and σ₀ is the background cross section.
At infinite dilution (σ₀ → ∞), this reduces to standard group_average.

Arguments:
- `energies`: PENDF energy grid
- `total_xs`: total cross section σ_t(E) on the grid
- `reaction_xs`: matrix of reaction cross sections (n_energies, n_reactions)
- `mt_list`: reaction MT numbers
- `group_bounds`: energy group boundaries
- `sigma0`: background cross section [barns]
- `weight_fn`: base weight function (default: 1/E)
"""
function group_average_shielded(energies::AbstractVector{<:Real},
                                total_xs::AbstractVector{<:Real},
                                reaction_xs::AbstractMatrix{<:Real},
                                mt_list::AbstractVector{<:Integer},
                                bounds::Union{AbstractVector{<:Real}, NTuple},
                                sigma0::Real;
                                weight_fn=weight_inv_e,
                                law::InterpolationLaw=LinLin)
    @assert length(total_xs) == length(energies)
    s0 = Float64(sigma0)

    # Build shielded weight: w(E) * sigma0 / (sigma_t(E) + sigma0)
    shielded_weight = let wf = weight_fn, st = total_xs, eg = energies, s = s0
        function(E::Real)
            # Find sigma_t by linear interpolation in the data
            idx = searchsortedlast(eg, E)
            idx = clamp(idx, 1, length(eg) - 1)
            t = (E - eg[idx]) / (eg[idx+1] - eg[idx])
            t = clamp(t, 0.0, 1.0)
            st_E = st[idx] * (1.0 - t) + st[idx+1] * t
            wf(E) * s / (st_E + s)
        end
    end

    group_average(energies, reaction_xs, mt_list, bounds;
                  weight_fn=shielded_weight, law=law)
end
