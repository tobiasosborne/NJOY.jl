# error_classifier.jl — Classify discrepancies between Julia and reference by error type
#
# Provides structured error analysis beyond simple relative error:
# identifies MISSING data, EXCESS peaks, MAGNITUDE errors, SHAPE errors,
# and GRID misses (resonance peaks the Julia grid missed entirely).

using Printf

const ABS_FLOOR_EC = 0.01  # barns — below this, XS is "near zero"

# =========================================================================
# Types
# =========================================================================

"""One classified error at a specific energy point."""
struct PointError
    energy::Float64
    julia_xs::Float64
    ref_xs::Float64
    rel_err::Float64
    error_type::Symbol    # :ok, :magnitude, :shape, :missing, :excess, :grid, :sign
end

"""Error classification for one reaction (MT) across all compared energies."""
struct ReactionErrorProfile
    mt::Int
    label::String
    n_compared::Int
    n_significant::Int        # where |ref_xs| >= ABS_FLOOR_EC

    # Classified error counts
    n_missing::Int            # julia near zero, ref significant
    n_excess::Int             # julia significant, ref near zero
    n_magnitude::Int          # both nonzero, differ >100%
    n_shape::Int              # both nonzero, differ 10-100%
    n_grid::Int               # resonance peak missed entirely

    # Summary statistics
    rms_rel_err::Float64
    max_rel_err::Float64
    worst_energy::Float64
    mean_abs_err::Float64

    # Top-N worst points
    worst_points::Vector{PointError}

    # Error regions: contiguous energy bands with dominant error type
    error_regions::Vector{@NamedTuple{e_lo::Float64, e_hi::Float64, type::Symbol, count::Int}}
end

# =========================================================================
# Point classification
# =========================================================================

"""
    classify_point(julia_xs, ref_xs; abs_floor) -> Symbol

Classify error type for a single energy point.
"""
function classify_point(julia_xs::Float64, ref_xs::Float64;
                        abs_floor::Float64=ABS_FLOOR_EC)
    ref_sig = abs(ref_xs) >= abs_floor
    jul_sig = abs(julia_xs) >= abs_floor

    # Both near zero → ok
    (!ref_sig && !jul_sig) && return :ok

    # Julia near zero, ref significant → missing data
    (!jul_sig && ref_sig) && return :missing

    # Julia significant, ref near zero → excess data
    (jul_sig && !ref_sig) && return :excess

    # Both significant → compute relative error
    rel = abs(julia_xs - ref_xs) / max(abs(ref_xs), 1e-30)
    rel <= 0.10 && return :ok
    rel <= 1.00 && return :shape
    return :magnitude
end

# =========================================================================
# Error profile builder
# =========================================================================

"""
    build_error_profile(our_e, our_xs, ref_e, ref_xs, mt, label; n_worst=20)

Build a full error profile for one reaction, comparing at reference energy points.
"""
function build_error_profile(our_e::Vector{Float64}, our_xs::Vector{Float64},
                             ref_e::Vector{Float64}, ref_xs::Vector{Float64},
                             mt::Int, label::String;
                             n_worst::Int=20,
                             abs_floor::Float64=ABS_FLOOR_EC)
    errors = PointError[]
    n_sig = 0
    sum_sq_rel = 0.0
    sum_abs = 0.0
    max_rel = 0.0
    worst_e = 0.0

    counts = Dict(:ok=>0, :missing=>0, :excess=>0, :magnitude=>0, :shape=>0, :grid=>0)

    for k in eachindex(ref_e)
        re = ref_e[k]
        (re <= 0 || !isfinite(re)) && continue

        xo = _interp_at_ec(our_e, our_xs, re)
        xr = ref_xs[k]
        etype = classify_point(xo, xr; abs_floor)

        rel = abs(xr) >= abs_floor ? abs(xo - xr) / max(abs(xr), 1e-30) : 0.0
        push!(errors, PointError(re, xo, xr, rel, etype))

        counts[etype] = get(counts, etype, 0) + 1

        if abs(xr) >= abs_floor
            n_sig += 1
            sum_sq_rel += rel^2
            sum_abs += abs(xo - xr)
            if rel > max_rel
                max_rel = rel
                worst_e = re
            end
        end
    end

    # Detect grid misses: ref has sharp peak in narrow window, Julia is flat
    n_grid = _detect_grid_misses!(errors, ref_e, ref_xs; abs_floor)
    counts[:grid] = n_grid

    # Sort by error for worst points
    sort!(errors; by=e -> e.rel_err, rev=true)
    worst = errors[1:min(n_worst, length(errors))]

    rms = n_sig > 0 ? sqrt(sum_sq_rel / n_sig) : 0.0
    mean_abs = n_sig > 0 ? sum_abs / n_sig : 0.0

    regions = identify_error_regions(errors)

    ReactionErrorProfile(
        mt, label, length(errors), n_sig,
        counts[:missing], counts[:excess], counts[:magnitude], counts[:shape], n_grid,
        rms, max_rel, worst_e, mean_abs,
        worst, regions
    )
end

# =========================================================================
# Grid miss detection
# =========================================================================

"""
Detect resonance peaks the Julia grid missed entirely.
A grid miss: reference has >10x variation in a narrow energy window
while Julia is essentially flat (< 2x variation at same points).
"""
function _detect_grid_misses!(errors::Vector{PointError},
                              ref_e::Vector{Float64}, ref_xs::Vector{Float64};
                              abs_floor::Float64=ABS_FLOOR_EC)
    n_grid = 0
    n = length(ref_e)
    n < 5 && return 0

    for i in 3:(n-2)
        abs(ref_xs[i]) < abs_floor && continue

        # Check if reference has a sharp peak: center > 10x neighbors
        left_avg = (abs(ref_xs[i-1]) + abs(ref_xs[i-2])) / 2
        right_avg = (abs(ref_xs[i+1]) + abs(ref_xs[i+2])) / 2
        neighbor_avg = max(left_avg, right_avg, abs_floor)

        if abs(ref_xs[i]) > 10 * neighbor_avg
            # Found a reference peak — check if Julia missed it
            idx = i  # errors[i] corresponds to ref point i
            if idx <= length(errors) && errors[idx].error_type in (:missing, :magnitude)
                errors[idx] = PointError(errors[idx].energy, errors[idx].julia_xs,
                                        errors[idx].ref_xs, errors[idx].rel_err, :grid)
                n_grid += 1
            end
        end
    end
    n_grid
end

# =========================================================================
# Error region identification
# =========================================================================

"""
    identify_error_regions(errors; merge_gap=0.1)

Group consecutive errors into contiguous energy regions by dominant error type.
`merge_gap` is the max relative energy gap for merging adjacent regions of same type.
"""
function identify_error_regions(errors::Vector{PointError};
                                merge_gap::Float64=0.1)
    # Sort by energy
    sorted = sort(filter(e -> e.error_type != :ok, errors); by=e -> e.energy)
    isempty(sorted) && return @NamedTuple{e_lo::Float64, e_hi::Float64, type::Symbol, count::Int}[]

    regions = @NamedTuple{e_lo::Float64, e_hi::Float64, type::Symbol, count::Int}[]
    cur_lo = sorted[1].energy
    cur_hi = sorted[1].energy
    cur_type = sorted[1].error_type
    cur_count = 1

    for i in 2:length(sorted)
        e = sorted[i]
        # Merge if same type and energies are close
        gap = cur_hi > 0 ? (e.energy - cur_hi) / cur_hi : Inf
        if e.error_type == cur_type && gap <= merge_gap
            cur_hi = e.energy
            cur_count += 1
        else
            push!(regions, (e_lo=cur_lo, e_hi=cur_hi, type=cur_type, count=cur_count))
            cur_lo = e.energy
            cur_hi = e.energy
            cur_type = e.error_type
            cur_count = 1
        end
    end
    push!(regions, (e_lo=cur_lo, e_hi=cur_hi, type=cur_type, count=cur_count))
    regions
end

# =========================================================================
# Formatting
# =========================================================================

"""Pretty-print a ReactionErrorProfile."""
function format_error_profile(prof::ReactionErrorProfile; indent::String="  ")
    buf = IOBuffer()
    @printf(buf, "%sMT=%-3d (%s): %d pts compared, %d significant\n",
            indent, prof.mt, prof.label, prof.n_compared, prof.n_significant)
    @printf(buf, "%s  RMS=%.2f%%  Max=%.2f%% at %.4e eV\n",
            indent, prof.rms_rel_err * 100, prof.max_rel_err * 100, prof.worst_energy)
    @printf(buf, "%s  Breakdown: %d ok, %d missing, %d excess, %d magnitude, %d shape, %d grid\n",
            indent,
            prof.n_compared - prof.n_missing - prof.n_excess - prof.n_magnitude - prof.n_shape - prof.n_grid,
            prof.n_missing, prof.n_excess, prof.n_magnitude, prof.n_shape, prof.n_grid)

    if !isempty(prof.error_regions)
        @printf(buf, "%s  Error regions:\n", indent)
        for r in prof.error_regions
            @printf(buf, "%s    %s in [%.3e, %.3e] eV (%d pts)\n",
                    indent, uppercase(string(r.type)), r.e_lo, r.e_hi, r.count)
        end
    end

    if !isempty(prof.worst_points) && prof.max_rel_err > 0.01
        n_show = min(5, length(prof.worst_points))
        @printf(buf, "%s  Top-%d worst points:\n", indent, n_show)
        @printf(buf, "%s    %14s  %14s  %14s  %10s  %s\n",
                indent, "E (eV)", "julia", "ref", "rel_err%", "type")
        for pe in prof.worst_points[1:n_show]
            @printf(buf, "%s    %14.6e  %14.6e  %14.6e  %10.2f  %s\n",
                    indent, pe.energy, pe.julia_xs, pe.ref_xs,
                    pe.rel_err * 100, pe.error_type)
        end
    end
    String(take!(buf))
end

# =========================================================================
# Interpolation helper (standalone, no dependency on reference_comparator)
# =========================================================================

"""Linear interpolation of xs at energy ep, given sorted energy grid e."""
function _interp_at_ec(e::Vector{Float64}, xs::Vector{Float64}, ep::Float64)
    n = length(e)
    n == 0 && return 0.0
    ep <= e[1] && return xs[1]
    ep >= e[n] && return xs[n]
    # Binary search
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >> 1
        e[mid] <= ep ? (lo = mid) : (hi = mid)
    end
    t = (ep - e[lo]) / max(e[hi] - e[lo], 1e-30)
    xs[lo] + t * (xs[hi] - xs[lo])
end
