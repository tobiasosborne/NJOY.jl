# Adaptive grid reconstruction -- generic higher-order function
#
# PROPOSAL B DESIGN: The adaptive grid is a GENERIC algorithm parameterized
# by any callable. It knows nothing about cross sections, ENDF, or nuclear
# physics. This enables reuse in BROADR, THERMR, or any other context
# requiring adaptive linearization.
#
# Key differences from a direct Fortran translation:
#   1. The function `f` is a type parameter -- any callable works.
#   2. Components are handled generically via NTuple{N,Float64}.
#   3. Pre-allocated workspace eliminates per-panel allocations.
#   4. Pure functions (round_sigfig) have no side effects.
#   5. The step-ratio guard and convergence tiers are cleanly separated.
#
# Correspondence to NJOY2016:
#   resxs (reconr.f90:2240-2569) -> adaptive_reconstruct
#   sigfig (util.f90:361-393)    -> round_sigfig

# ==========================================================================
# Configuration
# ==========================================================================

"""
    AdaptiveConfig

Configuration for adaptive grid refinement, matching NJOY2016's RECONR
tolerance parameters.

# Fields
- `err`: primary fractional tolerance
- `errmax`: relaxed tolerance (used with integral criterion)
- `errint`: integral error per grid point
- `max_depth`: maximum stack depth (NJOY uses 30)
- `thermal_threshold`: energy [eV] below which tolerance tightens
- `thermal_factor`: tightening divisor (NJOY uses 5.0 below 0.5 eV)
"""
struct AdaptiveConfig
    err::Float64
    errmax::Float64
    errint::Float64
    max_depth::Int
    thermal_threshold::Float64
    thermal_factor::Float64
    step_guard_limit::Float64
    force_boundaries::Vector{Float64}  # panels crossing these energies → forced convergence
end

function AdaptiveConfig(err::Real;
                        errmax::Real = 10 * err,
                        errint::Real = err / 20000,
                        max_depth::Int = 30,
                        thermal_threshold::Real = 0.4999,
                        thermal_factor::Real = 5.0,
                        step_guard_limit::Real = Inf,
                        force_boundaries::Vector{Float64} = Float64[])
    AdaptiveConfig(Float64(err), Float64(errmax), Float64(errint),
                   max_depth, Float64(thermal_threshold), Float64(thermal_factor),
                   Float64(step_guard_limit), force_boundaries)
end

# ==========================================================================
# Significant-figure rounding (pure function)
# ==========================================================================

"""
    round_sigfig(x, ndig, idig=0) -> Float64

Round `x` to `ndig` significant figures, then shade by `idig` in the last
significant digit. Matches NJOY2016 `sigfig` from util.f90.

- `idig = 0`: round to nearest
- `idig > 0`: shade up by `idig` units in the last significant digit
- `idig < 0`: shade down by `|idig|` units
"""
function round_sigfig(x::Real, ndig::Integer, idig::Integer = 0)
    x_f = Float64(x)
    x_f == 0.0 && return 0.0
    bias = 1.0000000000001  # tiny upward bias matching Fortran
    aa = log10(abs(x_f))
    ipwr = floor(Int, aa)
    ipwr = ndig - 1 - ipwr
    ii = round(Int, x_f * exp10(ipwr) + exp10(ndig - 11))
    if abs(ii) >= 10^ndig
        ii = div(ii, 10)
        ipwr -= 1
    end
    ii += Int(idig)
    return ii * exp10(-ipwr) * bias
end

# Efficient power-of-10 using a lookup or direct computation
@inline exp10(n::Int) = 10.0^n

# ==========================================================================
# Workspace (pre-allocated, avoids per-panel heap allocation)
# ==========================================================================

"""
    AdaptiveWorkspace{N}

Pre-allocated workspace for the adaptive stack. `N` is the number of
function components being tracked.
"""
struct AdaptiveWorkspace{N}
    stack_x::Vector{Float64}        # stack of energies
    stack_vals::Vector{NTuple{N, Float64}}  # stack of function values
end

function AdaptiveWorkspace{N}(max_depth::Int) where N
    AdaptiveWorkspace{N}(
        Vector{Float64}(undef, max_depth),
        Vector{NTuple{N, Float64}}(undef, max_depth)
    )
end

# ==========================================================================
# Core algorithm
# ==========================================================================

# Step-ratio guard (NJOY's estp=4.1)
const _STEP_RATIO = 4.1

"""
    adaptive_reconstruct(f, grid, config) -> (energies, values)

Adaptively refine `grid` until `f(x)` is linearly interpolable within
`config.err` relative tolerance for ALL components simultaneously.

# Arguments
- `f`: any callable `f(x::Float64)` returning either:
  - `NTuple{N, Float64}` (preferred, allocation-free)
  - `CrossSections` (converted to 4-tuple internally)
  - Any indexable container with `length` and `getindex`
- `grid`: sorted `Vector{Float64}` of initial energy nodes (>= 2 points)
- `config`: `AdaptiveConfig` with tolerance parameters

# Returns
`(energies::Vector{Float64}, values::Matrix{Float64})` where
`values[i, j]` is component `j` at `energies[i]`.

# Algorithm
Stack-based bisection matching NJOY2016's `resxs`:
1. For each pair of adjacent grid nodes, evaluate the midpoint.
2. Three-tier convergence test (primary, relaxed+integral, step-guard).
3. If converged, emit the upper stack point; otherwise bisect (push midpoint).
4. `round_sigfig` prevents infinite subdivision at floating-point limits.
"""
function adaptive_reconstruct(f, grid::AbstractVector{<:Real}, config::AdaptiveConfig)
    n_grid = length(grid)
    n_grid >= 2 || throw(ArgumentError("grid must have at least 2 points"))

    # Probe function to determine N (number of components)
    probe = _as_tuple(f(Float64(grid[1])))
    N = length(probe)

    # Pre-allocate workspace and output
    ws = AdaptiveWorkspace{N}(config.max_depth)
    out_e = Float64[]
    out_v = NTuple{N, Float64}[]
    sizehint!(out_e, max(4 * n_grid, 256))
    sizehint!(out_v, max(4 * n_grid, 256))

    n_added = 0
    n_integral = 0
    last_accepted_e = 0.0
    n_accepted = 0

    # Evaluate at first grid point
    egl = Float64(grid[1])
    val_left = probe

    for ig in 2:n_grid
        eg = Float64(grid[ig])
        val_right = _as_tuple(f(eg))

        n_added_panel, n_int_panel = _process_panel!(
            f, ws, out_e, out_v, config,
            egl, val_left, eg, val_right,
            last_accepted_e, n_accepted
        )
        n_added += n_added_panel
        n_integral += n_int_panel

        if !isempty(out_e)
            last_accepted_e = out_e[end]
            n_accepted = length(out_e)
        end

        egl = eg
        val_left = val_right
    end

    # Emit final point
    push!(out_e, egl)
    push!(out_v, val_left)

    # Convert to matrix
    n_out = length(out_e)
    values = Matrix{Float64}(undef, n_out, N)
    @inbounds for i in 1:n_out
        v = out_v[i]
        for j in 1:N
            values[i, j] = v[j]
        end
    end

    return (out_e, values)
end

"""
    _process_panel!(f, ws, out_e, out_v, config, e_lo, v_lo, e_hi, v_hi,
                    last_e, n_acc) -> (n_added, n_integral)

Process one panel [e_lo, e_hi] using stack-based bisection.
Emits points to `out_e` and `out_v` (excludes the final `e_hi` point,
which becomes the left endpoint of the next panel).
"""
function _process_panel!(f, ws::AdaptiveWorkspace{N}, out_e, out_v, config,
                         e_lo, v_lo, e_hi, v_hi,
                         last_accepted_e, n_accepted) where N
    sx = ws.stack_x
    sv = ws.stack_vals

    # Stack: index 1 = high energy end, index 2 = low energy end
    # Stack grows upward from index 2
    i = 2
    sx[2] = e_lo
    sv[2] = v_lo
    sx[1] = e_hi
    sv[1] = v_hi

    n_added = 0
    n_integral = 0
    local_last_e = last_accepted_e
    local_n_acc = n_accepted

    # Boundary-crossing forced convergence (reconr.f90:2353-2355).
    # Panels straddling resonance boundaries are forced to converge.
    _crosses_boundary = !isempty(config.force_boundaries)

    while i >= 2
        # Force convergence for panels crossing resonance boundaries
        if _crosses_boundary
            for eb in config.force_boundaries
                if sx[i] < eb && sx[i-1] > eb
                    @goto converged
                end
            end
        end

        xm = 0.5 * (sx[i] + sx[i-1])
        dx = sx[i-1] - sx[i]

        # --- Significant figures check ---
        ndig = (xm > 0.1 && xm < 1.0) ? 8 : 9

        if xm <= round_sigfig(sx[i], ndig, +1) || xm >= round_sigfig(sx[i-1], ndig, -1)
            # Cannot subdivide further -- forced convergence
            @goto converged
        end

        # Round midpoint
        if xm > round_sigfig(sx[i], 7, +1) && xm < round_sigfig(sx[i-1], 7, -1)
            xm = round_sigfig(xm, 7, 0)
        else
            xm = round_sigfig(xm, ndig, 0)
        end

        # Evaluate true function at midpoint
        mid_val = _as_tuple(f(xm))

        # Linear interpolation fractions
        fr2 = (xm - sx[i]) / dx
        fr1 = 1.0 - fr2

        # Effective tolerances (tighten near thermal)
        err_eff = config.err
        errmax_eff = config.errmax
        if sx[i-1] < config.thermal_threshold
            err_eff /= config.thermal_factor
            errmax_eff /= config.thermal_factor
        end

        # Compute deviations for all components
        all_primary = true
        any_exceeds_errmax = false
        @inbounds for j in 1:N
            sl = fr1 * sv[i][j] + fr2 * sv[i-1][j]
            dm = abs(mid_val[j] - sl)
            if dm > err_eff * abs(mid_val[j])
                all_primary = false
            end
            if dm > errmax_eff * abs(mid_val[j])
                any_exceeds_errmax = true
            end
        end

        if all_primary
            # Check step-ratio guard: prevent big jumps that miss peaks
            # Only apply when midpoint is below step_guard_limit (matching
            # Fortran's xm < eresr condition -- the guard is only needed in
            # the resolved resonance region).
            if local_n_acc > 3 && local_last_e > 0.0 && xm < config.step_guard_limit
                est = _STEP_RATIO * (sx[i] - local_last_e)
                if dx > est
                    @goto not_converged
                end
            end
            @goto converged

        elseif any_exceeds_errmax
            @goto not_converged

        else
            # Relaxed tolerance met -- apply integral criterion
            tsti = 2.0 * config.errint * xm / dx
            integral_ok = true
            @inbounds for j in 1:N
                sl = fr1 * sv[i][j] + fr2 * sv[i-1][j]
                dm = abs(mid_val[j] - sl)
                if dm >= tsti
                    integral_ok = false
                    break
                end
            end
            if integral_ok
                n_integral += 1
                @goto converged
            else
                @goto not_converged
            end
        end

        @label converged
        # Emit top of stack
        push!(out_e, sx[i])
        push!(out_v, sv[i])
        if abs(e_lo - sx[i]) > 1e-10 * e_lo
            n_added += 1
        end
        local_last_e = sx[i]
        local_n_acc += 1
        i -= 1
        continue

        @label not_converged
        # Push midpoint onto stack (bisect)
        i += 1
        if i > config.max_depth
            error("adaptive_reconstruct: stack depth $(config.max_depth) exceeded " *
                  "at E=$(xm). Increase max_depth or relax tolerance.")
        end
        sx[i] = sx[i-1]
        sv[i] = sv[i-1]
        sx[i-1] = xm
        sv[i-1] = mid_val
    end

    return (n_added, n_integral)
end

# ==========================================================================
# Type conversion helpers -- convert any return type to NTuple
# ==========================================================================

@inline _as_tuple(v::NTuple{N, Float64}) where N = v

@inline function _as_tuple(xs::CrossSections)
    (xs.total, xs.elastic, xs.fission, xs.capture)
end

# For generic indexable containers, build tuple dynamically
@inline function _as_tuple(v)
    N = length(v)
    ntuple(i -> Float64(v[i]), N)
end
