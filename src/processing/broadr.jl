# Doppler broadening pipeline (BROADR)
# Kernel math lives in sigma1.jl; this file does orchestration.
# broadn_grid matches Fortran broadn (broadr.f90:1256-1508) exactly:
#   node selection, convergence stack, sigfig rounding, thnmax handling.
# Correspondence: broadr->doppler_broaden, broadn->broadn_grid,
#                 bsigma->sigma1_at, thinb->thin_xs

# ==========================================================================
# broadn_grid: Fortran-faithful adaptive grid for Doppler broadening
# ==========================================================================
# Matches NJOY2016 broadn (broadr.f90:1256-1508).
#
# Julia has all data in memory so paging collapses to a single call.
# In Fortran bfile3, the first call to broadn has klow=nlow=mpage+1
# (NOT 1), so broadn takes the klow!=1 path (lines 1329-1346):
#   - Places the FIRST grid point directly into stack position 2
#   - Computes slope direction from first interval
#   - Evaluates broadened XS at first point
#   - Advances k to nlow+1 and enters node-selection loop at label 120
#
# The klow==1 path (goto 110) is only used for SUBSEQUENT page calls
# where the save-preserved state from the previous call is reused.

# Fortran constants (broadr.f90:1298-1314)
const _BROADN_NSTACK  = 12       # max convergence stack depth
const _BROADN_NMAX    = 10       # max original grid points between nodes
const _BROADN_THERM   = 0.0253   # thermal energy [eV]
const _BROADN_ESTP    = 4.1      # step-ratio guard factor
const _BROADN_STEP    = 2.01     # max energy ratio between nodes
const _BROADN_RMAX    = 3.0      # monotonicity ratio limit
const _BROADN_ERRMIN  = 1.0e-15  # absolute floor for XS significance
const _BROADN_TRANGE  = 0.4999   # thermal tightening threshold [eV]
const _BROADN_SSMALL  = 1.0e-6   # skip reaction if |sn/stot| < this
const _BROADN_SMALL   = 1.0e-9   # relative tolerance for sigfig comparison

"""
    broadn_grid(energies, xs_matrix, alpha, errthn, errmax, errint, thnmax;
                emtr=nothing) -> (out_e, out_xs)

Construct the Doppler-broadened adaptive grid exactly matching Fortran broadn
(broadr.f90:1256-1508).

Walks the original grid selecting "nodes" via slope-direction tracking, then
refines between nodes using a convergence stack (depth up to 12). Above
`thnmax`, original unbroadened points are copied directly.

# Arguments
- `energies`: sorted energy grid [eV], length N (positive, no duplicates)
- `xs_matrix`: cross sections (N x nreac), columns are reactions
- `alpha`: Doppler width parameter = AWR / (bk * T_eff)
- `errthn`: primary fractional tolerance (= tol from user input)
- `errmax`: relaxed tolerance (typically 10*errthn)
- `errint`: integral error criterion (typically errthn/20000)
- `thnmax`: upper energy limit for broadening [eV]
- `emtr`: optional threshold energies per reaction (length nreac); below
          the threshold, broadened XS is zeroed. Default: all zeros.

# Returns
`(out_e::Vector{Float64}, out_xs::Matrix{Float64})` where `out_xs` is
(n_out x nreac).
"""
function broadn_grid(energies::AbstractVector{<:Real},
                     xs_matrix::AbstractMatrix{<:Real},
                     alpha::Float64, errthn::Float64,
                     errmax_tol::Float64, errint_tol::Float64,
                     thnmax::Float64;
                     emtr::Union{Nothing, AbstractVector{<:Real}} = nothing)
    N = length(energies)
    n_pts, nreac = size(xs_matrix)
    @assert n_pts == N >= 2

    seg_e  = Float64.(energies)
    seg_xs = Float64.(xs_matrix)  # N x nreac

    # Threshold energies for zeroing below threshold (Fortran emtr, line 1441)
    thresh = emtr === nothing ? zeros(Float64, nreac) : Float64.(emtr)

    # Convert energies to velocity space (broadr.f90:1318-1321)
    vel = sqrt.(alpha .* seg_e)

    # Output accumulators
    out_e  = Float64[]
    out_xs = Vector{Vector{Float64}}()
    sizehint!(out_e, 2 * N)
    sizehint!(out_xs, 2 * N)

    # Convergence stack arrays (broadr.f90:1299-1301)
    # ss[level] is a Vector{Float64} of length nreac.
    stack_es = zeros(Float64, _BROADN_NSTACK)
    stack_ss = [zeros(Float64, nreac) for _ in 1:_BROADN_NSTACK]
    stack_ks = zeros(Int, _BROADN_NSTACK)
    stack_js = zeros(Int, _BROADN_NSTACK)

    sn_buf = zeros(Float64, nreac)   # midpoint evaluation workspace
    dl     = zeros(Float64, nreac)   # slope direction per reaction

    j      = 0     # output point counter (Fortran j)
    tt1old = 0.0   # energy of previously emitted point (broadr.f90:1327)

    # ==================================================================
    # Helpers
    # ==================================================================

    function eval_broadened!(out::Vector{Float64}, E_ev::Float64)
        for i in 1:nreac
            out[i] = sigma1_at(E_ev, seg_e, @view(seg_xs[:, i]), alpha)
        end
    end

    function emit!(e_ev::Float64, xs_vec::Vector{Float64})
        # Fortran lines 1439-1449: zero XS below threshold
        xs_out = Vector{Float64}(undef, nreac)
        for i in 1:nreac
            if e_ev > thresh[i] && tt1old >= thresh[i]
                xs_out[i] = xs_vec[i]
            else
                xs_out[i] = 0.0
            end
        end
        push!(out_e, e_ev)
        push!(out_xs, xs_out)
    end

    function update_slope_dir!(dl::Vector{Float64}, k::Int)
        # broadr.f90:1334-1339 and 1458-1464
        for i in 1:nreac
            dn = seg_xs[k+1, i] - seg_xs[k, i]
            if abs(dn) < abs(seg_xs[k, i]) / 1000
                dn = 0.0
            end
            dl[i] = dn < 0.0 ? -1.0 : 1.0
        end
    end

    function affected_flag(xs_broadened::Vector{Float64}, k::Int)::Int
        # broadr.f90:1341-1344
        for i in 1:nreac
            if abs(xs_broadened[i] - seg_xs[k, i]) > errthn * seg_xs[k, i]
                return 1
            end
        end
        return 0
    end

    # ==================================================================
    # SETUP: Place first grid point into stack position 2
    # (broadr.f90:1329-1346, the klow!=1 path)
    # ==================================================================
    k = 1                                                        # line 1323: k=nlow

    # Lines 1329-1346: set up first node in stack
    klast = k                                                    # line 1329
    stack_ks[2] = k                                              # line 1330
    et = vel[k]^2 / alpha                                        # line 1331
    stack_es[2] = round_sigfig(et, 7, 0)                         # line 1332
    xt_ev = stack_es[2]

    # Slope direction from first interval (lines 1334-1339)
    update_slope_dir!(dl, k)

    # Broadened XS at first node (line 1340)
    eval_broadened!(stack_ss[2], xt_ev)

    # Affected flag (lines 1341-1344)
    stack_js[2] = affected_flag(stack_ss[2], k)

    k += 1                                                       # line 1345: k=k+1

    # Now fall through to label 120 (goto 120, line 1346)

    # ==================================================================
    # MAIN LOOP: state machine matching Fortran goto structure
    # ==================================================================
    # States:
    #   :select_node  -> label 120 (node selection)
    #   :converge     -> label 140 (convergence stack)
    #   :copy_tail    -> label 190 (copy above thnmax)
    #   :done         -> return

    state = :select_node

    while state != :done
        if state == :select_node
            # ---- Node selection (label 120, lines 1351-1374) ----
            while true
                k > N && break

                et = vel[k]^2 / alpha                           # line 1352

                # Forced node: at/above thnmax or at end of grid
                (et >= thnmax) && break                          # line 1353
                (k >= N)      && break                           # line 1354

                # Skip: rounds to same 7-sigfig as prev node (line 1356 -> 126)
                test7 = round_sigfig(et, 7, 0)
                if abs(stack_es[2] - test7) < _BROADN_SMALL * test7
                    k += 1                                       # goto 126; k++; goto 120
                    continue
                end

                # Forced: too many skipped points (line 1357)
                (k > klast + _BROADN_NMAX) && break

                # Forced: energy jump too large (line 1358)
                (et > _BROADN_STEP * stack_es[2]) && break

                # Forced: at 3-sigfig boundary (line 1360)
                test3 = round_sigfig(et, 3, 0)
                (abs(et - test3) < _BROADN_SMALL * test3) && break

                # Forced: at thermal energy (line 1361)
                (abs(et - _BROADN_THERM) < _BROADN_SMALL * _BROADN_THERM) && break

                # Check slope direction change (lines 1362-1371)
                slope_changed = false
                for i in 1:nreac
                    dn = seg_xs[k+1, i] - seg_xs[k, i]
                    if abs(dn) < abs(seg_xs[k, i]) / 1000
                        dn = 0.0
                    end
                    dir = dn >= 0.0 ? 1.0 : -1.0
                    if dir != dl[i]
                        slope_changed = true
                        break
                    end
                end
                slope_changed && break                           # goto 130

                # All skip tests passed (goto 126; k++; goto 120)
                k += 1
            end

            # ---- Node found at index k (label 130, lines 1375-1383) ----
            et = vel[k]^2 / alpha
            stack_es[1] = round_sigfig(et, 7, 0)                # line 1376
            xt_ev = stack_es[1]
            eval_broadened!(stack_ss[1], xt_ev)                  # line 1378
            stack_ks[1] = k                                      # line 1379
            stack_js[1] = affected_flag(stack_ss[1], k)          # lines 1380-1383

            state = :converge

        elseif state == :converge
            # ---- Convergence stack (label 140, lines 1386-1488) ----
            is = 2
            while true
                converged = false

                # Test 1: stack full -> forced convergence (line 1388)
                if is >= _BROADN_NSTACK
                    converged = true
                end

                # Test 2: adjacent original points, both unaffected (lines 1389-1393)
                if !converged
                    if stack_ks[is-1] == stack_ks[is] + 1 &&
                       stack_js[is-1] == 0 && stack_js[is] == 0
                        converged = true
                    end
                end

                # Midpoint tests (lines 1395-1436)
                em = 0.0  # will be set by _broadn_test_midpoint if needed
                km = 0
                if !converged
                    em, km, converged = _broadn_test_midpoint(
                        is, nreac, stack_es, stack_ss, stack_ks, stack_js,
                        sn_buf, seg_e, seg_xs, alpha, errthn, errmax_tol,
                        errint_tol, j, out_e)
                end

                if converged
                    # ---- Emit top of stack (label 150, lines 1438-1451) ----
                    emit!(stack_es[is], stack_ss[is])
                    tt1old = stack_es[is]
                    j += 1

                    is -= 1
                    if is > 1
                        continue   # line 1451: goto 140
                    end

                    # Stack exhausted: promote node (lines 1453-1470)
                    stack_es[2] = stack_es[1]
                    stack_ks[2] = stack_ks[1]
                    stack_js[2] = stack_js[1]
                    stack_ss[2] .= stack_ss[1]

                    # Update slope direction (lines 1458-1464)
                    if k < N
                        for i in 1:nreac
                            dn = seg_xs[k+1, i] - seg_xs[k, i]
                            if abs(dn) < abs(seg_xs[k, i]) / 1000
                                dn = 0.0
                            end
                            dl[i] = dn >= 0.0 ? 1.0 : -1.0
                        end
                    end

                    # Termination checks
                    if k >= N                       # line 1466 -> label 180
                        state = :done
                        break
                    end
                    if stack_es[1] > thnmax         # line 1467 -> label 190
                        state = :copy_tail
                        break
                    end

                    klast = k       # line 1468
                    k += 1          # line 1469
                    state = :select_node   # goto 120
                    break

                else
                    # ---- Not converged: push midpoint (label 170, lines 1472-1488) ----
                    is += 1
                    stack_es[is]   = stack_es[is-1]
                    stack_es[is-1] = em
                    stack_ks[is]   = stack_ks[is-1]
                    stack_ks[is-1] = km
                    stack_js[is]   = stack_js[is-1]
                    # km==0 means interpolated midpoint, not original (line 1480)
                    stack_js[is-1] = (km == 0) ? 1 : 0
                    for i in 1:nreac
                        stack_ss[is][i]   = stack_ss[is-1][i]   # line 1482
                        stack_ss[is-1][i] = sn_buf[i]           # line 1483
                    end
                    if km != 0
                        for i in 1:nreac
                            if abs(sn_buf[i] - seg_xs[km, i]) >
                               errthn * seg_xs[km, i] + _BROADN_ERRMIN
                                stack_js[is-1] = 1
                            end
                        end
                    end
                    continue   # goto 140
                end
            end  # convergence while-loop

        elseif state == :copy_tail
            # ---- Copy above thnmax (label 190, lines 1495-1507) ----
            while k <= N
                et = vel[k]^2 / alpha
                xs_copy = Vector{Float64}(undef, nreac)
                for i in 1:nreac
                    xs_copy[i] = seg_xs[k, i]
                end
                push!(out_e, et)
                push!(out_xs, xs_copy)
                k += 1
            end
            state = :done
        end
    end  # state machine

    # Assemble output matrix
    n_out = length(out_e)
    result_xs = Matrix{Float64}(undef, n_out, nreac)
    for p in 1:n_out, i in 1:nreac
        result_xs[p, i] = out_xs[p][i]
    end
    return (out_e, result_xs)
end

"""
    _broadn_test_midpoint(...) -> (em, km, converged)

Compute and test the midpoint between stack levels `is` and `is-1`.
Returns the midpoint energy `em`, its original-grid index `km` (0 if
interpolated), and whether convergence was achieved.

Matches broadr.f90 lines 1395-1436.
"""
function _broadn_test_midpoint(is, nreac,
        stack_es, stack_ss, stack_ks, stack_js,
        sn_buf, seg_e, seg_xs, alpha,
        errthn, errmax_tol, errint_tol,
        j, out_e)

    # ---- Midpoint with sigfig rounding (lines 1395-1404) ----
    em = 0.5 * (stack_es[is-1] + stack_es[is])

    ndig = 9
    if em > 0.1 && em < 1.0                                    # line 1397
        ndig = 8
    end

    # Rounding: if midpoint is well within the 7-sigfig range of both
    # endpoints, round to 7 digits; otherwise round to ndig (line 1398-1404)
    if em > round_sigfig(stack_es[is], 7, +1)
        if em < round_sigfig(stack_es[is-1], 7, -1)
            em = round_sigfig(em, 7, 0)
        end
    else
        em = round_sigfig(em, ndig, 0)
    end

    km = 0                                                       # line 1405

    # Cannot subdivide further (lines 1406-1407)
    if em < round_sigfig(stack_es[is], ndig, +1) ||
       em > round_sigfig(stack_es[is-1], ndig, -1)
        return (em, km, true)   # forced convergence
    end

    # ---- Thermal tightening (lines 1408-1411) ----
    errt = errthn
    errm = errmax_tol
    if stack_es[is-1] < _BROADN_TRANGE
        errt /= 5
        errm /= 5
    end

    # ---- Evaluate broadened XS at midpoint (line 1413) ----
    for i in 1:nreac
        sn_buf[i] = sigma1_at(em, seg_e, @view(seg_xs[:, i]), alpha)
    end

    # ---- Linear interpolation fraction (lines 1414-1417) ----
    dx = stack_es[is-1] - stack_es[is]
    f  = (em - stack_es[is]) / dx

    # Fraction too close to 1 -> force refinement (line 1417: test=1-1/100)
    if f > 1.0 - 1.0/100
        return (em, km, false)
    end

    # ---- Sum of midpoint XS for significance test (line 1418-1421) ----
    stot = 0.0
    for i in 1:nreac
        stot += sn_buf[i]
    end

    # ---- Per-reaction convergence tests (lines 1422-1431) ----
    for i in 1:nreac
        # Skip negligible total (line 1423)
        stot < _BROADN_ERRMIN && continue
        # Skip negligible component (line 1424)
        abs(sn_buf[i] / stot) < _BROADN_SSMALL && continue

        # Monotonicity ratio (lines 1425-1426)
        if stack_ss[is-1][i] > _BROADN_RMAX * stack_ss[is][i]
            return (em, km, false)
        end
        if stack_ss[is-1][i] < stack_ss[is][i] / _BROADN_RMAX
            return (em, km, false)
        end

        # Linear interpolation error (lines 1427-1431)
        si = f * stack_ss[is-1][i] + (1 - f) * stack_ss[is][i]
        dy = abs(sn_buf[i] - si)

        # Primary tolerance (line 1429)
        dy <= errt * abs(sn_buf[i]) + _BROADN_ERRMIN && continue

        # Exceeds relaxed tolerance (line 1430)
        if dy > errm * abs(sn_buf[i]) + _BROADN_ERRMIN
            return (em, km, false)
        end

        # Integral criterion (line 1431)
        if dy * dx / 2 > errint_tol * em
            return (em, km, false)
        end
    end

    # ---- Step-size guard (lines 1435-1436) ----
    # tt(1) is last emitted energy; j is output count
    if j > 3 && !isempty(out_e)
        est = _BROADN_ESTP * (stack_es[is] - out_e[end])
        if dx > est
            return (em, km, false)
        end
    end

    return (em, km, true)  # converged
end

# ==========================================================================
# Public API: doppler_broaden using broadn_grid
# ==========================================================================

"""
    doppler_broaden(energies, xs, T, awr; tol=0.001, ...) -> (new_e, new_xs)
Doppler-broaden a single piecewise-linear cross section from 0K to temperature T.
Uses SIGMA1 exact kernel + broadn_grid adaptive refinement.
"""
function doppler_broaden(energies::AbstractVector{<:Real},
                         xs::AbstractVector{<:Real},
                         T::Real, awr::Real;
                         tol::Real = 0.001, do_thin::Bool = true,
                         errmax::Real = 10*tol, errint::Real = tol/20000,
                         thnmax::Real = Inf)
    n = length(energies)
    @assert n == length(xs) && n >= 2 && T > 0 && awr > 0
    @assert issorted(energies) "energies must be sorted"
    alpha = awr / (PhysicsConstants.bk * T)
    seg_e = Float64.(energies); seg_xs = reshape(Float64.(xs), n, 1)
    effective_thnmax = isfinite(thnmax) ? Float64(thnmax) : seg_e[end] * 1.01
    out_e, out_v = broadn_grid(seg_e, seg_xs, alpha,
                               Float64(tol), Float64(errmax), Float64(errint),
                               effective_thnmax)
    new_xs = vec(out_v)
    do_thin && ((out_e, new_xs) = thin_xs(out_e, new_xs; tol=Float64(tol),
                                           thnmax=effective_thnmax))
    return (out_e, new_xs)
end

"""
    doppler_broaden_multi(energies, xs_matrix, T, awr; ...) -> (new_e, new_xs)
Broaden multiple cross sections simultaneously on a shared adaptive grid.
`xs_matrix` is (N, n_reactions). Returns matrix output.
"""
function doppler_broaden_multi(energies::AbstractVector{<:Real},
                               xs_matrix::AbstractMatrix{<:Real},
                               T::Real, awr::Real;
                               tol::Real = 0.001, do_thin::Bool = true,
                               errmax::Real = 10*tol, errint::Real = tol/20000,
                               thnmax::Real = Inf)
    n_pts, n_reac = size(xs_matrix)
    @assert n_pts == length(energies) && T > 0
    alpha = awr / (PhysicsConstants.bk * T)
    seg_e = Float64.(energies)
    effective_thnmax = isfinite(thnmax) ? Float64(thnmax) : seg_e[end] * 1.01
    out_e, out_v = broadn_grid(seg_e, Float64.(xs_matrix), alpha,
                               Float64(tol), Float64(errmax), Float64(errint),
                               effective_thnmax)
    do_thin && ((out_e, out_v) = thin_xs(out_e, out_v; tol=Float64(tol),
                                          thnmax=effective_thnmax))
    return (out_e, out_v)
end

"""
    doppler_broaden(pendf::PointwiseMaterial, T_new; T_old=0.0, awr=1.0,
                    tol=0.001, thnmax=Inf) -> PointwiseMaterial
Broaden a PointwiseMaterial from T_old to T_new.
"""
function doppler_broaden(pendf::PointwiseMaterial, T_new::Float64;
                         T_old::Float64=0.0, awr::Float64=1.0,
                         tol::Float64=0.001, thnmax::Float64=Inf)
    T_eff = T_new - T_old
    T_eff == 0.0 && return pendf
    T_eff < 0.0 && error("doppler_broaden: T_new must be >= T_old")
    emax = isfinite(thnmax) && thnmax > 0.0 ? thnmax : pendf.energies[end]
    idx_broad = findall(e -> e <= emax, pendf.energies)
    idx_copy  = findall(e -> e >  emax, pendf.energies)
    (isempty(idx_broad) || length(idx_broad) < 2) && return pendf
    out_e, out_v = doppler_broaden_multi(
        pendf.energies[idx_broad], pendf.cross_sections[idx_broad, :],
        T_eff, awr; tol=tol, do_thin=true, thnmax=thnmax)
    if !isempty(idx_copy)
        out_e = vcat(out_e, pendf.energies[idx_copy])
        out_v = vcat(out_v, pendf.cross_sections[idx_copy, :])
    end
    return PointwiseMaterial(pendf.mat, out_e, out_v, copy(pendf.mt_list))
end

# ==========================================================================
# Thinning (matches NJOY2016 thinb, broadr.f90:1677-1760)
# ==========================================================================

"""
    thin_xs(energies, xs; tol=0.001, step_max=1.24, thnmax=Inf) -> (thinned_e, thinned_xs)
Remove redundant points where linear interpolation suffices within `tol`.
Above `thnmax`, points are never thinned (matching broadr.f90:1724).
"""
function thin_xs(energies::AbstractVector{<:Real}, xs::AbstractVector{<:Real};
                 tol::Real=0.001, step_max::Real=1.24, thnmax::Real=Inf)
    n = length(energies); @assert n == length(xs)
    n <= 2 && return (collect(Float64, energies), collect(Float64, xs))
    keep = _compute_thin_mask(energies, reshape(Float64.(xs), n, 1),
                              Float64(tol), Float64(step_max), Float64(thnmax))
    idx = findall(keep)
    return (Float64.(energies[idx]), Float64.(xs[idx]))
end

function thin_xs(energies::AbstractVector{<:Real}, xs::AbstractMatrix{<:Real};
                 tol::Real=0.001, step_max::Real=1.24, thnmax::Real=Inf)
    n = size(xs, 1); @assert n == length(energies)
    n <= 2 && return (collect(Float64, energies), Float64.(xs))
    keep = _compute_thin_mask(energies, Float64.(xs), Float64(tol), Float64(step_max),
                              Float64(thnmax))
    idx = findall(keep)
    return (Float64.(energies[idx]), Float64.(xs[idx, :]))
end

# Core thinning: test ALL reactions simultaneously for union grid preservation.
# Matches NJOY2016 thinb (broadr.f90:1677-1760).
# - tol*xs_matrix[j,r] WITHOUT abs() matches Fortran's errthn*s(i,k) -- negative
#   XS values effectively disable thinning at that point.
# - Reciprocal multiply matches Fortran's single-reciprocal division pattern.
# - thnmax guard: above thnmax, points are never thinned (broadr.f90:1724).
function _compute_thin_mask(energies, xs_matrix::Matrix{Float64},
                            tol::Float64, step_max::Float64,
                            thnmax::Float64)
    n, nr = size(xs_matrix)
    keep = falses(n); keep[1] = true; keep[n] = true
    last_kept = 1
    for i in 3:n
        can_thin = true
        e_lo = energies[last_kept]; denom = energies[i] - e_lo
        denom <= 0.0 && continue
        denom = 1.0 / denom                                        # reciprocal multiply (Fortran pattern)
        for j in (last_kept+1):(i-1)
            energies[j] >= step_max*e_lo && (can_thin = false; break)
            energies[j] * (1.0 + 1e-5) > thnmax && (can_thin = false; break)  # broadr.f90:1724
            for r in 1:nr
                slope = (xs_matrix[i,r] - xs_matrix[last_kept,r]) * denom  # reciprocal multiply
                sp = xs_matrix[last_kept,r] + slope*(energies[j] - e_lo)
                abs(sp - xs_matrix[j,r]) > tol*xs_matrix[j,r] && (can_thin = false; break)
            end
            can_thin || break
        end
        if !can_thin; keep[i-1] = true; last_kept = i-1; end
    end
    return keep
end

# ==========================================================================
# Grid enrichment (kept for backward compatibility with old API)
# ==========================================================================

function _prepare_grid(seg_e, seg_xs)
    grid = filter(e -> e > 0.0, seg_e); sort!(grid); unique!(grid)
    _enrich_broadr_grid(grid, seg_e, seg_xs)
end

function _enrich_broadr_grid(grid::Vector{Float64}, seg_e, seg_xs)
    enriched = copy(grid); n = length(seg_e)
    n <= 2 && return enriched
    for i in 2:(n-1)
        de_lo = seg_e[i] - seg_e[i-1]; de_hi = seg_e[i+1] - seg_e[i]
        (de_lo <= 0.0 || de_hi <= 0.0) && continue
        sl = (seg_xs[i] - seg_xs[i-1])/de_lo; sh = (seg_xs[i+1] - seg_xs[i])/de_hi
        ms = max(abs(sl), abs(sh))
        if ms > 0 && abs(sh - sl) > 0.5*ms
            push!(enriched, 0.5*(seg_e[i-1] + seg_e[i]))
            push!(enriched, 0.5*(seg_e[i] + seg_e[i+1]))
        end
    end
    sort!(enriched); unique!(enriched)
    return enriched
end
