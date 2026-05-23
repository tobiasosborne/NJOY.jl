# =========================================================================
# rescon — MF=32 → MF=33 resonance-parameter covariance sandwich
#
# Mirrors Fortran covout (errorr.f90:7465) → rescon (errorr.f90:8513-8819)
# + the rpxlc12 sensitivity build (errorr.f90:4399-4593).
#
# Algorithm (LRU=1 / LRF=3 / LCOMP=1, the U-238 JENDL path):
#   1. Read MF=32 → typed structs (mf32_reader.jl).
#   2. Read MF=2 → ReichMooreParameters (existing read_mf2).
#   3. For each MF=32 subsection (per range):
#      a. Build a dense pointwise E-grid that resolves each resonance
#         in [elr, ehr] (rpendf-equivalent — errorr.f90:5015-5089).
#      b. Match each MF=32 resonance (loopm) to its MF=2 (l, n) via
#         |ER/ER2 - 1| < 1e-6 && |AJ - AJ2| < 1e-4
#         (errorr.f90:4413-4443).
#      c. For each (loopm, loopn) ∈ NRB×MPAR perturbation:
#           · loopn=1 → ER  perturbed by ±0.0001·ER (gwidth=ER·0.0001)
#           · loopn≥2 → width perturbed by ±0.01·width (gwidth=w·0.01)
#         Two cross_section_rm evaluations per perturbation.
#         Central difference: dσ/dRP = (σ_+ - σ_-) / (2·gwidth).
#         Group-average via the iwt-aware Simpson formula from
#         `_rpxgrp_average` (rpxgrp port — errorr.f90:5227-5366), with
#         the per-E weight function selected by the user-deck `iwt`
#         (egtwtf, errorr.f90:10023-10110). T15 uses iwt=6 (thermal +
#         1/E + fission + fusion). Phase 72c fix: previously this
#         module assumed iwt=2 (flat) which biased C[1,1] by ~16%.
#         Stash sens[channel ∈ 1..4, p ∈ 1..npar, ig ∈ 1..ngn].
#   4. Sandwich (errorr.f90:4541-4593): for each (mt, mt2) ∈ RESCON_PAIRS:
#         absolute_cov[ig, ig'] = Σ_{i ≤ j} sens[ch,i,ig]·cov[i,j]·sens[ch,j,ig']
#                                + cross-doubled term when i ≠ j
#      then divide by (σ̄(ig)·σ̄(ig')) to get RELATIVE cov (the form
#      cov_matrices already stores after Phase 51 weighted-collapse).
#   5. Add into cov_matrices[(mt, mt2)] for the seven pairs.
#
# Channel mapping (errorr.f90:4550-4593):
#   sens[1,*,*] = total   → ctt (1,1)
#   sens[2,*,*] = elastic → cee/cef/ceg (2,2)/(2,18)/(2,102)
#   sens[3,*,*] = fission → cff/cef/cfg (18,18)/(2,18)/(18,102)
#   sens[4,*,*] = capture → cgg/ceg/cfg (102,102)/(2,102)/(18,102)
#
# Status: Phase 72c — T15 MT=102 row-1 canary GREEN.
# Acceptance: |jul_mt102[1,1] - 2.658914e-4| < 1e-7,  jul_mt102[1,9] < 0.
# =========================================================================

using LinearAlgebra: norm

# Seven (mt, mt2) reaction pairs that receive an RP-cov contribution.
# Mirrors errorr.f90:8531-8539 (rescon itp dispatch). Order matches the
# Fortran branch order so future sensitivity-builder code can dispatch
# by index 1..7 directly.
#   itp=1: (18,18)   cff triangular
#   itp=2: (18,102)  cfg full
#   itp=3: (102,102) cgg triangular  ← U-238 canary
#   itp=4: (2,2)     cee triangular
#   itp=5: (2,18)    cef full
#   itp=6: (2,102)   ceg full
#   itp=7: (1,1)     ctt triangular
const RESCON_PAIRS = (
    (18,  18 ),
    (18,  102),
    (102, 102),
    (2,   2  ),
    (2,   18 ),
    (2,   102),
    (1,   1  ),
)

# Channel index (1=total, 2=elastic, 3=fission, 4=capture) used as the
# `sens` first axis, per errorr.f90:4550-4593.
function _rescon_channel_pair(mt::Int, mt2::Int)
    canonical = mt <= mt2 ? (mt, mt2) : (mt2, mt)
    canonical == (1,   1  ) && return (1, 1)  # total/total
    canonical == (2,   2  ) && return (2, 2)  # elastic/elastic
    canonical == (2,   18 ) && return (2, 3)  # elastic/fission
    canonical == (2,   102) && return (2, 4)  # elastic/capture
    canonical == (18,  18 ) && return (3, 3)  # fission/fission
    canonical == (18,  102) && return (3, 4)  # fission/capture
    canonical == (102, 102) && return (4, 4)  # capture/capture
    return (0, 0)
end

"""
    rescon_pair_index(mt::Int, mt2::Int) -> Int

Map a `(mt, mt2)` pair to its rescon `itp` index 1..7, or 0 if the pair
does not receive an RP-cov contribution. Mirrors the seven-way `if`
chain at errorr.f90:8531-8539. Order is `mt ≤ mt2` (canonical).
"""
function rescon_pair_index(mt::Int, mt2::Int)
    canonical = mt <= mt2 ? (mt, mt2) : (mt2, mt)
    for (i, p) in enumerate(RESCON_PAIRS)
        canonical == p && return i
    end
    return 0
end

rescon_supports_pair(mt::Int, mt2::Int) = rescon_pair_index(mt, mt2) != 0

# -------------------------------------------------------------------------
# Pointwise E-grid for one resolved range (rpendf-equivalent)
# -------------------------------------------------------------------------

# Build a pointwise E-grid that resolves each resonance in [elr, ehr].
# Fortran's rpendf (errorr.f90:5015-5089) uses an adaptive multiplicative
# stepper with eskip1/2/3/4 factors keyed off proximity to each ER. We
# build a static union of:
#   - a base log-spaced grid from elr → ehr (~50 points/decade)
#   - dense linear refinement around each ER, from ER−30·Γ to ER+30·Γ
#     with ≥40 points (Γ = max(Gn+Gg, 1e-3·ER))
#   - all output-group boundaries within [elr, ehr] (clean integration
#     splits)
# Sorted, deduped, clipped to [elr, ehr].
function _build_pointwise_grid(elr::Float64, ehr::Float64,
                                resonances_E::AbstractVector{<:Real},
                                resonances_W::AbstractVector{<:Real},
                                egn::AbstractVector{<:Real})
    pts = Float64[]
    sizehint!(pts, 4000)

    # Base log grid (~200 pts/decade).
    elo = max(elr, 1e-5)
    nlog = max(200, ceil(Int, 200 * log10(ehr / elo)))
    for k in 0:nlog
        push!(pts, elo * (ehr / elo)^(k / nlog))
    end

    # Per-resonance refinement: 200 points spanning ±30Γ, with extra
    # density inside ±Γ (Lorentzian peak) via a tanh-stretched mapping.
    # ~33 pts/Γ in the wings, ~6 pts/(Γ/10) at the peak.
    for (Er, w) in zip(resonances_E, resonances_W)
        Er <= 0 && continue
        Γ = max(w, 1e-3 * abs(Er))
        e_lo = max(Er - 30 * Γ, elr)
        e_hi = min(Er + 30 * Γ, ehr)
        e_lo >= e_hi && continue
        npts = 200
        for k in 0:npts
            # Tanh-stretched: tau ∈ [-1, 1], mapped to [e_lo, e_hi]
            # with extra density at center (tau=0). atanh-inverse gives
            # uniform sampling in tanh space, dense at peak.
            tau = -1.0 + 2.0 * k / npts
            # Use a power mapping that concentrates near 0:
            #   E = Er + sign(tau) · (|tau|^p) · half_width
            # p > 1 → denser at center.
            half = 0.5 * (e_hi - e_lo)
            mid  = 0.5 * (e_hi + e_lo)
            p_pow = 2.0
            offset = sign(tau) * abs(tau)^p_pow * half
            push!(pts, mid + offset)
        end
    end

    # Output-group boundaries within range.
    for e in egn
        elr <= e <= ehr && push!(pts, Float64(e))
    end

    sort!(pts)
    unique!(pts)
    # Strip outside [elr, ehr] (clip).
    filter!(e -> elr <= e <= ehr, pts)
    isempty(pts) && error("_build_pointwise_grid: empty grid for [$elr, $ehr]")
    pts
end

# -------------------------------------------------------------------------
# Weight-aware Simpson group-average (rpxgrp port)
# -------------------------------------------------------------------------

"""
    _rpxgrp_average(grid, sig, egn, weight_fn) -> Vector{Float64}

Group-average a pointwise vector using the iwt-aware Simpson formula
from Fortran NJOY's `rpxgrp` (errorr.f90:5227-5366).

For each segment `(x2, x1)` fully inside the current group:

    z1     = (2·wt1 + wt2)·dx/6   # endpoint x1 coefficient
    z2     = (2·wt2 + wt1)·dx/6   # endpoint x2 coefficient
    de     = ½·(wt1 + wt2)·dx
    gsig  += y1·z1 + y2·z2        (per channel)
    sumde += de

At group boundaries: linear interpolation of both `wt` and `y` at the
crossing energy `ebb = egn(ig+1)`, then accumulate the partial segment
`[x2, ebb]` into the current group and finalize (`gsig /= sumde`).

The multi-group-spanning case (Fortran goto 1000 loop, lines 5323-5341)
fills intermediate groups via `gsig[ig] = yl*zl + yr*zr/sumde` — note
the Fortran precedence quirk: `yr*zr/sumde` is evaluated as
`(yr*zr)/sumde` first, then added to `yl*zl`. This looks like a Fortran
bug (the obvious intent is `(yl*zl + yr*zr)/sumde`) but we replicate it
verbatim to match canonical NJOY output (CLAUDE.md Rule 1, Law 2).

After the spanning loop, the segment tail `[egn(ig), x1]` is accumulated
into the new group without finalizing (Fortran 5343-5357).

For iwt=2 (`weight_fn(E) = 1`), the Simpson coefficients reduce to
`z1 = z2 = dx/3`, and the per-group integral is exactly the trapezoid
sum — algebraically equivalent to the previous `_group_average_flat`
implementation in the iwt=2 limit. For iwt=6 (T15 default), the
thermal/1/E/fission/fusion weight reshapes the integrand and changes
the group-averaged answer by ~16% on group 1 (capture peak energy
range), which was the root cause of the Phase 72b C[1,1] bias.

`grid` and `sig` must be the same length; `grid` must be sorted and
all entries must fall inside `[egn[1], egn[end]]`.
"""
function _rpxgrp_average(grid::AbstractVector{Float64},
                         sig::AbstractVector{Float64},
                         egn::AbstractVector{<:Real},
                         weight_fn)
    ngn = length(egn) - 1
    npts = length(grid)
    gsig = zeros(Float64, ngn)
    npts < 2 && return gsig

    # Locate the initial group containing grid[1] (Fortran lines
    # 5254-5264 do the same scan; we replace `goto 100 + i0=i0+1` with
    # explicit linear search and treat off-grid points as advancing the
    # scan, matching the Fortran `i0<ipoint` recovery path).
    i0 = 1
    ig = 0
    while i0 <= npts
        ig = _find_output_group(grid[i0], egn, ngn)
        ig > 0 && break
        i0 += 1
    end
    (ig == 0 || i0 >= npts) && return gsig

    sumde = 0.0
    x1 = grid[i0]
    wt1 = Float64(weight_fn(x1))

    @inbounds for i in (i0 + 1):npts
        x2 = x1
        wt2 = wt1
        y2 = sig[i - 1]

        x1 = grid[i]
        wt1 = Float64(weight_fn(x1))
        y1 = sig[i]
        x12 = x1 - x2
        x12 <= 0 && continue

        egnt  = Float64(egn[ig])
        egnt1 = Float64(egn[ig + 1])

        if x1 >= egnt && x1 <= egnt1
            # Case A — segment lies wholly inside the current group.
            # Ref: errorr.f90:5281-5300.
            de = 0.5 * (x1 - x2) * (wt1 + wt2)
            sumde += de
            z1 = (2.0 * wt1 + wt2) * x12 / 6.0
            z2 = (2.0 * wt2 + wt1) * x12 / 6.0
            if y1 != 0.0 || y2 != 0.0
                gsig[ig] += y1 * z1 + y2 * z2
            end
            if x1 == egnt1
                gsig[ig] = sumde != 0.0 ? gsig[ig] / sumde : 0.0
                ig += 1
                sumde = 0.0
                ig > ngn && break
            end
        elseif x1 > egnt1
            # Case B — segment crosses one or more group boundaries.
            # Ref: errorr.f90:5301-5358.
            ebb  = egnt1
            ebx  = ebb - x2
            wt12 = wt1 - wt2
            coef = ebx / x12
            wt3  = wt2 + wt12 * ebx / x12
            de   = ebx * (wt2 + wt3) / 2.0
            sumde += de
            z2 = (2.0 * wt2 + wt3) * ebx / 6.0
            z3 = (2.0 * wt3 + wt2) * ebx / 6.0
            y3 = y2 + (y1 - y2) * coef
            if y2 != 0.0 || y3 != 0.0
                gsig[ig] += y2 * z2 + y3 * z3
            end
            gsig[ig] = sumde != 0.0 ? gsig[ig] / sumde : 0.0

            # Walk through any wholly-spanned intermediate groups
            # (Fortran "1000 continue" goto loop, lines 5323-5341).
            while true
                ig += 1
                ig > ngn && break
                if x1 > Float64(egn[ig + 1])
                    xl = Float64(egn[ig])
                    xr = Float64(egn[ig + 1])
                    ebx_in = xr - xl
                    wtl = wt2 + wt12 * (xl - x2) / x12
                    wtr = wt2 + wt12 * (xr - x2) / x12
                    sumde = ebx_in * (wtl + wtr) / 2.0
                    zl = (2.0 * wtl + wtr) * ebx_in / 6.0
                    zr = (2.0 * wtr + wtl) * ebx_in / 6.0
                    yl = y2 + (y1 - y2) * (xl - x2) / (x1 - x2)
                    yr = y2 + (y1 - y2) * (xr - x2) / (x1 - x2)
                    if yl != 0.0 || yr != 0.0
                        # NOTE: Fortran precedence quirk (errorr.f90:5338):
                        #   gsig(j,ig) = yl*zl + yr*zr/sumde
                        # Parsed as `yl*zl + (yr*zr)/sumde`. Replicated
                        # verbatim to match canonical NJOY output;
                        # do NOT "fix" to `(yl*zl + yr*zr)/sumde`.
                        gsig[ig] = yl * zl + yr * zr / sumde
                    end
                    # Loop continues — try next spanned group.
                else
                    # Spanning ended in this group — accumulate the
                    # tail `[egn(ig), x1]` into it without finalizing.
                    # Ref: errorr.f90:5343-5357.
                    ebx_in = x1 - Float64(egn[ig])
                    de_tail = ebx_in * (wt3 + wt1) / 2.0
                    sumde = de_tail
                    z1 = (2.0 * wt3 + wt1) * ebx_in / 6.0
                    z3 = (2.0 * wt1 + wt3) * ebx_in / 6.0
                    y3_tail = y2 + (y1 - y2) * coef
                    if y3_tail != 0.0 || y1 != 0.0
                        gsig[ig] += y1 * z1 + y3_tail * z3
                    end
                    break
                end
            end
        end
    end

    # Finalize any group still mid-accumulation. (Fortran finalizes
    # only at exact boundary crossings; on the typical pointwise grid,
    # the last segment lands exactly on the topmost egn boundary, so
    # nothing remains. We add a guard for the rare case where the grid
    # ends mid-group, mirroring the implicit "leave gsig as ∫y·w" with
    # a `/sumde` normalize so the result is a proper average.)
    @inbounds if ig >= 1 && ig <= ngn && sumde != 0.0 && gsig[ig] != 0.0
        # Detect whether this group was already finalized: the Fortran
        # path only normalises on exact boundary crossing. If the last
        # x1 didn't equal egn(ig+1), `gsig[ig]` is still pre-normalised.
        # Heuristic: if we got here without an `x1 == egnt1` finalize
        # for the current `ig`, normalise now. This guards against the
        # "grid ends mid-group" edge case rpxgrp's input layout makes
        # impossible (its sig array is built to end on a group boundary).
        last_eb = Float64(egn[ig + 1])
        if grid[end] < last_eb
            gsig[ig] /= sumde
        end
    end

    return gsig
end

# Output-group index (1..ngn) containing energy `e`; 0 if outside.
function _find_output_group(e::Real, egn::AbstractVector{<:Real}, ngn::Int)
    @inbounds for ig in 1:ngn
        if Float64(egn[ig]) <= e < Float64(egn[ig + 1])
            return ig
        end
    end
    return 0
end

# -------------------------------------------------------------------------
# RM resonance-parameter cloning + perturbation
# -------------------------------------------------------------------------

# Match an MF=32 resonance (Er32, AJ32) to an MF=2 ReichMooreParameters
# (l_idx, n_idx). Mirrors errorr.f90:4413-4443 (rpxlc12 search loop).
# Tolerances per the Fortran: |Er/Er2 - 1| < 1e-6 and |AJ - AJ2| < 1e-4.
# Returns (l_idx, n_idx) or (0, 0) if no match.
function _match_mf32_to_mf2_rm(Er32::Float64, AJ32::Float64,
                                rm::ReichMooreParameters)
    Er32_abs = abs(Er32)
    AJ32_abs = abs(AJ32)
    @inbounds for il in 1:Int(rm.NLS)
        for n in eachindex(rm.Er[il])
            Er2 = rm.Er[il][n]
            AJ2 = rm.AJ[il][n]
            if Er32 * Er2 > 0
                rr = abs(Er32 / Er2 - 1)
                rr2 = abs(AJ32_abs - abs(AJ2))
                if rr < 1e-6 && rr2 < 1e-4
                    return (il, n)
                end
            end
        end
    end
    return (0, 0)
end

# Build a perturbed copy of `rm` where parameter `loopn` of resonance
# (l_idx, n_idx) is multiplied by `factor`. `loopn` is 1-indexed:
#   loopn=1 → ER  (errorr.f90:4448-4449 / 4480)
#   loopn=2 → GN  (errorr.f90:4459-4467, lrf=3 branch: il3=il2+loopn=il2+2)
#   loopn=3 → GG  (lrf=3: il3=il2+3)
#   loopn=4 → GFA (lrf=3: il3=il2+4)
#   loopn=5 → GFB (lrf=3: il3=il2+5)
function _perturb_rm(rm::ReichMooreParameters, l_idx::Int, n_idx::Int,
                     loopn::Int, factor::Float64)
    # Deep-copy the affected l-state vectors only.
    Er  = [copy(v) for v in rm.Er]
    Gn  = [copy(v) for v in rm.Gn]
    Gg  = [copy(v) for v in rm.Gg]
    Gfa = [copy(v) for v in rm.Gfa]
    Gfb = [copy(v) for v in rm.Gfb]

    if loopn == 1
        Er[l_idx][n_idx] *= factor
    elseif loopn == 2
        Gn[l_idx][n_idx] *= factor
    elseif loopn == 3
        Gg[l_idx][n_idx] *= factor
    elseif loopn == 4
        Gfa[l_idx][n_idx] *= factor
    elseif loopn == 5
        Gfb[l_idx][n_idx] *= factor
    else
        error("_perturb_rm: invalid loopn=$loopn (must be 1..5 for LRF=3)")
    end

    return ReichMooreParameters(rm.NLS, rm.SPI, rm.AP, rm.LAD,
                                rm.l_values, rm.AWRI, rm.APL,
                                Er, rm.AJ, Gn, Gg, Gfa, Gfb)
end

# Get the unperturbed value of parameter `loopn` for resonance (l, n).
# Used to compute `gwidth` for the central-difference denominator.
function _rm_param_value(rm::ReichMooreParameters, l_idx::Int, n_idx::Int,
                          loopn::Int)
    loopn == 1 && return rm.Er[l_idx][n_idx]
    loopn == 2 && return rm.Gn[l_idx][n_idx]
    loopn == 3 && return rm.Gg[l_idx][n_idx]
    loopn == 4 && return rm.Gfa[l_idx][n_idx]
    loopn == 5 && return rm.Gfb[l_idx][n_idx]
    error("_rm_param_value: invalid loopn=$loopn")
end

# Perturbation factors (errorr.f90:4448-4488, 4500-4505):
#   ER     : ±0.0001 (factor 0.9999 / 1.0001), gwidth = ER · 1e-4
#   widths : ±0.01   (factor 0.99   / 1.01  ), gwidth = w  · 1e-2
# Convention matches Fortran rpxlc12 verbatim — these are the
# perturbation magnitudes used as central-difference half-step.
const _RESCON_PERT_ER    = 1.0e-4
const _RESCON_PERT_WIDTH = 1.0e-2

function _rm_pert_factors(loopn::Int)
    p = loopn == 1 ? _RESCON_PERT_ER : _RESCON_PERT_WIDTH
    return (1.0 - p, 1.0 + p)
end

function _rm_gwidth(unperturbed::Float64, loopn::Int)
    p = loopn == 1 ? _RESCON_PERT_ER : _RESCON_PERT_WIDTH
    return unperturbed * p
end

# -------------------------------------------------------------------------
# Full evaluation: σ(E) on the grid for given RM parameters
# -------------------------------------------------------------------------

# Evaluate (total, elastic, fission, capture) at each E in `grid` for
# the given Reich-Moore parameter set + range. Returns 4 vectors of
# length(grid). Reuses the existing `cross_section_rm` evaluator.
function _eval_xs_grid_rm(grid::Vector{Float64},
                          rm::ReichMooreParameters,
                          range::ResonanceRange)
    n = length(grid)
    σtot = Vector{Float64}(undef, n)
    σel  = Vector{Float64}(undef, n)
    σfis = Vector{Float64}(undef, n)
    σcap = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        cs = cross_section_rm(grid[i], rm, range)
        σtot[i] = cs.total
        σel[i]  = cs.elastic
        σfis[i] = cs.fission
        σcap[i] = cs.capture
    end
    return σtot, σel, σfis, σcap
end

# -------------------------------------------------------------------------
# Subsection-level sensitivity build
# -------------------------------------------------------------------------

# For one MF=32 subsection (one NSRS block within one resolved range),
# build sens[channel, p, ig] via central-difference perturbation +
# group-averaging. Returns Array{Float64,3}(4, npar, ngn).
#
# Mirrors errorr.f90:4399-4523 (rpxlc12 inner perturbation loop):
#   for loopm = 1..nrb (resonance index)
#     for loopn = 1..mpar (parameter index)
#       ... central-difference, group-average ...
#
# `mf2_range` must be the MF=2 ResonanceRange whose RM parameters cover
# [sub.elr, sub.ehr] — found by overlapping ranges in MF=2 vs MF=32.
function _build_subsection_sensitivities_rm(
    sub::MF32ResolvedSubsection,
    mf32_range::MF32ResolvedRange,
    mf2_range::ResonanceRange{ReichMooreParameters},
    egn::AbstractVector{<:Real},
    ngn::Int,
    weight_fn,
)
    rm = mf2_range.parameters
    npar = sub.mpar * sub.nrb

    # Estimate width per resonance from MF=2 (Gn + Gg) for grid sizing.
    # Use absolute MF=32 ER as proxy when matching fails (rare).
    Es = Float64[sub.params[r, 1] for r in 1:sub.nrb]
    Ws = zeros(Float64, sub.nrb)
    for r in 1:sub.nrb
        Er32 = sub.params[r, 1]
        AJ32 = sub.params[r, 2]
        (l, n) = _match_mf32_to_mf2_rm(Er32, AJ32, rm)
        if l != 0
            Ws[r] = abs(rm.Gn[l][n]) + abs(rm.Gg[l][n])
        else
            Ws[r] = 0.025  # generic default ~25 meV
        end
    end

    grid = _build_pointwise_grid(mf32_range.elr, mf32_range.ehr, Es, Ws, egn)

    # σ̄_unperturbed per channel (for sandwich denominator and as the
    # central-difference reference in case we ever switch to forward).
    σtot0, σel0, σfis0, σcap0 = _eval_xs_grid_rm(grid, rm, mf2_range)
    σ̄_tot0 = _rpxgrp_average(grid, σtot0, egn, weight_fn)
    σ̄_el0  = _rpxgrp_average(grid, σel0,  egn, weight_fn)
    σ̄_fis0 = _rpxgrp_average(grid, σfis0, egn, weight_fn)
    σ̄_cap0 = _rpxgrp_average(grid, σcap0, egn, weight_fn)

    sens = zeros(Float64, 4, npar, ngn)

    p = 0
    for loopm in 1:sub.nrb
        Er32 = sub.params[loopm, 1]
        AJ32 = sub.params[loopm, 2]
        (l, n) = _match_mf32_to_mf2_rm(Er32, AJ32, rm)
        if l == 0
            # Unmatched MF=32 resonance — Fortran calls error(); we
            # warn and zero-fill its sensitivity entries (better than
            # crashing on slightly-mismatched MF=2/MF=32 pairs).
            @warn "rescon: MF=32 resonance #$loopm (ER=$Er32, AJ=$AJ32) \
                  has no matching MF=2 RP — zeroing sens entries."
            p += sub.mpar
            continue
        end

        for loopn in 1:sub.mpar
            p += 1
            unperturbed = _rm_param_value(rm, l, n, loopn)
            gwidth = _rm_gwidth(unperturbed, loopn)
            if gwidth == 0.0
                # Width is identically zero (e.g. fission = 0 for U-238)
                # — Fortran's `if (gwidth.ne.zero)` zeroes sens here too
                # (errorr.f90:4515-4520).
                continue
            end
            (fminus, fplus) = _rm_pert_factors(loopn)

            rm_minus = _perturb_rm(rm, l, n, loopn, fminus)
            rm_plus  = _perturb_rm(rm, l, n, loopn, fplus)

            σm = _eval_xs_grid_rm(grid, rm_minus, mf2_range)
            σp = _eval_xs_grid_rm(grid, rm_plus,  mf2_range)

            # Central difference per-point: (σ_+ - σ_-) / (2·gwidth).
            inv2g = 1.0 / (2.0 * gwidth)
            dtot = @. (σp[1] - σm[1]) * inv2g
            del  = @. (σp[2] - σm[2]) * inv2g
            dfis = @. (σp[3] - σm[3]) * inv2g
            dcap = @. (σp[4] - σm[4]) * inv2g

            ḡ_tot = _rpxgrp_average(grid, dtot, egn, weight_fn)
            ḡ_el  = _rpxgrp_average(grid, del,  egn, weight_fn)
            ḡ_fis = _rpxgrp_average(grid, dfis, egn, weight_fn)
            ḡ_cap = _rpxgrp_average(grid, dcap, egn, weight_fn)

            @inbounds for ig in 1:ngn
                sens[1, p, ig] = ḡ_tot[ig]
                sens[2, p, ig] = ḡ_el[ig]
                sens[3, p, ig] = ḡ_fis[ig]
                sens[4, p, ig] = ḡ_cap[ig]
            end
        end
    end

    # Trim sens to the resolved-range group window [iest, ieed].
    # Ref: errorr.f90:3093-3108 (iest/ieed derivation in resprx) + 4509
    # (sensitivity is only WRITTEN inside [iest, ieed]; outside, sens
    # stays at its initial zero from `sens = zeros(...)`). Without this,
    # _rpxgrp_average returns small but non-zero values for groups
    # outside the resolved range (FP noise from boundary handling),
    # which the sandwich amplifies across npar² parameters into
    # spurious rel-cov entries that survive the 1e-20 row-trim
    # threshold. The MF=33 self-cov rows for groups well above the
    # resolved range then bloat to span ~all output groups.
    iest, ieed = _resonance_group_window(mf32_range.elr, mf32_range.ehr, egn)
    if iest == 0 || ieed == 0 || iest > ieed
        fill!(sens, 0.0)
    else
        @inbounds for ig in 1:ngn
            (iest <= ig <= ieed) && continue
            for ch in 1:4, p in 1:npar
                sens[ch, p, ig] = 0.0
            end
        end
    end

    return sens, σ̄_tot0, σ̄_el0, σ̄_fis0, σ̄_cap0
end

"""
    _resonance_group_window(elr, ehr, egn) -> (iest::Int, ieed::Int)

Find the output-group window [iest, ieed] covering the resolved range
[elr, ehr], matching Fortran's iest/ieed derivation at errorr.f90:3093-3108.

- iest = first group whose upper bound exceeds elr
- ieed = last group whose upper bound is ≥ ehr
- if elg (lower bound of iest's group) ≥ ehr the range collapses to ieed=0

Returns (0, 0) if no group satisfies the criteria.
"""
function _resonance_group_window(elr::Real, ehr::Real,
                                 egn::AbstractVector{<:Real})
    ngn = length(egn) - 1
    iest = 0; ieed = 0
    elg = Float64(elr); ehg = Float64(ehr)
    ip1 = false; ip2 = false
    @inbounds for i in 2:(ngn+1)
        ee = Float64(egn[i])
        if !ip1 && ee > elr
            ip1 = true; elg = Float64(egn[i-1]); iest = i - 1
        end
        if !ip2 && ee >= ehr
            ip2 = true; ehg = ee; ieed = i - 1
        end
    end
    if elg >= ehr
        ieed = 0
    end
    return (iest, ieed)
end

# -------------------------------------------------------------------------
# Sandwich (errorr.f90:4541-4593, plus the σ̄·σ̄ relative conversion)
# -------------------------------------------------------------------------

# Apply the J·cov_RP·J^T sandwich for one (mt, mt2) pair to
# `cov_matrices`. Produces RELATIVE cov by dividing by σ̄·σ̄.
#
# Fortran rpxlc12 stores absolute cov in cff/cgg/.../ctt at lines
# 4550-4593. The σ̄·σ̄ relative conversion happens in covout (writer
# side) before it lands in MF=33. Julia's cov_matrices already holds
# RELATIVE cov post-Phase-51, so we do the conversion here.
#
# `sens_ch` is the channel index (1..4) per the rescon dispatch.
function _apply_sandwich!(
    cov_matrices::Dict{Tuple{Int,Int},Matrix{Float64}},
    mt::Int, mt2::Int,
    sens::Array{Float64,3},
    cov_RP::Matrix{Float64},
    σ̄_row::Vector{Float64}, σ̄_col::Vector{Float64},
    ngn::Int,
)
    (ch_row, ch_col) = _rescon_channel_pair(mt, mt2)
    ch_row == 0 && return  # not a rescon-supported pair

    npar = size(cov_RP, 1)
    canonical_mt, canonical_mt2 = mt <= mt2 ? (mt, mt2) : (mt2, mt)

    # Build absolute-cov contribution into a fresh ngn×ngn buffer.
    abs_cov = zeros(Float64, ngn, ngn)

    if canonical_mt == canonical_mt2
        # Self-pair: triangular walk + symmetric fill (Fortran 4541-4566).
        @inbounds for ig in 1:ngn
            for ig2 in ig:ngn
                acc = 0.0
                for i in 1:npar
                    si_row = sens[ch_row, i, ig]
                    si_col = sens[ch_col, i, ig2]
                    for j in i:npar
                        c = cov_RP[i, j]
                        c == 0.0 && continue
                        sj_row = sens[ch_row, j, ig]
                        sj_col = sens[ch_col, j, ig2]
                        acc += c * si_row * sj_col
                        if i != j
                            acc += c * sj_row * si_col
                        end
                    end
                end
                abs_cov[ig, ig2] = acc
                if ig != ig2
                    abs_cov[ig2, ig] = acc
                end
            end
        end
    else
        # Cross-pair: full ngn×ngn (Fortran 4571-4593). Note the
        # asymmetry — channel is row=mt, col=mt2 in canonical order.
        @inbounds for ig in 1:ngn
            for ig2 in 1:ngn
                acc = 0.0
                for i in 1:npar
                    si_row = sens[ch_row, i, ig]
                    si_col = sens[ch_col, i, ig2]
                    for j in i:npar
                        c = cov_RP[i, j]
                        c == 0.0 && continue
                        sj_row = sens[ch_row, j, ig]
                        sj_col = sens[ch_col, j, ig2]
                        acc += c * si_row * sj_col
                        if i != j
                            acc += c * sj_row * si_col
                        end
                    end
                end
                abs_cov[ig, ig2] = acc
            end
        end
    end

    # Divide by σ̄(g)·σ̄(g') to get RELATIVE cov. Skip cells where
    # either σ̄ is zero (sub-threshold groups stay zero).
    rel_cov = zeros(Float64, ngn, ngn)
    @inbounds for ig in 1:ngn, ig2 in 1:ngn
        denom = σ̄_row[ig] * σ̄_col[ig2]
        rel_cov[ig, ig2] = denom > 0 ? abs_cov[ig, ig2] / denom : 0.0
    end

    key = (canonical_mt, canonical_mt2)
    if haskey(cov_matrices, key)
        cov_matrices[key] .+= rel_cov
    else
        cov_matrices[key] = rel_cov
    end
    return nothing
end

# =========================================================================
# URR rescon (LRU=2) — Phase 75
#
# Mirrors Fortran rpxunr (errorr.f90:4785-4985) + ggunr1
# (errorr.f90:6800-6905). For each MF=32 LRU=2 range:
#   1. Flatten per-J params (D, AJ, GNO, GG, GF, GX) from MF=32.
#   2. Pull DOFs (mu, nu, lamda) per J-state from MF=2 (URR2Data via
#      rdumrd2's amur table, errorr.f90:5196-5213).
#   3. Build a dense URR XS grid via _build_urr_xs_grid (E-step ×1.015
#      with elr/ehr/ehg breakpoints — errorr.f90:4904-4930).
#   4. Group-average via _rpxgrp_average to get gsig_base[ch][ig].
#   5. For each (J, param-within-J) — npar = mpar × ΣNJS total —
#      perturb that param by 1.01×, recompute group-avg, build
#      sens = 100·(gsig_p - gsig_0)/orig (forward diff). Trim to
#      [iest, ieed] per errorr.f90:3093-3108.
#   6. Convert MF=32 URR cov from relative-relative to absolute by
#      multiplying by b_i·b_j (errorr.f90:4957-4964: cov(i,j)=a(l3)·bb
#      with bb=b(l2+i)·b(l2+j)).
#   7. Reuse _apply_sandwich! for each of the seven RESCON_PAIRS.
# =========================================================================

# MPAR-to-parameter-offset table (offsets index `params_per_j[k]` which
# stores [D, AJ, GNO, GG, GF, GX]). Mirrors Fortran rpxunr perturbation
# offsets at errorr.f90:4865-4897 — b(inow+1) D, b(inow+3) GNO,
# b(inow+4) GG, then either b(inow+5) GF (lfw=1) or b(inow+6) GX (lfw=0)
# for mpar=4, and finally b(inow+6) GX for mpar=5.
function _urr_param_offsets(mpar::Int, lfw::Int)
    mpar == 1 && return (1,)
    mpar == 2 && return (1, 3)
    mpar == 3 && return (1, 3, 4)          # U-238 JENDL
    mpar == 4 && return lfw == 1 ? (1, 3, 4, 5) : (1, 3, 4, 6)
    mpar == 5 && return (1, 3, 4, 5, 6)
    error("URR rescon: unsupported MPAR=$mpar (Fortran rpxunr supports 1..5).")
end

# Flatten MF=32 LRU=2 L-states into per-J vectors of mutable params.
# Returns (params, l_per_j, awri_per_j). params[k] is a 6-vector
# (D, AJ, GNO, GG, GF, GX) that can be mutated during the perturbation
# loop.
function _flatten_urr_params(urr::MF32UnresolvedRange)
    params = Vector{Float64}[]
    l_per_j = Int[]
    awri_per_j = Float64[]
    for ls in urr.l_states
        for js in ls.j_states
            push!(params, [js.D, js.AJ, js.GNO, js.GG, js.GF, js.GX])
            push!(l_per_j, ls.l)
            push!(awri_per_j, ls.awri)
        end
    end
    return params, l_per_j, awri_per_j
end

# DOFs (mu, nu, lamda) per J-state in (L, J) order — mirrors Fortran
# rdumrd2's `amur(1:3, nlru2) = (AMUN, AMUF, AMUX)` at errorr.f90:5206-5208.
function _urr_dofs_from_mf2(rng::ResonanceRange)
    p = rng.parameters
    if p isa URR2Data
        return [(round(Int, s.AMUN), round(Int, s.AMUF), round(Int, s.AMUX))
                for s in p.sequences]
    elseif p isa URRData
        return [(round(Int, s.AMUN), s.MUF, 0) for s in p.sequences]
    else
        error("_urr_dofs_from_mf2: unsupported URR formalism \
              $(typeof(p)) (need URR2Data or URRData).")
    end
end

# (L, J) tuple per sequence in MF=2 URR walk order. Used to verify
# that MF=2 sequences[] and MF=32 L-states/J-states list the same
# (L, J) pairs in the same order — the dofs/sandwich indexing assumes
# this 1:1 alignment.
function _urr_lj_pairs_from_mf2(rng::ResonanceRange)
    p = rng.parameters
    if p isa URR2Data || p isa URRData
        return [(s.l, Float64(s.J)) for s in p.sequences]
    else
        error("_urr_lj_pairs_from_mf2: unsupported URR formalism \
              $(typeof(p)).")
    end
end

# Evaluate URR XS at one E. Mirrors Fortran ggunr1 (errorr.f90:6800-6905).
# Returns (sig_tot, sig_el, sig_fis, sig_cap) in barns.
function _ggunr1(E::Float64,
                 params::Vector{Vector{Float64}},
                 l_per_j::Vector{Int},
                 awri_per_j::Vector{Float64},
                 dofs::Vector{NTuple{3,Int}},
                 spi::Float64, ap::Float64)
    Cp = PhysicsConstants
    cwaven = sqrt(2 * Cp.amassn * Cp.amu * Cp.ev) * 1e-12 / Cp.hbar
    sig_e = 0.0; sig_f = 0.0; sig_c = 0.0; spot = 0.0
    prev_l = -1
    @inbounds for k in eachindex(params)
        p = params[k]
        D, AJ, GNO, GG, GF, GX = p[1], p[2], p[3], p[4], p[5], p[6]
        ll  = l_per_j[k]
        awri = awri_per_j[k]
        rat  = awri / (awri + 1)
        aw   = awri * Cp.amassn
        ra   = _RC1 * aw^_THIRD_URR + _RC2
        cnst = (2 * Cp.pi^2) / (cwaven * rat)^2
        mu, nu, lamda = dofs[k]
        gj = (2 * AJ + 1) / (4 * spi + 2)
        e2v = sqrt(E)
        kwv = rat * e2v * cwaven
        rho = kwv * ra
        rhoc = kwv * ap
        vl, ps = _unfac(ll, rho, rhoc, Float64(mu))
        vl *= e2v
        if ll != prev_l
            spot += 4 * Cp.pi * (2 * ll + 1) * (sin(ps) / kwv)^2
            prev_l = ll
        end
        gnx = GNO * vl
        diff = GX
        den = E * D
        temp = cnst * gj * gnx / den
        terg = temp * GG
        ters = temp * gnx
        terf = temp * GF
        gs  = _gnrl(gnx, GF, GG, mu, nu, lamda, diff, 1)
        gc  = _gnrl(gnx, GF, GG, mu, nu, lamda, diff, 2)
        gff = _gnrl(gnx, GF, GG, mu, nu, lamda, diff, 3)
        gc  *= terg
        gff *= terf
        gs  *= ters
        # Interference correction (Fortran ggunr1 lines 6887-6890).
        add = cnst * gj * 2 * gnx * sin(ps)^2 / den
        gs -= add
        sig_e += gs
        sig_f += gff
        sig_c += gc
    end
    sig_e += spot
    return (sig_e + sig_f + sig_c, sig_e, sig_f, sig_c)
end

# Dense URR XS grid for one MF=32 URR range. Mirrors the E-stepping in
# Fortran rpxunr (errorr.f90:4904-4930): start at elg, step ×1.015 with
# forced breakpoints at elr / ehr / ehg. Outside [elr, ehr], XS = 0.
function _build_urr_xs_grid(urr::MF32UnresolvedRange,
                            params::Vector{Vector{Float64}},
                            l_per_j::Vector{Int},
                            awri_per_j::Vector{Float64},
                            dofs::Vector{NTuple{3,Int}},
                            elg::Float64, ehg::Float64)
    elr = urr.elr; ehr = urr.ehr
    spi = urr.spi; ap = urr.ap
    e_pts = Float64[]
    σtot = Float64[]; σel = Float64[]; σfis = Float64[]; σcap = Float64[]
    e1 = elg
    while true
        e1 > ehg && (e1 = ehg)
        if elr <= e1 <= ehr
            (st, se, sf, sc) = _ggunr1(e1, params, l_per_j, awri_per_j,
                                        dofs, spi, ap)
        else
            st = se = sf = sc = 0.0
        end
        push!(e_pts, e1); push!(σtot, st); push!(σel, se)
        push!(σfis, sf); push!(σcap, sc)
        e1 >= ehg && break
        e1_new = e1 * 1.015
        ebc = e1_new / 1.015
        ebc < elr && e1_new > elr && (e1_new = elr)
        ebc < ehr && e1_new > ehr && (e1_new = ehr)
        e1 = e1_new
    end
    return e_pts, σtot, σel, σfis, σcap
end

# Sensitivity build + sandwich for one MF=32 URR range. Returns the
# number of (mt, mt2) pair·subsection contributions applied.
function _apply_rescon_urr_range!(
    cov_matrices::Dict{Tuple{Int,Int},Matrix{Float64}},
    urr::MF32UnresolvedRange,
    mf2_urr_range::ResonanceRange,
    egn::AbstractVector{<:Real},
    group_xs::Dict{Int,Vector{Float64}},
    weight_fn,
)
    ngn = length(egn) - 1
    iest, ieed = _resonance_group_window(urr.elr, urr.ehr, egn)
    (iest == 0 || ieed == 0 || iest > ieed) && return 0
    elg = Float64(egn[iest])
    ehg = Float64(egn[ieed + 1])

    params, l_per_j, awri_per_j = _flatten_urr_params(urr)
    nJ = length(params)
    mpar = urr.mpar
    npar = mpar * nJ
    npar == size(urr.cov_RP, 1) || error(
        "rescon URR: NPAR=$npar (= MPAR·nJ = $mpar·$nJ) ≠ cov_RP size \
         $(size(urr.cov_RP)).")

    dofs = _urr_dofs_from_mf2(mf2_urr_range)
    length(dofs) == nJ || error(
        "rescon URR: MF=32 has $nJ J-states but MF=2 URR has \
         $(length(dofs)) — (L, J) walk count mismatch.")

    # Defensive (L, J) tuple alignment check — Phase 75 follow-up.
    # The dofs and sens arrays index J-states in MF=32's L-major walk
    # order; the sandwich relies on MF=2's `sequences[k]` referring to
    # the same (L, J) as MF=32's k-th flattened J-state. ENDF-6 §32.2
    # does not formally guarantee this — verified empirically on
    # U-238 JENDL. Error loudly if a future tape violates the assumption
    # so we don't silently apply mismatched DOFs.
    mf2_lj = _urr_lj_pairs_from_mf2(mf2_urr_range)
    for k in 1:nJ
        mf32_l = l_per_j[k]
        mf32_j = params[k][2]
        (mf2_l, mf2_j) = mf2_lj[k]
        if mf32_l != mf2_l || !isapprox(mf32_j, mf2_j; atol=1e-4)
            error("rescon URR: (L, J) mismatch at sequence $k — \
                   MF=32 has (L=$mf32_l, J=$mf32_j) but MF=2 has \
                   (L=$mf2_l, J=$mf2_j). The 1:1 walk assumption is \
                   violated; rescon would apply wrong DOFs/sens. \
                   Need a (L, J) lookup-based pairing instead.")
        end
    end

    offsets = _urr_param_offsets(mpar, Int(mf2_urr_range.LFW))

    # Base XS over the dense URR grid.
    e_pts, σt0, σe0, σf0, σc0 = _build_urr_xs_grid(
        urr, params, l_per_j, awri_per_j, dofs, elg, ehg)
    g_tot0 = _rpxgrp_average(e_pts, σt0, egn, weight_fn)
    g_el0  = _rpxgrp_average(e_pts, σe0, egn, weight_fn)
    g_fis0 = _rpxgrp_average(e_pts, σf0, egn, weight_fn)
    g_cap0 = _rpxgrp_average(e_pts, σc0, egn, weight_fn)

    sens = zeros(Float64, 4, npar, ngn)
    sigma_orig = zeros(Float64, npar)

    # Walk (J, param-in-J) in J-major order — matches cov_RP layout per
    # Fortran rpxunr (each J advances by mpar params before moving on).
    p_idx = 0
    for k in 1:nJ
        for q in 1:mpar
            p_idx += 1
            off = offsets[q]
            orig = params[k][off]
            sigma_orig[p_idx] = orig
            # Skip identically-zero params (e.g. GF=0 for U-238
            # sub-threshold fission) — Fortran skips via the
            # `if (b(l0+1+i).ne.zero)` guard at errorr.f90:4942.
            orig == 0.0 && continue

            params[k][off] = orig * 1.01
            _, σt, σe, σf, σc = _build_urr_xs_grid(
                urr, params, l_per_j, awri_per_j, dofs, elg, ehg)
            g_tot = _rpxgrp_average(e_pts, σt, egn, weight_fn)
            g_el  = _rpxgrp_average(e_pts, σe, egn, weight_fn)
            g_fis = _rpxgrp_average(e_pts, σf, egn, weight_fn)
            g_cap = _rpxgrp_average(e_pts, σc, egn, weight_fn)
            params[k][off] = orig

            # sens = (Δσ̄)/(0.01·orig) = forward-diff dσ̄/dp. Match
            # Fortran rpxunr's `s = 100*s/b(l0+1+i)` at errorr.f90:4944
            # (the `100·s/orig` form is algebraically identical to
            # `s/(0.01·orig)`). The 1e-10 guard at line 4945 zeroes
            # sub-threshold sensitivity components.
            inv = 100.0 / orig
            @inbounds for ig in iest:ieed
                s1 = (g_tot[ig] - g_tot0[ig]) * inv
                s2 = (g_el[ig]  - g_el0[ig])  * inv
                s3 = (g_fis[ig] - g_fis0[ig]) * inv
                s4 = (g_cap[ig] - g_cap0[ig]) * inv
                abs(s1) >= 1e-10 && (sens[1, p_idx, ig] = s1)
                abs(s2) >= 1e-10 && (sens[2, p_idx, ig] = s2)
                abs(s3) >= 1e-10 && (sens[3, p_idx, ig] = s3)
                abs(s4) >= 1e-10 && (sens[4, p_idx, ig] = s4)
            end
        end
    end

    # Convert MF=32 URR cov from relative-relative (cov_REL[i,j]) to
    # absolute (cov_ABS[i,j] = cov_REL · b_i · b_j) — Fortran rpxunr
    # lines 4957-4964 (`bb=b(l2+i)*b(l2+j); cov(i,j)=a(l3)*bb`).
    cov_abs = Matrix{Float64}(undef, npar, npar)
    @inbounds for i in 1:npar, j in 1:npar
        cov_abs[i, j] = urr.cov_RP[i, j] * sigma_orig[i] * sigma_orig[j]
    end

    # Sandwich for each RESCON_PAIR. For the σ̄·σ̄ relative-cov
    # denominator, prefer the GENDF-averaged σ from group_xs; fall back
    # to the URR-local σ̄ when group_xs lacks the channel-MT.
    σ̄_local = (g_tot0, g_el0, g_fis0, g_cap0)
    ch_to_mt = (1, 2, 18, 102)
    σ̄_by_ch = ntuple(ch -> begin
        mt_for_ch = ch_to_mt[ch]
        haskey(group_xs, mt_for_ch) ? group_xs[mt_for_ch] : σ̄_local[ch]
    end, 4)

    n_applied = 0
    for (mt, mt2) in RESCON_PAIRS
        (ch_row, ch_col) = _rescon_channel_pair(mt, mt2)
        σ̄_row = σ̄_by_ch[ch_row]
        σ̄_col = σ̄_by_ch[ch_col]
        _apply_sandwich!(cov_matrices, mt, mt2,
                         sens, cov_abs, σ̄_row, σ̄_col, ngn)
        n_applied += 1
    end
    return n_applied
end

# -------------------------------------------------------------------------
# Top-level: apply_rescon!
# -------------------------------------------------------------------------

"""
    apply_rescon!(cov_matrices, endf_path, mat, egn, group_xs; iwt=2)

Add MF=32 (resonance-parameter) covariance contributions to the seven
relevant `(mt, mt2)` entries of `cov_matrices`. Mirrors Fortran covout
(errorr.f90:7465) → rescon (errorr.f90:8513-8819) + the rpxlc12
sensitivity build (errorr.f90:4399-4593).

Reads MF=32 + MF=2 from `endf_path`, builds finite-difference
sensitivity Jacobians by perturbing each resonance parameter, and
applies the sandwich to update `cov_matrices` in place.

The `iwt` keyword threads the user-deck weight function (errorr card 2
`iwt`) through the group-averaging integrator (`_rpxgrp_average`),
mirroring Fortran `egtwtf` (errorr.f90:10023-10110). For T15 U-238
JENDL, `iwt=6` (thermal + 1/E + fission + fusion) is the correct
weight; Phase 72b mistakenly used `iwt=2` (flat), biasing the
sensitivities by ~16% on group 1.

Returns the parsed `MF32Data` (or `nothing` if MF=32 is absent).
"""
function apply_rescon!(
    cov_matrices::Dict{Tuple{Int,Int},Matrix{Float64}},
    endf_path::AbstractString,
    mat::Integer,
    egn::AbstractVector{<:Real},
    group_xs::Dict{Int,Vector{Float64}};
    iwt::Int = 2,
)
    weight_fn = get_weight_function(iwt)
    open(endf_path, "r") do io
        find_section(io, 32, 151; target_mat=Int(mat)) || return nothing
    end

    data = read_mf32(endf_path, mat)

    # Read MF=2 once for cross-resonance lookup.
    mf2 = open(endf_path, "r") do io
        find_section(io, 2, 151; target_mat=Int(mat)) || return nothing
        read_mf2(io)
    end
    mf2 === nothing && error(
        "apply_rescon!: MF=32 present but MF=2 missing for MAT=$mat — \
         RP perturbation requires resonance parameters from MF=2.")

    ngn = length(egn) - 1
    n_pairs_applied = 0

    for iso in data.isotopes
        # Find matching MF=2 isotope.
        mf2_iso_idx = findfirst(i -> isapprox(i.ZAI, iso.zai; atol=1e-3),
                                mf2.isotopes)
        mf2_iso_idx === nothing && continue
        mf2_iso = mf2.isotopes[mf2_iso_idx]

        for rng in iso.resolved_ranges
            # Find matching MF=2 range by [EL, EH] overlap.
            mf2_range_idx = findfirst(r ->
                isapprox(Float64(r.EL), rng.elr; atol=1e-3) &&
                isapprox(Float64(r.EH), rng.ehr; atol=1e-3),
                mf2_iso.ranges)
            if mf2_range_idx === nothing
                @warn "apply_rescon!: MF=32 range [$(rng.elr), $(rng.ehr)] \
                      has no matching MF=2 range — skipping."
                continue
            end
            mf2_range = mf2_iso.ranges[mf2_range_idx]
            if !(mf2_range.parameters isa ReichMooreParameters)
                @warn "apply_rescon!: MF=2 range is not Reich-Moore (LRF=3) \
                      — only RM port is implemented; skipping."
                continue
            end
            rm_range = mf2_range::ResonanceRange{ReichMooreParameters}

            for sub in rng.subsections
                # Build sensitivities for this subsection.
                sens, σ̄_local_tot, σ̄_local_el, σ̄_local_fis, σ̄_local_cap =
                    _build_subsection_sensitivities_rm(sub, rng, rm_range, egn, ngn,
                                                       weight_fn)

                # Apply sandwich for each of the seven pairs. For the
                # σ̄·σ̄ relative-cov denominator, prefer Fortran-equivalent
                # group-averaged σ from GENDF (already computed upstream
                # as `group_xs`) — cov_matrices stores RELATIVE cov in
                # the same convention as group_xs. Fall back to the
                # locally-built σ̄ when group_xs lacks an entry for the
                # channel-MT (e.g. no fission for U-238 → no MT=18).
                #
                # Channel-to-MT mapping: total=MT1, elastic=MT2,
                # fission=MT18, capture=MT102.
                σ̄_local = (σ̄_local_tot, σ̄_local_el,
                            σ̄_local_fis, σ̄_local_cap)
                ch_to_mt = (1, 2, 18, 102)
                σ̄_by_ch = ntuple(ch -> begin
                    mt_for_ch = ch_to_mt[ch]
                    haskey(group_xs, mt_for_ch) ?
                        group_xs[mt_for_ch] : σ̄_local[ch]
                end, 4)

                for (mt, mt2) in RESCON_PAIRS
                    (ch_row, ch_col) = _rescon_channel_pair(mt, mt2)
                    σ̄_row = σ̄_by_ch[ch_row]
                    σ̄_col = σ̄_by_ch[ch_col]
                    _apply_sandwich!(cov_matrices, mt, mt2,
                                     sens, sub.cov, σ̄_row, σ̄_col, ngn)
                    n_pairs_applied += 1
                end
            end
        end

        # URR (LRU=2) rescon — Phase 75 port. Each URR range maps to one
        # MF=2 LRU=2 range (either LRF=1 → URRData or LRF=2 → URR2Data).
        # The MF=32 LRR/LRF labelling differs from MF=2 (MF=32 LRU=2
        # always reports LRF=1 in its range head — see ENDF-6 §32.2),
        # so match purely by [EL, EH].
        for urr in iso.unresolved_ranges
            mf2_range_idx = findfirst(r ->
                Int(r.LRU) == 2 &&
                isapprox(Float64(r.EL), urr.elr; atol=1.0) &&
                isapprox(Float64(r.EH), urr.ehr; atol=1.0),
                mf2_iso.ranges)
            if mf2_range_idx === nothing
                @warn "apply_rescon!: MF=32 URR range \
                      [$(urr.elr), $(urr.ehr)] has no matching MF=2 \
                      URR range — skipping."
                continue
            end
            mf2_urr_range = mf2_iso.ranges[mf2_range_idx]
            n_applied = _apply_rescon_urr_range!(
                cov_matrices, urr, mf2_urr_range, egn, group_xs, weight_fn)
            n_pairs_applied += n_applied
        end
    end

    @info "apply_rescon!: MAT=$mat — applied RP-cov sandwich to \
          $n_pairs_applied (mt, mt2)·subsection contributions."

    return data
end
