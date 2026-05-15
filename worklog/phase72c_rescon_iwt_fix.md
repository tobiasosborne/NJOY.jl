# Phase 72c — rescon iwt=6 weight fix

**Date**: 2026-05-16
**Branch**: master
**Status**: GREEN — T15 MT=102 row-1 canary C[1,1] now within 1e-7 of reference.

## Summary (one line)

Replaced the iwt=2 (flat) trapezoidal group-averager in `rescon.jl` with a
faithful port of Fortran `rpxgrp` (errorr.f90:5227-5366) that consumes the
user-deck weight function via `get_weight_function(iwt)`. T15 MT=102 self-cov
C[1,1] moved from **3.580e-4 → 2.658408e-4** (vs reference 2.658914e-4,
**diff = -5.06e-8, ~0.019% relative**). The `@test_broken` for the 1e-7
canary in `test/validation/test_errorr_covcal_lb5.jl` now flips to `@test`.

## Bug

Phase 72b's `apply_rescon!` was written with the explicit assumption that
T15 used iwt=2 (constant flux weight). The rescon.jl:166-169 comment said
"T15 iwt=2 standard … per-E flux weight is 1.0".

**That was wrong.** Reading `njoy-reference/tests/15/input` card line 34:

```
 9237 3 6 1 1 /        ! mat ign iwt iprint irelco
                ^
                iwt=6 — thermal + 1/E + fission + fusion (Maxwell-like)
```

T15 uses iwt=6. The mismatch caused `_group_average_flat` to integrate
σ(E) with a unit weight instead of the iwt=6 weight `wtf(E)` from
Fortran `egtwtf`. On group 1 of the LANL-30 structure (1e-5 → 0.414 eV,
spanning the Maxwellian peak), this biased σ̄ uniformly by ~16%. The
sensitivity Jacobian inherited the bias; sandwiched against itself,
C[1,1] grew by 1.16² ≈ 1.346× — matching the observed ratio (Julia
3.580e-4 vs reference 2.659e-4).

## Diagnostic chain (from the orchestrator's brief)

1. **σ̄ bisection**: Phase 72b sandwich `C[1,1] / |C[1,9]|` cancellation
   pattern matched Fortran to 0.08% (1.1711 vs 1.1720) — pure scale bias,
   not differential. Ruled out finite-difference / grid-density / search-
   tolerance differences.
2. **Pointwise vs group-averaged bisection**: σ̄ for group 1 from
   `_group_average_flat` was 1.16× too large vs Fortran's `gsig(4,1)`.
3. **Input deck audit**: T15 `errorr` card has iwt=6, not the iwt=2
   Phase 72b assumed.
4. **Fortran ground truth**: `rpxgrp` (errorr.f90:5227-5366) calls
   `egtwtf(x,enext,idis,lord,wt)` per pointwise energy. For iwt=6 the
   weight is the four-region Maxwell/1-E/fission/fusion function
   (errorr.f90:10135-10156).
5. **Julia weight function check**: `get_weight_function(6)` exists in
   `src/processing/weight_functions.jl` (function `thermal_fission_fusion`)
   and matches the Fortran iwt=6 algorithmic spec.

## Fix

### Files changed

- `src/processing/rescon.jl`
  - Replaced `_group_average_flat` (split-segment trapezoidal, unit weight)
    with `_rpxgrp_average(grid, sig, egn, weight_fn)` — a line-by-line port
    of Fortran `rpxgrp` (errorr.f90:5227-5366). Per-segment Simpson
    coefficients `z1=(2·wt1+wt2)·dx/6`, `z2=(2·wt2+wt1)·dx/6`; boundary
    crossing uses linear interpolation of both `wt` and `y`. Multi-group
    spanning loop replicates Fortran's precedence quirk on line 5338
    (`gsig[ig] = yl*zl + yr*zr/sumde`, parsed as `yl*zl + (yr*zr)/sumde`
    not `(yl*zl+yr*zr)/sumde`).
  - Threaded `weight_fn` through `_build_subsection_sensitivities_rm`.
  - Updated `apply_rescon!` to accept `iwt::Int` kwarg (default 2) and
    materialize `weight_fn = get_weight_function(iwt)`.
  - Updated the file-header doctring to document the iwt-aware path.

- `src/orchestration/modules/errorr.jl`
  - At the `apply_rescon!` call site, pass `iwt=params.iwt` so the user
    deck's iwt threads through.

- `test/validation/test_errorr_covcal_lb5.jl`
  - Flipped `@test_broken abs(jul_mt102[1,1] - 2.658914e-4) < 1e-7`
    to `@test ... < 1e-7`.

### Before / after

| Quantity | Phase 72b | Phase 72c | Reference |
|----------|-----------|-----------|-----------|
| C[1,1]   | 3.580e-4 (37% high) | **2.658408e-4** | 2.658914e-4 |
| C[1,9] sign | negative | negative | negative |
| Row-1 nonzero col count | 10 | 10 | 10 |
| Sign pattern (10 nonzero cols) | matches | matches | — |

`diff = -5.06e-8`, `|relative| = 1.9e-4`. Within the 1e-7 absolute
tolerance for the canary.

## Test results

```
errorr covcal LB=5 XS·flux-weighted collapse (T15 MT=77)  | Pass: 4/4
errorr covcal Bug A — NK cross-pair stubs (T15 MT=77)     | Pass: 5/5
errorr covcal Phase 71 RED — MT=102 row-1 (rescon canary) | Pass: 41/41
T22 leapr (light-water)                                   | BIT_IDENTICAL 4636/4636
```

## Surprises

- The fix is "purely scale" on C[1,1] but the diff is tighter than expected
  (−5e-8 absolute, vs the ~1e-9 reference precision). The residual 0.019%
  is plausibly a combination of:
  - Julia's tanh-stretched pointwise grid vs Fortran's adaptive `rpendf`
    stepper — different finite-difference accuracy at the σ extrema.
  - Subtle constant difference in `weight_functions.jl`: Fortran iwt=6
    uses `wt6b=1.57855e-3` (5 digits); Julia uses `1.578551e-3` (6 digits).
    At E=bb=0.108 the Julia `cc = wt6b·exp(2)/bb²` evaluates to 1.000003
    where Fortran (iwt=6 branch) hard-codes `cc=1`. Fractional effect
    ~3e-6 in the thermal region only.
  None of these are critical for the canary tolerance; flagged for a
  future bit-faithful FP-grind if/when 1e-9 becomes the bar.

- The Fortran rpxgrp multi-group-spanning case has a precedence-quirk
  `gsig(j,ig) = yl*zl + yr*zr/sumde` on line 5338 that looks like a
  Fortran bug (obvious intent is `(yl*zl + yr*zr)/sumde`). On T15 this
  path is not exercised (all NSRS resonances are inside group 1 of
  LANL-30), so it doesn't affect the canary. Replicated verbatim per
  CLAUDE.md Rule 1 / Law 2 anyway, in case future tests do exercise it.

## Follow-up

- (P3) Bit-faithful FP-grind for the residual −5e-8 C[1,1] gap.
- (P3) Verify `weight_functions.jl` iwt=6 `cc=1` discrepancy against
  Fortran egtwtf precisely — likely a 6-digit constant copy-paste error.
- (P2) URR range cov port (`read_mf32` currently skips MAT=9237 range 11
  per the existing @info log line).
