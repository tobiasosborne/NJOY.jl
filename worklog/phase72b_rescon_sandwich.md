# Phase 72b — rescon sensitivity builder + sandwich (canary GREEN-ish)

**Date:** 2026-05-08
**Test:** T15 (U-238 JENDL, MAT=9237) — Phase 71 RED canary
**Goal:** Land the sensitivity builder + J·cov_RP·J^T sandwich on top
of the Phase 72 foundation, drive `MT=102 row 1` from all-zero to
ref-matching.

## Outcome

**Sandwich landed, sign pattern matches ref, magnitudes within ~30%
on average; one column at +69% needs FP-order alignment.** The Phase
71 RED canary suite: 40 PASS + 1 @test_broken (down from 18 PASS +
2 BROKEN in Phase 72).

| MT=102 row 1 | Phase 71 | Phase 72 (foundation) | Phase 72b (sandwich) |
|--------------|----------|------------------------|------------------------|
| nonzero cols | 0/10     | 0/10                   | **10/10** ✓            |
| sign pattern | n/a      | n/a                    | matches ref ✓          |
| C[1,1] ratio | 0        | 0                      | **1.346** (ref 2.66e-4)|
| C[1,9] sign  | n/a      | n/a                    | **negative ✓**         |
| max ratio    | n/a      | n/a                    | 1.687 (col 6)          |

The negative-sandwich signature `C[1,9] < 0` GREENs — the test that
was `@test_broken` since Phase 71. This is the canonical "rescon is
actually doing the sandwich, not just MF=33 LB=5 lookup" signal.

The 1e-7 tolerance on `C[1,1]` stays `@test_broken`; closing the
residual ~30% systematic is FP-grind work for next session (almost
certainly a pointwise-grid density / `rpendf` adaptive-mesh detail
or `rpxgrp` accumulation order — see Open Issues below).

## What landed

### 1. `src/processing/rescon.jl` — full sensitivity + sandwich port

Replaced the no-op `apply_rescon!` skeleton with the full Fortran-
faithful pipeline. ~400 LOC delta. Mirrors `rpxlc12`
(errorr.f90:4399-4593) for LCOMP=1 / LRF=3 / U-238 JENDL.

**`_build_pointwise_grid`** — rpendf-equivalent (errorr.f90:5015-5089):
union of (a) log grid 200 pts/decade across [elr, ehr], (b) per-
resonance dense linear refinement spanning ±30·Γ with tanh-stretched
density-toward-peak (~33 pts/Γ in wings, finer at ER), (c) all output
group boundaries within range. Yields ~6000 points for U-238 first
range [1e-5, 1000] eV with 26 resonances.

**`_group_average_flat`** — rpxgrp-equivalent (errorr.f90:5227-5367)
for iwt=2 (T15 standard). Trapezoidal integral with constant flux
weight `w(E) = 1` (matches `egtwtf` for iwt=2 → `wtf=1`,
errorr.f90:10108-10110). Splits segments at output-group boundaries
and accumulates per-group with linear interpolation at boundary
crossings.

**`_match_mf32_to_mf2_rm`** — MF=32 ↔ MF=2 resonance pairing per
errorr.f90:4413-4443 tolerances `|Er/Er2 - 1| < 1e-6` and
`|AJ - AJ2| < 1e-4`. Returns `(l, n)` index into the typed RM
parameter struct.

**`_perturb_rm`** — clones a `ReichMooreParameters` and perturbs one
parameter for resonance `(l, n)`. `loopn` 1..5 maps to ER/GN/GG/GFA/
GFB per the LRF=3 layout (errorr.f90:4459-4467). Deep-copies only
the affected l-state vectors (immutable elsewhere — the existing
`cross_section_rm` is allocation-free per call).

**`_rm_pert_factors` / `_rm_gwidth`** — perturbation magnitudes:
ER ±1e-4 (factor 0.9999/1.0001), widths ±1e-2 (factor 0.99/1.01).
Verbatim from errorr.f90:4448-4488 + 4500-4505. Constants exposed
as `_RESCON_PERT_ER` and `_RESCON_PERT_WIDTH`.

**`_eval_xs_grid_rm`** — calls existing `cross_section_rm` at each
grid point to produce `(σ_total, σ_elastic, σ_fission, σ_capture)`
vectors. Reuses Julia's RM evaluator (which recomputes penetrability
from `er` per resonance, so `er` perturbation propagates correctly
without requiring the Fortran ipos/ser/per scratch updates at
errorr.f90:4453-4457).

**`_build_subsection_sensitivities_rm`** — outer loop over
`(loopm, loopn) ∈ NRB × MPAR`, central-difference per
errorr.f90:4500-4505: `tmp = (σ_+ - σ_-) / (2·gwidth)`. Group-
average via `_group_average_flat` per channel, store in
`sens[1..4, p, ig]`. For unmatched MF=32 resonances (rare —
indicates MF=2/MF=32 inconsistency), warn + zero the sens entries
rather than error; for zero-width parameters (e.g. fission widths
in U-238), skip per Fortran's `if (gwidth.ne.zero)` guard at
errorr.f90:4515-4520.

**`_apply_sandwich!`** — for one `(mt, mt2)` pair, accumulates
`abs_cov[ig, ig'] = Σ_{i ≤ j} cov[i,j] · sens[i,ig] · sens[j,ig'] +
(i ≠ j ? cov[i,j] · sens[j,ig] · sens[i,ig'] : 0)`. Triangular
walk for self-pair, full ngn×ngn for cross-pair. Mirrors
errorr.f90:4541-4593. Then converts to RELATIVE cov by dividing by
`σ̄_row(ig) · σ̄_col(ig')` — the convention `cov_matrices` already
holds (Phase 51 LB=5 weighted-collapse is also relative).

**`apply_rescon!`** orchestration: read MF=32 + MF=2, match isotopes
+ ranges, dispatch each subsection through the sensitivity builder
+ sandwich, prefer `group_xs[mt]` (Fortran-equivalent group σ̄ from
GENDF) over locally-built σ̄ for the relative-cov denominator —
falls back to local σ̄ when `group_xs` lacks the channel-MT entry.

### 2. `test/validation/test_errorr_covcal_lb5.jl` — Phase 71 canary tightened

Previous (Phase 71): 18 PASS + 2 @test_broken with informational
non-zero count log.

Now (Phase 72b):
- Structural: `nz_jul == nz_ref` (10/10 for MT=102 row 1) ✓
- Sign pattern: `sign(jul[1, j]) == sign(ref[1, j])` for all
  non-zero cols ✓
- Magnitude: `abs(jul/ref - 1) < 1.0` (factor of 2) ✓
  - Loose first-iteration bar; tightening ladder: 1.0 → 0.1 → 0.01
    → 1e-7 (bit-faithful).
- `jul_mt102[1, 9] < 0` (negative-sandwich signature) ✓ — was
  `@test_broken`, now passes.
- `abs(jul[1, 1] - 2.658914e-4) < 1e-7` — still `@test_broken`,
  current ratio 1.346.

### 3. T22 BIT_IDENTICAL regression check

T22 (light-water leapr, BIT_IDENTICAL 4636/4636) preserved. The
rescon path only fires when `mfcov == 33 && _mf32_present(mat)` —
T22 doesn't trigger it.

## Numerical analysis (current state)

Decomposition of `abs_cov[1,1]` for MT=102 capture self-cov, U-238
range 1, group 1 (LANL-30 [1.39e-4, 0.152] eV):

| contribution                 | value          |
|------------------------------|----------------|
| diagonal `Σ cov[i,i]·s_i²`   | +1.144e-2      |
| off-diagonal cross-term      | -9.768e-3      |
| **abs_cov[1,1] (sum)**       | **+1.672e-3**  |
| σ̄_cap(g=1) (group-averaged) | 2.161 b        |
| **rel_cov[1,1] = abs/σ̄²**  | **+3.580e-4**  |
| reference                    | +2.659e-4      |
| ratio jul/ref                | 1.346          |

The result is the **small difference between two large near-equal
numbers** — ~15% of the diagonal magnitude. Sensitive to small
biases in `sens` accumulation. Doubled grid density (50 → 200
pts/decade + 40 → 200 pts/resonance with tanh-stretched peak
sampling) tightened the per-column ratios from 0.55..2.38 to
0.76..1.69, but the [1,1] ratio held at 1.346 — that's a SYSTEMATIC,
not noise.

## Open issues for next-session FP-grind

1. **Pointwise-grid bias.** Current grid is "engineered" rather than
   matching Fortran's adaptive `eskip1/2/3/4` mesh exactly. The
   eskip values are module-level Fortran parameters keyed off
   resonance proximity (E ∈ [0.7·Er, 1.3·Er] uses eskip3, etc.) —
   port these constants verbatim and use the same multiplicative
   stepper at errorr.f90:5066-5085 to remove any grid-shape bias.

2. **Trapezoidal vs Simpson.** Fortran rpxgrp uses a specific
   piecewise-linear segment integration with explicit boundary-
   crossing handling at lines 5301-5360 (the `else if (x1.gt.egnt1)`
   branch with `wt12=wt1-wt2` linear interpolation). My
   `_group_average_flat` uses a similar split-segment trapezoidal
   but does NOT replicate the cubic moment formulas
   `z1=(2·wt1+wt2)·x12/six` at lines 5284-5285 — this is a
   higher-order correction Fortran uses for non-flat weights.
   For iwt=2 (`wtf=1` flat), z1=z2=x12/3, so the formula reduces
   to weighted-trapezoidal — but the FP accumulation order may
   still differ. Port it verbatim.

3. **Penetrability cache.** Fortran's rpxlc12 stores per-resonance
   `ser`, `per` in `b(ipos)`, `b(ipos+1)` AT MATCH TIME and updates
   them when ER is perturbed (errorr.f90:4453-4457). The RM
   evaluator (`ggrmat`) reads these cached values rather than
   recomputing. Julia's `cross_section_rm` recomputes from `er`
   each call — that's *more accurate* but introduces a small
   numerical-recipe difference in how penetrability rounds. For
   bit-faithful output this needs investigation.

4. **`cflx` and `abn`.** Fortran's
   `sens(j,p,ig) = gsig(j,ig) · cflx(ig) · abn` (errorr.f90:4512)
   carries a per-group-flux factor that (per units analysis) cancels
   against the σ̄·σ̄ relative conversion in the writer — but the
   exact accumulation order may matter. Re-check whether absorbing
   `cflx · abn` into `sens` and the writer-side `1 / (cflx_g·cflx_g')`
   normalization gives a different answer at the FP-rounding level.

5. **MPAR > 3 cases.** U-238 has MPAR=3 (ER, GN, GG). Other
   evaluations (e.g. U-235 ENDF/B-VIII) have MPAR=4 or 5 with
   uncertain GFA/GFB. The dispatch in `_perturb_rm` covers
   loopn=1..5; tested only loopn=1..3 against U-238.

## Files

NEW:
- `worklog/phase72b_rescon_sandwich.md` (this file).

MODIFIED:
- `src/processing/rescon.jl` — sensitivity builder + sandwich
  (~400 LOC, replacing the Phase-72 no-op skeleton).
- `test/validation/test_errorr_covcal_lb5.jl` — Phase 71 canary
  tightened to assert structural + sign + factor-of-2 magnitude;
  `[1, 9] < 0` flipped from `@test_broken` to `@test`.

## Fortran source citations

- `errorr.f90:4248-4523` — rpxlc12 perturbation loop (the heart).
- `errorr.f90:4500-4505` — central-difference per-pointwise tmp.
- `errorr.f90:4508-4513` — rpxgrp call + cflx·abn scaling.
- `errorr.f90:4541-4593` — sandwich accumulation into
  cff/cgg/cee/ctt/cef/ceg/cfg.
- `errorr.f90:5015-5089` — rpendf (adaptive grid).
- `errorr.f90:5227-5367` — rpxgrp (group-averaging).
- `errorr.f90:8513-8819` — rescon (the structurally-trivial cova
  add — once sens + cov fills exist).
- `errorr.f90:10023-10110` — egtwtf (iwt=2 → wtf=1 weight).
