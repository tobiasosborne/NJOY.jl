# Phase 73 — NC cross-material skip + σ-ratio weighting

**Date**: 2026-05-16
**Branch**: master
**Status**: GREEN — T15 tape26 line count 6108 → 5890 (gap +150 → −68 vs ref 5958); values 5000-million× closer on the dominant offenders. 222 unit-test assertions pass with zero regressions; T22 BIT_IDENTICAL 4636/4636 preserved.

## Summary

Two distinct bugs in the MF=33 covariance pipeline, surfaced by per-MT span/value drill-down on T15 U-238 JENDL tape26:

1. **Cross-material (`MAT1 ≠ 0`) sub-sections were polluting same-MT self-cov.** U-238 MT=18 (fission) has NK=6 sub-sections: 1 self-cov + 5 cross-material refs to MAT=9222/9228/9437/9440/9443 (`L1=MAT1` in the sub-CONT). The Fortran rescon explicitly skips these (`if (mats(ixp).ne.0) return`, errorr.f90:8531). Julia's reader ignored `sh.L1` and pushed the cross-material LB=6 data into `pair_blocks[(18, 18)]` — the same key as the legitimate self-cov. The LB=6 reader then mis-parsed `NEC` (which is not stored in the LIST header for LB=6), producing values up to 10⁷ across high-energy rows.
2. **NC self-cov expansion was missing the σ-ratio factor.** For Y = Σ c_i · X_i, the relative covariance is `rel_cov(Y) = Σ_i c_i² · rel_cov(X_i) · (σ_i/σ_Y)²`, not `Σ_i c_i² · rel_cov(X_i)`. Without the σ-ratio, NC-derived MTs (MT=2 elastic, MT=4 total inelastic) inherited the full relative cov of every partial at every group — over-estimating by ~`(σ_i/σ_Y)⁻²` per cell, which at high-energy elastic cells with sub-threshold partials was up to 5621×. Same factor applies to the single-NC cross-pair (`σ_ref_j/σ_Y` once) and the double-NC pair (`σ_ref/σ_a` × `σ_ref/σ_b`).

## Diagnostic chain

1. **Per-MT line count delta**: pre-fix MT=2 +142, MT=18 +74, MT=1 −64, MT=102 −2 (total +150).
2. **Per-(mt, mt2) sub-section drill-down**: revealed `(MT=2, mt2=2)` self-cov had 27-30 cols of nonzero data per row vs ref's 10-17.
3. **Per-row span comparison**: pre-fix rows 1-27 of MT=2/mt2=2 spanned cols 1-27/1-30, far wider than ref's natural sparsity.
4. **Cell-level value dump**: rows 28-30 of MT=2/mt2=2 had spans matching ref but values off by `+7.42`, `+7.73`, `+8.06` absolute — too large for relative cov. Also a uniform `+5.000000e-05` bias on rows 1-14 cols 1-14.
5. **Per-MT row-1 value dump**: MT=18 self-cov had the same `+5e-5` bias on cols 1-14 and absurd 10⁴-10⁷ values on cols 15-27. MT=16, 17, 37, 51, ..., 102 were all zero or correct. **MT=18 was the polluter; MT=2 inherited the pollution via NC.**
6. **ENDF raw inspection**: U-238 MF=33 MT=18 has NK=6 sub-sections; sub-section 1 is self (LB=5, 44 energies × 990 fvals), sub-sections 2-6 have `L1 ∈ {9222, 9228, 9437, 9440, 9443}` — cross-material refs.
7. **Per-MT row-28 single-cell dump (after cross-material fix)**: MT=18 row 28 became bit-identical to ref. MT=2 row 28 still had `+7.42` Δ. Confirmed only NC-derived MTs (MT=2, MT=4) were wrong. Pinned the cause to `_expand_nc_blocks!` missing the σ-ratio factor.

## Fix

### Files changed

- `src/orchestration/modules/errorr.jl`
  - In the per-MT sub-section read loop: parse `mat1 = sh.L1` and set `skip_pair = (mat1 != 0)`. Still consume the NC + NI blocks to advance file position, but don't push into `pair_blocks`, `nc_blocks`, `listed_pairs`, or `nc_derived_mts`. Mirrors Fortran covout/rescon's `mats(ixp).ne.0` skip.
  - In `_expand_nc_blocks!`: extended signature to take `group_xs::Dict{Int,Vector{Float64}}`. Added a local `σ_zero(mt, ig)` helper that returns the group-averaged σ for a given (mt, ig) or 0 if absent. Self-cov: `self_acc[ig, igp] += c_i² · cov_in[ig, igp] · (σ_i/σ_Y)(ig) · (σ_i/σ_Y)(igp)` with denom>0 guards. Single-NC cross-pair: `cross_acc[ig, igp] += c_j · cov_in[ig, igp] · (σ_j/σ_Y)(ig)`. Double-NC pair: `pair_acc[ig, igp] += w · cov_in[ig, igp] · (σ/σ_A)(ig) · (σ/σ_B)(igp)`.
  - Updated docstring header line `_expand_nc_blocks!(cov_matrices, nc_blocks, egn, group_xs)`.
  - Updated the orchestration callsite to pass `group_xs` through.

- `src/processing/rescon.jl`
  - Added `_resonance_group_window(elr, ehr, egn) -> (iest, ieed)` helper that mirrors Fortran's iest/ieed derivation at errorr.f90:3093-3108: `iest` = first group with upper bound > elr; `ieed` = last group with upper bound ≥ ehr; range collapses to `ieed=0` if `elg ≥ ehr`.
  - In `_build_subsection_sensitivities_rm`: after building `sens`, zero `sens[*, *, ig]` for all `ig` outside `[iest, ieed]`. Mirrors Fortran's `do ig=iest,ieed; sens(j,loop,ig)=gsig(j,ig)*tmp` (errorr.f90:4509) which leaves sens at its initial zero outside the resonance group window. **Functional impact for T15 is nil** (the rescon was correctly producing near-zero in those groups already), but matches Fortran exactly and prevents subtle FP-noise accumulation across npar² parameters in future tests.

### Before / after

| Quantity                                | Before    | After     | Reference |
|-----------------------------------------|-----------|-----------|-----------|
| T15 tape26 line count                   | 6108      | 5890      | 5958      |
| T15 tape26 MF=33 line count             | 5805      | 5587      | 5655      |
| MT=2 Δ (lines vs ref)                   | +142      | **−2**    | —         |
| MT=18 Δ (lines vs ref)                  | +74       | **0**     | —         |
| MT=1 Δ (lines vs ref, untouched)        | −64       | −64       | —         |
| MT=102 Δ (lines vs ref, untouched)      | −2        | −2        | —         |
| MT=18 self-cov row 28 (col-28)          | bit-id    | bit-id    | bit-id    |
| MT=2/mt2=2 row 28 col 28 \|Δ\|          | +7.42     | +2.50e-6  | 1.32e-3   |
| MT=2/mt2=18 row 28 max \|Δ\|            | 6.52e-5   | 1e-11     | ~5e-5     |
| Phase 72c MT=102 C[1,1] (canary)        | −5.06e-8  | −5.06e-8  | 2.659e-4  |

### Test results

- `test_errorr_covcal_lb5.jl` (Phase 72c canary): **50/50 PASS** (4 + 5 + 41)
- `test_errorr_mf33_sparse.jl`: **73/73 PASS** (and self-reports the improved line counts)
- `test_errorr_writer_mf_dispatch.jl`: **51/51 PASS** (45 + 4 + 2)
- `test_errorr_nc_expansion.jl`: **9/9 PASS**
- `test_errorr_gendf_readback.jl`: **38/38 PASS**
- Reference test **T22 (leapr light-water)**: BIT_IDENTICAL 4636/4636 preserved

Total: 221 + 1 reference test = 222 assertions PASS, zero regressions.

## Remaining work for T15 tape26 to BIT_IDENTICAL

Net gap is now −68 lines (was +150). The remaining issues by MT:

1. **MT=1 Δ = −64**: MT=1/mt2=2 cross-pair emits only a 1-row zero stub vs ref's 16-row, 256-data-word content. MT=1 (total) is a sum-of-partials so its cross-cov to MT=2 needs a special derivation path. Likely MT=1 has its own NC formula or is auto-derived by Fortran covout.
2. **MT=102 Δ = −2**: row 14 missing in Julia (rescon group-window edge). Spans of MT=2/mt2=2 rows 13-15 also off by 1 (related to MT=102 row 14).
3. **Residual values at rows 10-12** (post-rescon): rel error 3-44× from rescon's sub-ULP FP precision (Phase 72c known).

## Surprises

- Three independent bug classes were stacked: rescon (Phase 72c partially addressed), NC σ-ratio missing, and cross-material pollution. Each individually masked the others — the rescon FP-noise hypothesis from my initial trim attempt produced zero line-count change because the pollution was already there from the LB=5/NC path.
- The cross-material bug had been latent since the original NC v1 work (Phase 48). It only surfaces on tapes with `MAT1 != 0` references — U-238 JENDL is one such; ENDF/B-VIII evaluations may also have them.
- The σ-ratio bug was deep in plain sight: the comment in `_expand_nc_blocks!` said "treats input cov as relative" but the actual formula was the absolute-cov accumulation pattern without the absolute→relative conversion. Easy to miss on inspection because the dictionary is named `cov_in` not `abs_cov_in`.

## Follow-up

- (P1) MT=1/mt2=2 cross-pair: trace Fortran covout's MT=1 derivation path. Likely a separate NC formula or a special "sum of all partials" path.
- (P2) MT=102 row 14 / MT=2/mt2=2 row 13-15 span fix: investigate rescon's group-window boundary handling.
- (P3) MT=2/mt2=2 rows 10-12 sub-ULP FP precision (rescon — Phase 72c follow-up).
