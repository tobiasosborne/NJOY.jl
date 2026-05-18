# Phase 74 — NC cross-pair to lower-MT (MT=1/mt2=2 closure)

**Date**: 2026-05-18
**Branch**: master
**Status**: GREEN — T15 tape26 line count 5890 → 5952 (gap −68 → **−6** vs ref 5958); MT=1/mt2=2 sub-section bit-identical on geometry (65 lines / 16 rows, matches ref). Phase 72c/73 canaries preserved. Zero regressions across 5 errorr test files + T22 reference test.

## Summary

Phase 73 closed cross-material pollution + the missing NC σ-ratio. The residual −68 line gap on T15 tape26 was dominated by MT=1/mt2=2 emitting a 2-line zero stub vs ref's 65-line populated cross-pair.

Root cause: `_expand_nc_blocks!` single-NC cross-pair pass filtered `ref_j > mt || continue`, skipping MT=2's NC reference to MT=1 (since 1 < 2). U-238 JENDL's MF=33 MT=2 carries an NC sub-subsection with formula `MT=2 = +1·MT=1 − Σ_partials c_i·MT_i` over `[1e-5, 2e7]` eV (J33U238 lines 16555-16568). By covariance symmetry, the cross-pair `Cov(MT=2, MT=1)` is equivalent to `Cov(MT=1, MT=2)^T`, which the writer keys on `cov_matrices[(1, 2)]` (min-MT-first convention).

## Fortran reference

Covout's sandwich loop (errorr.f90:7406-7438) walks the full input covariance tape for every output (ix, ixp). With MT=1 evaluated (akxy diag = 1) and MT=2 derived (akxy diag = 0 since MT=2 has its own NC formula), the loop accumulates:

```
cova(igp, ig) += akxy(iy, MT1, k) · akxy(iyp, MT2, kp) · cov(jgp)
```

For the input (iy=MT=1, iyp=MT=1) self-cov record, this picks up `akxy(MT1, MT1, k) = 1` and `akxy(MT1, MT2, kp) = +1` (from MT=2's NC formula), giving `+1 · Cov_abs(MT=1, MT=1)`.

In relative-cov terms with min-first storage convention:
```
rel_cov(MT=1, MT=2)[ig, igp] = c · rel_cov(MT=1, MT=1)[ig, igp] · σ(MT=1, igp) / σ(MT=2, igp)
```
(σ-ratio applies to the **column** when the derived MT is the column side — the transpose of the row-side σ-ratio that the existing `ref_j > mt` branch uses).

## Fix

`src/orchestration/modules/errorr.jl` — `_expand_nc_blocks!` single-NC cross-pair pass split into two storage branches keyed by `ref_j` vs `mt`:

- **`ref_j > mt`** (existing): store at `(mt, ref_j)` with row σ-ratio `σ(ref_j, ig) / σ(mt, ig)`.
- **`ref_j < mt`** (new): store at `(ref_j, mt)` with column σ-ratio `σ(ref_j, igp) / σ(mt, igp)`.

Two `cross_acc` dicts (`_hi` / `_lo`) collect per-`ref_j` contributions; at the end the `_hi` branch overwrites `cov_matrices[(mt, ref_j)]` (no other writer touches this key) while the `_lo` branch accumulates `+=` into `cov_matrices[(ref_j, mt)]` (multiple NC-derived MTs could collide; matches Fortran's additive `cova` semantics).

## Files changed

- `src/orchestration/modules/errorr.jl` — `_expand_nc_blocks!` cross-pair split (+47, −12 LOC including the formula comment block).
- `test/validation/test_errorr_mf33_sparse.jl` — new TDD assertions for MT=1/mt2=2 sub-section row count (`>= 12`) and total tape26 line count (`abs(total - ref) <= 15`).

## Before / after

| Quantity                                | Before    | After     | Reference |
|-----------------------------------------|-----------|-----------|-----------|
| T15 tape26 line count                   | 5890      | **5952**  | 5958      |
| Net gap vs ref                          | **−68**   | **−6**    | 0         |
| T15 tape26 MF=33 line count             | 5587      | 5649      | 5655      |
| MT=1/mt2=2 sub-section lines            | 3 (stub)  | **65**    | 65        |
| MT=1/mt2=2 row records                  | 1 (stub)  | **16**    | 16        |
| MT=2 line-count Δ                       | −2        | −2        | —         |
| MT=4 line-count Δ                       | 0         | 0         | —         |
| MT=102 line-count Δ                     | −2        | −2        | —         |
| Phase 72c MT=102 C[1,1] canary          | preserved | preserved | 2.659e-4  |
| T22 BIT_IDENTICAL                       | 4636/4636 | 4636/4636 | —         |

## Test results

- `test_errorr_mf33_sparse.jl`: **75/75 PASS** (was 73/73 — 2 new MT=1/mt2=2 assertions, both pass)
- `test_errorr_covcal_lb5.jl` (Phase 72c canary): **50/50 PASS** (4 + 5 + 41)
- `test_errorr_writer_mf_dispatch.jl`: **51/51 PASS** (45 + 4 + 2)
- `test_errorr_nc_expansion.jl`: **9/9 PASS**
- `test_errorr_gendf_readback.jl`: **38/38 PASS**
- Reference test **T22 (leapr light-water)**: BIT_IDENTICAL 4636/4636 preserved @ rtol=1e-9

Total: 223 unit assertions PASS + T22 reference test, zero regressions.

## Remaining work for T15 tape26 to BIT_IDENTICAL

Net gap is now −6 lines. Remaining per-MT issues:

1. **MT=102 Δ = −2**: row 14 missing in Julia (rescon group-window edge — Phase 73 follow-up #2).
2. **MT=2 Δ = −2**: cells around rows 13-15 (likely same rescon edge interaction).
3. **MT=2/mt2=2 rows 10-12 sub-ULP FP precision** (Phase 72c follow-up #3, sub-line-count effect — does not affect line count, only data values).

The −6 final-line gap is plausibly the union of items 1 and 2 (2+2 = 4) plus a couple of edge cases in single rows being one data-line longer/shorter than ref. Likely within reach for a Phase 75 pass once the rescon group-window edge is fixed.

## Surprises

- The `ref_j > mt` filter was an artefact of the original v1 NC implementation (Phase 48), which assumed cross-pair storage strictly in `(mt, ref_j)` order without considering that *most* NC formulas reference a mix of higher and lower MTs. MT=2's `+1·MT=1` term is the textbook NC pattern (elastic = total − rest); the same fix will benefit any evaluation where a derived MT references the total or another lower-numbered evaluated MT.
- The fix is symmetric and one-directional in terms of storage: we never need to write to both `(a, b)` and `(b, a)` because the writer only reads `(min, max)`. The column σ-ratio derivation falls out cleanly from `rel_cov(A, B)[i, j] = rel_cov(B, A)[j, i]`.

## Follow-up

- (P1) MT=102 row 14 / MT=2 row 13-15 span: investigate rescon's `_resonance_group_window` boundary handling (Phase 73 follow-up #2).
- (P3) MT=2/mt2=2 rows 10-12 sub-ULP FP precision (Phase 72c follow-up #3).
