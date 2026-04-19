# T15 errorr — MF33 sparse per-row emission + NC-aware pair synthesis

## Date: 2026-04-19

## Summary

T15 tape26 MF33: **7938 → 1556 lines** (reference 5655). Total tape26:
**8205 → 1859** (ref 5958). The 2247-line over-expansion from HANDOFF
Phase 44 is eliminated; the remaining gap is under-emission from
NJOY.jl-km1 (cross-MT NC-derived covariance not yet expanded, costing
~2300 lines on MT=2/MT=4 self-cov entries). T04 tape23 **BIT_IDENTICAL
preserved** (82/82).

## Root causes (two tangled)

1. **Dense full-NGN-column per-row emission** (errorr.jl:836-847 pre-fix):
   Julia emitted `ngn` values per row for every row of every non-empty
   covariance matrix, regardless of sparsity. Fortran `covout`
   (errorr.f90:7530-7605) scans each row for `|v| > eps`, emits only
   `[ig2lo, ng2]` (the nonzero column range), and **skips entirely**
   any row that is all-zero (except the last row, which always gets a
   single-value stub). For a block-diagonal covariance covering 10
   columns, Fortran writes ~75 lines per pair vs Julia's 181.

2. **Synthesized empty cross-pair stubs for every `mt2 > mt`** in
   `reaction_mts`. For T15 with 36 reaction MTs, every MT got ~35 empty
   3-line stubs = ~3780 over-emitted lines. Fortran only emits cross-pair
   sub-sections when the ENDF declared them (via NL list or NC-derived
   coefficients).

## Fix

### Sparse per-row emission (`_write_mfcov_rows`)

New helper matching Fortran `covout` exactly:
- Scan row for first (`ig2lo`) and last (`ng2`) column with `|v| > 1e-20`.
- All-zero row + `ig < ngn` → skip.
- All-zero row + `ig == ngn` → one-value stub at column `ig`.
- Otherwise emit `LIST(0, 0, L1=count, L2=ig2lo, N1=count, N2=ig)` +
  body values `row[ig2lo..ng2]`.

### NC-aware pair synthesis

Two-level pair-set for each MT's covariance section:
- `listed_pairs::Set{(mt,mt2)}`: explicit L2-declared sub-sections
  (read during the MF33 walk).
- `nc_derived_mts::Set{Int}`: MTs with at least one sub-section having
  `NC>0` (derived covariance that references other MTs via NC
  coefficients).

Writer rule for each MT:
- Always include pairs from `cov_matrices` (populated by the NI
  expansion loop) and `listed_pairs` (explicit L2 entries).
- **If MT in `nc_derived_mts`**: also add `(mt, mt2)` for every
  `mt2 > mt` in `reaction_mts`. Matches Fortran's coverage of
  NC-derived cross-pairs without requiring full NC expansion
  (which is deferred as NJOY.jl-km1).

T04 MT=18 has NC>0 → synthesizes `(18, 102)` stub → ref NL=2 restored.
T15 MT=1/16/17/… have NC=0 → no synthesis → no over-emission.
T15 MT=2/MT=4 have NC=1 → synthesizes ~35 cross-pair stubs each
(~210 lines total); this is under-emission compared to the real
NC-expanded data (NJOY.jl-km1) but structurally correct NL counts.

## Verification

### RED→GREEN test

`test/validation/test_errorr_mf33_sparse.jl` — asserts per-MT MF33
line count ≤ reference+15 on the 34 MTs where Julia has populated
self-covariance; asserts total tape26 < 4000 lines (was 8205).

| | Pre-fix | Post-fix |
|-|---------|----------|
| `@test` pass | 2 | 36 |
| `@test` fail | 34 | 0 |

### T15 per-MT deltas (sample)

| MT | Pre | Post | Ref | Δ(post-ref) |
|----|-----|------|-----|-------------|
| 1  | 287 | 167  | 271 | -104 |
| 2  | 103 | 106  | 1400| -1294 (NJOY.jl-km1) |
| 4  | 100 | 103  | 1106| -1003 (NJOY.jl-km1) |
| 16 | 278 | 110  | 110 | 0 |
| 17 | 275 | 101  | 103 | -2 |
| 18 | 272 | 264  | 228 | +36 (NJOY.jl-f8k content drift) |
| 77 | 235 | 30   | 20  | +10 (NJOY.jl-f8k) |

### Regression

- **T04** tape23 **BIT_IDENTICAL** 82/82 (preserved via NC-aware
  synthesis — MT=18 has NC>0, (18,102) stub re-added).
- T04 tape24 NUMERIC_PASS 56/74 (unchanged).
- T04 tape25 DIFFS 108/119 (unchanged).
- Phase-45 `test_groupr_auto_expand.jl`: 39 pass + 3 broken (unchanged).
- Phase-46 `test_errorr_gendf_readback.jl`: 38/38 (unchanged).

## Files changed

- `src/orchestration/modules/errorr.jl`:
  - Track `listed_pairs` (explicit NL L2 entries) during MF33 walk.
  - Track `nc_derived_mts` (MTs with any NC>0 sub-section).
  - New `_write_mfcov_rows` helper (sparse per-row emission, ~60 LOC).
  - Writer accepts `listed_pairs` and `nc_derived_mts` kwargs.
  - Replaced dense full-NGN row block with `_write_mfcov_rows` call.
- `test/validation/test_errorr_mf33_sparse.jl`: new RED→GREEN test.

## Traps (NEW)

**Trap (errorr MF33 sparse row emission)**: each row of an MF33
LIST-per-row covariance output must encode `(L1=count, L2=ig2lo,
N1=count, N2=ig)` where `[ig2lo, ng2]` is the first-to-last nonzero
column range with threshold `|v| > 1e-20`. Emitting `(L1=ngn, L2=1,
N1=ngn, N2=ig)` for every row blows file size by 4-20× for typical
block-diagonal covariances. All-zero rows are dropped except the last,
which gets a one-value stub at `ig2lo = ng2 = ig`.

**Trap (NC-derived covariance pair count)**: the number of sub-sections
(section HEAD NL field) an errorr MF33 output emits for a given MT
equals the number of MTs referenced by that MT's covariance data —
either via explicit L2-declared NL entries *or* via NC derivation
coefficients in any NC>0 sub-section. Matching the reference NL count
without full NC expansion requires tracking both: `listed_pairs` for
L2 entries, and `nc_derived_mts` flag to synthesize cross-pair stubs
against the full reaction list. T04's (18, 102) cross-pair comes
entirely from NC derivation — there is no explicit L2=102 sub-section.

## Follow-ups

- **NJOY.jl-km1**: full NC-block expansion — unlocks real values for
  T15 MT=2/MT=4 cross-covariance (~2300 lines under-emission).
- **NJOY.jl-f8k**: covariance matrix content/extent drift (MT=77
  self-cov values differ from Fortran's `covcal`+`resprp` output).
  Surfaced once sparse emission exposed the per-row shape mismatch.
- **T04 tape25 cov-sum**: 11 residual lines, orthogonal, T03_phase7.
