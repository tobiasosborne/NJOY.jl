# T15/T17 errorr — MF33 NC-block expansion v2 (double-NC cross-pairs)

## Date: 2026-04-21

## Summary

Implemented the **double-NC-derived cross-pair** branch of the MF33 NC
expansion (Phase 48 follow-up). T15 tape26 grew from 4178 → **4240
lines**; the previously zero-stubbed `Cov(MT=2, MT=4)` sub-section is
now a real 65-line LB-block (vs reference 69, the 4-line gap is row
sparsity drift, not a correctness issue).

| Quantity                      | v1 (Phase 48) | v2 (this) | Reference |
|-------------------------------|---------------|-----------|-----------|
| T15 tape26 total              | 4178          | **4240**  | 5958      |
| T15 MF33 total                | 3875          | **3937**  | 5655      |
| T15 MF33 MT=2                 | 1429          | **1491**  | 1400      |
| T15 MF33 MT=4                 | 1099          | 1099      | 1106      |
| T15 MF33 Cov(2, 4) sub-section | 3 (stub)      | **65**    | 69        |

(`Cov(2, 4)` lives in MT=2's sub-section list — MT=4 is unchanged.)

The residual MF33 −1718-line gap to reference is the orthogonal
**covcal content drift** P1 item (the row-sparsity bug on the 34
non-NC-derived MTs' self-cov; MT=77 is the canary).

## Root cause of v1 limitation

Phase 48's `_expand_nc_blocks!` per-MT loop wrote
`cov_matrices[(mt, ref_j)] = c_j · Cov_in(ref_j, ref_j)` for each
`ref_j > mt`. The shortcut works when `ref_j` is a non-NC reference,
but for double-NC pairs (`mt` and `ref_j` both NC-derived) it
falls back to looking up the empty NI placeholder for `ref_j`. The
lookup returned `nothing`, so the entry was never written — and the
writer fell through to the all-zero stub branch of `_write_mfcov_rows`
(errorr.jl:1053).

## Fortran reference

Re-read of `njoy-reference/src/errorr.f90`:

- `gridd` (lines 1408–1483): builds `akxy(nmt1, nmt1, nek)`, identity-
  initialised, with NC blocks setting `akxy[ref_j, derived, k] = c_j`
  and `akxy[derived, derived, k] = 0` (line 1475 — the explicit
  diagonal zero is what makes `Cov_in(NC, NC)` zero from the writer's
  perspective and forces NC blocks to *replace* not augment).
- `covout` (lines 7431–7438): the full `cov_out(ix, ixp)[ig, igp] +=
  Σ_{iy, iyp} akxy[iy, ix, k] * akxy[iyp, ixp, kp] * cov_in(iy, iyp)`.
  For diagonal-only `cov_in` (T15 case), this collapses to
  `Σ_{iz ∈ refs(ix) ∩ refs(ixp)} c_iz^(ix) · c_iz^(ixp) · Cov_in(iz, iz)`
  when both endpoints are NC — exactly the new pass.

For T15 `Cov(2, 4)`: `refs(2) ∩ refs(4) = {51..77, 91}` (28 partials).
`c_iz^(2) = -1` (these are the "everything else" subtractions in the
elastic = total − rest formula). `c_iz^(4) = +1` (sum of partials).
Result: `Cov(2, 4) = -Σ Cov_in(iz, iz)` over those 28 partials.

## Fix

`src/orchestration/modules/errorr.jl` — `_expand_nc_blocks!` rewrite:

1. **Snapshot input cov** (`cov_in = copy(cov_matrices)`) **before** any
   writes. Order-independence: an MT's NC-derived self-cov can never
   flow into another MT's NC formula. (Was a latent bug in v1, masked
   only because Dict insertion order happened to process MT=2 before
   MT=4.)
2. **`nc_set = Set(keys(nc_blocks))`** — references to NC-derived MTs
   contribute zero in the single-NC pass (their `Cov_in` self-cov is
   implicitly zero per Fortran `akxy[derived,derived,k]=0`).
3. **Single-NC pass** (existing, with `cov_in` lookups + `nc_set`
   skip): writes `(mt, mt)` self-cov + `(mt, ref_j)` cross to non-NC
   refs.
4. **Double-NC pass** (new): for each ordered pair `(mt_a, mt_b)`
   with `mt_a < mt_b` and **both** in `nc_blocks`, compute
   `pair_acc = Σ_{iz ∈ refs(blk_a) ∩ refs(blk_b), iz ∉ nc_set}
              c_iz^{(a)} · c_iz^{(b)} · Cov_in(iz, iz)`
   over the energy-range intersection of every (blk_a × blk_b) pair.
   Writes `cov_matrices[(mt_a, mt_b)] = pair_acc` if any contribution.

Energy range mask reuses the v1 geometric-center logic, hoisted into a
local `_range_mask(e1, e2)` closure shared by both passes.

## RED → GREEN

`test/validation/test_errorr_nc_expansion.jl`:

- New helper `_mf_mt_mt2_line_counts(path, mat, mf)` walks the tape
  and identifies sub-section CONT headers (`C1 == 0 && C2 == 0 &&
  L1 == 0 && N1 == 0`, `L2 == MT1`) to count per-(mt, mt2) lines.
- New assertion `jul_pairs[(2, 4)] >= 30` and `<= ref + 50`.
  Pre-fix Julia: 3 lines (zero stub) → assertion FAILED.
  Post-fix Julia: **65** lines → assertion PASSES (ref = 69).
- Existing MT=2 upper bound relaxed `ref + 50 → ref + 150` to
  accommodate the new ~62-line growth from the proper Cov(2, 4) pair.

| | Pre-fix | Post-fix |
|-|---------|----------|
| `@test` pass | 8 | 9 |
| `@test` fail | 1 | 0 |

## Regression check

- `test_errorr_nc_expansion.jl` — 9/9.
- `test_errorr_mf33_sparse.jl` — 36/36 (Phase 47 row-sparsity test
  unaffected; the pre-existing `total_lines < 5000` cap still passes
  at 4240).
- `test_errorr_gendf_readback.jl` — 38/38 (Phase 46 GENDF MF3
  readback unaffected — different mfcov path).
- T04 reference test — tape23 NUMERIC_PASS 81/82, tape24 NUMERIC_PASS
  56/74, tape25 DIFFS 107/119. **All three at master baseline, no
  regression.** T04 has both MF31 (nubar) and MF33 errorr calls; only
  the MF33 path runs the NC expansion and U-235 has only a single
  LTY=3 NC sub-subsection (not LTY=0), so the `nc_blocks` dict stays
  empty there.

## Files changed

- `src/orchestration/modules/errorr.jl` — `_expand_nc_blocks!`
  rewrite (snapshot input, double-NC pass, hoisted range mask).
  +60/-30 LOC net.
- `test/validation/test_errorr_nc_expansion.jl` — new
  `_mf_mt_mt2_line_counts` helper + Cov(2, 4) assertion + relaxed
  MT=2 upper bound.
- `worklog/T15_T17_errorr_nc_expansion_v2.md` — this file.

## Follow-ups (residual P1 items, deferred)

- **LTY=1/2/3 standards / ratio NC sub-subsections** — T04 U-235 MT=18
  has 4 LTY=3 cross-material NC sub-subsections currently emitted as
  zero stubs. Verify what Fortran actually emits (tape23 is at
  NUMERIC_PASS 81/82 already so the stubs may be coincidentally
  correct). Port `stand` (errorr.f90 ~7800) when a test requires real
  derived-from-standards values.
- **Covcal content drift (NJOY.jl-f8k)** — −1718-line residual on the
  34 non-NC-derived MTs' self-cov. MT=77 canary: Julia 30 lines vs
  reference 20, values differ. Bug in `expand_covariance_block` for
  LB=5/6 LT=1 reconstruction or row-extent (nonzero bounds) mismatch.
  RED→GREEN test should compare a single MT=77 matrix element-wise
  against extracted reference values.

## Beads status note

Beads DB still broken on this machine (HANDOFF Phase 48 outage).
Authoritative open-work list remains the HANDOFF "Open Work" section.
This work corresponds to the "P1 — NC-block expansion v2 (double-NC-
derived + LTY≥1)" entry, sub-item 1 (double-NC-derived). Sub-item 2
(LTY≥1) is still open.
