# T15/T17 errorr — MF33 NC-block (derived covariance) expansion

## Date: 2026-04-20

## Summary

Implemented the LTY=0 (linear-combination) NC block expansion in
`errorr_module`. T15 tape26 grew from **1859 → 4178 lines** (+2319,
target was ~+2300). MT=2 self+cross-cov: **106 → 1429** (vs ref 1400).
MT=4: **103 → 1099** (vs ref 1106). Total MF33 contribution: **1556 →
3875** (vs ref 5655). Closes NJOY.jl-km1 v1; remaining gap (~1780
lines on the non-NC-derived MTs) is the orthogonal NJOY.jl-f8k
covcal content drift.

T04 reference test status preserved (NUMERIC_PASS 81/82 on tape23 was
already the master baseline before this change — HANDOFF Phase 47's
"BIT_IDENTICAL" claim was over-stated; my change does not alter that
status).

## Root cause

`src/orchestration/modules/errorr.jl:131-134` (pre-fix) consumed NC
sub-subsections with a one-line `read_list(io)` per NC entry — but an
NC sub-subsection is **two** ENDF records (a CONT giving LTY, then a
LIST with the formula body). The pre-fix code mis-parsed the file on
NC≥1 sub-sections; for sub-section #1 NI was processed correctly so
the only effect on existing tests was occasional misalignment for
later sub-sections that no test depended on.

The deeper issue: the actual derivation (`Cov(MT, MT) = Σᵢⱼ cᵢ cⱼ
Cov(refᵢ, refⱼ)`) was simply never computed. For T15/T17 U-238 every
MT=2 / MT=4 cross-pair came out as a 3-line zero stub.

## Fortran reference

`njoy-reference/src/errorr.f90`:

- `gridd` (lines 1091-1483) — read NC LISTs, build `akxy(nmt, nmt,
  nek)` coefficient array. Identity initialization (line 1429-1430);
  per-NC-block fill `akxy[ref_j, derived, k] = c_j` over the energy
  bin range (lines 1473-1475); explicit diagonal zeroing
  `akxy[derived, derived, k] = 0` (line 1475) so the explicit NI
  self-cov is *replaced* — not summed — by the NC formula.
- `covout` (lines 7431-7438) — accumulate
  `cov_out[ig,igp] += Σ_{iy,iyp} akxy[iy,ix,k] * akxy[iyp,ixp,kp] *
   cov_in[jg,jgp]`. This is C^T Σ C in matrix form.
- `gridd` errors hard if a referenced MT (XMTI) is not in the MF33 MT
  list (lines 1456-1466). No fallback to MF31/MF35.

## T15/T17 U-238 NC inventory

`njoy-reference/tests/resources/J33U238` (MAT=9237, JENDL-3.3):

- **Only MT=2 and MT=4 have NC sub-subsections.** Both use one
  LTY=0 block over the full energy range [1e-5, 2e7].
- **MT=4** (trivial sum): NCI=28, all coeffs +1, references
  MT51-77 + MT91 (`σ_n,n' = Σ partials`).
- **MT=2** (mixed sign): NCI=34, references MT=1 (+1) and 33 others
  (-1) (`σ_elastic = σ_total - everything else`).
- All referenced MTs have NC=0 / NI=1 (pure self-covariance).
- Reference T15 tape26: 5958 lines total; MT=2 contributes 1400
  (NK=35 sub-sections), MT=4 contributes 1106 (NK=34).

T17 = same NC structure + 2 cross-material pairs each (MAT 9228 U-235
MT=18, MAT 9437 Pu-239 MT=18) for MT=2 and MT=4.

## Fix

### `_read_nc_subsection(io)` — new helper

Reads CONT (LTY in L2) + LIST (E1, E2 in C1/C2, NCI in N2, body
`(c_j, MT_j)` pairs). Returns `nothing` for LTY ∈ {1,2,3}
(standards/ratio variants — common in T04 U-235 cross-material
sub-sections, deferred). Returns `NCSubSubsection` for LTY=0.

### Read-loop change

Replaces the discarding `try; read_list(io); catch; break; end` with
a parsing call that captures LTY=0 blocks into a per-MT
`Dict{Int, Vector{NCSubSubsection}}`.

### `_expand_nc_blocks!(cov_matrices, nc_blocks, egn)` — new helper

For each MT with LTY=0 NC blocks, populates:

- `cov_matrices[(mt, mt)]` = self-cov, **replacing** any explicit NI
  placeholder (mirrors Fortran `akxy[derived,derived,k]=0`).
- `cov_matrices[(mt, ref_j)]` for each ref_j > mt = `c_j *
  Cov(ref_j, ref_j)` (the writer's pair iterator only emits pairs
  with `mt2 > mt`, matching Fortran convention).

Simplifying assumption: **input cross-MT covariances are zero** —
true for U-238 JENDL because every referenced MT has NC=0/NI=1
self-cov only. The full Fortran double sum
`Σ_{iy,iyp} akxy[iy,ix,k] * akxy[iyp,ixp,kp] * Cov_in(iy,iyp)`
collapses to `Σ_i c_i² Cov(ref_i, ref_i)` (self) and `c_j Cov(ref_j,
ref_j)` (cross). Group mask uses geometric center
`√(egn[g] * egn[g+1])` against `[E1, E2]` — for full-range NC blocks
the mask is all-true.

### v1 limitation

Cross-pairs where **both** endpoints are NC-derived (e.g. T15
Cov(2, 4): `-Σ_{iy ∈ {51..77,91}} Cov(iy, iy)`) are not computed.
Cost: ~40 lines short of reference per such pair. T15 has exactly
one (Cov(2, 4)). Filed as follow-up.

## Verification

### RED→GREEN test

`test/validation/test_errorr_nc_expansion.jl`:

| | Pre-fix | Post-fix |
|-|---------|----------|
| `@test` pass | 3 | 7 |
| `@test` fail | 4 | 0 |

| Quantity | Pre-fix | Post-fix | Reference |
|----------|---------|----------|-----------|
| T15 MF33 MT=2 lines | 106 | **1429** | 1400 |
| T15 MF33 MT=4 lines | 103 | **1099** | 1106 |
| T15 MF33 total      | 1556 | **3875** | 5655 |
| T15 tape26 total    | 1859 | **4178** | 5958 |

### Regression

- **T04** tape23 NUMERIC_PASS 81/82 (unchanged from master baseline —
  the line-14 ±1 ULP at 8 sigfigs was already present pre-fix; HANDOFF
  Phase 47's "BIT_IDENTICAL" claim was inaccurate).
- T04 tape24 NUMERIC_PASS 56/74 (unchanged).
- T04 tape25 DIFFS 107/119 (unchanged).
- Phase-46 `test_errorr_gendf_readback.jl`: 38/38 (unchanged).
- Phase-47 `test_errorr_mf33_sparse.jl`: 36/36 after one stale
  cap (`total_lines < 4000` → `< 5000`) was relaxed to permit NC's
  legitimate +2300-line growth.

## Files changed

- `src/orchestration/modules/errorr.jl`:
  - Read loop captures LTY=0 NC blocks into `nc_blocks`.
  - New `NCSubSubsection` struct.
  - New `_read_nc_subsection(io)` parser (CONT + LIST, two records).
  - New `_expand_nc_blocks!` post-processing pass (≈70 LOC).
  - One `_expand_nc_blocks!` call between read loop and writer.
- `test/validation/test_errorr_nc_expansion.jl`: new RED→GREEN test
  (~110 LOC).
- `test/validation/test_errorr_mf33_sparse.jl`: relaxed
  `total_lines < 4000` → `< 5000` to accommodate NC growth.

## Traps (NEW)

**Trap (NC sub-subsection = CONT + LIST)**: each NC entry counted in
the sub-section header's `N1` field is **two** ENDF records — a
1-line CONT giving LTY, then a LIST with E1/E2 in C1/C2 and a body of
NPL = 2·NCI floats interleaved as (c_j, XMTI_j). A single
`read_list` per entry mis-aligns the file pointer when LTY ≠ "no
NC".

**Trap (NC LIST head fields)**: E1 and E2 live in `lst.C1` and
`lst.C2` (LIST CONT header positions), **not** in `lst.data[1:2]`.
The body is purely the (c_j, XMTI_j) pairs. Reading E1/E2 from `data`
makes the parser overrun by one pair and silently fail.

**Trap (NC overrides NI on the diagonal)**: Fortran `akxy[derived,
derived, k] = 0` (errorr.f90:1475) means the NC formula **replaces**
any explicit NI self-covariance for an NC-derived MT. Treating NC as
additive (NI + NC) double-counts. The Julia post-pass uses the same
replace-not-add semantics for the (mt, mt) diagonal entry.

## Follow-ups

- **Double-NC-derived cross-pairs** — T15 Cov(2, 4) = `-Σ Cov(iy,
  iy)` over iy in MT2_refs ∩ MT4_refs. ~40 lines per such pair; one
  pair on T15. Expand `_expand_nc_blocks!` to do the full double
  sum once a test demands it.
- **LTY=1/2/3 (standards / ratio)** — T04 U-235 MT=18 has 4
  LTY=3 cross-material NC sub-subsections currently emitted as
  zero stubs (matches existing reference behavior). Port Fortran
  `stand` (errorr.f90 around line 7800) when a test needs real
  derived-from-standards covariance values.
- **NJOY.jl-f8k** — the −1780-line residual gap to reference T15
  tape26 5958 is content drift on the *other* 34 MTs' self-cov
  (covcal/resprp pipeline at errorr.f90:7170-7188), independent of
  NC expansion.

## Beads status note (2026-04-20)

Beads DB is broken on this machine — orphaned IDs referenced here
(`NJOY.jl-km1`, `NJOY.jl-f8k`) no longer resolve via `bd show`. The
authoritative open-work list is `HANDOFF.md` "Open Work" section.
Mapping for this worklog:

- `NJOY.jl-km1` → "NC-block expansion v1" (closed here) + "NC-block
  expansion v2 (double-NC-derived + LTY≥1)" (P1 in Open Work).
- `NJOY.jl-f8k` → "Covcal content drift" (P1 in Open Work).
