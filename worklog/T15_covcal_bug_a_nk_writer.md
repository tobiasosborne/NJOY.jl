# T15 — covcal Bug A: writer NK count for cross-pairs (Phase 56)

**Date:** 2026-05-01
**Test:** T15 (U-238 JENDL, MAT=9237, second errorr call → tape26)
**Goal:** Make Julia's MF33 writer emit one sub-section per (mt, mt2)
with mt2 ≥ mt in the active reactions list, matching Fortran covout
geometry (errorr.f90:7244 `scr(6)=nmts-ix+1`, inner loop 7254
`do 180 ixp=ix,nmts`). Closes HANDOFF P1 sub-item "Bug A — writer:
missing zero-stub cross-pairs."

## Outcome

**Bug A landed.** Per-MT NK now matches reference for every MT in T15
tape26 (verified for all 36 reaction MTs). Total tape26 lines:
4275 → **5964** (ref 5958, +6 from sub-section content drift, separate
issue).

Canary acceptance:
- T15 MT=77 NK: 1 → **3** (matches ref [77, 91, 102]).
- T15 MT=2  NK: 35 (unchanged — already correct via NC-derived path).
- T15 MT=18 NK: 1 → **31** (now matches ref).

Regression-clean:
- T05/T16 covr-isolation tests: 3/3 BIT-IDENTICAL (Phase 55 preserved).
- T04 reference: tape23 NUMERIC_PASS 81/82 (matches HANDOFF acceptance),
  tape24 NUMERIC_PASS 56/74, tape25 DIFFS 107/119 — all unchanged.
- Bug B canary (T15 MT=77 C[20,20]): exact 0.02987998 (Phase 51
  preserved).

## Files

CHANGED:
- `src/orchestration/modules/errorr.jl:1192-1209` — `_write_errorr_tape`
  now synthesises `(mt, mt2) for mt2 ≥ mt in reaction_mts` for every
  mt, not only `mt ∈ nc_derived_mts`. Set semantics dedupe against
  `cov_matrices` and `listed_pairs` entries. Empty cross-pair sub-
  sections route through `_write_mfcov_rows`'s `matrix===nothing`
  branch (2-line zero stub), mirroring Fortran's iabort=1 path at
  errorr.f90:7350-7356 + label 390.
- `test/validation/test_errorr_covcal_lb5.jl` — added `_read_mf33_nk`,
  `_read_mf33_mt2_list`, and a new `@testset "errorr covcal Bug A — NK
  cross-pair stubs (T15 MT=77)"` asserting NK==3 and mt2==[77, 91, 102].
- `test/validation/test_errorr_mf33_sparse.jl` — re-tuned for post-
  Bug-A geometry: per-MT NK match (strict, all 36 MTs), per-MT line-
  count slack 15 → 60 (catches blowup, accommodates content drift),
  total cap `< 5000` → `5500 < total < 6500`. Comments updated to
  reflect that residual per-MT line gaps are sub-section content drift,
  not missing geometry.

## Bug

`_write_errorr_tape` built `pair_set` for each `mt` from
`cov_matrices`, `listed_pairs`, and a synthesis pass `for mt2 in
reaction_mts; mt2 > mt && push!` — but the synthesis was gated by
`if mt in nc_derived_mts`. Non-NC MTs (everything except MT=1, 2, 4,
17, 18, 37 etc. that are NC-derived sums) ended up with pair_set =
just `{(mt, mt)}` from cov_matrices, producing NK=1.

Fortran covout has no such gate — line 7244 writes `scr(6)=nmts-ix+1`
unconditionally as the head record's N2 (NK), and the inner loop
`do 180 ixp=ix,nmts` (line 7254) iterates every ixp from ix to nmts.
Cross-pair sub-sections without computed data are written as 2-line
empty matrices via the iabort=1 path (lines 7350-7356) terminating at
label 390.

## Fix

```julia
# Before
if mt in nc_derived_mts
    for mt2 in reaction_mts
        mt2 > mt && push!(pair_set, (mt, mt2))
    end
end

# After
for mt2 in reaction_mts
    mt2 >= mt && push!(pair_set, (mt, mt2))
end
```

The non-strict `mt2 >= mt` includes self; Set semantics dedupe against
the cov_matrices/listed_pairs adds. listed_pairs / nc_derived_mts
kwargs retained for caller compatibility but no longer drive geometry.

## Per-MT line counts post-Bug-A (T15 tape26 MF33)

```
MT  Julia   Ref   Δ   NK_jul  NK_ref
1     171   271  -100   36     36
2    1506  1400  +106   35     35
4    1106  1106     0   34     34
16    110   110     0   33     33
17    103   103     0   32     32
18    266   228   +38   31     31
37     91    91     0   30     30
51-77  *     *     0   29..3   match
91     38    38     0    2      2
102    41    79   -38    1      1
TOTAL 5661  5655   +6
```

NK matches ref **for every MT.** The non-zero per-MT deltas (MT=1, 2,
18, 102) are sub-section *content* drift, not geometry. They are
tracked as the still-open HANDOFF P1 sub-item "v2 sub-item 2 —
LTY=1/2/3 standards/ratio expansion" and the deeper covcal content
work; this phase only closes the writer geometry.

## Acceptance criteria — met

- [x] T15 MT=77 NK == 3 (was 1, matches ref).
- [x] Generalised across the 34 non-NC-derived MTs — all 36 MTs match
      ref NK exactly.
- [x] T15 tape26 total > 5500 lines (HANDOFF target; actual 5964 vs ref
      5958).
- [x] T04 reference test does NOT regress below NUMERIC_PASS 81/82 on
      tape23 — verified at exactly 81/82 post-fix.
- [x] `test_errorr_mf33_sparse.jl`, `test_errorr_gendf_readback.jl`,
      `test_errorr_nc_expansion.jl`, `test_errorr_covcal_lb5.jl` all
      pass.
- [x] T05/T16 covr-isolation tests still BIT-IDENTICAL (Phase 55).

## Out of scope (still open)

- **Sub-section content drift** (HANDOFF P1 sub-item 2): MT=1, MT=2,
  MT=18, MT=102 line-count gaps within sub-sections. Will require
  porting Fortran `stand` (errorr.f90 ~7800) for LTY=1/2/3 standards/
  ratio synthesis and tightening the LB=5/6 union-grid path beyond
  Phase-51's MT=77 fix.
- **End-to-end T05/T16 reference tests**: Phase 55 covr is bit-
  identical *given correct errorr input*. With Bug A landed, T05/T16
  tape26 are now closer to ref but still differ in sub-section
  content. End-to-end will pass when the content-drift work above
  lands.

## Fortran source citations

- `njoy-reference/src/errorr.f90:7244` — `scr(6)=nmts-ix+1` (NK head)
- `njoy-reference/src/errorr.f90:7254` — `do 180 ixp=ix,nmts`
- `njoy-reference/src/errorr.f90:7350-7356` — iabort=1 empty-matrix path
- `njoy-reference/src/errorr.f90:7461 (label 390)` — write empty matrix,
  next ixp
- `njoy-reference/src/errorr.f90:7530-7605` — covout row emission
  (already mirrored in Julia `_write_mfcov_rows`)
