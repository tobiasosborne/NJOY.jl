# Phases 64-65 — powr lib=1 multi-T + gamff self-shielding

**Date:** 2026-05-03 (continuation of Phase 60-63 powr lib=1 port)
**Goal:** Round out powr lib=1 by porting the multi-temperature loop and
the `gamff` self-shielding f-factor accumulation. Bit-identical to Fortran
NJOY's tape50 on a Fe-56 multi-T multi-σ₀ oracle.

## Outcome

| Phase | Test                                            | Status                                |
|-------|-------------------------------------------------|----------------------------------------|
| 64+65 | Fe-56 ntemp=3, nsigz=4, GAM-I 68 group          | **BIT-IDENTICAL on first run** (177/177 lines, 14337/14337 bytes) |

All five powr lib=1 standalone tests pass clean:
- carbon abs-only (Phase 60), carbon + MF6 elastic (Phase 61), carbon + MF6
  elastic + inelastic (Phase 62), U-235 fission (Phase 63), Fe-56 multi-T
  multi-σ₀ self-shielding (Phase 64-65, this commit).

T01 regression-clean at NUMERIC_PASS 32812/32962 @ 1e-5 (unchanged).

## Why Phases 64 + 65 are coupled

The decision flow inside `fast()` (powr.f90:528-641) shows that the entire
single-temperature output path (header + absorption + fission/nu/chi + all
matrix blocks) is unaffected by `ntp`. **Multi-temperature output enters
only at lines 614-641**, gated on `nff != 0 .and. iwr != 0`. `nff = 1` is
hardcoded at line 414; `iwr` is computed inside `gamll` and is non-zero
only when `nsigz > 1 .and. iff = 1`. The block emits self-shielding
f-factors as a 3D array [nsigz, ntp, ngmax] — **the f-factors are what
varies per-temperature, not the cross sections**.

So splitting Phase 64 (multi-temp loop) from Phase 65 (`gamff`) was an
artifact of the original worklog grouping; in practice the smallest
falsifiable unit is the combined feature, exercised by an oracle with both
ntemp > 1 AND nsigz > 1 AND iff = 1.

## What landed

### `src/orchestration/modules/powr.jl` (~460 LOC delta, +376/-84)

Three substantive changes:

1. **`_powr_read_gendf_for_fast` rewrite** — now walks ALL temperatures on
   the GENDF (was: filtered by `temp ≈ rtemp`), captures the FULL per-σ₀
   XS vector per (mt, ig) (was: dropped silently — see "Latent bug" below),
   and returns `tmpr[]`, `sigz[]`, `rtemp_idx`. Mirrors gamll's discovery
   order: rtemp must be the FIRST temp on the tape (else error per
   powr.f90:847). MF=6 is still single-temp — gamff doesn't touch MF=6 and
   gamxs's `skiprz(-2) + findf(matd,1,451)` epilogue at powr.f90:1395-1397
   re-positions to the reference temp.

2. **`_powr_pack_ffactors`** (new, ~110 LOC) — direct port of `gamff`
   (powr.f90:1401-1542). Per-temperature absorption (MT 102-107) + fission
   (MT 18..21, 38) accumulation in a flat `[1..nsigz, 1..ntp, 1..ngnd]`
   column-major array. Final pass converts raw ratios to log-transformed
   f-factors `log(σ_z) - 2*log(ratio)` and computes `iwr = ngmax - iglo + 1`.
   Fortran's `mf3mt1` sentinel is collapsed: rebuild per-temp cflux from
   each temp's MF=3/MT=1 if present, else fall back to the gamxs reference
   cflux (matches the case where MT=1 is absent and the sentinel keeps the
   inverted ref cflux).

3. **`_powr_write_ffactor_block`** (new, ~25 LOC) — emits the powr.f90:614-641
   block: `(4i6) nsigz ntp misc jwf` (misc=1 hardcoded), sigz vector, tmpr
   vector, abs f-factors over groups [iglo, ngmax] × ntp × nsigz, fission
   f-factors (only when jwf != 0), and the dummy `iwr*10`-element zero
   amisc array.

The orchestrator `_powr_fast_one_material` is updated to:
- Use `g.mf3[g.rtemp_idx]` for the existing single-temp paths.
- Extract `xs_vec[1]` (= σ₀=∞ = the reference σ₀) for the absorption /
  fission accumulators.
- Compute the gamff path BEFORE writing the (6i10) header line, so the
  header `iwr` field uses the gamff-computed `ngmax - iglo + 1` value
  (matches Fortran's tape50 output for our Fe-56 oracle: iwr=43).

### `test/validation/test_powr_lib1_selfshielding.jl` (new, 130 LOC)

12-assertion bit-identical test against the Fe-56 oracle. Mirror of the
Phase 60 test structure plus assertions on the new params fields (`nsgz=4`).
Includes the regen recipe in the footer.

### `test/validation/oracle_cache/powr_lib1_selfshielding/` (gitignored)

Fe-56 JEFF-3.3 multi-T multi-σ₀ oracle:
- `tape20`: Fe-56 ENDF (copied from `njoy-reference/tests/resources/n-026_Fe_056-JEFF3.3.endf`).
- `tape21..tape24`: PENDF chain (reconr, broadr×3, unresr×3, groupr×3).
- `tape50`: Fortran NJOY reference (177 lines, 14337 bytes).

Note: Fe-56 has no unresolved-resonance parameters; `unresr` reports
"copy as is to nout". This is expected — the σ₀-weighted flux that
produces per-σ₀ XS differences in the resolved range is computed by
`groupr`'s narrow-resonance treatment, not unresr.

## Latent bug fixed alongside

The previous `_powr_read_gendf_for_fast` had a hardcoded XS index
`xs_p0 = data[nl + 1]` (powr.jl:506 in the old code). This happened to
match the correct multi-σ₀ index `data[nl*nz + 1]` ONLY when nz=1.
Carbon (Phase 60-62) and U-235 (Phase 63) all used nz=1 GENDFs, so the
bug was latent. With Fe-56 nz=4, Julia was reading flux values where it
expected XS, producing flat ~0.247-0.323 values across the entire
absorption block (vs. Fortran's resonance-structured 0.07-0.7). Fixed by
deriving the index from Fortran gamxs: `data[nl*(jz + nz - 1) + 1]` for
σ-zero jz=1..nz at Legendre order l=1.

## Three Fortran quirks worth pinning

1. **`gamll` ntp final state**: after gamll, `ntp = 1 + count_additional`
   where `tmpr(1)` is the reference and `tmpr(2..ntp)` are the additional
   temps. Set by powr.f90:894-896 (`ntp=0; if (rtemp.gt.0) ntp=1; tmpr(1)=temp`)
   then incremented by powr.f90:849-855 for each additional temp found.
2. **`gamff` reads ALL ntp temps including the reference** — gamxs's
   `skiprz(-2) + findf(matd,1,451)` at powr.f90:1395-1397 re-positions back
   to the reference temp's MF=1/MT=451, so gamff's `itp=1..ntp` loop covers
   the full set, not just the additional temps. The reference-temp
   absorption/fission ratios `a(loc) / σ_ref` therefore round to exactly 1.0
   at jz=izref, giving the canonical `log(σ_z) - 2*log(1) = log(σ_z)`
   value at the izref/ref-temp slot.
3. **Fission contributes to BOTH locsf AND locabs in gamff** (powr.f90:1497-1498),
   not just locsf. This mirrors the gamxs reference-temp behavior and is
   needed to match the absorption f-factor when MT=18 is present.

## Phase B residual items

| Item                              | Where in Fortran      | Trigger                  | Status |
|-----------------------------------|-----------------------|---------------------------|--------|
| Multi-temperature loop            | gamll lines 845-855   | rtemp > 0, ntp > 1       | **DONE (Phase 64)** |
| `gamff` self-shielding factors    | gamff (lines 1401-1542)| iwr = 1                  | **DONE (Phase 65)** |
| Delayed-neutron spectra           | nfs > 0 + MF=5/MT=455 | nfs > 0                  | TODO (Phase 66 candidate) |
| iread=1 (matd<0) abs-direct read  | fast lines 533-538    | matd < 0                 | TODO (low priority — dead in practice) |

Lib=2 (LIBRAR thermal) and lib=3 (CLIB cpm) remain at Phase A scaffold.

## Test commands

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/test_powr_lib1_standalone.jl     # Phase 60
julia --project=. test/validation/test_powr_lib1_elastic.jl        # Phase 61
julia --project=. test/validation/test_powr_lib1_matrices.jl       # Phase 62
julia --project=. test/validation/test_powr_lib1_fission.jl        # Phase 63
julia --project=. test/validation/test_powr_lib1_selfshielding.jl  # Phase 64-65 (NEW)
julia --project=. test/validation/reference_test.jl 1              # T01 regression
```

## Status: powr lib=1 ~95% covered

After Phase 64-65, powr lib=1 handles the full carbon (abs-only, +MF6
matrices), U-235 (fission), and Fe-56 (multi-T multi-σ₀ self-shielding)
code paths bit-identical. Remaining: delayed-neutron spectra (single phase,
needs an oracle with `nfs > 0`) and the (effectively dead) matd<0 corner.
lib=2 and lib=3 await dedicated phases.
