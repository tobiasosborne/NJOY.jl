# T38 — purr module wiring (Phase 1: structural match)

## Date: 2026-04-24

## Outcome

**T38 (Kr-83, purr MT152+MT153) → STRUCTURAL progress.** Output tape34
goes from 2849/3480 lines (no MT152/153) → 3642/3480 lines with MT152
and MT153 sections present in both DICT and MF2 body. First 196 lines
match; remaining line-count gap is pre-existing upstream issues, not
purr wiring.

**T21/T31 (URR reader mode-12 failure) → graceful fallback.** The
unresr URR reader fails for LRF=2 (mode 12) URR materials with a
BoundsError; purr_module catches the exception and falls back to
PENDF copy-through, so downstream modules (acer) don't crash.

**T10 (Rb-85, ntemp=3 × nsigz=7) → multi-block output.** Produces
23 energies × 3 temps × 7 sigma-zeros in 0.03s. No crashes.

**T22 (leapr BIT_IDENTICAL) — no regression.** 4636/4636 @ 1e-9.
**T02 (Pu-238 / ACE) — no regression.** tape28 NUMERIC_PASS 12519/13873.

## What's in scope

1. `src/processing/purr_writer.jl` (new, 340 LOC) — PENDF splicer
   + MT152/MT153 emitters matching Fortran purr.f90:316-576:
   - `_purr_mt152_nc` / `_purr_mt153_nc`: Fortran's exact
     `2+(NPL-1)÷6` DICT NC formulas (purr.f90:365-368).
   - `_write_mt152_purr`: HEAD(ZA,AWR,LSSF,0,0,INTUNR) + LIST
     (TEMZ,0,5,NSIGZ,NPL,NE) + [sigz] + per-energy (E, sigu 5×nsigz).
   - `_write_mt153_purr`: HEAD(ZA,AWR,IINEL,IABSO,INTUNR,NBIN) +
     LIST(TEMZ,0,LSSF,0,NPL,NE) + per-energy (E, tabl 5·nbin, heat nbin).
   - `_write_purr_pendf`: line-by-line state machine splice.
     DICT rewriting canonicalises order: always skip upstream MT152/153
     entries and insert a fresh MT152 then MT153 after MT151 (reconr
     currently seeds a bogus MT152 with NC=40; we replace it).

2. `src/orchestration/modules/purr.jl` (stub → real, ~120 LOC):
   - Reuses `_read_urr_for_unresr` and `_read_unresr_backgrounds` from
     unresr.jl.
   - Loops temperatures, builds `PurrBlock` with zero-stub probability
     tables (see below), writes output PENDF.
   - Exception-handler fallback to PENDF copy when URR reader fails.

3. `test/validation/test_purr_writer.jl` — 15 structural unit tests:
   - DICT NC formulas validated against T21 (15/264) and T38 (38/727)
     reference values.
   - MT152/MT153 record layout counts match Fortran.

## What's deferred (Phase 11)

**Numeric purr values.** `generate_ptable` in src/processing/purr.jl is
a "Proposal B" algorithm: it uses a 900·dmin window and per-sample
resonance convolution. On Kr-83 T38 it generated 23 billion pool
allocations and was still running after 10 minutes. For the wiring
pass we stubbed `_purr_generate_ptable` to return uniform 1/nbin
probabilities and zero cross sections; line counts match but values
are zero. Path to bit-identical values: port Fortran `unresx` (ladder
construction + infinite-dilution XS) and `unrest` (10000-history
Monte Carlo sampling) directly. Same shape as the T22→T80 leapr
contin port (filed as NJOY_jl-63f).

**iinel/iabso competition flags.** Reference T38 MT153 HEAD shows
`iinel=4, iabso=107` from the purr.f90:1165-1192 heuristic (multi-MT
inelastic → 4, single absorption MT → mtc). Julia currently defaults
to −1/−1; only the HEAD line differs — line counts unaffected.

**URR reader mode-12.** `_read_urr_subsection` in
src/orchestration/modules/unresr.jl:299 crashes on LRF=2 URR
evaluations (Fe-58 T21, Pu-240 T31 and many others). Blocks the
numeric MT152/153 path for those materials even after Phase 11.
Purr falls back gracefully; fixing properly requires extending the
URR reader to handle LRF=2. Out of purr scope.

**Pre-existing MF1 HEAD format divergence.** NJOY.jl's reconr emits
a compact 2-line MF1 HEAD (LRP=3 LFI=1 ... nmod=33); Fortran NJOY
emits 4 CONT lines + NWD hollerith description. Affects every PENDF
comparison. First diff on every purr test is this, not MT152/153.

**`_read_urr_for_unresr` grid subdivision drift.** Julia produces
39 URR energies for Kr-83 vs Fortran's 36 — accounts for the 162-line
Julia overshoot on T38. Pre-existing in unresr's `_build_unresr_grid`.

## RED → GREEN (T38)

### RED
`julia --project=. test/validation/reference_test.jl 38` pre-session:
```
T38 tape34 STRUCTURAL_FAIL 1/2849 lines — no MT152/153, 631-line deficit.
```

### GREEN (after wiring)
```
T38 tape34 STRUCTURAL_FAIL 196/3480 lines — MT152/153 emitted,
162-line overshoot from upstream URR grid drift.
```
Delta: +195 matching lines; reversed deficit into structural overshoot.
Full MT152+MT153 body now present in both DICT and MF2.

## Files touched

- `src/NJOY.jl` — added `include("processing/purr_writer.jl")`
  after `leapr_writer.jl` (both depend on params from input_parser).
- `src/orchestration/modules/purr.jl` — rewritten from 48-line
  stub to 126-line real module.
- `src/processing/purr_writer.jl` — new (344 LOC).
- `test/validation/test_purr_writer.jl` — new (60 LOC, 15 tests).

## Reference

- Fortran purr.f90: `njoy-reference/src/purr.f90:1-1400`
- Key subroutines:
  - main entry: `:61-679`
  - rdheat (MT301/302/318/402 on PENDF): `:1237-1289`
  - unresx (ladder + ∞-dilution XS): `:1291-1538`
  - unrest (MC sampling, 10000 histories): `:1540-1716`
  - DICT rewrite: `:316-400`
  - MT152 emission: `:415-470`
  - MT153 emission: `:472-576`
- Reference tapes:
  - `njoy-reference/tests/38/referenceTape34` (Kr-83, LSSF=0, 36 energies)
  - `njoy-reference/tests/21/referenceTape24` (Fe-58, LSSF=1, 13 energies)
- unresr writer template: `src/orchestration/modules/unresr.jl:539-785`
