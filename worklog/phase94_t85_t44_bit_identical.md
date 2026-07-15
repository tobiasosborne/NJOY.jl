# Phase 94 — T85 and T44 bit-identical

Date: 2026-07-15

This phase used targeted oracle runs only, at the user's request. The latest
full-suite baseline remains the Phase 93 sweep in `reports/REFERENCE_SWEEP.md`.

## Outcomes

- **T85: NUMERIC_PASS → BIT_IDENTICAL**, tape50 **7,809/7,809** at `1e-9`.
- **T44: DIFFS → BIT_IDENTICAL**, tape35 **322/322** at `1e-9`.
- The existing T84 direct-RECONR regression remained byte-identical
  (**1,115/1,115**).

## T85 — L=2 unresolved phase shift

Bead: `NJOY_jl-c33` (closed).

The red oracle assertion covered tape50 lines 82-88. Total and elastic were
high in the final printed digit while fission and capture were exact.

An initial research report blamed the order in which `csunr2` adds each
L-state's potential scattering. Implementing that order produced **no output
change**, so it was reverted. Direct comparison with the exact Fortran helper
found the real defect in `_unfac`:

- Fortran `reconr.f90:4488-4490` uses `rho` for L=2 penetrability, but
  `rhoc*rhoc` in the L=2 phase-shift denominator.
- Julia used `rho^2` in both places.

Changing only the phase-shift denominator to `3 - rhoc^2` made the exact
seven-line assertion pass and moved the complete official T85 classifier to
BIT_IDENTICAL.

## T44 — BROADR preserves the incoming PENDF header

Bead: `NJOY_jl-cgy` (closed).

The independent candidate audit identified T44 as the cleanest second target:
Julia wrote 319 records against the oracle's 322. The red assertion covered
tape35 lines 1-13 and showed the entire difference:

- Fortran output had four ENDF-6 MF1/MT451 CONT records and RECONR's canonical
  blank Hollerith description.
- Julia synthesized two CONT records, hard-coded `LRP=3/LFI=1`, and read NWD
  from the wrong record, losing the blank description.

Per `broadr.f90:638-680`, BROADR now copies the incoming PENDF metadata records
verbatim, rewrites only the temperature and corrected dictionary size in the
version-dependent control record, and extracts descriptions after that record.
The focused 13-line assertion and the complete official T44 classifier pass
byte-for-byte.

## Validation

All Julia runs were serial and preceded by clearing the NJOY precompile cache.

- `test/validation/test_reconr_pendf_serialization.jl`: T85 lines 82-88 pass;
  T84 lines 1-1115 pass.
- `test/validation/reference_test.jl 85`: BIT_IDENTICAL 7,809/7,809.
- `test/validation/test_broadr_header_serialization.jl`: T44 lines 1-13 pass.
- `test/validation/reference_test.jl 44`: BIT_IDENTICAL 322/322.

No full sweep was run. Based on the Phase 93 sweep plus these two official
targeted conversions, the repository has 13 known bit-identical reference
tests; this is not presented as a refreshed full-suite baseline.

## Next

1. Re-run T45 when returning to `NJOY_jl-1kf`; the shared BROADR header fix
   should remove its three missing metadata records, but its TPID and numerical
   residuals remain separate until proven otherwise.
2. `NJOY_jl-sc5` / T80 is the strongest next numerical grind (47 residual
   lines); re-run T33 immediately after it because the same LEAPR `contin`
   path may clear its eight residual records too.
3. Continue the P1 T83 unresolved-resonance lane in `NJOY_jl-fod`.
