# Phase 95 — T33 LEAPR endout bit-identical

Date: 2026-07-15

This phase used targeted reference tests only. The latest full-suite baseline
remains Phase 93 in `reports/REFERENCE_SWEEP.md`.

## Outcome

- **T33: NUMERIC_PASS → BIT_IDENTICAL**.
  - tape24: **53,151/53,151** at `1e-9`
  - tape34: **53,151/53,151** at `1e-9`
- T09 remains BIT_IDENTICAL (**1,830/1,830**).
- T22 remains BIT_IDENTICAL (**4,636/4,636**).
- T80 improved from **91,406/91,453** to **91,452/91,453** and now passes
  the first-round `1e-7` floor.

## Oracle and diagnosis

Bead: `NJOY_jl-4lw` (closed).

Three independent read-only audits inspected the Fortran writer, Julia LEAPR
path, and fixed-column residuals while the parent ran all Julia work serially.
The Phase 93 sweep reported eight tolerance failures, but raw byte comparison
found ten original records: seven S(α,β) final-digit boundaries and three
tape34 AWR-format records. Two tiny S-value differences were hidden by the
harness absolute tolerance.

The exact red assertion covered:

- tape24 lines 5886, 9668, 15683
- tape34 lines 2, 76, 78, 2850, 4388, 6469, 13584

The full classifier then exposed one further record on line 53147 of each tape
after canonical extended `a11` formatting was restored; those lines received a
second red assertion before the final change.

## Root cause and fix

This was output serialization, not `contin`/`discre!` computational drift.

1. Fortran `endout` applies biased `sigfig(value,7,0)` (six digits below
   `small=1e-9`) before the SMIN cutoff and `a11`
   (`leapr.f90:3354-3417`). Julia emitted the raw transformed S value.
2. Julia forced `extended=false` for every LEAPR field. Canonical Fortran
   `a11` (`endf.f90:882-981`) permits the nine-digit fixed form, required for
   tape34 AWR `15.8575110`. Once S values are sigfig-rounded first, their
   trailing-zero rule naturally retains the reference scientific form.
3. Fortran sigfig-rounds both effective-temperature TAB1 axes before `a11`
   (`leapr.f90:3588-3615`). Julia wrote raw `tempf`, producing fixed forms such
   as `861.768224` instead of `8.617682+2`.

`src/processing/leapr_writer.jl` now reproduces those three writer semantics.
No LEAPR physics or accumulation loop changed.

## Validation

Every Julia invocation was serial and preceded by clearing the NJOY precompile
cache.

- Focused exact serialization regression: 4/4 pass.
- Official T33: both tapes BIT_IDENTICAL, 106,302 total records.
- Official T09 and T22: remain BIT_IDENTICAL.
- Official T80: one residual record remains, now NUMERIC_PASS at `1e-7`.

No full suite was run.

## Next

`NJOY_jl-sc5` has been corrected from the stale phonon-FP premise. T80's sole
remaining record is tape24 line 19: the MF7 B-list writes raw
`therm*beta_max=5.00000105`, while Fortran uses
`sigfig(therm*beta(nbeta),7,0)` and emits `5.000001+0`. Add the exact line-19
red assertion before applying that separate one-field fix.
