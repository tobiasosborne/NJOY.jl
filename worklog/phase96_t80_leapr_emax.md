# Phase 96 — T80 LEAPR B(4) bit-identical

Date: 2026-08-22

This phase used targeted oracle runs only. The latest full-suite baseline
remains Phase 93 in `reports/REFERENCE_SWEEP.md`.

## Outcome

- **T80: NUMERIC_PASS → BIT_IDENTICAL**, tape24 **91,453/91,453** at `1e-9`.
- T09 remains BIT_IDENTICAL (**1,830/1,830**).
- T22 remains BIT_IDENTICAL (**4,636/4,636**).
- T33 remains BIT_IDENTICAL on tape24 and tape34 (**53,151/53,151 each**).

## Oracle and diagnosis

Bead: `NJOY_jl-sc5` (closed).

The official red T80 run reproduced one failing record, tape24 line 19:

- Fortran: `5.000001+0`
- Julia: `5.00000105`

The exact red regression names that tape and line and compares the complete
80-character record byte-for-byte.

Direct source inspection and an independent read-only Fortran audit agreed on
the serialization chain. `leapr.f90:3291-3323` (`endout`) assigns
`scr(10)=sigfig(therm*beta(nbeta),7,0)` before `listio`; `endf.f90:170-205`
and `838-981` then route the already-rounded B array through `lineio` and
canonical `a11`. ENDF-6 Formats Manual Section 7.4 defines B(4) as the
principal-scatterer EMAX field.

Julia's comment described that behavior, but `_write_mf7_mt4` stored the raw
`therm*beta_max` product. This was writer drift, not phonon-loop or scattering
law accumulation drift.

## Fix

`src/processing/leapr_writer.jl` now applies
`round_sigfig(therm * beta_max, 7, 0)` before placing B(4) in the LIST payload.
The biased sigfig result triggers Fortran's canonical scientific `a11` form.

## Validation

Every Julia invocation was serial and preceded by clearing the NJOY precompile
cache.

- Exact T80 line-19 regression: red on the raw product, then 2/2 pass.
- Official T80: BIT_IDENTICAL, 91,453/91,453.
- Official T09: BIT_IDENTICAL, 1,830/1,830.
- Official T22: BIT_IDENTICAL, 4,636/4,636.
- Official T33: both tapes BIT_IDENTICAL, 53,151/53,151 each.

No full sweep was run. Combining the Phase 93 sweep with targeted Phases
94-96 gives 15 known bit-identical reference tests.

## Next

Continue the P1 `NJOY_jl-fod` T83 unresolved-resonance grind at tape50 line
103, then re-evaluate the remaining immediate metadata and thermal residual
beads against fresh targeted oracle runs.
