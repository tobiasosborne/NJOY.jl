# Phase 97 — T83 RML total bit-identical

Date: 2026-08-22

This phase used targeted oracle and regression runs. The latest full-suite
baseline remains Phase 93 in `reports/REFERENCE_SWEEP.md`.

## Outcome

- **T83: DIFFS → BIT_IDENTICAL**, tape50 **71,931/71,931** at `1e-9`.
- The generated tape and oracle are also raw-byte identical, with matching
  SHA-256 `1eb70210a7ff8d3ace84ee1fef60de5a10eece84e510e1abd80ff550d61206ed`.
- T81, T84, and T85 remain BIT_IDENTICAL.
- T20 remains structurally incomplete, but its strict tape23 matches improve
  from the Phase 93 sweep's 1 to 1,039 of 28,549 compared records.

## Oracle-first re-scope

Bead: `NJOY_jl-fod` (closed).

The bead's original `csunr2`/MT152 premise was stale after Phase 94 changed
the L=2 phase-shift denominator from `rho^2` to Fortran's `rhoc^2`. A fresh
official run and three independent read-only audits established:

- MF2/MT152 lines 101-130 were already byte-identical.
- Trial and oracle had exactly 71,931 records and identical energy grids.
- The only residual was MF3/MT1 lines 134-11717, covering all 34,752 resolved
  RML points below the 5,000 eV unresolved-range boundary.
- Every Julia MT1 value was low by the separately exact MT107 contribution.

The exact red regression preserves lines 101-130 as an MT152 guard and names
tape50 line 134 as the first failing resolved MT1 record.

## Root cause and fix

Fortran `emerge` initializes `sn=0`, optionally reads the smooth MF3 value,
then unconditionally adds `res(1+itype)` before sigfig rounding and redundant
total accumulation (`reconr.f90:4791-4809,4832-4893`).

Julia's `merge_background_legacy` returned early whenever the smooth MF3
background was zero. The RML MT800/MT801 reaction addition occurred after
that guard, so those exact reaction sections were omitted only from MT1.

The guard now skips a zero smooth background only when no resonance reaction
exists for that MT. No formulas, grids, or writer paths changed.

## Validation

Every Julia invocation was serial and every run started from a cleared NJOY
precompile cache.

- Focused T83 regression: MT152 guard passed; line 134 failed red, then all
  3 assertions passed green.
- Official T83: BIT_IDENTICAL, 71,931/71,931.
- Raw full-tape `cmp`: exact.
- Canonical resonance cohort T01/T02/T08/T27/T34/T45/T46 completed with no
  status regression; T01 remains NUMERIC_PASS at 32,858/32,962 strict lines.
- RML cohort: T81 45,372/45,372, T84 1,115/1,115, and T85 7,809/7,809 remain
  BIT_IDENTICAL; T20's strict overlap improved.

No full sweep was run. Combining the Phase 93 sweep with targeted Phases
94-97 gives 16 known bit-identical reference tests.

## Next

Re-run and re-scope `NJOY_jl-1kf` from the current T45 output: Phase 94 already
closed its three-record MF1 structural gap, so the old bead premise is stale;
the current tape is 7,188/7,188 records with a blank-TPID and numerical grind.
