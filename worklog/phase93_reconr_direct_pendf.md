# Phase 93 — Direct RECONR PENDF serialization

Date: 2026-07-15

Beads: `NJOY_jl-gkw` and `NJOY_jl-rhm` closed; residuals split to
`NJOY_jl-fod` (T83) and `NJOY_jl-c33` (T85).

## Outcome

- **T81: DIFFS → BIT_IDENTICAL**, 45,372/45,372 lines.
- **T84: DIFFS → BIT_IDENTICAL**, 1,115/1,115 lines.
- **T85: DIFFS → NUMERIC_PASS at `1e-5`**, 7 differing lines remain.
- T83 remains DIFFS, but its three-record structural offset and MT152 channel
  packing are fixed; the first residual is now a `csunr2` value difference.
- Full 86-test sweep: **11 BIT_IDENTICAL, 5 NUMERIC_PASS, 68 DIFFS,
  1 NO_REFERENCE, 1 TIMEOUT, 0 CRASH**. This is +2 BI, +1 numeric, and -3
  DIFFS from Phase 92 with all prior BI tests preserved.

## Oracle-first diagnosis

The focused T84 regression compared every byte of `referenceTape100`. Before
the fix Julia produced 1,112 records against 1,115 and first differed at MF1
line 3. The official T81 harness later isolated its only post-serializer
residual at tape30 line 50: MF2 scattering radius `0.0` versus `0.71`.

Ground truth:

- `reconr.f90:419-453`: `ruina` turns `ncards=0` into one blank Hollerith
  description card.
- `reconr.f90:5028-5174`: full ENDF-6 MF1/MT451 header, dictionary, and
  per-section NC values.
- `reconr.f90:5230-5355`: redundant sections receive `L2=99`; ordinary
  sections retain source L2/LR.
- `endf.f90:1064-1170,1280-1355`: positive `npend` carries sequence numbers
  across copied sections/files, while computed redundant sections still use
  `SEND=99999` and restart the counter.
- `samm.f90:664,1093-1127`: RML overwrites the default scattering radius with
  each elastic-channel `RDEFF`; the final value is written to MF2.
- `reconr.f90:1628-1735`: MT152 keeps total, elastic, fission, and capture
  independent; the sixth per-energy payload value copies total.

Two independent read-only agents audited the Fortran storage/write order and
the Julia/oracle record layouts. The main agent was the sole Julia runner.

## Implementation

- `ReconrParams` now retains whether `npend` was positive before normalizing
  the tape unit.
- Direct single-material RECONR reads the source `Mf1HeaderInfo` and supplies it
  to the writer.
- The legacy body writer now emits the complete ENDF-6 MF1 header, canonical
  blank description, per-section dictionary counts, source ZA/AWR, faithful
  redundant L2 and source LR, and the correct coded/binary sequence behavior.
- RML MF2 output selects the final elastic-channel `RDEFF` as Fortran does.
- URR table builders no longer replace the independently accumulated total
  with `elastic+fission+capture`; MT152 serializes `total_copy` correctly.

## Verification

- Focused T84 byte oracle: 3/3 assertions pass, 1,115/1,115 exact records.
- Official reference tests:
  - T81: BIT_IDENTICAL 45,372/45,372.
  - T84: BIT_IDENTICAL 1,115/1,115.
  - T83: DIFFS 60,904/71,931; first residual MF2/MT152 line 103.
  - T85: NUMERIC_PASS 7,802/7,809 at `1e-5`; first residual line 82.
- Full sweep completed in 3,361.1 seconds (56.0 minutes), preserving all nine
  prior BI tests and the prior numeric-pass cohort. T17 remains the sole known
  300-second ERRORR timeout.

## Next work

1. `NJOY_jl-fod`: trace T83 `csunr2` accumulation/state and its downstream
   grid propagation.
2. `NJOY_jl-c33`: eliminate T85's six-node MT152 last-digit residual and reach
   the `1e-7` first-round floor.
3. Keep `NJOY_jl-1kf` separate: T45's final-assembly MF1 ownership still leaves
   tape40 at 7,185 versus 7,188 records.
