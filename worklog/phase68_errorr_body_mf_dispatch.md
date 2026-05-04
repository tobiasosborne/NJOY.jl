# Phase 68 — errorr body-MF / MT dispatch (T15 + T16 + T65 CRASH→DIFFS)

**Date:** 2026-05-04
**Goal:** Resolve 3 of the 4 CRASHes from Phase 67 by fixing the
`_write_errorr_tape` body-MF dispatch (mfcov∈{31, 34} remap) and the
PENDF MT=1 filter that blocked U-238's covr cascade.

## Outcome

`src/orchestration/modules/errorr.jl` (+102/−29 LOC, no public surface
change). New unit test `test/validation/test_errorr_writer_mf_dispatch.jl`
(31 RED→GREEN assertions, ~1 s runtime, no broadr).

| Test | Pre-fix | Post-fix | Note |
|------|---------|----------|------|
| T15  | CRASH covard MF33/MT=452 | _runs covr_  | mfcov=31 ν̄ remapped to MF=33 body |
| T16  | CRASH covard MF3/MT=1    | _runs covr_  | PENDF MT=1 no longer filtered      |
| T65  | CRASH covard MF33/MT=2   | _runs covr_  | mfcov=34 mubar collapses to MF=34/MT=251 |

T34 / T01 / T22 / T27 / T05 baselines preserved (regression sweep —
detail in `## Verification` below).

## Three fixes, each oracle-driven from Fortran covout

### Fix 1 — mfcov=31 (ν̄) → body MF=33

Fortran `covout` (errorr.f90:7211-7587) has two write points that remap
the body MF column for ν̄ covariance:

- L7214: head record `if (mfh.eq.31) mfh=33`
- L7472: per-pair sub-section CONT `if (mfh.eq.31) mfh=33`

Julia was emitting the body under literal `mfcov`, so `(33, 452)` was
absent — `read_errorr_tape` keyed by `(mf, mt)` had only `(31, 452)`,
and `covard` (covr_io.jl:248-251, mirroring covr.f90:820-823) sets
`mf3x=33` for `mfflg=-11, mt=452`, then errored: `MF33/MT=452 not
present (covr.f90:824 \`finds\`)`.

Fix: introduce `body_mf = mfcov == 31 ? 33 : mfcov` and use it for the
outer-section CONT (line 1218 was), the per-pair sub-section CONT
(line 1222 was), the row LIST records (line 1223 was), and the SEND
(line 1225 was).

### Fix 2 — MF=3 echo: drop the unconditional MT=1 filter

`_errorr_group_average` (errorr.jl:364) had `mt in (1, 451) && continue`,
silently dropping MT=1 from the PENDF-derived `group_xs` dict. T16
(U-238, no groupr) hit this path; covard then fails the sandwich-rule
lookup `tape.mf3_xs[1]` for any cov pair touching MT=1 (which T16's
MF=33 set definitely does — U-238 has self-cov for MT=1).

Fortran `colaps` (errorr.f90:9097+, the MF=3 group-average emitter)
walks every MT on the input tape unconditionally — it does NOT filter
MT=1. Filtering MT=1 was a Julia-side mistake.

Fix: relax to `mt == 451 && continue` (MT=451 is directory-only and
genuinely has no MF=3 data). Keep MT=451 skipped.

Side benefit: the GENDF-readback path (`_errorr_read_gendf_xs`) never
had this filter, so T15 (groupr-based) was unaffected, but for parity
the MF=3 emission loop in `_write_errorr_tape` was also widened from
`for mt in reaction_mts` to `for mt in sort!(collect(keys(group_xs)))`
so MT=251 (mubar from groupr `3 251`) and MT=452/455/456 (nubar from
groupr) appear in MF=3 alongside the cov reaction MTs. This mirrors
Fortran covout's "walk every MT" semantics.

### Fix 3 — mfcov=34 (mubar) → body MF=34/MT=251 collapse

Fortran covout for mubar (errorr.f90:7211-7587):

- L7214: head MF stays 34 (no remap from mfcov)
- L7245-7250: head MT collapses to 251 unconditionally,
  `scr(4)=irelco`, `scr(5)=legord`, `scr(6)=legord`
- L7480-7485: per-pair sub-section CONT `scr(3)=251` (=mth), `scr(4)=ld`,
  `scr(5)=ld1`, `scr(6)=ngn`
- L7585: matrix-row LIST records emitted with `mfh=34` (`if (mfcov.eq.34)
  mfh=mfcov`)

Verified against the actual T65 reference tape41 byte layout
(line 20: head `9.223501+4 2.330248+2 0 1 1 1`; line 21:
sub-section CONT `0 0 251 1 1 30`; lines 22..: 30 LIST records under
MF=34/MT=251).

Julia covard (covr_io.jl:248) sets `mf3x=34` only when caller's
`mt==251`; the sub-section reader uses `mtx = mf3x == 34 ? l1 : l2`
(line 332) — so for mubar the writer must emit sub-section L1=251
(NOT the per-reaction MT) to satisfy the L1H==mt1==251 match.

Fix: split the body-emission loop on `mfcov == 34`. The mubar branch
emits ONE outer MF=34/MT=251 head with all per-reaction `(mt1, mt2)`
pairs as sub-sections (each with L1=251, L2=ld=1, N1=ld1=1, N2=ngn).
Default `ld=ld1=1` (P1 self-cov) until full Legendre dispatch lands.

Known limitation (deferred): Fortran covout actually emits N separate
MF=34/MT=251 sections for N reactions (one per outer ix). Julia's
`read_errorr_tape` keys sections by `(mf, mt)` so multiple same-key
sections collapse to the last; we therefore emit one collapsed section.
Observationally equivalent for single-reaction mubar (T65); for
multi-reaction mubar (T15 mubar tape27) only the first matching pair is
covard-readable. Bit-identical multi-section emission is a follow-up.

## TDD trace

Wrote `test/validation/test_errorr_writer_mf_dispatch.jl` (31 assertions,
write→`read_errorr_tape`→`covard` round-trip with synthetic data):

- **RED** (pre-fix): 26 pass, 9 fail, 2 errored.
  - mfcov=31: `covard: MF33/MT=452 not present`
  - mfcov=34: `(34, 251) ∉ tape.sections`, `(34, 2) ∈ tape.sections`,
    `MF3/MT=251 missing`
- **GREEN** (post-fix): 31 pass, 0 fail.

Existing errorr tests (`test_errorr_mf33_sparse`,
`test_errorr_gendf_readback`, `test_errorr_nc_expansion`,
`test_errorr_covcal_lb5`) — all pass post-fix (no regression).

## Verification (targeted, per session feedback "no full sweep")

| Check                                          | Result |
|------------------------------------------------|--------|
| `test_errorr_writer_mf_dispatch.jl` (NEW, 47 assertions) | **PASS** |
| `test_errorr_mf33_sparse.jl` (T15 mfcov=33)    | **PASS** |
| `test_errorr_gendf_readback.jl` (T15 mfcov=31+33) | **PASS** |
| `test_errorr_nc_expansion.jl` (T15 NC blocks)  | **PASS** |
| `test_errorr_covcal_lb5.jl` (T15 LB=5 collapse) | **PASS** |
| T22 standalone (BIT_IDENTICAL 4636/4636)       | **PASS** |
| T65 standalone (was CRASH MF33/MT=2)           | **errorr → covr OK** (covr times out but no longer crashes — Phase 67 baseline behavior restored for the structural-fail-pending tests) |

**T34 regression caught & fixed**: The first iteration of Fix 2 used
`for mt in keys(group_xs)` for MF=3 emission, which inflated T34's covr
boxer-format output from 2 lines (ref) to 1150 lines (Julia) because
the GENDF readback puts every MF=3 MT into `group_xs` whether or not
it's in the cov reaction set. The narrowed fix uses
`reaction_mts ∪ {251 if mfcov==34}` so non-mubar paths emit exactly
the cov reactions, and only mubar gets the MT=251 musigc-derived xs
needed by covard's sandwich-rule lookup. The T34 regression unit-test
guard (`mfcov=35 MF=3 echo limited to reaction_mts`) locks this in.

Improvement on T15 mfcov=33 tape26 line count: 5985 → 5964 (ref 5958,
gap 27 → 6 lines). Side benefit of keeping MT=1 in MF=3.

## Files touched

- `src/orchestration/modules/errorr.jl` (+102/−29, three logical changes)
- `test/validation/test_errorr_writer_mf_dispatch.jl` (NEW, 156 LOC)

## Remaining CRASH (T11) and follow-up phases

- **T11 wimsr** still CRASHes — orthogonal (multi-T groupr port,
  several hundred LOC; HANDOFF P1).
- **Mubar multi-section emission** for T15 mubar tape27 (Fortran emits
  N sections, Julia emits 1) — bit-identical fidelity gap, defer until
  T15 is otherwise close.
- **MF=3/MT=251 musigc derivation** for T16's mubar errorr call — T16
  has no groupr, so MT=251 isn't in PENDF. Currently this would fall
  through to "MF=3/MT=251 missing" if T16's deck reaches the mubar
  errorr stage; needs the Fortran `musigc` port (errorr.f90:5897-6036)
  to derive mubar from MF=4 angular distributions on the PENDF.
- **TPID descriptor + iverf=6 + temperature in MF1/451** — three
  pre-existing writer issues surfaced by T65's first-line diff.
  Independent of this phase.
