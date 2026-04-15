# Phase 11 — Dispatch `covr` in `run_njoy`

## Date: 2026-04-15

## Summary

Wired `:covr` into the pipeline dispatcher per T10's priority-1
recommendation. The low-level covariance → correlation conversion and
boxer/CSV formatters live in `src/processing/covr.jl`; the missing piece
was the orchestration wrapper that (a) parses the covr card from the
input deck and (b) writes an output tape that downstream `viewr` can open.

Stub only — same shape as T10's `purr`/`leapr` stubs. Current `covr_module`
just `touch`es the output unit. Full plot-tape emission (the large block
of plotr-format commands that Fortran `covr.f90` writes for viewr to
render as PostScript) is deferred.

## Files

**Added:**

- `src/orchestration/modules/covr.jl` — stub `covr_module(tapes, params)`.
- `test/validation/test_phase11_covr_dispatch.jl` — RED→GREEN assertion.

**Modified:**

- `src/NJOY.jl` — `include("orchestration/modules/covr.jl")`.
- `src/orchestration/pipeline.jl` — `:covr` branch in the dispatch loop.
- `src/orchestration/input_parser.jl` — `CovrParams` struct + `parse_covr`.

## Behavior

`CovrParams` holds `nin, npend, nout` — the three tape unit numbers from
card 1. Other cards (nplt, reaction-pair lists for T16) are ignored by
the stub but live in `mc.raw_cards` for a future real implementation.

`covr_module` writes an empty file at `resolve(tapes, nout)` and
registers it. If `nout == 0` the stub is a no-op with a warning.

## Verification (RED → GREEN)

Ran the reference-test framework pre- and post-fix on T05 (fastest
covr-blocked test; 38 s vs T16's 20 min):

| State   | T05 outcome                                                    | tape34 state           |
|---------|----------------------------------------------------------------|------------------------|
| RED     | `CRASH: SystemError: opening file ".../tape34"`                | `MISSING`              |
| GREEN   | `0/3 tapes pass @ rtol=1e-09` (pipeline runs end-to-end)       | `STRUCTURAL_FAIL 0/78796` |

Test suite: `test/validation/test_phase11_covr_dispatch.jl` — 3/3 pass.

### Correction to T10 worklog

T10 §Recommendations item 1 claimed `covr` dispatch would unblock
**T05, T06, T16**. Investigation shows T06 is not a covr test at all —
its input deck is `plotr`, not `covr`. T06's crash is
`SystemError: opening file ".../tape31"` caused by the pipeline's
`plotr: skipped` branch (`src/orchestration/pipeline.jl:223`). T06 needs
a separate `plotr` dispatch (or for `viewr` to tolerate a missing input
tape), tracked as a follow-up.

## Sweep impact — estimated

Expected movement in next full sweep:

| Test | Pre-P11 | Expected Post-P11 | Note                                  |
|------|---------|-------------------|---------------------------------------|
| T05  | CRASH   | DIFFS             | verified above                        |
| T16  | CRASH   | DIFFS             | same structural fix; ~20 min runtime  |
| T06  | CRASH   | CRASH             | plotr, not covr — unchanged           |

So ~2 tests move CRASH → DIFFS. Full sweep not re-run (takes 2–3 h) —
next session or after the next batch of dispatches.

## Missing for bit-identical covr output

A real `covr_module` needs to reproduce Fortran `covr.f90`'s behavior:

1. Read MF33 (or MF34/MF35) covariance sections from the errorr tape at
   `params.nin`.
2. For each requested reaction pair (listed in later cards of the covr
   deck — absent from `CovrParams` today), build the correlation
   matrix via `covariance_to_correlation` (already ported).
3. Emit plotr-format plot commands: axes, heatmap, labels, contours.
   The format is the same as plotr's output — a series of ASCII
   command cards that viewr interprets.
4. Write those cards as tape N text so `viewr` can render PostScript.

The covariance math is done; the plot-tape serializer is the remaining
work. For T05 the expected tape34 is 78,796 lines, mostly coordinate
data for the correlation heatmap.

## Recommendations (priority order, post-Phase-11)

### Immediate

1. **T15/T17 `BoundsError`** — investigate empty-vector access after
   INT=0 fix let the reader past the TAB1 record (T10 §Recommendations
   item 2). One-liner in MF3 reader or broadr thinning is likely.

2. **T20 moder tape-unit-0 sentinel** — `moder_module` must skip card
   entries with unit 0. Simple guard.

3. **`plotr` dispatch (stub)** — same shape as this commit. Unblocks
   T06. `src/orchestration/pipeline.jl:223` currently `@info`s "skipped".
   A `plotr_module` stub that writes an empty tape at the output unit
   lets viewr proceed or fail with a clearer DIFFS classification.

### Medium-term

4. **Real `covr` plot-tape serializer** — unblocks bit-identical T05,
   T16 tape34/36/37.

5. **Real `leapr`** / **MF7/MT2 lattice reader** — biggest remaining
   modules (T10 §Missing for bit-identical).

## Framework health

Same pattern as Phase 10: one small dispatch + parser + test = one
verifiable CRASH → DIFFS transition per test. Per-fix cost stays low
because the reference-test framework does the oracle comparison for
free.

## Commit

- `src/orchestration/modules/covr.jl`: new.
- `src/NJOY.jl`, `src/orchestration/pipeline.jl`: include + dispatch.
- `src/orchestration/input_parser.jl`: `CovrParams`/`parse_covr`.
- `test/validation/test_phase11_covr_dispatch.jl`: new RED→GREEN suite.
- `HANDOFF.md`, this worklog.
