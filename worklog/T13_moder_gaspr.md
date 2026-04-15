# Phase 13 — moder extract-mode stub + gaspr dispatch stub

## Date: 2026-04-15

## Summary

Tape-unit plumbing batch targeting the 7 "SystemError opening tapeNN"
CRASHes (T12/T18/T24/T27/T34/T47/T65) from HANDOFF's crash taxonomy.

Reality check: **5 of 7 were already past CRASH** before this phase,
having moved to DIFFS as cascading wins from Phase 10–12 stubs (covr,
plotr). Only T24 still had a structural CRASH. T65 was a TIMEOUT, not
CRASH.

Fixes landed are general-purpose rather than T24-specific:

1. **moder extract-mode stub** — detect `1 ≤ |nin| ≤ 19` on card 1 and
   read the real input tape from card 3; cp to output. Matches Fortran
   moder's iopt=1 (extract material from ENDF tape) behaviour at the
   stub level. T24 was the test that surfaced this gap.
2. **gaspr dispatch stub** — same shape as purr/leapr/covr stubs. cp
   `npendf_in → npendf_out`. T12 and T24 both have `gaspr` calls;
   previously the pipeline warned "Module gaspr not yet implemented"
   and left the output tape MISSING.

## Verification

### Pre-session sweep (CRASH class → reality on this machine)

| Test | HANDOFF says | Reality now                                       |
|------|--------------|---------------------------------------------------|
| T12  | CRASH        | **DIFFS** (0/2, 41 s) — was already past CRASH    |
| T18  | CRASH        | **DIFFS** (0/4, 43 s) — was already past CRASH    |
| T24  | CRASH        | **CRASH** tape21 missing (moder iopt=1 unhandled) |
| T27  | CRASH        | **DIFFS** (0/6, 284 s)                            |
| T34  | CRASH        | **DIFFS** (0/2, 266 s)                            |
| T47  | CRASH        | **DIFFS** (0/8, 288 s)                            |
| T65  | CRASH        | **TIMEOUT** — errorr @ U-235 MF34 slow path       |

So: covr stub (P11) and plotr stub (P12) cascaded through to T12/T18/
T27/T34/T47 because those decks run covariance → plot-tape → viewr
chains.

### Post-fix assertions

`test/validation/test_phase13_moder_gaspr.jl` — 2/2 pass in 6m 23s.

- **T24**: `!occursin("tape21", run_error)` — moder extract-mode now
  writes tape21 correctly; crash (if any) has moved downstream to a
  different missing tape (tape26 from heatr's plot output — a separate
  gap tracked as `NJOY.jl-327`).
- **T12**: `run_ok` — gaspr stub produces tape22 (previously MISSING).

Phase 11 (3/3) and Phase 12 (8/8) regression suites still pass.

## T24 — not a clean CRASH → DIFFS

T24's deck is a 12-module stress test. After moder + gaspr fixes, it
reaches thermr (85 s) then crashes at viewr on tape26 (heatr's plot
output, never written). Additional gaps to close before T24 goes
DIFFS:

1. heatr's 4th card-1 int (plot tape unit) is silently dropped
   (`HeatrParams` has no `nplot` field). Bead `NJOY.jl-327`.
2. thermr 85 s runtime for Pu-239 @ emax=10 eV — perf, similar class
   as `NJOY.jl-326`.
3. acer iopt=7 (photoatomic) — still a stub; T24's second acer call
   wants this.

Minimum fix for T24 is the heatr plot-tape stub (bead NJOY.jl-327).
Deferred out of this phase to keep scope tight.

## Files

**Added:**

- `src/orchestration/modules/gaspr.jl` — stub `gaspr_module`.
- `test/validation/test_phase13_moder_gaspr.jl` — RED→GREEN suite.

**Modified:**

- `src/NJOY.jl` — include `gaspr.jl`.
- `src/orchestration/modules/moder.jl` — extract-mode branch +
  `_moder_extract_stub!` helper.
- `src/orchestration/pipeline.jl` — `:gaspr` dispatch.

## Sweep impact — estimated

Net change in CRASH bucket from this phase alone: T24 stays CRASH (crashes at a different tape later in the chain). T12 was already DIFFS.

The cascading Phase-10/11/12 effect on the 5 "already DIFFS" tests is
the real visible change — the HANDOFF's 16-CRASH taxonomy is now
≈9 CRASHes (pre-Phase-11 16, minus Phase-11 T05/T16, minus Phase-12
T06/T20/T43, minus cascading T12/T18/T27/T34/T47, plus shifts T15/T17
→ TIMEOUT). T24 and T65 remain structurally blocked.

## Recommendations (priority order, post-P13)

### Immediate (small, high-yield)

1. **heatr plot-tape stub** (bead `NJOY.jl-327`) — add `nplot::Int` to
   `HeatrParams` and touch the plot tape. Unblocks T24 (the next
   downstream crash).
2. **acer iopt=7 stub → touch** — currently the acer module already
   writes empty stub tapes for unsupported iopts, but T24's second
   acer call (iopt=7 photoatomic) is fine. This is "soft": no crash,
   just empty-file DIFFS.
3. **Remaining pre-P13 CRASHes** (non-cascading): T09 (leapr thin
   stub), T60 (Fe-nat IRDFF-II MT=1 missing in broadr).

### Medium-term

4. **Real plotr output** — biggest remaining dispatch singleton.
5. **Perf beads** `NJOY.jl-326` (broadr U-238), T65 errorr perf, T24
   thermr perf — profiling session needed.
6. **Full leapr MF7 emission** — Phase 10 outstanding.

## Commit

- `src/orchestration/modules/{moder,gaspr}.jl`, `src/NJOY.jl`,
  `src/orchestration/pipeline.jl`: extract stub + gaspr dispatch.
- `test/validation/test_phase13_moder_gaspr.jl`: new RED→GREEN suite.
- `HANDOFF.md`, this worklog.
