# Phase 67 — Full 84-test sweep + errorr/covr `mfflg` fix + wimsr loud diagnosis

**Date:** 2026-05-03 (continuation of Phase 66)
**Goal:** Establish the post-Phase-66 reference-test baseline. Triage the
recently-touched modules (Phases 50-66: gaspr / mixr / resxsr / powr /
wimsr / covr / leapr / errorr) before the long sweep. Fix any
quick-to-resolve regressions surfaced. Ship the full 84-test report.

## Outcome

Full 84-test sweep (1500s/test timeout, 15s heartbeat): **97.2 min**.

| Status            | Count | Apr-18 baseline | Δ        |
|-------------------|-------|-----------------|----------|
| `BIT_IDENTICAL`   | 2     | 1               | **+1**   |
| `NUMERIC_PASS`    | 2     | 1               | **+1**   |
| `DIFFS`           | 75    | 48              | +27      |
| `STRUCTURAL_FAIL` | 0     | 0               | =        |
| `MISSING_TAPE`    | 0     | 17              | **−17**  |
| `NO_REFERENCE`    | 1     | 1               | =        |
| `CRASH`           | 4     | 16              | **−12**  |
| `TIMEOUT`         | 0     | 0               | =        |

`MISSING_TAPE` 17→0 and `CRASH` 16→4 are the headline numbers — the cumulative
landings from Phases 50–66 + this phase's two fixes drove every previously
crashing or tape-less test through the pipeline. T01 (NUMERIC_PASS @1e-5,
32812/32962 lines) and T22 (BIT_IDENTICAL 4636/4636) baselines preserved.
T80 added as second NUMERIC_PASS (Phase 57); T03 added as second
BIT_IDENTICAL (photoatomic ASCII passthrough).

## Probe sequence (pre-sweep triage)

Per session feedback (`feedback_eager_flush.md`): every long-running command
must produce timestamped, eagerly-flushed progress so a hang surfaces in
seconds. The harness already does this (heartbeat Timer with `flush(stdout)`
and 300s `InterruptException` schedule); **verified before any sweep launch.**

- **Probe A** (90s budget): `using NJOY` smoke. Precompile 8.4s, total 9.2s. ✓
- **Probe B** (240s budget): T01 alone via `reference_test.jl 1`. 73s end-to-end,
  matches HANDOFF baseline exactly. ✓
- **Probe C** (~10 min budget): 7-test batch covering reconr / leapr / RM / FP-class
  diff tests (T01/T08/T19/T22/T34/T45/T80). 7.7 min wallclock. T22 + T80 + T01
  match baselines; surfaced **T34 CRASH** (`covard: MF33/MT=18 not present`).
- **Probe C2** (~15 min budget): 5-test batch covering Phase 50-66 module surface
  (T05 covr, T11 wimsr, T13 gaspr, T16 covr-multi, T27 covr-pu). 14.8 min
  wallclock. Surfaced **T11 BoundsError** in wimsr and **T27 CRASH** (same
  family as T34).

## Two fixes landed (commit `ac2a00e`)

### Fix 1: `_write_errorr_tape` mfflg dispatch

**Bug**: `src/orchestration/modules/errorr.jl:1157` hardcoded `N1=-11` in
the MF1/MT451 head regardless of `mfcov`. Fortran has THREE distinct writer
paths (errorr.f90:5927 `-11`, 6098 `-12`, 7823–7824 `-11/-14`); Julia has one.
Result: covr's covard reads `mfflg = -11`, dispatches `mf3x = 33`, fails to
find MF=35 cov on what is actually an MF=35 tape → `MF33/MT=18 not present`.

**Fix** (5 lines):
```julia
mfflg = mfcov == 35 ? -12 : (mfcov == 40 ? -14 : -11)
_write_cont_line(io, za, awr, 5, 0, mfflg, 0, mat, 1, 451, seq); seq += 1
```

**Impact**: T34 + T27 CRASH → DIFFS (run end-to-end through covr+viewr).
Each saved ~225–248s of broadr+groupr re-run.

### Fix 2: `wimsr_extract_resint` Rule-6 loud error

**Bug**: `src/processing/wimsr_resint.jl:275` sliced `xs_data.tempr[1:ires]`
with no length check. T11 has `ires=3` from the deck but the GENDF tape
contains only 1 temperature → `BoundsError [1:3] on 1-element vector`.
Root cause is upstream: Julia groupr (`src/orchestration/modules/groupr.jl:31`)
uses only `params.temperatures[1]` regardless of how many broadr produced,
so multi-T decks (T11 Pu-238 at 300/900/2100K) never get a multi-T GENDF.

**Fix** (Rule 6 — replace BoundsError with self-explaining `error()` naming
the exact upstream callsite). Real fix is the multi-T groupr port (separate
phase — added to HANDOFF P1).

**Impact**: T11 opaque BoundsError → diagnosed CRASH naming the gap. Test
still fails (groupr port out of scope this session) but the failure
is now actionable from the report alone.

## Verified on 5-test mini-batch (572s) before sweep

| Test | Pre-fix              | Post-fix                                  |
|------|----------------------|-------------------------------------------|
| T01  | NUMERIC_PASS 32812/32962 | NUMERIC_PASS 32812/32962 (preserved) |
| T22  | BIT_IDENTICAL 4636/4636 | BIT_IDENTICAL 4636/4636 (preserved)   |
| T34  | CRASH covard MF33/MT=18 | DIFFS — runs through covr+viewr       |
| T11  | CRASH BoundsError [1:3] | CRASH self-explaining error (Rule 6)  |
| T27  | CRASH covard MF33/MT=18 | DIFFS all 6 tapes run                 |

## 4 remaining CRASHes — all in errorr/covr/wimsr handoff

Each is well-localized and independent of the others. Listed here for the
next phase's pickup, not to relitigate now.

| Test | Crash signature                                       | Root cause (one-line)                                                               |
|------|-------------------------------------------------------|--------------------------------------------------------------------------------------|
| T11  | `wimsr_extract_resint: GENDF has 1 temp but ires=3`   | Upstream — Julia groupr is single-temperature only.                                  |
| T15  | `covard: MF33/MT=452 not present`                     | errorr writes ν̄ cov under literal `mfcov=31`; Fortran convention writes under MF=33. |
| T16  | `covard: MF3/MT=1 missing for MAT=9237`               | errorr output tape lacks MT=1 (total) in **MF=3** (cross-section, not cov).          |
| T65  | `covard: MF33/MT=2 not present`                       | errorr writes mubar cov under literal `mfcov=34`; Fortran writes under MF=33 + special-cases MT=251. |

The errorr→covr "output MF dispatch" is the same family of bug (Fix 1) but
covers only mfcov=35/40 separation; mfcov=31 (ν̄) and mfcov=34 (mubar) need
a *body* dispatch on the cov sections, not just the tape mfflg sentinel.
Approximate scope: 30–60 LOC in `_write_errorr_tape` + verify with T15,
T16, T65, T34 standalone re-runs. T11 wimsr is a much larger groupr port
(several hundred LOC).

## Harness observations worth pinning

1. **Hard timeout is at-next-yield, not preemptive.** Tight CPU loops in
   broadr (T16 broadr ran 506s on its own) ignore the scheduled
   `InterruptException` until they yield. The 300s default in the sweep
   script *post-hoc* marks the test `:TIMEOUT` but doesn't kill it; the
   1500s override (this phase) lets long-broadr tests complete naturally.
   Without the override every U-238/U-235/Pu-239 broadr-bound test would
   appear as TIMEOUT in the report despite being correct.

2. **First-test JIT amortizes.** First test in a Julia process pays ~60s of
   JIT/load before any module runs (heartbeat correctly catches it as
   `… idle for 60.8s last: initializing`). Subsequent tests reuse JIT — Pu-239
   reconr that took 31s on T27 took 14s on T34 in the same session.

3. **Per-module progress is module-boundary-only.** Inside broadr / sigma1
   / leapr there is no `_progress_start!` instrumentation; the heartbeat
   shows "broadr for Xs last: broadr starting" for the entire module
   duration. Acceptable for sweep — the 1500s timeout catches genuine hangs
   — but worth knowing if a future debug session needs fine granularity.

## Test commands

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
# single-test reproducible smoke
julia --project=. test/validation/reference_test.jl 1

# small-batch with the bumped timeout
julia --project=. -e '
include("test/validation/sweep_reference_tests.jl")
run_sweep([1, 22, 34, 11, 27]; timeout_hard_sec=1500.0, heartbeat_sec=15.0)
'

# full 84-test sweep (97 min on this hardware)
julia --project=. -e '
include("test/validation/sweep_reference_tests.jl")
run_sweep(all_tests(); timeout_hard_sec=1500.0, heartbeat_sec=15.0)
'
```

## Session memory addendum

Saved `feedback_eager_flush.md` per-user demand: every long-running
command MUST eager-flush + heartbeat + aggressive timeout; smoke run before
any 90+ minute sweep. Filed as a hard rule for all future sessions —
silent multi-hour stares are a recurring failure mode and the cost is
non-negotiable.

## Next phase candidates (priority-ordered)

- **P1 covr/errorr body MF dispatch** — fix T15 + T65 + T16 covr-cascade
  crashes by remapping mfcov∈{31, 34} cov sections to MF=33 in the output
  tape body (mirror Fortran's `nmtcov` / `covout` convention). One pass,
  ~60 LOC, three tests CRASH→DIFFS.
- **P1 multi-temperature groupr** — port the outer-temperature loop in
  Fortran groupr.f90 so multi-T decks (T11 Pu-238 prototype) get the
  correct number of GENDF temperature blocks. Larger — several hundred
  LOC. Unblocks T11 + likely cleaner numerics for T19/T67-T74.
- **P2 covr volume reduction** — T05/T07/T13/T82 emit 2–9× more lines than
  reference. Phase 55 covr port is BIT_IDENT on isolation fixtures but
  full-pipeline produces too much. Likely interacts with errorr's missing
  LTY=1/2/3 standards/ratio (HANDOFF P1 sub-item 2).
- **P2 leapr Phase B** (T80 BIT_IDENTICAL) — `naint` gating to close 47
  residual lines (HANDOFF P2; covered in Phase 57 worklog).
- **P3 broadr U-238 perf** — T15/T17 ran 498s/944s; HANDOFF P2 says Fortran
  does it in ~30s. Profile + hot-loop devectorize for ~20× speedup.
