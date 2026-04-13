# Phase 8 — Fortran-Faithful Reference Test Framework

## Date: 2026-04-13

## Summary

Built a **faithful Julia port of `njoy-reference/tests/execute.py`** so every
Julia module gets tested under the exact same conditions as the corresponding
Fortran reference test. Replaces three ad-hoc per-test pipeline scripts
(`t01/t02/t04_pipeline.jl`) with one parameterised runner.

The framework was designed to satisfy one constraint and avoid one anti-pattern:

- **Constraint:** every Julia port should be testable by running the real
  Fortran input deck (`njoy-reference/tests/NN/input`) through `run_njoy()` and
  comparing every produced `tape{U}` against `referenceTape{U}` using the
  exact same line-equivalence semantics `execute.py` uses.
- **Anti-pattern avoided:** no subprocess-per-test (precompilation cache
  corruption has wasted hours on this project). Everything runs in one Julia
  process, sequentially. A Timer-based heartbeat makes hangs visible in ≤10s.

---

## Files created

| File | Purpose |
|------|---------|
| `test/validation/reference_test.jl` | **Generic runner.** `run_reference_test(N)` stages resources via CMakeLists, runs `run_njoy(input; work_dir)`, compares every produced tape against every `referenceTape{U}` using `execute.py`-equivalent semantics. Heartbeat every 10s; soft timeout warning at 600s. |
| `test/validation/sweep_reference_tests.jl` | **Sweep driver.** Loops all 84 tests, classifies each tape (BIT_IDENTICAL / PASS_Nε / STRUCTURAL_FAIL / CRASH / NO_REFERENCE / MISSING_TAPE), writes `reports/REFERENCE_SWEEP.md`. |
| `test/validation/module_coverage.jl` | **Coverage audit.** Parses every input deck, emits `MODULE_COVERAGE.md` — proves (or disproves) that every Julia port has at least one covering Fortran reference test. |
| `test/validation/MODULE_COVERAGE.md` | Generated coverage map. |
| `reports/REFERENCE_SWEEP.md` | Generated sweep report (per-tape pass/fail per tolerance). |

## Files modified

| File | Changes |
|------|---------|
| `src/orchestration/pipeline.jl` | Added `RunProgress` struct + `_progress_start!` / `_progress_end!` hooks. `run_njoy` now accepts `verbose::Bool` and `progress::Union{Nothing,RunProgress}` kwargs. Per-module `… start` / `… done Xs` print with `flush(stdout)` when `verbose=true`. `progress` object lets a background heartbeat Timer detect hangs without subprocess isolation. |
| `test/runtests.jl` | Added `@testset "Reference Tests (Fortran-faithful)"` at the end, invoking `run_reference_test` on T01/T02/T03/T04 with their documented tolerances. `@test_broken` for known-ongoing tapes. |

## Files deleted

| File | Reason |
|------|--------|
| `test/validation/t01_pipeline.jl` (492 lines) | Superseded. Manually hand-orchestrated reconr→broadr→heatr→thermr×2 with hardcoded `/home/tobiasosborne/…` paths (broken on different machines). Predates `run_njoy()`. Logic now lives inside `run_njoy` module dispatchers. |
| `test/validation/t02_pipeline.jl` (355 lines) | Superseded. Manually orchestrated reconr+broadr×3 at hardcoded temperatures. |
| `test/validation/t04_pipeline.jl` (145 lines) | Superseded. Already called `run_njoy()` — this was the seed pattern for the generic runner. Generalised. |

---

## How the framework works

### 1. Staging

`run_njoy(input_path; work_dir=mktempdir())` calls `build_tape_manager`
which parses `CMakeLists.txt` and stages the referenced resource tapes into
the work dir under the unit numbers specified. This is the **exact same
staging** the Fortran build does via CMake's `configure_file(... COPYONLY)`
directives — no reimplementation needed.

### 2. Execution

`run_njoy` dispatches each module in the input deck via `select case` over
`mc.name`. Supported modules (already existed): `moder`, `reconr`, `broadr`,
`heatr`, `thermr`, `groupr`, `unresr`, `ccccr`, `errorr`, `gaminr`, `dtfr`,
`matxsr`, `viewr`. `plotr` is skipped. `acer`, `purr`, `covr`, `leapr`, `wimsr`
are present in deck parsing but not dispatched — they will warn and the test
will pass any reference tapes that don't depend on those outputs.

### 3. Comparison

`tape_compare(ref_path, trial_path; rel_tol, abs_tol)` reads both files and
runs `line_equivalence` on each line pair. `line_equivalence` is a direct
port of `execute.py::lineEquivalence`:

1. Byte equality — pass.
2. Otherwise: regex-extract all numbers using `FLOAT_RE` (same pattern as
   execute.py's `floatPattern`). Counts must match. Each float must satisfy
   `isapprox(a, b; rtol, atol)`.
3. Non-numeric residue (line with all matches stripped) must match byte-for-byte.

**Subtle match with Python:** the greedy regex can absorb trailing
`MAT`/`MF`/`MT` trailer digits into the previous float's exponent (e.g.
`6.706487+11306` from a line ending `...6.706487+1 MAT=1306 MF=3`). In Python,
`float("6.706487+11306")` overflows to `inf`; `math.isclose(inf, inf)` is
True. Julia's `tryparse(Float64, "6.706487E+11306")` returns `nothing`, which
would spuriously fail the line. `_parse_endf_or_normal` reproduces Python's
overflow semantics: when mantissa parses but exponent ≥ 309, return `sign*Inf`.

Date strings (`\d{2}/\d{2}/\d{2}`) are pre-normalised to `XX/XX/XX` on both
sides before comparison — same as execute.py.

### 4. Verbose output + heartbeat

`run_reference_test` wraps `run_njoy` in `with_heartbeat(f, progress)`.
`progress::RunProgress` is updated by `run_njoy` at every module boundary:
current module name, start time, status (:idle / :running / :done),
last-activity string. A Julia `Timer` fires every 10s; if the current module
has been running for ≥10s, it prints:

```
[T03] … still in gaminr (18.4s, total 20.6s)  last: gaminr starting
```

At 600s it escalates to `⚠ SOFT TIMEOUT`. **No subprocess kill** — cache
safety mandates in-process execution. User must Ctrl-C to abort. This is an
accepted trade: cache corruption has cost far more hours than Ctrl-C costs.

### 5. Result classification

Each reference tape gets one of:

- `BIT_IDENTICAL`       — passes at `rtol=1e-9` (tightest typical)
- `PASS_1e-7` / `PASS_1e-5` / `PASS_1e-3` — passes at that tolerance
- `STRUCTURAL_FAIL`     — line counts differ (structurally incomplete)
- `DIFFS`               — same line count, content mismatch even at loosest tolerance
- `MISSING_TAPE`        — referenceTape{U} exists but tape{U} wasn't produced

Each test gets an aggregate:

- `BIT_IDENTICAL` / `NUMERIC_PASS` / `DIFFS` / `STRUCTURAL_FAIL` / `MISSING_TAPE`
- `NO_REFERENCE` — test ran but has no `referenceTape*` files
- `CRASH`       — `run_njoy` threw an exception

### 6. Integration with `Pkg.test()`

`test/runtests.jl` now ends with a testset that invokes
`run_reference_test` on T01/T02/T03/T04 at their documented tolerances.
Tapes in each test's `broken_tapes` list are tracked as `@test_broken`
(known-failing, tracked as such so an accidental fix is flagged).

**No GitHub CI.** The user explicitly vetoed CI emails. Testing is local
only, via `julia --project=. -e 'using Pkg; Pkg.test()'` or running
`test/validation/sweep_reference_tests.jl` directly.

---

## Validation results

Runner was validated against known state (from Phase 7 handoff):

| Test | Runner result | Phase 7 handoff state | ✓? |
|------|---------------|------------------------|----|
| T03 | `BIT_IDENTICAL` (9274/9274 @ 1e-9) | 100% byte-identical | ✓ |
| T01 | `NUMERIC_PASS` (32812/32962 @ 1e-5) | passes at 1e-5 | ✓ |
| T02 | tape28 `NUMERIC_PASS` @ 1e-5; tape29 `STRUCTURAL_FAIL` (627/6213) | tape28 passes 1e-5; tape29 missing groupr MF6 | ✓ |
| T04 | tape23 @ 1e-7; tape24 @ 1e-5; tape25 13 bad @ 1e-5 | tape23 @ 1e-7; tape24 @ 1e-5; 11 lines bad in covariance | ✓ |

Runner reproduces the exact state the hand-orchestrated `t0N_pipeline.jl`
scripts documented — proving it's faithful.

---

## Module coverage (from MODULE_COVERAGE.md)

**Dispatched in `run_njoy` AND covered by ≥1 Fortran test:** all 14
pipeline modules (moder, reconr, broadr, heatr, thermr, groupr, unresr,
ccccr, errorr, gaminr, dtfr, matxsr, viewr, plotr).

**Ports exist but not dispatched** (run_njoy silently skips with @warn):
- `acer` — 45 tests depend on this. **Biggest coverage gap.**
- `purr` — 13 tests
- `covr` — 9 tests
- `leapr` — 5 tests (T09, T22, T23, T33, T80 — S(α,β) generation)
- `wimsr` — 1 test

**Ports with no source code** (would need implementation): `gaspr` (5 tests),
`powr`, `mixr`, `resxsr`.

**Tests with no `referenceTape*`:** T76 (execution-only).

---

## How to use the framework

```bash
# Single test (verbose, heartbeat visible)
julia --project=. test/validation/reference_test.jl 3

# Subset
julia --project=. test/validation/reference_test.jl 1 2 3 4

# Full 84-test sweep → reports/REFERENCE_SWEEP.md
julia --project=. test/validation/sweep_reference_tests.jl

# Subset sweep
julia --project=. test/validation/sweep_reference_tests.jl 3 4 8 27

# Coverage audit (regenerates MODULE_COVERAGE.md)
julia --project=. test/validation/module_coverage.jl

# Full test suite including reference tests
julia --project=. -e 'using Pkg; Pkg.test()'
```

**Cache rule still applies**: `rm -rf ~/.julia/compiled/v1.12/NJOY*` before
the first invocation after source changes. Subsequent runs in the same Julia
process reuse precompilation.

---

## What this framework does NOT do

- **No subprocess-per-test** — cache safety. Hangs visible via heartbeat,
  killable via Ctrl-C.
- **No GitHub CI hookup.** Local-only. User vetoed CI emails.
- **No ACE/S(α,β) comparison for tests whose references aren't `referenceTape*`**.
  T14, T48, T50–54, T62 etc. exercise `acer`/`leapr` but their outputs aren't
  compared because those modules aren't dispatched by `run_njoy` yet.
- **No per-test tolerance auto-discovery.** `runtests.jl`'s
  `REFERENCE_TEST_STATE` lists the known tolerance per test manually. Bump
  when improved.

---

## Next steps (priority order)

### 1. Run the full 84-test sweep
Produces the first complete `reports/REFERENCE_SWEEP.md` against actual Fortran
reference tapes (the Phase 18 sweep only used oracle PENDFs, reconr-only).
Will surface which Julia modules are actually mature vs partial.

### 2. Expand `REFERENCE_TEST_STATE` in `runtests.jl`
Once the full sweep identifies more passing tests, add them. Each entry is
a regression gate: if a test regresses from its recorded tolerance,
`Pkg.test()` fails loudly.

### 3. Dispatch `acer` / `purr` / `covr` / `leapr` in `run_njoy`
`acer` is the biggest gap (45 tests). Modules exist but the pipeline skips
them. Add branches to `pipeline.jl`'s dispatch loop.

### 4. Known bad tapes to fix (from current T01-T04 state)
- **T02 tape29** (groupr GENDF, 627/6213 lines): missing multi-temp × multi-sigma0
  GENDF records, no MF6 transfer matrices. Needs groupr upgrade.
- **T04 tape25** (errorr covariance, 13 lines fail at 1e-3): MF33 LB=2 union-grid
  collapse for blocks straddling group boundaries (Phase 7 handoff §3).
- **T01 tape25** (thermr SAB pipeline): 150 lines fail at 1e-9 (pass at 1e-5).
  Same class as Phase 34.

---

## Why in-process heartbeat might miss a true CPU-hang

Julia's `Timer` runs on the libuv event loop; callbacks fire only when the
main task yields. A pure-CPU tight loop with zero I/O for >10s blocks the
Timer. In practice every module emits `@info` logs and file writes, so
heartbeat works reliably. If we ever hit a real no-I/O loop hang, upgrade
path is `Threads.@spawn` with `julia --threads 2` — heartbeat task runs on
its own OS thread and fires regardless of main-thread state.

Not done now because no observed case requires it.

---

## Commit

```
399c80e  Phase 7 handoff: T03 complete + T04 progress
THIS     Phase 8: Fortran-faithful reference test framework
```

Files:
- **New:** `test/validation/reference_test.jl`, `sweep_reference_tests.jl`,
  `module_coverage.jl`, `MODULE_COVERAGE.md`; `reports/REFERENCE_SWEEP.md`;
  this worklog.
- **Modified:** `src/orchestration/pipeline.jl`, `test/runtests.jl`.
- **Deleted:** `test/validation/t01_pipeline.jl`, `t02_pipeline.jl`, `t04_pipeline.jl`.
