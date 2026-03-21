# NJOY.jl Session Handoff

## What is this project

NJOY.jl is a Julia port of NJOY2016 — the standard nuclear data processing system used worldwide for reactor physics, criticality safety, and radiation transport. The original is 119,613 lines of Fortran 90. Our Julia version is ~14,000 lines.

The goal: produce bit-compatible PENDF/ACE output that passes all 85 of NJOY's own reference test problems. The port must be idiomatic Julia — composable, differentiable, no global state — not a transliteration.

**Repo:** https://github.com/tobiasosborne/NJOY.jl (GPL-3.0)

## What exists right now

### Source code (49 files, ~14,000 lines)
All 23 NJOY2016 processing modules are implemented:

| Category | Modules | Location |
|----------|---------|----------|
| Foundation | ENDF I/O, types, constants, interpolation | src/endf/ |
| Resonances | SLBW, MLBW, Reich-Moore, SAMMY/LRF=7, Faddeeva | src/resonances/ |
| Processing | RECONR, BROADR, HEATR, THERMR, UNRESR, PURR, GROUPR, ERRORR, GASPR, LEAPR, GAMINR, MIXR, MODER, RESXSR, COVR | src/processing/ |
| Formats | ACER, CCCCR, MATXSR, WIMSR, DTFR, POWR | src/formats/ |
| Visualization | plotr/viewr replacement (ASCII + PostScript + Makie recipe) | src/visualization/ |

### Tests
- **Unit tests:** ~13,360 assertions in test/runtests.jl (all pass)
- **Integration tests:** test/integration_tests.jl (3 tests, 1 pre-existing failure: Sr-88 SAMMY)
- **Validation pipeline:** test/validation/ — full infrastructure for all 85 NJOY tests
- **Diagnostic runner:** test/validation/diagnose_test.jl — verbose per-test diagnostics

### Reference data (gitignored, available locally)
```bash
git clone https://github.com/njoy/NJOY2016.git ./njoy-reference
git clone https://github.com/njoy/ENDFtk.git ./endftk-reference
```
These contain the Fortran source, 85 test problems with reference tapes, and 52 ENDF evaluation files. **Already cloned locally — do not re-clone.**

**CRITICAL:** The 85-test suite has ENDF evaluations baked in (the .endf files referenced by each test's CMakeLists.txt → `njoy-reference/tests/resources/`). Always use these exact files, never substitute newer versions from NNDC — evaluation revisions change resonance parameters.

## What was done this session

### Bug Fix 1: `parse_endf_float` (COMPLETED)
**File:** `src/endf/io.jl:12-24`

Fixed two bugs:
1. **Spaces in exponents:** `"2.00300e+ 3"` → `filter(!isspace, s)` strips all whitespace, not just leading/trailing
2. **Compact notation via regex:** Replaced manual backward loop with `replace(t, r"([\d.])([+-])(\d)" => s"\1e\2\3")`

Tests added in `test/runtests.jl` for spaces, D/d notation, multi-digit exponents. All 13,360 unit tests pass.

### Bug Fix 2: Physics Constants CGS→SI (COMPLETED — THE BIG ONE)
**File:** `src/constants.jl`

**Root cause of 10,000x resonance cross-section errors.** The `PhysicsConstants` module used CGS units (ergs, grams, cm/s) but ALL formulas ported from NJOY2016 Fortran assume SI (Joules, kg, m/s). This made `cwaven_constant()` return 0.002197 instead of 0.2197 (100x too small), making `pifac = π/k²` 10,000x too large, inflating every resonance cross section by ~10,000x.

**What changed:**
- `ev`: 1.602e-12 (erg) → 1.602e-19 (J)
- `clight`: 3e10 (cm/s) → 3e8 (m/s)
- `amu`: derived in grams → 1.661e-27 (kg)
- `hbar`: derived in erg·s → 1.055e-34 (J·s)
- `finstri`: replaced formula with exact CODATA value 137.036

**Downstream fixes:**
- `src/processing/reconr_evaluator.jl:554`: Removed `clight/100` (was converting cm/s→m/s, no longer needed)
- `src/resonances/sammy.jl:66`: Same `clight/100` removal

**Verification:** Test 02 (U-234, MAT=1050) capture XS at 5.23 eV went from 6,537 barns → 0.65 barns (matches physics). Worst error against reference dropped from 8,127,229% → 713%.

### Infrastructure: Validation Pipeline Module Chaining (COMPLETED)
**Files:** `test/validation/test_executor.jl`, `input_parser.jl`, `njoy_test_runner.jl`, `run_all.jl`

The validation pipeline now chains modules per the NJOY input deck instead of only running `reconr()`. Supports: `reconr → broadr → heatr → thermr → gaspr → acer` with graceful degradation for `unresr`, `purr`, `groupr`, `errorr`, `covr` (pass-through).

Key changes:
- `ExecutionResult` now has single `pendf` field (was separate `reconr_output`/`broadr_output`)
- AWR extracted from ENDF file MF1/MT451 and passed through pipeline (was hardcoded 1.0)
- Added `NJOYTestCase.awr` field
- Added `reconr_to_pendf()`, `pendf_to_result()`, `ensure_pendf()` conversion helpers
- Added `execute_heatr`, `execute_thermr`, `execute_gaspr`, `execute_acer` functions
- Added input parsers: `ThermrParams`, `AcerParams`, `GasprParams`, `UnresrParams`, `PurrParams`

### Infrastructure: Diagnostic Runner (NEW)
**File:** `test/validation/diagnose_test.jl`

Verbose per-test diagnostic: shows parsed input deck params, step-by-step execution with timing, sample XS at key energies, detailed comparison against all reference tapes with top-10 worst points.

```bash
julia --project=. test/validation/diagnose_test.jl 2   # diagnose test 02
julia --project=. test/validation/diagnose_test.jl 7   # diagnose test 07
```

## Current state: HONEST ASSESSMENT (updated)

| What works | What doesn't |
|------------|-------------|
| All 13,360 unit tests pass | ~0/85 NJOY tests fully pass yet |
| Constants now SI (was CGS — 10,000x bug fixed) | MF3 backgrounds missing above resonance range |
| Pipeline chains reconr→broadr→heatr→thermr→gaspr→acer | Elastic/fission/capture go to zero above ~400 eV |
| AWR correctly extracted and passed through | Resonance grid misses some narrow peaks |
| Broadr sigma1 kernel verified correct vs Fortran | UNRESR/PURR/GROUPR not fully executed (pass-through) |
| Diagnostic runner gives detailed per-test output | ReferenceTape29-type (GROUPR multigroup) miscompared as PENDF |

## CRITICAL: What to do next (priority order)

### Step 1: Fix MF3 background merging above resonance range (HIGHEST PRIORITY)
**This is the #1 blocker for most tests.** After reconr, cross sections go to zero above ~400 eV because MF3 background data isn't being merged outside the resolved resonance range.

**Symptoms (from Test 02 diagnostic):**
- Elastic XS = 1e-8 above 400 eV (should be ~30 b from potential scattering)
- Fission XS = 0 in keV range (should be ~1.5 b from MF3/MT18)
- Capture XS = 0 in keV range (should be ~2 b from MF3/MT102)
- Total XS = 0 at 1000 eV (should be ~16 b)

**Root cause analysis (from deep agent investigation):**
The bug is in `src/processing/reconr_evaluator.jl` and `src/processing/reconr.jl`. Three interacting issues:

1. **Overlap region suppression** (`reconr_evaluator.jl:112-118`): The `merge_background!` function skips MT=2,18,102 backgrounds for energies in `[eresr, eresh)`. This is meant for the resolved/unresolved overlap region, but if `eresr ≈ eresh` (no unresolved range), it may suppress ALL backgrounds.

2. **Grid extension outside resonance range** (`reconr.jl:133-152`): Points outside `[eresl, eresh]` are added from MF3 breakpoints and decade multiples (1,2,5 per decade). These get zero resonance XS, then `merge_background!` is supposed to fill in MF3 data. But the overlap check in point #1 may prevent this.

3. **Inconsistent skip lists**: `reconr.jl:136` skips MT=(1,3,101) for grid extension, but `reconr_evaluator.jl:94` skips MT=(1,3,4,101,27). These should be consistent.

**How to debug:**
```bash
julia --project=. test/validation/diagnose_test.jl 2
```
Look at the RECONR output table — XS goes to zero at 1000 eV. Then check:
- What are `eresl`, `eresh`, `eresr` for Test 02's material?
- Does MF3/MT18 have data in the keV range?
- Is `merge_background!` being called for those energies?

**How to fix:** Compare against `njoy-reference/src/reconr.f90`, specifically:
- The `emerge` subroutine (handles background merging)
- The `panel` loop structure (iterates over energy panels)
- How backgrounds are added outside the resonance range

### Step 2: Fix resonance grid resolution
**Symptoms:** At ~191 eV in Test 02, there's a resonance the reference resolves (peak ~700 b) but our grid completely misses (we get 0.2 b). The adaptive grid isn't placing enough points near narrow resonances.

**Files:** `src/processing/adaptive_grid.jl`, `src/processing/reconr_grid.jl`

**Investigation notes from agent:** Issue in `reconr_grid.jl:149` — half-width nodes clamped to 7 sig-figs for resonances where `E_r >> half_width`. For a resonance at 191 eV with width 0.07 eV, the half-width offset rounds to the resonance energy itself, losing the flank node entirely.

### Step 3: Skip GROUPR reference tapes in comparator
**Quick fix:** ReferenceTape29-type files are GROUPR multigroup output, not PENDF pointwise data. The comparator reads them as PENDF and gets garbage. Either detect and skip them, or only compare against tapes matching the pattern of RECONR/BROADR output.

### Step 4: Check thermr.jl Bragg scattering units
**File:** `src/processing/thermr.jl:197-204`

The `build_bragg_data` function uses `amne = amassn * amu` and `econ = ev * 8 * amne / hbar²`. With the CGS→SI switch, `econ` changed by 1e4. If lattice parameters (`a`, `c`) are passed in Angstroms or some other unit, the formula may need a compensating factor. This only affects LEAPR/THERMR crystalline scattering (Tests 22, 23, 33, 80), not the main RECONR→BROADR chain.

### Step 5: Run the SAMMY reviewer (beads: NJOY.jl-eux)
Per the project rules, the SAMMY/LRF=7 implementation needs a skeptical reviewer (line-by-line diff against samm.f90).

### Step 6: Iterate on remaining tests
After Steps 1-3, run `julia --project=. test/validation/run_all.jl` or use the diagnostic runner on individual tests. Each failure will expose a real bug.

## MANDATORY RULES — READ BEFORE DOING ANYTHING

### The 3+1 Agent Workflow (NON-NEGOTIABLE)

For every module port, physics implementation, or significant architectural decision:

1. **DUAL PROPOSALS:** Spawn TWO independent subagents. They work WITHOUT seeing each other's output. Each produces a complete, working implementation.

2. **ORCHESTRATOR COMPARISON:** Compare both proposals against the Fortran reference code at njoy-reference/src/. Run tests. Select the better approach OR synthesize a hybrid. Document the reasoning.

3. **SKEPTICAL REVIEWER:** Spawn a reviewer agent that:
   - Diffs Julia code against Fortran line-by-line
   - Checks edge cases (zero widths, negative XS, energy boundaries)
   - Verifies physical invariants (sum rules, detailed balance)
   - Reviews for Julia idiom violations (type instability, allocations, global state)
   - Produces a REVIEW_REPORT.md with PASS/CONDITIONAL PASS/FAIL verdict

4. **NO CODE MERGED WITHOUT REVIEWER PASS.**

5. **DECISION LOG:** Every decision in reports/decisions/NNN-topic.md with: both proposals summarized, selection rationale, reviewer verdict.

### Critical Operational Rules

6. **NEVER run parallel Julia test processes.** Julia's precompilation cache (~/.julia/compiled/) is NOT safe for concurrent access. **Always:** kill all Julia processes, clear cache (`rm -rf ~/.julia/compiled/v1.12/NJOY*`), then run exactly ONE test.

7. **Use beads (bd) for issue tracking.** Run `bd ready` to see open work. Run `bd close <id> --reason "..."` when done.

8. **After every wave/phase, restate the 3+1 rules.** LLM attention drifts. Repeat the rules explicitly before starting new work.

9. **Read the Fortran before writing Julia.** The reference is at njoy-reference/src/. When the Julia disagrees with NJOY output, the Julia is wrong until proven otherwise.

10. **Use the bundled ENDF evaluations.** The test suite's .endf files are in `njoy-reference/tests/resources/`. Never use external ENDF files — evaluation revisions change resonance parameters and produce false failures.

11. **Write elegant, idiomatic Julia.** No Fortran transliterations. Use multiple dispatch, broadcasting, clear type annotations, functional patterns.

## How to run things

```bash
# Run unit tests (takes ~90 seconds, 13,360 assertions)
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e 'using Pkg; Pkg.test()'

# Run full validation pipeline (all 85 tests, ~5 minutes)
julia --project=. test/validation/run_all.jl

# Run validation on specific tests
julia --project=. test/validation/run_all.jl --tests 1,2,7,12,84

# Diagnose a single test (verbose output)
julia --project=. test/validation/diagnose_test.jl 2

# Quick smoke test
julia --project=. -e '
using NJOY
endf = "njoy-reference/tests/resources/n-001_H_002-ENDF8.0.endf"
result = reconr(endf; mat=128, err=0.001)
println("H-2: $(length(result.energies)) points, thermal=$(round(result.total[1], digits=2)) b")
'

# Check beads issues
bd ready
bd list
```

## Key reference files

| File | Purpose |
|------|---------|
| src/constants.jl | Physics constants — NOW SI (was CGS, caused 10,000x bug) |
| src/endf/io.jl | ENDF float parser (fixed: spaces, compact notation) |
| src/resonances/slbw.jl | SLBW formalism + cwaven_constant() |
| src/processing/reconr.jl | Main RECONR pipeline |
| src/processing/reconr_evaluator.jl | Background merging — CONTAINS KNOWN BUG (Step 1) |
| src/processing/broadr.jl | Sigma1 Doppler kernel (verified correct) |
| test/validation/test_executor.jl | Pipeline execution + module chaining |
| test/validation/diagnose_test.jl | Verbose per-test diagnostic runner |
| test/validation/input_parser.jl | NJOY input deck parser |
| test/validation/reference_comparator.jl | PENDF reference comparison |
| njoy-reference/src/reconr.f90 | RECONR Fortran (5747 lines) — THE authority |
| njoy-reference/src/broadr.f90 | BROADR Fortran (1947 lines) |
| njoy-reference/tests/ | 85 test problems with reference tapes |

## Test 02 detailed state (use as debugging reference)

**Material:** U-234 (ZA=92234, AWR=232.029, MAT=1050 in ENDF/B-IV notation)
**Chain:** moder → reconr → broadr(300,900,2100K) → moder → unresr → groupr → ccccr → moder → moder
**ENDF source:** `njoy-reference/tests/resources/t404` (tape20)
**Reference tapes:** referenceTape28 (PENDF via unresr), referenceTape29 (GROUPR — ignore)

**Current diagnostic results (after SI fix):**
- RECONR: 2249 pts, thermal total=594 b — reasonable
- BROADR: 1727 pts at T=300K, AWR=232.0 — runs correctly
- Comparison vs referenceTape28:
  - MT=1 total: worst 99.97% at 192 eV (missing resonance peak)
  - MT=2 elastic: 100% above 400 eV (MF3 background missing)
  - MT=18 fission: 100% in keV range (MF3 background missing)
  - MT=102 capture: worst 713% at 5.23 eV (much improved from 8,127,229%)

**21 resonances** at: -1.68, 5.19, 31.4, 46.4, 49.4, 78.3, 88.7, 95.3, 106.9, 112.1, 132.9, 145.9, 154.0, 179.0, 184.0, 191.0, 274.0, 295.0, 319.0, 357.0, 369.0 eV

## Lessons learned (save yourself time)

1. **Precompilation races are the enemy.** If tests fail with weird errors after agent work, it's almost certainly cache corruption. Clear and retry once.

2. **The ENDF float format is treacherous.** `1.234567+8` (no E), `1.23456-38` (2-digit exponent), spaces in exponents, blank fields — every edge case exists in real data. The parser must be bulletproof.

3. **CGS vs SI was a 10,000x bug.** `constants.jl` used CGS (ergs, grams, cm) while all Fortran formulas assume SI (Joules, kg, m). Fixed in this session. If you see cross sections that are orders of magnitude too high/low, check units first.

4. **Sum rule tests need `>=` not `==`.** After adding non-primary MTs to total (gaspr, inelastic channels), total > elastic + fission + capture. Use `total >= sum_of_parts - epsilon`.

5. **LRF=7 (SAMMY) is the hardest formalism.** Variable channel count, per-channel radii, eliminated capture channel, background R-matrix. The Fortran is 7169 lines for a reason.

6. **Reference tapes include ALL processing steps.** A "RECONR test" reference tape may actually include BROADR+HEATR+UNRESR effects if those modules appear in the input deck. Always check the input file to know what processing was applied.

7. **ReferenceTape29-type files are GROUPR output**, not PENDF. The comparator reads them as PENDF and gets millions of percent error. Ignore these until GROUPR comparison is implemented.

8. **Test 02 is MAT=1050 = U-234 (not Pu-238).** Old ENDF/B-IV MAT numbers are tape-sequential, not nuclide-derived. Always check ZA in the ENDF file header.

9. **The user has 64 threads available.** Use parallel agents for research. But NEVER run parallel Julia test processes (precomp cache corruption).

10. **The user wants elegant, idiomatic Julia.** No Fortran transliterations. Favor multiple dispatch, broadcasting, clean type annotations. File length limits are a soft guideline.
