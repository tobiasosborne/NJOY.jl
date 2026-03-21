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
- **Unit tests:** ~10,000+ assertions in test/runtests.jl
- **Integration tests:** test/integration_tests.jl (3 tests wired up)
- **Validation pipeline:** test/validation/ (infrastructure for all 85 NJOY tests)

### Documentation
- README.md, docs/api.md, docs/tutorial.md, docs/architecture.md

### Reports
- reports/RESEARCH_REPORT.md — deep analysis of NJOY2016 Fortran
- reports/decisions/001-011 — architectural decisions with rationale
- reports/REVIEW_WAVE{1-5}.md — skeptical reviewer reports
- reports/VALIDATION_PLAN.md, TEST_READINESS.md, EXTENSIONS_PLAN.md

### Reference data (gitignored, must be cloned)
```bash
git clone https://github.com/njoy/NJOY2016.git ./njoy-reference
git clone https://github.com/njoy/ENDFtk.git ./endftk-reference
```
These contain the Fortran source, 85 test problems with reference tapes, and 52 ENDF evaluation files.

## Current state: HONEST ASSESSMENT

**The architecture is done. The physics formulas match the Fortran (verified by 6 reviewer reports). But end-to-end validation is incomplete.**

| What works | What doesn't |
|------------|-------------|
| All modules load and export | Only 7/85 NJOY tests fully pass |
| Resonance XS evaluation matches Fortran | 11 tests blocked by ENDF float parsing bugs |
| SAMMY/LRF=7 fixed (Sr-88 thermal correct) | 29 tests need module chaining in pipeline |
| Adaptive grid, Doppler broadening work | ~25 tests not yet attempted (format issues) |
| Unit tests cover individual functions | No end-to-end RECONR→BROADR→HEATR→ACER test |

## CRITICAL: What to do next (priority order)

### Step 1: Fix parse_endf_float (beads: NJOY.jl-i1q)
**File:** src/endf/io.jl — the `parse_endf_float` function

Two bugs discovered by the validation pipeline:
1. `"2.00300e+ 3"` — space between `+` and `3` in exponent. Julia's `parse(Float64, ...)` rejects this. Fix: strip spaces from exponent before parsing.
2. `"2.530000-2"` — compact ENDF format (no `E` before exponent). This IS handled in `parse_endf_float` but some OTHER code path is calling `parse(Float64, ...)` directly instead of going through `parse_endf_float`. Find and fix all direct `parse(Float64, ...)` calls on ENDF data.

**Tests blocked:** 07, 09, 10, 11, 15, 16, 17 (and likely more)

### Step 2: Wire module chaining in test pipeline (beads: NJOY.jl-7zv)
**Files:** test/validation/test_executor.jl, test/validation/njoy_test_runner.jl

The validation pipeline currently only runs `reconr()` and compares against reference PENDF. But most reference tapes include effects from BROADR (Doppler broadening), HEATR, PURR, etc. The pipeline must chain:

```
reconr output → doppler_broaden(result, T) → compute_kerma(result, ...) → build_ace(result)
```

The NJOY input deck tells you what modules to run and in what order. The input parser (`test/validation/input_parser.jl`) already parses this. The executor (`test_executor.jl`) has stubs for chaining but only reconr→broadr is partially wired.

**Tests blocked:** ~29 tests where reference tape includes downstream processing

### Step 3: Run the SAMMY reviewer (beads: NJOY.jl-eux)
Per the project rules, the SAMMY/LRF=7 implementation needs a skeptical reviewer (line-by-line diff against samm.f90). Also sammy.jl is 641 lines (limit is 300) — needs splitting.

### Step 4: Iterate
After steps 1-2, run `julia --project test/validation/run_all.jl` (or run_all_tests.jl). Each failure will expose a real bug. Fix it, re-run, repeat. This is the grind phase.

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

6. **NEVER run parallel Julia test processes.** This was the #1 time-waster in this session. Julia's precompilation cache (~/.julia/compiled/) is NOT safe for concurrent access. When two `Pkg.test()` or `julia test/runtests.jl` processes run simultaneously, they corrupt each other's cache, causing cryptic errors that look like real test failures. **Always:** kill all Julia processes, clear cache (`rm -rf ~/.julia/compiled/v1.12/NJOY*`), then run exactly ONE test.

7. **Use beads (bd) for issue tracking.** Run `bd ready` to see open work. Run `bd close <id> --reason "..."` when done. See AGENTS.md for full bd workflow.

8. **Each source file must be ≤300 lines.** Split by physics concept. Currently sammy.jl (641), wimsr.jl (516), dtfr.jl (506), powr.jl (674) violate this.

9. **After every wave/phase, restate the 3+1 rules.** LLM attention drifts. Repeat the rules explicitly before starting new work.

10. **Read the Fortran before writing Julia.** The reference is at njoy-reference/src/. When the Julia disagrees with NJOY output, the Julia is wrong until proven otherwise.

## How to run things

```bash
# Run unit tests (takes ~60-90 seconds)
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e 'using Pkg; Pkg.test()'

# Run validation pipeline (takes ~2-3 minutes)
julia --project=. test/validation/run_all.jl

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
| prompt.md (now docs/design.md) | Original project specification with full requirements |
| reports/RESEARCH_REPORT.md | 928-line analysis of all NJOY2016 Fortran source |
| njoy-reference/src/reconr.f90 | RECONR Fortran (5747 lines) — the most important module |
| njoy-reference/src/samm.f90 | SAMMY Fortran (7169 lines) — complex R-matrix |
| njoy-reference/src/broadr.f90 | BROADR Fortran (1947 lines) — Doppler kernel |
| njoy-reference/tests/ | 85 test problems with reference output tapes |

## Lessons learned (save yourself time)

1. **Precompilation races are the enemy.** If tests fail with weird errors after agent work, it's almost certainly cache corruption. Clear and retry once.

2. **The ENDF float format is treacherous.** `1.234567+8` (no E), `1.23456-38` (2-digit exponent), spaces in exponents, blank fields — every edge case exists in real data. The parser must be bulletproof.

3. **Sum rule tests need `>=` not `==`.** After adding non-primary MTs to total (gaspr, inelastic channels), total > elastic + fission + capture. Use `total >= sum_of_parts - epsilon`.

4. **LRF=7 (SAMMY) is the hardest formalism.** Variable channel count, per-channel radii, eliminated capture channel, background R-matrix. The Fortran is 7169 lines for a reason.

5. **Reference tapes include ALL processing steps.** A "RECONR test" reference tape may actually include BROADR+HEATR effects if those modules appear in the input deck. Always check the input file to know what processing was applied.
