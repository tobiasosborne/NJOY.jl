# Decision 002: Wave 2 RECONR Architecture

**Date:** 2026-03-21
**Status:** PENDING REVIEWER

## Dual Proposals

### Proposal A (Faithful Translation)
- Direct translation of resxs, lunion, emerge from reconr.f90
- Monolithic sigma function with mode switch
- CrossSections struct-based return values
- Legacy reconr() returning NamedTuple

### Proposal B (Composable Julia Design)
- Generic `adaptive_reconstruct(f, grid, config)` higher-order function
- Closure-based evaluator via `build_evaluator(mf2)` — no global state
- `PointwiseMaterial` with matrix-based cross section storage
- Clean pipeline: `build_evaluator → build_grid → adaptive_reconstruct → merge_background!`
- `AdaptiveConfig` struct with sensible defaults
- Pre-allocated `AdaptiveWorkspace` to minimize inner-loop allocations

## Orchestrator Selection: Proposal B (with legacy compat from A)

**Rationale:**
1. The generic adaptive grid can be reused for BROADR and THERMR (both use similar bisection)
2. The closure-based evaluator eliminates global state — key project requirement
3. NTuple{N,Float64} return type enables type-stable compilation
4. Pre-allocated workspace reduces GC pressure in hot loops
5. The composable design enables `cross_section(E, material)` standalone use

**Compromise:** Both legacy (`reconr()` returning NamedTuple) and new (`reconstruct()` returning PointwiseMaterial) interfaces kept. File exceeds 300-line limit and should be split.

## Test Results
6244 tests pass including:
- Adaptive grid on Lorentzian and step functions
- Full H-2 reconstruction with reference tape comparison
- Physical invariants (sum rule, positivity, monotonic grid)

## Reviewer Verdict: CONDITIONAL PASS

Critical issues C1 (silent exception swallowing) and C2 (step-ratio guard eresr condition) fixed.
Remaining items tracked as beads issues: NJOY.jl-7jl, NJOY.jl-2dm, NJOY.jl-4n0, NJOY.jl-2c8.
See reports/REVIEW_WAVE2.md for full review.
