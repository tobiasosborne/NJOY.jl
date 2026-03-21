# Decision 003: Wave 3 BROADR Architecture

**Date:** 2026-03-21
**Status:** APPROVED

## Dual Proposals

Both proposers converged on the same architecture:
- `sigma1.jl` — pure mathematical kernel (F/H-functions, SIGMA1 integral)
- `broadr.jl` — pipeline reusing `adaptive_reconstruct` from Wave 2

### Key shared design decisions:
1. Kernel/pipeline separation (sigma1.jl is pure math, broadr.jl is orchestration)
2. Reuse of `adaptive_reconstruct` instead of reimplementing broadn's bisection
3. Three-pass sigma1: below segments, above segments, negative velocity
4. Cancellation avoidance via Taylor series fallback in H-functions
5. AD-compatible: no mutation in sigma1.jl

## Selection

Since proposals converged, no selection needed. The final code represents both.

## Reviewer Verdict: CONDITIONAL PASS → PASS (after fixes)

Must-fix issues resolved:
- h_taylor convergence: matched Fortran's aerr=1e30 behavior
- T_new==T_old: returns input unchanged instead of error
- h_taylor allocations: MVector replaces heap arrays
- Output clamping: matches Fortran sigmin=1e-15

## File sizes
- sigma1.jl: 229 lines (under 300 limit)
- broadr.jl: 156 lines (under 300 limit)
