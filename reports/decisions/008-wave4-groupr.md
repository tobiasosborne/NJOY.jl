# Decision 008: Wave 4 GROUPR Architecture

**Date:** 2026-03-21
**Status:** PENDING REVIEWER

## Dual Proposals
Both converged on the same 3-file structure:
- group_structures.jl (223 lines): 6 built-in structures as const NTuple
- weight_functions.jl (210 lines): 6 weight functions + tabulated
- groupr.jl (253 lines): group_integrate, group_average, group_average_shielded

## Key decisions
1. Weight functions are plain callables (any f(E)→Float64), not an enum
2. Exact panel_integral for piecewise-linear data (no quadrature)
3. Self-shielding as a separate concern wrapping the weight function
4. All files under 300-line limit
5. 12,590 lines of Fortran → 686 lines of Julia (18:1 reduction)
