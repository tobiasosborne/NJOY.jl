# Decision 006: Wave 4 UNRESR/PURR Architecture

**Date:** 2026-03-21
**Status:** PENDING REVIEWER

## Dual Proposals
Both converged on the same architecture:
- `unresr.jl` (255 lines): URRStatModel, Hwang quadrature, Bondarenko self-shielding
- `purr.jl` (255 lines): Chi-squared sampling via quantile table, Wigner spacing, ladder generation, probability tables

## Key decisions
1. Deterministic RNG (Xoshiro) for reproducible probability tables
2. Clean separation: statistical model → sampling → physics integration
3. Chi-squared quantile table matching NJOY's 20-point inverse CDF
4. Both files under 300-line limit
