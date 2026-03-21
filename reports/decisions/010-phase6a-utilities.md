# Decision 010: Phase 6a — Utility Modules (gaspr, mixr, resxsr, covr, gaminr)

**Date:** 2026-03-21
**Status:** PENDING REVIEWER

## Dual Proposals

### Pair 1: gaspr + mixr + resxsr
Both proposers converged on same designs:
- gaspr.jl (252 lines): MT203-207 gas production from partial reaction sums
- mixr.jl (199 lines): abundance-weighted material mixing on union grid
- resxsr.jl (267 lines): RESXS format output for thermal flux calculations

### Pair 2: covr + gaminr
Both proposers converged:
- covr.jl (283 lines): covariance→correlation conversion, boxer format, reuses CovarianceMatrix from errorr.jl
- gaminr.jl (282 lines): 9 built-in photon group structures, reuses group_integrate from groupr.jl

## Key decisions
- All 5 files under 300 lines
- Maximum reuse of existing infrastructure (errorr, groupr)
- Pure functions throughout
- 1,283 lines total for 5,958 lines Fortran (4.6:1 reduction)
