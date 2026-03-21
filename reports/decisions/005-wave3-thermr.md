# Decision 005: Wave 3 THERMR Architecture

**Date:** 2026-03-21
**Status:** PENDING REVIEWER

## Dual Proposals
Both converged on the same structure (292 lines each):
- Free gas model (analytical XS + differential kernel)
- S(alpha,beta) model (tabulated + kernel)
- Bragg edges (coherent elastic, step function)
- Incoherent elastic
- Standard 118-point thermal energy grid from NJOY2016

## Key design decisions
1. Detailed balance enforced by construction (exp(-beta/2) symmetrization)
2. Free gas formula corrected: sign in exp(-x²) term matters critically at low E
3. Both free gas and S(α,β) use the same kernel interface
4. All pure functions, AD-compatible

## Pending: Reviewer verdict
