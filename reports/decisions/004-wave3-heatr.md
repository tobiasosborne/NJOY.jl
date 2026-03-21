# Decision 004: Wave 3 HEATR Architecture

**Date:** 2026-03-21
**Status:** APPROVED

## Dual Proposals
- Proposal A: 421 lines (over limit), comprehensive but monolithic
- Proposal B: 290 lines (under limit), clean per-reaction separation, fixed aniso formula bug

## Selection: Proposal B
**Rationale:** Under 300-line limit, cleaner per-reaction function design, caught and fixed the anisotropic elastic heating formula (wrong outgoing energy coefficient).

## Reviewer Verdict: PASS (with minor notes)
Elastic heating formula confirmed matching Fortran. Lindhard damage, capture, fission all verified.
See reports/REVIEW_WAVE3_HEATR.md.
