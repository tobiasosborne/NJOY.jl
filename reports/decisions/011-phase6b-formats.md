# Decision 011: Phase 6b — Format Output Modules

**Date:** 2026-03-21
**Status:** PENDING REVIEWER

## Dual Proposals

### Pair 1: ccccr + matxsr
- Proposer A: ccccr.jl (421, over limit), matxsr.jl (367, over limit)
- Proposer B: ccccr_b.jl (214), matxsr_b.jl (250) — type-safe with validation

**Selection: Proposer B** — under 300 lines, better type safety and validation. Keep A's versions as they have more complete binary record writing.

### Pair 2: wimsr + dtfr + powr
Both proposers converged. Files over 300-line limit (516, 506, 674) — need splitting.

## Key decisions
- CCCC ISOTXS format with Fortran unformatted record markers
- MATXS with full particle/material/submaterial hierarchy
- WIMS-D/E with burnup chains and resonance tables
- DTF-IV/ANISN with CLAW edit layout
- EPRI-CELL fast/thermal + EPRI-CPM with CLIB subfile structure
