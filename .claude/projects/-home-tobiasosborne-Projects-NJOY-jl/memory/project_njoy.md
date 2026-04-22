---
name: NJOY.jl Port Status
description: Port of NJOY2016 Fortran nuclear data processing to Julia — scope, progress, and key decisions
type: project
---

Port of NJOY2016 (~100k lines of Fortran 90 code) to idiomatic Julia. Target: bit-compatible PENDF/ACE output, differentiable processing chain, ~24k lines Julia as of Phase 51 (code-only, per cloc).

**Why:** Enable differentiable nuclear data processing, uncertainty quantification, and compositional module development. NJOY2016 is the standard tool but monolithic Fortran.

**How to apply:** Follow the 5-wave porting order from prompt.md. Use 3+1 agent workflow (researcher → dual proposals → skeptical reviewer) for every major decision. No code merged without reviewer PASS.

**Status as of 2026-03-21:**
- Wave 0 (Research): COMPLETE — RESEARCH_REPORT.md (928 lines)
- Wave 1 (Foundation): COMPLETE — CONDITIONAL PASS from reviewer
  - Physics constants, ENDF I/O, types, interpolation, penetrability, Faddeeva, MF2 reader, XS evaluation (SLBW/MLBW/RM)
  - 4783 tests passing
  - Decision: reports/decisions/001-wave1-foundation.md
- Wave 2 (RECONR): CONDITIONAL PASS from reviewer — fixing C1/C2
  - Adaptive grid (generic bisection), RECONR pipeline, PENDF writer
  - 6244 tests passing
  - Remaining: LRU=1 end-to-end tests, lunion linearization, file splitting
  - Decision: reports/decisions/002-wave2-reconr.md
- Wave 3 (BROADR/HEATR/THERMR): NOT STARTED
- Reference repos: njoy-reference/, endftk-reference/
- Issue tracker: beads (bd) initialized, issues tracked
- Total Julia source: ~24 k code lines across 87 files (as of Phase 51, 2026-04-22)
