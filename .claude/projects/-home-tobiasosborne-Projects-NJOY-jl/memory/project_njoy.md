---
name: NJOY.jl Port Status
description: Port of NJOY2016 Fortran nuclear data processing to Julia — scope, progress, and key decisions
type: project
---

Port of NJOY2016 (120k lines Fortran 90) to idiomatic Julia. Target: bit-compatible PENDF/ACE output, differentiable processing chain, <20k lines Julia.

**Why:** Enable differentiable nuclear data processing, uncertainty quantification, and compositional module development. NJOY2016 is the standard tool but monolithic Fortran.

**How to apply:** Follow the 5-wave porting order from prompt.md. Use 3+1 agent workflow (researcher → dual proposals → skeptical reviewer) for every major decision. No code merged without reviewer PASS.

**Status as of 2026-03-21:**
- Wave 0 (Research): COMPLETE — RESEARCH_REPORT.md produced (928 lines)
- Wave 1 (Foundation): IN PROGRESS — Project skeleton created, 186 tests passing
  - Physics constants: done (matches phys.f90 exactly)
  - ENDF types (CONT/LIST/TAB1/TAB2): done
  - ENDF I/O (read/write): done
  - Interpolation (all 6 laws): done
  - Penetrability/shift/phase (l=0..4): done
  - Resonance type hierarchy: done
  - Need: test against real ENDF files, resonance readers (File 2 parsing)
- Wave 2 (RECONR): NOT STARTED
- Reference repos cloned: njoy-reference/, endftk-reference/
- ENDF-6 manual downloaded: endf-data/endf102.pdf
