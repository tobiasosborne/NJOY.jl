# Decision 007: Wave 4 ACER Architecture

**Date:** 2026-03-21
**Status:** PENDING REVIEWER

## Dual Proposals
Both created ace_types.jl + ace_writer.jl. Merged result has both flat (ACETable) and structured (ACENeutronTable) representations.

## Key decisions
1. NXS/JXS indices as named constants matching acefc.f90
2. Type 1 ASCII format with exact field widths (i20/1pe20.11)
3. ESZ block: 5 parallel arrays (energy, total, absorption, elastic, heating)
4. build_ace_from_pendf converts eV→MeV automatically

## Known issues
- ace_types.jl (406 lines) and ace_writer.jl (467 lines) exceed 300-line limit
- Both proposals wrote to same files, creating some redundancy
