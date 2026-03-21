# Decision 009: Wave 5 MODER + ERRORR

**Date:** 2026-03-21
**Status:** PENDING REVIEWER

## Dual Proposals

### MODER
Both converged: tape copy/extract utility (no binary ENDF in Julia).
- moder.jl (271 lines): moder_copy, read_endf_tape, write_endf_tape, extract_material, merge_tapes

### ERRORR
Proposal B added sandwich formula for differentiability story.
- errorr.jl (257 lines): read_mf33, LB=0/1/5/6 expansion, multigroup collapse, sandwich_covariance
- Key: covariance operations are pure linear algebra — compose with ForwardDiff

## File splitting (concurrent)
5 biggest files split into 13 smaller files, all ≤300 lines.
reconr.jl(684→4 files), faddeeva.jl(524→2), breit_wigner.jl(473→2), ace_writer.jl(467→2), ace_types.jl(406→2)
