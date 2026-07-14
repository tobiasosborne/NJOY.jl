# Phase 92 — staleness reconciliation and T45 photon totals

Date: 2026-07-14

## Scope

This was an orchestrated session with staleness handled before feature work.
Three independent read-only audits covered documentation, embedded Beads state,
and architecture. A separate 3+1 review of `NJOY_jl-2k4` covered Fortran
RECONR/BROADR/GASPR behavior, the T45 ENDF/oracle, current Julia state, and the
sole serial Julia runner.

## Staleness milestone

The embedded Beads database was behind tracked `.beads/issues.jsonl`; `bd
import` restored the tracked state before triage. Obsolete or completed issues
were closed/reopened/reframed, and six concrete documentation/architecture
issues were recorded. README, HANDOFF, architecture, design, tutorial, API,
and archive indices were refreshed in commit `09d0f84`. The old HANDOFF and
design text were archived byte-for-byte.

## `NJOY_jl-2k4`: corrected premise

The original “MF12/13/14 pass-through” diagnosis was wrong. Direct Fortran
review established:

- `reconr.f90:1771-2197` (`lunion`) accepts MF12 only for `LO=1` and accepts
  MF13, forces `NK=1`, keeps the first total TAB1, and adds its nodes to the
  shared union grid.
- `reconr.f90:4726-4960` (`emerge`) interpolates and rounds those totals on the
  final grid, trims the leading zero region, and emits one lin-lin region.
- MF14 is deliberately excluded.
- `broadr.f90:862-999` copies every non-MF3 section unchanged; GASPR likewise
  preserves the reconstructed photon sections.

## Red bar

Baseline `reference_test.jl 45` was `STRUCTURAL_FAIL`, 6634 generated lines
versus 7188, with no generated MF12/MF13. The new focused regression names the
oracle slices and compares columns 1-66:

- MF12/MT102: oracle lines 6639-6979
- MF13/MT4: oracle lines 6982-7104
- MF13/MT103: oracle lines 7106-7184

Before production edits it passed 2/5 assertions and failed exactly the three
empty photon-section comparisons. MF14 absence already passed.

## Implementation

- LRU=0 RECONR now reads MF12/MF13 before its early return, includes them in
  `lunion_grid`, and reconstructs the total sections with sigfig-7 values and
  Fortran threshold trimming.
- Legacy single- and multi-material PENDF writers emit the photon sections and
  exact MF1 directory counts.
- BROADR preserves arbitrary non-MF2/MF3 sections in tape order and includes
  them in regenerated directory entries.
- GASPR required no change.

## Evidence

Focused T45 regression: **7/7 PASS**, including exact photon section order and
SEND/FEND structure.

RECONR-only preserved output matched the oracle in columns 1-66:

| Section | NP | MF1 NC |
|---|---:|---:|
| MF12/MT102 | 1014 | 341 |
| MF13/MT4 | 360 | 123 |
| MF13/MT103 | 228 | 79 |

MF14 is absent, as required. The final full T45 pipeline contains the same
exact photon bodies and directory tuples.

Independent final review caught three general-path bugs before staging:

- negative-Q photon sections must be evaluated from lunion's adjusted scratch
  TAB1, not the original TAB1;
- HEAD ZA/AWR/MF12 LG and TAB1 L1/L2 metadata must survive reconstruction;
- MT460/501 and conditionally redundant MT522 must be filtered before both
  union-grid construction and output.

Each was fixed behind a red bar. The fast threshold/eligibility regression is
wired into `Pkg.test()` and passes 6/6. T84's direct oracle at
`referenceTape100:853` moved from 1/3 to **3/3 PASS**, preserving TAB1 L1/L2
`2/2` and matching the complete MF12/MT102 body in columns 1-66.

Serial regression canaries:

- T01 retained `NUMERIC_PASS` at `1e-5`.
- T02, T08, T27, T34, and T46 retained their known structural statuses and did
  not crash.
- T50, T52, T53, T61, and T62 remained `BIT_IDENTICAL` at `1e-9`.

No full 86-test sweep was run.

## Residual split

Full `reference_test.jl 45` improved to 7185 generated lines but remains
`STRUCTURAL_FAIL` against 7188. The first differences are TPID/MF1/MT451
metadata and the tape is short by three records. Photon processing is already
exact, so this was split into `NJOY_jl-1kf` rather than expanding `2k4`.

Final review also noted that BROADR passthrough directory entries default MOD
to zero. T45's relevant MOD values are zero and the section bodies are exact;
preserving nonzero incoming MOD is tracked separately as `NJOY_jl-6lg`.
