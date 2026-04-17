# T49 Zr-90 — URR upper-boundary sentinel fix

## Date: 2026-04-17

## Summary

T49 (Zr-90, MAT=4025, err=0.001) went from **41/46 → 44/46 MTs
bit-identical** (+3). Zero regressions across all 16 currently
bit-identical tests (T01, T02, T04, T07, T08, T18, T19, T20, T26,
T27, T34, T45, T46, T47, T55, T84).

## Root cause

Both Fortran `rdfil2` (reconr.f90:765-767) and Julia
`_add_mf2_nodes!` (reconr_grid.jl:92-93) add the URR upper-boundary
shaded pair `{sigfig(EH_URR,7,-1), sigfig(EH_URR,7,+1)}` to the
initial node list. For Zr-90's URR range [200000, 1780460] this is
`{1780459, 1780461}`.

**Fortran's lunion silently drops 1780461 via a sentinel mechanism.**
When the sub-panel loop at label 240 reads the LAST enode into `eg`,
the check at line 2051 `if (ig.ge.ngo) go to 250` routes control
past the `en=abs(eg)` branch, so 1780461 never becomes a panel
boundary and never reaches label 325's write to `inew`. Subsequent
MTs inherit this: after MT=2 finishes, 1780461 is not in `iold`
anymore.

**Julia's lunion_grid keeps 1780461** because Julia has no per-section
iold/inew swap; initial nodes stay in `grid` throughout, and there's
no equivalent of the `ig.ge.ngo` sentinel.

The extra grid point at 1780461 cascaded as an extra line in 5 MTs
(1, 2, 4, 91, 102) for T49.

## Why T27 Pu-239 (currently bit-identical) wasn't affected by the bug

For T27, MT=2 has an **explicit breakpoint** at `x[87]=30000.01`
(the evaluator pre-shaded the URR boundary). So 30000.01 enters the
grid via MT=2's own per-section processing, independent of the
enode sentinel.

Similarly for T02/T18/T19/T34/T47: MT=2 has a **duplicate pair**
`x[k]=x[k+1]=EH`, which Julia's lunion_grid shading logic
(reconr_grid.jl:394-403) converts to `sigfig(EH,7,+1)`.

T49 is the only URR test where MT=2 has *neither* a pre-shaded
breakpoint at `sigfig(EH,7,+1)` *nor* a duplicate pair at `EH`. Just
a single `x[78]=1780460`.

## Why HANDOFF §3 was wrong

HANDOFF.md §3 described the bug as "MT=2's x[78]=1780460 ends up in
Julia's grid" — predicting a Julia grid of `{1780459, 1780460,
1780461}`. Empirically, Julia's grid was
`{1780459, 1780461, 1780465}` (no 1780460). The 1780460 from MT=2
is correctly absorbed by Julia's existing label-220 coincidence
check at `reconr_grid.jl:415` (gated to `k == start_k`, and MT=2's
first breakpoint x[1]=1e-5 triggers it nowhere near 1.78 MeV anyway
— MT=2 handles 1780460 via its own panel processing and label-220
path). The real bug was the spurious 1780461, coming in via the
URR-boundary shading in `_add_mf2_nodes!`.

MT=51's threshold is `thrxx = sigfig(1780463.87, 7, +1) = 1780465`
— the "1780465" in Fortran's grid. This is added by Julia's
threshold replacement at `reconr_grid.jl:291-305` for MT=51's first
breakpoint, and is not the issue.

## Fortran diagnostic evidence

Patched `reconr.f90` `lunion` with `write(*,…)` statements at labels
210 (merge write), 220 (coincidence check), 222 (absorption), 240
(sub-panel entry), 300 (panel-too-small check), 325 (panel lower
bound write), 350 (final point), 410 (boundary dedup). Built and
ran against T49.

Key diagnostic outputs for MT=2 near 1.78 MeV:

```
L240E mt=2 eg=1780459 er=1780460 enl=1719580 ig=362 ngo=363
L240  mt=2 eg=1780459 er=1780460 enl=1719580              (en=|eg|, write path)
L300  mt=2 en=1780459 enl=1719580 -> F
L240E mt=2 eg=1780461 er=1780460 enl=1780459 ig=363 ngo=363  (sentinel!)
L300  mt=2 en=1780460 enl=1780459 -> F
L325W mt=2 ent=1780459.000000178
L240E mt=2 eg=1780461 er=1800000 enl=1780460 ig=363 ngo=363  (sentinel!)
L325W mt=2 ent=1780460.000000000
```

All subsequent `L240E` events for MT=2 have `ig >= ngo=363` and
`eg=1780461` (buffered), but control is routed to label 250 before
any write. 1780461 is never written.

`L410R` (label 410 removals) fires for `eg=1780460.000000000`
(removed by `eresm` boundary dedup at line 2212), but NOT for
1780461 (its distance from `eresh=1780460` is too large for the
`1e-9` relative tolerance).

## Fix

`src/processing/reconr_grid.jl`: new `_drop_unsupported_urr_plus_boundary!`
post-processor. For each LRU=2 range in MF2, computes
`sigfig(EH_URR, 7, +1)` and checks whether this value is in `grid`.
If it is, and NO MF3 section has either:

1. an explicit breakpoint at `sigfig(EH_URR, 7, +1)` (Pu-239 MT=2
   x[87]=30000.01 case), or
2. a duplicate pair `x[k]=x[k+1]=EH_URR` (Pu-238 / Cf-252 / Pu-241
   / Pu-240 case, which lunion_grid shades to
   `sigfig(EH_URR, 7, +1)` via the duplicate-pair handling at
   reconr_grid.jl:394-403),

then remove it. Mirrors Fortran's `ig.ge.ngo` sentinel: the URR+1
boundary node is only in the final grid if an MF3 section
contributes it via per-section processing.

`src/processing/reconr.jl`: invoke the post-processor right after
`lunion_grid` in both the main pipeline (around line 385).

## Verification

Regression matrix (all URR-containing + all bit-identical tests):

| Test | Material         | Pre-fix | Post-fix | Δ   |
|------|------------------|---------|----------|-----|
| T01  | C-nat            | 29/29   | 29/29    | =   |
| T02  | Pu-238           | 17/17   | 17/17    | =   |
| T04  | U-235 (err=0.10) | 27/27   | 27/27    | =   |
| T07  | U-235            | 27/27   | 27/27    | =   |
| T08  | Ni-61 RM         | 18/18   | 18/18    | =   |
| T18  | Cf-252           |  9/9    |  9/9     | =   |
| T19  | Pu-241           | 23/23   | 23/23    | =   |
| T20  | Cl-35 SAMMY      | 162/162 | 162/162  | =   |
| T26  | Pu-245           | 23/23   | 23/23    | =   |
| T27  | Pu-239           | 49/49   | 49/49    | =   |
| T34  | Pu-240           | 53/53   | 53/53    | =   |
| T45  | B-10             | 53/53   | 53/53    | =   |
| T46  | Fe-56 JEFF       | 73/73   | 73/73    | =   |
| T47  | Pu-239           | 49/49   | 49/49    | =   |
| **T49** | **Zr-90**     | **41/46** | **44/46** | **+3** |
| T55  | Fe-56 TENDL      | 61/61   | 61/61    | =   |
| T84  | H-2              |  4/4    |  4/4     | =   |

## Remaining T49 diffs

MT=1 and MT=2, one line each at E=110487.7 eV, ±1 in 7th sigfig:

- MT=1: 6.504111 (Julia) vs 6.504112 (Fortran)
- MT=2: 6.497532 (Julia) vs 6.497533 (Fortran)

Same energy, same ±1 direction. MT=102 is bit-identical, so this is
purely MT=2 elastic (MLBW) FP accumulation at one energy in the
resolved range. Same class as the label previously applied to T34's
3 MT=102 ±1 diffs before they were traced with gdb and fixed in
later phases — a candidate for a similar diagnostic session, but
not pursued here.

## Files changed

- `src/processing/reconr_grid.jl` — `_drop_unsupported_urr_plus_boundary!`
- `src/processing/reconr.jl` — call the post-processor after `lunion_grid`

## Trap (NEW)

**Trap N (URR sentinel)**: Fortran `lunion`'s sub-panel loop at label
240 uses `ig.ge.ngo` (inclusive), so the *last* enode in the sorted
enode list is never processed as a sub-panel boundary. It only
survives into the final grid if some MF3 section contributes it via
per-section processing (explicit breakpoint at the value, or
duplicate pair at the unshaded parent that triggers shading). For
URR upper-boundary shaded pairs `sigfig(EH,7,±1)`, the `+1` variant
is the last enode when there are no resonance peak nodes / URR
energy-table nodes above it. Julia's lunion_grid has no sentinel
analogue; a post-processor is required to match Fortran.
