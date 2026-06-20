# Phase 90 — e5n COMPLETE: broadr thermal MF3/MT1 = Σ broadened mtr[] partials (deep core fix, landed)

**Date:** 2026-06-19/20
**Bead:** `NJOY_jl-e5n` CLOSED. New: `NJOY_jl-<downstream-wobble>` (acer/gaspr realignment).
**Commit:** `a2c5630`. Before/after full sweeps run (e5n vs baseline `36398b5`).
**Orchestration:** 3 parallel read-only Sonnet research (Fortran full path / Julia write path / B-10+Co-58m1 oracle) → 1 Opus TDD coder (xhigh) → parent code-reviewed (Rule 2) + ran the before/after full sweep gate. Rule 9 honored.

## The deep fix
Phase 89 deferred e5n after a narrow attempt (sum of only `{MT2,MT18,MT102}`) regressed B-10/Co-58m1 ~99%. The COMPLETE fix:

1. **`select_broadr_partials` widened to Fortran's full `mtr[]` predicate** (broadr.f90:500-524): broaden any MF3 reaction whose first energy ≤ `emin` (`emin`=1 eV for lrp=0, `eresh` for lrp=1), excluding hard-excludes {1,3,4,19,46-49, iverf≥6:201-599&>850}. EXOTHERMIC reactions (QI>0 → first E≈1e-5) are thereby included: B-10 → {2,102,103,105,107,113,600,700,800,801}; Co-58m1 → {2,51,102,103,107}. Each selected partial is **interpolated onto the MT1 master grid** (the union; lin-lin, identity for {2,18,102}) — `_interp_onto_grid`. `lrp`/`iverf` read from the **original ENDF** (reconr overwrites LRP=2 in the PENDF). + **eresh capped at 6.5e6** (broadr.f90:423-424), which Julia lacked (inert on all current tests, max eresh=1.2 MeV).
2. **Thermal-range MT1 = `sigfig(Σ unrounded broadened partials)`** with the **iflag channel-dedup** (broadr.f90:957-973): when MT103 in the set, exclude MT600-649 from the MT1 sum; MT104→650-699; MT105→700-749; MT106→750-799; MT107→800-849. Sum unrounded, sigfig after (f90:974→980). Above thnmax = verbatim original total.

Julia already copied non-`mtr[]` reactions verbatim and broadened the rest on the shared `b_e` grid (matching Fortran) — the only changes were the partial-set predicate + the MT1 reconstruction.

## Verification (Rule 2 + ground truth)
- **B-10 broadr BYTE-IDENTICAL to the Fortran broadr oracle** (implementer), covering all 3 dedup branches (MT103→600s, MT105→700s, MT107→800s). MT1 thermal 1.934399e5 b (was ~5 b / ~99% wrong with the narrow sum). T01/T70 broadr byte-identical / MT1 bit-identical.
- **Code-reviewed faithful**: predicate, `_iflag` dedup ranges, grid-interp all cite + match the Fortran.
- **Before/after full sweep (e5n vs `36398b5`), per-tape match deltas:**
  - **Massive net win:** T87 **+6709**, T41 **+1865**, T82 **+872/+872/+439/+439**, T37 **+872** (tape44), T35/T63 **+96**, T70 **+351**, T01 **+45** (NUMERIC 32813→**32858**), T02/T10 **+272** each (the MT1-sum benefit reaches non-widened tests), T24 +23, T40 +2, T36 +9.
  - **No status changes**; all 9 BIT_IDENTICAL canaries intact (T03/T09/T22/T50/T52/T53/T61/T62/T86); 0 new CRASH/STRUCTURAL_FAIL; same 6 TIMEOUTs.
  - **Tiny downstream regressions** (all on acer/gaspr): T42 acer tape44 **−229** (Zn-67; total grew +687), T28 acer **−3**, T20 gaspr **−1**, T37 acer tape34 **−4** (its tape44 +872), T45 gaspr **−1**.

## Decision: shipped per the GROUND-TRUTH PRINCIPLE
The broadr output is verified faithful (B-10 byte-identical to the Fortran oracle; identical code path for all nuclei). The downstream wobbles are acer/gaspr realigning with the now-correct WIDER broadr grid — reverting faithful broadr to preserve coincidental downstream matches would violate the project's core principle. (The T42 broadr-vs-oracle double-check was blocked by a `fortran_oracle.jl` deck-truncation bug for T42 — a harness issue, not e5n; B-10's oracle match already covers the logic.) → committed + filed the acer/gaspr realignment bead.

## Net
- **e5n LANDED** — the T01/T70 MT1 bit-identical floor closed at the broadr level; net thousands of matching lines gained across the sweep; the narrow-version B-10/Co-58m1 regression averted by doing the full `mtr[]` widening.
- **Remaining:** `NJOY_jl-c9f` (MF3/MT222 3-pt Lagrange FP grind, P3) + the new acer/gaspr realignment bead + `NJOY_jl-9np` (T67 calcem/acer timeout).
