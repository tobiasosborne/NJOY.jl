# Phase 87 — thermr SAB: restore the `sigfig(emax,7,-1)` held point at the thermal cutoff

**Date:** 2026-06-15
**Beads:** `NJOY_jl-js4` (closed) — sub-work of `NJOY_jl-c3q`; spun off `NJOY_jl-nm2` (free-gas, P3) + a broadr-MT1-floor investigation bead (P2).
**Vehicle:** T70 (Al-27, thermr, nbin=20, iinc=2 S(α,β) + icoh=1 coherent elastic).
**Outcome:** **T70 tape60 MF3/MT221 + MF3/MT222 NP off-by-one (1484 → 1485) FIXED.** MF3/MT222 now byte-identical to `referenceTape60` (zero diffs); MF3/MT221 structural / `>1e-5` diffs gone (14 residual diffs all ≤4.8e-7, the Sigma1/calcem FP-grind class tracked under c3q). T01 NUMERIC_PASS 32812/32962 + T09 BIT_IDENTICAL 1830/1830 preserved.

## Orchestration shape (ultracode)

1. **RED capture + characterization** (1 background Julia, serial per Rule 9): ran T70 once into a fixed work_dir, then a **fixed-column field-level diff** (`/tmp/t70_diff2.py`) vs the reference. **Rule-2 catch:** the first harness used the reference test's `line_equivalence`/`FLOAT_RE` and reported "71767 numeric + 13598 structural diffs, maxrel 0.9" — an **artifact** of the greedy float regex merging the last 11-col data field with the trailing MAT number. The fixed-column parse is ground truth: only ~673 of 213800 lines differ (99.7% match).
2. **Research fan-out** (Workflow, 3 read-only sonnet agents, no Julia → Rule 9 safe): mapped the calcem value pipeline (Fortran σ/cosine, Fortran MF3 emit + ref decomposition, Julia thermr/writer). Produced conflicting root-cause claims for the leading MF3/MT1 diff; **R1's "missing thermr sigfig bias" was refuted by R2+R3 (2 vs 1)** — thermr copies MT=1 verbatim from broadr, so the MT1 grind is a *broadr* artifact (Target B, deferred).
3. **Single opus coding agent** (the only Julia process, background): oracle-driven RED→GREEN TDD + the regression set.
4. **Orchestrator independent verification** (Rule 2): re-read the applied diff at source, re-ran `/tmp/t70_diff2.py` on the produced tape60, and **independently re-ran T01 + T09** (not trusting the subagent's regression report).

## Root cause (confirmed against Fortran, LAW 2)

`_thermr_sab!` (`src/orchestration/modules/thermr.jl`) appended the thermal-cutoff grid as `[emax, sigfig(emax,7,+1), 2e7]`, **omitting `sigfig(emax,7,-1)`** (= 2.276999 for emax=2.277). The reference grid near the cutoff is a 4-point sequence with the thermal value **held** at both `sigfig(emax,7,-1)` and `emax`, stepping to zero at `sigfig(emax,7,+1)`, then `etop=2e7`:

- `thermr.f90:1194` (sigcoh) `enext=sigfig(elim,7,-1)` — the coherent grid's last point **is** `sigfig(emax,7,-1)`.
- `thermr.f90:393` (save-elastic) `enext=emax` — value held at emax (clamp).
- `thermr.f90:394` (save-elastic) `enext=up*emax` — step to zero above emax ≈ `sigfig(emax,7,+1)`.
- `thermr.f90:3211-3213` (tpend) `ex(1)=etop` appended at `ib=ne+1`. MT221 (ix=3) and MT222 (ix=4) are written from the **same `ex` array** → shared energy grid (matches ref: both NP=1485, identical energies).

The missing point → NP=1484 (vs 1485), a trailing-blank "structural" diff on the last data line, and a one-index misalignment over the last two grid points. The section **line count stayed 498** either way (2968 and 2970 values both `ceil` to 495 data lines), so the Phase-86 directory-NC (`3+div(np+1,3)` = 498 for both) remains correct — confirmed by zero MF1/MT451 diffs.

## Fix

`src/orchestration/modules/thermr.jl`: extracted the inline 4-line sentinel block into a documented helper `_append_emax_sentinels!(grid, emax)` (cites the four Fortran lines), and added the load-bearing line:

```julia
push!(grid, round_sigfig(emax, 7, -1))   # held thermal value; coh grid's last point (sigcoh, thermr.f90:1194)
```

The free-gas path (`_thermr_free_gas!`) and `intermediate_e` were left untouched — the free-gas path has the same *shape* but is UNCONFIRMED (no oracle diff yet) → bead `NJOY_jl-nm2` (LAW 1: no fix without a red bar).

## Tests

- **NEW** `test/validation/test_thermr_emax_sentinels.jl` — RED→GREEN unit test: asserts `sigfig(emax,7,-1)` survives in the grid, all four boundary points present, sorted, in order `-1 < 0 < +1 < etop`. RED before fix (`@test em1 in grid` fails), GREEN after (6/6).
- **Phase-86** `test/validation/test_thermr_mf3_nc.jl` re-run — 10/10 (directory-NC formula untouched by the np change).

## Verification (independent — orchestrator, not the subagent)

| Check | Result |
|---|---|
| T70 MF3/MT221 NP | 1484 → **1485 == ref** ✓ |
| T70 MF3/MT222 | **0 diffs** (was 4 incl. STRUCT) — tail byte-identical to ref ✓ |
| T70 MF3/MT221 residual | 14 diffs, all ≤4.8e-7 (FP-grind class, c3q) — structural / `>1e-5` gone ✓ |
| T70 tape60 total | 213800 == ref ✓ |
| T01 (free-gas + SAB chain, NUMERIC) | **32812/32962 — unchanged** ✓ |
| T09 (leapr→thermr, BIT_IDENTICAL) | **tape24 1830/1830 — unchanged** ✓ |

## Remaining on T70 tape60 (not this phase)

- **MF3/MT1** — 471 lines differ at ~1e-6 (max 9.86e-7, all `<1e-5`). MT=1 is copied verbatim from broadr; the diff is the **broadr Doppler/Sigma1 FP floor** (same class as the long-standing T01 NUMERIC-not-BI floor, Trap #14). Concrete lead from R2 (sum-of-broadened-partials vs broaden-the-total) → P2 investigation bead. HIGH regression risk (broadr feeds every test).
- **MF6/MT221** — 177 lines differ at ≤5.2e-6 (calcem σ/cosine FP grind, c3q).
- **MF3/MT221** — 14 lines ≤4.8e-7 (same FP class).
- **free-gas path** emax endpoint — bead `NJOY_jl-nm2`, needs its own oracle.

## Lesson

The reference test's regex-based `line_equivalence` is **unreliable for ENDF field-level triage** — its greedy float pattern merges the last data field with the trailing MAT number, wildly inflating both the diff count and the apparent magnitude. Always cross-check section-level diffs with a **fixed-column (6×11) ENDF parse** before believing a "huge value error" — here it turned a phantom "39% of the matrix off by up to 90%" into the real "one missing grid point + a ~1e-6 FP floor."
