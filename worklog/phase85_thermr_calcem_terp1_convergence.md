# Phase 85 — thermr calcem E′-adaptive convergence: `terp1` chord, not midpoint average

**Date:** 2026-06-14
**Beads:** `NJOY_jl-7dp` (closed), `NJOY_jl-h61` (new, open) — both under `NJOY_jl-c3q`
**Commit:** (this phase)
**Vehicle:** T70 (Al-27, thermr, nbin=20, iinc=2 S(α,β) + icoh=1 coherent elastic)
**Outcome:** **T70 tape60 MF6/MT221 +26-line structural overshoot ELIMINATED** — total 213826 → **213800 (== ref)**, MF6/MT221 185321 → **185295 (== ref)**, all 101 incident-energy blocks match; matching lines 125483 → **128064**. T01 NUMERIC unchanged (no regression).

## Orchestration shape (ultracode)

1. **Beads sync** — local `master` was 17 commits behind `origin/master`; the "stale beads" were the Dolt DB lagging the git-tracked JSONL. Fast-forwarded to `c0c910a`, `bd import` (62 issues + 3 memories).
2. **Research fan-out** (Workflow, 3 read-only sonnet agents, no Julia → Rule 9 safe): Fortran `thermr.f90` calcem rules / reference `tape60` exact counts / Julia current-state loops. Produced two competing root-cause hypotheses for the +26 (H1 last-seed zeroing, H2 `half_tol`).
3. **Single opus coding agent** (the only Julia process, background): oracle-driven TDD — RED + a decisive *per-incident-energy NW block diff* to pick the hypothesis empirically.
4. **Orchestrator verification** (Rule 2): re-read every Fortran claim at source; independently re-ran T70 + T01; **counted the actual section lengths myself** — which caught a mis-report (see below).

## Root cause (the real one — both research hypotheses were wrong)

calcem's E′-adaptive panel subdivision (Fortran `thermr.f90:2041-2068`) tests a freshly-evaluated `sigl` value at the panel midpoint against the **linear-interpolated reference** `ym`:

```fortran
xm = half*(x(i-1)+x(i)); xm = sigfig(xm,8,0)          ! 2051-2052 — SIGFIG-ROUNDED midpoint
call terp1(x(i),y(k,i),x(i-1),y(k,i-1),xm,ym,2)        ! 2059 — law 2 (lin-lin) CHORD at xm
```

`terp1` law 2 = `y1 + (x-x1)*(y2-y1)/(x2-x1)` (`endf.f90:1614`). Because `xm` is the **sigfig-rounded** midpoint, it is generally *not* the exact midpoint, so `terp1(xm) ≠ ½(y_lo+y_hi)`. Julia used the simple average `0.5*(stk[depth]+stk[depth-1])`, which shifted the add/skip convergence decision at many midpoints in the densely-refined high-incident-energy blocks → net **+7 E′-entries** spread (mixed sign) across blocks 78–94 → **+26 ENDF data lines** under `ceil(NW/6)`.

### Why the two research hypotheses were rejected (Rule 2)
- **H2 (`half_tol = 0.5*tol` is 2× too strict):** FALSE. Fortran `sigl` *itself* does `tol = half*tolin` (`thermr.f90:2694`), so Julia's `half_tol` is **faithful**.
- **H1 (last-seed zeroing):** FALSE for this symptom. Last-seed zeroing would give uniform +1 entries on affected blocks; the measured per-block deltas were mixed-sign and up to +4/−2 — the signature of a convergence-test divergence, not a trailing-point keep.

The decisive evidence was the **per-block NW diff** (10 blocks differ, mixed sign, multiples of NL=22) — it ruled out both hypotheses before any code changed.

## Fix

`src/processing/thermr.jl` — both calcem paths (S(α,β) ~L1220-1248, free-gas ~L1383-1397): replaced the two `0.5*(y_lo+y_hi)` midpoint averages (σ and each cosine) with a `terp1_lin` law-2 chord evaluated at `xm`. Cited `thermr.f90:2058` + `endf.f90:1614`. +28/−5 lines, one file.

## Verification (independent — orchestrator, not the subagent)

| Check | Result |
|---|---|
| T70 tape60 total | 213826 → **213800 == ref** ✓ |
| MF6/MT221 section (awk) | 185321 → **185295 == ref** ✓ |
| MF3/MT221, MF3/MT222 sections | **498 == 498** each ✓ |
| matching lines @1e-9 | 125483 → **128064** (59% → 59.9%) |
| T01 (free-gas thermr, NUMERIC) | **32812/32962 — unchanged, no regression** ✓ |

## What the orchestrator's independent count caught

The coding agent reported "MF3/MT221 already matched (498=498)" and the bead claimed a "grid off-by-one (497 vs 498)." Both were partly wrong. My own section count proved the **MF3 sections are correct (498 lines)**; the residual first-diff at reference_test lines 69-70 is the **MF1/MT451 directory NC value** (Julia writes 497 for MF3/MT221+222, ref says 498) — a directory-metadata off-by-one in the PENDF writer, **not** a Bragg-grid bug. The research agents (R1/R3) had hypothesized a `build_thermal_grid` first-past-emax mechanism that does not actually apply (R2 was right: "no off-by-one exists" in the grid). Reframed + tracked as **`NJOY_jl-h61`**.

## Remaining on T70 tape60 (not this phase)
- **`NJOY_jl-h61`** — MF1/MT451 directory NC = 497 vs 498 for MF3/MT221+222 (2 lines; directory writer).
- **value grind** (~40% of lines, the σ/cosine FP at >1e-9) — tracked under `NJOY_jl-c3q`.

## Lesson
"Off-by-one in a section" claims must be checked against the **actual section line count**, not a directory NC field or a single first-diff line — the two diverge, and they need opposite fixes (grid vs directory writer). Independent re-counting (Rule 2) is what separated the real win from the mis-framed residual.
