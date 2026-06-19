# Phase 89 — thermr build_thermal_grid (km2) + free-gas endpoint (nm2) + broadr MT1 (e5n confirmed/deferred)

**Date:** 2026-06-19
**Beads:** `NJOY_jl-km2` closed, `NJOY_jl-nm2` closed, `NJOY_jl-9np` new (T67 full-pipeline timeout), `NJOY_jl-e5n` confirmed + re-scoped (NOT landed).
**Commits:** `0807b41` (km2), `402117c` (nm2). e5n reverted (patch saved `/tmp/e5n_partial_sum.patch`).
**Orchestration:** parallel read-only Sonnet research (Explore) → serial Opus TDD coders → parent independently re-verified every result on a clean cache. Rule 9 honored (one NJOY Julia process at a time).

## km2 — `build_thermal_grid` under-inclusion + O(n²) perf (CLOSED, commit `0807b41`)

Surfaced by Phase 88b. `build_thermal_grid` (src/processing/thermr.jl, mirrors Fortran `coh` thermr.f90:748-922) dropped 18 input-grid points on T67 (D-in-7LiD, 3374 Bragg edges) — merged grid 7042 vs ref 7060 → MF3/MT221+222+223 = 2351 vs 2357 — and burned ~5.1e9 allocations (~13 min, 600s-killed).

**Two root causes (faithful Fortran port, not bandaids):**
1. **Termination was edge-count-based, not energy-based.** Fortran `coh` keeps fetching "edges" past the last physical Bragg edge — `sigcoh` returns an `enext ≥ emax` sentinel (sigcoh 1192-1195) — and processes the `[last_edge, emax]` interval via label-125 inclusion, terminating only when a popped point `> emax*(1+small)` (lines 870-875). Julia bailed at `bidx > length(bragg_e)`, never processing `[last_edge, emax]` → dropped the high-E input sequence (5.0…9.6875). Fix: append two emax sentinel edges + a `terminated` flag porting the emax cutoff.
2. **DEEPER (the index-2312 / `1.5` drop): `eps`=3e-5 used where Fortran uses `small`=1e-10.** Fortran uses `small` in the phase-1 (f90:797), label-125 inclusion (826, 830) and label-170 window-advance (853) tests; `eps` only in the convergence spacing test (839). Julia used `eps` blanket → input point 1.5 (9e-6 above edge 1.499986) was wrongly judged "not above stack bottom" and dropped. Fix: split the tolerances exactly as Fortran does.

**Perf:** `bragg_edges` (thermr.jl:409-414) is an exact prefix sum `scon/E·Σ_{τ²<E·econ} f_i`. Replaced the four internal O(n_edges) calls with a `cumsum(form_factor)` + `searchsortedfirst` local fast-path — **bit-identical accumulation order** (verified 200k samples). De-boxed `advance_window!`, added `sizehint!`. **~13 min → 2.39 s**, 5.1e9 allocs → 408 KB/call.

**Oracle (new):** `test/validation/test_thermr_build_grid_t67.jl` — isolated, calls `build_thermal_grid` on the real T67 inputs (tsl-DinLiD MF7/MT2 3374 edges + referenceTape60 MT2 grid), asserts merged+sentinels grid == referenceTape60 MF3/MT221 (7060 pts). **13/13, GREEN, 2.39 s** (re-verified clean cache). T70 (the original 568-edge validation case) no regression: 213800==ref, match count not decreased; T09 BIT_IDENTICAL preserved.

**Residual (separate bead `NJOY_jl-9np`):** the FULL T67 reference test still exceeds 580s, but the bottleneck has MOVED off `coh` onto the calcem inelastic S(α,β) matrix + two acer runs. Not a `build_thermal_grid` problem anymore.

## nm2 — thermr free-gas MF3 step-to-zero endpoint (CLOSED, commit `402117c`)

**Bead premise REFUTED, real bug found.** The hypothesis (a missing `sigfig(emax,7,-1)` point) is wrong — Fortran free-gas never calls `sigcoh` (the only source of that point), so its absence is correct. **Real defect:** the free-gas save-elastic loop (thermr.f90:383-394) steps to zero at `up*emax = 1.00001*emax` (param `up`, line 161), but Julia (`_thermr_free_gas!` modules/thermr.jl:122) used `round_sigfig(emax,7,1)`. Ref T01 tape25 MF3/MT221 endpoint = `1.200012` (=1.00001·1.2), Julia emitted `1.200001`. One-line fix → **T01 tape25 32812→32813**, no regression. Path-specific: the SAB *coherent* path correctly steps at `sigfig(emax,7,+1)` via sigcoh — left untouched.

## e5n — broadr thermal-range MF3/MT1 (CONFIRMED against Fortran, DEFERRED — would regress B-10)

**Hypothesis CONFIRMED (LAW 2):** Fortran broadr label-275 (broadr.f90:938-980) reconstructs thermal-range (E≤thnmax) MF3/MT1 as `sigfig(Σ broadened partials over mtr[], − iflag channel-dedup, + URR)`; above thnmax = verbatim original total. Julia broadens the TOTAL directly (`sigma1_at`) — wrong because the Doppler kernel is nonlinear.

**Partial fix worked for light nuclei but was REVERTED:** summing Julia's broadened `{MT2,MT18,MT102}` (broadr.jl:114-145; sum unrounded partials then sigfig, matching f90:974→980) gave **T70 MF3/MT1 BIT-IDENTICAL** (was a 471-line diff) and **T01 32813→32858**. BUT a pre-commit at-risk analysis (read-only) found it **regresses B-10 (T45) and Co-58m1 (T82) by ~99%**: Julia's `select_broadr_partials` = `{MT2,MT18,MT102}` is **narrower than Fortran's `mtr[]`**. EXOTHERMIC non-threshold reactions dominate the thermal total but aren't broadened by Julia — B-10 **MT107 (n,α)** ≈ 3840 b @ 2200 m/s (~99% of thermal MT1), Co-58m1 **MT103 (n,p)** + **MT51** (exothermic de-excitation). Per the user's "don't commit if the sweep regresses anything", reverted (patch `/tmp/e5n_partial_sum.patch`).

**Complete fix (re-scoped in bead, dedicated session):** widen `select_broadr_partials` (auto_params.jl:69-92) to Fortran's full `mtr[]` membership (broadr.f90:497-560) + apply the `iflag` channel-dedup (945-977) + URR getunx add-back. This also broadens+writes more partial MF3 sections (matching Fortran) → broad regression surface → needs full sweep. At-risk set already mapped: among the 51 broadr tests, only T45+T82 have exothermic non-`{2,18,102}` reactions below thnmax (no MF2/MT152 URR in any bundled resource). Lesson recorded via `bd remember` (`broadr-mt1-thermal-reconstruction-e5n-...`).

## Net
- **+1 structural milestone (km2)**, **+1 endpoint fix (nm2, T01 +1 line)**, **e5n confirmed + precisely scoped** (regression averted before commit). No regressions shipped.
- **Remaining cluster:** `NJOY_jl-c9f` (MF3/MT222 3-pt Lagrange FP grind, P3) and the complete `NJOY_jl-e5n` (broadr partial-set widening, P2).
