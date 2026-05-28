# Phase 77 — broadr O(n²)→O(n), premise corrections, fresh baseline, T20 fix, ACER T50 bit-identical

Date: 2026-05-28. Orchestrated multi-milestone session (one Julia process at a
time, per Rule 9). All work committed + pushed as it landed.

## Milestones landed

### 1. broadr U-238 performance — 1039 s → 40 s (~26×), byte-identical
- **Root cause:** `sigma1_at` (sigma1.jl:146) rebuilt the full length-`nseg`
  velocity array `v = map(e->sqrt(alpha*e), seg_e)` on *every* call.
  `broadn_grid` calls `sigma1_at` O(nseg) times → the Doppler kernel was
  O(nseg²) in time and allocations (20k-pt grid ≈ 1.8 GB of velocity arrays).
- **Fix:** `broadn_grid` already computes `vel = sqrt.(alpha .* seg_e)` once
  (broadr.jl:79). Thread that precomputed `vel` into `sigma1_at`; update the 3
  call sites (broadr.jl:106/417, orchestration/modules/broadr.jl:96 — precompute
  once before its total-xs loop) + the unit test. Mirrors Fortran which converts
  E→v once per page (broadr.f90:1318-1322).
- **FP-safety:** `sqrt.(alpha.*seg_e)` is bit-identical to the old per-call map
  (verified over 20k pts). T08 broadr tape23 md5 unchanged; T01 tape25
  32812/32962. Kernel integral logic untouched (`v = vel` alias).
- **Follow-up:** `h_all` used `ntuple(5) do k` → heap-allocated closure per
  interval (~239 B, O(nseg²) total). `ntuple(Val(5))` unrolls at compile time
  (bit-identical; matches Val(50)/Val(62) convention). Landed separately.
- **Impact:** real U-238 broadr now 40 s (verified). **Unblocked T15/T16/T65
  CRASH→DIFFS** — the Phase 71-76 covcal output is finally measurable
  (T15 tape26 5089/5958, tape35 3724/3812). Beads NJOY_jl-58o(.1/.2), kb8.

### 2. leapr `naint` — misdiagnosis corrected (NON-bug)
- HANDOFF claimed a missing `naint` SCT-gating in `generate_sab`. Verified vs
  Fortran: `naint` is hardcoded to 1 (leapr.f90:292-294, never read from input),
  so `iprt=mod(j-1,naint)+1=1` always — the gate is permanently open and Julia's
  unconditional SCT replacement already matches. Fixed the misleading comment +
  HANDOFF. Real T80 residual is phonon-loop FP order (bead NJOY_jl-sc5).
- Bead NJOY_jl-a6q closed not-a-bug.

### 3. Fresh full sweep (P0) — post-Phase-76 + post-broadr baseline
- 80 min. **2 BIT_IDENTICAL, 2 NUMERIC_PASS, 78 DIFFS, 1 CRASH, 1 NO_REFERENCE.**
  Crashes 4→1 vs the stale Phase-67 sweep (T11/T15/T16/T65 now complete).
  `reports/REFERENCE_SWEEP.md` regenerated. Bead NJOY_jl-0fx closed.

### 4. T20 LRF=7 crash — regression fixed (crashes 1→0)
- The fresh sweep exposed T20 (Cl-35 RML LRF=7) crashing in `read_mf32` — a
  regression from the Phase 72-76 rescon work (was DIFFS pre-rescon).
- **Fix:** `read_mf32` now `@warn`(maxlog=1)+skips unported MF32 sub-formats
  (LRF=4/7, LCOMP≠1, NRO≠0, NLRS≠0, ISR=1+LRF=3, LRU≠1) instead of `error()`,
  with a loop guard (variable-length bodies can't be skipped past). Genuine
  structural/record-length violations still error (Rule 6). `apply_rescon!` is a
  clean no-op for zero-usable-range data. INTERIM — proper fix is porting
  `rpxsamm` (errorr.f90:3252) = bead NJOY_jl-4pz.
- T20 CRASH→DIFFS; supported LRF=3+URR byte-unchanged (mf32_reader 83/83,
  covcal_lb5 50/50, tape26=5958, MT=102 C[1,1] canary OK). **Crashes 0/84.**
  Bead NJOY_jl-958 closed.

### 5. ACER T50 — BIT_IDENTICAL (143/163 → 163/163)
- **Root cause (a real ~1e-6 path diff, NOT an FP floor):** `acecpe` reads the
  *sigfig-rounded* ACE ESZ columns, not the raw PENDF. Julia fed raw values.
- **Fix (faithful):** `ace_charged.jl` — `xelas` interp + `signi` use the
  sigfig-7 elastic column (acefc.f90:5497/5499/6659); `_log_log_interp` rewritten
  to mirror `terp1` int=5 exactly (endf.f90:1627) with the acecpe bracket search
  (6653-6657) and ratio-then-log order. `acer.jl` (charged path only) builds the
  ESZ total as sigfig9(Σ sigfig7(reaction_xs)) over MT=2/MT≥5 (5497-5505,6304)
  instead of raw MT=1, so the post-loop total cancels cleanly.
- T50 tape34 163/163 @1e-9 (independently re-verified). No regression: T01
  32812/32962, T02 12519/13873 (change is charged-path-only). Bead NJOY_jl-cnh.1
  closed; epic NJOY_jl-cnh stays open for the siblings.

## Meta-lesson (recorded)
Three consecutive beads (leapr/ACER/Bragg), all derived from a drifted HANDOFF,
had **stale/wrong premises** that "Fortran before Julia" verification corrected
*before* any implementation: leapr naint was a phantom; ACER "unlocks 8 siblings"
was overstated (siblings need separate MF6 work); Bragg "replace the hardcoded
table" was backwards (keep lat=1/2/3 built-ins, ADD the lat=10 MF7/MT2 path,
sigcoh thermr.f90:1016). Always verify a bead's premise against the .f90 first.

## Net state after this session
- **3 BIT_IDENTICAL (T03, T22, T50), 0 CRASH** (the committed REFERENCE_SWEEP.md
  predates the T20 fix + ACER win — next sweep will show 3/0).
- broadr no longer the T15/T17 bottleneck (40 s).

## Remaining backlog (beads) — all multi-session grinds
- T15 covcal sub-section content fidelity (NJOY_jl-zso, now measurable/unblocked).
- ACER siblings (cnh.2/.3): T51-54/T62/T71 each need MF6 multi-MT / LTP<12 work.
- Bragg lat=10 MF7/MT2 read path (NJOY_jl-fvv) — unblocks T25/T67-70/T74.
- groupr MT=251/252 + empty-MT skip (NJOY_jl-01v).
- leapr Phase B phonon FP-order (NJOY_jl-sc5); T49 MLBW ±1 (NJOY_jl-4k2);
  purr fidelity (NJOY_jl-gvs); rpxsamm port (NJOY_jl-4pz); colaps (NJOY_jl-h2d);
  per-tape DIFFS grind (NJOY_jl-tbe); T04/T65/plotr (P4).
