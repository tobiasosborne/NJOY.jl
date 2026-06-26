# Phase 91 — Orchestrated bead grind: czw, 4k2 (_THIRD), xu3→k1v (lunion mixed-law shading)

**Date:** 2026-06-26
**Mode:** Orchestrated (read-only Sonnet research in parallel; serial Opus implementers; orchestrator commits/pushes per milestone).
**Beads:** `NJOY_jl-czw` (closed), `NJOY_jl-4k2` (closed), `NJOY_jl-xu3` (closed/refuted), `NJOY_jl-k1v` (created+closed), `NJOY_jl-2k4` (created/open), `NJOY_jl-c9f` (researched/deferred).
**Commits:** `f3c68c4` (czw), `895795e` (4k2), `abc1e93` (k1v).

---

## Headline

Three milestones landed; two are **systematic reconr faithfulness fixes**:

1. **czw** — removed dead `thermr_coh_ne` plumbing (8 sites/4 files); provably output-invariant.
2. **4k2** — channel-radius constant `_THIRD = 1.0/3.0` → Fortran-truncated `0.333333333` (CLAUDE.md Trap #3). Fixes T49 MLBW elastic ±1 ULP at E=110487.7 eV; makes Julia `rho`/elastic **bit-identical** to the Fortran oracle. `channel_radius` feeds SLBW/MLBW/RM/SAMMY (all NAPS=0 reconstruction).
3. **k1v** — `lunion_grid` over-shaded **every** breakpoint into `{sigfig(-1),sigfig(+1)}` pairs when a section had **any** histogram region; gated on per-breakpoint `_inta_for_point(tab,k)==Histogram` (Fortran `inta==1`, reconr.f90:2051). Fixes T42 (Zn-67) reconr over-refinement: MT1 grid **52307→51211 = exact Fortran** (the +1096-point low-E cluster gone).

**Methodology note (the through-line of the phase):** three confident hypotheses were **refuted by instrument-first oracle work** before any code changed —
- 4k2: the l=2-penetrability FP-order hypothesis was a verified no-op (perturbation probe → l=2 has zero effect at the divergent E); the real cause was `_THIRD` (candidate the research had flagged "likely latent").
- xu3: the "broadr_module is missing `thinb`" hypothesis was refuted by Fortran source (`bfile3` calls `broadn` OR `thinb` either/or; broadening never thins).
- k1v: the "step guard / peak-triplet hw" hypothesis was refuted by an instrumented trace (step guard + arithmetic midpoint are identical both sides; the seed was lunion shading). The "geometric midpoint" claim in one research analysis was refuted by direct read of reconr.f90:2358.

Every fix was confirmed against the Fortran oracle and regression-gated before commit. "Nothing is irreducible" held again — both 4k2 and k1v were real bugs, not cross-compiler floors (FMA explicitly ruled out for 4k2 via `-ffp-contract=off`).

---

## NJOY_jl-czw — remove vestigial `thermr_coh_ne` plumbing (commit `f3c68c4`)

Dead since h61 (Phase 86) changed the thermal MF3 directory-NC to use the actual emitted np (`3+div(np+1,3)`). The value was still set in `_recompute_thermr_mf6!` and threaded `ctx → final_assembly! → write_full_pendf`, where it was silently ignored (never read in any body). Removed all 8 sites: `RunContext` field + positional constructor default (pipeline.jl), `final_assembly!` call + signature/docstring + `write_full_pendf` pass-through (moder.jl), the writer kwarg (pendf_writer.jl), a stale `copy_with_modifications` docstring (pendf_io.jl).

The one footgun — the **positional `RunContext` constructor** — was handled carefully (removing only the `0` aligned to `thermr_coh_ne`, trailing fields verified unshifted). Provably output-invariant: grep clean; T70 unchanged (DIFFS, tape60 213800==ref); T01 unchanged (NUMERIC_PASS 32858/32962).

---

## NJOY_jl-4k2 — channel-radius `_THIRD` = Fortran-truncated `0.333333333` (commit `895795e`)

T49 (Zr-90) MLBW elastic differed by ±1 in the 7th sigfig at E=110487.7 eV (MT1 6.504111 vs 6.504112; MT2 6.497532 vs 6.497533); MT102 bit-identical.

**Disproof first (Law 1):** the research's primary hypothesis — l=2 penetrability/shift denominator FP-order (`den=9.0+3.0*r2+r4` vs Fortran `3*r2+r4+9`) — was applied and found to be a **verified no-op** (reconr PENDF byte-identical). A perturbation probe (×1.001 kick to the l=2 den, 13 orders > 1 ULP) moved **nothing** at E=110487.7 → no d-wave contribution there; the contributing resonances are l=0/l=1.

**Instrumented trace (the gdb-trace the bead always wanted):** dumped `csmlbw` intermediates Fortran (-O0) vs Julia. First divergence was `rho = k·ra`: `k` bit-identical, so the error is in `ra = _RC1·aw^_THIRD + _RC2`. Julia used `_THIRD = 1.0/3.0` (full precision); NJOY2016 uses `third = .333333333e0_kr` (TRUNCATED, reconr.f90:2877). **CLAUDE.md Trap #3.** The URR path already used the truncated `_THIRD_URR`; the resolved path (slbw.jl) was simply missed.

After the fix Julia's `rho`/elastic are bit-identical to Fortran: T49 MT2 6.497532→6.497533, MT1→6.504112; MF3 cols 1-66 now differ in only the unrelated HEAD-record `0 0` vs `0 99` metadata. **FMA ruled out** (`-ffp-contract=off` still yields 6.497533; instrumentation at -O3 perturbs the knife-edge value, so the true Fortran trace was captured at -O0). Zero regressions: T01/T02/T08/T27/T34/T45/T46.

---

## NJOY_jl-xu3 — REFUTED (closed); surfaced k1v + 2k4

Premise was "e5n's wider broadr grid introduced fixable acer/gaspr **consumption** regressions (T42 −229, T28 −3, T20/T45 −1, T37 −4)." Oracle probe of T42 + T45 refuted it:

- **T42**: the broadr-grid diff is **inherited from reconr** (Julia 52307 vs Fortran 51211 MT1 pts; 11-pt cluster at 1.308719e-5 eV), independent of 4k2 (`_THIRD` A/B test) and e5n (reconr is upstream of broadr). → filed `NJOY_jl-k1v` (fixed this phase).
- **T45**: broadr MF3 is **byte-perfect**; the tape40 failure is Julia dropping MF12/MF13/MF14 photon production wholesale (−554 lines). → filed `NJOY_jl-2k4` (open).

Both original candidates (acer `ie_start`, gaspr `ngas`) refuted. e5n is correct; the apparent "regressions" are metric wobbles on tapes already STRUCTURAL_FAIL for these pre-existing reasons (e5n actually improved the MF3 grids).

---

## NJOY_jl-k1v — lunion mixed-law histogram shading (commit `abc1e93`)

The 1.308720e-5 cluster is seeded by the **initial grid**, not the bisection. `lunion_grid` shaded **every** interior breakpoint into `{sigfig(E,7,-1), sigfig(E,7,+1)}` whenever the section contained **any** histogram region (`is_histogram = any(law==Histogram)`). Zn-67 MT107 (n,α) has mixed law `NBT/INT = 70/5, 99/1, 190/5` (log-log/histogram/log-log) with an explicit **log-log** breakpoint at 1.308720e-5; shading it into `{1.308719e-5, 1.308721e-5}` (Δ≈2e-11) seeded a tiny panel that `adaptive_reconstruct`'s (correct) step-ratio guard then over-refined into ~11 clustered points — +1096 on the full T42 grid, cascading into broadr/acer/gaspr line counts.

**Fortran ground truth:** lunion shades only when the breakpoint's **own region** interpolation `inta==1` (histogram): `if (inta.ne.1) go to 255` (reconr.f90:2051; `inta` set at 2019-2023/2076-2080). **Fix:** gate the shading branch on a new `_inta_for_point(tab,k)==Histogram` helper (the analogue of `inta`, distinct from `_law_for_interval` which uses `idx<NBT[r]`). No-op for pure-histogram and non-histogram sections; only removes spurious shading of non-histogram breakpoints in mixed-law sections.

**Regression (full BI cohort, mandatory for a core-reconr change):** T42 MT1 52307→51211 (cluster gone); all 9 BIT_IDENTICAL canaries (T03/09/22/50/52/53/61/62/86) stay BIT_IDENTICAL; 6 resonance canaries (T01/02/08/27/34/49) unchanged. Surgical: 1 file, +34/−9.

The instrumented trace also definitively resolved the research conflict: midpoint is **arithmetic** both sides (reconr.f90:2358 vs adaptive_grid.jl:251); the step guard (`est=estp·(x(i)-res(1))`, estp=4.1) is **identical** both sides — neither was the bug.

---

## NJOY_jl-2k4 — MF12/MF13/MF14 photon-production pass-through (OPEN, filed)

T45 (B-10) source ENDF carries MF12(41)/MF13(1059)/MF14(6); Fortran reference tape40 keeps MF12(342)+MF13(204); Julia drops them wholesale (−554 lines). Site: reconr/moder section pass-through. Julia has the readers (`read_mf12_lo1_sections`, `read_mf13_sections`) but the output path doesn't emit them for this chain. P3; touches BI-sensitive paths → regression-gate the acer BI cohort.

---

## NJOY_jl-c9f — researched, DEFERRED

Two read-only agents confirmed Julia `_lagrange_eval` (modules/thermr.jl:357) is **bit-identical** to Fortran `terp` (thermr.f90:1522-1536) — the bead's "Lagrange FP order" premise is wrong. Real suspect: `exp(-2·c2)` in `incoh_elastic_xs` (libm), or gfortran -O3 FMA. Per the 4k2 lesson it may still be a real upstream/constant bug, so it warrants the same instrumented-trace treatment (probe T69, which — unlike T74 — does not time out). Deferred to a focused session; bead retains the research.

---

## Full sweep (post-Phase-91 baseline)

Full 86-test sweep at `abc1e93` (93.2 min). Counts: **9 BIT_IDENTICAL · 4 NUMERIC_PASS · 72 DIFFS · 0 STRUCTURAL_FAIL · 1 NO_REFERENCE · 0 CRASH · 0 TIMEOUT** (vs Phase-90 baseline 9 / 4 / 66 / 0 / 1 / 0 / **6**).

Runtime-normalized diff vs the committed Phase-90 baseline — **only 7 tests changed tape results; ZERO regressions:**
- **T42 (the k1v win)**: tape44 Julia line count **106093 → 61721** (−44,372 spurious over-refinement lines removed), matched **8219 → 32296 (×4)**; tape34 37→38. The +1096-point reconr over-refinement cascading downstream into acer is now gone — direct, attributable confirmation of the lunion-shading fix.
- **6 ex-TIMEOUT → DIFFS** (now gradeable): T17 (612.9s), T25, T65, T67 (tape60 106610/463386), T68, T74. These were known *contention*-driven (not hangs; Phase-84 note); they completed this run (the resolution is sweep-environment/time, not attributable to the fixes — T17 still runs ~10 min).
- Every other test byte-identical; all 9 BI + 4 NUMERIC intact.

Interpretation: the two systematic fixes are correct and faithful, but their *measurable* line-count impact concentrates on T42 — `_THIRD` only crosses a 7-sigfig boundary in rare NAPS=0 cases (T49, whose only tape is the acer-stub tape71), and the lunion mixed-law over-shading was prominently triggered by Zn-67 MT107's log-log/histogram/log-log law. No test flipped status, but T42's quality jumped sharply and the regression surface is clean.

---

## Next work (handoff)

- **`NJOY_jl-2k4`** (P3) — MF12/13/14 photon-production pass-through (T45).
- **`NJOY_jl-c9f`** (P3) — thermr iel MT222 trace (probe T69; likely a real upstream bug or a documented libm floor).
- Re-check T28 (Pu-241) and other resolved-resonance tests' line counts post-k1v (k1v is systematic; the sweep will show which improved).
- The T49 HEAD-record `0 0` vs `0 99` metadata diff is a separate, pre-existing issue (noted in the 4k2 commit).
