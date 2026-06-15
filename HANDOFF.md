# NJOY.jl Session Handoff

> Historical detail — full session history, the Phase-18 sweep dump, the superseded Phase-33 'Current State', and dated session logs — has been moved to `handoff_archive/` (see `handoff_archive/INDEX.md`). This file holds only the live state.

---


## What is this project

NJOY.jl is a Julia port of NJOY2016 — the standard nuclear data processing system used worldwide for reactor physics, criticality safety, and radiation transport. The original is ~120 k lines of Fortran 90 code (`njoy-reference/src/`, 39 files; ~100 k code-only). Our Julia version is ~33 k lines of code (`src/`, 105 files, cloc code-only).

The goal: produce **bit-identical** PENDF/ACE output matching NJOY's own reference test problems (86 test directories in the cloned `develop` tree). The port must be idiomatic Julia — composable, differentiable, no global state — not a transliteration.

**Repo:** https://github.com/tobiasosborne/NJOY.jl (GPL-3.0)

---

## MANDATORY RULES — READ THESE FIRST

### Rule 1: Pursue Bit Agreement, Accept 1e-7

**Stretch goal (1e-9)**: Byte-for-byte identical MF3 output (columns 1-66 of each data line) with the Fortran NJOY2016 reference. This IS achievable — 19 RECONR tests already pass bit-identical. Pursue relentlessly for all modules.

**First-round acceptance (1e-7)**: Values agreeing within ±1 in the last digit of 7-sigfig ENDF format. This is the cross-compiler precision floor — even Fortran-to-Fortran fails at 1e-9 across architectures (25% failure on ifort, 30% on ARM64; see `reports/ACCEPTANCE_CRITERIA.md` for maintainer quotes).

**Structural match is non-negotiable at ANY tolerance**: Same line counts, same sections, same energy grid sizes. The NJOY maintainers explicitly flag line count differences as the real concern.

See **[reports/ACCEPTANCE_CRITERIA.md](reports/ACCEPTANCE_CRITERIA.md)** for the full tolerance hierarchy with verbatim quotes from NJOY maintainers (Wim Haeck, Jeremy Conlin) and published inter-code comparison standards.

### Rule 2: Fortran is Canonical Truth

If a unit test breaks after a change that makes output match Fortran, **fix the test** — the Fortran is ground truth. Never "fix" code to match a test if the Fortran disagrees.

### Rule 3: Read the Fortran Before Writing Julia

The Fortran source in `njoy-reference/src/` is the authoritative reference. Every formula, constant, rounding step, and accumulation order matters. Read the specific subroutine before implementing or fixing Julia code. Copy nothing blindly from the HANDOFF — verify against the source.

### Rule 4: No Parallel Julia Processes

Precompilation cache corruption is real and wastes hours. **Always** run:
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
```
before every test. One Julia process at a time. Weird runtime errors after agent work = cache corruption, not a code bug. **Use subagents for research (Fortran reading, ENDF parsing, code exploration) but NEVER for running Julia tests in parallel.**

### Rule 5: Use Bundled ENDF Files Only

All test data lives in `njoy-reference/tests/resources/`. Never download or use external ENDF files.

### Rule 6: Be Skeptical of Everything

Verify claims from this HANDOFF, from agents, from comments in the code. The codebase has had multiple agents working on it. Previous "fixes" have been wrong. Previous HANDOFF analysis has been wrong (e.g., "coincident breakpoint shading" was misidentified — the real cause was MF=12 histogram shading). Read the Fortran, test the claim, then proceed.

### Rule 7: Idiomatic Julia

No Fortran transliterations. Use multiple dispatch, broadcasting, clean types. But when the Fortran does something specific (like rounding to 8 sigfigs before adding potential scattering), match it exactly — just express it idiomatically.

### Rule 8: Use Subagents Wisely

Launch multiple **research** subagents in parallel for maximum efficiency:
- **Explore agents** to read Fortran source, search for patterns, understand data structures
- **Research agents** to analyze ENDF files, extract breakpoints, check data formats
- **Comparison agents** to run a SINGLE Julia diagnostic (grid extraction, MF3 comparison)

**NEVER** launch multiple agents that run Julia processes — Rule 4 applies. Only ONE agent should execute Julia at a time. Others should do read-only research.

---

## Beads Status (2026-05-28) — **RESTORED (fresh, empty)**

**Beads was wiped and re-initialized on 2026-05-28.** The previously broken
setup (external Dolt `sql-server` mode that lost its `NJOY_jl` database on a
server restart) was deleted entirely — `.beads/` removed from git and disk —
and `bd` was updated **v1.0.0 → v1.0.4**, then re-init'd with `bd init`.

The new database runs in **embedded Dolt mode** (in-process engine, no
external server — this removes the exact failure mode that broke the old DB).
Issue prefix is now `NJOY_jl`; IDs look like `NJOY_jl-<hash>` (e.g.
`NJOY_jl-a3f2dd`). The DB **starts empty**.

**Implication**: every historical `NJOY.jl-xxx` bead ID referenced in this
HANDOFF or any `worklog/*.md` is gone for good — those IDs will never resolve.
The work they tracked is re-captured in full in the "Open Work" section below,
which remains the authoritative work list. Re-file open items as fresh beads
(`bd create …`) as you pick them up.

**Sync model**: no Dolt remote is configured. Durable cross-machine state is
`.beads/issues.jsonl` (auto-exported, git-tracked); beads' git hooks
(`core.hooksPath=.beads/hooks`) sync Dolt + JSONL on commit/push, so a normal
`git push` carries beads state. There is no `bd dolt push` target.

Session rules:
- `bd ready` / `bd show` / `bd create` / `bd remember` all work again — use them.
- `bd remember` (beads' store) is separate from the harness auto-memory
  (`MEMORY.md`); both coexist (see CLAUDE.md §"Issue tracking: Beads").

---

## Recent Phases (Authoritative — Bead-Independent)

Beads-DB outage means the canonical phase log is here + per-phase worklogs
under `worklog/T*.md`. Most-recent first.

| Phase | Date       | Topic | Outcome | Worklog |
|-------|------------|-------|---------|---------|
| 87    | 2026-06-15 | thermr SAB: restore `sigfig(emax,7,-1)` held point at thermal cutoff (MF3/MT221+222 NP 1484→1485) — bead `NJOY_jl-js4` | **T70 tape60 MF3/MT221+MT222 NP off-by-one FIXED.** MF3/MT222 now byte-identical to ref (**0 diffs**, was 4 incl. structural); MF3/MT221 structural + `>1e-5` diffs gone (14 residual ≤4.8e-7, FP-grind/c3q). Total 213800==ref. **Root cause:** `_thermr_sab!` appended cutoff grid `[emax, sigfig(emax,7,+1), 2e7]`, omitting `sigfig(emax,7,-1)` (held thermal value = coh grid's last point, sigcoh thermr.f90:1194; value also held at emax via clamp 393, steps to 0 at up*emax 394, etop in tpend 3211-3213; MT221 ix=3 + MT222 ix=4 share the `ex` grid). Fix = extracted `_append_emax_sentinels!` helper + the one line `push!(grid, round_sigfig(emax,7,-1))`. **Section line count (498) + Phase-86 directory NC unchanged** (2968/2970 vals both ceil to 495 data lines; `3+div(np+1,3)`=498 both). New RED→GREEN unit test `test_thermr_emax_sentinels.jl` 6/6; Phase-86 `test_thermr_mf3_nc.jl` 10/10. **T01 NUMERIC_PASS 32812/32962 + T09 BIT_IDENTICAL 1830/1830 preserved (independently re-run).** **Rule-2 catch:** the reference test's regex `line_equivalence` mis-reported "39% of the matrix off by 90%" (greedy float regex merges last data field w/ MAT number) — a fixed-column ENDF parse showed only ~673/213800 lines differ. Free-gas path UNCONFIRMED → bead `NJOY_jl-nm2`; MF3/MT1 ~1e-6 grind = broadr Doppler floor (Target B bead, P2). Orchestrated: RED capture + 3 sonnet research (Workflow) + 1 opus TDD agent; serial Julia per Rule 9. | `worklog/phase87_thermr_emax_sentinel.md` |
| 86    | 2026-06-14 | thermr T70 MF1/MT451 directory-NC off-by-one (MF3/MT221+222: 497→498) — bead `NJOY_jl-h61` | **T70 tape60 directory lines 69-70 now match ref (NC=498); first-diff moved line 69→223 (value-grind, c3q); total 213800==ref, match 128064→128066 (+2).** Orchestrated (3+1 read-only research fan-out sonnet + 1 opus TDD agent; serial Julia per Rule 9). Root cause: the directory-NC formula substituted a stale pre-sentinel `thermr_coh_ne` (~1482) for the actual emitted np (1485). **Fortran ground truth (LAW 2) OVERTURNED the bead's own `3+cld(np,3)` suggestion:** thermr `tpend` uses `nc=3+(ne+2)/3` with `ne=NP_written-1` (one `etop=20e6` sentinel is in the written NP via `scr(6)=ne+1` but not in `ne`; thermr.f90:3087/3198/3211) → faithful form `3+div(np+1,3)`, which UNDERCOUNTS the true record count by 1 when `np≡1 (mod 3)` — reproduced per the GROUND-TRUTH PRINCIPLE, NOT "fixed" to ceil. New helpers `_thermr_mf3_dir_nc`/`_reconr_mf3_dir_nc` (reconr else-branch keeps `3+cld`, reconr.f90:5115-5116). New unit test `test_thermr_mf3_nc.jl` 10/10 (incl. the `np≡1 mod 3` case). T01 NUMERIC_PASS 32812/32962 + T68 directory NC=129 (ref-tape+unit verified) preserved. `thermr_coh_ne` plumbing now vestigial → cleanup bead `NJOY_jl-czw`. | `worklog/phase86_h61_thermr_mf3_dir_nc.md` |
| 85    | 2026-06-14 | thermr calcem MF6/MT221 +26-line overshoot — `terp1` chord vs midpoint average | **T70 tape60 +26 structural OVERSHOOT eliminated; total 213826→213800 (==ref), MF6/MT221 185321→185295 (==ref), all 101 blocks match, match 125483→128064; T01 NUMERIC 32812/32962 unchanged.** Root cause: calcem's E′-adaptive convergence (thermr.f90:2041-2068) tests `sigl(xm)` against the reference value `ym=terp1(x(i),y(k,i),x(i-1),y(k,i-1),xm,2)` — a **law-2 lin-lin chord at the SIGFIG-ROUNDED `xm`** (xm=sigfig(½(x(i-1)+x(i)),8,0), thermr.f90:2051-2052; terp1 endf.f90:1614). Julia used the simple average `½(y_lo+y_hi)`; since `xm` is sigfig-rounded it is not the exact midpoint, so the average shifted add/skip decisions in densely-refined high-E blocks → +7 E′-entries (mixed sign, blocks 78-94) → +26 data lines. Fixed both S(α,β)+free-gas calcem paths in `thermr.jl` with a `terp1_lin` helper. **Both research hypotheses were WRONG** (Rule 2): H2 `half_tol=0.5·tol` is FAITHFUL (sigl itself halves: `tol=half*tolin`, thermr.f90:2694); H1 last-seed zeroing would give uniform +1, but the per-block NW diff was mixed-sign — that diff decided it. **Orchestrator independent re-count caught a mis-report:** the inherited "MF3/MT221 off-by-one (497 vs 498)" is NOT a grid bug — the MF3 sections ARE 498; the residual first-diff (reference_test lines 69-70) is a **MF1/MT451 directory-NC** off-by-one (Julia writes NC=497) → reframed as bead `NJOY_jl-h61`. Beads: `7dp` closed, `h61` new (both under `c3q`). | `worklog/phase85_thermr_calcem_terp1_convergence.md` |
| 84    | 2026-06-09 | Oracle realign (ghost ac5adf5→pin 2c64dfb) + thermr MF6/MT221 cosine-width fix | **T70 tape60 112818→213826 lines (ref 213800), match 22%→59%; T01 NUMERIC + T09 BI preserved.** (1) **Oracle was a ghost**: `njoy-reference` HEAD `ac5adf5` (rebased-away 2026-05-30 baseline), not the pin. `setup_reference.sh` failed (`reference is not a tree`) because a plain `git fetch origin` can't pull a pin rebased off all branch tips — **fixed** by adding a direct `git fetch origin <SHA>` (bead `NJOY_jl-rby` closed). (2) **Rule-2 catch**: `c3q`'s "MF3 grid over-produced 213800 vs 112817" premise was INVERTED — at the pin, ref tape60=213800, **Julia=112818 (under-production)**; the bulk is the **MF6/MT221 thermal scattering MATRIX**, not an MF3 grid. (3) **Root cause** (empirically exact): MF6/MT221 stored 10 values/secondary (E′,σ,8 cosines) hardcoded instead of `nbin+2`; for `nbin=20`, 10/22=0.4545 = observed 0.455 deficit. Same 101 incident × 460 secondary energies — only the per-row width was wrong. (4) **Fix**: `thermr.jl` new `_mf6_entry` helper (`[E′,σ,μ₁…μ_nbin]`, len nbin+2) + `pendf_writer.jl` derives `nl=length(entries[1])` for LIST + directory-NC (Fortran `thermr.f90:1630/2230`). `nbin=8` byte-identical → T01/T09/T11 protected by construction. Corrects the thermal matrix for **all 10 `nbin=20` tests** (T24/25/32/49/67/68/69/70/72/74; the sweep "TIMEOUTs" T67/68/74 were the 300s cap+contention, not hangs — T68=265s, T70=143s). c3q stays in_progress (MF3 off-by-one, +26 MF6 overshoot, value grind, tape55/71/76 via `acer iopt=2` stub `p9q`). | `worklog/phase84_thermr_mf6_cosine_width.md` |
| 83    | 2026-06-04 | Orchestrated bead session: o2l (T72 crash) + lho (T09 leapr discrete oscillators) + 3rp close | **+1 BIT_IDENTICAL (T09), +1 NUMERIC_PASS (T33), T72 CRASH→DIFFS.** (1) **o2l**: T72 aplots `log10(-1e-6)` crash root-caused to `read_mf6_incident_energies` mis-parsing MF6 **LAW=7** (nested TAB2 — outer NE × inner NMU TAB1) as flat LAW=1, harvesting μ-cosines (−1.0,−0.9) as incident energies → negative ESZ grid. Fixed LAW=7/6 dispatch (`readers.jl`, Ref acefc.f90:6363-6366,6455-6462 + ENDF-6 §6.2.7); `_ascll` left as faithful Fortran port. T72 CRASH→DIFFS; aplots BI cohort T50/52/53/61/62 + T02/T08 verified clean (commit `c3bfd75`). Follow-up bead NJOY_jl-i23 (T72 ESZ structural gap 1711 vs 50748). (2) **lho**: T09 leapr **DIFFS→BIT_IDENTICAL 1830/1830** + **T33 DIFFS→NUMERIC_PASS** (75%@1e-9→99.996%@1e-5). Secondary scatterer was a RED HERRING (b7=1 free-gas → leapr.f90:398 skips secondary loop; primary-only correct). Real bug: **discrete oscillators never wired** — `leapr_module` called `generate_sab` without oscillators and ran no `discre` analog, silently dropping nd>0 deltas (T09: 2 oscillators, free-gas trans → off up to 565% in the alpha tail). Ported faithful `discre!` (leapr.f90:1320-1661) + `_bfact_leapr` (A&S-polynomial Bessel, leapr.f90:1663-1796) on the column-major ssm[β,α,T] layout; narrowed the misleading nss warning to b7≤0. T22 BI + T80 NUMERIC preserved; T23 (BeO nss=1/b7=0) DIFFS→STRUCTURAL_FAIL (oscillators apply but b7=0 SCT-secondary merge still unported → folded into NJOY_jl-gm1) (commit `5132e5b`). (3) **3rp** closed (suffix 80c→00c already in `0da9fa2`). (4) **c3q part 1** (commit `cf8e66a`): T70 thermr **MF1/MT451 faithful ENDF-6/5 header** — `write_full_pendf` emitted a corrupt 3-record header (LRP=0/NFOR=0/missing AWI-EMAX-NSUB-NVER); new `read_mf1_header_info` + `_write_full_mf1` hdr path reproduce it (reconr.f90:5028-5067). T70 records 1-8 byte-identical, T01 header fixed (NUMERIC preserved), T03/T09 BI preserved, T02/T08 untouched. **Triage correction:** T70 grid over-production (Julia 213800 vs ref 112818, ~1.9×) is REAL — the 2026-06-02 "premise stale" note predates the oracle pin alignment; that's c3q part 2 (calcem inelastic + lat=10 coherent grid, DEEP, open). Orchestration: Opus coders (think-hard/ultrathink), parallel read-only Sonnet researchers (3+1), serial Julia, verified-before-trust (parent re-ran T01/T02/T03/T09/T22/T23/T33/T70). | `worklog/phase83_o2l_lho_orchestrated.md` |
| 80    | 2026-05-31 | ACER aplots remaining blocks (heating + aplotr threshold + aploxp) — T53+T62 full BIT_IDENTICAL | **2 NEW FULL BIT_IDENTICAL: T53, T62.** Orchestrated session (Opus coders, parallel Sonnet researchers, serial Julia, verified-before-trust). Continues Phase 79: the scaffold's `AplotsNotPortedError` guards were already at the correct page-sequence slots, so this replaced guards with byte-faithful emission. Ported (all `acefc.f90`): (1) **aplotr threshold-reaction page** (16970-17093) — T62's sole block + T53 page 5; `1 0 2 1/` literal axes (no xtag/ytag), `_ascle` linear, MT filter excludes fission/inelastic/sub-threshold; (2) **log-log heating** (16191-16226) — `hmin=1e-10` floor, `ymin<ymax/scale` (scale=1e6), `_ascll`, `4 0 2 1/`; (3) **lin-lin heating** (16557-16595) — e≥0.2 threshold, always-write-last; (4) **aploxp** (new `_aploxp`, 18787-19460, gated `nxs[NXS_NTYPE]>0`) — C1 particle-heating-contributions (`xss[hpd+1+naa+i]`), C2 recoil heating (esz_heat−Σparticle), C3 production XS (`xss[hpd+1+i]`, photons/5), **C5 3D angular pages reuse `_aplof4`'s perspective body** (factored into `_aplof4_perspective!`; subtitle `angular distribution for <len_trim(name)> <particle-word>` e.g. `(d,p*0) triton`). **T62 tape33 5460/5460 + T53 tape33 17236/17236 cmp-clean** (independently re-verified). Fortran surprises: (a) **`1p,e14.6` drops the `E` for 3-digit exponents** (`1.155100-117` not `...E-117`) — fixed `_e14_6`, untouched for 2-digit (T50/T52 BI pages never hit this); (b) C5 subtitle strips name trailing blanks (unlike `_aplof4` elastic); (c) incident-particle `name(2:2)` swap is caller-side (`_aplots_name_izai`); (d) mtname table extended 50→500 entries from `acecm.f90:31-138`. **Zero regressions:** T50 432/432, T52 6042/6042, T61 49211/49211 BIT_IDENTICAL preserved. **T54** tape33 now structurally complete (11340/11340 lines, was empty stub) — DIFFS trace 1-ULP to its tape34 ACE values (bead NJOY_jl-53h); flips to FULL BI once that FP word lands. **T51/T71** aplots generalizes (no crash) but tape33 wrong-content because their ACE/tape34 is structurally wrong (pre-existing, ACE-level). Unported aplots sub-paths (MT=444 damage, URR, aplonu, aplodn, aplopf, aplof4 equiprobable, aploxp law-4/44/61 emission) stay loud-guarded; no current test hits them. | `worklog/phase80_acer_aplots_t53_t62.md` |
| 79    | 2026-05-30 | ACER aplots viewr plot tape (iopt=7 tape33) — T50+T52 full BIT_IDENTICAL | **2 NEW full BIT_IDENTICAL: T50, T52.** Orchestrated session (Opus coder, 3 parallel Sonnet researchers; serial Julia; verified-before-trust). Ported Fortran `aplots` (acefc.f90:15212-16830) + delegates aplotr/aplopf/aplof4 (16832-17428) + ascll/ascle axis scaling (19541-19724). New `src/formats/ace_reader.jl` (`read_ace_ascii` inverts the Type-1 ASCII writer, 13992-14013) feeds new `src/formats/ace_aplots.jl`. Wired into `_acer_iopt7`, gated on the ZAID class letter c/h/o/r/s/a per `acefix` (14198-14206) so thermal ('t') emits no plot. **T50 tape33 432/432 + T52 tape33 6042/6042 BYTE-IDENTICAL** (independently re-verified cache-nuked). Fortran surprises: (1) legend-tag coords use `e12.4` WITHOUT `1p` (the `0.ddddE±dd` form) unlike every other plot field's `1p` → custom formatter; (2) aplof4 elastic name is `character(10)='elastic'`="elastic   " (3 trailing blanks, all 10 emitted); (3) T50 has ntr=0/gpd=0/iurpt=0 so only the two principal-XS pages + the elastic angular-distribution page fire — all guarded blocks correctly skip. **Zero regressions:** T61 BIT_IDENTICAL preserved (thermal gate), T62/T53 tape34+tape35 BIT_IDENTICAL preserved — their richer tape33 blocks (heating / aplotr threshold / aploxp particle-production) raise a dedicated `AplotsNotPortedError` caught to leave the prior empty stub + still write tape35; any OTHER exception aborts loud (Rule 6). Remaining blocks → bead NJOY_jl-zsh (P2). | `worklog/phase79_acer_aplots.md` |
| 78    | 2026-05-29 | ACER charged-particle arc (ptlegc/coul, NTR>0 reactions, iopt=7 rename+xsdir, acelcp particle-type) + thermr lat=10 | **4 NEW BIT_IDENTICAL: T52, T62, T61, T53.** Orchestrated subagent session (Opus coders, Sonnet scopers; serial Julia; verified-before-commit; 5 commits pushed). **(1) ptlegc/coul** (acefc.f90:7980-8201) for LTP=1/2 Coulomb-Legendre charged elastic → **T52 tape34 BI 3986/3986** (root cause: `nt=nint(c(6))` is LIST CONT N2=NL, not data). **(2) thermr lat=10**: new `read_mf7_mt2` (MF7/MT2 LTHR=1 coherent elastic, sigcoh transform thermr.f90:1150-1166, icoh=10·lthr dispatch) + **latent bug fix** (`read_mf7_mt4` read natom from B(4)=emax not the input card, thermr.f90:1670-1676; T01 byte-identical). T70 lat=10 grid over-produces (213800 vs 112817)→bead c3q; iel/LTHR=2→69j. **(3) charged NTR>0 path** → **T62 tape34 BI 7221/7221** (5 deep bugs: MF3 INT=6 Coulomb interp terp1 i=6; disappear=Σ sigfig7 reaction xs acefc.f90:5648; heating zeroed izai≤2004; LQR=QI/1e6; a11 e20.11 drops E for \|exp\|≥100). **(4) iopt=7** suffix rename (newsuff acecm.f90:619) + class-aware xsdir (thermal aceth.f90:2422 vs neutron acefc.f90:13013) → **T61 BI** + **T50/T52/T62 tape35 xsdir 0/1→1/1 BI**. **(5) acelcp particle-type** (new ace_lcp.jl/ace_lcp_build.jl; PTYPE/NTRO/PLOCT + per-type HPD/SIG/AND/DLW/YH, MF6 LAW=2/3/4 + MT=2 recoil; acefc.f90:9121-11056) → **T53 tape34 BI 12030/12030** + **T54 STRUCTURAL_FAIL→DIFFS** (7327/7327 structure, 6383 match, tape35 BI). 2 prereq fixes (charged union grid unionx acefc.f90:1541-1693 →140 pts; SIG threshold ie_start acefc.f90:5594-5628), gated to charged path. Unported MF6 LAWs (1/6/7/Kalbach/phase-space) @warn+skip → **T51 LAW=6 crash fixed** (MISSING→1197). **Zero regressions:** T01 NUMERIC_PASS 32812/32962, T02 tape28 NUMERIC_PASS, T08 neutron, T50/T62/T22 BI all preserved/verified. **`reports/REFERENCE_SWEEP.md` is now STALE — run a fresh sweep.** Next: aplots viewr plot tape (tape33, bead cnh.4) blocks FULL pass on all 4 charged tests. | `worklog/phase78_acer_charged_arc.md` |
| 77    | 2026-05-28 | broadr O(n²)→O(n) perf + leapr/ACER/Bragg premise corrections + fresh sweep baseline | **broadr U-238 1039s→40s (~26×), byte-identical.** Root cause: `sigma1_at` rebuilt the full velocity array (`v=map(e->sqrt(alpha*e),seg_e)`) every call → O(nseg²); `broadn_grid` already had `vel` at line 79. Fix threads precomputed `vel` through `sigma1_at` + all 3 call sites; `ntuple(Val(5))` in `h_all` kills a per-interval closure alloc. Verified byte-identical (T08 tape23 md5 unchanged; T01 tape25 32812/32962). **This unblocked T15/T16/T65 CRASH→DIFFS** (broadr no longer the bottleneck; Phase 71-76 covcal output now measurable: T15 tape26 5089/5958, tape35 3724/3812). **Fresh full sweep (post-Phase-76, 80 min): 2 BIT_IDENTICAL (T03,T22), 2 NUMERIC_PASS (T01,T80), 78 DIFFS, 1 CRASH (T20), 1 NO_REFERENCE.** Crashes 4→1 vs stale Phase-67 sweep. **Premise corrections (3 stale-HANDOFF beads caught by Fortran verification):** (a) leapr `naint` gating is a NON-bug — `naint` hardcoded to 1 (leapr.f90:292-294, never read), Julia already matches; (b) ACER "8 siblings" overstated — T50 is 143/163 BI with 20 *real* path-diff lines (~1e-6), siblings T51-54/T62/T71 each need separate MF6 work; (c) Bragg bead reframed — keep hardcoded lat=1/2/3 (faithful), ADD the lat=10 MF7/MT2 read path (sigcoh thermr.f90:1016). **Regression found + fixed:** T20 CRASH (`read_mf32` aborted on LRF=7, introduced by Phase 72-76 rescon) → `read_mf32` now @warn+skips unported MF32 sub-formats (LRF=4/7, LCOMP≠1, NRO≠0, NLRS≠0, ISR=1+LRF=3, LRU≠1); genuine structural violations still error (Rule 6). T20 CRASH→DIFFS; supported LRF=3+URR byte-unchanged (mf32_reader 83/83, covcal_lb5 50/50). **Crashes now 0/84.** Proper rpxsamm port = bead NJOY_jl-4pz. **ACER T50 BIT_IDENTICAL:** tape34 143/163 → **163/163 @1e-9** — acecpe reads sigfig-rounded ESZ columns (xelas/signi from sigfig7 elastic; total=sigfig9(Σ sigfig7(reaction)) over MT=2/MT≥5; `_log_log_interp` rewritten to mirror terp1 int=5). Charged-path-only change; T01/T02 ACE unchanged. **NOTE:** `reports/REFERENCE_SWEEP.md` was generated BEFORE the T20 fix + ACER win, so its headline (2 BIT_IDENTICAL, 1 CRASH) is now **3 BIT_IDENTICAL (T03,T22,T50), 0 CRASH** — next sweep will reflect this. ACER siblings T51-54/T62/T71 are NOT free riders (each needs separate MF6 work). | `worklog/phase77_broadr_acer_sweep.md` |
| 74    | 2026-05-18 | NC cross-pair to lower-MT (T15 tape26 line gap −68 → **−6**, MT=1/mt2=2 bit-identical on geometry) | **Single-NC cross-pair pass now handles `ref_j < mt`.** The Phase 48 `_expand_nc_blocks!` filter `ref_j > mt || continue` skipped MT=2's NC reference to MT=1 (since 1 < 2). U-238 JENDL's MF=33 MT=2 carries `MT=2 = +1·MT=1 − Σ_partials c_i·MT_i` (J33U238 lines 16555-16568). By cov symmetry the cross-pair `Cov(MT=1, MT=2) = +1 · Cov(MT=1, MT=1) · σ(MT=1, igp)/σ(MT=2, igp)` (σ-ratio on the **column** side — the transpose of the existing row-side ratio). Fix: split into `cross_acc_hi` (`ref_j > mt`, store at `(mt, ref_j)` with row σ-ratio) + `cross_acc_lo` (`ref_j < mt`, store at `(ref_j, mt)` with column σ-ratio). Lower branch uses `+=` accumulation (matches Fortran `cova` semantics across iy/iyp passes). Mirrors errorr.f90:7406-7438 sandwich loop (`akxy(iy, ix, k) · akxy(iyp, ixp, kp) · cov(jgp)`) with the asymmetric storage convention `cov_matrices[(min, max)]`. Results: T15 tape26 **5890 → 5952 lines** (gap **−68 → −6** vs ref 5958); MT=1/mt2=2 sub-section **3 lines (stub) → 65 lines (16 rows × 3 data lines)** — bit-identical to ref on geometry. Phase 72c MT=102 C[1,1] canary preserved; T22 BIT_IDENTICAL 4636/4636 preserved. **Tests: 223 assertions PASS, zero regressions** across mf33_sparse 75/75 (2 new MT=1/mt2=2 assertions), covcal_lb5 50/50, writer_mf_dispatch 51/51, nc_expansion 9/9, gendf_readback 38/38. Remaining −6 line gap: MT=102 −2 (rescon row 14 edge), MT=2 −2 (rescon row 13-15 span), MT=1 ~−2 (single-row spans). | `worklog/phase74_mt1_nc_lower_ref.md` |
| 73    | 2026-05-16 | T15 MF=33 cross-material skip + NC σ-ratio weighting (T15 tape26 line gap +150 → −68) | **Two stacked bugs in the MF=33 covariance pipeline.** (1) Cross-material sub-section pollution: U-238 MT=18 (fission) has NK=6 sub-sections — 1 self + 5 cross-MAT refs (`L1=MAT1 ∈ {9222, 9228, 9437, 9440, 9443}`). Julia's reader ignored `sh.L1` and pushed cross-material LB=6 data into `pair_blocks[(18,18)]`, polluting MT=18 self with up to 10⁷ junk values via LB=6 NEC mis-parse. Fix: detect `mat1 = sh.L1 != 0`, skip pushes (still read records to advance file position) — mirrors Fortran rescon `if (mats(ixp).ne.0) return` (errorr.f90:8531). (2) NC σ-ratio missing: `_expand_nc_blocks!` did `self_acc += c_i² · rel_cov(ref_i)` instead of `c_i² · rel_cov(ref_i) · (σ_i/σ_Y)²`. At high-energy elastic cells (σ_Y=σ_elastic ~3b, σ_i=σ_inelastic_level ~0.1b), the missing factor caused 5621× over-estimate. Fix: thread `group_xs` into `_expand_nc_blocks!`, apply σ-ratio per cell with denom>0 guards. Self gets `(σ_i/σ_Y)²`, single-NC cross gets `σ_ref_j/σ_Y` once, double-NC pair gets `(σ_iz/σ_a)·(σ_iz/σ_b)`. Also (P3 cleanup): `rescon._build_subsection_sensitivities_rm` now zeroes `sens[*, *, ig]` outside `[iest, ieed]` per Fortran's `do ig=iest,ieed` write-bound (errorr.f90:4509) — functional impact on T15 nil but matches Fortran. Results: T15 tape26 **6108 → 5890 lines** (gap +150 → **−68** vs ref 5958); MT=2 Δ **+142 → −2**, MT=18 Δ **+74 → 0**, MT=1 Δ unchanged at −64 (cross-MT=1/2 separate bug), MT=102 Δ unchanged at −2 (row 14 separate). MT=2/mt2=2 row 28 col 28 \|Δ\| **+7.42 → +2.50e-6** (3M× improvement). MT=2/mt2=18 rows 15-30 now match within 1e-11. Phase 72c MT=102 C[1,1] canary preserved. **Tests: 222 assertions PASS, zero regressions**: covcal_lb5 50/50, mf33_sparse 73/73, writer_mf_dispatch 51/51, nc_expansion 9/9, gendf_readback 38/38; T22 BIT_IDENTICAL 4636/4636. | `worklog/phase73_nc_xmat_and_sigma_ratio.md` |
| 72c   | 2026-05-16 | rescon iwt=6 weight fix — Fortran `rpxgrp` + `egtwtf` port closes T15 MT=102 C[1,1] | **Phase 71 RED canary GREEN.** `_group_average_flat` (iwt=2 trapezoidal) replaced with `_rpxgrp_average` (rpxgrp port — errorr.f90:5227-5366), iwt-aware via `get_weight_function(iwt)` (egtwtf — errorr.f90:10023-10110). T15 input deck line 34 `9237 3 6 1 1 /` uses `iwt=6` (thermal + 1/E + fission + fusion); Phase 72b mistakenly assumed `iwt=2` (flat), biasing group-1 σ̄ by ~16% and C[1,1] = σ̄² × cov by ~35%. C[1,1] **3.580e-4 → 2.658408e-4** (ref 2.658914e-4, diff −5.06e-8, ~0.019% rel — within 1e-7 canary tolerance). `@test_broken` flipped to `@test`. `apply_rescon!` gained `iwt::Int` kwarg; errorr orchestrator threads `params.iwt` through. Fortran precedence quirk on errorr.f90:5338 (`gsig=yl*zl + yr*zr/sumde` parsed as `yl*zl + (yr*zr)/sumde`, not `(yl*zl+yr*zr)/sumde`) replicated verbatim. T15 MT=102 row-1 41/41 PASS (was 40 PASS + 1 BROKEN); T22 BIT_IDENTICAL 4636/4636 preserved. Residual −5e-8 plausibly Julia tanh grid vs Fortran rpendf eskip-mesh + a 6-digit-vs-5-digit `wt6b` constant in weight_functions.jl (3e-6 fractional, flagged for follow-up). | `worklog/phase72c_rescon_iwt_fix.md` |
| 72b   | 2026-05-08 | rescon sensitivity builder + sandwich (MT=102 row 1: all-zero → 10/10 nonzero, sign-matches ref, magnitudes within factor of 2) | **Phase 71 canary partially flipped GREEN.** `apply_rescon!` no-op skeleton replaced with full Fortran-faithful pipeline (~400 LOC delta in rescon.jl) mirroring `rpxlc12` (errorr.f90:4399-4593) + rpendf-equivalent grid (5015-5089) + rpxgrp-equivalent group-average (5227-5367) for LRU=1/LRF=3/LCOMP=1. Builds dense pointwise grid (200 pts/decade + tanh-stretched per-resonance refinement spanning ±30·Γ + output-group boundaries), evaluates `cross_section_rm` at each point for unperturbed + ±-perturbed RM parameters (ER ±1e-4, widths ±1e-2 — verbatim Fortran factors), central-differences pointwise dσ/dRP, group-averages with iwt=2 flat weight, sandwiches via triangular-or-full walk over MF=32 cov per (mt, mt2) pair, divides by σ̄·σ̄ for relative cov. T15 MT=102 row 1: jul=10/10 nonzero (was 0/10), sign-pattern matches ref for every nonzero col, ratios 0.76..1.69 (most within 35%), C[1,9]<0 negative-sandwich signature confirmed (flipped from @test_broken to @test). 1e-7 bound on C[1,1] still @test_broken (ratio 1.346 — needs FP-grind alignment with rpendf adaptive eskip mesh + rpxgrp accumulation order). T22 BIT_IDENTICAL preserved (4636/4636). T15 errorr tape26 6108 lines (was 5964; ref 5958). Test suite: 40 PASS + 1 BROKEN canary, 5 PASS Bug A NK, 57 PASS reader. Worklog: `worklog/phase72b_rescon_sandwich.md`. | `worklog/phase72b_rescon_sandwich.md` |
| 72    | 2026-05-08 | rescon foundation — MF=32 reader + apply_rescon! skeleton (LRU=1/LRF=1,2,3/LCOMP=1) | **Foundation landed; canary stays @test_broken by design.** New `read_mf32` (213 LOC) parses U-238 JENDL into typed structs (10 resolved ranges, 317 resonances, 951 uncertain parameters); URR (LRU=2) skipped cleanly. New `apply_rescon!` (rescon.jl, 86 LOC) parses + validates + logs but does NOT yet compute sensitivities. Wired into `errorr_module` replacing the Phase-71 `@warn`. `test_mf32_reader.jl` (NEW, 57 PASS assertions) exercises the reader against actual J33U238 ENDF data — top-level header + range 1 (1e-5..1000 eV, RM/LCOMP=1, NRB=26, MPAR=3) + cov symmetry/non-negative-diagonal across all 10 ranges + range 2 sanity. Phase 71 RED canary suite: 18 PASS + 2 BROKEN preserved. T15 errorr tape26 5964 lines unchanged. Next session: port `rpendf` (errorr.f90:5015-5089), `rpxgrp` (5227-5367), the rpxlc12 perturbation loop (4399-4523), and the sandwich fill (4541-4593) — then GREEN the `MT=102 C[1,1] ≈ 2.658914e-4` canary. | `worklog/phase72_mf32_reader.md` |
| 71    | 2026-05-07 | T15 covcal MT=102 row-1 diagnosis — rescon (MF=32→MF=33) is the missing piece | **Diagnosis only — no functional code change.** Matrix-level dump of T15 MT=102 self-cov reveals rows 1..14 zero in Julia, populated in ref (with negative entries proving sandwich-rule output, not MF=33 LB=5 lookup). Fortran covout (errorr.f90:7465) calls rescon (errorr.f90:8513-8819) to add MF=32 resonance-parameter cov contributions to seven (mt, mt2) pairs: (1,1), (2,2), (2,18), (2,102), (18,18), (18,102), (102,102). Disproves HANDOFF P1's "Bug-B siblings (LB=1/LB=2 midpoint sampling)" hypothesis as the dominant cause — the LB=5 weighted-collapse path (Phase 51) produces *exact* values for output rows that fall in covered input bins (MT=102 rows 16..30 are byte-identical between Julia and ref). The drift is missing rescon. Landed: MF=32 detection scaffold + loud-warn at the would-be rescon call site (errorr.jl, no functional change), Phase 71 RED testset in `test_errorr_covcal_lb5.jl` asserting `MT=102 C[1,1] ≈ 2.658914e-4` (currently `@test_broken`), HANDOFF P1 Covcal section updated to point at rescon. Multi-session port estimate: ~1500-2500 LOC across 8+ Fortran subroutines (resprx, rpxlc0/12/2, rpxsamm, rpxunr, rpendf, grpav4, rescon). T15/T22/T01 unit-test regression-clean. | `worklog/phase71_rescon_diagnosis.md` |
| 70    | 2026-05-06 | musigc MT=251 derivation + per-mfcov MF=3 echo restriction (T15+T16+T65 CRASH→runs, sweep CRASH count 4→0) | See worklog. | `worklog/phase70_musigc_mt251_derivation.md` |
| 69    | 2026-05-05 | multi-temperature groupr port (T11 CRASH→DIFFS, last CRASH from Phase 67 cleared) | **All 4 Phase-67 CRASHes now resolved.** T11 (Pu-238 WIMS, ntemp=3 at 300/900/2100K) was crashing in `wimsr_extract_resint` because Julia groupr only emitted one temperature, so wimsr's GENDF reader (`_wimsr_read_gendf_metadata`, wimsr_xsecs.jl:188-248) saw `tempr` of length 1 and threw `BoundsError` on `tempr[1:ires]` for ires=3. Three layered fixes: (1) `_pendf_material_temperature` + `extract_mf3_at_temperature` added to `pendf_io.jl` (+75 LOC) — heuristic locates TEMP CONT in MF1/MT451 by walking lines 2..5 picking last (L2=0, N1>0, 0≤C1<1e6) match; integer-fields-parse-first guard skips description Hollerith; handles ENDF-IV (TEMP at mf1_lines[2], T11) and ENDF-V/VI (mf1_lines[3] or [4], T01); errors loudly if no match listing available temperatures. (2) `groupr_module` (groupr.jl) refactored to outer-T loop mirroring Fortran groupr.f90:479-943 — `temps = copy(params.temperatures)`; per-T `extract_mf3_at_temperature` + per-T mt_results buffer indexed `[ti][mti]`. (3) `_write_groupr_tape` signature changed to take `temps::Vector` + `sigz_list::Vector` + `per_temp_results::Vector{Vector{MTRes}}` — emits per-T MF=1/MT=451 (HEAD `(za, awr, 0, nz, -1, 1)` + LIST `(tempin, 0, ngn, 0, NW, 0)` ← C1=tempin per Fortran 551, was 0.0) + per-MT MF=3 (head with `nz_mf3=1` simplification + per-(g, ig2) records `(tempin, 0, 3, 1, 3, g)`) + MEND between temps + single TEND. The `nz=1` MF=3 simplification keeps body length ≥ wimsr's `_bidx(i, iz, ig2, nl, nz)` requirement; full nsigz×ngn block-matrix output (URR self-shielding via Fortran `genflx`) is a follow-up. Side fix: removed Phase-58a leftover `_collect_gaspr!(ctx, tapes, params)` call inside the wimsr branch in pipeline.jl:262 — gaspr collector is typed `GasprParams`, the wimsr branch passes `WimsrParams`, MethodError for every wimsr-using test. Verified single-test: T01 NUMERIC_PASS 32812/32962 preserved, T22 BIT_IDENTICAL 4636/4636 preserved, T11 CRASH→DIFFS (tape27 13/455 STRUCTURAL_FAIL ref 1169; tape28 0/10379 STRUCTURAL_FAIL ref 15722). Remaining T11 work to BIT_IDENTICAL: multi-σ₀ MF=3 emission (~40 LOC), Fortran-faithful MF=1 LIST body (titles+sigz+egn+egg), MF=6 transfer matrix port (substantive), MT=251/252 derivation (HANDOFF P2 cdy), moder tape28 TPID. | `worklog/phase69_groupr_multi_T.md` |
| 68    | 2026-05-04 | errorr body-MF dispatch (mfcov∈{31, 34} remap) + MT=1 PENDF unfilter (T15 + T16 + T65 CRASH→runs) | **3 of 4 CRASHes resolved.** `_write_errorr_tape` (errorr.jl) now mirrors Fortran `covout` (errorr.f90:7211-7587): mfcov=31 (ν̄) → body MF=33 (L7214/7472 remap), mfcov=34 (mubar) → ONE outer MF=34/MT=251 with sub-section CONT carrying L1=251/L2=ld/N1=ld1 (L7245-7250, 7480-7485, 7585), mfcov=33/35/40 unchanged. Mubar fix collapses per-reaction pairs into one section (Fortran emits N sections; Julia reader keys by (mf, mt) so multi-section emission needs reader-side multi-section support, deferred). MF=3 echo: removed `mt in (1, 451) && continue` → `mt == 451` only (PENDF path, errorr.jl:364) — Fortran `colaps` (errorr.f90:9097+) walks every MT incl. MT=1; T16 crash was Julia-side over-filtering. MF=3 emission iterates `reaction_mts ∪ {251 if mfcov==34}` (NOT all of group_xs — first-iteration `keys(group_xs)` regressed T34 covr boxer-format output 2→1150 lines because GENDF readback puts every MF=3 MT into group_xs). 47-assertion unit test landed (`test/validation/test_errorr_writer_mf_dispatch.jl`) covering all 5 mfcov branches + T34 + T16 regression guards. T15 mfcov=33 tape26 dropped 5985→5964 lines (ref 5958, gap 27→6). T22 BIT_IDENTICAL preserved; T15/T16/T65 errorr no longer crashes covr (was: covard MF33/MT=452, MF3/MT=1, MF33/MT=2). T11 wimsr CRASH unchanged (orthogonal — multi-T groupr port). Deferred: mubar multi-section emission for T15 mubar tape27, MF=3/MT=251 musigc derivation for T16's mubar errorr call, TPID descriptor + iverf=6 + temperature in MF1/451. | `worklog/phase68_errorr_body_mf_dispatch.md` |
| 67    | 2026-05-03 | full 84-test sweep + errorr/covr `mfflg` fix + wimsr loud diagnosis | See worklog. | `worklog/phase67_full_sweep_post_phase66.md` |
| 66    | 2026-05-03 | powr lib=1 delayed-neutron spectra (MF=5/MT=455, U-235 nfs=6 oracle) | **BIT-IDENTICAL on first run.** 131-line / 10611-byte tape50 byte-for-byte. Three landings in `src/orchestration/modules/powr.jl` (~140 LOC delta): (1) Reader extended with `mf5_455`/`mf3_455` containers; gamll's `nfs += nl` accumulator now correctly reads NL from the MF=5/MT=455 HEAD's L1 (was previously a broken stub reading NG2 from the LIST CONT). (2) Pointer layout corrected: `locdla = locchi + ngnd` between locchi and locnus, sized `ngnd*nfs` (Phase 64 worklog had this wrong). (3) Delayed-chi accumulation mirrors gamxs.f90:1226-1262 — MF=3/MT=455 branch sets dnu/dnorm/jgdnu (defensive — not exercised by our oracle), MF=5/MT=455 branch accumulates `dnorm`-weighted prompt-chi + `dnu`-weighted delayed-chi + the `jgc=ngnd` OVERWRITE from `a[jgt+locb]` (quirk: `locb=6` makes this read from the absorption region, copying capture XS into the last cell of each delayed-chi block). (4) Final-pass normalization mirrors gamxs.f90:1366-1393 with the `sumd[il]*rnorm` replacement at ig=ngnd-1. (5) Writer emits `nfs` delayed-chi blocks (each `(i6,i2,10a4)` header with `i2=ifs` + ngnd-cell chi data) after the prompt chi. Loud-TODO at `nfs>0` removed. All 6 powr lib=1 standalone tests pass (carbon abs-only, +MF6 elastic, +inelastic, U-235 fission, Fe-56 self-shielding, U-235 delayed neutrons). T01 regression-clean at NUMERIC_PASS 32812/32962 @ 1e-5. **powr lib=1 now ~98% covered**; only the dead `matd<0`/iread=1 corner remains. | `worklog/phase66_powr_lib1_delayed.md` |
| 64+65 | 2026-05-03 | powr lib=1 multi-temperature loop + gamff self-shielding f-factors (Fe-56 ntemp=3 nsigz=4 oracle) | **BIT-IDENTICAL on first run.** 177×80-char tape50 byte-for-byte (14337/14337 bytes). Combined Phase 64 + 65 because the multi-temp output only enters via the `gamff` f-factor block at powr.f90:614-641 (gated `nff!=0 && iwr!=0` ≡ `nsigz>1 && iff=1`). Three substantive landings: (1) `_powr_read_gendf_for_fast` rewritten to walk all temps, return `tmpr[]`/`sigz[]`/`rtemp_idx`, capture full per-σ₀ XS vectors per (mt, ig); (2) `_powr_pack_ffactors` (~110 LOC) directly mirrors `gamff` (powr.f90:1401-1542) — per-temp absorption (MT 102-107) + fission (MT 18-21,38) accumulation in flat `[1..nsigz, 1..ntp, 1..ngnd]` column-major array, log-formula final pass `log(σ_z) - 2·log(ratio)`, `iwr = ngmax - iglo + 1`; (3) `_powr_write_ffactor_block` emits `(4i6) nsigz/ntp/misc/jwf` + sigz + tmpr + abs f-factors over groups [iglo, ngmax] + dummy `iwr*10` zeros. Latent bug fixed: previous reader's `xs = data[nl+1]` happened to match the correct multi-σ₀ index `data[nl*nz+1]` only when nz=1 (carbon, U-235); Fe-56 nz=4 exposed it as Julia reading flux values where it expected XS. Five powr lib=1 standalone tests now pass (carbon abs-only, +MF6 elastic, +inelastic, U-235 fission, Fe-56 self-shielding). T01 regression-clean at NUMERIC_PASS 32812/32962 @ 1e-5. | `worklog/phase64_65_powr_lib1_selfshielding.md` |
| 63    | 2026-05-02 | powr lib=1 fission path (MF=3/MT=18 sigma_f + MF=6/MT=18 chi/nu, U-235 oracle) | **BIT-IDENTICAL on U-235.** 53 × 80-char tape50 byte-for-byte. Switched powr to a unified `a[]` array with Fortran-style pointers (locab0/locsf0/locchi/locnus) so the gamxs spectrum-branch quirk (`jg2c = ngnd - ig2c - 1`, sign mismatch vs matrix branch's `+1`) replicates exactly — Fortran silently corrupts the last 2 sigf cells via `locchi-1/-2`, and we reproduce it. Fission MF=3 (MT=18 latch + redundant MT∈{19,20,21,38} skip), MF=6/MT=18 (cspc setup at ig=0 + spectrum-distribution at ig>0 with ig2lo=0), nu = nus/sigf, chi /= cnorm. Bug fix in GENDF reader: MF=6 ig=0 records were being filtered (relaxed to ig≥0). Output writer emits chi block first (kscr equivalent) then sigf+nu before matrices. | (commit ea6da24) |
| 61+62 | 2026-05-02 | powr lib=1 MF=6 elastic + inelastic + n2n matrices | **BIT-IDENTICAL on carbon** for both: 84 lines (elastic) and 113 lines (elastic+inelastic). Factored matrix accumulation into `_powr_pack_matrix(kind ∈ {:elastic, :inelastic, :n2n})`. Elastic uses Legendre P0+P1 with sigma-zero index izref and 3×/cfluxp1 multipliers for P1 (gamxs.f90:1311-1318). Inelastic + n2n use only first sigma-zero P0; n2n divides accumulated value by 2. Output emission order matches fast() lines 572-613: elastic P0, elastic P1, inelastic, n2n. | (commit 96cb744) |
| 60    | 2026-05-02 | powr lib=1 (fast/GAMTAP) carbon-only — replace Phase A scaffold with real implementation | **BIT-IDENTICAL on first try** against Fortran NJOY's tape50 (carbon 68-group GAM-I, full reconr→broadr→groupr→powr chain). 16/16 lines × 80 chars match exactly. parse_powr extended to consume cards 3-5 for lib=1 (per-material matd/rtemp/iff/nsgz/izref + word + fsn). Module orchestrates GENDF walk → cflux normalization → absorption accumulation (MT ∈ [102, 150]) → 4 header lines + 12 data lines via `_powr_e125` / `_powr_pad80` helpers. Multi-temp / MF6 matrices / fission / self-shielding / matd<0 paths all error LOUDLY with named TODOs (Rule 6). T01 regression-clean. Next: Phase 61+ tackles MF6 elastic matrix (Fe-56 oracle), fissionable path (U-235 oracle), then lib=2 / lib=3. | `worklog/phase60_powr_lib1_carbon.md` |
| 59    | 2026-05-02 | mixr + resxsr + powr scaffold — close 23/23 module wiring gap | **22/23 fully ported, 1/23 (powr) Phase A scaffold.** mixr (~280 LOC) BIT-IDENTICAL on two standalone oracles (carbon weight=1.0 + self-mix 0.5+0.5; 1059/1059 lines). resxsr (~340 LOC) BIT-IDENTICAL on first try (carbon @296K, eps=0.001; CCCC binary 7968/7968 bytes). powr (~50 LOC scaffold) parses cards 1-2, errors loudly per Rule 6 with mode-named TODO. Quirks pinned: mixr TPID is 66 blanks (not description); mixr writes no SEND between MT=451 and FEND; mixr dictionary entries use `(22x, 4i11)` not zero-floats; resxsr description goes nowhere (forced to blanks); resxsr thinner intentionally drops uncommitted stack contents at end of grid; resxsr stores `(a6)` material names left-justified in 8-byte slots. T01 regression-clean at NUMERIC_PASS 32812/32962. Next: powr Phases B-D (fast/therm/cpm). | `worklog/phase59_mixr_resxsr_powr_wiring.md` |
| 58a-d | 2026-05-02 | wimsr orchestration + xsecs + resint + multi-temp + p1flx (T11 Pu-238 WIMS) | **130/1169 lines bit-identical (11.1%) — was MISSING_TAPE.** Four-phase port of NJOY's wimsr module (wimsr.f90:48-2147, ~2150 LOC Fortran → ~1100 LOC Julia). Phases: 58a scaffold (commit `6eb9d7d`, 0→19); 58b xsecs port (`a93e5e7`, 19→55) — multi-temp dispatcher, MF=3 + MF=6 routing, ip1opt diagonal correction, l1/l2 zero-gate for 270 block, MT=2 thermal-skip, group-bound reverse; 58c resint port + writer concat fix (`d2bc264`, 55→103) — Goldstein-Cohen normalization for resonance integrals, transport+absorption emitted as one stream per Fortran nscr3 layout; 58d p1flx + nu-bar from fission matrix (`de90641`, 103→130) — MT=1 dispatch populates p1flx for csp1 weighting, snu derived from MF=6/MT=18 fission spectrum integral when MT=452 is absent (T11 case). All 1039 remaining diffs are ULP-class (8th sigfig in e15.8) — sub-7-sigfig precision floor. Phase 58e (next session): port `parse_endf_float` to Fortran-faithful precision, match cross-MT scat accumulation order, port resint flux/elastic from MF=3/MT=1. T01 regression-clean at NUMERIC_PASS 32812/32962 @1e-5. T22 / leapr / reconr untouched. | `worklog/T11_phase58_wimsr_complete.md` |
| 57    | 2026-05-01 | leapr `contin` Phase A — lat/sc/arat α/β rescaling | **T80 DIFFS → NUMERIC_PASS @1e-5.** `generate_sab` was missing Fortran `contin`'s `sc = therm/tev` rescaling for `lat=1` (leapr.f90:492-493 + 500/506). For T80 (lat=1, T=343 K, sc≈0.856) every α/β fed into `_terpt` and the SCT formula was off by ~14% in scale. Added `lat::Int=0, arat::Float64=1.0` kwargs to `generate_sab`; threaded through from `leapr_module`. Result: T80 tape24 75.7% (69231/91453, rtol=0) → **NUMERIC_PASS 99.95% (91406/91453 @1e-5)** — only 47 lines remain off, all ±1 in 7th sigfig. T22 (lat=0, sc=1 → unaffected) preserved BIT_IDENTICAL 4636/4636. All 4 leapr unit-test files pass. Phase B follow-up: SCT-replacement `naint` gating (leapr.f90:584-642 wraps `ssm(k,j)=ssct` in `iprt = mod(j-1, naint)+1 == 1`; Julia replaces all unconditionally). | `worklog/T80_leapr_contin_phase_a_lat_sc.md` |
| 56    | 2026-05-01 | T15 errorr writer NK = nmts−ix+1 (covcal Bug A) | **CLOSED.** `_write_errorr_tape` now synthesises one sub-section per (mt, mt2) with mt2 ≥ mt for **every** mt, not just `mt ∈ nc_derived_mts`. Empty cross-pairs route through `_write_mfcov_rows`'s `matrix===nothing` branch (2-line zero stub), mirroring Fortran covout iabort=1 path (errorr.f90:7350-7356, label 390). T15 MT=77 NK 1 → **3** (matches ref [77, 91, 102]); MT=18 NK 1 → 31; all 36 reaction MTs match ref NK exactly. T15 tape26 4275 → **5964 lines** (ref 5958, +6 from sub-section content drift). Residual per-MT line gaps (MT=1 −100, MT=2 +106, MT=18 +38, MT=102 −38) are sub-section *content* drift, not geometry — folded into HANDOFF P1 sub-item 2 (LTY=1/2/3 standards/ratio). Re-tuned `test_errorr_mf33_sparse.jl`: per-MT NK strict-match (all 36 MTs), per-MT slack 15 → 60, total cap `< 5000` → `5500 < total < 6500`. Regression-clean: T05/T16 covr-isolation 3/3 BIT-IDENTICAL, T04 tape23 NUMERIC_PASS 81/82 unchanged, T15 MT=77 C[20,20] = 0.02987998 exact (Phase 51 Bug B preserved). | `worklog/T15_covcal_bug_a_nk_writer.md` |
| 55    | 2026-04-30 | covr stub → full Fortran-faithful port (covr.f90 2249 lines, 18 subroutines) | **3/3 ISOLATION TAPES BIT-IDENTICAL.** Replaced empty-tape `covr_module` stub with the complete port: `covard` reader (errorr → group structure, MF3 cross sections, MF33/34/35/40 subsections with sparse-row LIST format), `corr` (cov→corr + ismall/izero flags), `truncg` (zero-xs lower-end truncation), `plotit` + `matshd` + `level` + `patlev` (3-frame plot-tape emission with connected-region shading), `smilab` + `matmes` + `elem` + `mtno` (Z + isotope + reaction labels with viewr `#EH`/`<x>` markup), `press` + `setfor` (boxer-format library output, plot-mode + library-mode dispatch). Card-1/2/2'/2a/3a/2b/3b/3c/4 parser handles `//` multi-card-per-line shorthand for consecutive defaults via new `ModuleCall.raw_lines` field. Bugs landed (each oracle-driven, each cited): title 80-char rpad, MT-name table row 13 (`]g<)` for n,γ), `(mt NNN)` fallback width, `plotit(x(ixmin:),...)` slice (rsdx/rsdy/xs_x/ys_y all ixmin-truncated), spurious NC sub-subsection skip in `_read_mf33_subsection` (Fortran covard never reads CONT.N1 — consumes only NI LIST records). T05 referenceTape34 (1 839 463 B), T16 referenceTape36 (270 922 B), T16 referenceTape37 (23 768 B): all byte-for-byte. T01 NUMERIC_PASS + T22 BIT-IDENTICAL regression-clean. Full-pipeline T05/T16 still diff because Julia's errorr (HANDOFF P1) produces sparse covariance — covr is correct end-to-end once errorr's NK-stub + LB=5 union-grid path lands. | `worklog/T05_T16_covr_full_port.md` |
| 54    | 2026-04-28 | gaspr stub → real MT=203/207 splice (T13/T45) | **STANDALONE BIT-IDENTICAL.** Wired the gaspr orchestration: read PENDF, accumulate gas production via existing `accumulate_gas` (multiplicity table for MT=11, 22-45, 103-117, 154-200), build new MT=203..207 sections on MT=1's grid (back-up-by-1 per channel, sigfig(7,0) values), splice into tape, update MF1/MT451 directory (NXC bump + MT-sorted MOD=1 entries). Standalone test (Fortran-oracle `after_heatr.pendf` → Julia gaspr): bit-identical to `after_gaspr.pendf` (25430/25430 lines). Pipeline integration via new `_collect_gaspr!` feeds MT=203..207 into `ctx.extra_mf3` so `final_assembly` places them at correct MT positions (T13 tape28 16878→18584; T45 tape40 now contains MT=203/204/205/207 with reasonable NC counts). Surfaced + fixed pre-existing pendf_io round-trip bugs: TPID had MAT=0 (spec is MAT=1), `write_pendf_tape` doubled SEND for MF1 + duplicate FEND between MFs, per-section sequence numbers didn't restart at 1, `read_pendf` retained SEND boundary records inside section data (now stripped on read). Deferred (open work): MOD=1 flag preservation in `write_full_pendf` for new sections, MF6/MT5 emission spectra, MT=51-91 LR=22..36 yields, multi-temperature loop. T22 BIT_IDENTICAL preserved; T01 NUMERIC_PASS @1e-5 unchanged. | (commits c77276f, 2a92c19) |
| 53    | 2026-04-27 | T50 ACER phases 5-6.2 (unionx step + LAW=5→LAW=14 + sigfig polish) | **NUMERIC_PASS @1e-5; 143/163 BIT-IDENTICAL.** Phase 5 ports the rest of unionx (MF6 anchor c-overwrite quirk, step-1.2 ratio enforcement off-by-one drop, gety1 dedup at aceout) — T50 NES 32→37, all ESZ values bit-identical. Phase 6 ports `acecpe`: Coulomb σ_C analytic formula (LIDP=1 identical-particle interference), per-μ trapezoidal (μ, pdf, cdf) build with adaptive midpoint subdivision when σ_C grows >2×, log-log interpolation of integrated yys onto ESZ for Coulomb-corrected total/elastic. Phase 6.1 applies sigfig(7) to μ/pdf and sigfig(9) to cdf at write time and stops recomputing absorption_xs (Fortran acecpe leaves it untouched). T50 tape34 status STRUCTURAL_FAIL → **NUMERIC_PASS @1e-5**, lines 59→**163** (matches ref), exact bit-identical lines 35→**143/163 (88%)**. T01 32812/32962 unchanged; T02 12519/13873 unchanged. Remaining 20 lines: ±1 in 7th sigfig from FP order in trapezoidal cumm and signow log-log interp — focused FP debug session. Neutron-path MF6 anchor union (originally added in Phase 4) restored after Phase 5 inadvertently dropped it. | `worklog/T50_acer_phase1_scaffolding.md` |
| 52    | 2026-04-24 | T50 ACER promotion, phases 1-4 (α+He-4) | **PARTIAL.** AcerParams card-1 parser fixed to 5 slots (every acer test was missing `ngend`); header BIT-IDENTICAL modulo date via MF1/MT451 reader + NSUB→letter + AWR override + Fortran-faithful hz/hd/hm alignment + f11.0 trailing-dot; JXS layout corrected (MTR/LQR/TYR/LSIG/SIG always populated; LDLW/DLW sentinel; END=length); MF6 incident-energy union added to ESZ grid (NES 29→32). T50 tape34 16/163 exact. T01/T02 regression-clean. Superseded by Phase 53. | `worklog/T50_acer_phase1_scaffolding.md` |
| 51    | 2026-04-22 | T15 covcal LB=5 σ·flx-weighted collapse (NJOY.jl-f8k Bug B) | **FIX LANDED.** MT=77 self-cov C[20,20] 0.04 → 0.02987998 (exact match to ref). `multigroup_covariance` extended with xs_row/xs_col/ugrid kwargs; orchestration refactored to collect blocks per (mt, mt2) pair and route LB=5/6 through union-grid σ·flx-weighted collapse. For iwt=2, flux = bin width (matches GENDF). T04 MF31 unchanged. | `worklog/T15_covcal_lb5_weighted.md` |
| 50    | 2026-04-21 | T15 covcal MT=77 diagnosis (NJOY.jl-f8k) | ROOT CAUSE PINNED — Bug A (writer NK count) + Bug B (midpoint sampling vs XS·flux-weighted union-grid expansion). Fortran trace in `/tmp/t15_fortran_diag/stdout.log`. **No code changes** — implementation deferred. | `worklog/T15_covcal_mt77_diagnosis.md` |
| 49    | 2026-04-21 | T15 errorr MF33 NC v2 sub-item 1 (double-NC cross-pairs) | Cov(2, 4) sub-section 3 → 65 lines (ref 69). T15 tape26 4178 → 4240. T04 baseline preserved. | `worklog/T15_T17_errorr_nc_expansion_v2.md` |
| 48    | 2026-04-20 | T15 errorr MF33 NC v1 (LTY=0 single-NC) | Cov(MT=2, *) and Cov(MT=4, *) cross-pairs to non-NC refs computed. T15 tape26 1859 → 4178. | `worklog/T15_T17_errorr_nc_expansion.md` |
| 47    | 2026-04-20 | errorr MF33 sparse per-row emission | T15 tape26 8205 → 1859 via Fortran-style ng2==0 row skip. | `worklog/T15_errorr_mf33_sparse.md` |
| 46    | 2026-04-19 | errorr GENDF MF3 readback | T15 tape26 MF3 0 → 36 MTs. | `worklog/T15_errorr_gendf_readback.md` |
| 45    | 2026-04-19 | groupr auto-expand reaction MTs from PENDF | T15 tape91 MF3 7 → 39 MTs (ref 40). | `worklog/T15_groupr_auto_expand.md` |
| 44    | 2026-04-18 | T15/T17 errorr MF33 size cap | Avoid 30M-line MF33 OOM by per-row sparse emission. | `worklog/T15_T17_errorr_size.md` |

Older phases (10-43) are recorded in commit messages and individual worklogs.

## Open Work (Authoritative — Bead-Independent)

Each entry below is fully self-contained: scope, file locations, Fortran
references, acceptance criteria, and notes. Retired bead IDs are listed
for cross-reference with older worklogs but are not required to pick up
the work.

### P1 — NC-block expansion v2: LTY=1/2/3 standards / ratio (sub-item 2 only)

- **Retired bead**: `NJOY.jl-km1` v2 sub-item 1 (double-NC-derived
  cross-pairs) landed in Phase 49, 2026-04-21 — see
  `worklog/T15_T17_errorr_nc_expansion_v2.md`. T15 Cov(2, 4) now
  emits a real 65-line LB-block (vs reference 69) instead of the
  3-line zero stub; tape26 grew 4178 → 4240. Sub-item 2 below remains.
- **Scope (sub-item 2 only)**: **LTY=1/2/3 standards / ratio blocks**:
  `_read_nc_subsection` currently returns `nothing` for these. Fortran
  handles them via `stand` around errorr.f90:7800 (NDS-standard
  covariance + ratio synthesis). T04 U-235 MT=18 has 4 LTY=3
  cross-material sub-subsections currently emitted as zero stubs
  (matches existing reference behaviour by coincidence — verify what
  Fortran actually emits before treating as correct).
- **Where**: `src/orchestration/modules/errorr.jl`
  (`_read_nc_subsection`); Fortran `njoy-reference/src/errorr.f90`
  `stand` (~7800).
- **Acceptance**: T04 reference test stays at NUMERIC_PASS 81/82 on
  tape23 (current baseline); when a test surfaces real LTY=3 demand,
  ported `stand` produces non-zero blocks matching reference.
- **Notes**: not blocking any current sweep failure. Pick this up
  when a test needs real derived-from-standards covariance values.

### P1 — Covcal content drift (NJOY.jl-f8k) — BUGS A + B LANDED; SUB-SECTION CONTENT DRIFT OPEN

- **Retired bead**: `NJOY.jl-f8k`.
- **Status (2026-05-01, post Phase 56)**: **Bug A fixed** — writer NK
  geometry now matches Fortran covout for every MT. T15 MT=77 NK
  1 → 3, MT=18 NK 1 → 31, all 36 reaction MTs match ref. T15 tape26
  4275 → 5964 lines (ref 5958). Worklog: `worklog/T15_covcal_bug_a_nk_writer.md`.
- **Status (2026-04-22, post Phase 51)**: **Bug B fixed** —
  `multigroup_covariance` extended with xs/flx weighting; orchestration
  refactored to route LB=5/6 through union-grid σ·flx-weighted collapse
  matching Fortran covcal. T15 MT=77 C[20, 20] = 0.02987998 (exact
  match to reference). Full walkthrough in
  `worklog/T15_covcal_lb5_weighted.md`.

- **Status (2026-05-18, post Phase 74)**: **T15 tape26 line-count
  gap closed from −68 to −6 (jul 5952 vs ref 5958).** `_expand_nc_blocks!`
  single-NC pass extended to handle `ref_j < mt` (store at `(ref_j, mt)`
  with column-side σ-ratio). MT=1/mt2=2 sub-section now matches ref on
  geometry (65 lines / 16 rows × 3 data lines each). Mirrors Fortran
  covout's sandwich-loop accumulation with the storage transpose
  derived from cov symmetry. Phase 72c canary preserved; T22
  BIT_IDENTICAL 4636/4636 preserved. 223 unit assertions PASS, zero
  regressions. Worklog: `worklog/phase74_mt1_nc_lower_ref.md`.

- **Status (2026-05-19, post Phase 75)**: **T15 tape26 line-count
  EXACT MATCH (jul 5958 = ref 5958, gap −6 → 0).** URR (LRU=2) rescon
  port. `read_mf32` now captures `MF32UnresolvedRange` (L-states with
  per-J (D, AJ, GNO, GG, GF, GX) bodies + relative cov matrix). New
  `_apply_rescon_urr_range!` in rescon.jl ports Fortran ggunr1
  (errorr.f90:6800-6905) + rpxunr (errorr.f90:4785-4985): forward-diff
  perturbation 1.01× per (J, param), ×1.015 E-step grid with
  elr/ehr/ehg breakpoints, cov_rel → cov_abs via b_i·b_j, sandwich
  through existing `_apply_sandwich!`. Sandwich contributions logged
  70 → 77 (one URR range × 7 RESCON_PAIRS). All MF=33 sub-section row
  layouts match ref byte-for-byte at the LIST-header level. Phase 72c
  MT=102 row-1 canary preserved (10/10 nonzero cols). Phase 73/74 NC
  expansion canaries preserved. T22 BIT_IDENTICAL 4636/4636 preserved.
  **306/306 unit assertions PASS** across 6 errorr/MF32 test files,
  zero regressions. Worklog: `worklog/phase75_urr_rescon.md`.

  The Phase 74 worklog framed the residual −6 lines as a
  `_resonance_group_window` boundary bug. That diagnosis was wrong —
  `_resonance_group_window` matches errorr.f90:3093-3108 exactly. The
  actual gap was URR rescon being unported (deferred since Phase 72).
  Per-row diff of (1,1)/(2,2)/(2,102)/(102,102) showed cols 12-15 of
  rows 13-15 zero in Julia where ref had nonzero values — exactly the
  cells the URR range [10 keV, 150 keV] / LANL groups 13..15 contributes
  to. Lesson: trust the Fortran source over the worklog when they
  disagree (CLAUDE.md Rule 2 fired).

  Residual sub-line-count work: MT=2/mt2=2 rows 10-12 sub-ULP FP
  precision (~5e-8 relative on MT=102 C[1,1]) plausibly traces to
  Julia's tanh-stretched perturbation grid vs Fortran's adaptive rpendf
  eskip mesh + a 6-digit wt6b constant in weight_functions.jl
  (3e-6 thermal effect). Deferred P3 — sub-line-count, doesn't affect
  tape26 line totals.

- **Status (2026-05-16, post Phase 73)**: **T15 tape26 line-count
  gap closed from +150 to −68 (jul 5890 vs ref 5958).** Two stacked
  bugs in the MF=33 path resolved: (1) cross-material sub-section
  pollution in U-238 MT=18 (5 MAT1-bearing sub-sections were
  contributing junk LB=6 data to the self-cov via the same
  `(18,18)` key); (2) `_expand_nc_blocks!` was missing the
  `(σ_i/σ_Y)²` per-cell σ-ratio in the absolute→relative cov
  conversion for NC-derived MTs (MT=2 elastic, MT=4 total inelastic).
  MT=2 Δ +142 → −2; MT=18 Δ +74 → 0; row-28 cell value error
  +7.42 → +2.50e-6 (3M× improvement). MT=1 Δ −64 and MT=102 Δ −2
  unchanged (separate cross-MT and rescon-row-14 issues, deferred).
  All Phase 72c canary tests preserved. Worklog:
  `worklog/phase73_nc_xmat_and_sigma_ratio.md`.

- **Status (2026-05-16, post Phase 72c)**: **T15 MT=102 row-1 canary
  GREEN at 1e-7.** Root cause of the Phase 72b ~35% systematic was an
  iwt mismatch: rescon assumed iwt=2 (flat weight), but T15's deck
  specifies `iwt=6` (thermal + 1/E + fission + fusion, errorr card 2
  line 34: `9237 3 6 1 1`). Fix replaced `_group_average_flat` with
  `_rpxgrp_average` (faithful Fortran rpxgrp port — errorr.f90:5227-
  5366) consuming `get_weight_function(iwt)` (egtwtf port —
  errorr.f90:10023-10110). C[1,1]: 3.580e-4 → **2.658408e-4** (ref
  2.658914e-4, diff −5.06e-8, ~0.019% relative). `@test_broken`
  flipped to `@test`. T22 BIT_IDENTICAL preserved. Worklog:
  `worklog/phase72c_rescon_iwt_fix.md`. Residual −5e-8 is plausibly
  Julia's tanh-stretched perturbation grid vs Fortran's adaptive
  rpendf eskip mesh + a 6-digit-vs-5-digit `wt6b` constant in
  `weight_functions.jl` (3e-6 fractional thermal effect) — flagged
  P3 for a future bit-faithful FP-grind.

- **Status (2026-05-08, post Phase 72b)**: **Sandwich landed — MT=102
  row 1 0/10 nonzero → 10/10 nonzero, signs match ref, magnitudes
  within factor of 2 (most within 35%); 1e-7 bound on `[1,1]` still
  pending FP-grind.** `apply_rescon!` now does the full pipeline:
  pointwise grid → 78 perturbation pairs (×2 for central diff) →
  `cross_section_rm` at each point → group-averaged sens →
  triangular/full sandwich over the 7 RESCON_PAIRS → divide by σ̄·σ̄
  for relative cov. Phase 71 canary suite: 40 PASS + 1 @test_broken.
  T22 BIT_IDENTICAL preserved. Phase 72c (2026-05-16) closed this
  out — see status block above. Worklog:
  `worklog/phase72b_rescon_sandwich.md`.

- **Status (2026-05-08, post Phase 72)**: **Foundation landed — MF=32
  reader + apply_rescon! skeleton.** `src/processing/mf32_reader.jl`
  parses U-238 JENDL (LRU=1/LRF=3/LCOMP=1) into typed Julia structs
  (10 resolved ranges, 317 resonances, 951 uncertain parameters);
  URR (LRU=2) range skipped cleanly. `src/processing/rescon.jl`
  defines `apply_rescon!` and the seven `RESCON_PAIRS` dispatch
  constant. Wired into `errorr_module` replacing the Phase-71
  `@warn`. 57 reader unit-test assertions PASS against the actual
  J33U238 ENDF tape (`test_mf32_reader.jl`). Phase 71 RED canary
  stays @test_broken (sensitivity builder not yet ported).
  Worklog: `worklog/phase72_mf32_reader.md`. Next session: port
  rpendf + rpxgrp + perturbation loop + sandwich fill to flip the
  canary GREEN.

- **Status (2026-05-07, post Phase 71)**: **Diagnosis revised — rescon
  is the dominant missing piece, NOT Bug-B siblings.** Matrix-level
  diff of T15 MT=102 self-cov reveals rows 1..14 zero in Julia,
  populated in ref (with negative entries proving sandwich-rule output).
  Fortran's `covout` (errorr.f90:7465) calls `rescon` (errorr.f90:8513)
  to add MF=32 RP-cov contributions to seven (mt, mt2) pairs:
  (1,1), (2,2), (2,18), (2,102), (18,18), (18,102), (102,102). Julia
  has zero MF=32 → MF=33 propagation. U-238 JENDL has an 8092-line
  MF=32 section. MF=32 detection scaffold + loud-warn landed in
  `errorr_module`; no functional rescon port yet. Worklog:
  `worklog/phase71_rescon_diagnosis.md`. RED canary in
  `test/validation/test_errorr_covcal_lb5.jl` (Phase 71 testset).

- **Scope (remaining — sub-section content drift)**: Per-MT line-count
  gaps within sub-sections after Bug A: MT=1 −100, MT=2 +106, MT=18
  +38, MT=102 −38, net +6 vs ref (Julia 5964 vs 5958). NK is correct;
  the gap is in *what* each sub-section emits. Three open paths,
  ordered by impact:
  - **rescon (MF=32 → MF=33 RP-cov sandwich)** — multi-session port
    (~1500-2500 LOC across 8+ Fortran subroutines: resprx, rpxlc0,
    rpxlc12, rpxlc2, rpxsamm, rpxunr, rpendf, grpav4, rescon).
    Closes the dominant share of MT=1 (−100), MT=18 (+38 mixed),
    MT=102 (−38), and contributes to MT=2 (+106 mixed). The rows-1..14
    drift in MT=102 self-cov is purely this.
  - **NC expansion v2 sub-item 2 — LTY=1/2/3 standards/ratio**
    (HANDOFF P1 sub-item 2 above) — Fortran's `stand` (errorr.f90:2939)
    synthesises NDS-standard covariance + ratio data; Julia returns
    `nothing`. Affects MT=2 NC-cross-pairs (the ~+106 overflow).
  - **Bug-B siblings (LB=1/LB=2 union-grid collapse)** — `expand_lb1`
    / `expand_lb2` use midpoint-sampling on the output grid, bypassing
    union-grid + flux/xs-weighted collapse. Smallest blast radius;
    likely contributes the small (e.g. MT=91 −6 lines) drifts.
    Originally hypothesised as the dominant cause; matrix dump in
    Phase 71 disproved this — the LB=5 weighted-collapse path
    produces *exact* values for output rows that fall in covered
    input bins (e.g. MT=102 rows 16..30 byte-identical).

#### Bug B — `expand_lb5_symmetric` midpoint sampling (large)

- **Symptom (canary)**: T15 MAT=9237 MF=33 MT=77 self-cov, output
  group 20 (LANL [1.353e6, 1.738e6]):
  - **Reference**: `2.987998e-2`
  - **Julia**: `4.000000e-2` (exactly the input bin's verbatim value)
  - The "round" Julia values across the entire matrix (4.0e-2,
    3.8e-2, 3.6e-2, etc.) are direct samples from MT=77's LB=5
    block with NO group averaging.
- **Where (bug)**: `src/processing/errorr.jl:70-82`
  `expand_lb5_symmetric` does midpoint-sampling: for each output group
  pair (g_i, g_j) it copies the single input bin value at the
  midpoint via `_find_bin(0.5*(egrid[i]+egrid[i+1]), ek)`. No
  group averaging at all. Same issue in `expand_lb5_full` (LB=6).
  `expand_lb1`/`expand_lb2` likely have analogous midpoint-only bugs
  but smaller blast radius (LB=1/2 are diagonal, simpler structure).
- **Where (caller)**: `src/orchestration/modules/errorr.jl:168` calls
  `expand_covariance_block(block, egn)` directly with the OUTPUT
  group grid `egn`, bypassing the union-grid + flux/XS-weighted
  collapse path. The `multigroup_covariance` path in
  `src/processing/errorr.jl:137` exists with a union-grid collapse
  but is unused by the orchestration.

- **Reference algorithm (Fortran covcal, errorr.f90:1770-2417)**:
  Build union grid `un = sort(unique(egn ∪ all_input_bin_edges))`.
  For each (jg, jh) of union groups (jh ≤ jg, symmetric):
  1. Find input bin `k = bin(un[jg])`, `l = bin(un[jh])`.
  2. Look up cov: `LB5_fvals[k, l]` (upper-tri index for LT=1).
  3. Accumulate `cov[jh] += LB5_fvals[k, l] * sig(jg) * sig1(jh)`
     (errorr.f90:2233 — XS·XS weighted at union-group level).
  4. After per-row accumulation, write
     `b[jg, jh] = cov[jh] * flx(jg) * flx(jh)` (errorr.f90:2341 —
     also flux-weighted).
  In `covout` (errorr.f90:7431-7438), collapse to output by summation
  over (jg→ig, jh→igp). Final relative cov (per the writer's reverse
  conversion):
  ```
  C_out_rel[ig, igp] = (Σ b[jg,jh]) / (sig_out[ig]·sig_out[igp]·flx_out[ig]·flx_out[igp])
  ```
  For piecewise-constant XS within an output group this collapses to
  pure flux-weighted averaging of the input fvals over (ig × igp).

- **Trace evidence (Phase 50)**: 990 DIAG lines in
  `/tmp/t15_fortran_diag/stdout.log` (regenerable — patch
  `errorr.f90:2230` per the diagnosis worklog). For LANL group 20
  there are **9 union groups** (jg ∈ 56..64) covering 1.35e6 → 1.738e6,
  with `sig` varying 7.46e-5 → 8.13e-3 across them. Bin 1
  (`[1e-5, 1.4e6]`) is sub-threshold (zero LB-block contribution) and
  drags the average from 4.0e-2 down to the reference 2.988e-2.

- **Implementation plan (next session)**:
  1. **Union grid**: in `errorr_module`, compute
     `ugrid = sort(unique([egn..., (b.energies for b in all_blocks)..., ...]))`.
  2. **Per-union-group XS** `sig_u[g]`: assign piecewise-constant from
     `group_xs[mt][ig_containing(g)]` as a first pass; later iteration
     can interpolate from the PENDF for higher fidelity.
  3. **Per-union-group flux** `flx_u[g]`: first pass use lethargy
     `log(ugrid[g+1]/ugrid[g])` or bin width `ugrid[g+1] - ugrid[g]`;
     real flux from groupr's GENDF tape (the `flx` array) in a
     follow-up.
  4. **Sandwich on union grid**:
     ```julia
     b[g, h] = LB5_fvals[k(g), k(h)] * sig_u[g] * sig_u[h] * flx_u[g] * flx_u[h]
     ```
  5. **Collapse to output**:
     ```julia
     for ig in 1:ngn, igp in 1:ngn
       num = sum(b[g,h] for g→ig, h→igp)
       den = (Σ sig_u[g]·flx_u[g] for g→ig) · (Σ sig_u[h]·flx_u[h] for h→igp)
       C_out[ig, igp] = num / den
     end
     ```
  6. **Bug A in same pass**: writer auto-synthesis of zero-stub
     cross-pairs for non-NC MTs.

- **TDD hook (RED → GREEN)**:
  Extend `test/validation/test_errorr_mf33_sparse.jl` with element-wise
  assertions on the MT=77 canary cell:
  ```julia
  @test abs(jul_mt77[20, 20] - 2.987998e-2) < 1e-6
  @test jul_pair_lines[(77, 91)] >= 2     # zero-stub presence
  @test jul_pair_lines[(77, 102)] >= 2
  ```
  Plus relax/tighten the global gap assertion as expected line counts
  shift (target T15 tape26 total > 5500 lines).

- **Acceptance**:
  - MT=77 self-cov C[20,20] within 1e-6 of reference 2.987998e-2.
  - MT=77 NK = 3 sub-sections (self + 2 zero stubs).
  - MT=77 total line count drops from 30 → ~20.
  - T15 tape26 total > 5500 lines (vs current 4240, ref 5958).
  - T04 reference test does NOT regress below NUMERIC_PASS 81/82 on
    tape23 (the MF31 path uses the same expand pipeline; LB=2 nubar
    blocks are smaller and may already be close, but verify).
  - `test_errorr_mf33_sparse.jl` and `test_errorr_gendf_readback.jl`
    must still pass; the per-MT line-count bands in the sparse test
    will need re-tuning as line counts fall toward reference.

- **Estimated cost**: ~half day for first iteration (piecewise-constant
  XS + lethargy flux); another ~half day for fidelity work (real
  per-bin XS + GENDF flux read).

- **Risk / known scope-creep**:
  - LB=6 (asymmetric) needs the same union-grid treatment.
  - LB=1/2 (diagonal) probably needs the same averaging — currently
    `expand_lb1` and `expand_lb2` also bin-sample at midpoint, so they
    likely have the same class of bug for fine-output / coarse-input
    cases. Worth checking once LB=5 lands.
  - Bug A change in the writer affects every MT, not just MT=77 —
    will likely shift the T15 tape26 line count by hundreds of lines
    via zero-stub additions to each non-NC MT's sub-section list.

- **Where to look in Fortran (subroutines to read first)**:
  - `covcal` (errorr.f90:1770-2417) — main loop
  - `covcal` LB=5 branch (errorr.f90:2208-2235) — the per-pair lookup
  - `covcal` per-row write (errorr.f90:2336-2353) — flux multiplication
  - `covout` (errorr.f90:7018-7787) — collapse from union grid to
    output groups; pair iteration (lines ~7393-7440)
  - `resprp` (errorr.f90:8009-8511) — NOT relevant for MT=77 (that's
    resonance-parameter covariance, used only when MF=32 present).

### P2 — groupr: skip empty-MT, derive MT=251/252

- **Retired beads**: `NJOY.jl-5oi` (empty-MT skip), `NJOY.jl-cdy`
  (MT=251/252 derivation).
- **Scope**:
  1. **5oi**: U-238 MT=37 emits all-zero groups (threshold 17.82 MeV
     above LANL-30 top 17.0 MeV). Fortran skips such MTs; Julia
     emits them, producing an extra 3-line block per MT. Surfaced
     from T15 tape91 `MF=3: 3 → 39 MTs (ref 40)` mismatch
     (`worklog/T15_groupr_auto_expand.md`).
  2. **cdy**: MT=251 (mubar, average cosine of scattering) and
     MT=252 (xi, average logarithmic energy decrement) are derived
     by Fortran groupr from MF=4/6 angular distributions. Julia
     emits them only if already present on the PENDF. Needs a
     dedicated derivation path modelled on the Fortran `gendir` /
     `convr` data flow.
- **Where**: `src/orchestration/modules/groupr.jl` (`_groupr_expand_auto`,
  plus a new MT=251/252 builder); Fortran `njoy-reference/src/groupr.f90`
  around `nextr` (1087-1123) and MT=251/252 derivation (search for
  `mt251`, `mt252`).
- **Acceptance**: T15 tape91 MF=3 MT count exactly 40 (matches
  reference), no all-zero groups in Julia output, MT=251 and MT=252
  values within 1e-7 of reference.
- **Notes**: `worklog/T15_groupr_auto_expand.md` has the full
  context. The Phase 45 auto-expand work is the upstream entry
  point. Scope-creep risk: MT=251/252 derivation pulls in the MF=6
  reader which is partially implemented.

### P2 — Cross-ign flux-weighted collapse (colaps)

- **Retired bead**: none (never filed).
- **Scope**: `errorr_module` currently handles only same-ign (groupr
  ign == errorr ign) GENDF readback. For cross-ign (coarser output
  grid than input), Fortran `colaps` (errorr.f90:9255-9283) does a
  flux-weighted sum of σ over fine groups falling inside each coarse
  group. No Julia port yet.
- **Where**: Build `_errorr_colaps` helper in
  `src/orchestration/modules/errorr.jl`; read flux from input
  GENDF's MF=3/MT=0 or compute from user_egn.
- **Acceptance**: A test case where input GENDF uses one group
  structure and errorr's `ign` is a coarser one produces
  flux-weighted σ matching reference.
- **Notes**: NOT blocking T15/T16/T17 (same ign). File a trigger
  when a sweep failure points here.

### P2 — T49 MLBW ±1 at E=110487.7 eV (MT=1, MT=2)

- **Retired bead**: none.
- **Scope**: T49 RECONR passes 44/46 MTs PERFECT; MT=1 and MT=2 have
  ±1 in last sigfig at a single energy E=110487.7 eV. Looks like the
  same "irreducible FP precision" class that Phases 12-14 proved is
  always a real bug (see T34 Frobenius-Schur fix).
- **Where**: `src/resonances/mlbw.jl` `cross_section_mlbw`; Fortran
  `njoy-reference/src/reconr.f90` `csmlbw` (~2988-3023 area).
- **Acceptance**: T49 reports 46/46 PERFECT, no ±1 at
  E=110487.7 eV.
- **Notes**: Use gdb to trace Fortran csmlbw intermediate values at
  that energy (same recipe used to diagnose T34). Patch Julia to
  match the exact intermediate-rounding order.

### P3 — T04 tape25 12-line gap (sub-precision FP, same class as T34)

- **Retired bead**: none.
- **Status (2026-05-23, Phase 76c diagnosis)**: HANDOFF's previous
  "MF31 LB=2 union-grid collapse" diagnosis is **wrong** — verified
  by tracing `_errorr_second_call`. The 12 differing lines (107/119
  match at no-tolerance; 119/119 pass at 1e-7) are MF=3/MT=452 nubar
  values, not MF=33 covariance.

  The ENDF U-235 t511 MF=1/MT=452 evaluation has ν(E) = constant
  2.4367 across [1e-5, 25000] eV (3 breakpoints with identical
  value, INT=2). Yet Fortran's GENDF tape24 emits ν=2.436701 for
  LANL group 10 (1235-3350 eV) while ν=2.436700 for the other groups
  — a 1-unit-in-7th-digit shift from sigfig(7) rounding at the FP
  boundary. Julia's groupr emits 2.436700 for all groups (rounds
  the same way as the integer 2.4367). The user-grid nubar collapse
  then propagates this 1e-6 fractional difference to user groups 3-7
  of tape25.
- **Where**: `src/orchestration/modules/groupr.jl`
  `_groupr_nubar_records` — uses σ_f·flux-weighted integration with
  1% sub-stepping. Already matches Fortran's groupr panel + getwtf
  to first order. The remaining ULP shift comes from FP accumulation
  order in the inner loop — same class as T34 capture `csrmat`
  Frobenius-Schur accumulation (Phase 14 gdb-confirmed irreducible).
- **Acceptance**: T04 tape25 119/119 at no-tolerance would require
  matching Fortran's exact `panel` integration order in `_groupr_
  nubar_records`. Needs a gdb diagnostic dump of Fortran's per-bracket
  integrand values vs Julia's at LANL group 10.
- **Notes**: Tape25 already passes at 1e-7 (the project's first-round
  acceptance bar per CLAUDE.md Rule 1). Demoted from P2 to P3.
  Orthogonal to NC-block / LB=5 / rescon work.

### P2 — broadr U-238 JENDL performance (TIMEOUT)

- **Retired bead**: `NJOY.jl-326`.
- **Scope**: T15 and T17 full-pipeline reference tests TIMEOUT in
  broadr. Latest measurement (Phase 48 post-fix run): broadr
  completes in 1039 s, exceeding the 300 s soft cap and the 1800 s
  hard timeout (the framework actually gives ~1168 s wallclock, but
  the reference-test harness enforces its own limits).
- **Where**: `src/orchestration/modules/broadr.jl`; Fortran
  `njoy-reference/src/broadr.f90`. The slow case is sigma1 Doppler
  broadening over U-238's dense URR tabulation.
- **Acceptance**: T15 broadr under 300 s, tapes 22/23 produced.
- **Notes**: Profile first with `@profile` / flamegraph; don't
  refactor blind. The Fortran version of this exact case runs in
  ~30 s, so a 30× speedup is plausible from hot-loop
  devectorization or preallocation wins.

### P3 — T65 errorr performance (U-235 MF34)

- **Retired bead**: none (noted under `NJOY.jl-326` but orthogonal).
- **Scope**: T65 errorr MF=34 (angular distribution cov) is slow
  enough to be on the sweep watchlist. Not yet profiled.
- **Where**: `src/orchestration/modules/errorr.jl` MF=34 path
  (grep `mfcov.*34`).
- **Acceptance**: T65 under sweep default timeout with real output.
- **Notes**: Worth re-measuring first — may have been subsumed by
  Phase 44's 30M-line MF=33 fix.

### P1 — ACER promotion (T50 α+He-4 and 8 sibling tests)

- **Retired bead**: none.
- **Status (Phase 53, 2026-04-27)**: parser + iopt=7 passthrough + header
  bit-identical + JXS layout + total-XS-from-partials (P1-3) + full
  `unionx` port (MF6 anchor merge with c-overwrite quirk + step-1.2
  ratio enforcement + gety1 dedup) (P5) + `acecpe` port (LAW=5 → LAW=14
  Coulomb+nuclear, identical-particle interference, adaptive subdivision
  at Mott singularity, log-log interpolated Coulomb-corrected
  total/elastic) (P6). T50 `tape34` status STRUCTURAL_FAIL →
  **NUMERIC_PASS @1e-5** with line count exactly matching ref (163);
  103/163 lines pass at 1e-5. T01/T02 unchanged. See
  `worklog/T50_acer_phase1_scaffolding.md` for full detail.
- **Scope**: promote `acer_module` from MF3-only stub to real ACE
  generator passing T50 at 1e-9. Unlocks 8 sibling tests (T14, T48,
  T51-T54, T62, T71 — all RAN_OK today without oracle comparison).
- **Next work items**:
  1. **Tighten T50 from 1e-5 to 1e-9** (60 sub-1e-5 lines remain): FP
     rounding in per-μ Coulomb evaluation, trapezoidal cumm
     accumulation order, log-log interp on yys. Method: `write(*,...)`
     diagnostics in Fortran `acecpe` to capture intermediate sigc /
     signi / cumm / signow values; match Julia's order step-by-step.
  2. **Sibling tests T51-T54, T62, T71**: verify NSUB→letter mappings
     (proton+H-2, deuteron+H-2, etc.), exercise additional reaction
     types (non-elastic MTs in MF6, MT=600+ for proton emission).
     Each test surfaces feature gaps to port piecewise.
  3. **iopt=7 aplots port**: real viewr plot-tape generator, not empty
     stub. Lower priority — unblocks viewr testing only.
- **Where**: `src/orchestration/modules/acer.jl`,
  `src/formats/ace_builder.jl`, `src/formats/ace_charged.jl` (NEW),
  `src/endf/readers.jl`.
- **Acceptance**: T50 `tape34` 163/163 exact at 1e-9 via execute.py;
  no regression on T01 ACE or any currently-passing test.
- **Estimated cost**: ~½ session for the FP-grind to 1e-9 (Phase 7);
  ~1-2 sessions per sibling test as new feature gaps surface.

### P3 — ~55 DIFFS cases — per-tape bisection (Grind Method)

- **Retired bead**: none.
- **Scope**: After Phases 10-48, the post-sweep state is roughly
  ~64 DIFFS, 17 MISSING_TAPE, 0 CRASH, and 2 NUMERIC_PASS + 1
  BIT_IDENTICAL on the 84-test suite. Most DIFFS are small last-digit
  or structural issues per the Grind Method.
- **Where**: Each test's section in `reports/REFERENCE_SWEEP.md`.
  Use the Grind Method (HANDOFF "The Grind Method" section + Rule 0
  oracle-driven TDD).
- **Acceptance**: each DIFFS → NUMERIC_PASS or better.
- **Notes**: Re-run a fresh full sweep first to get the post-Phase-48
  baseline. Expected ~90-100 min.

### P3 — MF7/MT2 Bragg reader

- **Retired bead**: none.
- **Scope**: thermr currently uses a hardcoded `BRAGG_LATTICE_PARAMS`
  table for the small set of materials it supports. For general
  ENDF-6 thermal evaluations (T25, T67-T74) the per-material Bragg
  lattice parameters should be read directly from the evaluation's
  MF=7 MT=2 section.
- **Where**: `src/processing/thermr.jl` (or wherever the hardcoded
  constants live — grep `BRAGG_LATTICE_PARAMS`); Fortran
  `njoy-reference/src/thermr.f90` `bragg` subroutine.
- **Acceptance**: BRAGG_LATTICE_PARAMS constant removed; T25, T67,
  T68, T69, T70, T74 all use real MF7/MT2 reads.
- **Notes**: landed as a stub in Phases 12-14 gating; unblocks a
  cluster of 6-7 tests when done.

### P2 — leapr `contin` Phase B: T80 BIT_IDENTICAL @1e-9 (47 lines remain)

- **Retired bead**: none.
- **Status (2026-05-01, post Phase 57)**: Phase A landed (lat/sc/arat
  α/β rescaling). T80 went DIFFS 75.7% → **NUMERIC_PASS 99.95% @1e-5**
  (91406/91453 lines). T22 BIT_IDENTICAL preserved. Worklog:
  `worklog/T80_leapr_contin_phase_a_lat_sc.md`. Three follow-up items
  to close the residual 47 lines:
  1. **SCT-replacement `naint` gating — NON-ISSUE (verified 2026-05-28).**
     The earlier claim here was WRONG. `naint` is hardcoded to 1 in Fortran
     (leapr.f90:292-294: `naint=0; if(naint.eq.0) naint=1`) and is **never
     read from input**. So `iprt = mod(j-1,naint)+1 = mod(j-1,1)+1 = 1` for
     every α — the gate at leapr.f90:585-587 is always open. Julia's
     unconditional SCT replacement (leapr.jl:251-257) already matches Fortran
     exactly; there is no naint kwarg to thread. (`naint` is NOT a card-3
     field — that was HANDOFF drift.) Bead `NJOY_jl-a6q` closed as not-a-bug.
  2. **Phonon-loop FP-order alignment.** `_convol!` and the in-place
     `xa[j] += log(al*f0/n)` accumulation (leapr.f90:533) may differ
     in IEEE-754 order. Diagnose with `write(*,...)` traces in Fortran
     `contin` at a single (j, k) point — recipe in CLAUDE.md.
  3. **Moment-check loop structure.** Fortran combines diagnostics
     printout and SCT replacement in one outer-α loop (lines 581-642).
     Once `naint` gating lands, the Julia structure should mirror it
     exactly so additional accumulation-order matching becomes trivial.
- **Where**: `src/processing/leapr.jl` `generate_sab`;
  `src/orchestration/modules/leapr.jl` callsite; `LeaprParams.naint`.
- **Acceptance**: T80 tape24 BIT_IDENTICAL 91453/91453 @ rtol=1e-9.
  T22 stays BIT_IDENTICAL 4636/4636.

### P3 — Real plotr / purr output (covr/leapr/gaspr now wired)

- **Retired bead**: none.
- **Status**:
  - **covr**: Phase 55 landed — full Fortran-faithful port,
    T05/T16 isolation tapes 3/3 BIT_IDENTICAL.
  - **gaspr**: Phase 54 landed — MT=203/207 splice, T13 standalone
    BIT_IDENTICAL.
  - **leapr**: Phase 57 landed — T22 BIT_IDENTICAL, T80 NUMERIC_PASS
    @1e-5 (Phase B above closes the remaining gap).
  - **purr**: Phase 38 wired MT=152/153 structurally; T38 STRUCTURAL_FAIL
    (3642 lines vs ref 3480). Real value-fidelity work remains.
  - **plotr**: Phase 11 dispatch stub only; pure plot-tape, no
    downstream chain — lowest priority.
- **Scope (remaining)**: purr value-fidelity (T28, T34-partial, T35-42,
  T63, T72) + plotr real output (T06).
- **Acceptance**: per-module real output matches reference within
  1e-5 for at least one test each.
- **Notes**: leapr P2 above is the immediate priority; purr/plotr can
  follow.

### Known completed (retired bead IDs, do not re-do)

| Retired ID | Scope | Resolved in |
|------------|-------|-------------|
| `NJOY.jl-327` | heatr plot-tape stub | Phase 14 |
| `NJOY.jl-km1` | NC-block expansion v1 | Phase 48 (2026-04-20) |
| `NJOY.jl-km1` | NC-block expansion v2 sub-item 1 (double-NC-derived cross-pairs) | Phase 49 (2026-04-21) |

---

## Test Infrastructure

### Fortran-faithful reference test framework (Phase 8 — 2026-04-13)

**Primary test harness.** For each Fortran reference test NN, runs the
**real input deck** `njoy-reference/tests/NN/input` through `run_njoy(…)` and
compares every produced `tape{U}` against `referenceTape{U}` using line-
equivalence semantics identical to `njoy-reference/tests/execute.py`.

- `test/validation/reference_test.jl` — `run_reference_test(N)` generic runner.
- `test/validation/sweep_reference_tests.jl` — loops all 84, writes `reports/REFERENCE_SWEEP.md`.
- `test/validation/module_coverage.jl` — writes `MODULE_COVERAGE.md` (which Fortran tests cover each Julia port).
- `test/runtests.jl` — `@testset "Reference Tests (Fortran-faithful)"` at end, invoked by `Pkg.test()`.

Verbose output + 10s heartbeat (`[Txx] … still in <module> (Xs) last: …`) means
hangs are visible immediately. Single process, sequential — cache safety.
Hard timeout at 1800s per test (throws `InterruptException` into the task; may
overshoot if task has no yield points, but flags `:TIMEOUT` cleanly).

See `worklog/T08_reference_test_framework.md` for full design notes.

```bash
julia --project=. test/validation/reference_test.jl 3        # single test
julia --project=. test/validation/sweep_reference_tests.jl   # all 84
julia --project=. -e 'using Pkg; Pkg.test()'                 # full suite
```

### Sweep baseline (2026-04-13, 2h 13m, full 84 tests)

First clean sweep after Phase 8 framework landed. Superseded by the
Phase 10 sweep below — kept here for Δ tracking.

| Status           | Count |
|------------------|-------|
| `BIT_IDENTICAL`  | 1     |
| `NUMERIC_PASS`   | 1     |
| `DIFFS`          | 12    |
| `MISSING_TAPE`   | 14    |
| `NO_REFERENCE`   | 1     |
| `CRASH`          | 55    |

### Post-Phase-10 sweep (2026-04-14, 2h 43m, full 84 tests) — **CURRENT**

`reports/REFERENCE_SWEEP.md`. Numbers are from the sweep run immediately
after Phase 10's five dispatch/fix items landed:

| Status           | Pre-P10 | Post-P10 | Δ      |
|------------------|---------|----------|--------|
| `BIT_IDENTICAL`  | 1       | 1        | =      |
| `NUMERIC_PASS`   | 1       | 1        | =      |
| `DIFFS`          | 12      | **48**   | **+36**|
| `MISSING_TAPE`   | 14      | 17       | +3     |
| `NO_REFERENCE`   | 1       | 1        | =      |
| `CRASH`          | **55**  | **16**   | **−39**|

**39 tests moved CRASH → DIFFS.** Pipeline runs end-to-end on 68 of 84.

**Crash taxonomy (16 remaining)** — `/reports/REFERENCE_SWEEP.md` has
per-test detail:

| Count | Exception                                                       | Root cause                                                                    |
|-------|-----------------------------------------------------------------|-------------------------------------------------------------------------------|
| 9     | `SystemError: opening file ".../tapeNN"`                        | Undispatched upstream modules: `covr` (T05/T06/T16), plus separate tape-unit plumbing issues in T12/T18/T24/T27/T34/T47/T65. |
| 2     | `BoundsError: attempt to access 0-element Vector{Float64}`      | **NEW failure mode uncovered by the INT=0 fix** — T15/T17 JENDL U-238 now gets past the reader and trips a downstream empty-vector indexing bug. |
| 1     | `MF7/MT4 not found for MAT=101`                                 | T09: leapr stub writes empty file; thermr tries to read real MF7. Upgrade path: real leapr output. |
| 1     | `tape unit 0 is invalid`                                        | T20 moder card uses `0` as a "no tape" sentinel; moder_module needs to skip it. |
| 1     | `broadr: MT=1 (total) not found on input PENDF`                 | T60 Fe-nat IRDFF-II (MF10-only dosimetry).                                    |
| 1     | `InexactError: Int64(NaN)`                                      | T43 broadr at T=0 — NaN→Int conversion.                                       |
| 1     | leapr chain downstream effects                                  | — (covered under T09)                                                         |

**Phase-10-fixed crashes** (CRASH → DIFFS or better):
- heatr `h` (6 tests): T08, T13, T21, T26, T49, T79.
- acer dispatch (7 tests): T14, T50-T54, T62, T71 moved CRASH → DIFFS/MISSING_TAPE.
- purr dispatch (13 tests): T28, T34(partial), T35-T42, T63, T72.
- leapr dispatch (3 tests): T22, T23, T80.
- Bragg gating + graceful degrade (7 tests): T25, T32, T67, T68, T69, T70, T74.
- INT=0 (2 tests): T15/T17 progress past reader; new downstream bug revealed.

**Priority follow-ups**: historical post-Phase-14 list superseded by
the "Open Work" section at the top of this file. Items 1-8 of the old
list all landed in Phases 11-14 (covr, plotr, errorr 999-mode stub,
broadr T=0, moder iopt=1, heatr nplot, broadr MT=1 fallback, thermr
free-gas fallback). The remaining items are folded into the current
Open Work entries: the fresh sweep re-run is pending; TIMEOUTs tie
into "broadr U-238 JENDL performance" and "T65 errorr performance";
"Real plotr / covr / leapr / purr output" is its own P3 entry; the
MF7/MT2 Bragg reader is unfiled (add it to Open Work if a test needs
it); the DIFFS bisection is "~55 DIFFS cases" P3.

See `worklog/T10_phase10_batch_dispatch.md`, `T11_covr_dispatch.md`,
`T12_small_batch.md`, `T13_moder_gaspr.md`, `T14_drain_crash.md` for
implementation detail on the completed items.

### Legacy oracle system (superseded by above for cross-module tests, still useful for reconr-only grind)
Each test's Fortran reference output is cached in `test/validation/oracle_cache/testNN/`. The `diagnose_harness.jl` script generates these by running the Fortran NJOY binary with truncated input decks. Each oracle directory contains:
- `after_reconr.pendf` — ASCII PENDF after reconr (the main comparison target)
- `after_broadr.pendf` — ASCII PENDF after broadr (if the test chain includes broadr)
- `run_reconr/input` — truncated input deck
- `run_reconr/tape20` — ASCII ENDF input tape (preferred for Julia)
- `run_reconr/tape21` or `tape22` — may be binary (from moder) — avoid with Julia

### Quick comparison script
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e '
using NJOY
r = reconr("test/validation/oracle_cache/testNN/run_reconr/tape20"; mat=MAT, err=ERR)
write_pendf_file("/tmp/testNN.pendf", r; mat=MAT, err=ERR)
# ... then parse and compare MF3 columns 1-66 ...
'
```

### Tape selection
- **Always prefer tape20** (ASCII). Some tests use `moder` to convert ASCII→binary before reconr.
- If only binary tapes exist (tape21/tape22), you need Julia's `moder()` function first, or find the original ASCII source in `njoy-reference/tests/resources/`.
- The reconr input deck tells you which tape reconr reads (the first number after `reconr`). Negative = binary.

### Current oracle coverage
- **30 tests have oracle caches**: T01-04, T07-13, T15-21, T24-27, T30, T34, T45-47, T49, T55-58, T60, T63-65, T84
- **16 tests need oracles**: T24\*, T28, T29, T31, T32, T35-44, T63\* (\*have cache dirs but no run_reconr)
- **11 tests skip** (no RECONR in chain): T05, T14, T50-54, T59, T61, T62

---

## The Grind Method

This is how we achieve bit-identical output. It works. Don't skip steps.

### Per-test workflow

1. **Pick a test** (ordered by complexity: LRU=0 → SLBW → MLBW → Reich-Moore → SAMMY)
2. **Generate the Fortran oracle**: `julia --project=. test/validation/diagnose_harness.jl <test_number>`. This runs the Fortran binary and caches per-module PENDF output in `test/validation/oracle_cache/test<NN>/`
3. **Run Julia RECONR** and write PENDF:
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e '
using NJOY
result = reconr("njoy-reference/tests/resources/<file>"; mat=<MAT>, err=<ERR>)
write_pendf_file("/tmp/julia_test.pendf", result; mat=<MAT>, err=<ERR>)
'
```
4. **Compare MF3 data** (columns 1-66, ignoring sequence numbers in columns 76-80):
```bash
julia --project=. -e '
using NJOY
function parse_all_mf3(fn)
    result = Dict{Int, Vector{String}}()
    lines = readlines(fn); idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = NJOY._parse_int(p[71:72]); mt = NJOY._parse_int(p[73:75])
        mat = NJOY._parse_int(p[67:70])
        if mf == 3 && mt > 0 && mat > 0
            idx += 3; data = String[]
            while idx <= length(lines)
                length(lines[idx]) < 75 && (idx += 1; continue)
                pd = rpad(lines[idx], 80)
                mfc = NJOY._parse_int(pd[71:72]); mtc = NJOY._parse_int(pd[73:75])
                (mfc != 3 || mtc != mt) && break
                push!(data, pd[1:66]); idx += 1
            end; result[mt] = data
        else; idx += 1; end
    end; return result
end
function compare(name, jf, ff)
    j = parse_all_mf3(jf); f = parse_all_mf3(ff)
    mts = sort(collect(union(keys(j), keys(f)))); np = 0
    for mt in mts
        jd = get(j, mt, String[]); fd = get(f, mt, String[])
        if jd == fd; np += 1
        elseif isempty(jd); println("$name MT=$(lpad(mt,3)): MISSING")
        else
            d = findfirst(i -> i > length(jd) || i > length(fd) || jd[i] != fd[i],
                          1:max(length(jd),length(fd)))
            println("$name MT=$(lpad(mt,3)): DIFF $(length(jd))v$(length(fd)) @$d")
            d !== nothing && d <= min(length(jd),length(fd)) &&
                println("  J: $(jd[d])\n  F: $(fd[d])")
        end
    end
    println("$name: $np/$(length(mts)) PERFECT")
end
compare("T07", "/tmp/julia_test.pendf",
        "test/validation/oracle_cache/test07/after_reconr.pendf")
'
```
5. **Find the first byte that differs**. Trace to root cause. Read the Fortran. Fix.
6. **Repeat** until all MTs show PERFECT.

### Debugging a diff

When you see a diff, classify it:
- **Grid diff** (different energies): issue is in `lunion_grid`, `adaptive_reconstruct`, or the initial grid construction
- **XS diff at same energy** (±1 in last digit): issue is in sigma evaluation, sigfig rounding, float accumulation order, or threshold interpolation
- **Large XS diff**: missing feature (e.g., unresolved resonance evaluation was completely missing)
- **Format diff** (same value, different string): issue is in `format_endf_float` (the `a11` equivalent)
- **Line count diff** (same data, shifted lines): extra or missing grid points from shading or threshold handling

For ±1 diffs, common root causes (all encountered and fixed in this project):
- Missing intermediate `round_sigfig(x, 8)` call before addition (SLBW elastic)
- Wrong float accumulation order (non-associativity of IEEE 754)
- Testing total instead of partials in adaptive convergence
- Hardcoded linear interpolation instead of using MF3's actual law (LinLog, LogLog)
- Missing boundary-crossing forced convergence in adaptive reconstruction
- MF3 backgrounds added rounded instead of unrounded to resonance values
- Missing `round_sigfig(x, 7)` on primary channel output values
- URR table missing rdfil2 boundary nodes (sigfig(EL, 7, ±1) shading)
- MF=12 histogram sections not processed (creates shaded pairs at shared breakpoints)

### Effective subagent strategy (proven in Phase 6)

For complex debugging, launch 3+ research subagents **in parallel**:
1. **Fortran source reader** — have an Explore agent deeply read the relevant Fortran subroutine
2. **ENDF data inspector** — have an agent examine the raw ENDF file (grep for specific energies, list breakpoints, check interpolation laws, find MF=12/13/23 sections)
3. **Grid comparator** — have ONE agent run Julia and compare grids/MF3 output

This 3+1 pattern (3 read-only researchers + 1 Julia runner) was how the MF=12 breakthrough was found. The data inspector discovered that MF=12 MT=102 had histogram interpolation at the exact mystery energies — something that manual Fortran code tracing missed for hours.

---

## Reference Oracle Pin

The Fortran oracle (`njoy-reference/`, a per-machine git clone of NJOY2016) is **pinned** so "bit-identical" validation is measured against a frozen commit instead of upstream `develop` (which rebases). The pin lives in `REFERENCE_PIN` (repo root):

- **Pinned SHA:** `2c64dfb3d7cd57bea9fa8e36b995eee4a9d88d58` (upstream `develop` tip at pin time, commit dated 2026-03-18, "updating release notes"). Pinned 2026-06-01.
- **Align your local oracle:** `bash scripts/setup_reference.sh` (idempotent — clones if absent, else fetch + checkout the pin, then verifies HEAD).
- **Preflight:** `test/validation/reference_pin.jl::check_reference_pin()` runs at the top of both `reference_test.jl` and `sweep_reference_tests.jl`. It only `@warn`s on drift/missing oracle (never aborts) and points at the setup script.
- **ac5adf5 is a ghost:** the prior recorded baseline short-SHA `ac5adf5` is **not reachable** in current upstream history (rebased away) and is unreproducible — it is superseded by the pin above.
- **Drift from the old laptop oracle** (`3c9873d`, Dec 2025) to this pin: only 5 existing reference tapes changed, each a single MF1/MT451 record-3 field-2 (EMAX) header going `0.000000+0` → a real value (tests 09, 22, 23, 33, 80); plus two brand-new tests **86, 87** were added.

---

## Current State (2026-06-15 — T70 tape60 MF3/MT221+222 emax point restored (Phase 87); 9 BIT_IDENTICAL, 4 NUMERIC_PASS)

**2026-06-15 (Phase 87, bead `NJOY_jl-js4` closed):** T70 tape60's **MF3/MT221 + MF3/MT222 NP off-by-one is fixed** (1484 → 1485). `_thermr_sab!` (`src/orchestration/modules/thermr.jl`) was appending the thermal-cutoff grid as `[emax, sigfig(emax,7,+1), 2e7]`, omitting `sigfig(emax,7,-1)` (=2.276999 for emax=2.277) — the shaded-DOWN held point that the Fortran emits as the coherent grid's last point (sigcoh thermr.f90:1194), with the value also held at emax (clamp, 393), stepping to zero at `up*emax` (394) and `etop=2e7` (tpend 3211-3213); MT221 (ix=3) and MT222 (ix=4) share the `ex` energy grid. **MF3/MT222 is now byte-identical to `referenceTape60` (0 diffs, was 4 incl. structural); MF3/MT221's structural + `>1e-5` diffs are gone** (14 residual diffs all ≤4.8e-7, the calcem/Sigma1 FP-grind class under c3q). The fix = a new `_append_emax_sentinels!` helper (cites the 4 Fortran lines) + the one line `push!(grid, round_sigfig(emax,7,-1))`. **Section line count stays 498 and the Phase-86 directory NC stays 498** (2968/2970 vals both `ceil` to 495 data lines; `3+div(np+1,3)`=498 both) — confirmed by zero MF1/MT451 diffs. New RED→GREEN unit test `test_thermr_emax_sentinels.jl` (6/6); Phase-86 `test_thermr_mf3_nc.jl` 10/10. **T01 NUMERIC_PASS 32812/32962 + T09 BIT_IDENTICAL 1830/1830 preserved (independently re-run, not just subagent-reported).** **Rule-2 catch:** the reference test's regex `line_equivalence` mis-reported the grind as "71767 numeric + 13598 structural, maxrel 0.9" — an artifact of the greedy `FLOAT_RE` merging the last 11-col field with the trailing MAT number; a fixed-column (6×11) ENDF parse showed only ~673/213800 lines actually differ (99.7% match). **T70 tape60 remaining:** (a) MF3/MT1 471 lines at ~1e-6 = the **broadr Doppler/Sigma1 FP floor** (MT=1 copied verbatim from broadr; same class as the T01 NUMERIC-not-BI floor, Trap #14) — concrete lead (sum-of-broadened-partials vs broaden-the-total) tracked in a new P2 investigation bead, HIGH regression risk; (b) MF6/MT221 177 lines ≤5.2e-6 + MF3/MT221 14 lines ≤4.8e-7 = calcem σ/cosine FP grind (c3q); (c) the free-gas path emax endpoint has the same shape but is UNCONFIRMED → bead `NJOY_jl-nm2` (LAW 1, needs its own oracle).


**2026-06-14 (Phase 86, bead `NJOY_jl-h61` closed):** T70 tape60's **MF1/MT451 directory-NC off-by-one is fixed** — the MF3/MT221 and MF3/MT222 directory entries now read NC=498 (were 497), matching `referenceTape60` lines 69-70. T70 tape60's first structural diff moved from **line 69 → line 223** (now purely the c3q σ/cosine value-grind: `2.169385` vs `2.169386`, ±1 ulp); total 213800==ref, matching lines 128064→**128066** (+2 = exactly the two fixed directory lines). The directory-NC formula no longer uses the stale `thermr_coh_ne`; new `_thermr_mf3_dir_nc(np) = 3 + div(np+1, 3)` mirrors Fortran thermr `tpend` (`nc=3+(ne+2)/3` with `ne=NP_written-1` for the appended `etop=20e6` sentinel — thermr.f90:3087/3198/3211). **Fortran ground truth (LAW 2) OVERTURNED the bead's own `3+cld(np,3)` suggestion:** the faithful form UNDERCOUNTS the true record count by 1 when `np≡1 (mod 3)` — reproduced per the GROUND-TRUTH PRINCIPLE, NOT "fixed" to ceil (the reconr else-branch keeps `3+cld`, reconr.f90:5115-5116). New unit test `test_thermr_mf3_nc.jl` 10/10 (incl. the `np≡1 mod 3` quirk). T01 NUMERIC 32812/32962 + T68 directory NC=129 (ref-tape + unit verified) preserved. The now-vestigial `thermr_coh_ne` plumbing removal is tracked by `NJOY_jl-czw` (P3). Orchestrated: 3+1 read-only research fan-out (sonnet) + 1 opus TDD agent, serial Julia per Rule 9.

**2026-06-14 (Phase 85, bead `NJOY_jl-7dp` closed):** T70 thermr tape60's **MF6/MT221 +26-line structural overshoot is eliminated** — total now **213800 == ref**, every section line count matches exactly (MF3/MT221=498, MF3/MT222=498, MF6/MT221=185295, all 101 incident-energy blocks), matching lines 125483→128064. Root cause: calcem's E′-adaptive convergence used a simple midpoint average for the reference `ym` instead of Fortran's `terp1` law-2 lin-lin chord at the **sigfig-rounded** `xm` (thermr.f90:2051-2059 + endf.f90:1614); fixed both S(α,β)+free-gas paths in `thermr.jl`. T01 NUMERIC 32812/32962 unchanged (no regression). Both pre-session hypotheses (`half_tol`, last-seed zeroing) were wrong — a per-block NW diff decided it; `half_tol=0.5·tol` is faithful (sigl halves `tolin`, thermr.f90:2694). **T70 tape60 remaining:** (a) ~~MF1/MT451 directory-NC off-by-one for MF3/MT221+222~~ **[FIXED — Phase 86, bead `NJOY_jl-h61` closed: NC 497→498 via faithful `3+div(np+1,3)`]**; (b) the ~40% σ/cosine value grind @>1e-9 (under `c3q`; T70 tape60 first-diff now at line 223). *(Phase 84, 2026-06-09: the prior thermr MF6/MT221 cosine-width fix, 10→nbin+2 words/secondary, that took T70 112818→213826 lines; see Recent Phases table.)*

**2026-06-04 orchestrated bead session (Phase 83):** three beads landed, each Fortran-root-caused and verified.
- **NJOY_jl-o2l** — T72 **CRASH→DIFFS**. The aplots `log10(-1e-6)` crash was NOT an aplots bug: `read_mf6_incident_energies` mis-parsed MF6 **LAW=7** (nested TAB2) and harvested μ-cosine values (−1.0,−0.9) as incident energies → negative ESZ grid. Fixed LAW=7/6 dispatch (commit `c3bfd75`). `_ascll` left as a faithful Fortran port. Aplots BI cohort T50/52/53/61/62 preserved; T02/T08 unaffected (no acer step). Follow-up: T72 structural gap (bead NJOY_jl-i23).
- **NJOY_jl-lho** — T09 leapr **DIFFS→BIT_IDENTICAL 1830/1830**, and bonus **T33 DIFFS→NUMERIC_PASS** (75%@1e-9 → 99.996%@1e-5). The secondary scatterer was a red herring (b7=1 free-gas → leapr.f90:398 skips the secondary loop; primary-only is correct). Real bug: **discrete oscillators were never wired into the pipeline** (`generate_sab` called without oscillators; dormant helper was not a faithful port). Ported `discre!` + `_bfact_leapr` faithfully (commit `5132e5b`). T22 BI + T80 NUMERIC preserved. T23 (BeO, nss=1/b7=0) shifted DIFFS→STRUCTURAL_FAIL — oscillators now apply but its b7=0 SCT-secondary merge is still unported (folded into NJOY_jl-gm1).
- **NJOY_jl-3rp** — closed (ACE suffix 80c→00c fix already in commit `0da9fa2`).
- **NJOY_jl-c3q (part 1)** — T70 thermr **MF1/MT451 header now faithful** (commit `cf8e66a`). `write_full_pendf` emitted a corrupt 3-record header (LRP=0, NFOR=0, missing the AWI/EMAX/NSUB/NVER record); now `read_mf1_header_info` + `_write_full_mf1`'s `hdr` path reproduce the iverf-conditioned ENDF-6/5 header (reconr.f90:5028-5067). T70 records 1-8 byte-identical; **T01's header also fixed** (NUMERIC_PASS 32812/32962 preserved); T03/T09 BI preserved; T02/T08 untouched (legacy path). **Triage correction:** the 2026-06-02 note that T70's grid "premise was stale (ref=213800)" is WRONG — after the oracle pin alignment the current ref is **112818**, so the grid over-production (Julia 213800, ~1.9×) is REAL and is c3q's remaining (part 2) blocker: inelastic calcem MF3 grid (MT=221) + lat=10 coherent (MT=223/229) + the `sigfig(emax,7,-1)` nudge. **DEEP — open.**

**Prior baseline (2026-06-02 — fresh sweep vs pin 2c64dfb: 8 BIT_IDENTICAL, 3 NUMERIC_PASS)**

Validated against NJOY2016 pinned at `2c64dfb` (see "## Reference Oracle Pin"). **2026-06-02: the local `njoy-reference/` clone was STALE at `3c9873d` (Dec-2025 laptop oracle, EMAX=0 references) and was checked out to the pin `2c64dfb` this session** — this recovered T22 to BIT_IDENTICAL (its EMAX reference is now the nonzero pin value the committed leapr fix emits) and added tests 86/87. `njoy-reference/` is git-ignored, so alignment is a per-machine step (`git -C njoy-reference checkout 2c64dfb` or `bash scripts/setup_reference.sh`); the reference_test preflight prints `reference oracle at pin` when aligned.

**Fresh full sweep (2026-06-02, 84.8 min, `reports/REFERENCE_SWEEP.md`): 8 BIT_IDENTICAL, 3 NUMERIC_PASS, 67 DIFFS, 1 CRASH, 6 TIMEOUT, 1 NO_REFERENCE.** Up from the prior stale report's 2 BI. **Caveat — the sweep ran under CPU contention** (orphaned sweep-worker julia processes; T01 took 145s vs its usual 83s, ~1.7×): the **6 TIMEOUTs (T17/T65 broadr, T25 acer, T67/T68/T74 moder) are inflated by contention + known pre-existing perf** (broadr U-238/U-235; thermr lat=10 Bragg grid c3q) — none in code touched this session; a clean re-sweep should reduce them. The **1 CRASH (T72)** is a **pre-existing aplots bug, not a regression**: `_ascll` (ace_aplots.jl:129) calls `log10(-1e-6)` on a neutron acer-iopt7 plot curve — surfaced by the first full sweep since the Phase 79-80 aplots port (new bead filed). The BI/NUMERIC cohort below was each independently re-verified in targeted cache-nuked runs.

**Full test bit-identical** (rtol 1e-9, every produced tape):
- T03 (moder→reconr, photoatomic) — tape37 9274/9274
- T61 (acer iopt=7 thermal) — tape71 49211/49211 + tape72
- T50 (acer α+He-4) — tape33 432/432 + tape34 163/163 + tape35 — **aplots plot tape ported (Phase 79)**
- T52 (acer p+H-1) — tape33 6042/6042 + tape34 3986/3986 + tape35 — **aplots ported (Phase 79)**
- T53 (acer d+H-2) — tape33 17236/17236 + tape34 12030/12030 + tape35 — **aplots heating + aploxp ported (Phase 80)**
- T62 (acer d+He-3) — tape33 5460/5460 + tape34 7221/7221 + tape35 — **aplots aplotr threshold ported (Phase 80)**
- T22 (leapr S(α,β)) — tape20 4636/4636 — **EMAX fix (NJOY_jl-ixb); recovered to BI when the oracle was aligned to the pin 2026-06-02**
- **T09 (leapr S(α,β), H-in-H₂O) — tape24 1830/1830 — NEW 2026-06-04 (Phase 83, NJOY_jl-lho, commit 5132e5b).** Discrete-oscillator port (`discre!`/`_bfact_leapr`).
- **T86 (groupr ign=1 + card9a extended format, Hf-177 TENDL2021) — tape25 52/52 — NEW 2026-06-02 (Phase 82, beads NJOY_jl-57z + ow6, commits 675fdb5 + 2375874).** Required: ign=1 user group read + card9a/compound-mfd parse; a faithful iwt∈{2,3} panel Lobatto-2 quadrature port (the old `group_integrate` linearized 1/E → g1 flux 10× low, never validated for bit-identical GENDF); reconr MF10 (radioactive nuclide production) reconstruction onto the PENDF (gated on MF10 presence → no-op elsewhere); groupr MF10 activation averaging emitting MF3/MTd with izam/fzam in C2 + igzero threshold-group skip.

**aplots generalizes but blocked elsewhere** (Phase 80):
- T54 (acer p+H-3) — **tape34 DIFFS → NUMERIC_PASS (2026-06-02, commits 648c6cc + 8c622c7, bead NJOY_jl-53h).** Two real bugs fixed: (1) `_terpa_scr2` below-range (`ace_lcp_build.jl:601`: `e<=Ei(1)&&return Vi(1)` → `e<Ei(1)&&return 0.0`, Fortran `terpa` label 170, endf.f90:1812-1817) — zeroed the triton MT=2 recoil heating[1]; (2) AND-block LC locators written as floats (`ace_lcp_build.jl:502` missing `isint=true`; Fortran writes i20 integers, acefc.f90:9520 + acecm `typen` iflag=1) — fixed 52 locator words. tape34 1e-9 6383→6398, NUMERIC_PASS 7326/7327 (only the wildcarded date). T53/T50/T52/T62 FULL BI + T02/T08 preserved. **But T54 the TEST is still DIFFS — the HANDOFF's old "tape33 flips once the recoil word lands" was WRONG on two counts (locator bug + FP residual).** tape33 DIFFS 9919/11340: all diffs are the charged heating/elastic ~7e-7 FP-floor (T01/T34 class), with 4 outliers >1e-5 (max 1.84e-5, lines 3109/3115/3130/3132, E=3.1–6.6 MeV) in the **'recoil heating' PLOT curve** — aplots computes it as `esz_heat − Σ particle_heating`, so subtractive cancellation amplifies the underlying 7e-7 ACE-heating residual to 1.8e-5. Same residual = tape34's 949 lines off @1e-9. **Remaining (bead 53h, deep FP grind):** to flip T54→NUMERIC_PASS, drop those 4 plot outliers <1e-5 (match Fortran aplots `esz_heat−Σ` order in `ace_aplots.jl`, or reduce the 7e-7 charged-heating residual); to flip→BI, kill the 7e-7 residual. See `worklog/T54_recoil_heating_terpa_diagnosis.md` + `worklog/phase81_t54_tape34_numeric_pass.md`.
- T51 (acer d+H-2 .10h) / T71 (acer 78184.10o) — aplots no longer crashes, but tape33 is wrong-content because their ACE (tape34) is itself STRUCTURAL_FAIL (T51 unported MF6 LAW=6; T71 over-produces). ACE-level work, separate from aplots.

**Numeric pass** (rtol 1e-5):
- T01 (reconr→broadr→heatr→thermr→groupr) — tape25 32812/32962
- T80 (leapr S(α,β)) — tape24 91405/91453
- T54 (acer p+H-3) — tape34 7326/7327 (Phase 81)
- **T33 (leapr S(α,β)) — tape24 53149/53151, tape34 53145/53151 — NEW 2026-06-04 (Phase 83, NJOY_jl-lho, commit 5132e5b).** Discrete-oscillator port; was DIFFS 75%@1e-9.

Notes:
- `reports/REFERENCE_SWEEP.md` predates Phase 79 and shows T50/T52 as DIFFS — they are now **full BIT_IDENTICAL** (tape33 aplots ported + verified via targeted cache-nuked `reference_test.jl` runs); the next full sweep will reflect this. A test counts bit-identical only if EVERY produced tape passes, so always read the per-tape detail (T53/T62 still carry a deferred `tape33`, bead NJOY_jl-zsh).
- T22 (leapr) is now **BIT_IDENTICAL 4636/4636** — the MF1/MT451 record-3 EMAX was hardcoded 0; fixed to `round_sigfig(.0253·β_max, 7, 0)` per `leapr.f90:3083` (upstream PR #372). Bead NJOY_jl-ixb closed. The same one-line fix corrects the EMAX line in T09/T23/T33/T80, which carry separate value diffs (T80 FP-order → sc5; T09 MF7/MT4 S(α,β) → NJOY_jl-lho).
- RECONR is the most mature module (~19 materials bit-identical on the produced PENDF across all resonance formalisms).
- `reports/REFERENCE_SWEEP.md` is the fresh sweep; README.md was rewritten this session to reflect this baseline.

---
## Architecture

### Key source files

| File | What it does | Fortran equivalent |
|------|-------------|-------------------|
| `src/processing/reconr.jl` | Top-level RECONR pipeline. `reconr()` returns NamedTuple, `reconstruct()` returns PointwiseMaterial. | reconr.f90 main |
| `src/processing/reconr_grid.jl` | Grid construction: `lunion_grid()` (union grid builder matching Fortran lunion) | reconr.f90 lunion |
| `src/processing/reconr_evaluator.jl` | `merge_background_legacy()` adds MF3 backgrounds. `sigma_mf2()` evaluates resonance XS. | reconr.f90 emerge, sigma |
| `src/processing/pendf_writer.jl` | PENDF output. `_get_legacy_section()` handles thresholds, redundant sums. | reconr.f90 emerge/recout |
| `src/processing/adaptive_grid.jl` | Generic adaptive linearization: `adaptive_reconstruct()`, `round_sigfig()` | reconr.f90 resxs/panel |
| `src/resonances/slbw.jl` | Single-Level Breit-Wigner cross section evaluation | reconr.f90 csslbw |
| `src/resonances/mlbw.jl` | Multi-Level Breit-Wigner + Reich-Moore dispatch | reconr.f90 csmlbw, csrmat |
| `src/resonances/reich_moore.jl` | Reich-Moore R-matrix cross section evaluation | reconr.f90 csrmat |
| `src/resonances/unresolved.jl` | Unresolved resonance XS: `_csunr1`, `_csunr2`, `_gnrl`, table builder + interpolator | reconr.f90 csunr1, csunr2, gnrl, genunr, sigunr |
| `src/resonances/reader.jl` | MF2 reader: SLBW, MLBW, RM, SAMMY, URR (LRU=2 modes 11+12) | reconr.f90 rdfil2, rdf2u1, rdf2u2 |
| `src/processing/reconr_types.jl` | MF3 + MF12 readers: `read_mf3_sections`, `read_mf12_lo1_sections` | reconr.f90 lunion (MF processing) |
| `src/endf/io.jl` | ENDF I/O: `format_endf_float` (with 9-sigfig `a11` extension), `parse_endf_float` | endf.f90 a11, lineio |
| `src/constants.jl` | Physics constants in CGS matching `phys.f90` | phys.f90 |

### Fortran subroutines you'll need to read

| Subroutine | reconr.f90 lines | Purpose |
|-----------|-----------------|---------|
| `lunion` | 1771-2238 | Union energy grid from MF3+MF12 sections. Key: decade points (2112-2127), histogram shading (2050-2064), singularity shading (1930-1937, 2024-2033, 2093-2096), step ratio (2136-2141) |
| `resxs` | 2240-2569 | Adaptive reconstruction in resonance range. Key: thermal tightening (2390), step guard estp=4.1 (2419), sigfig midpoint ndig=8/9 rounding (2360-2372) |
| `sigma` | 2571-2667 | Dispatch to SLBW/MLBW/RM + unresolved (sigunr at 2659-2665) |
| `csslbw` | 2669-2856 | SLBW cross section formula. Key: sigfig(sigp(2),8,0) at line 2845 |
| `csrmat` | 2858-3023 | Reich-Moore R-matrix cross section formula |
| `emerge` | 4646-4982 | Merge grids + evaluate + write output. Key: itype dispatch (4756-4770), MF3 bg suppression in URR range (4800), sigfig(sn,7,0) at line 4832, ith pseudo-threshold tracking (4834) |
| `recout` | 4984-5441 | Write PENDF. Key: redundant reactions at line 5308 |
| `genunr` | 1628-1735 | Build unresolved XS table. Dispatches: mode=11→csunr1, mode=12→csunr2 (1682-1683) |
| `sigunr` | 1737-1769 | Interpolate from unresolved table |
| `rdfil2` | 700-881 | Read all MF2 ranges. Key: eunr boundary nodes (759-775), eunr sort+dedup (856-869) |
| `rdf2u1` | 1312-1425 | Read LRU=2 unresolved data (mode=11, energy-dependent fission widths) |
| `rdf2u2` | 1428-1540 | Read LRU=2 unresolved data (mode=12, all widths energy-dependent). 6-word stride per energy, gap filling with egridu at 1/1000 threshold |
| `csunr2` | 4079-4317 | Unresolved XS for mode=12. ALL parameters interpolated via terp1 at each energy. Clamps small values <1e-8 to 0. Final XS interpolation uses per-J INT (typically 2=linear) |
| `rdf2bw` | 884-1310 | Read SLBW/MLBW resonance parameters. Adds peak nodes to enode with sigfig rounding |
| `a11` | endf.f90:882-981 | Float → 11-char ENDF format |
| `sigfig` | util.f90:361-393 | Round to N sigfigs with shading. Key: bias=1.0000000000001 at line 390 |

### Pipeline flow (reconr.jl for LRU=1 + URR materials)

```
1. read_mf2 → MF2Data (includes resolved + unresolved ranges)
2. read_mf3_sections → Vector{MF3Section}
3. read_mf12_lo1_sections → Vector{MF3Section} (photon multiplicities)
4. build_unresolved_table → URRTable (csunr1/csunr2 + MF3 bkg at egridu nodes)
   - Includes rdfil2 boundary nodes: sigfig(EL,7,±1), sigfig(EH,7,±1)
   - Sort, dedup, overlap marking (negative energies for E < eresr)
5. _add_mf2_nodes! → MF2 peak/width nodes + URR table energies
6. lunion_grid → bg_grid (union of MF3+MF12 panels + MF2 nodes + URR nodes)
   - Shades duplicate MF3 breakpoints to sigfig(x,7,±1)
   - Shades histogram interior breakpoints to sigfig(x,7,±1) pairs
   - elim = min(0.99e6, eresr) — decade/ratio checks only below elim
   - Processes sections cumulatively (each sees previous grid)
7. filter to [eresl, eresh) → res_grid
8. adaptive_reconstruct(xs_partials, res_grid) → res_energies
   - xs_partials returns (elastic, fission, capture) — NOT total
   - force_boundaries=[eresr] prevents subdivision across resolved/unresolved boundary
   - Thermal tightening: err/5 below 0.4999 eV
   - Step guard: estp=4.1
   - Midpoint rounding: ndig=9 (or 8 for 0.1<xm<1), sigfig comparison check
9. merge bg_outside + res_energies → all_energies
10. for each energy: sigma_mf2 (resolved) + eval_unresolved (URR table) → res_xs
11. merge_background_legacy(all_energies, res_xs, mf3_sections) → final XS
    - Primary channels (MT=2,18,102): add UNROUNDED MF3 bg to resonance
    - Non-primary channels: round MF3 bg to 7 sigfigs before accumulating
    - Round each primary channel to 7 sigfigs (matching emerge line 4832)
    - Total = round_sigfig(other_bg + sigfig(el) + sigfig(fi) + sigfig(ca), 7)
12. write_pendf_file → PENDF output
```

---

## Traps and Lessons (from 6 debugging sessions)

### Critical traps

1. **Precomp cache corruption** — `rm -rf ~/.julia/compiled/v1.12/NJOY*` before EVERY run. Non-negotiable.

2. **`sigfig` bias** — `round_sigfig` multiplies by 1.0000000000001. This creates near-duplicates that `unique!` misses. Use `_dedup_tol!` instead.

3. **Constants are hardcoded in Fortran** — `ehigh = 20e6`, `elow = 1e-5`, `elim = 0.99e6`, `emax = 19e6`, `third = 0.333333333` (truncated, NOT 1/3). Never read these from MF2.

4. **Multi-material tapes** — `find_section(io, 2, 151; target_mat=MAT)` must filter by MAT. Without it, you get another material's data.

5. **CGS units** — `constants.jl` uses ergs, grams, cm/s matching `phys.f90`. A previous session changed to SI and broke everything. Don't touch.

6. **MF2 has MULTIPLE ranges** — Pu-238 has LRU=1 (resolved, 1-200 eV) AND LRU=2 (unresolved, 200-10000 eV). U-235 has SLBW (1-82 eV) + URR/LRF=2 (82-25000 eV). The reader must parse BOTH. `eresh = max(EH_resolved, EH_unresolved)`, `eresr = EH_resolved_only`.

### Subtle Fortran behaviors that must be matched

7. **SLBW sigfig rounding** — `sigfig(sigp(2),8,0)` and `sigfig(spot,8,0)` BEFORE adding potential scattering to elastic. Without this, values at 7th-sigfig boundaries round wrong (±1 in last digit). See reconr.f90:2845-2847.

8. **Adaptive convergence tests partials only** — The Fortran tests j=1..nsig-1 = (elastic, fission, capture), NOT total. Testing total makes convergence stricter, producing extra grid points. See reconr.f90:2394.

9. **Boundary-crossing forced convergence** — Panels straddling eresr (resolved/unresolved boundary) are forced to converge without subdivision. See reconr.f90:2353-2355. Also applies at eresu and eresm boundaries. Implemented via `force_boundaries` field in AdaptiveConfig.

10. **Threshold interpolation uses MF3's own law** — Near thresholds, the first breakpoint is modified to (thrxx, 0.0) and the section's own interpolation law (LinLog, LogLog, etc.) is used — NOT hardcoded linear. Implemented via `_threshold_interp()`.

11. **Total = sum of sigfig'd sections, NOT component sum** — Fortran `lunion` skips MT=1 (line 1882), so emerge never processes MT=1. The total is accumulated from all non-redundant sections, each `sigfig(sn,7,0)`'d. Then `recout` applies a final `sigfig(total,7,0)`.

12. **ENDF float format has two modes** — The Fortran `a11` uses 9-sigfig fixed-point for values in (0.1, 1e7) with genuine precision, otherwise 7-sigfig scientific. Falls back to scientific when trailing zeros indicate only 7 sigfigs. Only applied to energy values (odd positions), not XS values.

13. **URR table needs rdfil2 boundary nodes** — sigfig(EL,7,±1) and sigfig(EH,7,±1) must be in the URR energy grid. Without this, sigunr interpolates at boundary energies instead of returning exact table values.

14. **lunion skips MT=251-300** (except 261) and MT=1, MT=3, MT=101. Also conditionally skips MT=4 (if mtr4>0), MT=103-107 (if redundant), MT=18 (if mtr18>0).

15. **Duplicate MF3 breakpoints** — Encode discontinuities. Shaded to sigfig(x,7,-1) and sigfig(x,7,+1). Implemented in lunion_grid.

16. **elim = min(0.99e6, eresr)** — NOT a fixed constant. For U-235 (eresr=82), elim=82. For Pu-238 (eresr=200), elim=200. For Ni-61 (eresr=70000), elim=70000. Decade forcing and ratio checks only below elim.

17. **MF=12 histogram shading in lunion** — The Fortran lunion processes MF=12 LO=1 sections alongside MF=3. Histogram interpolation (INT=1) triggers the two-pass shading mechanism (reconr.f90:2050-2064): each interior breakpoint E → {sigfig(E,7,-1), sigfig(E,7,+1)}. For Ni-61, MF=12 MT=102 has INT=1; for U-235, all MF=12 sections use INT=2 (linear, no shading). **Check each material's MF=12 interpolation law!**

18. **Redundant MT=4 pseudo-threshold** — Fortran recout tracks first nonzero index (`ith`) and starts MT=4 from the energy before. Without this, output includes leading zeros. Implemented as pseudo-threshold skip in `_get_legacy_section`.

19. **Mode=12 csunr2 clamps small widths** — GX and GF values below 1e-8 are clamped to exactly 0 (reconr.f90:4243-4244). Mode=11 doesn't do this.

20. **Mode=12 csunr2 potential scattering condition** — Uses `if (j.le.1)` (Fortran line 4249) vs mode=11's `if (j.eq.1)`. The `.le.` means potential scattering is added for j=0 AND j=1, not just j=1.

21. **Mode=12 final XS interpolation uses per-J INT** — The INT variable in csunr2 is overwritten by each J-state (line 4219). The final terp1 call at line 4312 uses the LAST J-state's INT, typically 2 (linear). Mode=11 hardcodes log-log (intlaw=5).

22. **MT=18 is redundant when MT=19 exists** — Fortran `anlyzd` sets `mtr18=1` when MT=19 is found (line 557-561). Then `lunion` skips MT=18 (line 1893), `emerge` never sees it, and `recout` outputs MT=18 as a computed sum of MT=19+20+21+38. Without this, fission XS is exactly 2x at low energies where MT=18 = MT=19.

23. **SLBW/MLBW half-width uses GT from ENDF, not sum of partials** — Fortran `rdf2bw` uses `hw = res(jnow+2)/2` where `res(jnow+2)` is GT (total width stored directly in ENDF). The Julia was using `(|GN|+|GG|+|GF|)/2`. These differ because GT may include competitive width and floating-point rounding. Fixed by adding GT field to SLBWParameters/MLBWParameters.

24. **MF=13 processed in lunion alongside MF=3** — Fortran lunion processes MF=3, MF=10, MF=12, MF=13, MF=23 sections (line 1868). **WARNING**: Fortran line 1880 (`scr(5)=1`) is DEAD CODE — `tab1io` at line 1913 overwrites `scr(5)` with the actual NR. MF=13 does NOT force histogram interpolation; it uses the ENDF's own interpolation law (INT=2 for U-235). Non-MF=3 sections bypass the MT skip checks (line 1881: `if (mfh.ne.3) go to 180`).

25. **Initial vs mid-data discontinuity shading differs** — Fortran lunion's initial discontinuity check (label 207, line 1979) uses `sigfig(er,7,0)` (round without nudge). The mid-data discontinuity check (label 270, line 2029) uses `sigfig(er,7,-1)` (nudge down). These produce DIFFERENT values. For U-235 at 1.09e6 eV, MF=12 has a mid-data duplicate → 1089999, while MF=13 has an initial duplicate → 1090000. Julia must distinguish `k == start_k` (initial) from `k > start_k` (mid-data) in the breakpoint insertion loop.

26. **Primary fission channel is MT=19, not MT=18+19+20+21** — In `merge_background_legacy`, only MT=18 or MT=19 should be added to the fission accumulator (whichever is the primary). MT=20, 21, 38 are non-primary: their backgrounds go to `other_bg` with sigfig rounding, matching Fortran emerge `itype` dispatch (only MT=2, MT=18/19, MT=102 get `itype` assignments).

27. **Coincidence shading comparison uses sigfig(|eg|,7,-1), NOT |eg|** — Fortran lunion label 220 (line 1996) compares `(er - sigfig(|eg|,7,-1)) > 1e-8*er`. The comparison is against the shaded-DOWN version of the grid point, not the grid point itself. This means coincidence only fires when a section's first breakpoint matches a previously-shaded grid point (where `sigfig(eg,7,+1)` then `sigfig(sigfig(eg,7,+1),7,-1) ≈ eg_original`). For two sections sharing the exact same unshaded breakpoint, the check does NOT fire because the sigfig(7,-1) subtraction creates a gap of ~1e-6*er which is much larger than the 1e-8*er threshold.

28. **Frobenius-Schur matrix inversion for Reich-Moore** — Fortran `csrmat` (reconr.f90:3503-3607) inverts the complex (I+R+iS) matrix using the Frobenius-Schur decomposition into real operations: `frobns` → `thrinv` (symmetric matrix inverse via Gaussian elimination with internal I-D transform) → `abcmat` (3x3 matrix multiply). Julia must use this exact algorithm (`_frobns`, `_thrinv!`, `_abcmat` in `reich_moore.jl`) instead of `inv(SMatrix{3,3})` to match Fortran's intermediate FP rounding. Key detail: `thrinv(D)` computes `D^{-1}` (not `(I-D)^{-1}`) despite internally negating D and adding identity — the Gaussian elimination undoes the transform.

---

## Immediate Next Steps — PRIORITY ORDER

**Use the Fortran debugger** — the NJOY2016 binary is compiled with debug symbols at `njoy-reference/build/njoy`. Patch `reconr.f90`/`samm.f90` with `write(*,...)` diagnostics and recompile for speed. Restore clean source with `cd njoy-reference && git checkout -- src/`.

### 0. PRIORITY: Generate Oracle Caches for All 84 Tests

**Why this matters**: Only 30 of 84 tests have oracle caches. Without oracles, we can't verify correctness. Many "untested" tests may already work perfectly — we just don't know.

**The oracle system**: `test/validation/diagnose_harness.jl` generates per-module Fortran reference output. For each test, it runs the Fortran NJOY binary with a truncated input deck (stopping after RECONR or BROADR) and caches the ASCII PENDF output. Oracles are stored in `test/validation/oracle_cache/testNN/`.

**How to generate an oracle**:
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/diagnose_harness.jl <test_number>
```
This creates `test/validation/oracle_cache/testNN/after_reconr.pendf` and `run_reconr/` with the input deck and tapes.

**Tests that need oracles** (16 tests, no oracle cache at all):
T24, T28, T29, T31, T32, T35, T36, T37, T38, T39, T40, T41, T42, T43, T44, T63

**Tests with oracles but not yet compared** (newly discovered in Phase 17):
| Test | MAT | err | Material | Status | Tape Issue |
|------|-----|-----|----------|--------|------------|
| T17 | 9237 | 0.001 | U-238 JENDL | No tape20 | Needs `moder` to create ASCII tape |
| T21 | 2637 | 0.001 | Fe-58 ENDF/B-8 | **54/79** | Use tape20. Grid shortfall (34k vs 50k) — adaptive reconstruction density in dense RM region (262 resonances, err=0.001) |
| T65 | 9228 | 0.001 | U-235 ENDF/B-8 | **42/87** | Use tape20. XS precision diffs, likely URR boundary class |

**Binary tape issue**: Some tests use `moder` to convert ASCII→binary before reconr reads. The oracle caches store BOTH tape20 (ASCII original) and tape21 (binary from moder). Julia's reconr reads ASCII only. **Always use tape20** (ASCII) when available. If only tape21 (binary) exists, you need to run `moder` first or find the ASCII source.

Affected tests: T16 (tape20 works), T17 (no tape20 — needs moder), T58/T60/T64/T65 (tape20 works)

**Photonuclear tests (WILL CRASH — expected)**:
T03, T56, T57, T58, T64 — these use photon-induced ENDF files with NO MF2/MT151. reconr currently requires MF2. These need a new code path (MF=23 processing) and are NOT expected to pass. Skip them for now but note the error.

**Dosimetry test (T60 — produces 0 points)**:
Fe-nat IRDFF-II has no MF=3 sections, only MF=10 and MF=40. reconr produces 0 energy points. Needs investigation — the Fortran somehow produces 1 MT.

**How to run a comparison after generating oracle**:
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e '
using NJOY
# Read oracle input to find MAT and err (see input deck in oracle_cache/testNN/run_reconr/input)
r = reconr("test/validation/oracle_cache/testNN/run_reconr/tape20"; mat=MAT, err=ERR)
write_pendf_file("/tmp/testNN.pendf", r; mat=MAT, err=ERR)
# Then compare MF3 columns 1-66
'
```

**What to record for each test**: Run status (OK/ERROR/CRASH), grid size (Julia vs Fortran), MTs PERFECT/total, error classification (grid diff, XS ±1, missing MTs, crash reason).

**Goal**: Build a complete test matrix showing which tests run, which pass, and what the failure mode is for each. This lets us identify which code paths are brittle and prioritize fixes.

### 1. RECONR: T20 (Cl-35 RML/SAMMY, 158/162) — NEAR COMPLETE

**Status**: 158/162 PERFECT. Grid: 10730 pts (exact match). Only 4 MTs with ±1 FP diffs.

**Phase 17 fixes (this session — supersede Phase 16 grid/MT=103 analysis)**:

**Bug 7 — FIXED (SAMMY peak/width nodes)**: `_add_peak_nodes!` had NO dispatch for `SAMMYParameters` — the fallback did nothing, so SAMMY resonance peak/width nodes were never added to the initial grid. Added `_add_peak_nodes!(nodes, params::SAMMYParameters, el, eh)` matching Fortran `rdsammy` (samm.f90:1151-1177): for each resonance, adds {Er-hw, Er, Er+hw} with hw = gamgam/2 + Σ|gamma_j|/2. Fixed 675-point grid discrepancy.

**Bug 8 — FIXED (per-section AWR for thresholds)**: Cl-35 ENDF has MF2 AWR=34.66850 but MF3 AWR=34.66845. Julia passed `mf2.AWR` to `lunion_grid`, but Fortran reads AWR from each MF3 section's HEAD record (`awrx=c2h/awin` at reconr.f90:1911). Added `awr` field to `MF3Section` struct, populated from HEAD record. Now `lunion_grid` and PENDF writer use per-section AWR. Fixed 10 threshold energies rounding to wrong sigfig boundaries. **Confirmed via gdb**: Fortran `thrx=1.317612e6` for MT=807 (using MF3 AWR=34.66845), Julia was computing `thrx=1.317611e6` (using MF2 AWR=34.66850).

**Bug 9 — FIXED (reaction XS in redundant sums)**: MT=103-107 redundant sums computed only MF3 backgrounds, missing SAMMY `reaction_xs` contributions. MT=103 showed zeros instead of 24.1 barns at thermal. Added `reaction_xs[smt]` to each partial's contribution in the sum. Fixed 10 threshold sections' pseudo-threshold start points.

**Previous bugs (still fixed)**: Bug 1 (2x doubling), Bug 2 (Coulomb l=0), Bug 3 (XQ indexing), Bug 4 (reaction channel plumbing), Bug 5 (l>0 Coulomb), Bug 6 (4-channel convergence)

**Remaining issues for T20 (4 MTs)**:
1. **MT=1 (total)**: 2 differing lines at E≈1.085e6 and E≈1.198e6 — ±1 in 7th sigfig. Resonance XS precision.
2. **MT=2 (elastic)**: 1 differing line at E≈1.085e6 — ±1. Likely from hard-sphere vs Coulomb phase shift (see below).
3. **MT=600 (proton)**: 1973 differing lines — all ±1 in 7th sigfig. IEEE 754 non-associativity in R-matrix accumulation over ~200 resonances (same class as T34 irreducible diffs). Each `|XXXX_ij|²` computation has ~1e-13 intermediate differences that cascade through sigfig rounding boundaries.
4. **MT=103 (proton total)**: Same 1973 lines as MT=600 (MT=103 = sum of MT=600-649; only MT=600 contributes below 1 MeV).

**Remaining known correctness issue**: Coulomb phase shift (Julia uses hard-sphere `_sammy_sinsix` for all channels; Fortran uses `pghcou` for Coulomb channels). Impact: entrance-channel elastic only (MT=2). The proton channel is exit-only, so its phase shift is unused. Would require implementing Coulomb phase shift extraction from `_coulomb_pen_shift` and making two `pghcou` calls when rho≠rhof.

**Trap 37 (FIXED)**: Per-section AWR — ENDF files can have different AWR values in MF2 vs MF3 HEAD records. Always use AWR from the MF3 section being processed for threshold computation, not the MF2 AWR.

**Trap 38 (FIXED)**: SAMMY resonances need peak/width nodes — unlike BW/RM which have explicit `_add_peak_nodes!` dispatches, SAMMY previously fell through to the do-nothing default. Without these seed nodes, the adaptive reconstruction missed narrow resonances entirely.

### 2. RECONR: Fix T34 ±1 FP diffs (52/53) — CONFIRMED HARD

**Status**: 52/53. 3 remaining ±1 diffs in MT=102 (capture) at E=630.04, 2089.07, 4526.46 eV.

**Phase 14 gdb confirmation**: Fortran values traced via diagnostic prints in csrmat (reconr.f90:3496-3499, after pifac multiplication). Results:
- E=630.04: Fortran cap=0.096262384**997**, Julia cap=0.096262384**999** (diff=+2.3e-12)
- E=2089.07: Fortran cap=0.041812484**989**, Julia cap=0.041812484**999** (diff=+1.0e-11)
- E=4526.46: Fortran cap=0.004613320**4999**2, Julia cap=0.004613320**4998** (diff=-1.3e-13)

All values within 1e-11 of the 0.5 boundary at 7 sigfigs. Both codes use identical Frobenius-Schur algorithm (`_frobns`, `_thrinv!`, `_abcmat`). The difference is purely from IEEE 754 non-associativity in the accumulation of `gg4*a1*a1` over 437 l=0 fissile resonances.

**Possible approach (not yet tried)**: Match the Fortran's exact loop order for the R-matrix accumulation. The Fortran accumulates `r(1,1) += gg4*a1*a1` sequentially over all resonances with matching J-value. Julia does the same loop but the IEEE 754 intermediates may differ due to compiler optimizations. Try: (a) force sequential accumulation with `@fastmath false`, (b) match Fortran's `per=res(in+1)` precomputed penetrability instead of on-the-fly computation (the `in` index in csrmat line 3340 reads precomputed values from `rdf2bw` line 1002).

### 3. RECONR: T49 (Zr-90) — URR sentinel FIXED, 44/46 (2 FP diffs remain)

**Status (2026-04-17)**: **44/46** after URR upper-boundary sentinel fix
(was 41/46). Remaining 2 MTs (MT=1, MT=2) differ by ±1 in 7th sigfig at
E=110487.7 eV — MLBW elastic evaluation FP class. No regressions across 16
other bit-identical tests.

**Root cause (corrected — HANDOFF's prior Phase 14 analysis was wrong)**:
`_add_mf2_nodes!` added `sigfig(EH_URR,7,+1) = 1780461` as the last enode.
In Fortran lunion, the last enode is a *sentinel* — label 240's check
`ig.ge.ngo` at reconr.f90:2051 routes control past `en=abs(eg)`, so the
last enode never becomes a sub-panel boundary. It only reaches `inew` if
some MF3 section contributes it via per-section processing (e.g. Pu-239
T27 has MT=2 x[87]=30000.01; Pu-238/Cf-252/Pu-240 have duplicate EH
breakpoints that lunion_grid shades). Zr-90 has neither — MT=2 has only
`x[78]=1780460` — so Fortran drops 1780461 but Julia kept it.

HANDOFF's prior "MT=2's 1780460 gets absorbed into the shaded pair"
diagnosis was wrong: Julia already lacks 1780460 (only Fortran writes it,
then removes it at label 410 via `eresm` dedup). The actual bug was the
spurious 1780461.

**Fix (landed)**: `_drop_unsupported_urr_plus_boundary!` in
`src/processing/reconr_grid.jl`, invoked by `reconr.jl` after
`lunion_grid`. For each LRU=2 range, if `sigfig(EH,7,+1)` is in grid AND
no MF3 section has (a) an explicit breakpoint at that energy OR (b) a
duplicate pair at `EH` that lunion_grid's shading would convert into
`sigfig(EH,7,+1)`, remove it.

Worklog: `worklog/T49_urr_sentinel.md`.

**Remaining**: MT=1 (6.504111 vs 6.504112) and MT=2 (6.497532 vs 6.497533)
at E=110487.7 eV. Candidate for Fortran diagnostic on MLBW elastic eval
(same approach that found T34's "irreducible" ±1 diffs were real bugs).

### 4. RECONR: Fix T46 MT=1 total (72/73) — SUMMATION ORDER

**Status**: 72/73. MT=1 (total) differs by ±1-3 at 2 high-energy lines.

All individual MTs are BIT-IDENTICAL at the diff energies. The total differs because Julia computes `total = sigfig(other_bg + elastic + fission + capture, 7)` while Fortran accumulates each sigfig'd section in tape order (emerge line 4893: `tot(2)=tot(2)+sn`). IEEE 754 non-associativity of the summation causes ±1.

**Fix approach**: Rewrite the total computation in `merge_background_legacy` to accumulate in section order (iterate `mf3_sections`, add each sigfig'd value to the total). Previous attempt regressed because it re-evaluated MF3 backgrounds in the total loop without matching all the same conditions (threshold, LSSF suppression). The correct implementation must use the ALREADY-COMPUTED sigfig'd values for each section, not re-interpolate. Store each section's contribution in a vector, then sum in order.

### 5. RECONR: Fix T04/T07 ±1 FP diffs (24/27) — SAME CLASS AS T34

**Status**: 24/27 each. 3 MTs (MT=18/19/102) differ by ±1 at E=82.00001 eV (URR boundary).

Same class as T34: `_gnrl` Gauss-Laguerre quadrature accumulates 100 terms. The raw values are at the sigfig rounding boundary. Not investigated with gdb yet — the same diagnostic approach as T34 applies.

### 6. RECONR: Fix T15 ±1 FP diffs (32/36)

**Status**: 32/36. 4 MTs (1,2,18,102) with ±1 diffs, 59 total differing lines. Mostly at 7th sigfig, plus one energy value diff at E=9238.377 (9-sigfig format: Julia produces 9238.37706, Fortran 9238.37707).

The energy diff suggests a midpoint rounding issue in the adaptive reconstruction — the `sigfig(xm, ndig, 0)` at reconr.f90:2369-2371 may differ slightly from Julia's `round_sigfig`.

### 7. Grind BROADR to bit-identical

BROADR is fully implemented in `src/processing/broadr.jl` and `src/processing/sigma1.jl`. BROADR oracles exist for 13 tests. Apply the Grind Method.

### 8. Fix remaining RECONR test errors

- **T03, T56-58, T64**: photoatomic/photonuclear files need MF=23 processing in lunion (no MF2/MT151)
- **T60** (Fe-nat IRDFF): dosimetry file with no MF=3 (only MF=10 + MF=40), needs new code path
- **T21** (Fe-58): adaptive reconstruction density shortfall in dense RM region. 34k vs 50k grid points. Root cause unclear — peak nodes are present but adaptive reconstruction converges earlier than Fortran's
- **T65** (U-235 ENDF/B-8): 42/87, 43 XS-only diffs + 2 grid diffs. Likely URR boundary precision class

---

## How to Run T01

The T01 pipeline test script is committed at `test/validation/t01_pipeline.jl`. Run:
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/t01_pipeline.jl
```

This produces `/tmp/t01_tape25.pendf` and compares section-by-section and line-by-line against `njoy-reference/tests/01/referenceTape25` (32,962 lines).

**What the test script does** (read it before modifying!):
1. **reconr**: Runs reconr on tape20 → gets 0K pointwise XS (r.energies, r.elastic, r.capture, r.total)
2. **broadr**: Runs `broadn_grid(r.energies, hcat(r.elastic, r.capture), alpha, ...)` with PARTIALS ONLY (Trap 53). Then broadens total via sigma1_at on the same grid.
3. **heatr**: Reads MF12/MT=102 gamma data from ENDF, computes KERMA via `compute_kerma` with `gamma_data` dict.
4. **thermr1 (free gas)**: MT=221 = broadened elastic below emax (Trap 50). MF6/MT=221 via `calcem_free_gas` (returns `(esi, xsi, records)` — xsi used for MF6 normalization).
5. **thermr2 (SAB)**: Reads SAB from t322. Builds thermal grid via `build_thermal_grid` using BROADENED elastic grid (NOT reconr grid — broadr thins 17 points). Computes MT=229 via `calcem`, MT=230 via `build_bragg_data`. Grid has 571 points (569 from coh + emax + sigfig(emax,7,+1) + 2e7 sentinel). Calcem xsi interpolated to thermal grid via ORDER-5 Lagrangian with `terp_lagrange`.
6. **write**: `write_full_pendf` with override_mf3, extra_mf3, MF6 records (with xsi normalization via `mf6_xsi` and `mf6_emax` parameters), MF12/13 passthrough, description text, thermr_mts=Set([221,229,230]).

**Key parameters** (must match `test/validation/oracle_cache/test01/run_*/input`):
- reconr: tape20, mat=1306, err=0.005
- broadr: alpha=AWR/(bk*296), tol=0.005, thnmax=4.81207e6
- heatr: Z=6 (carbon), gamma_data from MF12/MT=102
- thermr1: emax=1.2, nbin=8, tol=0.05, natom=1, iinc=1 (free gas), **sigma_b=((A+1)/A)^2** (Trap 58), **awr=A passed to sigl** (Trap 73)
- thermr2: emax=1.2, nbin=8, tol=0.05, natom=1, iinc=2 (SAB), MAT_sab=1065
- bragg: a=2.4573e-8, c=6.7e-8, sigma_coh=5.50, A_mass=12.011, natom=1, DW=2.1997, lat=1
- MF6 normalization: sigma divided by xsi(ie) in _write_mf6_section (Trap 74)

## How to Use gdb / Fortran Diagnostics (PROVEN EFFECTIVE)

This project has found 8+ major bugs using Fortran diagnostic prints. The approach:

```bash
# 1. Patch the Fortran source with write(*,*) diagnostics
#    ALWAYS use write(*,*) (list-directed) — NOT formatted writes
#    Filter condition: if (iinc.eq.2.and.ie.eq.36) to limit output
#    Example: write(*,*) 'TAG',j,x(i),y(1,i)
vi njoy-reference/src/thermr.f90  # or broadr.f90, heatr.f90, reconr.f90

# 2. Rebuild (fast, ~5 seconds)
cd njoy-reference/build && cmake --build . --target njoy

# 3. Run with the appropriate oracle input
cd test/validation/oracle_cache/test01/run_thermr_2   # or run_broadr, run_heatr
../../build/njoy < input 2>&1 | grep 'TAG' > /tmp/fortran_diag.txt

# 4. Compare with Julia computation at the same energy
# Write a Julia script that prints the same intermediate values
rm -rf ~/.julia/compiled/v1.12/NJOY* && julia --project=. /tmp/diag.jl

# 5. ALWAYS restore clean source when done
cd njoy-reference && git checkout -- src/
cd build && cmake --build . --target njoy
```

**Oracle cache directories**: Each `test/validation/oracle_cache/test01/run_XXX/` has the exact Fortran input deck and tapes needed to reproduce that processing stage in isolation:
- `run_reconr/` — reconr input with tape20 (ASCII ENDF)
- `run_broadr/` — broadr input reading from reconr output
- `run_heatr/` — heatr input
- `run_thermr/` — first thermr (free gas)
- `run_thermr_2/` — second thermr (SAB from tape26=t322)

**Proven diagnostic patterns (from Phases 19-23)**:
1. **Side-by-side E' trace**: Patch Fortran calcem label 360 to print each accepted `(j, E', sigma)`. Write Julia script that does the same. Compare sequences to find first divergence. Found: missing E'=0 seed (Trap 56).
2. **Interpolation trace**: Patch Fortran tpend to print `terp` result at specific energies. Compare with Julia's interpolation. Found: terp uses order-5 Lagrangian, not linear (Trap 57).
3. **Kernel value trace**: Patch Fortran `sig` to print `(E, E', mu, alpha, beta, s, sig)`. Compare with Julia's `sab_kernel`. Found: bilinear vs biquadratic interpolation (Trap 54).
4. **Variable trace at specific function**: Patch broadr.f90 `broadn` at label 120 to print slope variables. Found: nreac=2 (partials only, Trap 53).
5. **Bragg edge trace**: Patch sigcoh to print `(i, tau_sq, ff, tau_sq*recon)` after sorting. Compare with Julia's `build_bragg_data`. Found: missing tsqx merge threshold + one-sided merge (Trap 59).
6. **Convergence decision trace**: Patch calcem do-350 loop to print `(k, xm, yt(k), ym, test2, 'REJECT'/'PASS')` for specific ie/iinc. Found: sigma test at xm=4.986e-5 makes identical borderline decisions (both reject cosine k=3 at |diff|=0.052>tol=0.05), but sigma at E'≈E differs by 57% due to amin floor (Trap 61).
7. **Parameter trace**: Print `sb`, `smz`, `az` in calcem initialization. Found: free gas sigma_b=((AWR+1)/AWR)^2, not A*sigma_free (Trap 58).

## Key Files for T01 Pipeline

| File | What it does | Key functions |
|------|-------------|---------------|
| `src/processing/broadr.jl` | Doppler broadening | `broadn_grid` (convergence stack + slope tracking) |
| `src/processing/sigma1.jl` | Doppler kernel | `sigma1_at` (h-function/f-function Voigt integral) |
| `src/processing/heatr.jl` | KERMA computation | `compute_kerma`, `photon_recoil_heating`, `photon_recoil_damage` |
| `src/processing/thermr.jl` | Thermal scattering | `calcem`, `calcem_free_gas`, `sigl_equiprobable`, `build_bragg_data`, `build_thermal_grid`, `read_mf7_mt4` |
| `src/processing/pendf_writer.jl` | PENDF output | `write_full_pendf`, `_write_mf6_section`, `_write_mf6_coherent_stub` |
| `src/processing/reconr.jl` | Resonance reconstruction | `reconr`, `reconstruct` |
| `test/validation/t01_pipeline.jl` | T01 test script | Full pipeline + comparison |

---

