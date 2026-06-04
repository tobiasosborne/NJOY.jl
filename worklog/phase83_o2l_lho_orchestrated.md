# Phase 83 — Orchestrated bead session: o2l (T72 crash) + lho (T09 leapr discrete oscillators)

**Date:** 2026-06-04
**Outcome:** +1 BIT_IDENTICAL (T09), +1 NUMERIC_PASS (T33), T72 CRASH→DIFFS. Two beads fixed (o2l, lho), one closed (3rp), two follow-ups filed/updated (i23, gm1).
**Commits:** `c3bfd75` (o2l), `5132e5b` (lho).
**Orchestration model:** parent (Opus) plans + monitors + commits; coding delegated to Opus subagents (think-hard / ultrathink); Fortran/Julia research delegated to parallel read-only Sonnet subagents (3+1 pattern); all Julia serial; verified-before-trust (parent independently re-ran T02, T09, T22, T23, T33).

---

## NJOY_jl-3rp — ACE suffix default 80c→00c (closed, no new code)

The 2-line fix (`input_parser.jl:298,326`, Fortran acer.f90:292 suff=0→`.00`) was already committed in `0da9fa2` and present in the working tree. Verified and closed. Correctness freebie; flips no test (T14 milestone is a separate deep MF6 LAW=1 task).

---

## NJOY_jl-o2l — T72 aplots `log10(-1e-6)` crash → DIFFS (commit `c3bfd75`)

### Symptom
T72 (MAT=425; reconr→broadr→heatr→groupr→acer, with a 2nd acer = iopt=7 plot pass) crashed with `DomainError(-1.0e-6)` in `_ascll` (`src/formats/ace_aplots.jl:129`, `log10(amin)`).

### Root cause (NOT in aplots)
`read_mf6_incident_energies` (`src/endf/readers.jl`) treated MF6 **LAW=7** (laboratory angle-energy, ENDF-6 §6.2.7) like the flat LAW=1/5 structure — read a CONT and skipped ⌈N1/6⌉ lines per record. LAW=7 is **nested**: outer TAB2 over NE incident energies, then for each E an *inner TAB2* over NMU cosines (its C2 = the incident energy) followed by NMU TAB1 records. The flat parse walked into the per-cosine TAB1 headers and harvested μ-cosine values (−1.0, −0.9, …) as bogus "incident energies." For T72 (Be-9 MT=16 n,2n, LAW=7) this seeded 9 negative points into the ACER `master_e`/ESZ energy grid; aplots then took `log10` of a negative energy.

`_ascll` is a faithful port of Fortran `ascll` (acefc.f90:19684-19724), which itself has no guard against `amin≤0` because valid ACE energy grids are always positive. So the fix belongs upstream, not in `_ascll`.

### Fix
Split the LAW dispatch in `read_mf6_incident_energies`:
- `law in (1,5)`: unchanged (outer TAB2 → NE LIST records, C2=E, skip ⌈NW/6⌉).
- `law == 7` (new): outer TAB2 → per E read inner TAB2 (push its C2), skip NMU TAB1.
- `law == 6` (new): single CONT, no incident-energy grid → emit none.
- `law == 2`: documented no-op (unexercised).

Ref: njoy-reference/src/acefc.f90:6363-6366,6455-6462 (acensd law==7 nmu loop); ENDF-6 §6.2.7.
ESZ grid now starts at 1.0e-11 MeV matching the reference.

### Verification / regression
- T72: CRASH → DIFFS (runs to completion). tape41 1711 vs ref 50748 — large structural gap remains (URR/probability-table expansion + missing MT pages) → **follow-up bead NJOY_jl-i23**.
- Aplots BI cohort: T50, T52, T53, T61, T62 all still BIT_IDENTICAL.
- T08 DIFFS confirmed pre-existing (stash test, identical MF1/MT451 first-diff).
- T02 unaffected — its deck has no `acer` step (moder→reconr→broadr→unresr→groupr→ccccr→moder→moder); `read_mf6_incident_energies` is called only from `acer.jl`. T02 tape29 6213-vs-1295 is a pre-existing groupr GENDF structural issue (baseline `2/1295`, reproduced exactly).

---

## NJOY_jl-lho — T09 leapr discrete oscillators → BIT_IDENTICAL (commit `5132e5b`)

### Symptom
T09 (leapr, H-in-H₂O, MAT=101, lat=1) tape24 DIFFS: 194/1830 lines differed; 4700/10392 S(α,β) fields >1e-3, max rel **565%**, growing with α, present even at the phonon peak (~4e-4). NOT FP-order.

### Diagnosis path (the value of LAW 2)
First hypothesis from the bead + a Julia warning ("nss=1 secondary scatterer not yet ported") was the **secondary scatterer**. The 3+1 Fortran research **disproved** this: Fortran leapr.f90:398 (`if (nss.eq.0.or.b7.gt.zero.or.isecs.gt.0) idone=1`) means for **b7=1 (free gas — T09's case)** the secondary loop never runs; the secondary is ENDF B-array metadata only and primary-only output IS correct. The warning was misleading.

Sharpened by comparing T09 to the passing T80 (NUMERIC_PASS @1e-5, also lat=1): T80 = diffusion translation + **0 discrete oscillators**; T09 = free-gas translation + **2 discrete oscillators**. Diff pattern (per-β-block, worst rows evenly spaced by 24 = same slot each block, α-growing) pointed squarely at the discrete-oscillator assembly.

### Root cause
The discrete oscillators were **never applied**: `leapr_module` called `generate_sab` without the oscillators and ran no `discre` analog, so any deck with nd>0 silently dropped its oscillator deltas. The dormant `_add_discrete_oscillators!` helper was also NOT a faithful port (exact `SpecialFunctions` Bessel vs NJOY's A&S-polynomial I₀/I₁ + reverse recursion; wrong interpolation path; missing line sort+cull and the twt≤0 elastic block).

### Fix (`src/processing/leapr.jl` + `src/orchestration/modules/leapr.jl`)
- `_bfact_leapr` — faithful port of `bfact` (leapr.f90:1663-1796): A&S-polynomial bessi0/bessi1, downward recursion normalized against I₁, y>1 scaled-Bessel branch folding +x into exponents.
- `discre!` — faithful port of `discre` (leapr.f90:1320-1661) on the column-major ssm[β,α,T] layout (matching accumulation order): oscillator setup, per-α delta-line accumulation, descending sort + small-weight cull, `sint` continuum convolution (`be=-betan-bes`, `wt=tbeta+twt`), twt≤0 elastic delta, `tempf=(tbeta+twt)*tempf+tsave`.
- Wiring: per-T `dwpix[itemp]=f0` (leapr.f90:715); call `discre!` after `trans!` when nd[itemp]>0 (Fortran order contin→trans→discre, leapr.f90:377-383). Narrowed the nss warning to the genuinely-unported `b7≤0` SCT-secondary case.

### Verification / regression (all cache-nuked, parent-verified)
- **T09 tape24 DIFFS → BIT_IDENTICAL 1830/1830** (ALL PASS @1e-9).
- **T33 tape24/tape34 DIFFS (75%@1e-9) → NUMERIC_PASS 53149/53151, 53145/53151 @1e-5** (bonus).
- T22 BIT_IDENTICAL 4636/4636 preserved (lat=0, nd=0).
- T80 NUMERIC_PASS preserved (diffusion, nd=0).
- T23 (BeO, nss=1/b7=0) DIFFS→STRUCTURAL_FAIL: oscillators now apply correctly (line counts for T09/T33 are exact, so `discre!` is not over-producing) but its **b7=0 SCT-secondary merge is still unported** → folded into **NJOY_jl-gm1** (second-pass arat=aws/awr + endout merge ssm=(sbs/sb)·ssm_sec+ssm_prin, leapr.f90:3018-3030; tempf1 secondary teff, leapr.f90:3578-3617). T23 was 0.02%-matching pre-fix.

### Lesson
The misleading "not yet ported" warning nearly sent the fix down the wrong path. Fortran-before-Julia (LAW 2) + the T09-vs-T80 component diff (which primary pieces does the failing test exercise that a passing test doesn't?) localized the real bug. The dormant-but-unwired helper is a recurring trap — check the call graph, not just the presence of a function.

---

## NJOY_jl-c3q (part 1) — T70 thermr MF1/MT451 faithful header (commit `cf8e66a`)

### Symptom
T70 (moder→moder→reconr→broadr→thermr→moder; compared tape = thermr output tape60) had a corrupt MF1/MT451 header.

### Root cause
`_write_full_mf1` (`src/processing/pendf_writer.jl`) wrote a hardcoded **3-record** header regardless of ENDF version: record 1 with LRP=0, a blank CONT (NFOR=0), then the TEMP/ERR control record. For ENDF-6 inputs the AWI/EMAX/NSUB/NVER record (record 3) was **entirely missing** and LRP/NFOR were wrong. Fortran THERMR copies records 1-3 of MF1/MT451 from its input PENDF verbatim (thermr.f90:3041-3061/3140-3154), which RECONR builds per reconr.f90:5028-5067 (iverf=6 → LRP=2, NFOR=6, record-3 AWI/EMAX/NSUB/NVER).

### Fix (minimal blast radius)
- New `read_mf1_header_info(endf_path, mat)` — reads the original ENDF MF1/MT451 records 1-3 with iverf detection mirroring reconr `ruina` (reconr.f90:219-256: N1≠0⇒v4, N2==0⇒v5, else v6). Returns `nothing` if absent → graceful fallback to the legacy 3-record path.
- `_write_full_mf1` gains `hdr` kwarg → emits the faithful iverf-conditioned header (ENDF-6: 4 records incl. LRP=2/NFOR=6/AWI/EMAX/NSUB/NVER/LDRV=1; ENDF-5: 3 records), each cited to reconr.f90. Self-NC dir entry `nwd+nxc+2` (ENDF-6 carry-forward, thermr.f90:3156 + reconr.f90:651,5072) vs `nwd+nxc` (ENDF-5/legacy).
- Plumbed via `RunContext.mf1_header` (read in the reconr branch, pipeline.jl) → `final_assembly!` (moder.jl) → `write_full_pendf` → `_write_full_mf1`.
- **`_write_legacy_mf1` (reconr-only T02/T08) untouched.**

### Verification
- T70 tape60 MF1/MT451 records 1-8 now **BYTE-IDENTICAL** to reference (was corrupt); 25562/112818 lines match; first whole-tape diff moved off the header to MF3/MT221.
- T01 NUMERIC_PASS 32812/32962 preserved; its ENDF-5 header now byte-identical.
- T03 BIT_IDENTICAL 9274/9274 preserved (pipeline header-read safe); T09 BIT_IDENTICAL 1830/1830 preserved.

### Triage correction (important)
The 2026-06-02 c3q triage claimed the grid-over-production premise was stale ("Julia 213800 == local ref 213800; 112817 from an older clone"). **That is WRONG.** After the oracle was aligned to pin `2c64dfb` (2026-06-02), the current T70 referenceTape60 = **112818** lines; Julia produces **213800** (~1.9× over). The ORIGINAL bead premise (grid over-production) is CORRECT and is c3q's remaining **part 2** blocker: the inelastic S(α,β) calcem MF3 grid (MT=221) + lat=10 coherent (MT=223/229) + the missing `sigfig(emax,7,-1)` nudge in `sigcoh`. DEEP — left open.

### Lesson
Re-verify a "stale premise" note against the *current* pinned oracle before trusting it (Rule 2). A bead's own triage can drift from the reference when the oracle changes underneath it.
