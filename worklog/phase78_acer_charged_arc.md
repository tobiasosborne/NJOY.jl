# Phase 78 — ACER charged-particle arc + thermr lat=10 (2026-05-29)

Orchestrated, subagent-driven session. Coding delegated to Opus subagents
(read Fortran → oracle-driven TDD → verify); scoping/queries to Sonnet
subagents. All work serial (Rule 9), cache-cleared before every Julia run,
verified by the orchestrator before each commit (Rule 2), committed + pushed
per milestone (commit-cadence).

## Headline

**4 new BIT_IDENTICAL reference tests: T52, T62, T61, T53.** Plus T54
STRUCTURAL_FAIL → DIFFS (structure complete), 3 xsdir (tape35) flips, the
thermr lat=10 read path, and a latent `natom` bug fix. ~10 deep bugs
root-caused. Zero regressions across all guards (T01/T02/T08 neutron,
T50/T62/T22 BI).

Running BIT_IDENTICAL count (tape34/various): T03, T22, T50 (Phase 77) +
**T52, T62, T61, T53** (this session). `reports/REFERENCE_SWEEP.md` is now
stale — a fresh sweep is the recommended first step next session.

## Commits (master)

1. `edaeefe` ACER ptlegc/coul (LTP=1/2) → **T52 tape34 BIT_IDENTICAL** (3986/3986).
2. `2bc4cb8` thermr lat=10 coherent-elastic MF7/MT2 read path + latent `natom` bug fix.
3. `d394cd7` ACER charged NTR>0 reaction path → **T62 tape34 BIT_IDENTICAL** (7221/7221) (5 deep bugs).
4. `ce6723a` ACER iopt=7 suffix rename + class-aware xsdir → **T61 BIT_IDENTICAL** + T50/T52/T62 tape35 0/1→1/1 BI.
5. `e796bac` ACER acelcp particle-type blocks → **T53 tape34 BIT_IDENTICAL** (12030/12030) + T54 structure complete; T51 crash-hardened.

## Per-fix detail

### T52 — ptlegc/coul (bead NJOY_jl-4ee)
`acecpe_one_incident` (ace_charged.jl) only handled LTP=12 (tabulated, T50).
Ported `ptlegc`+`coul` (acefc.f90:7980-8201) for LTP<12 (Legendre nuclear-
amplitude / residual-XS Coulomb expansion); dispatched when `sub.ltp<12`.
`read_mf6_law5` keeps raw LIST `coeffs` for LTP<12. **Root cause of initial
under-refinement:** `nt = nint(c(6))` is the LIST CONT N2 (=NL) field, NOT a
data value — threaded as `sub.nl`.

### thermr lat=10 (bead NJOY_jl-fvv)
New `read_mf7_mt2` (thermr.jl): ENDF MF7/MT2 LTHR=1 coherent-elastic structure
factors, temperature-matched (`rdelas` thermr.f90:477-592), sigcoh lat=10
transform (E_i→τ²_i, S_cum→δS_i, scon=1; thermr.f90:1150-1166). LTHR dispatch
(icoh=10·lthr) in `_thermr_sab!` + `_recompute_thermr_mf6!`, hardcoded
lat=1/2/3 fallback kept. **Latent bug:** `read_mf7_mt4` read `natom` from B(4)
(=2.277 elastic emax for Al → round=2 → /2 bound-XS error); Fortran calcem
takes natom from the input card (thermr.f90:1670-1676). T01 byte-identical
(natom no-op for H-in-H2O). T70 lat=10 elastic now present but MF3 grid
over-produces (213800 vs 112817) → bead NJOY_jl-c3q. iel/LTHR=2 → NJOY_jl-69j.

### T62 — charged NTR>0 reaction path (bead NJOY_jl-vk9)
61 value lines off (line count already exact). Five deep bugs:
1. MF3 interp law dropped for charged reactions (MT=600 INT=6 Coulomb
   penetrability: terp1 i=6 endf.f90:1634). New `_parse_mf3_tab1`/
   `extract_mf3_tab1_all` carry NBT/INT + QM/QI; acer samples native law.
2. Disappear column = Σ sigfig-7 reaction xs (mt≤200 or ≥600) (acelod
   acefc.f90:5648), not max(total−elastic,0).
3. Heating zeroed for izai≤2004 (acefc.f90:6663).
4. LQR = QI/1e6 MeV.
5. a11 e20.11 drops 'E' for |exp|≥100 (acecm.f90:791) — inert for neutron.

### T61 — iopt=7 suffix rename + xsdir (bead NJOY_jl-2h6)
iopt=7 was a verbatim copy. Apply card-2 `suff` (acer.f90:295) via `newsuff`
port (acecm.f90:619). Class-aware xsdir: thermal `i2,' 1 ',i9` (aceth.f90:2422)
vs neutron `i2,i4,1x,i8` (acefc.f90:13013). **Bonus:** T50/T52/T62 tape35 xsdir
0/1→1/1 BI (shared format fix).

### T53 + T54 — acelcp particle-type blocks (bead NJOY_jl-sta)
New `ace_lcp.jl` + `ace_lcp_build.jl`: PTYPE/NTRO/PLOCT + per-type HPD/MTRH/
TYRH/LSIGH/SIGH/LANDH/ANDH/LDLWH/DLWH/YH (acelcp acefc.f90:9121-11056 + first()
230-800), MF6 LAW=2/3/4 + MT=2 elastic recoil + LAW=4 partner back-read.
`particle_production` field on `ACENeutronTable`, appended in build_xss
(NTYPE=0 → basic tape byte-identical, the gate keeping T62/T50/T52 BI).
ptleg2/pttab2 reused from ace_charged.jl. **Two prerequisite bugs (gated to
charged path):** charged union grid (unionx acefc.f90:1541-1693, all MF3 grids
∪ MF6 anchors ∪ step-1.2; T54 125→140 pts), SIG threshold ie_start walk
(acefc.f90:5594-5628, includes the threshold zero). Unported LAWs
(1/6/7/Kalbach/phase-space/MF4) @warn+skip (Phase-77 read_mf32 pattern) →
**T51 LAW=6 no longer crashes** (was MISSING → 1197 lines). Resolved the spin=1
ptlegc point-count concern (bead NJOY_jl-fu6 — the union-grid fix). T54
residual: ~2000 words @7e-7 (pre-existing Coulomb elastic 1-ULP, T01/T34 class)
+ 1 triton recoil heating[1] word >1e-5 → bead NJOY_jl-53h.

## Beads
Filed/updated: 4ee(✓) fvv(✓) vk9(✓) 2h6(✓) sta(✓) fu6(✓); plus 3rp, aia, k2v,
69j, c3q, 6si, 53h, 5tm, fu6 (history). Closed 6, opened the accurate
follow-on backlog.

## Recommended next steps
1. **Fresh full sweep** — `reports/REFERENCE_SWEEP.md` is stale (predates this
   session's 4 BI flips + T54 structure + xsdir flips).
2. **aplots viewr plot tape (tape33, bead NJOY_jl-cnh.4)** — the ONLY remaining
   blocker to FULL pass on T50/T52/T62 (and T54 after its FP grind). 4 potential
   full-test flips from one feature. Complex viewr port.
3. T54 FP residual (53h): the 1 heating word flips tape34 → NUMERIC_PASS@1e-5.
4. acelcp unported LAWs (5tm) for T51/T14/T71; T51 value bug (6si).
5. groupr MT=251/252 + empty-skip (01v) — diversify; T15 tape91=40.
