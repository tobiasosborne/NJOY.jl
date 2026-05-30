### Phase 31: SAB kernel sabflg fix + inelastic disbar stepping — 1072→725 failures at 1e-9

**WARNING TO NEXT AGENT: READ THIS ENTIRE SECTION AND ALL TRAPS BEFORE WRITING ANY CODE.**

**2 bugs fixed, 2 dead-end investigations documented. Rule 6 saved hours.**

**Bug 1 — FIXED (SAB kernel sabflg/SCT path mismatch in sab_kernel, thermr.jl)**:

**Root cause**: When `_interp_sab` returned sabflg (-225), `sab_kernel` unconditionally fell through to the SCT (Short Collision Time) approximation. But in the Fortran `sig` function (thermr.f90:2514-2612), there are TWO paths that lead to sabflg:

1. **sabflg corner check** (lines 2556-2559): If ANY of the 4 corner cells `sab[ia,ib]`, `sab[ia+1,ib]`, `sab[ia,ib+1]`, `sab[ia+1,ib+1]` is ≤ sabflg → jump to label 170 (SCT). This ONLY fires when `a*az >= test2 (30) OR b >= test2 (30)` (line 2555).

2. **terpq returns sabflg** (line 2569): When `a*az < 30 AND b < 30`, the corner check is SKIPPED (go to label 155). Terpq runs, but may return sabflg if the SAB table values at those positions are all -225. The Fortran then computes `exp(sabflg - bb/2) ≈ exp(-215) ≈ 0`, which is below sigmin (1e-10) → sig=0.

Julia's `sab_kernel` treated BOTH cases the same → SCT fallback. But for case 2, SCT gives ~1e-5 barn (nonzero) instead of ~0. This caused Julia's `sigl_equiprobable` to integrate a nonzero kernel tail at mu > -0.453 where the Fortran correctly returns 0, producing 3.6% sigma overestimates at high-beta secondary energies.

**How it was found**: 3+1 agent pattern. Patched Fortran `sig` function with `write(*,*)` at both the terpq return (line 2574) and the SCT return (line 2611). Compared side-by-side with Julia's `sab_kernel` at E=0.56 eV, E'=0.049 eV (beta=20.2). The Fortran SCT values MATCHED Julia exactly (2.235e-4 at mu=-1). But the Fortran sigl total sigma was 4.247e-5 while Julia got 4.402e-5 (3.6% higher). The extra came from Julia's nonzero kernel at mu > -0.453 where `a*az < 30` causes the Fortran to take the terpq path (returning 0) while Julia took SCT (returning nonzero).

**Fix**: In `sab_kernel` (thermr.jl line 322), added a check: when `_interp_sab` returns sabflg AND `at*awr < 30 && bt < 30` (the skip condition was true), use the terpq result (`exp(sabflg - bb/2) → 0`) instead of the SCT fallback.

**Impact**: MF6/MT=229 diffs: **511 → 50** (90% reduction). Total 1e-9 failures: **1072 → 742** (31% reduction).

**Trap 96 (NEW — FIXED)**: The Fortran `sig` function has TWO paths to sabflg: (1) corner check → SCT, (2) terpq → sabflg value. These give DIFFERENT results. Julia must distinguish them. The key: when `a*az < 30 AND b < 30`, the sabflg came from terpq, NOT from the corner check. Use `exp(sabflg - bb/2)` (≈0), NOT SCT. The condition `a*az = alpha_lat * AWR` crosses 30 at a specific mu value, creating a sharp kernel discontinuity that both Julia and Fortran must match exactly.

---

**Bug 2 — FIXED (inelastic disbar stepping state machine in heatr.jl)**:

**Root cause**: Julia's `compute_kerma` computed discrete inelastic damage (MT=51-90) by calling `_disbar_damage_fl` directly at EVERY query energy. The Fortran's `disbar` subroutine (heatr.f90:1829-2013) uses a **stepping state machine** with `save` variables: it evaluates the 64-point GL quadrature at 1.1x-stepped bracket points (clamped to MF4 grid boundaries via `enext` from `hgtfle`), then LINEARLY INTERPOLATES to the query energy. This applies to ALL discrete MTs, not just elastic — the Fortran initializes disbar (`e=0` call) once per MT, then calls it for each energy with state persisting.

**How it was found**: 3+1 agent pattern. Agent 1 deeply read Fortran disbar, confirmed the stepping state machine applies to ALL MTs (disbar is called at line 1190 for all `icon=1` reactions). Agent 2 read `hgtfle`, confirmed it properly zero-pads (disproving the stale-fl hypothesis). Agent 3 confirmed Julia only uses the state machine for elastic. Direct comparison showed Julia elastic damage MATCHED Fortran exactly at E=10.448 MeV (both have bracket [10.0, 10.5] with identical dame values), while the TOTAL MT=444 differed by 93.8 eV-barn (0.28%), proving the error comes from inelastic MTs.

**Fix**: Generalized `build_disbar_damage_vector` with a `thresh` parameter for inelastic MTs. For `thresh > 0`: `r = 0` when `thresh >= e*(1-small)`, `r = sqrt(1-thresh/e)` above threshold (matching Fortran lines 1957-1963). In `compute_kerma`, precompute damage vectors for all discrete MTs with MF4 data using this state machine, stored in `_inel_dame_vecs`. The per-energy loop uses precomputed values (`_inel_dame_vecs[mt][ie] * sigma`) instead of direct evaluation.

**Impact**: MT=444 exact matches: **36 → 52** (+16). 1e-9 failures: **742 → 725** (-17). No regression at any tolerance.

**Trap 97 (NEW — FIXED)**: Fortran disbar uses the SAME stepping state machine (1.1x bracket + hgtfle enext clamping + linear interpolation) for ALL discrete MTs, not just elastic. Each MT gets its own disbar instance (initialized at e=0 per MT, save variables persist per-MT). Julia must precompute damage vectors for MT=51-90 the same way as MT=2.

**Trap 98 (NEW — FIXED)**: The threshold guard for inelastic r must match Fortran EXACTLY: `thresh >= e*(1-small)` → `r=0` (lines 1958-1959). Using `e < thresh*(1-small)` has different boundary behavior and causes 2.6% errors near threshold. Also: when `r=0`, the recoil energy `e2 = e*(1+g^2)*afact/arat = 0` for all mu, so `dame = 0` naturally — no separate threshold skip needed.

---

**Dead-end investigation 1 — Fortran hgtfle stale fl coefficients (WRONG HYPOTHESIS)**:

Investigated whether Fortran `hgtfle` retains stale Legendre coefficients from previous bracket evaluations when NL decreases. Diagnostic showed fl[6]=-0.002168 at E=3.0 MeV (which has NL=4), matching the E=2.98 MeV value (NL=6). Hypothesis: hgtfle doesn't zero-pad shorter records.

**DISPROVED by Agent 2**: hgtfle (heatr.f90 lines 4218-4406) at label 130 (lines 4363-4381) explicitly does:
```fortran
nlmax=max(nlo,nhi)
do i=1,nle
   if (i.le.nlmax) then
      call terp1(elo,flo(i),ehi,fhi(i),e,fle(i),innt)
   else
      fle(i)=0  ! ZERO-PAD above nlmax
   endif
enddo
nle=nlmax
```

**Rule 6 confirmed**: The stale fl values seen in the diagnostic came from the fl(65) array in disbar itself (save variable), NOT from hgtfle's output. The hgtfle function properly zero-pads. The stale values in disbar's fl array are overwritten by hgtfle's output. The diagnostic printed BEFORE hgtfle was called, showing the pre-call state.

Implementing the stale-fl hypothesis introduced marginal regressions (+10 at 1e-4, +12 at 1e-3) and was reverted.

**Dead-end investigation 2 — Direct inelastic disbar without threshold (WRONG THRESHOLD)**:

First attempt at the inelastic state machine used `e < thresh * (1.0 - small)` as the threshold guard. This caused 2.6% errors at E=4.98 MeV (near MT=51 threshold at 4.81 MeV). The bracket [el, en] straddled the threshold, and the wrong guard allowed nonzero damage at el < thresh. Fixed by using the Fortran's exact condition: `thresh >= e * (1.0 - small)`.

---

**T01 results after Phase 31:**

```
Structural: 41/41 sections, 32962/32962 lines — EXACT MATCH
rel_tol=1e-9: 725 fail (was 1072 at start of session)
rel_tol=1e-7: 725 fail (was 1072)
rel_tol=1e-4: 191 fail (was 226)
rel_tol=1e-3:  74 fail (was 98)
Total data exact: 1658/2386 (69.5%) — was 1641 (68.8%)
```

**Per-section diff breakdown (828 total line diffs, down from 1290):**

| Section | Diffs | Class | Fixability |
|---------|-------|-------|------------|
| MT=444 | 271 | Damage energy | 0.28% worst. Remaining: evaporation (MT=91) adaptive convergence + elastic stale-fl (minor) |
| MT=301 | 242 | KERMA | 0.27% worst. Cascades from MT=1 sigma1 ULP |
| MF3/MT=229 | 161 | calcem XS | 0.002% worst. terp_lagrange interpolation |
| MT=1 | 67 | sigma1 ULP | IRREDUCIBLE. Broadening FP accumulation order |
| MF6/MT=229 | 50 | sigl FP | Mostly ±1 ULP at high IEs |
| MF6/MT=221 | 28 | Free gas FP | ±1 ULP |
| MF1/MT=451 | 5 | Directory NC | Mechanical fix: `div(N,6)` vs `ceil(N/6)` |
| MT=221 | 2 | emax boundary | Minor |
| MT=2 | 1 | sigma1 ULP | IRREDUCIBLE |
| MF2/MT=151 | 1 | EH value | Trivial |

**Session improvement: -347 tolerance failures at 1e-9 (32% reduction), -24 at 1e-3 (24% reduction).**

**Files changed**:
- `src/processing/thermr.jl` — sabflg/SCT path fix in `sab_kernel` (Trap 96)
- `src/processing/heatr.jl` — inelastic disbar stepping in `build_disbar_damage_vector` + `compute_kerma` precompute (Traps 97-98)
- `test/validation/t01_pipeline.jl` — path fix (old `/home/tobias/` → `/home/tobiasosborne/`)

**How to verify Phase 31 fixes:**
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/t01_pipeline.jl
# Expected: 41/41 sections, 32962 lines, 1658/2386 data exact
# Then run tolerance test (see Python script in HANDOFF section above)
# Expected: 725 at 1e-9, 191 at 1e-4, 74 at 1e-3
```

**Trap 99 (CORRECTED)**: `nld` in Fortran disbar is a LOCAL variable (line 1840), NOT in the save statement. BUT the stale-fl hypothesis for `nld` is different from the stale-fl issue in `hgtfle`'s `flo` array. See Phase 32 Trap 102 for the REAL stale-fl bug in hgtfle.

**Trap 100 (CRITICAL WARNING)**: The `test/validation/t01_pipeline.jl` script has HARDCODED absolute paths. These were updated from `/home/tobias/` to `/home/tobiasosborne/` in this session. If running on a different machine, ALL 5 path references must be updated. Consider using `@__DIR__` or relative paths in a future cleanup.

**Trap 101 (CRITICAL WARNING)**: When implementing the inelastic disbar state machine, the threshold condition MUST be `thresh >= e * (1.0 - small)` (matching Fortran line 1958), NOT `e < thresh * (1.0 - small)`. These are NOT equivalent near the boundary due to floating-point rounding. The wrong condition caused 2.6% errors at MT=51 threshold (E=4.98 MeV). Always read the Fortran (Rule 3) and match the exact comparison direction.

---

### Phase 32: HEATR deep grind — 5 bugs fixed, 725→652 at 1e-9

**WARNING TO NEXT AGENT: THIS SESSION CORRECTED MULTIPLE PREVIOUS WRONG CLAIMS. READ EVERY TRAP CAREFULLY. RULE 6 IS NOT OPTIONAL.**

**Previous claim "stale-fl hypothesis definitively dead" (Phase 31) was WRONG.** Phase 31 investigated whether `nld` (number of Legendre coefficients) was stale — correctly concluded it is NOT stale (local variable). But the REAL stale-fl issue is in `hgtfle`'s `flo` SAVE array: the SLIDE at label 120 (line 4333) only copies `nhi` values from `fhi` to `flo`, leaving higher-order `flo` coefficients STALE from previous brackets. This is a DIFFERENT mechanism than the `nld` stale hypothesis. Confirmed via gdb.

**Previous claim "MT=301 cascades from sigma1 ULP" (Phase 31) was WRONG.** The sigma1 ULP is ~2e-6 relative, but MT=301 worst error was 0.27%. The REAL root cause was the discrete inelastic ebar computation using direct evaluation instead of Fortran's stepped cn interpolation (see Bug 4 below).

**5 bugs found and fixed (all via 3+1 agent pattern + Fortran gdb diagnostics):**

1. **MF2/MT=151 EH value (FIXED — Trap 102)**: Julia wrote EH from the ENDF range data (e.g. 1e5 for C-nat LRU=0). Fortran reconr recout uses `eresr`: resolved range EH for materials with LRU=1, `ehigh=2e7` for LRU=0-only. Fix: check for resolved ranges, use `ehigh=2e7` when all ranges are LRU=0. Verified across T01/T02/T08/T34/T45 oracles.

2. **MF1/MT=451 directory NC values (FIXED — Traps 103-104)**: Five NC values differed. Two root causes:
   - **MF3 thermr sections**: Fortran tpend NC formula is `3 + (ne+2)/3` (Fortran integer division) where `ne` is the coh/elastic grid count BEFORE tpend adds sentinel points. For T01 thermr2: ne=569 (from coh), NP=571 (after sentinels). `(569+2)/3 = 571/3 = 190` in Fortran → NC=193. Julia was computing `3 + cld(571, 3) = 194`. Fix: pass `thermr_coh_ne` parameter to `write_full_pendf`.
   - **MF6 sections**: Fortran `ncds` undercounts TAB1 yield by 1 (assumes packed interp+data on 1 line, but output has 2 separate lines). Fix: subtract 1 from MF6 NC.
   - **MF6 stubs**: Fortran `ncdse=3` (same undercount). Fix: hardcode NC=3.

3. **hgtfle stale-fl in damage GL integration (FIXED — Trap 105)**: The Fortran `hgtfle` subroutine (heatr.f90 lines 4223-4411) maintains SAVE arrays `flo(65)` and `fhi(65)`. The SLIDE at label 120 (line 4333): `do i=1,nhi; flo(i)=fhi(i); enddo` only copies `nhi` values. When NL decreases at a grid boundary (e.g. NL=7→5 at E=3.0 MeV for C-12), `flo[6]` retains its stale value from the previous bracket. This stale coefficient participates in the subsequent GL integration via `nlmax = max(nlo, nhi)`. C-12 MF4/MT=2 has 4 NL decrease transitions (E=2.06, 3.0, 5.0, 6.8 MeV). Julia's `_interp_legendre` correctly zero-pads, giving fl[6]=0 instead of the stale value. Fix: implemented persistent hgtfle state via closure in `build_disbar_damage_vector`, matching the Fortran SLIDE partial-copy behavior.

   **Verified via gdb**: Patched Fortran heatr.f90 disbar with `write(*,*)` at the GL integration output. At E=3.0 MeV: Fortran fl[6]=-0.00217 (stale from 2.98 MeV bracket), Julia now matches.

   **Impact**: MT=444 exact matches +4 (237→233), 1e-9 failures −2.

   **CRITICAL: Phase 31's claim "stale-fl hypothesis definitively dead" was about a DIFFERENT variable (nld) than what was actually stale (flo array). Rule 6 applies: be skeptical of previous analysis.**

4. **Non-broadened MF3 sigma interpolation for heatr (FIXED — Trap 106)**: The Julia pipeline interpolated raw ENDF MF3 TAB1 data for non-broadened MTs (51-91, 107). The Fortran heatr reads from the PENDF (reconr output), which has threshold-adjusted energies (`thrxx = sigfig(thrx, 7, +1)`) and pseudo-threshold zeros. For MT=51 at threshold E=4.812123 MeV: raw MF3 interpolation gave sigma=2.59e-5, reconr PENDF has sigma=0.0 (threshold shift creates a zero at the adjusted threshold). Fix: use `_get_legacy_section` (reconr-processed data) for non-broadened MT cross sections.

   **Found via 3+1 agent pattern**: Agent 1 patched Fortran nheat, found sigma=0 at threshold. Agent 2 ran Julia diagnostic, found sigma=2.59e-5. Agent 3 compared side-by-side.

   **Impact**: MT=444 exact matches +2 (233→233), 1e-9 failures −2.

5. **Discrete inelastic ebar: stepped cn interpolation (FIXED — Trap 107, MAJOR)**: Julia's `compute_kerma` computed `discrete_inelastic_ebar(E, Q, A; mu_bar)` DIRECTLY at every broadened grid energy. The Fortran `disbar` computes `cn = (1 + 2*b*wbar + b²)*afact` at 1.1x-stepped bracket endpoints (using hgtfle MF4 interpolation) and LINEARLY INTERPOLATES cn to each query energy via `terp1`. The cn fraction varies with energy through `r = sqrt(1-thresh/E)`, so direct evaluation at each energy differs from the bracket interpolation by up to 0.27%.

   **Found via 3+1 agent pattern**: Agent 1 patched Fortran nheat to print per-MT heating at E=11.5 MeV. Agent 2 ran Julia diagnostic. Agent 3 read the Julia code. Comparison revealed MT=51 ebar diff of 7400 eV (0.12%), MT=53 ebar diff of 19000 eV (1.9%), all from the different ebar computation methods.

   Fix: `build_disbar_damage_vector` now returns BOTH `dame_out` and `ebar_out` vectors. The cn computation uses the same hgtfle stale-fl state and bracket stepping as damage. `compute_kerma` uses `_inel_ebar_vecs[mt][ie]` for discrete inelastic MTs (51-90) instead of direct `discrete_inelastic_ebar`.

   **Impact**: MT=301 worst error halved (2.71e-3 → 1.33e-3). MT=301 exact matches doubled (65 → 130 out of 307). 1e-9 failures −63. **This was the single largest fix.**

**T01 results after Phase 32:**

```
Structural: 41/41 sections, 32962/32962 lines — EXACT MATCH
rel_tol=1e-9: 652 fail (was 725 at start, was 1072 at Phase 31 start)
rel_tol=1e-7: 652 fail
rel_tol=1e-5: 235 fail (was 311)
rel_tol=1e-4: 148 fail (was 191)
rel_tol=1e-3:  52 fail (was 74)
Total data exact: 1728/2386 (72.4%) — was 1658 (69.5%)
```

**Per-section diff breakdown (657 total diffs, down from 828):**

| Section | Diffs | Worst | Change | Root cause |
|---------|-------|-------|--------|------------|
| MT=444 | 233 | 2.77e-3 | −4 | Stale-fl fix helped at NL transitions |
| MT=301 | 164 | 1.33e-3 | −63 | Stepped ebar fixed the ebar computation |
| MF3/MT=229 | 159 | 5.27e-4 | 0 | calcem XS interpolation (unchanged) |
| MT=1 | 46 | 2.31e-6 | 0 | sigma1 ULP (IRREDUCIBLE) |
| MF6/MT=229 | 36 | 2.11e+0 | 0 | 2 near-zero kernel + ±1 ULP cosines |
| MF6/MT=221 | 17 | 7.71e-7 | 0 | Free gas cosine ULP |
| MT=221 | 1 | 9.17e-6 | 0 | emax boundary |
| MT=2 | 1 | 4.81e-7 | 0 | sigma1 ULP (IRREDUCIBLE) |

**Files changed**:
- `src/processing/pendf_writer.jl` — MF2 EH, MF1 NC directory (Traps 102-104)
- `src/processing/heatr.jl` — hgtfle stale-fl, stepped ebar (Traps 105, 107)
- `test/validation/t01_pipeline.jl` — reconr MF3 for heatr, thermr_coh_ne param (Trap 106)

**How to verify Phase 32 fixes:**
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/t01_pipeline.jl
# Expected: 41/41 sections, 32962 lines, 1728/2386 data exact
# Then run tolerance test:
# Expected: 652 at 1e-9, 148 at 1e-4, 52 at 1e-3
```

**Trap 102 (NEW — FIXED)**: MF2/MT=151 EH value: Fortran reconr recout uses `eresr` — resolved range EH for materials with LRU=1 resonances, `ehigh=2e7` for LRU=0-only materials. Julia was using `maximum(r.EH for r in iso.ranges)` which gives the ENDF range EH (can be < 2e7 for LRU=0). Verified: T01 (LRU=0) needs EH=2e7; T02 (LRU=1+2, resolved EH=200) needs EH=200; T08 (LRU=1, EH=70000) needs EH=70000.

**Trap 103 (NEW — FIXED)**: Fortran tpend MF3 NC formula is `3 + (ne+2)/3` where `ne` is the value passed to tpend (the coh/elastic grid count BEFORE tpend adds sentinel points). For free gas (ne=NP): NC = `3 + (NP+2)/3`. For coh grid (ne=NP−2): NC = `3 + NP/3`. Fortran integer division means NC is 1 less than actual line count when NP is not divisible by 3 AND sentinels were added. This is a cosmetic bug in the Fortran tpend that must be reproduced.

**Trap 104 (NEW — FIXED)**: Fortran `ncds` undercounts MF6 TAB1 yield by 1. The `tab1io` with nw=12 (NR=1, NP=2) writes 3 lines (CONT + interp + data), but ncds adds only 2 (line 1950). The output has 3 lines, but the directory NC says 2 fewer. For MF6 stubs, `ncdse=3` (same undercount: actual 4, NC says 3). Must reproduce to match.

**Trap 105 (NEW — FIXED, CORRECTS Phase 31)**: The hgtfle stale-fl issue IS REAL — Phase 31 was wrong to declare it "definitively dead." The Phase 31 investigation checked whether `nld` (Legendre order count) was stale — it is NOT (local variable). But the `flo(65)` ARRAY in hgtfle IS stale: the SLIDE at label 120 only copies `nhi` values, leaving higher indices untouched. For C-12 at E=3.0 MeV (NL=7→5 transition): flo[6] retains the 2.98 MeV value (-0.00217) instead of 0. This causes 0.016% error in elastic damage at that energy. The fix uses a persistent closure state in `build_disbar_damage_vector` matching the Fortran SAVE arrays.

**Trap 106 (NEW — FIXED)**: Non-broadened MT cross sections for heatr MUST use the reconr-processed data (`_get_legacy_section`), NOT raw ENDF MF3 interpolation. The reconr output has threshold-adjusted energies (`thrxx = sigfig(thrx, 7, +1)`) and pseudo-threshold zeros. At MT=51 threshold (4.812123 MeV): reconr PENDF has sigma=0, raw MF3 gives sigma=2.59e-5. The Fortran heatr reads from the PENDF.

**Trap 107 (NEW — FIXED, MAJOR)**: The Fortran `disbar` computes ebar via the cn fraction at 1.1x-stepped bracket endpoints, then INTERPOLATES cn to query energies via `terp1`. Julia was computing `discrete_inelastic_ebar` DIRECTLY at every energy. The cn fraction `(1 + 2*b*wbar + b²)*afact` varies with energy through `r = sqrt(1-thresh/E)`, so direct evaluation differs from bracket interpolation by up to 0.27%. Fix: `build_disbar_damage_vector` now returns both `dame_out` and `ebar_out`, and `compute_kerma` uses the stepped ebar for discrete inelastic MTs.

**Trap 108 (CRITICAL WARNING)**: Previous agents' diagnostic scripts often used APPROXIMATIONS (e.g., `wbar = 0.0  # approximate`) that DON'T match the actual pipeline computation. ALWAYS verify that diagnostic results match the pipeline output BEFORE drawing conclusions. The Phase 32 Julia agent computed h=3.623e5 for MT=53 at 11.5 MeV with wbar=0, but the actual pipeline (with MF4 mu_bar) gives a different value. This led to initially wrong ebar diff estimates.

---

### Phase 33: HEATR bracket stepping + below-threshold skip — 4 bugs, 652→578 at 1e-9, 52→2 at 1e-3

**4 bugs found and fixed, all via the 3+1 agent pattern + Fortran gdb diagnostics.**

**Bug 1 — FIXED (MT=91 conbar damage stepping — Trap 109, MAJOR)**: Julia's `compute_kerma` called `evaporation_damage(E, ...)` DIRECTLY at every broadened grid energy. The Fortran `conbar` (heatr.f90:2273-2290) evaluates `anadam` at 1.5x-stepped bracket endpoints and LINEARLY INTERPOLATES damage to query energies. The ebar is computed directly via `anabar` at the query energy (line 2267) — NOT stepped.

   **Found via 3+1 agents**: Agent 1 patched Fortran nheat to print per-MT dame at E=10.448 MeV. Agent 2 ran Julia diagnostic. Comparison: MT=91 Julia dame=36393 vs Fortran dame=23982 (+51.8%), contributing +92.8 eV to MT=444 — the dominant error source.

   Fix: Added `build_conbar_damage_vector` with step=1.5 bracket stepping matching conbar. KEY: only damage is stepped; ebar uses direct `evaporation_ebar` (Fortran anabar computes analytically at each energy, no bracket interpolation).

   **Impact**: MT=444 worst: 0.277% → 0.049%. 1e-3 failures: 52 → 22.

**Bug 2 — FIXED (MT=54-62 isotropic disbar ebar — Trap 110)**: For MTs WITHOUT MF4 data, Fortran disbar still runs the 1.1x bracket stepping with isotropic coefficients (fl(1)=1, fl(2)=0, imiss=1, enext=etop). Julia was falling through to direct `discrete_inelastic_ebar` evaluation. At E=16.55 MeV: MT=54-62 ebar diffs of 6k-35k eV each, totaling ~5000 eV heating difference.

   Fix: Generate synthetic isotropic MF4 data `([1e-5, 2e7], [[1.0, 0.0], [1.0, 0.0]])` for all MTs 51-90 that lack MF4. Pass to `build_disbar_damage_vector` for both damage and ebar stepping.

   **Impact**: 1e-3: 22 → 6. 1e-4: 91 → 88.

**Bug 3 — FIXED (MT=107 capdam stepping — Trap 111)**: Fortran `capdam` (heatr.f90:1792-1826) uses 1.1x bracket stepping with SAVE variables for charged-particle damage (MT=103-107), same pattern as disbar. Julia's `capdam_particle` evaluated directly. At E=9.56 MeV: MT=107 Julia dame=20033 vs Fortran dame=19964 (+68.8 eV).

   Fix: Added `build_capdam_damage_vector` with step=1.1 matching capdam state machine.

   **Impact**: MT=444 worst: 0.049% → 0.037%. 1e-9: 652 → 651.

**Bug 4 — FIXED (Below-threshold bracket skip — Trap 112, MAJOR)**: Fortran `disbar` is only called at energies where sigma > 0 (above threshold). Julia's `build_disbar_damage_vector` iterated over ALL energies including below-threshold ones, building 1.1x-stepped brackets that STRADDLED the threshold. When a bracket endpoint fell above threshold (r>0, cn >> afact) but the interpolation target was below threshold, the linearly-interpolated cn was wrong by up to 4x.

   **Found via 3+1 agents**: Traced Fortran disbar bracket state for MT=63 at E=16.55 MeV. Fortran bracket: [0, 15990000] with el=0 from initialization (disbar never called below thresh). Julia bracket: [15400000, 16940000] straddling thresh=15990000, with corrupted cn at right endpoint.

   Fix: Skip energies below threshold in the stepping loop: `if thresh > 0 && ee < thresh * (1.0 - small); continue; end`

   **Impact**: MT=301 worst: 0.114% → 1.0e-6 (1000x improvement). MT=301 exact: 130/307 → 198/307. 1e-9: 651 → 578 (-73). 1e-3: 6 → 2.

**MF6/MT=229 cosine investigation (2 remaining 1e-3 failures — DEEPLY INVESTIGATED, NOT FIXED)**:

   The 2 remaining failures are equi-probable cosines at E=0.01820 eV, E'=0.59286167 eV where sigma=3.79719e-10 (physically zero). The kernel values at this E/E' straddle the sigmin=1e-10 cutoff:

   - Fortran kernel at mu=0.995: sig=1.006e-10 (0.6% above sigmin → NONZERO after cutoff)
   - Julia kernel at mu=0.995: sig=9.939e-11 (0.6% below sigmin → ZERO after cutoff)

   Root cause: Julia's `_terpq` quadratic interpolation returns log S(α,β) = -18.1958 at alpha=1.343, beta=22.714. Fortran's terpq returns -18.1923. Diff = 0.0035 in log space. This shifts exp(s - bb/2) by ~0.35%, which at the sigmin cliff edge flips the kernel from above to below. The CDF then has one fewer nonzero point, completely changing the equi-probable bin boundaries.

   Verified: terpq formulas are IDENTICAL. Stencil indices (ia=5, ib=73) are IDENTICAL. The 0.0035 difference is from IEEE 754 intermediate precision — the Fortran compiler may use 80-bit x87 extended precision registers while Julia uses 64-bit SSE. Two-step computation `b=...; b=b-...` vs single-expression `bq = ... - ...` was tested with no effect.

   **This is genuinely at the FP precision floor.** The physical impact is zero (sigma=3.8e-10). Both codes produce the same total sigma. Only the cosine distribution differs.

**Trap 109 (NEW — FIXED)**: Fortran conbar evaluates anadam damage at 1.5x-stepped bracket endpoints (step=1.5 at line 2077), then linearly interpolates via terp1 (line 2290). Ebar is computed directly via anabar at the query energy (line 2267) — NOT bracket-stepped. Julia was computing both directly. The damage stepping accounts for ~93 eV systematic excess at E=10.4 MeV.

**Trap 110 (NEW — FIXED)**: MTs without MF4 data still need disbar bracket stepping. Fortran disbar with imiss=1: fl(1)=1, fl(2)=0, enext=etop. Julia was falling through to direct evaluation. For isotropic MTs, generate synthetic MF4 data covering [1e-5, 2e7] with coefficients [1.0, 0.0].

**Trap 111 (NEW — FIXED)**: Fortran capdam uses 1.1x bracket stepping with SAVE variables (en, damn, el, daml) for charged-particle damage (MT=103-107). Same pattern as disbar. Julia's capdam_particle evaluated directly.

**Trap 112 (NEW — FIXED, MAJOR)**: Fortran disbar is ONLY CALLED at energies where sigma > 0 (above threshold). The SAVE variables stay at initialization values (en=0, cn=0) until the first above-threshold call. Julia's build_disbar_damage_vector iterated over all energies, advancing the bracket through below-threshold points. When the bracket straddled the threshold, the right-endpoint cn (computed above threshold) corrupted the interpolation for below-threshold query energies. Fix: skip below-threshold energies entirely.

**Trap 113 (INVESTIGATED, NOT FIXABLE)**: The 2 remaining MF6/MT=229 failures at 1e-3 are from terpq FP precision at the sigmin=1e-10 boundary. At E=0.01820 eV, E'=0.5929 eV: Julia's terpq returns -18.1958 vs Fortran's -18.1923 (diff=0.0035 in log S). This flips the kernel from 9.94e-11 (below sigmin) to 1.01e-10 (above) at mu=0.995, changing the equi-probable cosine CDF. The formulas and stencil indices are identical — the difference is from IEEE 754 intermediate precision (compiler-level). Physical impact: zero (sigma=3.8e-10).

**T01 results after Phase 33:**

```
Structural: 41/41 sections, 32962/32962 lines — EXACT MATCH
rel_tol=1e-9: 578 fail (was 652 at start of session)
rel_tol=1e-7: 578 fail
rel_tol=1e-5: 144 fail (was 235)
rel_tol=1e-4:  47 fail (was 148)
rel_tol=1e-3:   2 fail (was 52) — only irreducible MF6/MT=229 cosines
Total data exact: 1799/2386 (75.4%) — was 1728 (72.4%)
```

**Per-section diff breakdown (578 total at 1e-7):**

| Section | Fails | Worst | Change from Phase 32 | Root cause |
|---------|-------|-------|---------------------|------------|
| MT=444 | 227 | 3.69e-4 | −6 | Bracket stepping residual + ±1 ULP |
| MF3/MT=229 | 159 | 5.27e-4 | 0 | calcem XS interpolation (unchanged) |
| MT=301 | 96 | 9.99e-7 | −68 | **Now sub-ULP!** Only sigma1 cascade |
| MT=1 | 46 | 2.31e-6 | 0 | sigma1 broadening ULP |
| MF6/MT=229 | 32 | 3.39e-1 | −4 | sigl FP + 2 near-zero kernel |
| MF6/MT=221 | 16 | 7.71e-7 | −1 | Free gas cosine ULP |
| MT=2 | 1 | 4.81e-7 | 0 | sigma1 ±1 ULP |
| MT=221 | 1 | 9.17e-6 | 0 | emax boundary |

**Files changed**:
- `src/processing/heatr.jl` — build_conbar_damage_vector, build_capdam_damage_vector, compute_kerma precomputation + MT=91/103-107 wiring, below-threshold skip in build_disbar_damage_vector
- (thermr.jl terpq two-step change tested and reverted — no effect)

**How to verify Phase 33 fixes:**
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/t01_pipeline.jl
# Expected: 41/41 sections, 32962 lines, 2017/2386 data exact
# Then run tolerance test (see Python script in HANDOFF):
# Expected: 355 at 1e-9, 0 at 1e-5, 0 at 1e-4, 0 at 1e-3
```

---

### Phase 34: T01 passes at 1e-5 — 7 bugs fixed, 578→355 at 1e-9, 47→0 at 1e-4, 2→0 at 1e-3

**7 bugs found and fixed, all via the 3+1 agent pattern + Fortran gdb diagnostics.**

**Bug 1 — FIXED (sigl_equiprobable phase 2 bin-finding inverted condition — Trap 114, MAJOR)**:

**Root cause**: In `sigl_equiprobable` phase 2 CDF inversion (thermr.jl line 863), when the kernel value `yl` at the panel left edge was below `sigmin` (e.g., yl=0 at a zero-kernel-to-nonzero transition), Julia used the LINEAR formula `xn = xl + (fract-cum_sum)/max(yl,sigmin)` which divides by sigmin=1e-32 → xn≈1e21 (garbage). The Fortran sigl (label 170→175) sends `yl < sigmin` to the QUADRATIC root-finding path.

**How it was found**: 3+1 agent pattern. Patched Fortran `sig` function to print kernel values at E=0.0182, E'=0.5929 (the failing MF6/MT=229 point). Found ALL kernel values within 7% of sigmin=1e-10. Then patched Fortran `sigl` to print phase 1 total sigma and phase 2 bin boundaries. Compared side-by-side with Julia — kernel values matched to 1e-12 at same mu points, total sigma matched, but bin boundaries diverged. Traced the CDF inversion code and found the inverted condition.

**Fix**: Changed condition from `if yl < sigmin || ...` (→LINEAR) to `if yl >= sigmin && ...` (→LINEAR), sending `yl < sigmin` to the ELSE (QUADRATIC) path.

**Impact**: 1e-3 failures 2→0. 1e-9 failures 578→576.

**Trap 114 (NEW — FIXED, CORRECTS Trap 113)**: The "2 irreducible MF6/MT=229 cosine failures" from Phase 33 were NOT irreducible. They were caused by an inverted condition in the CDF bin-finding code. When yl=0 (zero kernel at panel left edge), the Fortran uses quadratic root-finding (label 175) while Julia used the linear formula. The linear formula divides by yl≈0 → garbage xn → corrupted cosine bins. **EVERY "irreducible" label so far has been wrong. NOTHING IS IRREDUCIBLE.**

---

**Bug 2 — FIXED (MT=229 one-step vs two-step interpolation — Trap 115, MAJOR)**:

**Root cause**: Julia interpolated calcem XS (94 points) directly onto the 571-point thermal grid via order-5 Lagrangian. The Fortran does TWO steps: (1) calcem label 610 interpolates xsi from 94 calcem points onto the broadened elastic grid (~145 points) via `terp(esi,xsi,94,enow,5)`, (2) coh label 190 interpolates from the ~145-point intermediate grid onto new Bragg edge energies via `terp(x,z,5,ej,4)` (5-point sliding window, ORDER 4). The one-step approach uses distant calcem stencil points; the two-step uses closer intermediate grid points.

**How it was found**: Patched Fortran calcem to print esi/xsi (94 points) — confirmed IDENTICAL to Julia's calcem output. Patched Fortran terp to print all 145 intermediate interpolation calls — confirmed only 145 calls, not 571. Realized the merge to 571 points happens in coh, not calcem.

**Fix**: Implemented two-step interpolation: (1) terp_lagrange(esi, xsi, broadened_energy, 5) for each broadened energy → intermediate grid, (2) streaming coh window interpolation from intermediate → thermal grid.

**Trap 115 (NEW — FIXED)**: Fortran calcem label 610 reads from iold (broadened elastic grid, ~145 points below emax), calls `terp(esi,xsi,94,enow,5)` at each. Then coh merges Bragg edges using a 5-point sliding window. Julia must match this two-step process, NOT interpolate directly from the 94-point calcem grid.

---

**Bug 3 — FIXED (coh interpolation order nlt1=nlt-1=4 — Trap 116)**:

**Root cause**: Fortran coh (line 783) sets `nlt1=nlt-1`. With nlt=5 (window size), nlt1=4 (interpolation order). The coh call `terp(x,z,nlt,ej,nlt1)` does ORDER 4 Lagrangian (4 of 5 window points), not order 5. Julia was using order 5.

**How it was found**: Patched Fortran coh label 190 to print nlt, nlt1, window contents, and interpolation result at E≈0.09381. Output showed `nlt=5, nlt1=4`. Computed exact 5-point Lagrangian in Julia — result didn't match. Computed 4-point — matched.

**Fix**: Changed terp_lagrange order from 5 to 4 for the coh interpolation step.

**Trap 116 (NEW — FIXED)**: Fortran coh uses `nlt1=nlt-1=4` for the interpolation order, NOT nlt=5. Near the grid end, nlt and nlt1 decrease further (lines 864-865: `if (iex.eq.ne) nlt=nlt-1; nlt1=nlt1-1`). The Julia coh interpolation must use a streaming window that tracks nlt/nlt1 reductions.

---

**Bug 4 — FIXED (coh streaming window with boundary order reduction — Trap 117)**:

**Root cause**: When the coh sliding window reaches the last input grid point (iex=ne), Fortran decreases nlt and nlt1 by 1. At E≈1.09 eV (near emax=1.2), this reduces the interpolation from order 4 to order 3, using only 4 points instead of 5. Julia's terp_lagrange on the full intermediate grid always used order 4, giving different stencils at boundaries.

**Fix**: Implemented `coh_interp_streaming` function that maintains a sliding window matching coh's exact advancement logic (label 170) and order-reduction behavior.

**Trap 117 (NEW — FIXED)**: At the grid boundary, Fortran coh's window uses FEWER points and LOWER order. The 5-point window with order 4 becomes a 4-point window with order 3 when iex reaches ne. The interpolation at high energies near emax is systematically affected.

---

**Bug 5 — FIXED (last intermediate point XS=0 — Trap 118)**:

**Root cause**: Fortran calcem label 610 (line 2460): `if (ie.eq.ne) xs=0` — the last intermediate grid point gets XS forced to zero. Julia didn't do this, so the last point had a nonzero interpolated calcem XS.

**Fix**: `intermediate_xsi[end] = 0.0` after computing the intermediate grid.

**Trap 118 (NEW — FIXED)**: Fortran calcem label 610 forces xs=0 at the last grid point (ie==ne). This creates a zero-XS sentinel at emax that affects the coh window interpolation near the grid boundary.

---

**Bug 6 — FIXED (strict emax cutoff for MT=229)**:

**Root cause**: At E=sigfig(emax,7,+1)=1.200001, the Fortran writes XS=0 for MT=229. Julia's streaming window interpolated a nonzero value because the check `e > emax*(1+small)` was too loose (1.200001 < 1.200036).

**Fix**: Changed to strict `e > emax_thermr` cutoff.

---

**Bug 7 — FIXED (capdam below-threshold bracket corruption — Trap 119, MAJOR)**:

**Root cause**: Julia's `build_capdam_damage_vector` advanced the 1.1x bracket stepping at EVERY broadened grid energy, including those below the MT=107 threshold (6.18 MeV). The Fortran capdam is only called at energies where sigma > 0 (nheat line 1374: `if (e.lt.thresh) go to 290`). This caused Julia to make ~730 bracket advances (from E=1e-5 to 9.36e6) while Fortran made ~30 (from E=6.2e6 to 9.36e6). The geometric 1.1x sequences diverged completely: Julia bracket at 9.36 MeV was [9.35, 10.29] while Fortran was [9.04, 9.94]. The different brackets produced 0.04% systematic damage errors.

**How it was found**: Patched Fortran nheat to print per-MT damage at E=9.36 MeV. Comparison showed MT=107 as the dominant error source (12.8 eV-barn, 0.21%). Patched Fortran capdam to print bracket endpoints — found [9.039, 9.943]. Julia diagnostic showed [9.351, 10.286]. Traced the bracket evolution from initialization: Julia started at E=1e-5 (first broadened energy), Fortran started at E=6.2e6 (first energy above threshold). The below-threshold energies corrupted Julia's bracket sequence.

**Fix**: Added threshold guard in `build_capdam_damage_vector`: `if thresh_cap > 0 && ee < thresh_cap * (1-small); continue; end`.

**Impact**: 1e-4 failures 8→0. 1e-5 failures 40→0. 1e-9 failures 418→355. Data exact 82.1%→84.5%.

**Trap 119 (NEW — FIXED, MAJOR)**: Fortran capdam is ONLY CALLED at energies where sigma > 0 (above threshold). The SAVE bracket state starts from the first above-threshold energy. Julia's build_capdam_damage_vector iterated over ALL broadened energies, advancing the bracket at below-threshold energies. This corrupted the entire geometric bracket sequence. Fix: skip energies below threshold in the bracket loop. **This is the same class of bug as Trap 112 (disbar below-threshold bracket skip) but for capdam instead of disbar.**

---

**T01 results after Phase 34:**

```
Structural: 41/41 sections, 32962/32962 lines — EXACT MATCH
rel_tol=1e-9: 355 fail (was 578 at start of session)
rel_tol=1e-7: 355 fail
rel_tol=1e-5:   0 fail (was 144) ← PASSES ✓
rel_tol=1e-4:   0 fail (was 47)  ← PASSES ✓
rel_tol=1e-3:   0 fail (was 2)   ← PASSES ✓
Total data exact: 2017/2386 (84.5%) — was 1799 (75.4%)
```

**Per-section diff breakdown (355 total at 1e-9):**

| Section | Fails | Worst | Root cause |
|---------|-------|-------|------------|
| MT=444 | 164 | 4.4e-6 | ±1 ULP cascade through damage×sigma |
| MT=301 | 96 | 1.0e-6 | Cascades from MT=1 sigma1 ULP |
| MT=1 | 46 | 2.3e-6 | sigma1 broadening FP accumulation |
| MF6/MT=229 | 30 | 6.0e-7 | sigl FP accumulation at high IEs |
| MF6/MT=221 | 16 | 7.7e-7 | Free gas cosine FP |
| MF3/MT=229 | 1 | 4.2e-6 | coh window at grid edge |
| MF3/MT=221 | 1 | 9.2e-6 | emax boundary |
| MF3/MT=2 | 1 | 4.8e-7 | sigma1 ±1 ULP |

**Files changed**:
- `src/processing/thermr.jl` — sigl_equiprobable phase 2 bin-finding condition fix (Bug 1)
- `src/processing/heatr.jl` — capdam below-threshold bracket guard (Bug 7)
- `test/validation/t01_pipeline.jl` — two-step interpolation, coh streaming window, order 4, last-point zero, emax cutoff (Bugs 2-6)

**How to verify Phase 34 fixes:**
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/t01_pipeline.jl
# Expected: 41/41 sections, 32962 lines, 2017/2386 data exact
# Then run tolerance test:
# Expected: 355 at 1e-9, 0 at 1e-5, 0 at 1e-4, 0 at 1e-3
```

---

### Remaining work — PRIORITIZED for next agent

**CRITICAL RULES FOR THE NEXT AGENT:**
1. **FOLLOW THE METHOD**: First diff → 3+1 agents → gdb diagnostic → fix → retest → repeat. IT WORKS. Every shortcut fails. Phase 34 found 7 bugs with this method. Phase 33 found 4. Phase 32 found 5. The method is PROVEN over 16+ bugs.
2. **Rule 6 IS NOT OPTIONAL**: Verify EVERYTHING. Phase 34 proved "irreducible" labels WRONG: the 2 "irreducible" 1e-3 failures were a simple inverted condition (Trap 114). Phase 33's "irreducible MT=301 cascades from sigma1" was wrong — the real bug was below-threshold bracket straddling. READ THE FORTRAN.
3. **NOTHING IS IRREDUCIBLE**: Every diff found so far was a real bug. The pattern "values are at the FP precision floor" has been wrong EVERY TIME — there was always a logic error or missing computation step.
4. **3+1 pattern**: Launch 3 RESEARCH subagents in parallel. Keep 1 slot for the Julia runner. NEVER run Julia in parallel (Rule 4).
5. **Per-MT decomposition**: For damage/heating diffs, decompose by MT first (patch Fortran nheat to print per-MT dame at the failing energy). The dominant contributing MT reveals the bug. Phase 34 found MT=107 capdam as the dominant source — missing threshold guard.

**T01 PASSES at 1e-5. To pass at 1e-7 (= 1e-9), need to fix 355 lines. All remaining diffs are > 1e-7 but < 1e-5.**

**Per-section failure breakdown (355 lines at 1e-9)**:

| Section | Diffs | Worst | Root cause | Approach |
|---------|-------|-------|------------|----------|
| MT=444 | 164 | 4.4e-6 | ±1 ULP cascade through damage×sigma | Per-MT decomposition at worst energy, trace GL quadrature |
| MT=301 | 96 | 1.0e-6 | Cascades from MT=1 sigma1 | Fix MT=1 first, MT=301 follows |
| MT=1 | 46 | 2.3e-6 | sigma1 broadening FP accumulation | Match h-function/f-function intermediate rounding in sigma1_at |
| MF6/MT=229 | 30 | 6.0e-7 | sigl phase 1 FP accumulation at high IEs | Trace sigl adaptive stack at worst IE |
| MF6/MT=221 | 16 | 7.7e-7 | Free gas cosine FP | Same as MF6/MT=229 but free_gas_kernel |
| MF3/MT=229 | 1 | 4.2e-6 | coh window at grid edge | Streaming window boundary behavior |
| MF3/MT=221 | 1 | 9.2e-6 | emax boundary | Fortran tpend emax handling |
| MF3/MT=2 | 1 | 4.8e-7 | sigma1 ±1 ULP | Same class as MT=1 |

**Priority 1 — MT=444 (164 diffs, worst 4.4e-6)**:
All diffs are now < 5e-6. The dominant source was capdam below-threshold bracket corruption (fixed in Phase 34). Remaining diffs are ±1 ULP cascaded through damage×sigma products. Per-MT decomposition at the worst energy (E≈9.36 MeV) showed MT=2/51/52 damage PERFECT, MT=107 now close after capdam fix. The residual is cumulative ±1 ULP across many MTs.

**Approach**: These may genuinely be at the FP floor now (4.4e-6 worst = ~1 ULP at 7 sigfigs). But Rule 6 says verify. Trace per-MT damage decomposition at the worst-error energy. If ALL per-MT contributions match to < 1e-6, the total ±1 ULP comes from summation order (match Fortran's MT iteration order in compute_kerma).

**Priority 2 — MT=1/MT=301/MT=2 (143 diffs, worst 2.3e-6)**:
MT=1 broadened total XS differs by ±1 ULP from sigma1_at Doppler broadening kernel. These cascade to MT=301 (96 diffs) through KERMA = h × σ. The sigma1 kernel computes an integral involving h-function and f-function terms. If the intermediate accumulation order differs between Julia and Fortran, the result differs by ±1 ULP.

**Approach**: Patch Fortran broadr sigma1 (bsigma function) to print intermediate h/f values at one of the 46 failing energies. Compare with Julia's sigma1_at. The bsigma function has an integration loop with multiple terms — match the EXACT loop order and intermediate rounding.

**Priority 3 — MF6/MT=229 + MF6/MT=221 (46 diffs, worst 7.7e-7)**:
Equi-probable cosine values differ by ±1 ULP at high incident energies. The sigl_equiprobable adaptive linearization and CDF inversion accumulate FP errors through the angular integration. These are right at the 1e-7 boundary.

**Approach**: Trace Fortran sigl phase 1 linearization at one of the worst IEs. Print each (mu, sig) pair in the adaptive stack. Compare with Julia. Find the first intermediate value that diverges.

**Priority 4 — Generate oracle caches for untested RECONR tests**:
31 tests RAN_OK without oracles. Many share materials with BIT-IDENTICAL tests.

**Priority 5 — Grind BROADR to bit-identical on more tests**.

---

### Phase 35: T02 broadr pipeline test script (PARTIAL — NOT A PROPER PORT)

**What was done**: Created `test/validation/t02_pipeline.jl` — a manually-wired test script that runs reconr→broadr(3 temperatures) for T02 (Pu-238, MAT=1050) and compares against the `after_broadr.pendf` oracle.

**T02 pipeline results**:
- RECONR: 17/17 MTs BIT-IDENTICAL, 3567 pts exact
- BROADR grids match Fortran exactly: 2925/2592/2418 pts at 300K/900K/2100K
- MT=2 (elastic): PERFECT at all 3 temperatures
- 13 non-broadened MTs: ALL PERFECT at all 3 temperatures (39/39)
- MT=18 (fission): 99.3-99.6% (±1 ULP at URR boundary ~70-82 eV)
- MT=102 (capture): 98.4-99.6% (same class)
- MT=1 (total): 78-81% (sigma1 FP accumulation cascade — same class as T01)
- Byte-identical: 12,552/13,133 data lines (95.6%)
- Tolerance: 0 failures at 1e-5, 581 at 1e-9

**T02 broadr parameters** (from Fortran output, verified):
- thnmax = 200.0 eV (resolved resonance range upper limit)
- nreac = 3: broadened partials are MT=2 (elastic), MT=18 (fission), MT=102 (capture)
- Sequential broadening: 0K→300K (T_eff=300, alpha=9135), 300K→900K (T_eff=600, alpha=4568), 900K→2100K (T_eff=1200, alpha=2284)
- Total (MT=1) broadened separately via sigma1_at on same grid (Trap 53)
- All broadened XS rounded to 7 sigfigs via round_sigfig(x,7,0) before format_endf_float

**Trap 120 (NEW — FIXED)**: Fortran broadr rounds all broadened cross section values to 7 sigfigs via sigfig(x,7,0) at broadr.f90 line 980 before writing through a11. This creates trailing zeros in the 9-sigfig fixed format, causing a11 to fall back to 7-sigfig scientific notation. Without this rounding, Julia produces 9-sigfig fixed format (e.g., `31584.1069`) while Fortran produces 7-sigfig scientific (`3.158411+4`). Non-broadened MTs pass through from reconr with original 9-sigfig precision unaffected.

**Files created**:
- `test/validation/t02_pipeline.jl` — T02 pipeline test script

**CRITICAL WARNING — THIS IS NOT A PROPER PORT**:

The T02 pipeline script (like T01's) is a **manually-wired test harness**, NOT a faithful port of Fortran NJOY. Problems:

1. **Hardcoded parameters**: thnmax=200.0, nreac=3, broadened MTs=[2,18,102] are all hardcoded from Fortran diagnostic output. A proper broadr() must auto-compute these from the ENDF/PENDF tapes.

2. **No input deck parsing**: The script doesn't read an NJOY input deck. It hardcodes mat=1050, err=0.005, temps=[300,900,2100].

3. **No PENDF file output**: The script compares data values against the oracle but doesn't produce a complete multi-temperature PENDF file.

4. **No moder execution**: The Fortran T02 chain starts with `moder 20 -21/` (ASCII→binary). The script skips this entirely.

5. **Comparison is against oracle, not reference tape**: The official T02 test compares against `referenceTape28` (unresr output) and `referenceTape29` (groupr output). This script only compares the broadr stage against `after_broadr.pendf`.

6. **Helper functions in test script**: `_format_tab1_data`, `_parse_broadr_oracle`, `_compare_mf3` live in the test script, not the library.

**What's needed for a true T02 port**: A proper broadr() top-level function, proper moder handling, input deck parsing, and full pipeline execution (reconr→broadr→unresr→groupr) producing referenceTape28 and referenceTape29. The kernel algorithms (broadn_grid, sigma1_at) are proven correct — the gap is module-level orchestration.

**The project owner has stated the requirement clearly: NJOY.jl must be a 100% faithful drop-in Julia replacement for ALL 23 Fortran NJOY modules, passing ALL 84 reference tests. No module is out of scope. Every module in the input deck must execute. The current test scripts with hardcoded parameters are not acceptable as a final deliverable — they are validation tools only.**

**How to run T02 pipeline**:
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/t02_pipeline.jl
# Expected: 17/17 reconr PERFECT, 95.6% broadr data match, 0 failures at 1e-5
```

---

### Phase 36: Module-Level Orchestration — `run_njoy()` with Tape-Based Architecture

**Goal**: Build a proper `run_njoy(input_path)` function that reads any Fortran NJOY input deck and executes the module chain automatically, with no hardcoded parameters. Tape-based module communication matching Fortran architecture.

**T01 result**: `run_njoy("njoy-reference/tests/01/input")` passes at 1e-5 with 0 failures. 32,962/32,962 lines. 212 diffs at 1e-9 (sigma1 FP accumulation — same class as hand-wired pipeline).

**Design principle**: Each Julia module is a **drop-in replacement** for the Fortran module. It reads input tapes, writes output tapes. No shared mutable state between modules. The tape files are the communication mechanism — exactly like Fortran NJOY. The `run_njoy()` dispatcher is trivial (matching Fortran main.f90's select-case).

**Architecture**:
```
run_njoy(input_path; work_dir)
  ├─ parse_njoy_input(input_path) → NJOYInputDeck (ModuleCall objects)
  ├─ build_tape_manager(work_dir) → TapeManager (unit→path mapping)
  ├─ for each ModuleCall:
  │   ├─ moder_module(tapes, mc) → copy/register tape files
  │   ├─ reconr_module(tapes, params) → reconr() + write_pendf_file()
  │   ├─ broadr_module(tapes, params) → read PENDF, broadn_grid, write PENDF
  │   ├─ heatr_module(tapes, params) → read PENDF+ENDF, compute_kerma, write PENDF
  │   ├─ thermr_module(tapes, params) → read PENDF+SAB, calcem, bragg, write PENDF
  │   └─ groupr: stub (GENDF not validated for T01)
  └─ final_assembly!() → write_full_pendf() with all accumulated data
```

**New files created** (all under `src/orchestration/`):

| File | Lines | Purpose |
|------|-------|---------|
| `src/orchestration/types.jl` | 85 | TapeManager, PENDFTape/Material/Section structs |
| `src/orchestration/input_parser.jl` | 380 | Input deck tokenizer + module param parsers (moved from test/) |
| `src/orchestration/auto_params.jl` | 120 | compute_thnmax, select_broadr_partials, BRAGG_LATTICE_PARAMS |
| `src/orchestration/pendf_io.jl` | 280 | read_pendf, write_pendf_tape, extract_mf3, copy_with_modifications |
| `src/orchestration/pipeline.jl` | 300 | run_njoy() dispatcher + RunContext + collection helpers |
| `src/orchestration/modules/reconr.jl` | 25 | reconr_module: thin wrapper around reconr() |
| `src/orchestration/modules/broadr.jl` | 110 | broadr_module: auto thnmax, partials, sequential temps |
| `src/orchestration/modules/heatr.jl` | 170 | heatr_module: reads ENDF MF4/5/12/13, compute_kerma |
| `src/orchestration/modules/thermr.jl` | 230 | thermr_module: free gas + SAB + Bragg, calcem, MF6 |
| `src/orchestration/modules/moder.jl` | 80 | moder_module + final_assembly! via write_full_pendf |
| `src/endf/readers.jl` | 80 | read_mf12_gammas, read_mf5_evaporation |

**Modified files**:
- `src/NJOY.jl` — includes + exports for orchestration layer
- `src/processing/pendf_writer.jl` — L2=99 fix for redundant MTs in override_mf3 path

**Key auto-computed parameters** (no hardcoded material-specific values):

| Parameter | Source | Function |
|-----------|--------|----------|
| thnmax | Lowest inelastic threshold in MF3 sections | `compute_thnmax()` |
| Broadr partials | Scan for MT=2/18/102 in PENDF | `select_broadr_partials()` |
| Bragg lattice (a, c, σ_coh, A_mass) | Lookup table by MAT (matches Fortran sigcoh) | `lookup_bragg_params()` |
| Z (atomic number) | ZA / 1000 from PENDF HEAD | `extract_Z()` |
| MF12 gamma data | `read_mf12_gammas()` from ENDF | New reader in endf/readers.jl |
| MF5 evaporation (u, θ) | `read_mf5_evaporation()` from ENDF | New reader in endf/readers.jl |
| MF12/MF13 passthrough | Linearized onto reconr grid via `interpolate()` | `_linearize_mf_section()` |

**Bugs found and fixed during Phase 36**:

1. **MF12 NK interpretation (CRITICAL)**: `read_mf12_gammas` read NK-1 per-gamma subsections instead of NK. NK is the number of photon transitions, not total subsections. The total yield TAB1 is subsection 0; gammas are 1 through NK. Missing the 3rd gamma (1.2625 MeV) caused 2.4% systematic KERMA error. **Fix**: `for k in 1:nk` instead of `for k in 2:nk`.

2. **_interp_to_grid! threshold boundary**: Used `E <= src_e[1]` → zero at E == first energy. Should be `E < src_e[1]`. Caused zero KERMA at E=1e-5 eV.

3. **MF12/MF13 passthrough**: Raw ENDF lines (13+43) instead of reconr-linearized lines (348+139). Fortran emerge interpolates MF12/MF13 via gety1 onto the reconr grid. **Fix**: `_linearize_mf_section()` evaluates TabulatedFunction at each reconr energy, applies `round_sigfig(y, 7, 0)`, formats as ENDF TAB1 with proper sequence numbers.

4. **MF13 TAB1 L2 field**: Wrote 0 instead of 2 (LF=2 = tabulated distribution). The `read_mf13_sections` doesn't capture LF from the TAB1 record. **Fix**: Hardcode L2=2 for MF13 (always tabulated in reconr output).

5. **MF12/MF13 HEAD ZA/AWR**: Wrote 0.0 instead of material ZA/AWR. **Fix**: Pass ZA/AWR from reconr result.

6. **MF12/MF13 sequence numbers**: Wrote "    0" instead of incrementing 1,2,3... Tolerance test parses integers from sequence column. **Fix**: Track sequence counter in `_linearize_mf_section`.

7. **L2=99 for MT=4 in override_mf3 path**: `write_full_pendf` only set L2=99 for MT=1 in the broadened-override path, not for MT=4 or MT=103-107. **Fix**: Apply same redundant-detection logic as the reconr path. This fix is in `src/processing/pendf_writer.jl`.

8. **ThermrParams parser field ordering**: Card 2 fields were shifted — nbin was at position 3 (not ntemp). **Fix**: Corrected field mapping to match Fortran: matde, matdp, nbin, ntemp, iinc, icoh, iform, natom, mtref, iprint.

9. **Input parser '/' inside quotes**: The tokenizer stripped '/' from ALL cards including quoted strings like `'pendf tape for c-nat from endf/b tape 511'`. **Fix**: `_find_slash_outside_quotes()` skips '/' inside single/double quotes.

**Trap 121 (NEW — FIXED)**: MF12 NK=3 means 3 photon transitions = 3 per-gamma TAB1 subsections AFTER the total yield subsection. Total subsections = NK+1. The reader must iterate `1:nk` after skipping the total, not `2:nk`.

**Trap 122 (NEW — FIXED)**: Fortran emerge interpolates MF12/MF13 onto the reconr energy grid via gety1 and applies sigfig(sn,7,0) before writing. The PENDF MF12/MF13 is NOT a raw passthrough of the ENDF — it's a linearized TAB1 on the same grid as MF3. Julia must do the same: evaluate TabulatedFunction at each reconr energy, round to 7 sigfigs, format as ENDF lines.

**Trap 123 (NEW — FIXED)**: The `write_full_pendf` L2=99 logic for redundant MTs only applied in the reconr-path branch (lines 286-307), NOT the override_mf3 branch (lines 260-273). When broadr copies MT=4 through to its output PENDF and it ends up in override_mf3, the override branch used L2 from the ENDF (L2=1) instead of L2=99. Both branches must detect redundant MTs.

**How to run T01 via orchestrated pipeline**:
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e '
using NJOY
tapes = run_njoy("njoy-reference/tests/01/input"; work_dir="/tmp/t01_orch")
println("Output: ", resolve(tapes, 25))
println("Lines: ", countlines(resolve(tapes, 25)))
'
# Expected: 32962 lines, matches referenceTape25
# Tolerance test:
python3 -c "
import re; from math import isclose
fp=re.compile(r'([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)')
ref=open('njoy-reference/tests/01/referenceTape25').readlines()
trial=open('/tmp/t01_orch/tape25').readlines()
n=min(len(ref),len(trial))
for rtol in [1e-9,1e-7,1e-5,1e-4,1e-3]:
    f=0
    for i in range(n):
        r=ref[i]; t=trial[i]
        if r==t: continue
        rF=fp.findall(r); tF=fp.findall(t)
        if not rF:
            if r.rstrip()!=t.rstrip(): f+=1
            continue
        if len(rF)!=len(tF): f+=1; continue
        ok=True
        for a,b in zip(rF,tF):
            try: rv=float(a); tv=float(b)
            except: continue
            if not isclose(rv,tv,rel_tol=rtol,abs_tol=1e-10): ok=False; break
        if not ok: f+=1
    print(f'rel_tol={rtol:.0e}: {f} lines fail ({f/n*100:.1f}%)')
"
# Expected: 0 at 1e-5, 212 at 1e-9
```

**Remaining work — PRIORITIZED for next agent**:

### 1. Extend run_njoy to T02 and beyond

T02 (Pu-238, reconr→broadr with 3 temperatures) should work with the existing orchestration since `broadr_module` supports sequential multi-temperature broadening. Test:
```bash
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/02/input"; work_dir="/tmp/t02_orch")'
```
The T02 chain continues to unresr→groupr which are stubbed. The broadr stage should match the oracle.

### 2. Grind T01 from 1e-5 → 1e-9

212 remaining diffs at 1e-9. Breakdown:
- MT=444: ~100 diffs (damage ±1 ULP cascade)
- MT=301: ~50 diffs (KERMA cascade from MT=1)
- MT=1: ~40 diffs (sigma1 broadening FP accumulation)
- MF6/MT=229: ~15 diffs (sigl FP at high IEs)
- MF6/MT=221: ~5 diffs (free gas cosine FP)
- MF3/MT=229: 1 diff (coh window edge)
- MF3/MT=2: 1 diff (sigma1 ±1 ULP)

Use the 3+1 agent grind method. Rule 6: nothing is irreducible.

### 3. Implement remaining module stubs

- `groupr_module`: multigroup averaging (needed for T01 GENDF, T02 chain)
- `unresr_module`: unresolved self-shielding (needed for T02 chain)
- `acer_module`: ACE format output (T14, T48, T50-54)
- `gaspr_module`: gas production (T45)
- `leapr_module`: S(α,β) generation (T22, T23, T33)

### 4. Generate oracle caches for untested tests

31 reconr tests RAN_OK without oracles. `run_njoy` can now be used to test these.

### 5. Key files for the next agent

| File | What it does |
|------|-------------|
| `src/orchestration/pipeline.jl` | `run_njoy()` entry point + RunContext + helpers |
| `src/orchestration/modules/broadr.jl` | broadr_module with auto thnmax + partials |
| `src/orchestration/modules/heatr.jl` | heatr_module with ENDF readers |
| `src/orchestration/modules/thermr.jl` | thermr_module for free gas + SAB |
| `src/orchestration/auto_params.jl` | Parameter auto-computation |
| `src/orchestration/pendf_io.jl` | PENDF read/write/extract/modify |
| `src/orchestration/input_parser.jl` | Input deck tokenizer (moved from test/) |
| `src/processing/pendf_writer.jl` | `write_full_pendf` — the final PENDF assembly |
| `test/validation/t01_pipeline.jl` | Original hand-wired pipeline (still useful for comparison) |

---

### Phase 37: Multi-temp broadr, multi-material reconr, MF=23 photoatomic support

**Goal**: Make `run_njoy` work for T02 (multi-temp broadr) and T03 (photoatomic, multi-material reconr).

**9 bugs found and fixed:**

1. **Multi-temperature broadr output (CRITICAL)**: `broadr_module` only wrote the final temperature. Fortran writes ALL temperatures as separate MAT blocks. Added `write_broadr_pendf()` in `broadr.jl` that writes one complete MAT block per temperature (MF1+MF2+MF3+MEND). Each block has correct directory NC values, temperature in MF1/MT451 CONT and MF2/MT152.

2. **Sequential T_eff for broadening (CRITICAL)**: `broadr_module` used `alpha = awr / (bk * T)` (full temperature). Fortran uses `T_eff = T_new - T_old` for sequential broadening. Without this, 900K and 2100K grids differed (2540 vs 2592, 2307 vs 2418). Fix: `t_eff = temp - t_old`.

3. **thnmax from eresh (CRITICAL for T02)**: `resolve_thnmax(0)` computed lowest inelastic threshold (44 keV for Pu-238). Fortran broadr.f90 lines 427-431: when `thnmx=0`, `thnmax = min(eresh, lowest_threshold)` where `eresh` is from MF2/MT151. For Pu-238: `eresh=200 eV → thnmax=200`. Added `eresh` kwarg to `resolve_thnmax`, `_extract_eresh_from_pendf()` to read EH from MF2/MT151.

4. **Multi-material reconr (T03)**: `reconr_module` only processed one MAT. T03 has MAT=1 (hydrogen) + MAT=92 (uranium). Added `ReconrMatSpec` struct, multi-material loop in `reconr_module`, `_write_reconr_mat_block()` for appending MAT blocks to same PENDF.

5. **parse_reconr ncards handling (T02 regression fix)**: T02 input has `ncards=3` (3 description lines after err). The parser tried to parse description lines as MAT numbers, crashing. Fix: read `ncards` from card3 field 2, skip that many cards per material.

6. **MF=23 reader**: Added `read_mf23_sections()` in `reconr_types.jl` — identical format to MF=3 (HEAD + TAB1), returns `MF3Section` with `mf=23`.

7. **Photoatomic reconr path**: When no MF=3 but MF=23 exists (photoatomic), builds union grid from MF=23 breakpoints, linearizes to INT=2, evaluates each section with threshold zero at first energy (Fortran emerge line 4794), computes MT=501 as sum of partials excluding MT=515/517.

8. **MF-aware PENDF writer**: `_write_legacy_mf3()` now uses the section's actual `mf` field (3 or 23) instead of hardcoding 3. Directory entries also use correct MF via `mt_to_mf` dict passed through `_write_legacy_mf1()`.

9. **MT=501/460 skip in lunion_grid**: Both the pre-pass breakpoint loop and the main bisection loop now skip MT=501 (photoatomic total, redundant) and MT=460 (delayed photon total), matching Fortran reconr.f90 lines 1875-1877.

**T02 broadr results** (multi-temperature):
```
Grid sizes: 2925/2592/2418 — EXACT MATCH with oracle ✓
Data lines: 13133/13133 — exact count match ✓
Data exact: 86.0% (PENDF round-trip precision)
Tolerance:  0 failures at 1e-5 ✓
```

**T03 reconr results** (photoatomic, multi-material):
```
MF=23 data: 2112/2112 lines — exact count match ✓
Byte-identical: 99.5% (2102/2112)
Tolerance 1e-5: 2 failures (MAT=92 MT=516 grid ±1 ULP from lunion vs emerge)
MAT=1 MT=501: 100% PERFECT (225/225 lines)
```

**T03 remaining 2 failures at 1e-5** (investigated, root cause found):

The Fortran reconr **skips lunion entirely** for photoatomic materials (`lrp ≠ 1`). The grid comes from emerge's processing of MF=23 breakpoints directly — NOT from lunion's panel bisection. Julia uses lunion_grid + `_linearize_mf3!` which produces a similar but not identical grid. Confirmed via gdb: label 140 is NEVER reached for MF=23 sections in T03. The grid points at E≈1.024e6 (MAT=92 MT=516) differ by ±1 in 7th sigfig between Julia's `_linearize_mf3!` midpoints and Fortran's emerge evaluation grid.

**To fix the remaining 2 diffs**: need to understand how Fortran emerge builds the grid for photoatomic materials (when lunion is skipped). The Fortran emerge reads each section's TAB1 from the scratch tape and evaluates at the breakpoints of the UNION grid. The union grid for photoatomic is built by emerge itself (not lunion) — each section's breakpoints contribute to the common grid. The emerge grid construction for photoatomic needs to be matched exactly.

**Trap 124 (NEW — FIXED)**: Fortran broadr uses `T_eff = T_new - T_old` for sequential multi-temperature broadening, NOT the full temperature. `alpha = awr / (bk * T_eff)`. Without this, second and third temperature grids differ significantly.

**Trap 125 (NEW — FIXED)**: Fortran broadr.f90 lines 427-431: when `thnmx=0`, `thnmax = eresh` (resolved resonance range boundary from MF2), then reduced by lowest inelastic threshold. For Pu-238: eresh=200 eV, no thresholds below 200 → thnmax=200. For C-nat: eresh=2e7 (from reconr MF2/MT151), lowest threshold=4.81e6 → thnmax=4.81e6.

**Trap 126 (NEW — FIXED)**: Fortran reconr.f90 line 4794: `if (thresh.gt.one.and.abs(thresh-eg).lt.test*thresh) sn=0`. At the exact threshold energy (E=1000 eV for photoatomic sections), XS is forced to 0. This applies to ALL sections including MF=23 partials. Without this, the MT=501 total was 12.24 instead of 0 at E=1000.

**Trap 127 (NEW — INVESTIGATED)**: Fortran reconr SKIPS lunion entirely for photoatomic materials (`lrp ≠ 1` at line 304). The grid comes from emerge processing the MF=23 breakpoints directly. Julia currently uses lunion_grid + `_linearize_mf3!` which produces a similar but not identical grid. The remaining 2 diffs at 1e-5 (MAT=92 MT=516 at E≈1.024e6) are from this fundamental path difference.

**PENDF round-trip and SEND handling (Trap 128)**: `read_pendf` includes the SEND record in section lines. When copying non-broadened sections verbatim, the SEND must be stripped (via `_count_data_lines`) to avoid double SEND. Same applies to MF2 sections and directory NC computation.

**Files changed**:
- `src/orchestration/modules/broadr.jl` — multi-temp collection, T_eff, eresh, write_broadr_pendf + helpers
- `src/orchestration/modules/reconr.jl` — multi-material, _write_reconr_mat_block
- `src/orchestration/input_parser.jl` — ReconrMatSpec, parse_reconr with ncards/ngrid
- `src/orchestration/auto_params.jl` — resolve_thnmax with eresh kwarg
- `src/processing/reconr_types.jl` — read_mf23_sections
- `src/processing/reconr.jl` — photoatomic path, MF=23 in all_lunion_sections
- `src/processing/reconr_grid.jl` — MT=501/460 skips
- `src/processing/pendf_writer.jl` — MF-aware writer, MT=501 handling, _collect_reactions photoatomic

### Phase 38: T04 pipeline — ERRORR + GROUPR orchestration, 2/3 tapes BIT-IDENTICAL

**Goal**: Run T04's full chain (moder→reconr→errorr→groupr→errorr) in the completely modular faithful Fortran drop-in fashion, producing output tapes matching referenceTape23/24/25.

**T04 test structure** (from input deck):
- Chain: moder→reconr→errorr→groupr→errorr
- Input: tape20 (t511, U-235 ENDF/B-V)
- Comparison: tape23 vs referenceTape23, tape24 vs referenceTape24, tape25 vs referenceTape25
- MAT=1395, err=0.10 (10% tolerance for errorr test problem)
- ERRORR #1: MF33 covariance on 9-group structure (from MF33 energies intersected with [1.0, 1000.0])
- GROUPR: MT=452 (nubar) on LANL-30 (30 groups), ign=3, iwt=3 (1/E weight)
- ERRORR #2: MF31 nubar covariance on 7-group user structure [1, 10, 100, ..., 1e7], appended to tape23

**New modules implemented:**

| Module | File | Lines | Purpose |
|--------|------|-------|---------|
| ErrorrParams + parse_errorr | `src/orchestration/input_parser.jl` | +115 | Parse ERRORR input cards |
| GrouprParams + parse_groupr | `src/orchestration/input_parser.jl` | (above) | Parse GROUPR input cards |
| errorr_module | `src/orchestration/modules/errorr.jl` | 531 | Full ERRORR orchestration |
| groupr_module | `src/orchestration/modules/groupr.jl` | 297 | Full GROUPR orchestration |
| Pipeline dispatch | `src/orchestration/pipeline.jl` | +8 | Dispatch for errorr/groupr |
| T04 validation | `test/validation/t04_pipeline.jl` | 145 | Full pipeline test script |

**Results:**
| Tape | Lines | Status | Details |
|------|-------|--------|---------|
| **tape23** (errorr #1) | 82/82 | **BIT-IDENTICAL** | MF3 XS + MF33 covariance all match |
| **tape24** (groupr) | 74/74 | **BIT-IDENTICAL** | Nubar + flux + σf_avg all match |
| **tape25** (errorr #2) | 119/119 | 16 diffs | Copied tape23 matches; MT452 nubar cov needs full covcal |

**Key breakthrough: 1% energy stepping** — The Fortran errorr/groupr `egtwtf` function for iwt=3 (1/E weight) returns `enext = s101 * e` where `s101 = 1.01`. This creates dense sub-panels every 1% of energy within each PENDF panel. Both modules use trapezoidal integration on these sub-panels. Replicating this stepping (instead of piecewise-linear panel integration) produced **exact bit-identical match** on both tape23 and tape24.

**Bugs fixed:**

1. **read_mf33 NC bounds crash (Bug 14)**: NC sub-subsections with odd N1 values (ENDF/B-V quirk) caused BoundsError. Fixed with safe bounds checking.

2. **expand_lb1 off-by-one (Bug 15)**: For LB=0/1/2, energies and data have equal length (E,F pairs). The loop iterated to `length(fk)` but needed `min(length(fk), length(ek)) - 1`.

3. **Covariance double-counting (Bug 16)**: MF33/MT18 Sub 3/4 have NC+NI blocks — NI blocks from NC-derived sub-sections shouldn't be directly summed with standalone Sub 1 NI blocks. Fix: only expand NI blocks from sub-sections with NC=0.

4. **Cross-covariance placement (Bug 17)**: MT18/MT102 cross-covariance stored only in the lower-MT section (MT18). MT102 section has NL=1 (self only), not NL=2.

5. **MF3 XS format (Bug 18)**: ERRORR output uses `format_endf_float` (a11 format), not free-format `@sprintf`. Confirmed via Fortran endf.f90 lineio→a11 chain.

6. **Sequence numbers (Bug 19)**: Fortran uses continuous sequence numbers across MF1/MF3/MF33 sections. Julia was restarting at 1 for each MF.

7. **Integration method (Bug 20 — CRITICAL)**: Group-averaged XS used `group_integrate` with piecewise-linear σ/E → 0.3% error. Fortran's `epanel` uses trapezoidal with 1% stepping sub-panels. Fix: `_group_average_inv_e` with `s101 = 1.01` sub-stepping.

8. **GENDF record format (Bug 21)**: For ratio quantities (MT=452), the GENDF stores (flux, nubar, σf_avg) NOT (flux, nubar, ∫ν·σf·W dE). Fortran `displa` transforms: production/reaction = nubar, reaction/flux = σf_avg.

9. **parse_errorr ign=1 (Bug 22)**: ign=1 means "read user group structure" (same as ign<0). The second errorr call uses ign=1 with 7 user energy boundaries.

10. **FEND duplication (Bug 23)**: Second errorr call copied tape23 including FEND, then added another FEND. Fix: only add MEND after copy.

**tape25 remaining 16 diffs** — Two categories:
1. **Nubar XS (lines 89-90)**: 0.1% off in high-energy groups (regrouping from LANL-30 to 7-group)
2. **Nubar covariance (lines 96+)**: 20x off — requires full ERRORR `covcal` pipeline with union-grid expansion + LB=3 correlation handling + flux-weighted collapse. Direct LB=5 block expansion gives 2.4e-6 vs reference 4.8e-5.

**Trap 129 (NEW — FIXED)**: Fortran errorr/groupr `egtwtf` for iwt=3 returns `enext = 1.01*E`. Both modules use trapezoidal rule with these dense sub-panels for 1/E weighting. Without this, group-averaged XS differ by 0.02-0.3% — enough to fail at 1e-5 tolerance. Julia's `group_integrate` with piecewise-linear σ/E is analytically correct but gives DIFFERENT results from the Fortran's approximate quadrature on dense sub-panels.

**Trap 130 (NEW — FIXED)**: For GENDF ratio quantities (MT=251-253, MT=452/455/456), the Fortran `displa` subroutine transforms raw integrals to (flux, yield, cross_section). The panel accumulates 3 integrals: flux=∫W dE, production=∫ν·σf·W dE, reaction=∫σf·W dE. The displa transforms: ans(2) = production/reaction = nubar, ans(3) = reaction/flux = σf_avg.

**Trap 131 (NEW — INVESTIGATED)**: MF31/MT452 Sub 1 has 5 NI blocks (LB=5, LB=1, LB=2, LB=2, LB=2). The LB=5 block has the symmetric covariance matrix (~5e-5 values). The LB=1/2 blocks add large diagonal corrections (~4e-3). The Fortran ERRORR covcal processes all blocks via union-grid expansion + flux-weighted collapse, which properly distributes the corrections. Direct expansion onto the output grid gives wrong values (20x off for LB=5 alone, or 200x with all blocks summed).

**How to run T04:**
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/t04_pipeline.jl
```

---

### Phase 39: UNRESR + CCCCR orchestration — MT152 BIT-IDENTICAL, T02 end-to-end

**Goal**: Wire unresr and ccccr into `run_njoy()` so T02 runs its full chain (moder→reconr→broadr→moder→unresr→groupr→ccccr→moder→moder).

**9 bugs found and fixed (all via 3+1 agent pattern + Fortran gdb diagnostics):**

1. **ajku Faddeeva function (CRITICAL)**: Julia's `ajku` was a naive 40-point midpoint rule that gave ~2.4x wrong xj values. Fortran uses the complex probability integral (Faddeeva function) with precomputed 62×62 tables + 8-point Gaussian quadrature. Ported full algorithm: `_uw` (complex w function via asymptotic/Taylor series), `_UNRESR_W_TABLE` (precomputed table), `_quikw` (fast table lookup with 4-region rational approximation), `ajku` (Gaussian quadrature). **After fix: ajku matches Fortran to machine precision.**

2. **bondarenko_xs multi-sequence accumulation (CRITICAL)**: Julia computed `sigu += sigm * t_part_seq / (1 - ttj_seq)` per-sequence independently. Fortran accumulates `tl`, `tj`, `tk` across ALL sequences first, applies cross-sequence correction `(1 - xj_other)`, then divides by global `(1 - yj)` (lines 1152-1191). Rewrote to two-phase: phase 1 accumulates per-sequence arrays, phase 2 sums with cross-isotope correction matching Fortran exactly. **After fix: MT152 426/426 BIT-IDENTICAL.**

3. **sigbt missing from total**: `bondarenko_xs` line 428 only added `bkg[1]` to total XS, missing `spot + sint` (potential scattering + interference). Fortran line 1195: `sigu(1,is) = sigf(1,is,4) + sigbt` where `sigbt = sigbkg(1) + spot + sint`.

4. **URR energy grid construction**: ENDF has only 3 fission-width energy nodes (200, 500, 10000 eV). Fortran `rdunf2` synthesizes 23-point grid by subdividing with 78-point equi-lethargy reference grid (`egridu`) using `wide=1.26` threshold and `step=1.01` ratio. Implemented `_build_unresr_grid` + `_UNRESR_EGRIDU` constant. Grid matches oracle exactly.

5. **MF2 URR reader format**: Mode 11 (LFW=1) subsection header is a LIST record (C1=SPI, C2=AP, N1=NE, N2=NLS) with NE energy values in the body, not a CONT. Matched existing `_read_urr_lfw1` pattern from reader.jl.

6. **MF3 backgrounds from ENDF**: Fortran `rdunf3` reads MF3 backgrounds from the ENDF tape (nendf), not the PENDF. And skips the HEAD record before reading the TAB1. For Pu-238, MF3 in the URR range is 0 (ENDF convention: resonance XS comes from MF2, not MF3).

7. **Broadr 90-char lines**: `write_broadr_pendf` used `@printf "%32s%11d..."` producing 90-char lines instead of ENDF-standard 80. Fixed to `"%22s%11d..."` (22+44+14=80).

8. **Broadr extra MEND lines**: Broadr wrote `(0,0,0,0)` MEND records after each MF's FEND, creating 6 extra lines. Removed intermediate MENDs; kept only one MEND per temperature block.

9. **TPID and descriptions passthrough**: Added `title` and `descriptions` parameters to `reconr_module` → `write_pendf_file` → `_write_legacy_mf1`. Pipeline now passes reconr input deck description cards through to the PENDF.

**New files created:**
- `src/orchestration/modules/unresr.jl` (~700 lines): Full unresr_module with MF2 URR reader, MF3 background reader, 78-point grid construction, Bondarenko computation, MT152 ENDF format writer, PENDF copy-with-insert
- `src/orchestration/modules/ccccr.jl` (~90 lines): ccccr_module stub (writes placeholder CCCC files, doesn't crash)

**Files modified:**
- `src/processing/unresr.jl`: Replaced naive `ajku` with Faddeeva function port (~200 lines). Rewrote `bondarenko_xs` with two-phase accumulation + sigbt fix.
- `src/orchestration/pipeline.jl`: Added `:unresr` and `:ccccr` dispatch cases. Title/descriptions passthrough.
- `src/orchestration/modules/broadr.jl`: Fixed 90→80 char lines. Removed extra MEND lines.
- `src/orchestration/modules/reconr.jl`: Added title/descriptions parameters.
- `src/processing/pendf_writer.jl`: Added title/descriptions to `write_pendf` and `_write_legacy_mf1`.
- `src/NJOY.jl`: Includes for new modules.

**T02 results:**
```
tape28: 13870 lines (ref: 13873) — 3 lines short
MT152: 426/426 BIT-IDENTICAL (zero failures at any tolerance)
MF3 data: correct values, correct positions
```

**Remaining 3-line gap**: Julia's reconr doesn't write preliminary MF2/MT152 (the Fortran reconr does). This causes broadr temp blocks 2-3 to have NXC=19 instead of 20 (missing 1 directory entry per block = 3 lines). Fix requires adding MF2/MT152 to the reconr PENDF writer (`write_pendf` / `_write_legacy_mf1`).

**Trap 132 (NEW — FIXED)**: Julia's `ajku` was a naive 40-point midpoint rule — completely wrong algorithm. The Fortran uses the complex probability integral (Faddeeva function) via precomputed 62×62 tables with bilinear interpolation (for |z|²<36), rational approximation (36-144), simpler rational (144-10000), and asymptotic form (>10000). The 8-point Gaussian quadrature evaluates the integral with 4 `quikw` calls per point. Without this, Bondarenko XS had ~2.4x systematic error.

**Trap 133 (NEW — FIXED)**: `bondarenko_xs` must accumulate `tl`, `tj`, `tk` from ALL sequences first, then sum with cross-sequence correction `(1-xj_other)*abns(ks)` before dividing by global `(1-yj)`. Per-sequence `1/(1-ttj)` gives wrong results at low sigma0 (up to 174% error). The Fortran's two-phase approach (lines 1076-1191) is NOT equivalent to the per-sequence formula.

**Trap 134 (NEW — FIXED)**: Fortran `rdunf2` synthesizes the URR energy grid from a 78-point equi-lethargy reference grid (`egridu`). The ENDF provides only 3 fission-width energy nodes. Between consecutive nodes with ratio > 1.26, intermediate points from `egridu` are inserted using 1.01 step ratio. Result: 23-point grid for Pu-238. The sentinel is 1e6 (placed at END of sorted list via `ilist`), not 0 (at start).

**Trap 135 (NEW — FIXED)**: Broadr `write_broadr_pendf` produced 90-char lines via `@printf "%32s%11d..."`. ENDF standard is 80 chars. Fixed to `"%22s%11d..."` (22 blanks + 4×11 data + 4+2+3+5 trailer = 80). Without this, every line comparison fails because MAT/MF/MT fields are at wrong positions.

**Trap 136 (NEW)**: Julia's reconr doesn't write preliminary MF2/MT152 (infinite-dilution Bondarenko table, nsigz=1). The Fortran reconr does. This causes the broadr PENDF to have NXC=19 (no MT152 directory entry) while the Fortran has NXC=20. The unresr module fixes block 1's directory but blocks 2-3 remain 1 line short each.

**How to run T02:**
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/02/input"; work_dir="/tmp/t02_orch")'
```

---

### Phase 40: T02 tape28 PASSES at 1e-5 — 8 bugs fixed in reconr/broadr/unresr

**Goal**: Get T02's full chain (moder→reconr→broadr→unresr→groupr→ccccr→moder) passing the official tolerance test.

**T02 official test** compares tape28 (unresr PENDF, 13873 lines) AND tape29 (groupr GENDF, 6213 lines) at rel_tol=1e-9.

**Results after Phase 40:**
```
tape28: 13873/13873 lines — STRUCTURAL MATCH ✓
tape28 rel_tol=1e-5: 0 failures ← PASSES ✓
tape28 rel_tol=1e-9: 1836 failures (all ±1 ULP sigma1 broadening class)
tape29: 627/6213 lines — needs groupr MF3 multi-temp/sigma0 + MF6 transfer matrices
T01 tape25: no regression (32962 lines, 0 failures at 1e-5, 212 at 1e-9)
```

**8 bugs found and fixed (all via grind method + Fortran source reading):**

1. **parse_groupr weight function card skip (FIXED)**: For iwt=1/4/5, the parser didn't skip the weight function parameter card before reading the MT list. For T02 (iwt=4), card 6 (`.1 0.025 0.8208e06 1.4e06`) was parsed as MFD=0 → 0 MT requests. Fix: skip 1 card for iwt in (1,4,5). **Impact: groupr processes 17 MTs instead of 0.**

2. **ign=5 RRD 50-group structure (ADDED)**: T02 uses ign=5 (50 groups). Julia only had ign=2,3,4. Added `RRD_50` constant with 51 boundaries computed from Fortran gl5 lethargy array via `sigfig(1e7 * exp(-u), 7, 0)` matching groupr.f90 line 4471. Verified exact match against reference tape29 group boundaries.

3. **Pipeline final_assembly! overwrite (FIXED)**: `run_njoy` called `final_assembly!` unconditionally on the last moder output unit. For T02, this overwrote tape29 (groupr GENDF) with a reconr PENDF. Fix: only run final_assembly! when thermr/heatr data needs merging. **Impact: tape29 preserved as groupr output.**

4. **Unresr multi-block MT152 insertion (FIXED)**: `_write_unresr_pendf` only modified block 1 of the multi-temperature PENDF. Blocks 2-3 were copied without MT152. Complete rewrite as state machine: detects block boundaries via MF1/MT451 HEAD after MEND, modifies MF1 directory (incrementing NXC, inserting MT152 entry), inserts MT152 data after MT151 SEND. Matches Fortran unresr.f90 lines 138-321 (per-temperature loop). **Impact: all 3 blocks now have MT152.**

5. **Reconr preliminary MT152 (ADDED)**: Fortran reconr (lines 5203-5214) writes a preliminary MF2/MT152 section with nsigz=1 (infinite dilution) Bondarenko XS. Julia reconr didn't write it. Added to `write_pendf` via `_write_legacy_mf2`: CONT(ZA,AWR,LSSF,0,0,INTUNR) + LIST(temp,0,5,1,NCP,nunr) + data(sigma0=1e10, per-energy: E,total,elastic,fission,capture,elastic_copy). Directory: NC=3+nunr. Also added `urr_table` and `urr_lssf` to reconr result tuple. **Impact: tape28 MT152 NC=27 matches Fortran (was 141).**

6. **MF2/MT151 EL value (FIXED)**: Julia wrote EL from `iso.ranges[1].EL` (=1.0 eV for Pu-238 resolved range). Fortran reconr recout writes `elow=1e-5` (the ENDF minimum). Fix: hardcode EL=1e-5 matching Fortran. **Impact: -3 failures.**

7. **Broadr TAB1 interp format (FIXED)**: Julia wrote NR/NP as ENDF floats (`2.925000+3 2.000000+0`). Fortran tab1io writes as integers (`2925  2`). Fix: use `@printf "%11d%11d"` for interp record in broadr. **Impact: -12 failures.**

8. **Reconr MF3 HEAD L2 values (FIXED)**: `_write_legacy_mf3` hardcoded L2=99 for ALL MTs. Fortran reconr recout line 5331: `scr(4)=lfs` (level number from ENDF) for non-redundant MTs, L2=99 only for redundant MTs (1, 4, 103-107). Fix: read L2 from `sec.L2` per MT. **Impact: -27 failures.**

**Trap 137 (NEW — FIXED)**: parse_groupr must skip the weight function card for iwt=1/4/5 before reading the MT request list. Without this, the first MT card is parsed as MFD=0 → empty MT list.

**Trap 138 (NEW — FIXED)**: ign=5 is the RRD 50-group structure. Boundaries are computed from gl5 lethargy array: `E = sigfig(ezero * exp(-u), 7, 0)` where `ezero = 1e7`. The gl5 array has 51 values from 27.631 down to -0.6917.

**Trap 139 (NEW — FIXED)**: final_assembly! should only run when thermr/heatr data needs merging. For chains like T02 (reconr→broadr→unresr→groupr), each module writes its own complete tape. Running final_assembly! overwrites the groupr GENDF output.

**Trap 140 (NEW — FIXED)**: Fortran reconr writes preliminary MF2/MT152 with nsigz=1 (infinite dilution) for materials with URR data. This preliminary table is later replaced by unresr with multi-sigma0 data, but the directory NC=3+nunr is PRESERVED by unresr (Fortran unresr line 283: `j` only advances when `new>0`). Julia must write this preliminary MT152 in reconr output.

**Trap 141 (NEW — FIXED)**: Fortran tab1io writes the interpolation record (NBT, INT) as integers using `%11d` format, not ENDF floats. Julia's broadr was using format_endf_float. The difference (`2925  2` vs `2.925000+3 2.000000+0`) causes tolerance failures because the float values parse as different numbers.

**Trap 142 (NEW — FIXED)**: Fortran reconr recout uses `scr(4)=lfs` (level number from ENDF MF3 HEAD) for non-redundant MTs (line 5331) and `scr(4)=99` only for redundant MTs (lines 5234, 5266, 5348). `_write_legacy_mf3` was hardcoding L2=99 for ALL MTs.

**Remaining for tape28 (all ±1 ULP sigma1 class)**:
- 1836 failures at 1e-9, 0 at 1e-5
- All from Doppler broadening FP accumulation order (same class as T01 MT=1)

**Remaining for tape29 (groupr GENDF, 627/6213 lines)**:
- MF3: present (17 MTs × 50 groups) but single-temp single-sigma0 format. Need multi-temp (3) × multi-sigma0 (7) GENDF records.
- MF6: missing entirely. Need transfer matrices for elastic, inelastic, fission.
- MF1: needs ntemp=3, nsigz=7 in HEAD; sigma0 values in directory.
- iwt=4: weight function params not yet used (1/E weight used as default).

**Files changed**:
- `src/orchestration/input_parser.jl` — weight function card skip for iwt=1/4/5
- `src/processing/group_structures.jl` — RRD_50 constant + IGN_RRD50 enum
- `src/orchestration/modules/groupr.jl` — ign=5 dispatch
- `src/orchestration/pipeline.jl` — conditional final_assembly!
- `src/orchestration/modules/unresr.jl` — complete rewrite of _write_unresr_pendf as state machine
- `src/processing/reconr.jl` — urr_table + urr_lssf in result tuple
- `src/processing/pendf_writer.jl` — preliminary MT152, L2 from sec.L2, urr_table plumbing
- `src/orchestration/modules/broadr.jl` — integer TAB1 interp format
- `src/NJOY.jl` — RRD_50 + IGN_RRD50 exports

**How to verify Phase 40:**
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/02/input"; work_dir="/tmp/t02_test")'
# tape28: 13873 lines, 0 failures at 1e-5
# tape29: 627 lines (needs groupr work for 6213)
# T01 regression check:
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/01/input"; work_dir="/tmp/t01_test")'
# tape25: 32962 lines, 0 failures at 1e-5
```

---

### Phase 41: viewr+graph.f90 faithful port — 100% BIT-IDENTICAL PostScript

**Goal**: Complete, faithful, idiomatic Julia port of Fortran viewr.f90 (1669 lines) + graph.f90 (2899 lines) — the NJOY PostScript rendering engine.

**Result**: **9274/9274 lines BIT-IDENTICAL** (100.0%, zero diffs) when viewr reads the correct Fortran plot tape. The T03 referenceTape37 is reproduced byte-for-byte.

**Architecture**: New `src/viewr/` directory with 11 focused files. Single `mutable struct GraphState` carries all ~60 module-level variables from graph.f90 — no global mutable state. All functions take `gs::GraphState` as first argument.

**Files created** (4462 lines total, porting 4568 Fortran lines):

| File | Lines | Fortran equivalent | What it implements |
|------|-------|-------------------|-------------------|
| `src/viewr/constants.jl` | 152 | graph.f90 data decls | CW(128,3) font width tables (Times-Roman, Helvetica, Symbol), IBRGB/IFRGB/ISRGB color tables, page geometry constants |
| `src/viewr/types.jl` | 203 | graph.f90 + viewr.f90 module vars | `GraphState` (56 fields) + `ViewrState` (57 fields) mutable structs |
| `src/viewr/ps_primitives.jl` | 264 | graph.f90:2103-2374 | gplot!, gdone!, newp!, endp!, moveh!, drawh!, fillh!, gset!, gend!, cclip!, nclip!, ncurve! |
| `src/viewr/transforms.jl` | 377 | graph.f90:212-636 | transw, trans3, initp!, window!, endw!, init2!, frame2!, endfr!, vect2!, vect3!, poly2!, poly3! |
| `src/viewr/text.jl` | 651 | graph.f90:1634-2095,2795-2896 | DISSPLA markup engine (text3!, text2!, dchr!, txtlen, ssum, charw, stripv, iget, rget) |
| `src/viewr/axes.jl` | 763 | graph.f90:638-1216 | axis3! (485-line Fortran routine), axis2!, xscale/yscale/zscale, xinvrs/yinvrs |
| `src/viewr/curves.jl` | 343 | graph.f90:1218-1632 | grid2!/grid3!, curv2!, curv3!, hatch! |
| `src/viewr/symbols.jl` | 197 | graph.f90:2376-2793 | dsym! (25 symbol types), circle! |
| `src/viewr/plot_tape.jl` | 163 | viewr.f90 card reading | read_plot_reals, read_plot_string, _read_card_reals, _read_card_string |
| `src/viewr/set2d.jl` | 945 | viewr.f90:757-1524, graph.f90:487-569 | ascalv, set2d!, erbar!, legndb!, tagit!, init3!, set3d! |
| `src/viewr/viewr.jl` | 404 | viewr.f90:42-755 | viewr_render! — main plot tape loop |

**Files modified:**
- `src/NJOY.jl` — includes for 11 viewr engine files in dependency order
- `src/orchestration/modules/viewr.jl` — replaced 367-line ad-hoc stub with 32-line wrapper calling viewr_render!

**Key design decisions:**
- All PS format strings match Fortran exactly (`f6.3`, `f9.2`, `f8.3`, `f7.2`, `i4`, `i5`)
- Coordinate clamping to [-1000, 2000] PS points (prevents PostScript overflow)
- Landscape rotation in software, not PS `rotate`
- `moveh!` always strokes previous path before moveto (PS accumulation model)
- `drawh!` caches linewidth (`wlast`) and dash pattern (`ldash`) to suppress redundant PS commands
- Text rendered character-by-character with per-char `findfont/makefont/setfont/show`
- Log axis guards for non-positive data values (Fortran would produce NaN/garbage, Julia would crash)

**Verification**:
```bash
# Generate Fortran reference plot tape
cd /tmp && mkdir -p t03_fortran && cd t03_fortran
cp ~/Projects/NJOY.jl/njoy-reference/tests/03/input .
cp ~/Projects/NJOY.jl/njoy-reference/tests/resources/gam23 tape30
cp ~/Projects/NJOY.jl/njoy-reference/tests/resources/gam27 tape32
~/Projects/NJOY.jl/njoy-reference/build/njoy < input > output 2>&1

# Run Julia viewr on Fortran plot tape
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=~/Projects/NJOY.jl -e '
using NJOY
open("/tmp/t03_viewr_test.ps", "w") do io
    NJOY.viewr_render!(io, "/tmp/t03_fortran/tape36")
end
jl = readlines("/tmp/t03_viewr_test.ps")
rf = readlines("njoy-reference/tests/03/referenceTape37")
n = min(length(jl), length(rf))
m = count(i -> jl[i] == rf[i], 1:n)
println("$m / $n lines match ($(round(100m/n, digits=1))%)")
# Expected: 9274 / 9274 lines match (100.0%)
'
```

**T03 full pipeline** runs end-to-end: `run_njoy("njoy-reference/tests/03/input")` produces all 7 tapes (31, 33, 34, 35, 36, 37).

### Phase 42: dtfr rewrite — tape34 BIT-IDENTICAL, plot tape 73% match

**Result**: **tape34 (DTF): 340/340 BIT-IDENTICAL** when reading Fortran GENDF. **tape36 (plot): 1187/1628 lines match (72.9%)** with Fortran inputs — first 1182 lines identical, remaining diffs in uranium PENDF overlay thinning.

**Bugs fixed**:
1. **`_dpend` infinite loop (FIXED)**: `gety1_interp` returned `enext == e` when hitting tabulated energy point. Early return now advances to `energies_raw[2]`, while loop uses `<=` not `<`.
2. **Plot tape formatting (FIXED)**: Title trailing spaces, histogram data leading spaces, axis tag E format (`_fmt_e10_2` for Fortran 0p format).
3. **matxsr API break (FIXED)**: Migrated from old `gmat.sections`/`sec.data` to new `gmat.mf23`/`gmat.mf26`/`sec.sigma` + `gmat.flux`.

**Files changed**:
- `src/orchestration/modules/dtfr.jl` — complete rewrite: GENDF reader (MF23+MF26), DTF writer, plot tape writer, dpend fix, fmt_e10_2 helper, GendfMaterial flux field
- `src/orchestration/modules/matxsr.jl` — API migration

**Remaining T03 blockers (priority order)**:

1. **gaminr MT621 (photon heating KERMA)**: Julia writes raw XS (~0.6 b), Fortran writes KERMA (~4488 b·eV). Missing energy multiplication in `src/processing/gaminr.jl`. Compare Fortran `gaminr.f90` subroutine `gheat`. **This is the #1 blocker** — all tape34/36/37 value differences trace to this.

2. **matxsr transfer matrices**: Writes diagonal-only (131 lines vs 211). Need full banded matrix from `sec.transfer` array with proper `jband`/`ijj` computation.

3. **`_dpend` thinning**: Missing Fortran ns states 3-5. Produces 451 vs 535 points for uranium MT=501 overlay. See `dtfr.f90` lines 1241-1268.

4. **gaminr GENDF structure**: 688 vs 604 lines — dense MF26 format (zeros below threshold), wrong section ordering (all MF23 then all MF26 instead of interleaved per MT).

**How to verify**:
```bash
# Generate Fortran oracle
mkdir -p /tmp/t03_fortran && cd /tmp/t03_fortran
cp ~/Projects/NJOY.jl/njoy-reference/tests/03/input .
cp ~/Projects/NJOY.jl/njoy-reference/tests/resources/gam23 tape30
cp ~/Projects/NJOY.jl/njoy-reference/tests/resources/gam27 tape32
~/Projects/NJOY.jl/njoy-reference/build/njoy < input > output 2>&1

# Test dtfr in isolation with Fortran GENDF (isolates dtfr from gaminr bugs)
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e '
using NJOY
mats = NJOY._read_gendf_tape("/tmp/t03_fortran/tape33")
params = NJOY.DtfrParams(33, 34, 31, 36, 1, 1, 0, 5, 12, 4, 5, 16, 1, 0, 0, 0, 0,
  ["pheat"], [NJOY.DtfrEditSpec(1, 621, 1)], 0, 0,
  [NJOY.DtfrMaterial("h", 1, 1, 0.0), NJOY.DtfrMaterial("u", 92, 1, 0.0)])
buf = IOBuffer(); NJOY._write_dtf_output(buf, mats, params)
jl = split(String(take!(buf)), "\n"); ref = readlines("/tmp/t03_fortran/tape34")
n = min(length(jl), length(ref))
println("tape34: $(count(i -> jl[i] == ref[i], 1:n)) / $(length(ref)) BIT-IDENTICAL")
# Expected: 340 / 340 BIT-IDENTICAL
'
```

**Trap 143 (NEW)**: viewr is 100% correct. When T03 produces wrong PS, the bug is in dtfr or upstream (gaminr). Test viewr in isolation with Fortran plot tape.

**Trap 144 (NEW — FIXED)**: `_dpend` hangs when first PENDF energy equals `elow`. Fix: `gety1_interp` early return must return `energies_raw[2]` as `enext`.

**Trap 145 (NEW)**: Julia `%10.2E` → `1.20E+04`. Fortran `E10.2` (0p) → `0.12E+05`. Use `_fmt_e10_2()` for viewr tag values.

**Trap 146 (NEW)**: matxsr uses `gmat.mf23`/`gmat.mf26`/`gmat.flux` now, not old `gmat.sections`/`sec.data`.

**Trap 147 (NEW)**: Test dtfr with Fortran GENDF to isolate dtfr bugs from gaminr bugs. Julia gaminr MT621 is wrong → tape34 shows ~35% match but this is gaminr, not dtfr.

See `worklog/T03_phase2_handoff.md` for full session details.

### Phase 43: gaminr total heating + incoherent normalization + dpend state machine

**Goal**: Get T03 tape37 passing at 1e-5.

**5 bugs found and fixed:**

1. **MT=621 total heating (CRITICAL)**: Julia computed heating only from photoelectric (MT=602). Fortran accumulates from incoherent (MT=504), pair production (MT=516), AND photoelectric. The Fortran gtff sets `ff(1,2)=e` for MT=602, creating a heating column that dspla normalizes in-place before accumulation into `toth`. Added `toth` array and accumulation in `_write_gaminr_tape`. MT=621 group 1 heating went from 7.9 → 4488 eV·barn (matching Fortran). Trap 148.

2. **Incoherent scattering normalization (CRITICAL)**: Julia multiplied σ_MF23 × KN × S(q) = σ² (3-4x too small). Fortran normalizes the angular distribution by siginc (gtff lines 1456-1464), then gpanel multiplies by σ_MF23. Rewrote `_gaminr_incoherent_matrix!` to normalize by siginc before multiplying by σ_MF23. Also fixed: S(q) was squared for incoherent — MF27/MT=504 stores S(q,Z) not F(q), use linearly. Traps 149-150.

3. **GENDF section ordering**: Interleaved MF23/MF26 per MT matching Fortran mtlst/mflst arrays. No FEND between interleaved sections. Trap 151.

4. **MF26 format**: Added flux as position 1 in LIST records. Coherent (MT=502) now ng2=2. Below-threshold skip for MF23/MT=516. Traps 153-155.

5. **_dpend thinning state machine**: Implemented full 5-state machine (Fortran computed GOTO, states 0-5: sort and keep extremes). Plot tape 1627 → 1711 (structural match). Trap 152.

**Results:**
```
tape33 (GENDF): 688 → 632 (target: 604, diff: +28)
tape36 (plot):  1627 → 1711 (target: 1711) ✓ STRUCTURAL MATCH
tape37 (PS):    9038 → 9110 (target: 9274)
MT=621 g1:      7.9 → 4488 (target: 4488) ✓ EXACT
MT=621 g12:     59k → 463k (target: 413k, 12% off)
```

**Remaining for T03 1e-5** (see `worklog/T03_phase3_handoff.md`):
1. Port Fortran gtff incoherent integration (Lobatto over p' panels, gaminr.f90:1341-1464) to replace 20-pt midpoint. This fixes both the 12% heating error and the -14 line GENDF format diff.
2. Fix pair production MF26 format (NL=1 not NL=5, Fortran line 300)
3. Fix MAT=92 MT=621 (form factor loading)
4. Apply igzero skip to MF26 sections

**Trap 148-158**: See worklog/T03_phase3_handoff.md for full trap descriptions.

**Files changed**: `src/orchestration/modules/gaminr.jl` (heating accumulation, incoherent normalization, section ordering, MF26 writer), `src/orchestration/modules/dtfr.jl` (_dpend state machine)

### Phase 44: T15/T17 JENDL U-238 — CRASH bucket drained, 30M-line errorr blowup fixed

