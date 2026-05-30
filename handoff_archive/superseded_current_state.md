<!-- [SUPERSEDED 2026-05-30 — see HANDOFF.md Current State] -->

## Current State

### What WORKS — bit-identical

| Test | Material | Formalism | MTs | Grid | Status |
|------|----------|-----------|-----|------|--------|
| Test 84 | H-2 (deuterium) | LRU=0 | 4/4 | 769 pts exact | **BIT-IDENTICAL** |
| Test 01 | C-nat (carbon) | LRU=0 | 29/29 | 1033 pts exact | **BIT-IDENTICAL** |
| Test 02 | Pu-238 | SLBW + URR(mode=11) | 17/17 | 3567 pts exact | **BIT-IDENTICAL** |
| Test 08 | Ni-61 | Reich-Moore (LRF=3) | 18/18 | 4825 pts exact | **BIT-IDENTICAL** |
| Test 18 | Cf-252 | SLBW + URR(mode=12) | 9/9 | exact | **BIT-IDENTICAL** — NEW Phase 14 |
| Test 27 | Pu-239 | Reich-Moore | 49/49 | exact | **BIT-IDENTICAL** — NEW Phase 14 |
| Test 47 | Pu-239 | Reich-Moore | 49/49 | exact | **BIT-IDENTICAL** — NEW Phase 14 |
| Test 45 | B-10 | LRU=0 | 53/53 | 338 pts exact | **BIT-IDENTICAL** — NEW Phase 11 |
| Test 22 | H-1 (para-H₂ 20K) | leapr (phonon + cold-H + trans) | MF7/MT4 | 4636 pts exact | **BIT-IDENTICAL** — NEW 2026-04-23 (leapr wired) |
| Test 61 | H-in-H₂O thermal (acer iopt=7 rename) | ACE edit/rename | tape71 49211 + tape72 1 | exact | **BIT-IDENTICAL (full test)** — NEW Phase 78 |

### ACER charged-particle ACE — tape34 BIT_IDENTICAL (Phase 78)

T50 (α+He-4), T52 (p+H-1), T62 (d+He-3), T53 (d+H-2) all have **ACE
`tape34` + `tape35` (xsdir) BIT_IDENTICAL** as of Phase 78 (ptlegc/coul,
charged NTR>0 reactions, acelcp particle-type blocks). These tests are NOT
yet FULL-test passes — the ONLY remaining failing tape is `tape33`, the
aplots viewr plot tape (a stub: 0 lines vs ref ~5-11k), tracked in bead
NJOY_jl-cnh.4. Implementing that plot tape would flip all four to full
BIT_IDENTICAL/NUMERIC_PASS. T54 (p+H-3) tape34 is DIFFS (7327/7327
structure correct, 6383 match; residual = pre-existing Coulomb-elastic 7e-7
1-ULP noise + 1 triton recoil heating[1] word, bead NJOY_jl-53h).

### Test 80 (H_HF, 343K) — NUMERIC_PASS @1e-5 (post-Phase-57)

91453/91453 line count matches; all MF1/MF7 record layouts correct.
Phase 57 (lat/sc/arat plumbing in `generate_sab`) lifted the dominant
gap: T80 tape24 went **DIFFS 75.7% → NUMERIC_PASS 91406/91453 @1e-5
(99.95%)**. The remaining 47 lines differ by ±1 in 7th sigfig — Phase B
class (SCT-replacement `naint` gating + phonon-loop FP order). See
`worklog/T80_leapr_contin_phase_a_lat_sc.md` for the full Phase A
analysis and Phase B plan. Pre-Phase-57: 76.8% of meaningful S(α,β)
values bit-identical, 77.9% pass 1e-5; the gap was Fortran `contin`'s
missing `sc = therm/tev` rescaling for lat=1.

### In Progress — Test 07 (U-235, SLBW + URR mode=12)

**Test details**: MAT=1395, err=0.005, ENDF file `t511`, resolved range [1, 82] eV (SLBW, LRF=1), unresolved range [82, 25000] eV (LRF=2, mode=12). 27 MTs total.

**Status: 24/27 PERFECT** — grid now matches exactly (2315 pts). Only 3 MTs differ (MT=18/19/102), all by ±1 in last digit at a single energy (82.00001 eV).

**What's been fixed for Test 07 (Phases 8-9):**
1. **Fission doubling (FIXED)**: MT=18 MF3 background was being added alongside MT=19, exactly doubling fission XS. Now MT=18 is filtered from mf3_sections when MT=19 exists, matching Fortran `mtr18` flag (reconr.f90:557-561,1893). MT=18 PENDF output is now a computed redundant sum of MT=19+20+21+38.
2. **GT peak nodes (FIXED)**: `_add_bw_peaks!` now uses GT/2 (total width from ENDF) instead of (|GN|+|GG|+|GF|)/2. Confirmed 4 peak nodes differ for U-235 but sigfig rounding washes them out at ndig=5-6.
3. **merge_background_legacy (FIXED)**: MT=20/21/38 backgrounds no longer accumulated into primary fission channel — they go to `other_bg`, matching Fortran emerge `itype` dispatch.
4. **MF=13 reader (FIXED)**: Removed incorrect histogram interpolation forcing. Fortran `scr(5)=1` at line 1880 is DEAD CODE — overwritten by `tab1io`. MF=13 uses actual ENDF interpolation law (INT=2 for U-235). This eliminated 3 extra histogram-shaded grid points.
5. **MF3Section `mf` field (ADDED)**: lunion skip logic correctly bypasses MT redundancy checks for MF=12/13 sections, matching Fortran (reconr.f90:1881).
6. **Threshold replacement in lunion_grid (FIXED)**: Julia now replaces the first MF3 breakpoint with `thrxx = sigfig(thrx, 7, +1)` when `tab.x[1] < thrxx`, matching Fortran reconr.f90:1925-1936. Also fixes subsequent breakpoints that are now <= the new first energy (lines 1937-1943). Only applied to breakpoint insertion; panel bisection still uses original data to avoid cascading regression.
7. **Pseudo-threshold advancement (FIXED)**: For non-primary MTs, skip leading zero-XS breakpoints where both current and next XS < 1e-30, matching Fortran label 205 (lines 1968-1976). Previously only checked first two breakpoints.
8. **Initial vs mid-data duplicate shading (FIXED)**: Fortran uses different shading for initial discontinuities (label 207, line 1979: `sigfig(er,7,0)`) vs mid-data discontinuities (label 270, line 2029: `sigfig(er,7,-1)`). Julia now distinguishes: `k == start_k` → sigfig(e,7,0), otherwise → sigfig(e,7,-1). This resolved the {1089999, 1090000, 1090001} three-point pattern at the MF=13 discontinuity energy.

**Remaining issues for Test 07 (3 MTs):**
1. **±1 XS at URR boundary (MT=18/19/102)** — At 82.00001 eV (resolved/unresolved boundary): fission 3.380198e1 vs 3.380197e1, capture 1.981499e1 vs 1.981500e1. Root cause: floating-point accumulation differences in `_gnrl` Gauss-Laguerre quadrature (100 terms accumulated). The raw values are within 1 ULP of the 7-sigfig rounding boundary. The two-step interference correction split (matching Fortran reconr.f90:4271-4272) was applied but didn't resolve it. This is a cross-compiler precision issue, not a logic bug.

**Fixed in Phase 9 (no longer issues):**
- **MT=1 ±4 at 12.2 MeV (FIXED)**: Was caused by missing below-threshold skip in `merge_background_legacy`. MT=17 was interpolating a spurious ~3.98e-6 contribution at E=1.219910e7 (below its threshold). Added guard matching Fortran emerge line 4792: `if (thrx - e) > 1e-10 * thrx → skip`.

### In Progress — Test 34 (Pu-240, Reich-Moore + URR mode=11, LSSF=1)

**Test details**: MAT=9440, err=0.001, ENDF file `n-094_Pu_240-ENDF8.0.endf`, resolved range [1e-5, 5700] eV (RM, LRF=3, 437 fissile l=0 + 121 non-fissile l=1 resonances), unresolved range above 5700 eV.

**Status: 52/53 PERFECT** — 3 remaining ±1 diffs in MT=102 (capture) at E=630.04, 2089.07, 4526.46 eV.

**What's been fixed for Test 34 (Phase 12):**
1. **Coincidence shading bug (FIXED)**: `lunion_grid` used `abs(grid[gi] - e) <= 1e-8*e` (direct comparison) for cross-section coincidence check at label 220. Fortran (reconr.f90:1996) compares `(er - sigfig(|eg|,7,-1)) <= 1e-8*er` (against shaded-down grid point). For E=42980 (shared by MT=2, MT=102, MF=12/MT=4), Julia falsely detected coincidence and created spurious grid point 42980.01, cascading 389+ line shifts across 7 MTs. This bug was introduced by the Phase 11 coincidence shading commit.
2. **Frobenius-Schur matrix inversion (FIXED)**: `cross_section_rm` used `inv(SMatrix{3,3,ComplexF64})` (generic complex inverse). Fortran `csrmat` (reconr.f90:3503-3607) uses Frobenius-Schur method via `frobns`→`thrinv`→`abcmat`. Implemented exact Fortran algorithm (`_frobns`, `_thrinv!`, `_abcmat` in `reich_moore.jl`). Both compute `(I+R+iS)^{-1}` but with different intermediate FP rounding. Resolved 2 of 5 ±1 capture diffs.

**Remaining issues for Test 34 (3 MTs in MT=102):**
1. **±1 capture XS at 3 energies** — **CONFIRMED IRREDUCIBLE via gdb (Phase 14)**: Fortran values traced with diagnostic prints in csrmat. At E=630.04: Fortran cap=0.096262384997, Julia cap=0.096262384999 (diff=2.3e-12). At E=2089.07: diff=1.0e-11. At E=4526.46: diff=-1.3e-13. All within 1e-11 of the 0.5 boundary at 7 sigfigs. The same Frobenius-Schur algorithm is used in both; the difference is purely from IEEE 754 non-associativity across 437 resonance accumulations. No fix possible without matching exact intermediate rounding order.

### Not Yet Attempted

| Test | Material | Formalism | Notes |
|------|----------|-----------|-------|
| Test 15-17 | U-238 JENDL | Likely Reich-Moore | ENDF file not yet located in resources/ |

### Formalisms status

| Formalism | Reader | Evaluator | Tested |
|-----------|--------|-----------|--------|
| LRU=0 (no resonances) | N/A | N/A | Tests 84, 01 BIT-IDENTICAL |
| SLBW (LRF=1) | `_read_bw_params` | `cross_section_slbw` | Test 02 BIT-IDENTICAL |
| MLBW (LRF=2) | `_read_bw_params` | `cross_section_mlbw` | Not tested |
| Reich-Moore (LRF=3) | `_read_rm_params` | `cross_section_rm` | Test 08 BIT-IDENTICAL |
| SAMMY/RML (LRF=7) | `_read_sammy_params` | `build_rml_evaluator` | Not tested |
| URR mode=11 (LRU=2, LRF≤1, LFW=1) | `_read_urr_lfw1` | `_csunr1`/`_gnrl` | Test 02 BIT-IDENTICAL |
| URR mode=12 (LRU=2, LRF=2) | `_read_urr_lrf2` | `_csunr2`/`_gnrl` | Test 07 (runs, 0/27, grid diffs) |

---


---

## Immediate Next Steps — PRIORITY ORDER

### 0. CHECK: Run the pipeline FIRST to establish baseline

Before making ANY changes:
```bash
cd ~/Projects/NJOY.jl
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/t01_pipeline.jl
```

**Expected results (Phase 34, as of 2026-03-31)**:
```
tape25: 32962 (ref: 32962)       ← STRUCTURAL: EXACT MATCH ✓
MATCH: 41 / 41 sections          ← ALL sections match ✓
MT=102: DATA 307/307 PERFECT     hdr=3/3
MT=  4: DATA 135/135 PERFECT     hdr=3/3
MT= 51: DATA 135/135 PERFECT     hdr=3/3
MT= 91: DATA 88/88 PERFECT       hdr=3/3
MT=103: DATA 26/26 PERFECT       hdr=3/3
MT=230: DATA 191/191 PERFECT     hdr=3/3
MT=229: DATA 189/191 (99.0%)     hdr=3/3  ← 2 diffs from coh window FP at grid edge
MT=  2: DATA 306/307 (99.7%)     hdr=3/3  ← 1 line ±1 ULP sigma1
MT=221: DATA 47/49 (95.9%)       hdr=3/3  ← 2 lines emax boundary
MT=  1: DATA 240/307 (78.2%)     hdr=3/3  ← ±1 ULP sigma1 broadening
MT=301: DATA 198/307 (64.5%)     hdr=3/3  ← sub-ULP, cascades from MT=1
MT=444: DATA 119/307 (38.8%)     hdr=3/3  ← ±1 ULP cascaded through damage×sigma
Total data: 2017 / 2386 (84.5%)
```

**Tolerance test**:
```
rel_tol=1e-9: 355 lines fail (1.1%)
rel_tol=1e-7: 355 lines fail (1.1%)  ← same as 1e-9 (all diffs > 1e-7)
rel_tol=1e-5:   0 lines fail (0.0%)  ← PASSES at 1e-5 ✓
rel_tol=1e-4:   0 lines fail (0.0%)  ← PASSES at 1e-4 ✓
rel_tol=1e-3:   0 lines fail (0.0%)  ← PASSES at 1e-3 ✓
```

**T01 passes the official test at 1e-5 tolerance.** All remaining 355 diffs are below 1e-5 relative.

**How to run the tolerance test** (simulates execute.py):
```bash
python3 -c "
import re; from math import isclose
fp=re.compile(r'(([-+]?\d+)(\.\d+)?(([eE])?[-+]\d+)?)',re.VERBOSE)
def mf(D):
    try: return float(D[0])
    except: return float('{}{}E{}'.format(D[1],D[2],D[3])) if not D[-1] else None
ref=open('njoy-reference/tests/01/referenceTape25').readlines()
trial=open('/tmp/t01_tape25.pendf').readlines()
for rtol in [1e-9,1e-7,1e-5,1e-4,1e-3]:
    f=0
    for r,t in zip(ref,trial):
        if r==t: continue
        rF=fp.findall(r); tF=fp.findall(t)
        if not rF:
            if fp.sub('',r)!=fp.sub('',t): f+=1
            continue
        if len(rF)!=len(tF): f+=1; continue
        ok=True
        for a,b in zip(rF,tF):
            rv=mf(a); tv=mf(b)
            if rv is None or tv is None: continue
            if not isclose(rv,tv,rel_tol=rtol,abs_tol=1e-10): ok=False; break
        if not ok: f+=1
    print(f'rel_tol={rtol:.0e}: {f} lines fail ({f/len(ref)*100:.1f}%)')
"
```

### YOUR GOAL: Reduce T01 physics failures

T01 currently has **98 physics failures at 1e-3** and **226 at 1e-4**. MT=301 and MT=444 worst errors are both under 0.3% — the HEATR grind is essentially done. The remaining targets:

**Priority 1 — Fix 487 MF12/MF13 sequence number diffs (cosmetic, biggest line count)**:
The MF12/MT=102 and MF13/MT=51 passthrough reads from `after_reconr.pendf` oracle. The data (cols 1-75) is correct but sequence numbers (cols 76-80) differ. Fix: read from the original ENDF tape (`tape20`) instead, or reset sequence counters to start from 1 per section. This is in `t01_pipeline.jl` around line 284 where `after_reconr.pendf` lines are collected.

**Priority 2 — Fix 3 MF1/MT=451 directory NC values (mechanical)**:
Fortran `tpend` uses truncated integer division `div(N,6)` for line count in directory. Julia uses `ceil`. Fix: in `_write_full_mf1` (pendf_writer.jl), for thermr and MF6 sections, compute NC = HEAD(1) + CONT(1) + div(2*NR+2*NP, 6) for TAB1, using integer division not ceiling. Check divisibility: only subtract 1 when `mod(data_values, 6) != 0`.

**Priority 3 — Fix MF2/MT=151 EH value (trivial)**:
One line diff. Check what EH value the Fortran writes vs Julia.

**Priority 4 — Grind MT=444 remaining 30 diffs at 1e-3**:
All below 0.28%. Root cause: disbar stepping state machine doesn't perfectly reproduce Fortran bracket sequence. The `build_disbar_damage_vector` uses `step=1.1` and MF4 energy clamping, but the exact sequence of (el, daml, en, damn) brackets may diverge from Fortran when multiple query points fall in the same bracket. Use gdb to trace Fortran disbar at specific energies where Julia differs.

**Priority 5 — Grind MF6/MT=229 remaining 26 diffs at 1e-3**:
sigl FP accumulation order. Same class as Phase 27. Use gdb on Fortran sigl to trace intermediate values at specific IE/E' pairs.

**Priority 6 — Grind MT=301 remaining 38 diffs at 1e-3**:
All cascade from sigma1 broadening ±1 ULP in MT=1. These are irreducible without matching the exact FP accumulation order in the sigma1 Doppler broadening kernel. Same class as reconr T34 diffs.

**Priority 7 — Generate oracle caches for untested RECONR tests**:
31 tests RAN_OK without oracles. Many share materials with BIT-IDENTICAL tests.

**Priority 8 — Grind BROADR to bit-identical on more tests**.

### Remaining diagnosed failures — exact breakdown

**MF6/MT=229 near-zero kernel cosines (2 lines at 1e-7)**:
At IE=44 (incident E ≈ 0.01 eV), secondary energy E' ≈ 0.593 eV has sigma = 3.8e-10 (near the sigmin=1e-10 floor). The equi-probable cosines are completely different:

```
Julia:  -3.400e-1  5.744e-2  0.1821  0.3065  0.4306  0.5608  0.7529  9.216e-3
Fortran: -2.674e-1 -7.267e-2  5.862e-2  0.1833  0.3077  0.4317  0.5622  0.7560
```

**Root cause**: `sigl_equiprobable` integrates σ(μ) via adaptive linearization, then inverts the CDF to get equi-probable cosines. When sigma ≈ 0, the CDF is nearly flat, making the inversion numerically unstable. Any tiny FP difference in the kernel evaluation (even 1e-15) produces completely different cosines.

**How to fix**: Match the Fortran `sigl` computation at this exact E/E' pair. The Fortran `sig` function (thermr.f90:2514-2566) zeros kernel values below sigmin=1e-10. If Julia's `sab_kernel` returns a value epsilon above sigmin while Fortran returns epsilon below (or vice versa), the CDF shapes diverge. Use gdb to trace:
```bash
# Patch thermr.f90 sigl at line ~2660 with filter:
#   if (ie.eq.44.and.j.eq.???) write(*,*) 'SIGL',e,ep,mu,sig_val
# Compare with Julia's sigl_equiprobable at same E, E'
```

This is a borderline convergence decision. The sigma is so small that the cosines have no physical impact, but matching them requires exact FP agreement in the kernel at this E'.

**Category C: 4 MF6/MT=229 ±1 ULP cosines (lines 31109, 31841, 32199, 32201)**

At high incident energies (IEs 85-93, E > 0.7 eV), individual cosine values differ by ±1 in the last significant digit:

```
Line 31109: Julia -0.3598930 vs Fortran -0.3598931 (diff = 1e-7)
Line 31841: Julia -0.3692905 vs Fortran -0.3692906 (diff = 1e-7)
Line 32199: Julia -0.3471784 vs Fortran -0.3471790 (diff = 6e-7)
Line 32201: Julia  0.25425679 vs Fortran 0.25425661 (diff = 2e-7)
```

**Root cause**: `sigl_equiprobable` phases 1 and 2 accumulate σ(μ) via adaptive integration. The Julia and Fortran accumulate in the same order but IEEE 754 non-associativity in the intermediate FP operations produces ~1e-7 relative differences in the CDF, which shift the equi-probable cosine boundaries.

**How to fix**: Use gdb to trace the exact intermediate values in `sigl` at one of these IEs. The approach that has worked throughout this project:

1. Patch Fortran `sigl` (thermr.f90:2660-2872) with `write(*,*)` at the phase 1 linearization loop to print each accepted (μ, σ) pair
2. Print the phase 2 CDF inversion results
3. Compare side-by-side with Julia's `sigl_equiprobable` intermediate values
4. Find the first intermediate value that diverges and trace to a FP operation ordering difference

The most likely culprits:
- **seep computation**: Fortran line 2700 computes `seep = 1/sqrt(e*ep)` as a single divide. Julia might compute differently.
- **peak location**: The x_peak formula uses AWR, kT, and seep. Any intermediate rounding difference shifts the peak, changing the adaptive linearization.
- **CDF moment**: Uses `third = 0.333333333` (Fortran truncated constant). Matched in Phase 25 but verify.

**Expected difficulty**: Lines 31109 and 31841 have ~1e-7 relative error — right at the boundary. Lines 32199 and 32201 have larger errors (~6e-7) — these may have a systematic cause.

### 7. Grind remaining RECONR tests + BROADR

- **19 BIT-IDENTICAL RECONR tests** (T01-03,08-13,18-19,25-27,30,45,47,55,84)
- **31 tests RAN_OK without oracle** — many share materials with BIT-IDENTICAL tests
- BROADR is implemented (broadn_grid matches 919/919 points for T01)
- Apply the Grind Method to each test

---


---

## Date: 2026-04-18

**Goal**: Drain the last 2 CRASHes in the 84-test sweep (T15, T17, both
JENDL U-238 `BoundsError`).

**Fresh full sweep (93.9 min, ~42% faster than Apr-14's 163 min post-P10)**:

| Status | Post-P10 (Apr-14) | Apr-18 baseline | **Post-P44** | Δ |
|--------|-------------------|-----------------|--------------|----|
| BIT_IDENTICAL | 1 | 1 | 1 | = |
| NUMERIC_PASS  | 1 | 1 | 1 | = |
| DIFFS         | 48 | 62 | **64** | **+16** |
| MISSING_TAPE  | 17 | 17 | 17 | = |
| NO_REFERENCE  | 1 | 1 | 1 | = |
| **CRASH**     | **16** | **2** | **0** | **−16** |

**Zero structural crashes across all 84 tests.** The Apr-14 → Apr-18
baseline delta (+14 DIFFS) is the cumulative effect of T11-T14 stubs
and fallbacks finally registering in the sweep framework.

**Fix 1 — MT=455 LIST-skip in groupr `_read_nubar`** (commit `60f5bdc`)

Root cause: MF1/MT=455 (delayed ν̄) with LNU=2 has a LIST of NNF
precursor decay constants between the HEAD and the TAB1 (ENDF-6 §1.5;
Fortran `groupr.f90:6472-6483` label 110). `_read_nubar` read the TAB1
immediately, returning empty `(energies, values)`. Downstream
`_groupr_nubar_records` → `_interp_linlin` crashed on
`x_data[1]`. Fixed: `mt == 455 && read_list(io)` before `read_tab1`.
Gated — MT=452/456 unchanged.

The HANDOFF's prior description of this as "INT=0 fix exposing a
downstream empty-vector indexing bug" was incomplete. The INT=0 fix
(Phase 10) let T15/T17 run past the reader, and the real downstream
bug turned out to be in groupr, not some unnamed "downstream" component.

**Fix 2 — errorr output grid by ign** (commit `c5730f3`)

Root cause: `errorr_module` used the union of user_egn + all MFcov
breakpoints as the **output** group structure for every `ign`. For
T15/T17 (`ign == 3`, LANL-30, 30 groups expected) this produced 2305
groups → 2305×2305 covariance writes × 36 reactions = **30 253 210
lines, 2.5 GB tape26**.

Fortran `egngpn` (`errorr.f90:9716`): the union replaces egn only
when `ign == -1`. For `ign == 1`, egn = user_egn. For `ign >= 2`, egn
is the library structure (LANL-30 etc. via `gengpn`,
`groupr.f90:1557`).

Fix: dispatch grid construction on `params.ign` via new helper
`_errorr_output_grid`. Union path preserved only for `ign == -1`.

Result: T15 tape26 **30 253 210 → 7953 lines (3800× reduction)**. T17
similar. T04 (ign=-1 + ign=1) zero regression: 81/82 NUMERIC_PASS @
1e-7, 56/74 @ 1e-5, 107/119 DIFFS — identical to pre-fix baseline.

**Still left on T15/T17** — all three items below landed in
Phases 45-47 (2026-04-19). Summary:

1. ~~**Missing MF3 sections**~~ → **FIXED Phase 46**:
   `_errorr_read_gendf_xs` reads per-group sigma from the input GENDF
   when `npend==0 && abs(ngout)>0`. T15 tape26 MF3 0 → 36 MTs.
2. ~~**MF33 over-expansion**~~ → **FIXED Phase 47**: sparse per-row
   emission in `_write_mfcov_rows` + NC-aware pair synthesis.
   T15 tape26 8205 → 1859 lines (vs ref 5958).
3. ~~**Groupr `3 /` auto-expand**~~ → **FIXED Phase 45**: parser emits
   `-1000` sentinel on bare cards; `groupr_module` expands against
   PENDF MF=3 keys via `_nextr_filter`. T15 tape91 3 → 39 MTs.

Remaining T15 tape26 gap (~100 lines total after the three fixes):
- −2300 lines on MT=2/MT=4 self-cov from NC-derived cross-MT
  covariance not yet expanded (`NJOY.jl-km1`).
- Content drift in several MT self-cov matrices (`NJOY.jl-f8k`).

**Traps (NEW)**

- **Trap (MT=455 LIST)**: MF1/MT=455 LNU=2 has a LIST of NNF
  precursor decay constants between the HEAD and the TAB1. MT=452
  and MT=456 don't. Any reader handling all three must special-case
  MT=455.
- **Trap (errorr output grid by ign)**: `ign == -1` → union; `ign == 1`
  → user_egn; `ign >= 2` → library structure. The MFcov breakpoint
  union is NEVER the output grid for `ign >= 1`. Using it inflates
  LANL-30's 30-group write to 2305-group (~4000× file size).
- **Trap (groupr `3 /` auto-expand — FIXED Phase 45)**: groupr card
  `<mfd> /` with no mtd is a Fortran sentinel (`mtdp=-1000`) meaning
  "auto-process all MF=3 MTs passing `mt<=200 OR 203<=mt<=207 OR
  mt>300`" (groupr.f90:622 → nextr at 1087-1123). Thermal 201-202 and
  derived 208-300 must be named explicitly. Julia parser now emits the
  sentinel; `_groupr_expand_auto` in `groupr_module` expands against
  PENDF `keys(mf3)`. T15 tape91 MF=3: **3 → 39 MTs** (ref 40).

Worklogs: `worklog/T15_T17_mt455_crash.md`,
`worklog/T15_T17_errorr_size.md`,
`worklog/T15_groupr_auto_expand.md`.

**Files changed**: `src/orchestration/modules/groupr.jl` (1-line LIST
skip for MT=455), `src/orchestration/modules/errorr.jl` (ign dispatch
+ `_errorr_output_grid` helper).

### Phase 45: groupr `3 /` auto-expand — 3 → 39 MTs on T15 tape91

## Date: 2026-04-19

**Fix**: Parser distinguishes "missing mtd" (emit `-1000` sentinel)
from explicit `0` by token count. `groupr_module` expands sentinels
against the PENDF MF=3 dict, applying the Fortran nextr filter. Scoped
to `mfd=3`; other mfd auto-paths need `conver` port (deferred).

T04 regression: zero (tape23 BIT_IDENTICAL, tape24 NUMERIC_PASS
56/74, tape25 DIFFS 108/119 — identical to Phase-44 baseline). T01/T02
groupr decks use only explicit-mtd cards; untouched.

**Orthogonal follow-ups surfaced** (filed as beads):
- `NJOY.jl-cdy`: groupr MT=251/252 need MF=4/6 derivation.
- `NJOY.jl-5oi`: groupr must skip MTs with all-zero group-averaged XS
  (U-238 MT=37 threshold 17.82 MeV > LANL-30 top 17.0 MeV).

### Phase 46: errorr GENDF MF3 readback (same-ign `colaps`) — T15 tape26 MF3 0 → 36

## Date: 2026-04-19

**Fix**: `_errorr_read_gendf_xs` walks MF=3 on the input GENDF and
populates `group_xs` when `npend==0 && abs(ngout)>0 && iread==0`.
Position 2 of the per-group LIST body (sigma for standard MTs, nubar
for 452/455/456) gets extracted directly. Scoped to same-ign only
(groupr ign == errorr ign). Cross-ign flux-weighted collapse deferred.

T15 tape25 MF3: 0 → 3 MTs (nubar). T15 tape26 MF3: 0 → 36 MTs; line
count 7953 → 8205 (reference 5958; remaining +2247 is the known MF33
over-expansion).

T04 zero regression (tape23/24/25 identical to Phase-45 baseline).

Worklog: `worklog/T15_errorr_gendf_readback.md`.

### Phase 47: errorr MF33 sparse per-row emission + NC-aware stubs — T15 tape26 8205 → 1859

## Date: 2026-04-19

**Fix**: `_write_mfcov_rows` helper emits each row as a LIST with
`(L1=count, L2=ig2lo, N1=count, N2=ig)` where `[ig2lo, ng2]` is the
nonzero column range (threshold `|v|>1e-20`). All-zero rows dropped
except the last (stub at `ig2lo=ng2=ig`), matching Fortran covout at
errorr.f90:7530-7605. NC-aware pair synthesis: MTs with any NC>0
sub-section get cross-pair stubs for `mt2>mt` in reaction_mts —
preserves T04's (18,102) NC-derived pair without over-emitting for
T15 MTs with pure self-cov NI.

T15 tape26: **8205 → 1859 lines** (reference 5958). Remaining gap
(−2300 from MT=2/4 self-cov) is NJOY.jl-km1 (NC expansion for cross-MT
values). T04 tape23 **BIT_IDENTICAL 82/82 preserved**; phase-45/46
tests no regression.

Worklog: `worklog/T15_errorr_mf33_sparse.md`.

### Phase 48: errorr MF33 NC-block expansion (LTY=0) — T15 tape26 1859 → 4178

## Date: 2026-04-20

**Fix**: Implement Fortran NJOY's NC-derived covariance via the
collapsed `Σᵢ cᵢ² Cov(refᵢ, refᵢ)` (self) and `cⱼ Cov(refⱼ, refⱼ)`
(cross-pair) formulas — valid where input cross-MT covariances are
zero (true for U-238 JENDL). New `NCSubSubsection` struct, new
`_read_nc_subsection` parser (CONT + LIST = 2 records, not 1), new
`_expand_nc_blocks!` post-processing pass between read loop and
writer. The diagonal `(mt, mt)` entry is **replaced**, not summed,
matching Fortran `akxy[derived,derived,k]=0` (errorr.f90:1475).

**Results** (T15 tape26):

| Quantity | Pre-Phase-48 | Post-Phase-48 | Reference |
|----------|--------------|---------------|-----------|
| MF33 MT=2 lines | 106 | **1429** | 1400 |
| MF33 MT=4 lines | 103 | **1099** | 1106 |
| MF33 total | 1556 | **3875** | 5655 |
| Tape26 total | 1859 | **4178** | 5958 |

T04 tape23 unchanged (NUMERIC_PASS 81/82 — already at this status on
master; HANDOFF Phase 47's "BIT_IDENTICAL" claim was inaccurate);
phase-46/47 tests no regression (one stale `total < 4000` cap relaxed
to `< 5000` to admit NC's legitimate +2300-line growth).

**Limitation v1**: cross-pairs where both endpoints are NC-derived
(T15 Cov(2, 4) = `−Σ Cov(iy, iy)` over `iy ∈ MT2_refs ∩ MT4_refs`)
not yet computed. ~40 lines per pair; one such on T15.
LTY=1/2/3 NC sub-subsections (T04 U-235 cross-material standards)
still emit zero stubs — port `stand` (errorr.f90 ~7800) when needed.

Worklog: `worklog/T15_T17_errorr_nc_expansion.md`.

**Immediate next-step candidates**: see the authoritative "Open Work"
section at the top of this file. P1 items closest to the current
context are:

1. **NC-block expansion v2** — double-NC-derived cross-pairs + LTY=1/2/3
   standards. Closes T15 tape26 MT=2/MT=4 residual (~+40 lines);
   unlocks T04 U-235 cross-material covariance.
2. **Covcal content drift** — Julia covariance matrix values/extents
   differ from Fortran for non-NC-derived MTs. The −1780-line
   residual on T15 tape26 after Phase 48 is concentrated here.
