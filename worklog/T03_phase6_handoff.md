# T03 Pipeline — Phase 6 Handoff

## Date: 2026-04-11

## Goal

Get T03 tape37 (viewr PostScript) passing the official execute.py test at rel_tol=1e-5.

T03 chain: `reconr → gaminr → dtfr → matxsr → viewr`. Only **referenceTape37** (9274 lines) is compared.

## What Was Done This Session

### 1. Fortran gpanel rndoff/delta boundary nudges (FIX)

**Root cause**: Fortran gpanel (gaminr.f90:907-937) applies `rndoff=1.000002` to the lower boundary of the first sub-panel and `delta=0.999995` to the upper boundary of the last sub-panel within each group. Julia used raw PENDF boundaries. For group 1 (1e4–1e5), this caused a 1.5e-5 relative flux error: the delta nudge at 1e5 (where W=1.0) removes ∫W over a 0.5 eV gap, exactly matching the observed 0.5002 absolute difference.

**Fix**: Added `_GPANEL_RNDOFF=1.000002` and `_GPANEL_DELTA=0.999995` constants. Applied nudges to all 5 integration paths:
- `_gaminr_group_average` (MF23 σ averaging)
- `_gaminr_heating_average` (MF23 heating)
- `_gaminr_scatter_matrix` flux computation (MF26)
- `_gaminr_coherent_matrix!` energy quadrature (MF26/MT502)
- `_gaminr_incoherent_matrix!` energy quadrature (MF26/MT504)
- `_gaminr_pair_production_matrix!` (MF26/MT516)

**Impact**: All MF23 flux values now match exactly. ~200 tape33 lines improved.

### 2. Base flux precomputation for threshold MTs (FIX)

**Root cause**: Fortran gpanel computes the same flux for ALL MTs because it walks all breakpoints (PENDF + weight function + form factor). Julia computed flux per-MT from PENDF data. For MT=516 (pair production, threshold at 1.022 MeV), the PENDF starts at the threshold, so the flux integral for group 4 (1–2 MeV) missed the [1.0, 1.022] MeV gap — a 3.1% flux error cascading to 3.2% σ error.

**Fix**: Precompute `base_flux` from MT=501 data (full energy range) and pass it to all MTs via new `base_flux` keyword argument on `_gaminr_group_average`, `_gaminr_heating_average`, and `_gaminr_scatter_matrix`. The pair production function also uses the passed `flux` array (overridden by base_flux) instead of computing its own denominator.

**Impact**: MT=516 errors eliminated (3.2% → ULP). Both MF23 and MF26 sections.

### 3. Free-KN mechanism for incoherent scattering (NEW FEATURE)

**Root cause**: At high energies, the incoherent scattering function S(q,Z) → Z (all momentum transfers are large enough that the form factor saturates). Fortran exploits this: the first (lowest-Z) material is the "reference" — all groups computed via gpanel. Higher-Z materials reuse the reference results scaled by Z for groups above the KN transition energy `E > Z × ekn` (ekn=12.4 keV, gaminr.f90:126).

**Discovery method**: Patched Fortran gpanel caller loop (line 337) with diagnostics. Found MAT=92 (Z=92) only uses gpanel for groups 1–4 (elo < 92×12400 = 1.14 MeV). Groups 5–12 use saved results from MAT=1 (Z=1) scaled by 92. The `zref` variable starts at 101.0 (compiler-initialized), so MAT=1 (Z=1 ≤ 101) becomes the reference with znow=-1 (negative flag). MAT=92 (Z=92 > zref=1) uses the stored `akn` array.

**Key Fortran code (gaminr.f90:311-393)**:
```fortran
! Lines 311-316: set znow = Z, flag reference
if (mfd.eq.26.and.mtd.eq.504) then
   znow=nint(za/1000)
   if (znow.le.zref) then; zref=znow; znow=-znow; endif
endif

! Line 324: free-KN branch for high energies
if (znow.gt.zero.and.elo.ge.znow*ekn) go to 325

! Lines 357-361: retrieve scaled reference
ans(j,1,i) = akn(j,i,ig)
if (i.gt.1) ans(j,1,i) = ans(j,1,i) * znow

! Lines 376-384: save reference (BEFORE dspla normalization)
if (znow.lt.zero) then
   akn(j,i,ig) = ans(j,1,i)
   if (i.gt.1) akn(j,i,ig) = akn(j,i,ig) / abs(znow)
endif
```

**Important subtlety**: The akn is saved BEFORE dspla normalization (raw ans / Z), so when retrieved and scaled by Z_cur, dspla normalizes it correctly. In Julia, the scatter matrix is already normalized by flux, so the save/retrieve handles normalized values differently (divide by Z_ref on save, multiply by Z_cur on retrieve).

**Fix**: Added free-KN mechanism in `_write_gaminr_tape`:
- `akn_zref = Inf` (first material with Z ≤ zref becomes reference)
- Reference material: save all scat_matrix values / Z to `akn` array
- Higher-Z material: for groups where `egg[ig] >= Z * ekn`, replace scat_matrix with akn * Z

**Impact**: MF26/MT504 MAT=92 worst error from 24.7% → 0.054%.

### 4. Pair production heating integral (FIX)

**Root cause**: Julia computed pair production heating as `σ_avg × (E_mid - 2me)` — a midpoint energy approximation. Fortran's gtff MT=516 returns `ff(1,ng)=E` as the heating feed, so gpanel integrates `∫E×σ×W dE / ∫W dE` properly. For wide groups like group 12 (9–20 MeV), the midpoint approximation (E_mid=14.5 MeV) is inaccurate because σ varies significantly across the group.

**Fix**: Added `heat_num` accumulator in `_gaminr_pair_production_matrix!` to compute `∫E×σ×W dE`. Heating = `∫EσW/∫W - 2me × ∫σW/∫W`.

**Impact**: MT=621 heating errors eliminated (MAT=1: 2%→7.5e-6, MAT=92: 17%→gone from significant diffs).

## Current State

```
tape33 (GENDF):   294/604  48.7% byte-identical (was 42.2%)
tape34 (DTF):     294/340  86.5% byte-identical (was 40.6%)
tape36 (plot):   1683/1711 98.4% byte-identical (was 86.3%)
tape37 (PS):     9304/9274 32.4% byte-identical (was 32.4%)
```

### T01 regression: 32750/32962, 212 diffs — UNCHANGED ✓

### Remaining errors by section:

| Section | Worst error | Val diffs | Notes |
|---------|------------|-----------|-------|
| MF23 σ (all MTs) | ~2e-6 | ~50 | ULP from PENDF grid differences |
| MF26/MT502 MAT=1 | 8.1e-6 | 45 | Coherent — very close |
| MF26/MT502 MAT=92 | 1.1e-4 | 57 | Coherent energy quadrature |
| MF26/MT504 MAT=1 | 6.5e-4 | 270 | Incoherent energy quadrature |
| MF26/MT504 MAT=92 | 5.5e-4 | 325 | Incoherent (free-KN groups match, gpanel groups ~0.05%) |
| MT=621 heating | 7.5e-6 | 12 | Essentially matched |
| MF1 header | cosmetic | 1 | PENDF hash identifier |
| tape37 lines | 9304 vs 9274 | +30 | From residual value diffs → different viewr point counts |

### Why tape37 match % didn't improve despite massive tape34/36 improvement:

tape37 has 9304 lines vs 9274 reference (+30 lines). The extra lines come from value diffs in the GENDF (tape33) that cascade through dtfr to different Y-axis ranges → viewr renders different numbers of PostScript path points. With different line counts, many lines can't match positionally even when values are close. The tape37 test at rel_tol=1e-5 should be more forgiving of this.

## Remaining Blockers — Priority Order

### 1. MF26/MT504 incoherent: ~6e-4 systematic (MEDIUM IMPACT)

**Root cause**: Julia's energy quadrature uses PENDF panel boundaries only. Fortran gpanel also splits at weight function breakpoints (from gtflx) and form factor breakpoints (from gtff). For group 12 (9–20 MeV), the wt breakpoint at 1e7 creates a kink in W(E) that Julia's GL-6 quadrature doesn't capture perfectly within a single PENDF panel.

**How to diagnose**: Patch Fortran gpanel to print running `ans` totals after each panel for MAT=1 MT=504 group 12 (31 panels). Compare with Julia's per-panel contributions.

**How to fix**: Add wt breakpoint splitting to the MF26 energy quadrature in `_gaminr_coherent_matrix!` and `_gaminr_incoherent_matrix!`. Split PENDF panels at the 4 wt breakpoints [1e3, 1e5, 1e7, 3e7], same as `_gaminr_group_average` already does.

**Expected impact**: Would reduce 6e-4 to ~1e-5 range. Cascades through tape34/36/37.

### 2. MF26/MT502 coherent MAT=92: ~1.1e-4 (LOW IMPACT)

**Root cause**: Julia uses GL-6 energy quadrature for coherent, while Fortran uses GL-2 (trapezoidal) because gtff MT=502 returns nq=0. For PENDF panels with max ratio ~1.03, GL-2 is adequate. Julia's GL-6 is more accurate mathematically but produces slightly different results.

**How to fix**: Switch coherent energy quadrature from GL-6 to GL-2 (trapezoidal), matching Fortran. This is a simple change: replace the GL-6 loop with endpoint evaluation.

### 3. matxsr transfer matrices (tape35: 131/211)

**Root cause**: `src/orchestration/modules/matxsr.jl` writes diagonal-only P0 transfer matrices. Fortran writes full banded matrices.

**Impact on tape37**: None — dtfr reads GENDF directly. Low priority.

### 4. tape37 line count: 9304 vs 9274 (+30)

Not a separate bug. Residual value diffs → different Y-axis ranges → different point counts.

## Key Traps (NEW this session)

**Trap 170 (CRITICAL)**: Fortran `zref` variable in gaminr starts at 101.0 (NOT 0.0). This is because gfortran initializes uninitialized `real(kr)` locals to 101.0 in this context (likely from stack contents or compiler-specific behavior). The consequence: MAT=1 (Z=1) has `Z <= zref(101)` → TRUE → becomes the reference material with `znow=-1`. If you assume zref=0, the logic inverts and the free-KN mechanism breaks completely.

**Trap 171 (FIXED)**: Pair production heating must use `∫E×σ×W dE / ∫W dE`, NOT `σ_avg × E_mid`. The midpoint approximation has 2–17% error for wide energy groups. The Fortran gtff MT=516 returns `ff(1,ng)=E` as the heating feed function, and gpanel integrates it at each energy quadrature point.

**Trap 172 (FIXED)**: Fortran gpanel computes the SAME flux for ALL MTs regardless of threshold coverage. Julia must use precomputed base flux from MT=501 (total XS, full range). Otherwise, threshold MTs like MT=516 get a smaller flux denominator, producing 3.2% errors.

**Trap 173 (DIAGNOSED)**: The free-KN `akn` array stores raw (BEFORE dspla normalization) values divided by Z_ref. In Julia, the scatter matrix is already normalized by flux. The save/retrieve logic must account for this: save normalized/Z_ref, retrieve normalized×Z_cur.

## How to Run

```bash
# Generate Fortran oracle (if /tmp/t03_fortran/ doesn't exist)
mkdir -p /tmp/t03_fortran && cd /tmp/t03_fortran
cp ~/Projects/NJOY.jl/njoy-reference/tests/03/input .
cp ~/Projects/NJOY.jl/njoy-reference/tests/resources/gam23 tape30
cp ~/Projects/NJOY.jl/njoy-reference/tests/resources/gam27 tape32
~/Projects/NJOY.jl/njoy-reference/build/njoy < input > output 2>&1
# Produces: tape31(2180) tape33(604) tape34(340) tape35(211) tape36(1711) tape37(9274)

# Run Julia pipeline
cd ~/Projects/NJOY.jl
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/03/input"; work_dir="/tmp/t03_test")'

# Quick comparison
python3 -c "
for u in [33, 34, 36, 37]:
    r = open(f'/tmp/t03_fortran/tape{u}').readlines()
    t = open(f'/tmp/t03_test/tape{u}').readlines()
    mn = min(len(r), len(t))
    m = sum(1 for a,b in zip(r[:mn], t[:mn]) if a == b)
    print(f'tape{u}: {m}/{len(r)} ({100*m/len(r):.1f}%) julia={len(t)}')
"

# T01 regression check
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/01/input"; work_dir="/tmp/t01_test")' 2>&1 | tail -3
python3 -c "
ref=open('njoy-reference/tests/01/referenceTape25').readlines()
trial=open('/tmp/t01_test/tape25').readlines()
mn=min(len(ref),len(trial))
m=sum(1 for a,b in zip(ref[:mn],trial[:mn]) if a==b)
print(f'T01 tape25: {m}/{len(ref)} match, {mn-m} diffs')
"
```

## Files Changed

| File | What changed |
|------|-------------|
| `src/orchestration/modules/gaminr.jl` | +170/-29 lines total |
| | `_GPANEL_RNDOFF`, `_GPANEL_DELTA` constants (line ~715) |
| | `_gaminr_group_average` — rndoff/delta nudges + base_flux kwarg |
| | `_gaminr_heating_average` — rndoff/delta nudges + base_flux kwarg |
| | `_gaminr_scatter_matrix` — rndoff/delta on flux + base_flux kwarg |
| | `_gaminr_coherent_matrix!` — rndoff/delta on energy panels |
| | `_gaminr_incoherent_matrix!` — rndoff/delta on energy panels |
| | `_gaminr_pair_production_matrix!` — rndoff/delta + proper ∫EσW heating |
| | `_write_gaminr_tape` — base_flux precomputation + free-KN mechanism |

## Next Steps for Next Agent

### Step 1: Add wt breakpoint splitting to MF26 energy quadrature

The ~6e-4 incoherent error is from PENDF panels spanning the wt breakpoint at 1e7 (inside group 12). Add the same wt breakpoint splitting that `_gaminr_group_average` already does to `_gaminr_coherent_matrix!` and `_gaminr_incoherent_matrix!`.

### Step 2: Switch coherent energy quadrature from GL-6 to GL-2

Fortran uses GL-2 (trapezoidal) for coherent because gtff MT=502 returns nq=0 (→ nq=0+2=2). Julia's GL-6 produces slightly different results. Changing to GL-2 would match Fortran exactly for coherent.

### Step 3: Run execute.py test

After steps 1-2, the tape33 values should be close enough that tape37 line counts match. Run `python3 njoy-reference/tests/execute.py` to check if T03 passes at rel_tol=1e-5.
