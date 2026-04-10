# T03 Pipeline — Phase 5 Handoff

## Date: 2026-04-10

## Goal

Get T03 tape37 (viewr PostScript) passing the official execute.py test at rel_tol=1e-5.

T03 chain: `reconr → gaminr → dtfr → matxsr → viewr`. Only **referenceTape37** (9274 lines) is compared.

## What Was Done This Session

### 1. Coherent angular integration — x-space GL quadrature (CRITICAL FIX)

**Root cause**: Julia's coherent angular integration used 16-point midpoint rule in cos(θ). The Fortran `gtff` MT=502 (gaminr.f90:1259-1334) integrates in **momentum-transfer x-space** with GL-6/GL-10 quadrature, panels at form factor table breakpoints.

**Key difference**: The Fortran integrand is `(1+u²) × |F(x)|² × P_l(u) × |du/dx|` where `u = 1 - x²/(c1·E)²`, integrated over x from 0 to `xlim = √2·c1·E`. The Jacobian `fact × xnow` maps from x-space to the differential cross section. Panels are the MF27 form factor tabulation intervals.

**Fix**: Rewrote `_coherent_feed_at_energy!` to faithfully port the Fortran. GL-6 for nl≤4, GL-10 for nl>4. Walks through form factor breakpoints as panel boundaries. Returns normalized Legendre moments.

**Impact**: MF26/MT502 worst error from 55%→0.006% (MAT=1), 73%→0.01% (MAT=92).

### 2. Incoherent energy quadrature — GL-6 (STRUCTURAL FIX)

**Root cause**: Julia evaluated the angular distribution (p-space GL) at ONE midpoint energy per PENDF panel. The Fortran `gpanel` evaluates `gtff` at nq=6 GL quadrature points per sub-panel.

**Fix**: Extracted `_incoherent_feed_at_energy!` — the p-space GL integration at a single energy (faithful port of gtff MT=504, verified identical to 8+ digits against Fortran diagnostic). Then added GL-6 energy quadrature in the outer PENDF-panel loop, calling this function at each quadrature point.

**Verification**: Patched Fortran with diagnostic print at gtff normalization point. Julia and Fortran produce identical siginc, ebar, and normalized ff values at every test energy.

**Impact**: MF26/MT504 worst error from 320%→0.5% (MAT=1).

### 3. Reaction rate linear interpolation — σ×W product (MATCH FIX)

**Root cause**: Julia computed `σ_linear(eq) × W_exact(eq)` at GL points. Fortran gpanel linearly interpolates the PRODUCT `σ×W`:
```fortran
a = sig(iz,1) * flux(iz,il)     ! upper boundary σ × W
b = slst(iz,1) * flst(iz,il)    ! lower boundary σ × W
rr = b + (a-b) * t1             ! linear interp of product
```

**Verification**: Patched Fortran gpanel to print sub-panel boundaries. Confirmed Julia and Fortran use the SAME 17 PENDF panels for MAT=92 MT=504 group 5 (identical elo/ehi values). The PENDF grid has max_ratio ~1.03, finer than the 1.05 step cap from gtflx.

**Fix**: Changed Julia GL loop to compute `rr_lo = sa*wa`, `rr_hi = sb*wb`, then `rr = rr_lo + (rr_hi - rr_lo)*t1` at each GL point.

### 4. Diagnostics methodology established

**gtff diagnostic**: Patch after line 1464 (after normalization), print E/siginc/ebar/ff values for specific MAT/MT. ALWAYS use `open(99,file=...,position='append')`, declare all new variables, full rebuild with `rm -rf CMakeFiles libnjoy.so`.

**gpanel diagnostic**: Patch the caller loop (around line 337), save `elo_save=elo` before gpanel call, print panel boundaries.

**Key finding**: PENDF max_ratio 1.03 < 1.05 step cap → PENDF grid determines sub-panels, NOT the 1.05 step. The 1.05 step from gtflx never fires for this test case.

## Current State

```
tape33 (GENDF):  604/604   42.2% byte-identical
tape34 (DTF):    340/340   40.6% byte-identical (was 34.7%)
tape35 (MATXS):  131/211   diagonal-only
tape36 (plot):  1711/1711  86.3% byte-identical (was 79.8%)
tape37 (PS):    9304/9274  32.4% byte-identical (was 19.1%)
```

### Remaining errors by section:

| Section | Worst error | Notes |
|---------|------------|-------|
| MF23 flux | 6e-5 | Systematic across all sections, root cause TBD |
| MF26/MT502 MAT=1 | 6e-5 | Essentially exact |
| MF26/MT502 MAT=92 | 1.1e-4 | Very good |
| MF26/MT504 MAT=1 | 5.5e-3 (0.55%) | From energy quadrature residual |
| MF26/MT504 MAT=92 | 2.4e-1 (24%) | On tiny values only; systematic ~1.4% |
| MF26/MT516 | 3.2e-2 (3.2%) | Pair production — no GL energy quad yet |
| MT=621 MAT=1 | 2.0e-2 (2.0%) | Cascading from MF26/MT504 |
| MT=621 MAT=92 | 1.7e-1 (17%) | Cascading from MF26/MT504 |

### Regression checks:
- T01 tape25: 32750/32962, 212 diffs — UNCHANGED ✓

## Remaining Blockers — Priority Order

### 1. MF23 flux: ~6e-5 systematic error

**Root cause**: TBD. Both Fortran and Julia use trapezoidal on PENDF panels. Possible sources:
- Fortran `rndoff=1.000002` nudge at panel start (elo *= rndoff)
- Fortran `delta=0.999995` backing off from upper boundary (ehigh = ehi*delta)
- Different handling of panel boundary evaluation

**How to diagnose**: Patch Fortran gpanel to print flux accumulation per panel for a specific group. Compare with Julia's trapezoidal flux computation.

**Impact**: Fixes ~200 of the non-matching tape33 lines (every MF23 data line has flux as first value). Cascades to tape34/36/37.

### 2. MF26/MT504 MAT=92 incoherent: ~1.4% systematic error

**Root cause**: Despite using the same PENDF panels and matching angular distribution to 8+ digits, the integrated scatter matrix values differ by ~1.4% at high energies. The rr linear interpolation and flux consistency match Fortran's approach.

**Possible causes**:
- The `rndoff` nudge changes the effective panel boundaries slightly
- Different cumulative floating-point rounding across 17 panels
- The form factor table interpolation (_interp_ff stateless binary search vs Fortran terpa stateful forward-walk) could produce slightly different values at the PANEL BOUNDARIES where xnow is updated

**How to diagnose**: Patch Fortran to print the accumulated ans values after each panel for one specific source group. Compare running totals.

### 3. MF26/MT516 pair production: 3.2% error

**Root cause**: The pair production XS in MF23/MT=516 group 4 (1-2 MeV) has 3.2% error. This group spans the 1.022 MeV threshold. The per-energy threshold handling clips the integration start differently.

**How to diagnose**: Compare MF23/MT=516 values line by line. Check if the Fortran produces the same sigma for threshold group.

### 4. matxsr transfer matrices (tape35: 131/211)

**Root cause**: Diagonal-only output. Needs full banded matrices from GENDF scatter data. This doesn't affect tape37.

### 5. tape37 line count: 9304 vs 9274 (+30)

**Root cause**: The remaining value diffs cause slightly different Y-axis ranges in dtfr plot commands, which causes viewr to render more data points. Will converge as tape33 values improve.

## Key Traps (NEW this session)

**Trap 164 (FIXED)**: Coherent angular integration must be in x-space (momentum transfer), NOT cos(θ). The Fortran integrand is `(1+u²)|F(x)|² P_l(u) × fact × x` where `u = 1-x²/(c1·E)²`. Integration panels are the MF27 form factor tabulation intervals. Using cos(θ) as the integration variable produces 20-55% errors because F(x(μ)) varies rapidly near F(q) breakpoints.

**Trap 165 (DIAGNOSED)**: Fortran gpanel linearly interpolates the PRODUCT σ×W (reaction rate), not σ and W separately. `rr = b + (a-b)*t1` where `a = sig_hi*W_hi`, `b = sig_lo*W_lo`. This matters because W is log-log interpolated (iwt=3) but σ is lin-lin (PENDF). Within a 1.03-ratio PENDF panel, the product σ×W varies less than 3%, making linear interpolation adequate.

**Trap 166 (VERIFIED)**: The PENDF grid for this test case has max_ratio ~1.03, which is FINER than the Fortran's 1.05 step cap from gtflx. So the 1.05 step never fires — the PENDF breakpoints from gtsig always come first. Both Fortran and Julia use the same panel boundaries (verified by diagnostic).

**Trap 167 (DIAGNOSED)**: The ~6e-5 flux error is systematic and grows from 1.5e-5 (group 1) to 6e-5 (group 11). Not from the 1.05 step or wt breakpoints. Likely from Fortran's rndoff/delta nudges which slightly shift panel boundaries.

## How to Run

```bash
# Generate Fortran oracle (if /tmp/t03_fortran/ doesn't exist)
mkdir -p /tmp/t03_fortran && cd /tmp/t03_fortran
cp ~/Projects/NJOY.jl/njoy-reference/tests/03/input .
cp ~/Projects/NJOY.jl/njoy-reference/tests/resources/gam23 tape30
cp ~/Projects/NJOY.jl/njoy-reference/tests/resources/gam27 tape32
~/Projects/NJOY.jl/njoy-reference/build/njoy < input > output 2>&1

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
```

## How to Use gdb Diagnostics

**CRITICAL**: Always do FULL clean rebuild:
```bash
cd njoy-reference/build
rm -rf CMakeFiles libnjoy.so
cmake ..
cmake --build . --target njoy
```

**gtff diagnostic** — print angular distribution at specific energies:
```fortran
! After line 1464 (after ff normalization)
if (matd.eq.92.and.mtd.eq.504) then
   open(99,file='/tmp/gtff_diag.txt',position='append')
   write(99,'(a,1p,e15.8,a,e15.8,a,e15.8)') 'E=',e,' siginc=',siginc,' ebar=',ebar
   write(99,'(a,5(1p,e15.8))') '  ff(1:5,1)=',(ff(il,1),il=1,min(nl,5))
   close(99)
endif
```

**gpanel diagnostic** — print sub-panel boundaries:
```fortran
! In the caller loop (around line 337), DECLARE: integer::npanel_diag; real(kr)::elo_save
! Before gpanel call: elo_save=elo; npanel_diag=npanel_diag+1
! After gpanel call:
if (matd.eq.92.and.mtd.eq.504.and.ig.eq.5) then
   open(99,file='/tmp/gpanel_diag.txt',position='append')
   write(99,'(a,i4,a,1p,e15.8,a,e15.8)') 'panel=',npanel_diag,' elo=',elo_save,' ehi=',enext
   close(99)
endif
```

ALWAYS restore: `cd njoy-reference && git checkout -- src/` then full rebuild.

## Files Changed

| File | What changed |
|------|-------------|
| `src/orchestration/modules/gaminr.jl` | Coherent x-space GL, incoherent GL-6 energy quad, reaction rate σ×W interp, wt breakpoint splitting |

## Next Steps for Next Agent

1. **Fix MF23 flux** — implement rndoff/delta nudges matching Fortran gpanel. Patch Fortran to print flux per panel, compare with Julia.
2. **Fix MF26/MT504 MAT=92** — patch Fortran to print running ans totals per panel, find where Julia diverges.
3. **Fix MF26/MT516** — compare pair production threshold handling with Fortran gtff MT=516.
4. **Fix matxsr** — full banded transfer matrices (needed for tape35 but not tape37).
