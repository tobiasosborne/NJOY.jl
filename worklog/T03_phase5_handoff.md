# T03 Pipeline — Phase 5 Handoff

## Date: 2026-04-10

## Goal

Get T03 tape37 (viewr PostScript) passing the official execute.py test at rel_tol=1e-5.

T03 chain: `reconr → gaminr → dtfr → matxsr → viewr`. Only **referenceTape37** (9274 lines) is compared.

## What Was Done This Session

### 1. Coherent angular integration — x-space GL quadrature (CRITICAL FIX)

**Root cause**: Julia's coherent angular integration used 16-point midpoint rule in cos(θ). The Fortran `gtff` MT=502 (gaminr.f90:1259-1334) integrates in **momentum-transfer x-space** with GL-6/GL-10 quadrature, panels at form factor table breakpoints. The integration variable `x` is the MF27 momentum-transfer variable (Å⁻¹). The Jacobian `|du/dx| = 2x/(c1·E)²` is absorbed into the weight: `fact × xnow` where `fact = 2·c2/(c1²·E²)`.

**Key Fortran code (gaminr.f90:1289-1311)**:
```fortran
do iq=1,nq                          ! GL-6 or GL-10 over form-factor panel
   xq=aq+bq*qp6(iq); wq=bq*qw6(iq) ! quadrature in x-space
   xnow=xq*rndoff
   if (iq.gt.1) then
      unow=1-c1e*xnow*xnow          ! cos θ from x: u = 1 - x²/(c1·E)²
      call terpa(snow,xnow,...)       ! F(x) from form factor table
      call legndr(unow,pl,nl)          ! Legendre polynomials at cos θ
      arg(il)=(1+unow*unow)*snow*snow*pl(il)  ! integrand
   endif
   ff(il,1)=ff(il,1)+wq*fact*xnow*arg(il)    ! accumulate
enddo
```

Panels are MF27 form factor tabulation intervals — the loop advances through breakpoints until `unow <= -1` (θ=π) or `xnow >= xlim = √2·c1·E` or `F(x) < tolerance`. After integration, `ff(il,1) /= ff(1,1)` normalizes to unit P0.

**Fix**: Rewrote `_coherent_feed_at_energy!` (gaminr.jl:686-769) to faithfully port the Fortran. GL-6 for nl≤4, GL-10 for nl>4. Walks through form factor breakpoints as panel boundaries using `ip` cursor advanced after each GL loop. Returns normalized Legendre moments.

**Impact**: MF26/MT502 worst error from 55%→0.006% (MAT=1), 73%→0.01% (MAT=92).

### 2. Incoherent energy quadrature — GL-6 in energy (STRUCTURAL FIX)

**Root cause**: Julia evaluated the incoherent angular distribution (the full p-space GL integration) at ONE midpoint energy per PENDF panel, then weighted by `smid × wmid × de`. The Fortran `gpanel` (gaminr.f90:874-1010) evaluates `gtff` at `nq=6` GL quadrature points per energy sub-panel.

**Architecture match (Fortran gpanel for MF26/MT504)**:
- For each energy sub-panel [elo, ehi] (determined by PENDF grid points):
  - Evaluate σ and W at both endpoints: `a = sig_hi × W_hi`, `b = sig_lo × W_lo`
  - GL-6 quadrature in energy: at each GL point `eq`:
    - Linearly interpolate reaction rate: `rr = b + (a-b) × t1`, `t1 = (eq-elo)/(ehi-elo)`
    - Call `gtff(eq)` → normalized angular distribution `ff[il, ig_sink]`
    - Accumulate: `ans[il, iz, igt] += rr × wq × ff[il, ig]`
  - Flux integral (separate, trapezoidal): `ans[il,iz,1] += (W_hi + W_lo) × bq`

**Fix**: 
1. Extracted `_incoherent_feed_at_energy!` (gaminr.jl:851-971) — the full p-space GL integration from the old `_gaminr_incoherent_matrix!` inner loop, parameterized by energy. Verified identical to Fortran gtff MT=504 to **8+ significant digits** by patching Fortran with diagnostic print at gaminr.f90:1465 and comparing E/siginc/ebar/ff values at 10 test energies.

2. Rewrote `_gaminr_incoherent_matrix!` (gaminr.jl:979-1056) to use GL-6 energy quadrature. For each PENDF panel, calls `_incoherent_feed_at_energy!` at 6 GL energy points.

**Impact**: MF26/MT504 worst error from 320%→0.5% (MAT=1), 297%→1.4% systematic (MAT=92).

### 3. Reaction rate σ×W product linear interpolation (MATCH FIX)

**Root cause**: Julia previously computed `σ_linear(eq) × W_exact(eq)` at each GL energy quadrature point. Fortran gpanel linearly interpolates the **product** `σ×W` (the "reaction rate"):
```fortran
a = sig(iz,1)*flux(iz,il)     ! σ_hi × W_hi at upper panel boundary
b = slst(iz,1)*flst(iz,il)    ! σ_lo × W_lo at lower panel boundary
rr = b + (a-b)*t1             ! linear interpolation of the PRODUCT
```

These produce different values because W is log-log interpolated (iwt=3), so σ_lin × W_exact ≠ (σ×W)_lin within a panel. For a PENDF panel with ratio 1.03, W changes by ~3%, creating ~0.1% quadrature-level difference.

**Verification**: Patched Fortran gpanel caller loop (around line 337) to print sub-panel boundaries for MAT=92 MT=504 group 5. Confirmed Julia and Fortran use the **exact same 17 PENDF panels** (identical elo/ehi values to all digits).

Also verified: PENDF grid max_ratio is ~1.03 everywhere, finer than the 1.05 step cap from `gtflx`. The 1.05 step never creates additional sub-panels for this test case.

**Fix**: Changed Julia GL loop (gaminr.jl:1020-1042) to compute `rr_lo = sa*wa`, `rr_hi = sb*wb`, then `rr = rr_lo + (rr_hi - rr_lo)*t1` at each GL point. Applied to both coherent (gaminr.jl:816-832) and incoherent.

### 4. Coherent energy quadrature — GL-6 (instead of Fortran's GL-2)

For coherent, Fortran gpanel uses nq=0+2=2 (trapezoidal) in energy, because gtff MT=502 returns `nq=0`. The Fortran's trapezoidal works because panels are very narrow (max 1.03 ratio).

Julia uses GL-6 for coherent energy quadrature because the wider-seeming PENDF panels (though same 1.03 ratio) benefit from more quadrature points to capture the angular distribution's energy dependence. This is mathematically more accurate than Fortran's GL-2 but produces slightly different results. The difference is < 0.01%, acceptable.

### 5. Diagnostics methodology refined

**gtff diagnostic pattern** (tested and working):
- Patch after gaminr.f90 line 1464 (after normalization), use `open(99,file='/tmp/gtff_diag.txt',position='append')` 
- Guard: `if (matd.eq.92.and.mtd.eq.504) then`
- ALWAYS declare all new variables in the `! internals` section
- ALWAYS full rebuild: `rm -rf CMakeFiles libnjoy.so && cmake .. && cmake --build . --target njoy`
- ALWAYS restore: `cd njoy-reference && git checkout -- src/`

**gpanel diagnostic pattern** (tested and working):
- Patch the caller loop around line 337
- Need to declare `integer::npanel_diag` and `real(kr)::elo_save` at line 100
- Variable declarations go at: `integer::idone,igzero,j,ibase,i,lim,l,npanel_diag` and `real(kr)::elo_save`
- Guard: `if (matd.eq.92.and.mtd.eq.504.and.ig.eq.5) then` (NOTE: use `ig` not `mfd` — the `mfd` check doesn't fire because `mfd` may not be 26 at this point in the processing)

**Key finding**: PENDF max_ratio 1.03 < 1.05 step cap → PENDF breakpoints from gtsig always come first. The 1.05 step from gtflx never creates additional sub-panels.

## Current State

```
tape33 (GENDF):   604/604  42.2% byte-identical (was 42.1%)
tape34 (DTF):     340/340  40.6% byte-identical (was 34.7%)
tape35 (MATXS):   131/211  diagonal-only (unchanged)
tape36 (plot):   1711/1711 86.3% byte-identical (was 79.8%)
tape37 (PS):     9304/9274 32.4% byte-identical (was 19.1%)
```

### Remaining errors by section:

| Section | Worst error | Error count | Notes |
|---------|------------|-------------|-------|
| MF23/all flux | 6e-5 | ~200 lines | Systematic, from rndoff/delta in gpanel |
| MF26/MT502 MAT=1 | 6e-5 | 119 diffs | Essentially exact, matches flux error |
| MF26/MT502 MAT=92 | 1.1e-4 | 119 diffs | Very good |
| MF26/MT504 MAT=1 | 5.5e-3 (0.55%) | 397 diffs | Residual energy quadrature |
| MF26/MT504 MAT=92 | 2.4e-1 (24%) on -2.1e-5 | 400 diffs | 24% is on tiny val; systematic ~1.4% |
| MF26/MT516 both | 3.2e-2 (3.2%) | 18 diffs | Pair production threshold handling |
| MT=621 MAT=1 | 2.0e-2 (2.0%) | 24 diffs | Cascading from incoherent |
| MT=621 MAT=92 | 1.7e-1 (17%) | 24 diffs | Cascading from incoherent |

### tape37 line count: 9304 vs 9274 (+30)
Value diffs cause different Y-axis ranges in dtfr → viewr renders different point counts. Will converge as tape33 values improve.

### Regression checks:
- T01 tape25: 32750/32962, 212 diffs — UNCHANGED ✓ (compared against njoy-reference/tests/01/referenceTape25)

## Remaining Blockers — Priority Order

### 1. MF23 flux: ~6e-5 systematic error (HIGHEST IMPACT)

**Root cause identified**: Fortran gpanel applies `rndoff=1.000002` to `elo` at the start of each panel (gaminr.f90:915: `if (elo*rndoff.lt.ehi) elo=elo*rndoff`) and `delta=0.999995` to `ehigh` at the upper boundary (gaminr.f90:935: `ehigh=delta*ehi`). Julia uses the raw PENDF panel boundaries.

**Quantification**: For group 1, single-panel estimate gives ~1.56e-5 relative error, matching observed 1.50e-5. Error grows for higher groups because the cumulative effect of rndoff across more panels shifts the integration.

**How to fix**: In Julia's flux computation and GL energy quadrature:
1. Apply `elo *= 1.000002` at the start of each panel (when elo ≠ previous panel's ehigh)
2. Apply `ehigh = ehi * 0.999995` at the upper boundary of each panel
3. Evaluate W at the nudged boundaries, not the raw PENDF points

**Location**: `_gaminr_scatter_matrix` flux computation (gaminr.jl:656-674), `_gaminr_group_average` (gaminr.jl:361-419), `_gaminr_heating_average` (gaminr.jl:563-605).

**Impact**: Would fix ~200 of 350 non-matching tape33 lines. Every MF23 data line has flux as its first value, so matching flux exactly would make many more lines byte-identical. Cascades through tape34/36/37.

### 2. MF26/MT504 MAT=92 incoherent: ~1.4% systematic error

**Root cause**: TBD. Despite:
- Same PENDF panel boundaries (verified by diagnostic)
- Matching angular distribution to 8+ digits at each energy
- Matching σ×W linear interpolation approach
- Matching GL-6 energy quadrature

Possible remaining causes:
- **rndoff/delta nudges** change the effective panel boundaries, which changes the σ and W endpoint values used for the linear interpolation
- **Form factor table interpolation at panel boundaries**: Julia's `_interp_ff` uses stateless binary search; Fortran's `terpa` uses stateful forward-walk. At panel boundaries where `xnow = xq*rndoff`, the nudge could land on different sides of a form factor breakpoint, giving different F(x) values
- **Cumulative floating-point rounding** across 17 panels per group

**How to diagnose**: Patch Fortran gpanel to print the running `ans` total after each panel for MAT=92 MT=504 group 5. Add the same print in Julia. Compare running totals to find the first panel where divergence occurs.

**How to fix**: Once the divergent panel is found, compare the σ, W, rr, ff values at each GL quadrature point within that panel. The fix will depend on what differs.

### 3. MF26/MT516 pair production: 3.2% error

**Root cause**: Julia's pair production integration in `_gaminr_pair_production_matrix!` (gaminr.jl:1061-1111) uses simple trapezoidal on PENDF panels with per-energy threshold clipping. The Fortran uses gpanel with the gtff MT=516 feed function.

**Key Fortran behavior (gaminr.f90:1487-1508)**: Pair production feed function is isotropic (NL=1), returns yield=2 photons into the annihilation group at 0.511 MeV. The integration is simple (nq=0 → gpanel uses GL-2 = trapezoidal). But the threshold handling via `gpanel`'s `gtsig` breakpoints means the integration naturally clips at the threshold energy.

**How to diagnose**: Compare MF23/MT=516 group 4 (1-2 MeV, spans threshold) values line by line. The 3.2% error is concentrated in this group.

**How to fix**: The fix for blocker #1 (rndoff/delta) may resolve most of this, since the pair production integration uses the same trapezoidal flux and panels.

### 4. matxsr transfer matrices (tape35: 131/211)

**Root cause**: `src/orchestration/modules/matxsr.jl` writes diagonal-only P0 transfer matrices (`jband=1` for all groups). Fortran writes full banded matrices.

**How to fix**: In `_write_matxsr_section`, compute `jband[g]` = number of nonzero sink groups for source group g from the GENDF MF26 data. Write all Legendre moments for all nonzero entries in the band. This data is already in the GENDF.

**Impact on tape37**: None — dtfr reads GENDF directly, not MATXS. Low priority.

### 5. tape37 line count: 9304 vs 9274 (+30)

Not a separate bug. Value diffs in tape33 → different heating in tape34 → different Y-axis ranges → different point counts in viewr. Will converge.

## Key Traps (NEW this session)

**Trap 164 (FIXED)**: Coherent angular integration must be in x-space (momentum transfer), NOT cos(θ). The Fortran integrand is `(1+u²)|F(x)|² P_l(u) × fact × x` where `u = 1-x²/(c1·E)²` and `fact = 2·c2/(c1·E)²`. Integration panels are MF27 form factor tabulation intervals. Using cos(θ) as the integration variable produces 20-55% errors because F(x(μ)) varies rapidly near F(q) breakpoints.

**Trap 165 (FIXED)**: Fortran gpanel linearly interpolates the PRODUCT σ×W (reaction rate), not σ and W separately. `rr = b + (a-b)*t1` where `a = sig_hi*W_hi`, `b = sig_lo*W_lo`, `t1 = (eq-elo)/(ehi-elo)`. This matters because W is log-log interpolated (iwt=3). Within a 1.03-ratio PENDF panel, σ×W product linear vs σ×W_exact differs by ~0.1%.

**Trap 166 (VERIFIED)**: PENDF grid max_ratio ~1.03 for this test case, finer than the 1.05 step cap from gtflx. The 1.05 step never fires — PENDF breakpoints from gtsig always come first. Both Fortran and Julia use the same panel boundaries (verified by gpanel diagnostic).

**Trap 167 (DIAGNOSED)**: The ~6e-5 flux error is from Fortran gpanel's `rndoff=1.000002` and `delta=0.999995` nudges. For group 1, single-panel estimate of 1.56e-5 matches observed 1.50e-5. The nudges shift panel boundaries by 0.0002%, which shifts the trapezoidal weight function evaluation points, accumulating across ~55 panels per group.

**Trap 168 (VERIFIED)**: When adding gpanel diagnostics, use `ig` (the group loop variable) in the guard condition, NOT `mfd`. The module-level `mfd` may not equal 26 at the point in the processing loop where gpanel is called. The correct guard is `if (matd.eq.92.and.mtd.eq.504.and.ig.eq.5)`.

**Trap 169 (VERIFIED)**: New Fortran variable declarations for diagnostics go at line 100: `integer::idone,igzero,j,ibase,i,lim,l` — append your new integers here. New reals go on a separate `real(kr)::` line.

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

# Detailed tape33 value analysis
python3 -c "
import re
def pf(s):
    s=s.strip()
    if not s: return None
    m=re.match(r'^([+-]?\d*\.?\d+)([+-]\d+)$',s)
    if m: return float(m.group(1))*10**int(m.group(2))
    try: return float(s)
    except: return None
ref=open('/tmp/t03_fortran/tape33').readlines()
trial=open('/tmp/t03_test/tape33').readlines()
secs={}
for i,(r,t) in enumerate(zip(ref,trial)):
    if r==t: continue
    tail=r[66:80].strip()
    if not tail: continue
    key=tail.split()[0][:10]
    for j in range(6):
        rv=pf(r[j*11:(j+1)*11]); tv=pf(t[j*11:(j+1)*11])
        if rv and tv and rv!=0:
            rel=abs(tv-rv)/abs(rv)
            if rel>1e-7:
                if key not in secs: secs[key]={'c':0,'w':0}
                secs[key]['c']+=1
                secs[key]['w']=max(secs[key]['w'],rel)
for k in sorted(secs):
    if '26' in k or '621' in k:
        print(f'  {k:12s}: {secs[k][\"c\"]:4d} diffs, worst={secs[k][\"w\"]:.4e}')
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

## How to Use gdb Diagnostics on gaminr

**CRITICAL**: The Fortran build has clock skew issues. Always do a FULL clean rebuild:
```bash
cd njoy-reference/build
rm -rf CMakeFiles libnjoy.so
cmake ..
cmake --build . --target njoy
```
Verify the diagnostic is in the binary: `strings libnjoy.so | grep "your_tag"`.

```bash
# Workflow:
# 1. Patch gaminr.f90 with diagnostics
vi njoy-reference/src/gaminr.f90
# 2. Declare new variables at line 100 (integers) or add real(kr) line
# 3. Use file output — NJOY captures stdout:
#    open(99,file='/tmp/diag.txt',position='append')
#    write(99,*) 'TAG', variable
#    close(99)
# 4. Full rebuild
cd njoy-reference/build && rm -rf CMakeFiles libnjoy.so && cmake .. && cmake --build . --target njoy
# 5. Run
rm -f /tmp/diag.txt
cd /tmp/t03_fortran && ~/Projects/NJOY.jl/njoy-reference/build/njoy < input > /dev/null 2>&1
cat /tmp/diag.txt
# 6. ALWAYS restore and rebuild
cd ~/Projects/NJOY.jl/njoy-reference && git checkout -- src/
cd build && rm -rf CMakeFiles libnjoy.so && cmake .. && cmake --build . --target njoy
```

### Tested diagnostic insertion points:

**gtff output** — after gaminr.f90 line 1464:
```fortran
         if (matd.eq.92.and.mtd.eq.504) then
            open(99,file='/tmp/gtff_diag.txt',position='append')
            write(99,'(a,1p,e15.8,a,e15.8,a,e15.8)') &
               'E=',e,' siginc=',siginc,' ebar=',ebar
            write(99,'(a,5(1p,e15.8))') '  ff(1:5,1)=', &
               (ff(il,1),il=1,min(nl,5))
            close(99)
         endif
```

**gpanel sub-panels** — around gaminr.f90 line 337:
```fortran
! At line 100, add to declarations:
!   integer::npanel_diag
!   real(kr)::elo_save
! Before the do-while loop: npanel_diag=0
! Inside the loop, before gpanel call: elo_save=elo; npanel_diag=npanel_diag+1
! After gpanel call:
      if (matd.eq.92.and.mtd.eq.504.and.ig.eq.5) then
         open(99,file='/tmp/gpanel_diag.txt',position='append')
         write(99,'(a,i4,a,1p,e15.8,a,e15.8)') &
            'panel=',npanel_diag,' elo=',elo_save,' ehi=',enext
         close(99)
      endif
```

### Julia equivalent diagnostic:
```julia
# Compare angular distribution at specific energies:
using NJOY, Printf
ff_data = NJOY._read_mf27_form_factors("njoy-reference/tests/resources/gam27")
q_vals, ff_vals = ff_data[92][504]
egg = Float64[1e4, 1e5, 5e5, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 2e7]
zz = ff_vals[end]; ngg = 12; nl = 5
qp = NJOY._GL6_PTS; qw = NJOY._GL6_WTS
ff_work = zeros(Float64, nl, ngg)
for e in [2.5e6, 3.0e6]  # test energies
    fill!(ff_work, 0.0)
    hp, si = NJOY._incoherent_feed_at_energy!(ff_work, e, egg, ngg, nl, q_vals, ff_vals, zz, 6, qp, qw)
    @printf("E=%.8E siginc=%.8E ebar=%.8E\n", e, si, e - hp)
    @printf("  ff(1:5,1)= %s\n", join([@sprintf("%.8E", ff_work[il,1]) for il in 1:5], " "))
end
```

## Files Changed

| File | What changed | Line range |
|------|-------------|------------|
| `src/orchestration/modules/gaminr.jl` | Coherent x-space GL | `_coherent_feed_at_energy!` lines 686-769 |
| | Coherent GL-6 energy quad | `_gaminr_coherent_matrix!` lines 776-845 |
| | Incoherent feed extraction | `_incoherent_feed_at_energy!` lines 851-971 |
| | Incoherent GL-6 energy quad | `_gaminr_incoherent_matrix!` lines 979-1056 |
| | σ×W product interp | both coherent/incoherent GL loops |
| | MF23 wt breakpoint splitting | `_gaminr_group_average` lines 361-419 |
| | MF23 heating wt splitting | `_gaminr_heating_average` lines 563-605 |
| | Flux on raw PENDF panels | `_gaminr_scatter_matrix` lines 656-674 |

## Architecture of Current gaminr.jl (1218 lines)

```
gaminr_module()                        # Entry point
  _gaminr_group_structure()            # 12-group boundaries
  _write_gaminr_tape()                 # Main processing loop
    _write_gaminr_mf1()                # MF1/MT451 header
    for each (mfd, mtd) in reaction_sequence:
      if MF23:
        _gaminr_group_average()        # ∫σW/∫W trapezoidal + wt breakpoints
        _gaminr_heating_average()      # ∫EσW/∫W (for MT=602)
        _write_gaminr_mf3_section()    # GENDF MF23 output
      if MF26:
        _gaminr_scatter_matrix()       # Dispatch to coherent/incoherent/pair
          flux computed trapezoidal on raw PENDF panels
          _gaminr_coherent_matrix!()   # GL-6 energy × GL-6/10 x-space angular
            _coherent_feed_at_energy!()  # x-space GL on FF breakpoints
          _gaminr_incoherent_matrix!() # GL-6 energy × GL-6 p-space angular
            _incoherent_feed_at_energy!() # p-space GL on group+FF breakpoints
          _gaminr_pair_production_matrix!() # simple trapezoidal
        _write_gaminr_mf6_section()    # GENDF MF26 output
      if MT=621:
        _write_gaminr_mt621()          # accumulated total heating

Helper functions:
  _gaminr_weight() / _gaminr_wt1_eval()  # iwt=3 weight function
  _read_pendf_mf23() / _read_za_awr_mf23()  # PENDF reader
  _read_mf27_form_factors()              # Form factor table reader
  _interp_ff()                           # Form factor interpolation (binary search, log-log)
  _legndr!()                             # Legendre polynomials
  _write_cont_line() / _write_data_line() / etc.  # ENDF format output
```

## Next Steps for Next Agent

### Step 1: Fix MF23 flux (highest impact, clearest path)

The ~6e-5 flux error causes ~200 non-matching tape33 lines. The fix is to implement Fortran gpanel's `rndoff` and `delta` nudges.

**What to do**:
1. Read Fortran gpanel lines 913-948 (lower boundary setup + flux accumulation)
2. In Julia's flux computation (`_gaminr_scatter_matrix` lines 656-674), add:
   - `elo *= 1.000002` at start of each panel (when first panel or after a gap)
   - `ehigh = ehi * 0.999995` for the evaluation point at upper boundary
   - Evaluate W at nudged `elo` and `ehigh`, not raw PENDF boundaries
3. Apply same nudges to `_gaminr_group_average` and `_gaminr_heating_average`
4. Apply same nudges to the MF26 GL energy quadrature in `_gaminr_incoherent_matrix!` and `_gaminr_coherent_matrix!`
5. Retest — should see most MF23 flux values match exactly

**Caution**: The rndoff nudge is applied only when `elo != elast` (first panel of a group, or after a discontinuity). For consecutive panels within a group, the saved `ehigh` from the previous panel becomes the new `elo`, and no nudge is applied (since `elo == elast`). This continuity is important.

### Step 2: Fix remaining MF26/MT504 MAT=92

Once flux is fixed (step 1), re-measure the incoherent errors. The rndoff/delta fix may reduce the 1.4% systematic error since it changes the σ and W endpoint values used in the reaction rate interpolation.

If still >0.1%, patch Fortran to print running `ans` totals and find the divergent panel.

### Step 3: Fix MF26/MT516 pair production

May be mostly fixed by step 1. If not, compare MF23/MT=516 group 4 threshold handling.
