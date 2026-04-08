# T03 Pipeline — Phase 3 Handoff

## Date: 2026-04-08

## Goal

Get T03 tape37 (viewr PostScript) passing the official execute.py test at rel_tol=1e-5.

T03 chain: `reconr → gaminr → dtfr → matxsr → viewr`. Only **referenceTape37** (9274 lines) is compared.

## What Was Done This Session

### 1. gaminr MT=621 total heating accumulation (CRITICAL FIX)

**Root cause**: Julia computed MT=621 heating only from photoelectric (MT=602). The Fortran accumulates heating from ALL reactions: incoherent (MT=504), pair production (MT=516), AND photoelectric (MT=602).

**How it was found**: Grind method — first diff in plot tape at line 5 showed Y-axis max 0.83E+06 (Fortran) vs 0.79E+01 (Julia). Traced to MT=621 values. Patched Fortran gaminr.f90 to print toth accumulation via file output (unit 99 → `/tmp/gtff_diag.txt`). Compared intermediate values.

**Key Fortran code (gaminr.f90:388-395)**:
```fortran
if (ng2.ne.2.and.mtd.ne.502) then
   toth(l+1)=toth(l+1)+ans(1,1,ng2)*ans(1,1,1)
```
This fires for MT=504/MF26, MT=516/MF26, AND MT=602/MF23 (because MT=602 has ng=2 feed columns from gtff line 1218: `if (mtd.eq.602) ng=2; ff(1,2)=e`).

**Important**: `dspla` (lines 1066-1074) normalizes `ans(1,1,i) /= ans(1,1,1)` IN-PLACE before the heating accumulation. So `ans(1,1,ng2)` is already the heating KERMA (eV·barn), and `ans(1,1,1)` is still the raw flux. The toth accumulation stores `heating_kerma × flux`. At MT=621 output, `dspla` normalizes again: `toth(l+1)/toth(l)` = total heating KERMA.

**Fix**: Added `toth` array to `_write_gaminr_tape`. Each reaction accumulates its heating into toth. MT=621 outputs accumulated toth normalized by flux.

**File**: `src/orchestration/modules/gaminr.jl` — `_write_gaminr_tape`, `_write_gaminr_mt621`

### 2. Incoherent scattering matrix normalization (CRITICAL FIX)

**Root cause**: Julia's `_gaminr_incoherent_matrix!` multiplied `σ_MF23 × KN × S(q)`, giving σ² (double-counting). The Fortran normalizes the angular distribution by siginc (total KN integral) at gtff lines 1456-1464:
```fortran
ebar=ebar/siginc              ! average scattered photon energy
ff(1,ng+2)=(e-ebar)*siginc    ! heating = (E-<E'>) × σ_KN
do il=1,nl; do i=1,ng
   ff(il,i)=ff(il,i)/siginc   ! NORMALIZE angular fractions
enddo; enddo
```
Then gpanel multiplies by σ_MF23: `ans += σ_MF23 × W × ff_normalized × dE`.

**Fix**: Rewrote `_gaminr_incoherent_matrix!` to compute KN×S(q) angular distribution, normalize by siginc, then multiply by σ_MF23. Heating = σ_MF23 × (E - <E'>).

**Also fixed**: S(q) was squared (`ff_val * ff_val`). MF27/MT=504 stores the incoherent scattering function S(q,Z), not the coherent form factor F(q). Must be used linearly. Fortran gtff line 1427: `fact = snow * (...)` — snow = S(q), not S(q)².

**Impact**: MT=621 group 1 went from 7.9 → 4488 eV·barn (matching Fortran). High-energy groups still ~12% off due to crude angular integration (see remaining work).

### 3. GENDF section ordering (interleaved MF23/MF26)

**Root cause**: Julia wrote all MF23 sections first, then all MF26. The Fortran processes reactions in `mtlst/mflst` order (gaminr.f90:111-114):
```
mtlst = [501, 502, 502, 504, 504, 516, 516, 602, 621]
mflst = [ 23,  23,  26,  23,  26,  23,  26,  23,  23]
```
This interleaves MF23 and MF26 for each MT. No FEND records between interleaved sections — only SEND after each section, MEND after material.

**Fix**: Changed `_write_gaminr_tape` to process `reaction_sequence` in Fortran order. Removed per-MF FEND records.

### 4. MF=26 GENDF format fixes

- **Flux column**: MF26 LIST records include flux as position 1 (NL Legendre moments, all equal to flux). Previously missing.
- **Coherent (MT=502)**: ng2=2 (flux + self-scatter, no total/heating). Previously wrote ng2=3.
- **Below-threshold skip**: MF23/MT=516 now skips groups 1-3 (below 1.022 MeV pair production threshold), matching Fortran `igzero` logic.

### 5. `_dpend` thinning state machine (STRUCTURAL FIX)

**Root cause**: Julia only implemented _dpend thinning states 0-2. Fortran (dtfr.f90:1240-1268) has states 0-5 via computed GOTO: `go to (230,240,240,250,270),ns`. States 3-5 sort and keep extremes within a pixel column.

**Fix**: Implemented full state machine with `@goto` labels matching Fortran labels 230, 240, 250, 260, 270, 280.

**Impact**: Plot tape went from 1627 → **1711 lines** (exact structural match with Fortran).

## Current State

```
tape33 (GENDF):  632 lines (target: 604, diff: +28)
tape36 (plot):  1711 lines (target: 1711) ✓ STRUCTURAL MATCH
tape37 (PS):    9110 lines (target: 9274, diff: -164)
```

### MT=621 heating values (MAT=1):
| Group | Energy range | Fortran | Julia | Error |
|-------|-------------|---------|-------|-------|
| 1 | 10-100 keV | 4488 | 4488 | **0.0%** |
| 2 | 100-500 keV | 22806 | 22857 | 0.2% |
| 3 | 0.5-1 MeV | 70470 | 70005 | 0.7% |
| 4 | 1-2 MeV | 123909 | 124450 | 0.4% |
| 5 | 2-3 MeV | 178392 | 182290 | 2.2% |
| 8 | 5-6 MeV | 279822 | 297060 | 6.2% |
| 12 | 9-20 MeV | 412977 | 463343 | 12.2% |

Errors grow at high energies where pair production dominates and the crude angular integration matters most.

### GENDF line count breakdown (+28 lines):
| Section | Ref | Julia | Diff | Root cause |
|---------|-----|-------|------|-----------|
| MF23/MT516 | 19 | 25 | +6 | Below-threshold skip not applied to MF26 path |
| MF26/MT504 | 85 | 71 | -14 | Fewer nonzero downscatter groups (angular precision) |
| MF26/MT516 | 19 | 41 | +22 | Pair production format wrong |
| (×2 materials) | | | ×2 | Same pattern for MAT=1 and MAT=92 |

### Regression check:
- T01 tape25: NOT checked this session (should verify)
- T02 tape28: NOT checked this session

## Remaining Blockers — Priority Order

### 1. GENDF incoherent scattering precision (28 format lines + 12% heating error)

**The fundamental issue**: Julia uses 20-point midpoint rule over cos θ for the incoherent angular integration. Fortran uses Lobatto quadrature (6 or 10 points) over scattered photon momentum p', with panels defined by group energy boundaries AND MF27 form factor breakpoints. The Fortran approach is superior because:
- Integration variable p' naturally maps to sink groups (each p' range = one sink group)
- Form factor breakpoints define natural panel boundaries
- The KN×S(q) product varies smoothly in p' but oscillates in cos θ

**How to fix**: Port the Fortran `gtff` MT=504 integration (gaminr.f90:1341-1464):
1. Convert energy E to dimensionless `enow = c3*E`
2. Loop over panels in p' space defined by group boundaries + form factor breakpoints
3. Use 6-point Lobatto quadrature within each panel
4. At each quadrature point: compute cos θ from p', look up S(q) via terpa, compute KN factor
5. Accumulate into sink group determined by p' range
6. Normalize by siginc, set heating = (E - ebar)*siginc

This is ~120 lines of Fortran to port. The key helper functions needed:
- `c3, c4, c5` constants (already in gaminr.jl)
- `terpa` for form factor interpolation (can use existing `_interp_ff`)
- Group boundary → p' conversion: `p_boundary = c3 * E_boundary`

### 2. MAT=92 MT=621 missing

Julia doesn't write MT=621 for uranium because the form factor lookup fails. Check:
- Does `ff_data[92]` exist? (The gam27 file has both MAT=1 and MAT=92 data)
- The `_read_mf27_form_factors` reader may fail for MAT=92

### 3. Pair production MF26 format (+22 lines per material)

Julia writes MT=516/MF26 with ng2 values larger than Fortran. Check:
- Fortran pair production `gtff` at lines 1487-1508 — what ng2 does it produce?
- The pair production is isotropic (NL=1 not NL=5) — Fortran line 300: `if (mtd.eq.516) nl=1`

### 4. MF23/MT516 below-threshold skip for MF26

Julia applies the igzero skip to MF23 but not MF26. For pair production MF26, groups below 1.022 MeV should also be skipped. Add the same igzero logic to the MF26 writer.

### 5. PENDF extra sections (MF2/MT151, MF3/MT103)

Julia reconr writes MF2/MT151 and MF3/MT103 for photoatomic materials. Fortran doesn't — photoatomic PENDF has only MF1 + MF23. These add 492 lines to tape31. Impact: minor flux/XS differences in gaminr.

## Key Files

| File | What changed | What to change next |
|------|-------------|-------------------|
| `src/orchestration/modules/gaminr.jl` | Heating accumulation, MF26 writer, section ordering | Port gtff incoherent integration, fix pair production |
| `src/orchestration/modules/dtfr.jl` | `_dpend` state machine | Done for now |
| `src/processing/gaminr.jl` | Group structures, weight functions | Unchanged |
| `njoy-reference/src/gaminr.f90` | Reference source | Read lines 1341-1514 for gtff incoherent+pair |

## How to Run

```bash
# Generate Fortran oracle (if not cached)
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

# Compare tape37 (official test target)
python3 -c "
ref = open('/tmp/t03_fortran/tape37').readlines()
trial = open('/tmp/t03_test/tape37').readlines()
print(f'Lines: ref={len(ref)}, julia={len(trial)}')
if len(ref) != len(trial): print('STRUCTURAL FAIL')
"

# Compare GENDF sections
python3 -c "
ref = open('/tmp/t03_fortran/tape33').readlines()
trial = open('/tmp/t03_test/tape33').readlines()
print(f'GENDF: ref={len(ref)}, julia={len(trial)}')
"
```

## How to Use gdb Diagnostics on gaminr

```bash
# Patch Fortran gaminr.f90 with write diagnostics
# Use file output (unit 99) since NJOY captures stdout
vi njoy-reference/src/gaminr.f90
# Add: open(99,file='/tmp/diag.txt',position='append')
#      write(99,*) 'TAG', variable
#      close(99)
cd njoy-reference/build && cmake --build . --target njoy
rm -f /tmp/diag.txt
cd /tmp/t03_fortran && ~/Projects/NJOY.jl/njoy-reference/build/njoy < input > /dev/null 2>&1
cat /tmp/diag.txt
# ALWAYS restore: cd njoy-reference && git checkout -- src/
```

## Traps

**Trap 148 (NEW — FIXED)**: Fortran gtff sets `ff(1,2)=e` for MT=602 MF=23 (gaminr.f90:1226). This creates a heating column (ng=2 feed columns → ans has 3 columns). After dspla normalization, `ans(1,1,3)` = heating KERMA = ∫E×σ×W/∫W. The accumulation `toth += heating_kerma × flux` fires because `ng2=3 ≠ 2`. Julia previously computed heating only from `_gaminr_heating_average` and didn't accumulate it into toth.

**Trap 149 (NEW — FIXED)**: Fortran gtff normalizes the incoherent angular distribution by siginc (lines 1460-1462: `ff(il,i)=ff(il,i)/siginc`). Then gpanel multiplies by σ_MF23. Without normalization, the product is σ_MF23 × σ_KN ≈ σ² — giving 3-4x too-small values for hydrogen and much worse for high-Z. Julia's previous code multiplied `smid * dsig` = σ_MF23 × KN × S(q) without normalizing.

**Trap 150 (NEW — FIXED)**: MF27/MT=504 stores the incoherent scattering function S(q,Z), NOT the coherent form factor F(q). For incoherent scattering, use S(q) linearly in the KN formula (not S²). Fortran gtff line 1427: `fact = snow * (...)`. For coherent (MT=502), the squared form factor |F(q)|² IS correct.

**Trap 151 (NEW — FIXED)**: Fortran gaminr does NOT write FEND between interleaved MF=23/MF=26 sections. Only SEND (asend) after each section, MEND after material. The section ordering follows `mtlst/mflst` arrays.

**Trap 152 (NEW — FIXED)**: `_dpend` uses a 5-state machine (Fortran computed GOTO at line 1242). States 3-5 sort and keep extremes within a pixel column. Julia previously only had states 0-2, producing 84 fewer overlay data points.

**Trap 153 (NEW — FIXED)**: GENDF MF26 LIST records include flux as position 1 with NL Legendre moments (all equal to the group flux value). Julia previously omitted the flux position.

**Trap 154 (NEW — FIXED)**: For coherent scattering (MT=502), Fortran writes ng2=2 (flux + self-scatter). No total or heating columns. Julia was writing ng2=3.

**Trap 155 (NEW — FIXED)**: Fortran skips below-threshold groups in GENDF output via `igzero` flag (dspla sets it when finding nonzero values). MT=516 pair production threshold is 1.022 MeV — groups 1-3 are skipped.

**Trap 156 (NEW)**: Fortran pair production (MT=516 MF=26) uses NL=1 (isotropic, no angular dependence) per gaminr.f90 line 300: `if (mtd.eq.516) nl=1`. Julia may be using NL=5.

**Trap 157 (NEW)**: Fortran `dspla` for MF=23 normalizes ans values in-place (line 1070-1071: `ans(1,1,i)=result(i-1)` where `result=ans/flux`). This happens BEFORE the heating accumulation at line 392. The heating value in `ans(1,1,ng2)` is already flux-normalized when accumulated into toth.

**Trap 158 (NEW)**: For Fortran gdb diagnostics on gaminr, use file output (`open(99,file=...,position='append')`) not `write(*,*)`. NJOY captures stdout for its own output. Also cannot use stderr (unit 0) — NJOY may redirect it.
