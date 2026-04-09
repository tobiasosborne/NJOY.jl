# T03 Pipeline — Phase 4 Handoff

## Date: 2026-04-09

## Goal

Get T03 tape37 (viewr PostScript) passing the official execute.py test at rel_tol=1e-5.

T03 chain: `reconr → gaminr → dtfr → matxsr → viewr`. Only **referenceTape37** (9274 lines) is compared.

## What Was Done This Session

### 5 bugs found and fixed in `src/orchestration/modules/gaminr.jl`:

1. **KN formula 2x factor (CRITICAL)**: Julia's Klein-Nishina integrand had `0.5 * c2 * ratio^2 * (...)` — the `0.5` was spurious. The `c2 = r_e²/2` already includes the 1/2. The mu-space Jacobian `|dp/dmu|` accounts for the rest. Removing the `0.5` doubled siginc from 0.040 → 0.080 (Fortran: 0.083). Found by computing the Jacobian `fact * |dp/dmu|` vs `kn` and getting ratio=0.5 exactly.

2. **MT=516 pair production NL=1 (STRUCTURAL)**: Pair production is isotropic — Fortran forces `nl=1` at gaminr.f90:300 (`if (mtd.eq.516) nl=1`) even when `lord > 0`. Julia was using `nl=lord+1=5`, producing 5 Legendre moments per entry. Fix: pass `nl_mt = (mtd == 516) ? 1 : nl` to both `_gaminr_scatter_matrix` and `_write_gaminr_mf6_section`. Impact: -22 lines/material on MT=516 MF=26.

3. **Gauss-Lobatto p-space integration (CRITICAL)**: Complete rewrite of `_gaminr_incoherent_matrix!`. Old code used 20-point uniform midpoint in mu-space. New code ports Fortran `gtff` MT=504 (gaminr.f90:1338-1464): Gauss-Lobatto 6-point quadrature in scattered-photon p-space with adaptive panels at group boundaries + form factor breakpoints. Key implementation details:
   - `zz = ff_vals[end]` = Z (atomic number) for incoherent — used to decide table vs constant path
   - `xzz = c5*sqrt(enow/500)` compared against `zz`: for hydrogen (Z=1) at E>500 keV, uses constant S=1 (no table)
   - `ff_ip` cursor tracks position in form factor table (matches Fortran's `terpa` state)
   - Must advance `ff_ip` INSIDE the GL quadrature loop (line 1417 equivalent) to prevent stale `xnext` causing infinite panel loops
   - `arg[]` initialized to 0, reused across panels via GL endpoint sharing (iq=1 skips recomputation)
   - Panel count: 23-31 per energy point (confirmed via Fortran gdb diagnostic)
   - `_interp_ff` upgraded to binary search (O(log N)) for performance with 45-point tables

4. **igzero skip in MF=26 writer**: Below-threshold groups with all-zero scatter data were being written. Fortran uses `igzero` flag in `dspla` (line 1072) to skip them. Added `igzero_seen` flag in `_write_gaminr_mf6_section` — only writes when nonzero data found or at last group. Impact: -6 lines/material for MT=516.

5. **Per-energy pair production threshold**: Used `elo < 2*epair` (per-group skip) instead of per-energy-point threshold. Group 4 (1-2 MeV) contains the 1.022 MeV pair production threshold — the above-threshold portion has nonzero XS. Fix: skip only groups where `ehi <= 2*epair`, and within the group, clip panel start to threshold. Also strip "total" column for MT=516 (same as yield → redundant). Impact: +2 lines/material restored.

### New constants added:
- `_GAMINR_C4 = 0.0485262` (form factor x→q conversion)
- `_GAMINR_C5 = 20.60744` (dimensionless→x conversion)
- `_GAMINR_RNDOFF = 1.0000001` (breakpoint nudge)
- `_GAMINR_CLOSE = 0.99999` (backward scatter threshold)
- `_GL6_PTS/_GL6_WTS` (6-point Gauss-Lobatto)
- `_GL10_PTS/_GL10_WTS` (10-point Gauss-Lobatto)

### Fortran oracle generated:
Oracle at `/tmp/t03_fortran/` with tape31(2180), tape33(604), tape34(340), tape35(211), tape36(1711), tape37(9274). Generated from `njoy-reference/tests/03/input` + `gam23`(tape30) + `gam27`(tape32).

## Current State

```
tape31 (PENDF):  2672/2180  (+492)    0.0% match  — extra MF2/MF3 from reconr (known)
tape33 (GENDF):   604/604   MATCH    42.1% match  — value diffs only
tape34 (DTF):     340/340   MATCH    34.7% match  — tracks GENDF improvement
tape35 (MATXS):   131/211   (-80)     6.9% match  — diagonal-only transfer matrices
tape36 (plot):   1711/1711  MATCH    79.8% match  — tracks GENDF values
tape37 (PS):     9110/9274  (-164)   19.4% match  — official target
```

### MT=621 heating values (MAT=1, hydrogen) after fixes:
Values NOT re-measured after GL fix. The KN 2x fix was measured before GL and showed same heating errors as Phase 3. Need to re-measure — the GL integration should significantly improve high-energy heating accuracy.

### Regression checks:
- T01 tape25: 32962/32962 lines, 212 diffs at 1e-9 — UNCHANGED ✓
- T02: not checked (gaminr changes don't affect T02 pipeline)

## Remaining Blockers — Priority Order

### 1. GENDF value diffs (tape33: 58% lines differ)

**Root cause**: Julia's gpanel-equivalent integrates σ×W over PENDF panels using simple trapezoidal rule with midpoint evaluation. The Fortran `gpanel` (gaminr.f90:874-1010) uses Gauss-Lobatto quadrature in ENERGY space (separate from the angular GL quadrature in gtff). gpanel also breaks at weight function TAB1 breakpoints (from `gtflx`/`terpa`), creating finer sub-panels than Julia's PENDF-panel-based integration.

**How to diagnose**: Patch Fortran gpanel to print the number of sub-panels per group for MT=501 MF=23 MAT=1. Compare the energy grid resolution.

**How to fix**: Either (a) implement gpanel's energy quadrature matching the Fortran (Lobatto with breakpoints from gtsig+gtflx+gtff), or (b) add sub-panel splitting at weight function breakpoints. Approach (b) is simpler — the weight function for iwt=3 has only 4 breakpoints [1e3, 1e5, 1e7, 3e7], so adding panel breaks there would capture most of the difference.

**Expected impact**: Would improve the flux-weighted energies (currently ~0.5 eV off) and the group-averaged XS (currently ~0.1-1% off). This cascades to MT=621 heating and all downstream tapes.

### 2. matxsr transfer matrices (tape35: -80 lines)

**Root cause**: `src/orchestration/modules/matxsr.jl` writes diagonal-only P0 transfer matrices (jband=1 for all groups). The Fortran writes full banded matrices using the actual Compton downscatter bandwidth. For hydrogen incoherent scattering, a 5 MeV photon can scatter into groups 2-7, so the matrix has bandwidth 6, not 1.

**How to fix**: In `_write_matxsr_section`, compute `jband[g]` = number of nonzero sink groups for source group g, and `ijj[g]` = lowest sink group index. Write all Legendre moments for all nonzero entries in the band. This data is already computed in the GENDF (MF=26 sections) — matxsr just needs to read it properly.

**Expected impact**: tape35 from 131→~211 lines. Won't directly affect tape37 since dtfr reads GENDF, not MATXS.

### 3. tape31 PENDF extra lines (+492)

**Root cause**: Julia reconr writes MF2/MT151 and MF3/MT103 for photoatomic materials. Fortran doesn't — photoatomic PENDF has only MF1 + MF23. These extra sections don't affect gaminr (which only reads MF23), but they make the PENDF different.

**How to fix**: In the photoatomic reconr path, skip MF2/MF3 output. Only write MF1/MT451 + MF23 sections.

**Low priority** — doesn't affect tape37.

### 4. tape37 line count diff (-164)

**Root cause**: Value differences in tape33 → different heating values in tape34 → different Y-axis ranges in tape36 plot commands → viewr renders different number of data points. When tape33 values match, tape37 will match (viewr is 100% BIT-IDENTICAL on correct input).

**Not a separate bug** — will resolve as tape33 values improve.

## Key Traps (NEW this session)

**Trap 159 (FIXED)**: The KN formula `dσ/dΩ = (r_e²/2)(E'/E)²[(E'/E)+(E/E')-sin²θ]` when integrated over dμ in mu-space needs the Jacobian |dp/dmu|. The factor `c2 = r_e²/2 = 0.249467` already includes the 1/2 from the Thomson cross section. Julia's extra `0.5 *` factor produced a 2x systematic error in siginc.

**Trap 160 (FIXED)**: `zz` in Fortran gtff is the MF27 TAB1 C2 field = Z (atomic number). For incoherent scattering, `S(x→∞) = Z`. The `ff_vals[end]` of the table equals Z. The comparison `xzz <= zz` determines whether to use the form factor table or the constant S=Z approximation. For hydrogen (Z=1) at E>500 keV: constant. For uranium (Z=92) at all photon energies: table (since xzz < 5.6 < 92).

**Trap 161 (FIXED)**: Fortran `terpa` is STATEFUL — `ip,ir` cursors advance monotonically through the form factor table. Julia's `_interp_ff` is stateless. Must track `ff_ip` explicitly and advance it both in the inner breakpoint loop AND inside the GL quadrature loop (where `terpa` is called at line 1417). Without advancing `ff_ip` inside the GL loop, the same `xnext` breakpoint persists across outer panel iterations, creating infinite tiny panels.

**Trap 162 (FIXED)**: MT=516 pair production threshold check must be per-energy, not per-group. Group 4 (1-2 MeV) spans the 1.022 MeV threshold — the above-threshold portion has nonzero XS. Using `elo < 2*epair` skips the entire group. Fix: `ehi <= 2*epair` for group skip, then clip panel start to threshold within the group.

**Trap 163 (FIXED)**: MT=516 pair production in GENDF has ng2=2 (flux + yield), NOT ng2=3. The "total" column is redundant with the yield. Must pass `strip_total=true` for MT=516 as well as MT=504.

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

# Compare tapes
for u in 33 34 36 37; do
  echo "tape$u:"
  wc -l /tmp/t03_fortran/tape$u /tmp/t03_test/tape$u
done
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
# Patch gaminr.f90 with diagnostics
# Use file output (unit 99) since NJOY captures stdout
vi njoy-reference/src/gaminr.f90
# Add: open(99,file='/tmp/diag.txt',position='append')
#      write(99,*) 'TAG', variable
#      close(99)
# DECLARE all new variables (e.g., npanel) in the "! internals" section
cd njoy-reference/build && rm -rf CMakeFiles libnjoy.so && cmake .. && cmake --build . --target njoy
rm -f /tmp/diag.txt
cd /tmp/t03_fortran && ~/Projects/NJOY.jl/njoy-reference/build/njoy < input > /dev/null 2>&1
cat /tmp/diag.txt
# ALWAYS restore: cd njoy-reference && git checkout -- src/
# ALWAYS full rebuild after restore
```

## Files Changed

| File | What changed |
|------|-------------|
| `src/orchestration/modules/gaminr.jl` | KN fix, NL=1, GL integration, igzero, pair prod threshold, constants, binary search interp |

## Next Steps for Next Agent

1. **Improve GENDF energy values** — implement gpanel energy quadrature or add weight function breakpoints to the integration panels. This is the highest-impact fix for tape37.
2. **Fix matxsr transfer matrices** — full banded output from GENDF scatter data. Needed for tape35 but not tape37.
3. **Grind tape33 value diffs** — patch Fortran gpanel to print flux/XS at each sub-panel, compare with Julia's trapezoidal. Find the dominant error source.
