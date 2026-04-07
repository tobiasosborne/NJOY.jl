# T03 Pipeline — Phase 2 Handoff

## Date: 2026-04-07

## What Was Done This Session

### 1. Recovered from crashed session

Previous session crashed mid-Opus-subagent while rewriting `dtfr.jl`. The subagent had written ~950 lines of partial code (new GendfMF23Section/GendfMF26Section structs, rewritten GENDF reader with MF23+MF26 parsing, DTF writer with group reversal, plot tape writer with histograms + PENDF overlay + 3D plotnn). Regenerated Fortran oracle at `/tmp/t03_fortran/`.

### 2. Fixed `_dpend` infinite loop (CRITICAL BUG)

**Root cause**: `gety1_interp` closure in `_dpend` returned `enext == e` when `e` equaled a tabulated energy point, causing the plot tape writer to hang forever. Two fixes:
1. `energies_raw[idx + 1] < e` → `<= e` in the while loop (advance past equal points)
2. Early return `if e <= energies_raw[1]` now returns `energies_raw[2]` as `enext` (not `energies_raw[1]`)

**File**: `src/orchestration/modules/dtfr.jl` lines 752-774

### 3. Fixed plot tape formatting (3 issues)

1. **Title trailing spaces**: `"%-6s    %-6s    "` → `"%-6s    %-6s"` (removed 4 trailing spaces)
2. **Histogram data spacing**: `"  %13.4E%13.4E/\n"` → `"%13.4E%13.4E/\n"` (removed 2 leading spaces)
3. **Axis tag E format**: Added `_fmt_e10_2()` helper for Fortran 0p E format (mantissa in [0.1,1) like `0.12E+05` vs Julia's standard `1.20E+04`)

### 4. Added flux storage to GendfMaterial

Added `flux::Vector{Float64}` field to `GendfMaterial` struct and populated it from `dvals[1]` during MF23 parsing. Needed by matxsr for flux weight output.

### 5. Fixed matxsr API migration

Updated `src/orchestration/modules/matxsr.jl` to use new GENDF struct API:
- `gmat.sections` → `gmat.mf23` / `gmat.mf26`
- `sec.data[g]` (Tuple) → `sec.sigma[g]` (Float64) + `gmat.flux[g]`
- Scatter diagonal: `sec.data[g]` → `sec.transfer[1, g, g]` (P0 diagonal from MF26)

### 6. Verified dtfr BIT-IDENTICAL output

**tape34 (DTF): 340/340 lines BIT-IDENTICAL** when reading Fortran GENDF (tape33).

**tape36 (plot tape): 1187/1628 lines match (72.9%)** when reading Fortran GENDF + PENDF.
- First 1182 lines are BIT-IDENTICAL
- Remaining diffs are in uranium PENDF overlay (451 vs 535 points — `_dpend` thinning missing Fortran ns states 3-5)

## Current State

### Pipeline runs end-to-end

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/03/input"; work_dir="/tmp/t03_julia")'
```

Output: tape31 (2672 lines), tape33 (688), tape34 (340), tape35 (131), tape36 (1627), tape37 (9038)

### Fortran oracle

```bash
mkdir -p /tmp/t03_fortran && cd /tmp/t03_fortran
cp ~/Projects/NJOY.jl/njoy-reference/tests/03/input .
cp ~/Projects/NJOY.jl/njoy-reference/tests/resources/gam23 tape30
cp ~/Projects/NJOY.jl/njoy-reference/tests/resources/gam27 tape32
~/Projects/NJOY.jl/njoy-reference/build/njoy < input > output 2>&1
```

Output: tape31 (2180), tape33 (604), tape34 (340), tape35 (211), tape36 (1711), tape37 (9274)

### Tape comparison (Julia GENDF → Julia dtfr)

| Tape | Julia lines | Fortran lines | Match % | Root cause of diffs |
|------|-------------|---------------|---------|---------------------|
| tape34 | 340 | 340 | 34.7% | gaminr MT621 wrong values |
| tape35 | 131 | 211 | ~30% | matxsr incomplete (diagonal scatter) |
| tape36 | 1627 | 1711 | 52.0% | gaminr + dpend thinning |
| tape37 | 9038 | 9274 | ~50% | downstream of tape36 |

### Tape comparison (Fortran GENDF → Julia dtfr)

| Tape | Match % | Notes |
|------|---------|-------|
| tape34 | **100%** | BIT-IDENTICAL |
| tape36 | **72.9%** | First 1182 lines identical; uranium PENDF thinning differs |

## Remaining Work — Priority Order

### Priority 1: Fix gaminr MT621 (photon heating KERMA)

**The bug**: Julia's gaminr writes raw cross sections (~0.6 barns for H group 1) instead of heating KERMA (~4488 barns·eV). The Fortran multiplies by average energy deposited per interaction.

**Evidence**:
```
Fortran tape33 MF23/MT621 group 1: flux=33300.9, sigma=4488.4
Julia   tape33 MF23/MT621 group 1: flux=33301.4, sigma=0.631
```

**Where to look**: `src/processing/gaminr.jl` — the function that computes MT621 (and possibly MT602 which is pair production heating). Compare with Fortran `gaminr.f90` subroutine `gheat` or equivalent. The heating KERMA = sum over partial reactions of (XS × average deposited energy). For photons: photoelectric deposits nearly all energy, Compton deposits partial energy, pair production deposits (E - 2*m_e*c²).

**Impact**: Fixes tape34 (all pheat values), tape36 (pheat histogram + y-axis range), and cascades to tape37.

### Priority 2: Fix matxsr transfer matrices

**The bug**: matxsr writes diagonal-only scatter (131 lines). Fortran writes full banded transfer matrices (211 lines).

**Where to look**: `src/orchestration/modules/matxsr.jl` lines 236-270. Currently writes `jband=1` (diagonal only). Need to compute actual bandwidth from `sec.transfer` array and write full banded matrix data.

**Key**: For each source group `ig`, scan `sec.transfer[1, ig, :]` to find the range of nonzero sink groups. Write `jband[ig]` = number of nonzero groups, `ijj[ig]` = lowest sink group. Then write all values in the band.

### Priority 3: Fix `_dpend` thinning (84 missing points)

**The bug**: Julia's `_dpend` implements only ns states 0-2 of Fortran's 5-state thinning machine. This produces 451 instead of 535 points for uranium MT=501 PENDF overlay.

**Where to look**: `src/orchestration/modules/dtfr.jl` lines 800-818. Compare with Fortran `dtfr.f90` lines 1241-1268 (subroutine `dpend`). Need to implement states ns=3 (keep both), ns=4 (always keep), ns=5 (triangle sort for local extrema preservation).

**Impact**: 84 extra lines in tape36 → matches Fortran exactly for uranium.

### Priority 4: Fix gaminr GENDF line count (688 vs 604)

Julia's gaminr writes 688 lines vs Fortran's 604. The extra lines come from:
1. Dense MF26 format (writing zeros for all sink groups instead of sparse above-threshold)
2. Section ordering (Julia: all MF23 then all MF26; Fortran: interleaved per MT)
3. Extra groups written for MT=516 pair production below threshold

These structural differences don't affect dtfr correctness (dtfr reads all sections regardless of order) but do affect matxsr and any line-count validation.

### Priority 5: Fix gaminr integration accuracy

Even for MT=501 (total), Julia values differ at ~0.002% from Fortran (trapezoidal vs Gauss-Lobatto quadrature). This is below 1e-5 tolerance but prevents bit-identical output.

## Files Modified This Session

- `src/orchestration/modules/dtfr.jl` — GENDF reader rewrite, DTF writer, plot tape writer, _dpend fix, _fmt_e10_2 helper, GendfMaterial flux field
- `src/orchestration/modules/matxsr.jl` — API migration to new GendfMF23Section/GendfMF26Section

## Traps

**Trap 144 (NEW — FIXED)**: `_dpend` hangs forever when PENDF data's first energy point equals `elow`. The `gety1_interp` closure returns `enext = energies_raw[1] = e`, creating an infinite loop. Fix: early return must return `energies_raw[2]` as `enext`, and while loop must use `<=` not `<`.

**Trap 145 (NEW)**: Julia's `@printf "%10.2E"` produces `1.20E+04` (mantissa ≥ 1). Fortran's `E10.2` (without `1p`) produces `0.12E+05` (mantissa < 1). Use `_fmt_e10_2()` helper for viewr tag values. Axis limits use `1p` format which matches Julia's standard `%E`.

**Trap 146 (NEW)**: matxsr references `gmat.sections` — old API. After the GENDF struct change, use `gmat.mf23` for cross sections, `gmat.mf26` for transfer matrices, `gmat.flux` for flux weights. The old `sec.data[g]` tuple is gone; use `sec.sigma[g]` and `gmat.flux[g]`.

**Trap 147 (NEW)**: When testing dtfr in isolation, use Fortran GENDF (tape33) as input to separate dtfr bugs from gaminr bugs. Julia's gaminr produces wrong MT621 values, so comparing Julia-pipeline tape34 against Fortran tape34 shows ~35% match — but this is a gaminr bug, not dtfr.
