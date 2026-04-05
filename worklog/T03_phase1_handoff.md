# T03 Pipeline — Phase 1 Handoff

## Date: 2026-04-05

## What T03 Is

T03 is a **photon interaction library generation** test.
- Chain: reconr → gaminr → dtfr → matxsr → viewr
- Input: tape30 (gam23 = DLC-7E photon XS library), tape32 (gam27 = MF27 form factors)
- Materials: MAT=1 (hydrogen), MAT=92 (uranium)
- Official test compares **tape37 only** (viewr PostScript output) at rel_tol=1e-9
- Fortran oracle generated at `/tmp/t03_fortran/` with all intermediate tapes

## What Was Done This Session

### 1. Full codebase audit (Rule 4: be skeptical)

Discovered that the processing/format layer for ALL T03 modules already exists:
- `src/processing/gaminr.jl` — group averaging, 9 photon group structures, weight functions
- `src/formats/dtfr.jl` — DTFMaterial, DTFLayout, write_dtf
- `src/formats/matxsr.jl` — write_matxs ASCII+binary
- `src/visualization/backends.jl` — render_postscript

The gap was entirely in the **orchestration layer**: no parsers, no `*_module` functions, no pipeline dispatch.

### 2. Input parsers added to `src/orchestration/input_parser.jl`

| Struct | Parser | Cards parsed |
|--------|--------|-------------|
| `GaminrParams` | `parse_gaminr` | Card 1 (units), Card 2 (matb,igg,iwt,lord,iprint), Card 3 (title), Card 6-7 (material loop) |
| `DtfrParams` + `DtfrEditSpec` + `DtfrMaterial` | `parse_dtfr` | Card 1-8 (iedit=0 branch, edit specs, material list) |
| `MatxsrParams` + `MatxsrMaterial` | `parse_matxsr` | Card 1-10 (all MATXS control cards + material list) |
| `ViewrParams` | `parse_viewr` | Card 1 (infile, nps) |

Updated `AVAILABLE_MODULES` to include `:dtfr`, `:matxsr`, `:viewr`.
Removed `:viewr` from `SKIP_MODULES`.

### 3. Four orchestration modules created

| File | Status | Key functions |
|------|--------|---------------|
| `src/orchestration/modules/gaminr.jl` | **Substantially working** | `gaminr_module`, `_write_gaminr_tape`, `_gaminr_group_average`, `_gaminr_wt1_eval`, `_read_pendf_mf23`, `_read_mf27_form_factors`, coherent/incoherent/pair scattering matrices, MT=621 heating |
| `src/orchestration/modules/dtfr.jl` | **Structural stub** | `dtfr_module`, `_read_gendf_tape`, `_write_dtf_output`, `_write_dtfr_plot_tape` |
| `src/orchestration/modules/matxsr.jl` | **Structural stub** | `matxsr_module`, `_write_matxsr_tape` |
| `src/orchestration/modules/viewr.jl` | **Structural stub** | `viewr_module`, `_render_njoy_plot_tape`, plot tape parser, PS drawing primitives |

### 4. Pipeline dispatch wired in `src/orchestration/pipeline.jl`

Added `elseif` branches for `:gaminr`, `:dtfr`, `:matxsr`, `:viewr`.

### 5. Pipeline runs end-to-end without crashing

```
run_njoy("njoy-reference/tests/03/input")
```
Produces all 5 tapes (33, 34, 35, 36, 37).

## Critical Discovery: iwt=3 Weight Function

**iwt=3 is NOT 1/E.** It's a 4-point tabulated weight function with log-log (INT=5) interpolation.

From Fortran `gnwtf` (gaminr.f90:795-798), the `wt1` array encodes a TAB1:
```
NR=1, NP=4, NBT=4, INT=5 (log-log)
Data: (1e3 eV, 1e-4), (1e5 eV, 1.0), (1e7 eV, 1e-2), (3e7 eV, 1e-4)
```

For group 1 [1e4, 1e5]: W(E) = 1e-4 * (E/1e3)^2 (log-log interpolation).
- Analytical flux integral: ∫W dE = 1e-10 * (1e15 - 1e12)/3 ≈ 33333
- Fortran oracle: 33300.8990
- Julia (after fix): 33301.3992

Before the fix, Julia used simple 1/E giving flux = 2.3 — a 14,000x error.

Implemented as `_gaminr_wt1_eval()` with `_GAMINR_WT1_E` and `_GAMINR_WT1_W` constants.

## Current tape33 Comparison (gaminr GENDF)

### Structure: MATCH
Julia and Fortran produce the same MF/MT sections:
- MF1/MT451 (header) — both materials
- MF23/MT501,502,504,516,602,621 — both materials
- MF26/MT502,504,516 — both materials

### Line count: 688 Julia vs 604 Fortran (+84 lines)

Per-section differences:
| Section | Julia | Fortran | Issue |
|---------|-------|---------|-------|
| MF23/MT516 | 25 | 19 | Julia writes all 12 groups; Fortran only writes above-threshold groups |
| MF26/MT502 | 49 | 37 | Julia writes dense; Fortran uses sparse format (omits zero sink groups) |
| MF26/MT504 | 81 | 85 | Julia sparse format not quite matching Fortran's record packing |
| MF26/MT516 | 45 | 19 | Julia writes too many groups; needs threshold handling |

### Values: ~0.002% accuracy for MF23

| Group | Julia flux | Fortran flux | Ratio |
|-------|-----------|-------------|-------|
| 1 | 33301.3992 | 33300.8990 | 1.000015 |
| 2 | 160951.784 | 160951.084 | 1.000004 |
| 3 | 69318.8269 | 69318.1266 | 1.000010 |

The ~0.002% gap is from trapezoidal integration vs Fortran's 10-point Gauss-Lobatto quadrature.

### Section ordering: WRONG
Julia writes all MF23 then all MF26. Fortran interleaves: MF23/MT502 → MF26/MT502 → MF23/MT504 → MF26/MT504.

## Remaining Work for tape33 (Priority Order)

1. **Fix section ordering** — interleave MF23/MF26 per MT matching Fortran
2. **Fix MF26 sparse format** — omit zero sink groups, match Fortran CONT record ng2/iglo
3. **Fix MT=516 threshold** — only write groups above 1.022 MeV pair threshold
4. **Implement Lobatto quadrature** — replace trapezoidal with 10-point Gauss-Lobatto to match Fortran gpanel
5. **Match TPID** — Julia reconr writes "PENDF tape produced by NJOY.jl" vs Fortran's material-specific title

## Remaining Work for Full T03 Pipeline

### tape34 (dtfr DTF tables) — 340 vs 340 lines
Same line count but wrong values. Needs:
- Read GENDF tape33 properly (including MF26 scattering matrices)
- Fill DTF table columns correctly: pheat (MT=621), absorption, total (MT=501), in-group, scatter bands
- Match Fortran `e12.4` format exactly

### tape35 (matxsr MATXS) — 175 vs 211 lines
Needs:
- Read MF26 scattering matrices from tape33
- Write MATXS matrix blocks ("8d"/"9d" records)
- Match Fortran record format exactly

### tape36 (dtfr plot tape) — 10528 vs 1711 lines
Julia produces wrong plot set. Needs:
- Correct plot set: pheat, absorp, total, phot-phot 3D table (per material)
- NOT one plot per MF23 MT
- Read Fortran `dtfplt` subroutine to match exact card format
- PENDF overlay curves on 2D plots

### tape37 (viewr PostScript) — 10663 vs 9274 lines
Completely dependent on tape36. Also needs:
- Full PostScript renderer matching Fortran graph.f90
- Exact coordinate transforms (landscape, 576x756 pt page)
- Tick marks, axis labels, fonts, 3D scatter matrix view
- Field widths: Fortran uses `f9.2` for coordinates (9 chars, 2 decimal)

## Physical Constants (from gaminr.f90:1198-1202)

```
c1 = 57.03156e-6    # momentum transfer ↔ angle conversion
c2 = 0.249467       # r_e^2/2 (Thomson XS prefactor)
c3 = 1.95693e-6     # E/(m_e c^2) conversion
c4 = 0.0485262      # related constant
c5 = 20.60744       # related constant
epair = 0.511e6 eV  # electron rest mass
```

## Fortran Oracle Location

All Fortran reference tapes at `/tmp/t03_fortran/tape{31,33,34,35,36,37}`.
Generated by running: `cd /tmp/t03_fortran && njoy < input` with gam23→tape30, gam27→tape32.

## Files Changed

| File | Change |
|------|--------|
| `src/orchestration/input_parser.jl` | Added GaminrParams, DtfrParams, MatxsrParams, ViewrParams + parsers. Updated AVAILABLE_MODULES. |
| `src/orchestration/pipeline.jl` | Added dispatch for :gaminr, :dtfr, :matxsr, :viewr |
| `src/orchestration/modules/gaminr.jl` | **NEW** — full gaminr module with MF23 averaging, MF26 scattering, MT=621 heating, form factor reader |
| `src/orchestration/modules/dtfr.jl` | **NEW** — structural stub, GENDF reader, DTF writer, plot tape writer |
| `src/orchestration/modules/matxsr.jl` | **NEW** — structural stub, MATXS writer |
| `src/orchestration/modules/viewr.jl` | **NEW** — structural stub, plot tape parser, PS renderer skeleton |
| `src/NJOY.jl` | Added includes for 4 new module files |

## How to Run

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/03/input"; work_dir="/tmp/t03_test")'
```

## How to Compare Against Oracle

```bash
# Generate oracle (already done, at /tmp/t03_fortran/)
cd /tmp && mkdir -p t03_fortran && cd t03_fortran
cp njoy-reference/tests/resources/gam23 tape30
cp njoy-reference/tests/resources/gam27 tape32
njoy-reference/build/njoy < njoy-reference/tests/03/input

# Compare
diff /tmp/t03_test/tape33 /tmp/t03_fortran/tape33 | head -40
```

## Key Rules Followed

1. Read the Fortran before writing Julia — read gaminr.f90 gnwtf, gpanel, gtff, genggp
2. Fortran is ground truth — matched iwt=3 weight function to wt1 TAB1 encoding
3. Be skeptical — discovered existing Julia processing layer had wrong 1/E weight for iwt=3
4. One Julia process at a time — all subagents were read-only research
5. Grind method — diff → read Fortran → fix → retest → repeat
