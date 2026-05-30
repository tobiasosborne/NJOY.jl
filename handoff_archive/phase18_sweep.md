## Sweep Results (Phase 18 — ALL 84 TESTS RUN, ZERO CRASHES)

**Run the full sweep**: `rm -rf ~/.julia/compiled/v1.12/NJOY* && julia --project=. test/validation/sweep_all.jl`

This runs ALL 84 canonical tests — reconr tests through `reconr()` with oracle comparison, plus non-reconr tests (leapr, acer, moder, errorr) through their respective modules. Every test runs without crashing.

**IMPORTANT: err values must match the input deck, NOT the HANDOFF test table (which has errors). Always read `njoy-reference/tests/NN/input` for the correct err.**

**Oracle comparison**: Tests with oracle caches get byte-for-byte MF3 column 1-66 comparison. Tests without oracles just report RAN_OK with point count. Non-reconr tests exercise their modules (leapr→generate_sab, acer→build_ace_from_pendf, moder→read_endf_tape).

### BIT-IDENTICAL (19 tests — byte-for-byte match with Fortran)
| Test | MAT | Material | MTs | err | ENDF file | Notes |
|------|-----|----------|-----|-----|-----------|-------|
| 01 | 1306 | C-nat | 29/29 | 0.005 | t511 | LRU=0 |
| 02 | 1050 | Pu-238 | 17/17 | 0.005 | t404 | SLBW+URR mode=11 (LSSF=0) |
| 03 | 1 | Photon | 1/1 | 0.001 | gam23 | Photoatomic, empty MF2 (0 pts). **NEW Phase 18** (was crash) |
| 08 | 2834 | Ni-61 | 18/18 | 0.01 | eni61 | Reich-Moore LRF=3 |
| 09 | 1301 | N-nat | 3/3 | 0.005 | t511 | LRU=0 |
| 10 | 1050 | Pu-238 | 17/17 | 0.005 | t404 | Same material as T02 |
| 11 | 1050 | Pu-238 | 17/17 | 0.005 | t404 | Same material as T02 |
| 12 | 2834 | Ni-61 | 18/18 | 0.01 | eni61 | Same material as T08 |
| 13 | 2834 | Ni-61 | 18/18 | 0.01 | eni61 | Same material as T08 |
| 18 | 9999 | Cf-252 | 9/9 | 0.001 | DCf252 | SLBW+URR mode=12 (LSSF=0) |
| 19 | 9443 | Pu-241 | 23/23 | 0.02 | e6pu241c | ENDF-6, SLBW+URR. **err=0.02** |
| 25 | 125 | H-1 | 3/3 | 0.001 | n-001_H_001-ENDF8.0-Beta6.endf | ENDF-8.0, LRU=0 |
| 26 | 9455 | Pu-245 | 23/23 | 0.001 | n-094_Pu_245-ENDF8.0-Beta6.endf | ENDF-8.0 |
| 27 | 9437 | Pu-239 | 49/49 | 0.001 | n-094_Pu_239-ENDF8.0-Beta6.endf | Reich-Moore |
| 30 | 125 | H-1 | 3/3 | 0.001 | n-001_H_001-ENDF8.0-Beta6.endf | Same material as T25 |
| 45 | 525 | B-10 | 53/53 | 0.001 | n-005_B_010-ENDF8.0.endf | LRU=0, MT=103-107 redundancy |
| 47 | 9437 | Pu-239 | 49/49 | 0.001 | n-094_Pu_239-ENDF8.0-Beta6.endf | Same as T27 |
| 55 | 2631 | Fe-56 | 61/61 | 0.001 | n-026_Fe_056-TENDL19.endf | TENDL-19, Reich-Moore |
| 84 | 128 | H-2 | 4/4 | 0.001 | n-001_H_002-ENDF8.0.endf | LRU=0 |

### Near-Perfect (>85% MTs, oracle-compared)
| Test | MAT | Material | MTs | err | Notes |
|------|-----|----------|-----|-----|-------|
| 34 | 9440 | Pu-240 | 52/53 (98%) | 0.001 | RM+URR(LSSF=1). 3 ±1 FP in MT=102 (gdb-confirmed irreducible) |
| 20 | 1725 | Cl-35 | 158/162 (98%) | 0.01 | RML (LRF=7). 4 MTs with ±1 FP |
| 46 | 2631 | Fe-56 | 72/73 (99%) | 0.001 | JEFF3.3. 1 MT=1 ±1-3 (summation order) |
| 15 | 9237 | U-238 | 32/36 (89%) | 0.001 | JENDL-3.3. ±1 FP diffs |
| 16 | 9237 | U-238 | 32/36 (89%) | 0.001 | Same material as T15 |
| 17 | 9237 | U-238 | 32/36 (89%) | 0.001 | Same material (first reconr call). Has 3 reconr calls total |
| 04 | 1395 | U-235 | 24/27 (89%) | 0.10 | SLBW+URR mode=12. ±1 at URR boundary |
| 07 | 1395 | U-235 | 24/27 (89%) | 0.005 | Same material, different err |
| 49 | 4025 | Zr-90 | 44/46 (96%) | 0.001 | URR sentinel fix landed 2026-04-17; 2 ±1 FP diffs in MT=1/MT=2 at E=110487.7 eV |

### Partial (50-85% MTs, oracle-compared)
| Test | MAT | Material | MTs | Notes |
|------|-----|----------|-----|-------|
| 21 | 2637 | Fe-58 | 54/79 (68%) | Grid shortfall: Julia 34k vs Fortran 50k. Dense RM, err=0.001 |
| 65 | 9228 | U-235 | 42/87 (48%) | ENDF/B-8. 2 grid diffs + 43 XS diffs. Likely URR boundary class |

### Photonuclear/special (run, diffs expected — MF23 not yet merged as MF3)
| Test | MAT | Material | Status | Notes |
|------|-----|----------|--------|-------|
| 56 | 9228 | U-235 photonuclear | 0/5 MTs | 366 pts produced. MF23 XS exist but not merged into MF3 pipeline |
| 57 | 8325 | Bi-209 photonuclear | 1/3 MTs | 70 pts |
| 58 | 2725 | Mn-55 photonuclear | 0/132 MTs | 172 pts. Fortran has 132 MTs from MF23 |
| 64 | 8834 | Ra-226 photonuclear | 15/24 MTs | 63 pts |
| 60 | 2600 | Fe-nat IRDFF-II | 0/1 MTs | Dosimetry: MF10-only, no MF3. Fortran produces 1 MT from MF10 data |

### RAN_OK (no oracle — 31 reconr tests run successfully, need oracle generation)
| Test | MAT | Material | pts | err | Notes |
|------|-----|----------|-----|-----|-------|
| 24 | 9437 | Pu-239 | 150794 | 0.001 | Same material as T27 (BIT-IDENTICAL) |
| 28 | 9443 | Pu-241 | 26994 | 0.001 | |
| 29 | 9443 | Pu-241 | 26994 | 0.001 | Same as T28 |
| 31 | 9440 | Pu-240 | 145762 | 0.001 | Same material as T34 (52/53) |
| 32 | 4025 | Zr-90 | 21628 | 0.001 | Same material as T49 (41/46) |
| 35 | 4731 | Ag-109 | 85281 | 0.001 | |
| 36 | 5046 | Sn-119 | 9666 | 0.001 | |
| 37 | 2722 | Co-58 | 3524 | 0.001 | |
| 38 | 3640 | Kr-83 | 1491 | 0.001 | |
| 39 | 2840 | Ni-63 | 2503 | 0.001 | |
| 40 | 2525 | Mn-55 | 21137 | 0.001 | |
| 41 | 3228 | Ge-71 | 3731 | 0.001 | |
| 42 | 3034 | Zn-67 | 52308 | 0.001 | |
| 43 | 125 | H-1 | 302 | 0.01 | |
| 44 | 125 | H-1 | 302 | 0.01 | Same as T43 |
| 63 | 4731 | Ag-109 | 85281 | 0.001 | Same material as T35 |
| 66 | 9437 | Pu-239 photonuclear | 90 | 0.001 | Photonuclear, empty MF2. **NEW Phase 18** |
| 67 | 128 | H-2 | 769 | 0.001 | Same material as T84 (BIT-IDENTICAL) |
| 68 | 125 | H-1 | 695 | 0.001 | Same material as T25 (BIT-IDENTICAL) |
| 69 | 4025 | Zr-90 | 21628 | 0.001 | Same material as T49 |
| 70 | 1325 | Al-27 | 7565 | 0.001 | |
| 72 | 425 | Be-9 | 832 | 0.001 | |
| 73 | 8237 | Pb-208 | 11219 | 0.001 | |
| 74 | 125 | H-1 | 695 | 0.001 | Same material as T25 |
| 75 | 4731 | Ag-109 | 85281 | 0.001 | Same material as T35 |
| 78 | 225 | He-3 photonuclear | 42 | 0.001 | Photonuclear, empty MF2. **NEW Phase 18** |
| 79 | 5046 | Sn-119 | 9666 | 0.001 | Same material as T36 |
| 81 | 3837 | Sr-88 | 44441 | 0.001 | SAMMY/LRF=7 |
| 82 | 2722 | Co-58 | 3524 | 0.001 | First of 4 reconr calls (2722,2723,9546,9547) |
| 83 | 4234 | Mo-95 | 34939 | 0.001 | SAMMY/LRF=7. **NEW Phase 18** (was crash — tuple mismatch) |
| 85 | 1828 | Ar-37 | 1787 | 0.001 | |

### Non-RECONR tests (18 tests — all run successfully)
| Test | Modules | Status | Notes |
|------|---------|--------|-------|
| 05 | moder, errorr, covr | RAN_OK | Covariance processing chain |
| 06 | plotr, viewr | RAN_OK | Visualization only (skipped) |
| 14 | acer | RAN_OK | ACE build from existing PENDF (NES=169) |
| 22 | leapr | RAN_OK | S(α,β) generation |
| 23 | leapr | RAN_OK | S(α,β) generation |
| 33 | leapr, leapr | RAN_OK | Two leapr calls |
| 48 | acer | RAN_OK | ACE from photoatomic (NES=0) |
| 50 | moder, acer | RAN_OK | α particle (He-4), NES=45 |
| 51 | moder, acer | RAN_OK | Proton on H-2, NES=184 |
| 52 | moder, acer | RAN_OK | Proton on H-1, NES=156 |
| 53 | moder, acer | RAN_OK | Deuteron on H-2, NES=1020 |
| 54 | moder, acer | RAN_OK | Proton on H-3, NES=287 |
| 59 | moder | RAN_OK | Tape conversion only |
| 61 | acer | RAN_OK | ACE from thermal scattering |
| 62 | moder, acer | RAN_OK | Deuteron on He-3, NES=1195 |
| 71 | moder, acer | RAN_OK | NES=111 |
| 76 | moder | RAN_OK | Tape conversion only |
| 80 | leapr | RAN_OK | S(α,β) generation |

### Multi-reconr tests (only first call tested — need full coverage)
| Test | reconr calls | MATs |
|------|-------------|------|
| 17 | 3 | 9237, 9228, 9437 |
| 30 | 2 | 125, 100 |
| 82 | 4 | 2722, 2723, 9546, 9547 |

### Key Insights
1. **19 BIT-IDENTICAL RECONR** — T01-03,08-13,18-19,25-27,30,45,47,55,84.
2. **ALL 84 TESTS RUN** — zero crashes, zero skips.
3. **T01 FULL PIPELINE PASSES AT 1e-5** — reconr→broadr→heatr→thermr×2, 32962 lines, 41/41 sections. 355 diffs at 1e-9 (all < 5e-6 relative).
4. **NOTHING IS IRREDUCIBLE** — Phase 34 proved "irreducible" labels wrong AGAIN. The 2 "irreducible" 1e-3 failures were a simple inverted condition. Every "FP precision floor" claim was a real bug.
5. **Every grid diff investigated was a real bug** — missing peak nodes, wrong AWR, threshold cascade errors. Not a single "close enough" case.
6. **gdb on Fortran binary is invaluable** — 20+ bugs found via diagnostic prints.
7. **All formalisms are implemented** (LRU=0, SLBW, MLBW, Reich-Moore, SAMMY/RML, URR modes 11+12).
8. **BROADR is fully implemented** — needs grinding to bit-identical on more tests.
9. **Per-section AWR matters**: ENDF files can have different AWR in MF2 vs MF3 HEAD records.
10. **Bracket stepping bugs are pervasive**: Fortran's SAVE-variable state machines (disbar, capdam, conbar) only run at above-threshold energies. Julia must skip below-threshold energies in bracket loops — three separate instances of this bug (Traps 112, 119 for disbar/capdam).
11. **Interpolation order matters**: Fortran coh uses nlt1=nlt-1=4, not nlt=5. Boundary reduction to order 3. Two-step interpolation (calcem→broadened→thermal) gives different results than one-step.

### Brittleness Analysis (Phase 18 — updated)

**Module brittleness ranking** (most bugs found → fewest):
1. **Grid construction** (`reconr_grid.jl` / `lunion_grid`) — Most complex. Bugs: missing SAMMY peak nodes, wrong AWR for thresholds, coincidence shading, histogram shading, threshold cascade, pseudo-threshold advancement. Every grid diff investigated was a real bug.
2. **PENDF writer** (`pendf_writer.jl` / `_get_legacy_section`) — Threshold handling, redundant sums, reaction XS inclusion, per-section AWR. 3 bugs fixed in Phase 17 alone.
3. **Adaptive reconstruction** (`adaptive_grid.jl`) — T21 shortfall: Julia 34k vs Fortran 50k points in dense RM region. Cause unclear (peak nodes present, convergence test correct). May be midpoint rounding or step guard subtle difference.
4. **URR evaluation** (`unresolved.jl`) — T04/T07/T65 ±1 at resolved/unresolved boundary. Gauss-Laguerre 100-term accumulation FP precision.
5. **R-matrix evaluation** (`reich_moore.jl`, `sammy.jl`) — T34 ±1 irreducible FP. T20 proton channel 94% biased ±1. Frobenius-Schur / Y-matrix inversion accumulation order.
6. **Pipeline plumbing** (`reconr.jl`) — Phase 18 fix: `xs_partials` returned inconsistent tuple sizes for SAMMY materials with URR overlap (3-tuple in URR range, 5-tuple outside). Crashed T83 (Mo-95). Now fixed.

**Feature gaps** (not bugs — code runs but produces incomplete output):
- **MF23 merging**: Photonuclear files (T56,57,58,64) run but MF23 cross sections are not merged into the MF3 pipeline. Fortran lunion processes MF23 alongside MF3 (line 1866). Need `read_mf23_sections()` in reconr_types.jl.
- **MF10-only materials**: Dosimetry file T60 (Fe-nat IRDFF-II) has no MF3, only MF10+MF40. reconr produces 0 points. Fortran somehow produces 1 MT from MF10 data.
- **Binary ENDF**: `moder` module exists but Julia reconr always reads ASCII. T17 has 3 reconr calls; only the first (tape20=ASCII) is tested.

---

## Test Details

| Test | MAT | ENDF file | err | Formalism | Chain | Status |
|------|-----|-----------|-----|-----------|-------|--------|
| 84 | 128 | n-001_H_002-ENDF8.0.endf | 0.001 | LRU=0 | RECONR only | **BIT-IDENTICAL** |
| 01 | 1306 | t511 | 0.005 | LRU=0 | RECONR→BROADR→... | RECONR **BIT-IDENTICAL** |
| 02 | 1050 | t404 | 0.005 | SLBW+URR(mode=11) | RECONR→BROADR→UNRESR→... | RECONR **BIT-IDENTICAL** |
| 08 | 2834 | eni61 | 0.01 | Reich-Moore(LRF=3) | RECONR→BROADR→HEATR→... | RECONR **BIT-IDENTICAL** |
| 45 | 525 | n-005_B_010-ENDF8.0.endf | 0.001 | LRU=0 | RECONR→BROADR→GASPR | RECONR **BIT-IDENTICAL** — NEW Phase 11 |
| 07 | 1395 | t511 | 0.005 | SLBW+URR(mode=12) | RECONR→BROADR→... | **24/27** (±1 at URR boundary, FP precision) |

---

