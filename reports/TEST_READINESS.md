# NJOY2016 Test Suite Readiness for NJOY.jl

Generated: 2026-03-21

## Summary

| Status    | Count | Description                                          |
|-----------|-------|------------------------------------------------------|
| READY     | 24    | All modules available in NJOY.jl                     |
| PARTIAL   | 51    | Some modules available (usually only `viewr` missing) |
| NOT_READY | 5     | Requires `leapr` or `plotr`/`viewr` only             |
| EMPTY     | 4     | No NJOY modules in input (tests 21, 48, 75, 76)      |
| SKIPPED   | 1     | Test 77 does not exist in the suite                  |

NJOY.jl available processing modules: `reconr`, `broadr`, `heatr`, `thermr`, `unresr`, `purr`, `groupr`, `errorr`, `moder`, `acer` (neutron mode 1 only).

Modules NOT implemented: `plotr`, `viewr`, `leapr`, `gaminr`, `dtfr`, `ccccr`, `matxsr`, `wimsr`, `gaspr`, `mixr`, `powr`, `resxsr`.

---

## Full Classification (Tests 01-85, skipping 77)

### READY -- All modules available (24 tests)

| Test | Modules Used | ENDF Input | Nuclide | Description |
|------|-------------|------------|---------|-------------|
| 01 | moder, reconr, broadr, heatr, thermr, groupr | t511 (tape20), t322 (tape26) | C-nat (MAT 1306) | Full pipeline: reconr+broadr+heatr+thermr+groupr |
| 04 | moder, reconr, errorr, groupr | t511 (tape20) | U-235 (MAT 1395) | ERRORR covariance + GROUPR nubar |
| 07 | moder, reconr, broadr, heatr, groupr, acer | t511 (tape20) | U-235 (MAT 1395) | Full pipeline through ACER |
| 08 | moder, reconr, broadr, heatr, groupr, acer | eni61 (tape20) | Ni-61 (MAT 2834) | HEATR with 6 MT KERMA + GROUPR + ACER |
| 10 | moder, reconr, broadr, unresr, purr, acer | t404 (tape20) | Pu-238 (MAT 1050) | UNRESR + PURR probability tables + ACER |
| 17 | reconr, broadr, groupr, moder, errorr | J33U238, J33U235, J33Pu239 | U-238/U-235/Pu-239 | Multi-material GROUPR + cross-material ERRORR |
| 19 | moder, reconr, broadr, unresr, heatr, purr, acer | e6pu241c (tape20) | Pu-241 (MAT 9443) | Full actinide pipeline with URR |
| 25 | moder, reconr, broadr, heatr, thermr, acer | H-1 ENDF8+tsl-HinH2O | H-1 (MAT 125) | Thermal scattering ACER for H in H2O, 3 temps |
| 26 | moder, reconr, broadr, heatr | Pu-245 ENDF8 (tape20) | Pu-245 (MAT 9455) | HEATR without MF1/MT458 |
| 29 | moder, reconr, broadr, groupr | Pu-241 ENDF8 (tape20) | Pu-241 (MAT 9443) | GROUPR delayed chi (3 runs with different options) |
| 31 | moder, reconr, broadr, purr | Pu-240 ENDF8 (tape20) | Pu-240 (MAT 9440) | PURR probability tables |
| 32 | moder, reconr, broadr, thermr | Zr-90 ENDF8+tsl-ZrinZrH | Zr-90 (MAT 4025) | THERMR with thermal scattering law |
| 43 | moder, reconr, broadr | H-1 ENDF8 (tape20) | H-1 (MAT 125) | BROADR at T=0 K edge case |
| 44 | moder, reconr, broadr | H-1 ENDF8 (tape20) | H-1 (MAT 125) | BROADR at T=0.1 K edge case |
| 46 | moder, reconr, broadr, groupr, errorr | Fe-56 JEFF3.3 (tape20) | Fe-56 (MAT 2631) | 238-group GROUPR + ERRORR covariance |
| 49 | moder, reconr, broadr, heatr, thermr, acer | Zr-90 ENDF8+tsl-ZrinZrH | Zr-90 (MAT 4025) | Thermal ACER for Zr in ZrH |
| 60 | moder, reconr, broadr, groupr | Fe-nat IRDFF-II (tape20) | Fe-nat (MAT 2600) | GROUPR with 140-group custom structure (x2) |
| 73 | moder, reconr, broadr, groupr | Pb-208 TENDL2021 (tape20) | Pb-208 (MAT 8237) | GROUPR with custom 2-group structure |
| 79 | moder, reconr, broadr, heatr | Sn-119 ENDF8 (tape20) | Sn-119 (MAT 5046) | HEATR with local photon option |
| 81 | moder, reconr | Sr-88 ENDF8.1 (tape20) | Sr-88 (MAT 3837) | RECONR with resolved resonances |
| 82 | reconr, broadr, acer | Co-58, Co-58m1, Am-242, Am-242m1 ENDF8 | Multiple (MAT 2722/2723/9546/9547) | ACER izaopt tests for 4 nuclides |
| 83 | moder, reconr | Mo-95 beta (tape20) | Mo-95 (MAT 4234) | RECONR for beta evaluation |
| 84 | reconr | H-2 ENDF8 (tape20) | H-2 (MAT 128) | RECONR only, no resonances |
| 85 | moder, reconr | Ar-37 TENDL2023 (tape20) | Ar-37 (MAT 1828) | RECONR for TENDL evaluation |

### PARTIAL -- Some modules available (51 tests)

| Test | Modules Used | Missing | Notes |
|------|-------------|---------|-------|
| 02 | moder, reconr, broadr, unresr, groupr, ccccr | ccccr | CCCCR output format |
| 03 | reconr, gaminr, dtfr, matxsr, viewr | gaminr, dtfr, matxsr, viewr | Photon + transport |
| 05 | moder, errorr, viewr | viewr | ERRORR + VIEWR plot |
| 09 | moder, reconr, broadr, leapr, thermr | leapr | LEAPR thermal law generation |
| 11 | moder, reconr, broadr, unresr, thermr, groupr, wimsr | wimsr | WIMS library format |
| 12 | reconr, gaspr, plotr, viewr | gaspr, plotr, viewr | Gas production + plot |
| 13 | moder, reconr, broadr, heatr, gaspr, acer, viewr | gaspr, viewr | Gas production + plot |
| 14 | acer, viewr | viewr | Proton ACER + VIEWR |
| 15 | moder, reconr, broadr, groupr, errorr, viewr | viewr | ERRORR covariance plot |
| 16 | moder, reconr, broadr, errorr, viewr | viewr | ERRORR pointwise plot |
| 18 | moder, reconr, broadr, groupr, errorr, viewr | viewr | ERRORR MG covariance plot |
| 20 | errorr, reconr, broadr, groupr, viewr | viewr | ERRORR plot |
| 24 | moder, reconr, broadr, heatr, thermr, gaspr, acer, viewr | gaspr, viewr | Full pipeline + gas |
| 27 | moder, reconr, broadr, groupr, errorr, viewr | viewr | ERRORR covariance plot |
| 28 | moder, reconr, broadr, acer, viewr | viewr | ACER + VIEWR plot |
| 30 | moder, reconr, broadr, groupr, gaminr, matxsr | gaminr, matxsr | Photon multigroup + MATXS output |
| 34 | moder, reconr, broadr, groupr, errorr, viewr | viewr | ERRORR plot |
| 35 | moder, reconr, broadr, purr, acer, viewr | viewr | PURR lssf=0 + VIEWR |
| 36 | moder, reconr, broadr, purr, acer, viewr | viewr | PURR lssf=0 variants |
| 37 | moder, reconr, broadr, purr, acer, viewr | viewr | PURR lssf=0 variants |
| 38 | moder, reconr, broadr, purr, acer, viewr | viewr | PURR lssf=0 variants |
| 39 | moder, reconr, broadr, purr, acer, viewr | viewr | PURR lssf=1 |
| 40 | moder, reconr, broadr, purr, acer, viewr | viewr | PURR lssf=1 variants |
| 41 | moder, reconr, broadr, purr, acer, viewr | viewr | PURR lssf=1 variants |
| 42 | moder, reconr, broadr, purr, acer, viewr | viewr | PURR lssf=1 variants |
| 45 | moder, reconr, broadr, gaspr | gaspr | Gas production |
| 47 | moder, reconr, broadr, unresr, groupr, errorr, viewr | viewr | ERRORR plot |
| 50 | moder, acer, viewr | viewr | ACER deuteron MF6 MT2 |
| 51 | moder, acer, viewr | viewr | ACER deuteron MF6 MT2 |
| 52 | moder, acer, viewr | viewr | ACER deuteron MF6 MT2 |
| 53 | moder, acer, viewr | viewr | ACER deuteron MF6 MT2 |
| 54 | moder, acer, viewr | viewr | ACER deuteron MF6 MT2 |
| 55 | moder, reconr, broadr, acer, viewr | viewr | ACER + VIEWR |
| 56 | moder, reconr, acer, viewr | viewr | ACER + VIEWR |
| 57 | moder, reconr, acer, viewr | viewr | ACER + VIEWR |
| 58 | moder, reconr, acer, viewr | viewr | ACER + VIEWR |
| 59 | moder, acer(mode 4) | acer photoatomic mode | Photoatomic ACE (acer mode 4) |
| 61 | acer(mode 7) | acer check/rename mode | ACE check/rename (acer mode 7) |
| 62 | moder, acer, viewr | viewr | ACER + VIEWR |
| 63 | moder, reconr, broadr, purr, acer, viewr | viewr | PURR nunx=2 + VIEWR |
| 64 | moder, reconr, acer, viewr | viewr | ACER + VIEWR |
| 65 | moder, reconr, errorr, viewr | viewr | ERRORR + VIEWR |
| 66 | moder, reconr, acer, viewr | viewr | ACER + VIEWR |
| 67 | moder, reconr, broadr, thermr, acer, viewr | viewr | Thermal ACER + VIEWR |
| 68 | moder, reconr, broadr, thermr, acer, viewr | viewr | Thermal ACER + VIEWR |
| 69 | moder, reconr, broadr, thermr, acer, viewr | viewr | Thermal ACER + VIEWR |
| 70 | moder, reconr, broadr, thermr, acer, viewr | viewr | Thermal ACER + VIEWR |
| 71 | moder, acer, viewr | viewr | ACER deuteron + VIEWR |
| 72 | moder, reconr, broadr, heatr, thermr, gaspr, acer, viewr | gaspr, viewr | Full pipeline + gas |
| 74 | moder, reconr, broadr, thermr, acer, viewr | viewr | Thermal ACER + VIEWR |
| 78 | moder, reconr, acer, viewr | viewr | ACER + VIEWR |

### NOT_READY -- Missing critical modules (5 tests)

| Test | Modules Used | Missing | Notes |
|------|-------------|---------|-------|
| 06 | plotr, viewr | plotr, viewr | Pure plotting test |
| 22 | leapr | leapr | H in H2O thermal scattering law |
| 23 | leapr | leapr | D in D2O thermal scattering law |
| 33 | leapr | leapr | He-4 thermal scattering law |
| 80 | leapr | leapr | Thermal scattering law (Be) |

### EMPTY -- No processing modules (4 tests)

| Test | Notes |
|------|-------|
| 21 | Empty or comment-only input |
| 48 | Empty or comment-only input |
| 75 | Empty or comment-only input |
| 76 | Empty or comment-only input |

---

## Top 10 READY Tests -- Execution Plan

### 1. Test 84: RECONR Only -- H-2 (Simplest Case)

**Already implemented** in `test/integration_tests.jl`.

- **Pipeline**: `reconr` only
- **Input**: `resources/n-001_H_002-ENDF8.0.endf` (tape20)
- **MAT**: 128, tolerance 0.001
- **Reference**: `referenceTape100` (ASCII PENDF)
- **Comparison**: Parse referenceTape100 as PENDF, compare MF3/MT1,2,102 cross sections
- **Status**: Already passing. No resonances (LRU=0), simplest validation.

```julia
result = reconr(endf_file; mat=128, err=0.001)
# Compare result.total, result.elastic, result.capture against referenceTape100
```

### 2. Test 43: RECONR + BROADR -- H-1 at T=0 K

**Already implemented** in `test/integration_tests.jl`.

- **Pipeline**: `moder` -> `reconr` -> `broadr` (T=0.0 K)
- **Input**: `resources/n-001_H_001-ENDF8.0-Beta6.endf` (tape20)
- **MAT**: 125, tolerance 0.01
- **Reference**: `referenceTape35` (PENDF after broadening)
- **Action**: BROADR at T=0 is identity; validates pipeline composition.

```julia
result = reconr(endf_file; mat=125, err=0.01)
xs_matrix = hcat(result.total, result.elastic, result.fission, result.capture)
pendf = PointwiseMaterial(Int32(125), result.energies, xs_matrix, [1,2,18,102])
broadened = doppler_broaden(pendf, 0.0; T_old=0.0, awr=0.9992, tol=0.01)
# Compare broadened against referenceTape35
```

### 3. Test 81: MODER + RECONR -- Sr-88 (Resolved Resonances)

**Already implemented** in `test/integration_tests.jl`.

- **Pipeline**: `moder` -> `reconr`
- **Input**: `resources/n-038_Sr_088-ENDF8.1.endf` (tape20)
- **MAT**: 3837, tolerance 0.001
- **Reference**: `referenceTape30` (44,441 points expected)
- **Challenge**: RML (LRF=7) formalism. Current reconr uses fallback evaluator.

```julia
result = reconr(endf_file; mat=3837, err=0.001)
# Compare against referenceTape30 (44441-point PENDF)
```

### 4. Test 44: RECONR + BROADR -- H-1 at T=0.1 K

- **Pipeline**: `moder` -> `reconr` -> `broadr` (T=0.1 K)
- **Input**: `resources/n-001_H_001-ENDF8.0-Beta6.endf` (tape20)
- **MAT**: 125, tolerance 0.01, thnmax=-19500000.0
- **Reference**: `referenceTape35` (PENDF)
- **Action needed**: Same as Test 43 but with T=0.1 K broadening. Validates sigma1 kernel at near-zero temperature.

```julia
result = reconr(endf_file; mat=125, err=0.01)
xs_matrix = hcat(result.total, result.elastic, result.fission, result.capture)
pendf = PointwiseMaterial(Int32(125), result.energies, xs_matrix, [1,2,18,102])
broadened = doppler_broaden(pendf, 0.1; T_old=0.0, awr=0.9992, tol=0.01)
# Compare broadened against referenceTape35
```

### 5. Test 85: MODER + RECONR -- Ar-37 (TENDL)

- **Pipeline**: `moder` -> `reconr`
- **Input**: `resources/n-018_Ar_37-tendl2023.endf` (tape20)
- **MAT**: 1828, tolerance 0.001
- **Reference**: `referenceTape50` (ASCII PENDF)
- **Action needed**: Run reconr, write PENDF output, compare against reference.

```julia
result = reconr(endf_file; mat=1828, err=0.001)
# Write PENDF and compare against referenceTape50
```

### 6. Test 83: MODER + RECONR -- Mo-95 (Beta Evaluation)

- **Pipeline**: `moder` -> `reconr`
- **Input**: `resources/n-042_Mo_095-beta.endf` (tape20)
- **MAT**: 4234, tolerance 0.001
- **Reference**: `referenceTape50` (ASCII PENDF)
- **Action needed**: Validates reconr on a beta ENDF evaluation.

```julia
result = reconr(endf_file; mat=4234, err=0.001)
# Write PENDF and compare against referenceTape50
```

### 7. Test 26: MODER + RECONR + BROADR + HEATR -- Pu-245

- **Pipeline**: `moder` -> `reconr` -> `broadr` (293.6 K) -> `heatr` (7 MT KERMA)
- **Input**: `resources/n-094_Pu_245-ENDF8.0-Beta6.endf` (tape20)
- **MAT**: 9455, tolerance 0.001
- **Reference**: `referenceTape40` (PENDF with KERMA)
- **Action needed**: Full pipeline through HEATR. Tests HEATR without MF1/MT458 fission Q data.

```julia
result = reconr(endf_file; mat=9455, err=0.001)
xs_matrix = hcat(result.total, result.elastic, result.fission, result.capture)
pendf = PointwiseMaterial(Int32(9455), result.energies, xs_matrix, [1,2,18,102])
broadened = doppler_broaden(pendf, 293.6; T_old=0.0, awr=awr_pu245, tol=0.001)
kerma = compute_kerma(broadened; mts=[302,303,304,318,442,443,444])
# Compare against referenceTape40
```

### 8. Test 79: MODER + RECONR + BROADR + HEATR -- Sn-119

- **Pipeline**: `moder` -> `reconr` -> `broadr` (296 K) -> `heatr` (4 MT KERMA)
- **Input**: `resources/n-050_Sn_119-ENDF8.0.endf` (tape20)
- **MAT**: 5046, tolerance 0.001
- **Reference**: `referenceTape30` (PENDF with KERMA)
- **Action needed**: HEATR with local photon option (iprint=1).

```julia
result = reconr(endf_file; mat=5046, err=0.001)
# ... broadr at 296 K ...
# ... heatr with mts=[302, 402, 443, 444] ...
# Compare against referenceTape30
```

### 9. Test 31: MODER + RECONR + BROADR + PURR -- Pu-240

- **Pipeline**: `moder` -> `reconr` -> `broadr` (293.6 K) -> `purr` (1 temp, 1 sigma0, 20 bins, 4 ladders)
- **Input**: `resources/n-094_Pu_240-ENDF8.0.endf` (tape20)
- **MAT**: 9440, tolerance 0.001
- **Reference**: `referenceTape30` (PENDF with probability tables)
- **Action needed**: Validates PURR probability table generation.

```julia
result = reconr(endf_file; mat=9440, err=0.001)
# ... broadr at 293.6 K ...
# ... purr with 1 temp, 1 sigma0, 20 bins, 4 ladders ...
# Compare against referenceTape30
```

### 10. Test 73: MODER + RECONR + BROADR + GROUPR -- Pb-208

- **Pipeline**: `moder` -> `reconr` -> `broadr` (293.6 K) -> `groupr` (custom 2-group structure)
- **Input**: `resources/n-082_Pb_208-TENDL2021.endf` (tape20)
- **MAT**: 8237, tolerance 0.001
- **Reference**: `referenceTape30` (GENDF with group-averaged data)
- **Action needed**: Simplest GROUPR test (only 2 groups). Good first GROUPR validation.
- **Group boundaries**: `[1.0e-5, 1.8e8, 2.0e8]` (custom structure)

```julia
result = reconr(endf_file; mat=8237, err=0.001)
# ... broadr at 293.6 K ...
# ... groupr with 2-group custom structure, IWT=2, 10 Legendre orders ...
# Compare against referenceTape30
```

---

## Implementation Notes

### What is needed to run these tests end-to-end

1. **Input deck parser**: NJOY.jl does not currently parse NJOY2016 input decks. Each test must be manually translated into Julia API calls. The input deck format is module-name + parameter cards separated by `/`.

2. **PENDF writer for comparison**: The `write_pendf` / `write_pendf_file` functions exist but must produce output in the exact ENDF-6 format (66-char data fields, 11-char floats, MAT/MF/MT columns) to match reference tapes byte-for-byte. The NJOY2016 comparison uses `execute.py` which compares floats to 1e-9 relative tolerance.

3. **Pipeline composition**: Each test is a chain of module calls where intermediate tapes are passed between stages. In NJOY.jl, this maps to:
   - `reconr()` returns a `PointwiseMaterial`
   - `doppler_broaden()` takes and returns `PointwiseMaterial`
   - `compute_kerma()` takes `PointwiseMaterial`
   - `compute_thermal()` takes `PointwiseMaterial`
   - `group_average()` / `group_average_shielded()` takes pointwise data
   - `build_ace_from_pendf()` takes `PointwiseMaterial`
   - `process_covariance()` takes ENDF file path

4. **Reference tape format**: Reference tapes are in ENDF-6 text format. The existing `read_reference_pendf()` in `test/integration_tests.jl` parses MF3 sections. For GENDF (groupr output), ACE (acer output), and covariance (errorr output) tapes, additional parsers are needed.

### Priority recommendations

**Tier 1 -- Run immediately** (reconr/broadr only, minimal infrastructure):
- Tests 84, 43, 44, 81, 83, 85 (reconr and broadr pipeline, 3 already implemented)

**Tier 2 -- Run with HEATR integration** (add KERMA output comparison):
- Tests 26, 79, 08

**Tier 3 -- Run with GROUPR/ERRORR integration** (add GENDF/covariance parsers):
- Tests 73, 60, 29, 46, 04

**Tier 4 -- Full pipeline with ACER** (add ACE comparison logic):
- Tests 07, 08, 10, 19, 25, 49, 82

**Tier 5 -- THERMR/PURR integration** (specialized processing):
- Tests 31, 32, 01, 17

### Key observation about PARTIAL tests

39 of the 51 PARTIAL tests are missing **only** `viewr` (a PostScript plotting module). The `viewr` step is always the final step and produces only a `.ps` plot file -- it does not affect the nuclear data output tapes. All 39 have reference tapes for the pre-viewr output stages. These tests could be run by simply skipping the `viewr` step and comparing against the available reference tapes. This would effectively raise the total runnable count from **24 to 63 tests** (75% of the suite).

The remaining 12 PARTIAL tests need other unavailable modules:
- `gaspr`: tests 12, 13, 24, 45, 72
- `ccccr`: test 02
- `wimsr`: test 11
- `leapr`: test 09
- `gaminr`/`matxsr`: tests 03, 30
- `acer` photoatomic mode (iopt=4): test 59
- `acer` check/rename mode (iopt=7): test 61

---

## Quick Reference: All 84 Tests by Status

```
READY (24):   01, 04, 07, 08, 10, 17, 19, 25, 26, 29, 31, 32, 43, 44, 46,
              49, 60, 73, 79, 81, 82, 83, 84, 85
PARTIAL (51): 02, 03, 05, 09, 11, 12, 13, 14, 15, 16, 18, 20, 24, 27, 28,
              30, 34, 35, 36, 37, 38, 39, 40, 41, 42, 45, 47, 50, 51, 52,
              53, 54, 55, 56, 57, 58, 59, 61, 62, 63, 64, 65, 66, 67, 68,
              69, 70, 71, 72, 74, 78
NOT_READY (5): 06, 22, 23, 33, 80
EMPTY (4):     21, 48, 75, 76
```
