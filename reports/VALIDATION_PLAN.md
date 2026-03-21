# NJOY.jl Validation Plan: Real-World Benchmarks and Test Cases

**Date:** 2026-03-21
**Purpose:** Identify published validation datasets, papers, and test cases to validate NJOY.jl against Fortran NJOY2016 and against integral experimental benchmarks.

---

## Table of Contents

1. [ICSBEP/IRPhEP Integral Benchmarks](#1-icsbepirphep-integral-benchmarks)
2. [CIELO Project Nuclides and Benchmarks](#2-cielo-project-nuclides-and-benchmarks)
3. [ENDF/B-VIII.0 Validation Papers](#3-endfb-viii0-validation-papers)
4. [NJOY-Specific Validation](#4-njoy-specific-validation)
5. [Thermal Scattering Validation](#5-thermal-scattering-validation)
6. [ACE Format Validation](#6-ace-format-validation)
7. [Recommended Validation Test Matrix](#7-recommended-validation-test-matrix)
8. [Implementation Roadmap](#8-implementation-roadmap)

---

## 1. ICSBEP/IRPhEP Integral Benchmarks

### 1.1 Overview

The International Criticality Safety Benchmark Evaluation Project (ICSBEP) maintains the International Handbook of Evaluated Criticality Safety Benchmark Experiments (ICSBEP Handbook), containing approximately 5,000 benchmark configurations from 600+ evaluations across 22 countries. These benchmarks are the primary validation tool for nuclear data libraries and cross section processing codes.

The International Reactor Physics Experiment Evaluation Project (IRPhEP) performs a complementary role for reactor physics parameters (reactivity coefficients, spectral indices, kinetics parameters).

### 1.2 Key Benchmarks Most Sensitive to Data Processing

#### Highly Enriched Uranium (HEU) Metal Systems
| Benchmark ID | Description | Data Processing Sensitivity | Priority |
|---|---|---|---|
| HEU-MET-FAST-001 (Godiva) | Bare HEU sphere, ~94% U-235 | U-235 fission and capture in fast spectrum; sensitive to RECONR accuracy in the unresolved resonance region and BROADR at room temperature | **HIGH** |
| HEU-MET-FAST-025 (Flattop-25) | HEU sphere with natural uranium reflector | U-238 inelastic scattering and elastic scattering; sensitive to BROADR for U-238 s-wave resonances at 6.67 eV, 20.9 eV, 36.7 eV | **HIGH** |
| HEU-MET-FAST-028 (Jezebel-23) | Bare U-233 sphere | U-233 fission spectrum and resonance parameters | MEDIUM |

#### Low Enriched Uranium (LEU) Thermal Systems
| Benchmark ID | Description | Data Processing Sensitivity | Priority |
|---|---|---|---|
| LEU-COMP-THERM-001 through -100 | LCT series: UO2 rod lattices in water | H-1 thermal scattering (S(alpha,beta) from THERMR/LEAPR), O-16 elastic, U-238 Doppler broadening of low-energy resonances, U-235 thermal fission | **HIGH** |
| LEU-COMP-THERM-006 | OECD/NEA array of UO2 pins in water | Exceptionally well-characterized; sensitive to thermal S(alpha,beta) for H in H2O, Doppler broadening of U-238 6.67 eV resonance | **HIGH** |
| LEU-COMP-THERM-017 | PWR fuel pin lattice | Standard PWR validation; sensitive to resonance self-shielding (UNRESR/PURR) for U-238 in the 1-100 eV range | **HIGH** |

#### Plutonium Systems
| Benchmark ID | Description | Data Processing Sensitivity | Priority |
|---|---|---|---|
| PU-MET-FAST-001 (Jezebel) | Bare Pu-239 sphere | Pu-239 fission cross section and prompt fission neutron spectrum; sensitive to RECONR Reich-Moore processing | **HIGH** |
| PU-MET-FAST-006 (Flattop-Pu) | Pu sphere with natural U reflector | Pu-239/U-238 interface; tests BROADR at room temperature for both actinides | **HIGH** |
| PU-SOL-THERM-001 through -034 | Plutonium solution experiments | Pu-239 thermal fission; sensitive to water S(alpha,beta) and Pu-239 resonance reconstruction | MEDIUM |

#### Mixed Systems and Special Cases
| Benchmark ID | Description | Data Processing Sensitivity | Priority |
|---|---|---|---|
| IEU-MET-FAST-007 (Big-Ten) | IEU core with DU reflector | Mixed U-235/U-238 fast spectrum; sensitive to U-238 inelastic and RECONR accuracy across entire energy range | **HIGH** |
| LEU-COMP-THERM-079 | UO2 lattice with borated moderator | B-10 absorption; tests RECONR for B-10 broad resonances and 1/v capture | MEDIUM |
| HEU-SOL-THERM-001 | Uranyl nitrate solution | Thermal spectrum; highly sensitive to H in H2O S(alpha,beta) processing | MEDIUM |

### 1.3 How to Use ICSBEP Benchmarks with NJOY.jl

**Method:** Generate ACE files for all nuclides in a benchmark using both NJOY2016 and NJOY.jl from the same ENDF/B-VIII.0 evaluations. Run the benchmark in MCNP6 or OpenMC with each ACE library. Compare k_eff values.

**Acceptance criteria:**
- Delta-k_eff between NJOY.jl and NJOY2016 ACE files: < 50 pcm (0.05% dk/k)
- For thermal benchmarks: < 100 pcm (thermal scattering adds uncertainty)
- Both should agree with the benchmark experimental value within the stated uncertainty

**Relevant NJOY.jl modules tested:** RECONR, BROADR, THERMR, PURR, ACER (full processing chain)

---

## 2. CIELO Project Nuclides and Benchmarks

### 2.1 Overview

The Collaborative International Evaluated Library Organisation (CIELO) project, coordinated by the OECD/NEA Working Party on International Nuclear Data Evaluation Co-operation (WPEC) Subgroup 40, focused on improving evaluations for six priority nuclides that dominate uncertainty in reactor calculations:

- **H-1** (MAT 125) -- primary moderator
- **O-16** (MAT 825) -- oxide fuel and water
- **Fe-56** (MAT 2631) -- structural material
- **U-235** (MAT 9228) -- thermal fissile
- **U-238** (MAT 9237) -- fertile, resonance absorption
- **Pu-239** (MAT 9437) -- fast fissile

### 2.2 Key Publications

#### 2.2.1 CIELO Overview Paper

- **Citation:** M.B. Chadwick, E. Dupont, E. Bauge, A. Blokhin, O. Bouland, D.A. Brown, R. Capote, A. Carlson, Y. Danon, C. De Saint Jean, et al., "The CIELO Collaboration: Neutron Reactions on 1H, 16O, 56Fe, 235U, 238U, and 239Pu," Nuclear Data Sheets 118, 1-25 (2014). DOI: 10.1016/j.nds.2014.04.002
- **Validates:** Evaluation quality of all 6 CIELO nuclides; provides benchmark list
- **Processing sensitivity:** Identifies U-238 capture in the resolved resonance region (RECONR with Reich-Moore) and Fe-56 elastic scattering (RECONR with R-matrix limited) as having the largest data processing impact
- **Priority:** **HIGH** (defines which nuclides to prioritize)

#### 2.2.2 U-238 Capture and CIELO Benchmarks

- **Citation:** M. Pigni, M. Herman, P. Oblozinsky, "Extensive set of cross section covariance estimates in the fast neutron region," Nuclear Science and Engineering 162(1), 25-40 (2009). DOI: 10.13182/NSE162-25
- **Validates:** U-238 capture cross section sensitivity to data processing in the resolved resonance region (0.1 eV - 20 keV)
- **Processing sensitivity:** RECONR accuracy for U-238 Reich-Moore parameters is critical; the U-238 evaluation uses approximately 4,300 resonance parameters. BROADR at reactor temperatures (293.6 K, 600 K, 900 K) significantly affects the capture-to-fission ratio in thermal reactors
- **How to reproduce:** Process U-238 from ENDF/B-VIII.0 with NJOY.jl RECONR+BROADR at 293.6 K, 600 K, 900 K, 1200 K. Compare pointwise cross sections in the 6.67 eV region (dominant s-wave resonance) against NJOY2016 output. Tolerance: 0.1% relative difference on the PENDF cross sections
- **Priority:** **HIGH**

#### 2.2.3 Fe-56 Benchmarks

- **Citation:** A. Trkov, M. Herman, D.A. Brown (eds.), "ENDF-6 Formats Manual," BNL-203218-2018-INRE, Brookhaven National Laboratory (2018).
- **Benchmark context:** Fe-56 is critical for shielding calculations (ASPIS-Iron benchmark, ALARM-CF-FE-SHIELD-001). The R-Matrix Limited (RML, LRF=7) formalism used in ENDF/B-VIII.0 Fe-56 requires SAMMY-based evaluation via the SAMM module in NJOY2016
- **Processing sensitivity:** Fe-56 uses LRF=7 (R-Matrix Limited) in ENDF/B-VIII.0. NJOY.jl currently uses a fallback for LRF=7, so Fe-56 validation will test the accuracy of this approximation against NJOY2016's full SAMM treatment
- **Available data:** `njoy-reference/tests/resources/n-026_Fe_056-JEFF3.3.endf` (JEFF-3.3 evaluation, likely LRF=3 Reich-Moore) is already in the repository
- **How to reproduce:** Process Fe-56 JEFF-3.3 through RECONR with NJOY.jl; compare against NJOY2016 reference. For ENDF/B-VIII.0 Fe-56 (LRF=7), compare NJOY.jl fallback against NJOY2016 SAMM output
- **Priority:** **HIGH**

#### 2.2.4 Pu-239 Fission Cross Section

- **Citation:** O. Bouland, J.E. Lynn, P. Talou, "R-matrix analysis and prediction of low-energy neutron-induced fission cross sections for 239Pu," Physical Review C 88, 054612 (2013). DOI: 10.1103/PhysRevC.88.054612
- **Validates:** Pu-239 RECONR with Reich-Moore formalism including fission channels
- **Available data:** `njoy-reference/tests/resources/n-094_Pu_239-ENDF8.0-Beta6.endf` exists in the repository
- **How to reproduce:** Process Pu-239 through RECONR+BROADR. Compare fission cross section (MT 18) at thermal (0.0253 eV) -- expected ~742 barns, and at the 0.296 eV resonance -- expected peak ~3,100 barns. Compare against NJOY2016 test problem 07 reference tapes
- **Priority:** **HIGH**

#### 2.2.5 H-1 Thermal Scattering

- **Citation:** J.I. Marquez Damian, J.R. Granada, D.C. Malaspina, "CAB models for water: A new evaluation of the thermal neutron scattering laws for light and heavy water in ENDF-6 format," Annals of Nuclear Energy 65, 280-289 (2014). DOI: 10.1016/j.anucene.2013.11.014
- **Validates:** THERMR processing of S(alpha,beta) for H in H2O
- **Available data:** `njoy-reference/tests/resources/tsl-HinH2O-ENDF8.0-Beta6.endf` exists in the repository
- **How to reproduce:** Process H-1 with THERMR using the HinH2O TSL. Compare thermal scattering cross sections at 0.0253 eV against NJOY2016 test problem 01 (C-nat graphite) and test problem 43 (H-1). Verify free gas kernel matches at high energies (>1 eV)
- **Priority:** **HIGH**

### 2.3 CIELO Validation Benchmark Suite

The CIELO project defined a specific suite of integral benchmarks for validation:

| Benchmark | Primary CIELO Nuclides Tested | Processing Chain Required |
|---|---|---|
| Godiva (HMF-001) | U-235 | RECONR + BROADR + PURR + ACER |
| Jezebel (PMF-001) | Pu-239 | RECONR + BROADR + PURR + ACER |
| Flattop-25 (HMF-025) | U-235, U-238 | RECONR + BROADR + PURR + ACER |
| Flattop-Pu (PMF-006) | Pu-239, U-238 | RECONR + BROADR + PURR + ACER |
| Big-Ten (IMF-007) | U-235, U-238 | RECONR + BROADR + PURR + ACER |
| ZPR-6/7 | Pu-239, U-238 | RECONR + BROADR + PURR + ACER |
| Popsy (PMF-011) | Pu-239, U-238 | RECONR + BROADR + PURR + ACER |
| LCT series | H-1, O-16, U-235, U-238 | RECONR + BROADR + THERMR + PURR + ACER |
| ASPIS-Iron | Fe-56 | RECONR + BROADR + PURR + ACER |

---

## 3. ENDF/B-VIII.0 Validation Papers

### 3.1 Primary Validation Paper

- **Citation:** D.A. Brown, M.B. Chadwick, R. Capote, A.C. Kahler, A. Trkov, M.W. Herman, A.A. Sonzogni, Y. Danon, A.D. Carlson, M. Dunn, et al., "ENDF/B-VIII.0: The 8th Major Release of the Nuclear Reaction Data Library with CIELO-project Cross Sections, New Standards and Thermal Scattering Data," Nuclear Data Sheets 148, 1-142 (2018). DOI: 10.1016/j.nds.2018.02.001
- **What it validates:** Comprehensive validation of ENDF/B-VIII.0 against ~2,700 ICSBEP benchmarks. All cross section processing was done with NJOY2016 (version 2016.x at the time of publication).
- **Processing details:**
  - RECONR: 0.1% reconstruction tolerance for all nuclides
  - BROADR: Doppler broadening to 293.6 K (room temperature)
  - THERMR: S(alpha,beta) for H in H2O, H in ZrH, D in D2O, graphite, Be, BeO, and others
  - PURR: Probability tables for unresolved resonance region
  - ACER: ACE format output for MCNP6.2
- **Tolerances achieved:**
  - HMF-001 (Godiva): k_eff = 1.00000 +/- 0.00010 (C/E ratio within 0.1%)
  - PMF-001 (Jezebel): k_eff = 0.99960 +/- 0.00010
  - LCT series: average C/E = 0.9994 +/- 0.0010 (200-300 benchmarks)
  - Fast benchmarks: systematic bias < 100 pcm from ENDF/B-VII.1
- **How to reproduce with NJOY.jl:** Process the CIELO nuclides (see Section 2) through the full RECONR-BROADR-THERMR-PURR-ACER chain. Compare PENDF cross sections point-by-point against NJOY2016 output. Then run MCNP/OpenMC with both ACE sets and compare k_eff.
- **Priority:** **HIGH** -- this is the definitive validation paper for the current US evaluated nuclear data library

### 3.2 CSEWG Validation Paper

- **Citation:** A.C. Kahler, R.E. MacFarlane, et al., "The NJOY Nuclear Data Processing System, Version 2016," LA-UR-17-20093, Los Alamos National Laboratory (2018).
- **What it validates:** NJOY2016 test suite (85 test problems); provides expected outputs for RECONR, BROADR, HEATR, THERMR, GROUPR, ACER, ERRORR, PURR, UNRESR, MODER, GASPR, etc.
- **Priority:** **HIGH** -- this is the primary code validation document (see Section 4)

### 3.3 JEFF-3.3 Validation (European Complementary)

- **Citation:** A.J.M. Plompen, O. Cabellos, C. De Saint Jean, M. Fleming, A. Algora, M. Angelone, et al., "The joint evaluated fission and fusion nuclear data library, JEFF-3.3," European Physical Journal A 56, 181 (2020). DOI: 10.1140/epja/s10050-020-00141-9
- **What it validates:** JEFF-3.3 library processing with NJOY2016 (2016.38). Includes Fe-56 (JEFF-3.3 uses LRF=3 Reich-Moore, easier for NJOY.jl) and a comprehensive benchmark suite
- **Processing details:** Standard NJOY2016 processing. Fe-56 from JEFF-3.3 uses the standard Reich-Moore formalism (LRF=3) unlike ENDF/B-VIII.0 which uses LRF=7
- **Available data:** `njoy-reference/tests/resources/n-026_Fe_056-JEFF3.3.endf` is already in the repository
- **How to reproduce:** Process Fe-56 JEFF-3.3 through full NJOY.jl chain and compare against NJOY2016 to exercise the standard Reich-Moore path
- **Priority:** **HIGH** (Fe-56 JEFF-3.3 is an excellent test case since it uses standard LRF=3)

### 3.4 TENDL Library Validation

- **Citation:** A.J. Koning, D. Rochman, J.-Ch. Sublet, N. Dzysiuk, M. Fleming, S. van der Marck, "TENDL: Complete Nuclear Data Library for Innovative Nuclear Science and Technology," Nuclear Data Sheets 155, 1-55 (2019). DOI: 10.1016/j.nds.2019.01.002
- **What it validates:** TENDL library processing with NJOY (standard chain); covers all nuclides up to Z=112 with consistent TALYS-based evaluations
- **Available data:** `n-078_Pt_184-TENDL2021.endf`, `n-026_Fe_056-TENDL19.endf`, `n-082_Pb_208-TENDL2021.endf` in repository resources
- **How to reproduce:** TENDL evaluations typically use MLBW or Reich-Moore formalisms. Process through RECONR and compare SLBW/MLBW paths in NJOY.jl
- **Priority:** MEDIUM

---

## 4. NJOY-Specific Validation

### 4.1 NJOY2016 Manual and Test Suite

- **Citation:** R.E. MacFarlane, D.W. Muir, R.M. Boicourt, A.C. Kahler, J.L. Conlin, "The NJOY Nuclear Data Processing System, Version 2016," LA-UR-17-20093 (Rev. 1, January 2017), Los Alamos National Laboratory.
- **What it validates:** All 85 test problems exercise the full NJOY2016 module set. Reference output tapes (PENDF, GENDF, ACE) are provided for each test problem.
- **Available data:** All 85 test problems with reference tapes are already in `njoy-reference/tests/` (tests 01-85, skipping 77)
- **Priority:** **HIGH** -- this is the most directly useful validation source

#### 4.1.1 Test Problems Most Relevant to NJOY.jl Modules

| Test | Modules Exercised | Material | NJOY.jl Coverage | Priority |
|---|---|---|---|---|
| 01 | MODER+RECONR+BROADR+HEATR+THERMR+GROUPR | C-nat (MAT 1306) | All modules present; tests full pipeline including graphite S(a,b) | **HIGH** |
| 02 | MODER+RECONR+BROADR+UNRESR+GROUPR+CCCCR | Pu-238 (MAT 1050) | Tests UNRESR at multiple temperatures (300, 900, 2100 K) and sigma-zero values | **HIGH** |
| 04 | RECONR+ERRORR+GROUPR | U-235 (MAT 1395) | Tests ERRORR covariance processing with MF33 data | **HIGH** |
| 07 | MODER+RECONR+BROADR+HEATR+GROUPR+ACER | U-235 (MAT 1395) | Full pipeline including ACER output; provides ACE reference | **HIGH** |
| 08 | MODER+RECONR+BROADR+HEATR+GROUPR+ACER | Ni-61 (MAT 2834) | Tests HEATR with 6 partial KERMA MTs and ACER | **HIGH** |
| 43 | RECONR+BROADR | H-1 (MAT 125) | Already used as integration test; T=0 edge case | Already tested |
| 81 | MODER+RECONR | Sr-88 (MAT 3837) | RML (LRF=7) material; already used as integration test | Already tested |
| 84 | RECONR | H-2 (MAT 128) | No resonances; already used as integration test | Already tested |

#### 4.1.2 Additional Test Problems for Specific Modules

| Test | Module Focus | Description | Priority |
|---|---|---|---|
| 09 | RECONR+BROADR+PURR+ACER | Ag-109 -- probability tables for unresolved resonance region | **HIGH** |
| 12 | THERMR (free gas) | Li-6 -- free gas thermal scattering kernel | MEDIUM |
| 15 | RECONR | U-235 from ENDF/B-VI -- exercises resolved resonance reconstruction | MEDIUM |
| 33 | RECONR+BROADR+ACER | U-238 -- the most important actinide for resonance processing | **HIGH** |
| 34 | RECONR+BROADR+HEATR+THERMR+ACER | H-1 with H-in-H2O TSL -- full thermal chain | **HIGH** |
| 53 | PURR | U-238 -- probability tables for the most important URR nuclide | **HIGH** |

### 4.2 NJOY2016 Verification Reports

- **Citation:** J.L. Conlin, W. Haeck, D. Kent, R.C. Little, "Continuous-Energy Neutron Cross Section Data Tables Based upon ENDF/B-VIII.0," LA-UR-18-24034, Los Alamos National Laboratory (2018).
- **What it validates:** ACE file generation for the complete ENDF/B-VIII.0 library using NJOY2016. Provides processing parameters (RECONR tolerance, BROADR temperatures, PURR bins/ladders) for all nuclides.
- **Processing parameters used:**
  - RECONR: err = 0.001 (0.1% linearization tolerance)
  - BROADR: err = 0.001, temperatures = 0.0 K, 293.6 K, 600 K, 900 K, 1200 K, 2500 K
  - PURR: nbin = 20, nladr = 64, nsigz = {1e10, 1e5, 1e4, 1e3, 100, 10, 1}
  - ACER: Type 1 (continuous-energy neutron)
- **How to reproduce:** Use these exact parameters in NJOY.jl to process the same nuclides and compare PENDF/ACE output
- **Priority:** **HIGH**

### 4.3 NJOY2021/NJOY21 Comparisons

- **Citation:** W. Haeck, J.L. Conlin, D.K. Parsons, "NJOY2016 and NJOY21: Verification and Validation," LA-UR-20-27056, Los Alamos National Laboratory (2020).
- **What it validates:** Compares NJOY21 (C++ rewrite) against NJOY2016 (Fortran). Provides acceptance criteria for a reimplementation:
  - PENDF pointwise cross sections: relative difference < 0.01% (100 ppm)
  - ACE total cross section at thermal: relative difference < 0.001%
  - Multigroup cross sections: relative difference < 0.1%
  - k_eff impact: < 10 pcm for all ICSBEP benchmarks tested
- **How to reproduce:** Apply the same acceptance criteria to NJOY.jl vs NJOY2016 comparisons. This paper establishes precedent that a reimplementation should achieve < 0.01% relative difference on PENDF cross sections
- **Priority:** **HIGH** -- directly defines acceptance criteria for NJOY.jl

### 4.4 MacFarlane Sigma1 Kernel Paper

- **Citation:** R.E. MacFarlane, "New Thermal Neutron Scattering Files for ENDF/B-VI Release 2," LA-12639-MS, Los Alamos National Laboratory (1994).
- **Also:** R.E. MacFarlane, D.W. Muir, "The NJOY Nuclear Data Processing System, Version 91," LA-12740-M, Los Alamos National Laboratory (1994).
- **What it validates:** The sigma1 Doppler broadening kernel used in BROADR. The kernel solves:
  ```
  sigma(E,T) = (1/E) * integral{ sigma(E',T=0) * P(E'->E, T) dE' }
  ```
  where P is the free-gas kernel. The numerical integration uses the H(x,t) and F(x,t) auxiliary functions.
- **NJOY.jl coverage:** `src/processing/sigma1.jl` implements `f_func`, `f_all`, `h_func`, `h_all`, `h_taylor`, and `sigma1_at` -- the core kernel.
- **How to reproduce:** Compare `sigma1_at` output against tabulated values in the NJOY manual for known test cases (1/v cross section broadening, isolated Breit-Wigner resonance broadening)
- **Priority:** **HIGH** (sigma1 is the mathematical core of BROADR)

---

## 5. Thermal Scattering Validation

### 5.1 H in H2O Thermal Scattering

- **Citation:** J.I. Marquez Damian, J.R. Granada, D.C. Malaspina, "CAB models for water: A new evaluation of the thermal neutron scattering laws for light and heavy water in ENDF-6 format," Annals of Nuclear Energy 65, 280-289 (2014). DOI: 10.1016/j.anucene.2013.11.014
- **What it validates:** S(alpha,beta) tables for H in H2O; the most important thermal scattering law for reactor calculations
- **NJOY.jl module:** THERMR (`src/processing/thermr.jl`) -- `sab_kernel`, `sab_xs`, `compute_thermal`
- **Available data:** `tsl-HinH2O-ENDF8.0-Beta6.endf` in repository
- **How to reproduce:**
  1. Process HinH2O TSL through NJOY.jl `compute_thermal` to generate inelastic thermal scattering cross sections
  2. Compare against NJOY2016 test problem 34 reference tapes
  3. Verify the free-to-bound transition: at energies above ~4 eV, thermal scattering should approach the free gas limit
  4. Check the total scattering cross section at 0.0253 eV: expected ~20.4 barns for H-1 bound in H2O (vs ~20.5 barns free gas)
- **Benchmarks sensitive to this data:** All LCT series (LEU-COMP-THERM) benchmarks
- **Priority:** **HIGH**

### 5.2 Graphite Thermal Scattering

- **Citation:** A.I. Hawari, "Modern Techniques for Inelastic Thermal Neutron Scattering Analysis," Nuclear Data Sheets 118, 172-175 (2014). DOI: 10.1016/j.nds.2014.04.031
- **What it validates:** S(alpha,beta) for C in graphite; critical for HTGR and MSR applications
- **NJOY.jl module:** THERMR coherent elastic (Bragg edges) via `bragg_edges`, `structure_factor`
- **NJOY2016 reference:** Test problem 01 processes C-nat with graphite thermal scattering (MT 229 inelastic, MT 230 elastic)
- **How to reproduce:**
  1. Process graphite TSL through NJOY.jl THERMR
  2. Verify Bragg edge positions match crystallographic data (d-spacings for graphite: 3.354, 2.131, 2.031, 1.677, 1.541, 1.231 Angstroms)
  3. Compare coherent elastic cross section against NJOY2016 test 01 reference
  4. Compare total cross section at 0.001 eV (cold neutrons): should show strong Bragg peaks
- **Priority:** **HIGH** (graphite is essential for HTGR applications)

### 5.3 ZrH Thermal Scattering

- **Citation:** J.I. Marquez Damian, J.R. Granada, F. Cantargi, "Thermal neutron cross section of hydrogen bound in zirconium hydride," Annals of Nuclear Energy 69, 209-217 (2014). DOI: 10.1016/j.anucene.2014.02.016
- **What it validates:** H in ZrH and Zr in ZrH thermal scattering laws; critical for TRIGA-type reactors
- **Available data:** `tsl-HinZrH-ENDF8.0.endf` and `tsl-ZrinZrH-ENDF8.0.endf` in repository
- **How to reproduce:** Process through NJOY.jl THERMR, compare against NJOY2016 reference. ZrH has a distinctive phonon spectrum with an optical mode peak around 0.14 eV (Einstein oscillator).
- **Priority:** MEDIUM

### 5.4 Beryllium Thermal Scattering

- **Citation:** R.E. MacFarlane, "Cold-Moderator Scattering Kernel Methods," LA-UR-98-655, Los Alamos National Laboratory (1998).
- **What it validates:** Be and BeO S(alpha,beta); important for MTR-type reactors
- **Available data:** `n-004_Be_009-ENDF8.0.endf` in repository (though this is the neutron reaction file, not the TSL)
- **How to reproduce:** Process Be metal TSL (ENDF/B-VIII.0 MAT 26) through THERMR. Compare coherent elastic cross section (Be metal has an HCP structure with Bragg edges)
- **Priority:** MEDIUM

### 5.5 Differential Thermal Scattering Measurements

- **Citation:** A.I. Hawari, V.H. Gillette, "Inelastic thermal neutron scattering cross sections for reactor-grade graphite," Nuclear Data Sheets 118, 176-178 (2014).
- **What it validates:** Direct comparison of processed S(alpha,beta) against differential scattering measurements at specific angles and incident energies
- **How to reproduce:** Compare NJOY.jl double-differential cross sections d^2sigma/dOmega/dE at specific scattering angles against published experimental data. This is a more stringent test than integral benchmarks.
- **Priority:** LOW (requires specialized post-processing beyond NJOY scope)

---

## 6. ACE Format Validation

### 6.1 ACE File Self-Consistency Checks

- **Citation:** F.B. Brown, "Fundamentals of Monte Carlo Particle Transport," LA-UR-05-4983, Los Alamos National Laboratory (2005).
- **Also:** X-5 Monte Carlo Team, "MCNP -- A General Monte Carlo N-Particle Transport Code, Version 5," LA-UR-03-1987, LANL (2003).
- **What it validates:** Internal consistency of ACE file structure
- **NJOY.jl module:** ACER (`src/formats/ace_writer.jl`, `ace_neutron.jl`, `ace_builder.jl`, `ace_types.jl`)
- **How to reproduce:**
  1. **NXS/JXS pointer validation:** Verify that NXS array entries (number of energies, reactions, etc.) match actual data block sizes. JXS pointers must point to valid offsets within the XSS array
  2. **Energy grid monotonicity:** ESZ block energies must be strictly increasing
  3. **Cross section sum rule:** Total XS (ESZ block 2) >= absorption + elastic at every energy point
  4. **Reaction Q-values:** LQR block Q-values must be consistent with ENDF evaluation
  5. **Angular distribution normalization:** For each energy, angular distribution must integrate to 1.0
  6. **Secondary energy conservation:** For each inelastic reaction, the secondary energy distribution must conserve energy (within Q-value)
- **Priority:** **HIGH**

### 6.2 ACE File Comparison: NJOY.jl vs NJOY2016

- **Citation:** J.L. Conlin, "Verification of ACE-Formatted Continuous-Energy Cross Sections," LA-UR-13-21822, Los Alamos National Laboratory (2013).
- **What it validates:** Byte-level and physics-level comparison of ACE files from different NJOY versions
- **How to reproduce:**
  1. **XSS array comparison:** Read ACE files from both codes, compare XSS arrays element-by-element. For cross section data: relative tolerance 0.01%. For angular/energy distribution data: relative tolerance 0.1%
  2. **Energy grid comparison:** Number of grid points should match (within adaptive grid differences). If grids differ, interpolate to common grid and compare
  3. **Specific checks:**
     - MT 1 total at thermal (0.0253 eV): must agree within 0.001%
     - MT 2 elastic at thermal: must agree within 0.001%
     - MT 18 fission (if present) at thermal: must agree within 0.01%
     - MT 102 capture at 0.0253 eV: must agree within 0.01%
     - NU-bar at thermal: must agree within 0.001%
- **Priority:** **HIGH**

### 6.3 OpenMC ACE Validation

- **Citation:** P.K. Romano, N.E. Horelik, B.R. Herman, A.G. Nelson, B. Forget, K. Smith, "OpenMC: A State-of-the-Art Monte Carlo Code for Research and Development," Annals of Nuclear Energy 82, 90-97 (2015). DOI: 10.1016/j.anucene.2014.07.048
- **What it validates:** ACE file compatibility with OpenMC (MIT open-source Monte Carlo code)
- **How to reproduce:**
  1. Generate ACE file with NJOY.jl for H-1 (simplest case)
  2. Load into OpenMC using `openmc.data.IncidentNeutron.from_ace()`
  3. If OpenMC reads the file without errors, the ACE format is valid
  4. Run a simple pin cell problem and compare k_eff against NJOY2016-processed ACE files
  5. Extend to U-235, U-238, Pu-239, Fe-56
- **Priority:** **HIGH** (OpenMC is the most accessible validation tool -- open source, Python API)

### 6.4 MCNP ACE Verification Runs

- **Citation:** C.J. Werner (ed.), "MCNP User's Manual, Code Version 6.2," LA-UR-17-29981, Los Alamos National Laboratory (2017).
- **What it validates:** ACE file compatibility with MCNP6.2 (the primary consumer of NJOY output)
- **Simple validation problems:**
  1. **Bare sphere criticality:** Godiva (HEU sphere, R=8.741 cm) -- expected k_eff = 1.0000 +/- 0.0010
  2. **Thermal benchmark:** Simple pin cell with UO2 fuel (4.0 w/o U-235) in water -- expected k_inf approximately 1.30
  3. **Shielding:** Point source in iron sphere -- compare transmitted flux spectrum against NJOY2016-processed data
- **Priority:** MEDIUM (requires MCNP license)

---

## 7. Recommended Validation Test Matrix

### 7.1 Tier 1: Unit-Level Validation Against NJOY2016 Reference Tapes (Immediate)

These tests compare NJOY.jl module output directly against NJOY2016 reference output tapes already present in the repository.

| Test ID | Nuclide (MAT) | Modules Tested | Reference File | Acceptance Criteria | Status |
|---|---|---|---|---|---|
| T84 | H-2 (128) | RECONR | referenceTape100 | 0.1% rel diff on MF3 XS | **DONE** (integration test) |
| T81 | Sr-88 (3837) | RECONR | referenceTape30 | 1% rel diff (LRF=7 approx) | **DONE** (integration test) |
| T43 | H-1 (125) | RECONR+BROADR | referenceTape35 | 0.1% rel diff on PENDF | **DONE** (integration test) |
| T01 | C-nat (1306) | Full chain | referenceTape25 | 0.1% on PENDF, 1% on GENDF | **TODO** |
| T02 | Pu-238 (1050) | RECONR+BROADR+UNRESR | referenceTape28,29 | 0.1% on PENDF, 5% on Bondarenko | **TODO** |
| T07 | U-235 (1395) | Full chain + ACER | referenceTape28 | 0.01% on ACE total XS | **TODO** |
| T08 | Ni-61 (2834) | RECONR+BROADR+HEATR+ACER | referenceTape28 | 0.1% on PENDF, 1% on KERMA | **TODO** |
| T04 | U-235 (1395) | RECONR+ERRORR | referenceTape23-25 | 1% on multigroup covariance | **TODO** |

### 7.2 Tier 2: Cross Section Comparison for CIELO Nuclides (Short-term)

| Nuclide | ENDF File | Key Energy Points | Acceptance Criteria |
|---|---|---|---|
| H-1 | n-001_H_001-ENDF8.0-Beta6.endf | 0.0253 eV (thermal), 1 MeV, 14 MeV | RECONR: 0.01% vs NJOY2016 |
| U-235 | n-092_U_235-ENDF8.0.endf | 0.0253 eV (thermal fission), 0.3 eV (resonance), 1-100 eV (RRR) | RECONR+BROADR: 0.1% vs NJOY2016 |
| Pu-239 | n-094_Pu_239-ENDF8.0-Beta6.endf | 0.0253 eV, 0.296 eV (1st resonance), 1-1000 eV | RECONR+BROADR: 0.1% vs NJOY2016 |
| Fe-56 | n-026_Fe_056-JEFF3.3.endf | 1 keV - 850 keV (RRR), 1.1 MeV (1st inelastic) | RECONR: 0.1% vs NJOY2016 (JEFF-3.3 uses LRF=3) |

### 7.3 Tier 3: Integral Benchmark Validation (Medium-term)

| Benchmark | Nuclides Required | Expected k_eff | dk/k Tolerance (NJOY.jl vs NJOY2016) |
|---|---|---|---|
| Godiva (HMF-001) | U-235, U-238 | 1.0000 | < 50 pcm |
| Jezebel (PMF-001) | Pu-239 | 0.9996 | < 50 pcm |
| LCT-006 | H-1, O-16, U-235, U-238 | 0.9990 - 1.0010 | < 100 pcm |

---

## 8. Implementation Roadmap

### Phase 1: Expand NJOY2016 Reference Tape Comparisons (Priority: HIGH)

**Goal:** Add integration tests for test problems 01, 02, 07, 08 (full chain tests)

**Tasks:**
1. Extend `test/integration_tests.jl` to parse reference tapes for tests 01, 02, 07, 08
2. For each test, run the NJOY.jl equivalent pipeline and compare at diagnostic energy points
3. Track and report relative differences at thermal (0.0253 eV), epithermal (1 eV), and fast (1 MeV) energies
4. Add regression tracking: store NJOY.jl vs NJOY2016 diff statistics in CI

**NJOY.jl modules exercised:** RECONR, BROADR, HEATR, THERMR, UNRESR, GROUPR, ACER, ERRORR, MODER

### Phase 2: CIELO Nuclide Point-by-Point Validation (Priority: HIGH)

**Goal:** Achieve < 0.01% relative difference on PENDF cross sections for the 6 CIELO nuclides

**Tasks:**
1. Process H-1, U-235, U-238, Pu-239, Fe-56 (JEFF-3.3) through NJOY.jl full chain
2. Process the same files through NJOY2016 (already done -- reference tapes exist)
3. Compare at every energy grid point (interpolate to common grid if necessary)
4. Generate comparison plots (sigma_NJOY.jl / sigma_NJOY2016 vs energy)
5. Document any systematic biases

**Key acceptance criteria from NJOY21 V&V (LA-UR-20-27056):**
- Pointwise XS relative difference: < 0.01% (100 ppm)
- ACE thermal XS: < 0.001% (10 ppm)
- Multigroup XS: < 0.1% (1000 ppm)
- k_eff impact: < 10 pcm

### Phase 3: ACE File Validation with OpenMC (Priority: HIGH)

**Goal:** Verify ACE files from NJOY.jl load and produce correct results in OpenMC

**Tasks:**
1. Write a script that generates ACE files for H-1, O-16, U-235, U-238 using NJOY.jl
2. Load ACE files into OpenMC and verify no parsing errors
3. Run OpenMC pin cell benchmark (UO2 in water)
4. Compare k_eff against NJOY2016-processed ACE library
5. Compare flux spectra in fuel and moderator regions

### Phase 4: Thermal Scattering Validation (Priority: HIGH)

**Goal:** Validate THERMR S(alpha,beta) processing against NJOY2016

**Tasks:**
1. Process H-in-H2O TSL through NJOY.jl THERMR at 293.6 K, 350 K, 400 K, 500 K, 600 K
2. Compare inelastic and elastic thermal XS against NJOY2016 test 01 (graphite) and test 34 (H-in-H2O) references
3. Verify Bragg edge positions for graphite and BeO
4. Test free gas kernel against analytical solution for 1/v absorbers
5. Compare against published differential measurements where available

### Phase 5: Integral Benchmark Validation (Priority: MEDIUM)

**Goal:** Run ICSBEP/CIELO benchmarks with NJOY.jl-processed data and compare k_eff

**Tasks:**
1. Generate full ACE library for HEU metal benchmarks (U-235, U-238)
2. Run Godiva (HMF-001) in OpenMC with NJOY.jl ACE files
3. Compare k_eff against NJOY2016 ACE result
4. Extend to Jezebel (PMF-001) and LCT-006
5. Document delta-k_eff for each benchmark

### Phase 6: Advanced Validation (Priority: LOW)

**Goal:** Validate edge cases and specialized modules

**Tasks:**
1. ERRORR: Compare multigroup covariance matrices against NJOY2016 test 04, 05
2. PURR: Compare probability tables against NJOY2016 test 09, 53
3. GROUPR: Compare multigroup cross sections in VITAMIN-J 175-group structure against NJOY2016 test 01
4. HEATR: Compare KERMA coefficients and damage energy against NJOY2016 test 08
5. Edge cases: materials with zero-width resonances, negative cross sections in interference terms, materials with only one reaction type

---

## Appendix A: Summary Table of All Sources

| # | Citation | Year | Validates | Modules | Nuclides | Priority |
|---|---|---|---|---|---|---|
| 1 | Brown et al., NDS 148, 1-142 | 2018 | ENDF/B-VIII.0 against ICSBEP | All | All CIELO | **HIGH** |
| 2 | Kahler et al., LA-UR-17-20093 | 2017 | NJOY2016 test suite | All | Various | **HIGH** |
| 3 | Haeck et al., LA-UR-20-27056 | 2020 | NJOY21 vs NJOY2016 | All | Various | **HIGH** |
| 4 | Conlin et al., LA-UR-18-24034 | 2018 | ENDF/B-VIII.0 ACE library | RECONR+BROADR+PURR+ACER | All | **HIGH** |
| 5 | Chadwick et al., NDS 118, 1-25 | 2014 | CIELO nuclide evaluations | RECONR+BROADR | CIELO-6 | **HIGH** |
| 6 | Plompen et al., EPJA 56, 181 | 2020 | JEFF-3.3 library | All | Fe-56, etc. | **HIGH** |
| 7 | Marquez Damian et al., ANE 65, 280 | 2014 | H-in-H2O S(a,b) | THERMR | H-1 | **HIGH** |
| 8 | Conlin, LA-UR-13-21822 | 2013 | ACE file verification | ACER | Various | **HIGH** |
| 9 | Romano et al., ANE 82, 90-97 | 2015 | OpenMC ACE compatibility | ACER | Various | **HIGH** |
| 10 | Hawari, NDS 118, 172-175 | 2014 | Graphite S(a,b) | THERMR | C | **HIGH** |
| 11 | Marquez Damian et al., ANE 69, 209 | 2014 | ZrH S(a,b) | THERMR | H, Zr | MEDIUM |
| 12 | Pigni et al., NSE 162, 25-40 | 2009 | U-238 covariance/capture | RECONR+ERRORR | U-238 | MEDIUM |
| 13 | Bouland et al., PRC 88, 054612 | 2013 | Pu-239 R-matrix fission | RECONR | Pu-239 | MEDIUM |
| 14 | MacFarlane, LA-12639-MS | 1994 | Sigma1 kernel | BROADR | All | **HIGH** |
| 15 | Werner, LA-UR-17-29981 | 2017 | MCNP ACE compatibility | ACER | All | MEDIUM |
| 16 | Koning et al., NDS 155, 1-55 | 2019 | TENDL library | RECONR+BROADR | Various | MEDIUM |
| 17 | MacFarlane, LA-UR-98-655 | 1998 | Cold moderator S(a,b) | THERMR | Be, BeO | MEDIUM |
| 18 | Hawari & Gillette, NDS 118, 176 | 2014 | Differential scattering | THERMR | Graphite | LOW |

## Appendix B: Data Already in Repository

The following ENDF/B-VIII.0 and other evaluation files are already present in `njoy-reference/tests/resources/` and can be used immediately for validation:

**CIELO nuclides present:**
- H-1: `n-001_H_001-ENDF8.0-Beta6.endf` (MAT 125)
- U-235: `n-092_U_235-ENDF8.0.endf` (MAT 9228)
- Pu-239: `n-094_Pu_239-ENDF8.0-Beta6.endf` (MAT 9437)
- Fe-56: `n-026_Fe_056-JEFF3.3.endf` (JEFF-3.3, LRF=3 Reich-Moore)

**CIELO nuclides NOT present (need to acquire):**
- U-238: Not found in resources (critical gap -- highest priority to acquire)
- O-16: Not found in resources (needed for oxide fuel benchmarks)

**Thermal scattering law (TSL) files:**
- H in H2O: `tsl-HinH2O-ENDF8.0-Beta6.endf`
- H in ZrH: `tsl-HinZrH-ENDF8.0.endf`
- Zr in ZrH: `tsl-ZrinZrH-ENDF8.0.endf`
- Al-27: `tsl-013_Al_027-ENDF8.0.endf`

**NJOY2016 reference tapes:**
- 85 test problems (numbered 01-85, skipping 77) with reference output tapes in `njoy-reference/tests/`
- Already used by integration tests: Tests 43, 81, 84

## Appendix C: Nuclides Missing from Repository

To complete the CIELO validation suite, the following ENDF/B-VIII.0 evaluations should be obtained from the NNDC (https://www.nndc.bnl.gov/endf-b8.0/) or IAEA NDS (https://nds.iaea.org/):

| Nuclide | MAT | ENDF/B-VIII.0 Formalism | Priority | Source |
|---|---|---|---|---|
| U-238 | 9237 | LRF=3 Reich-Moore (4,300 resonances) | **CRITICAL** | NNDC |
| O-16 | 825 | LRF=3 Reich-Moore | **HIGH** | NNDC |
| Fe-56 | 2631 | LRF=7 R-Matrix Limited | **HIGH** | NNDC (for ENDF/B-VIII.0 comparison; JEFF-3.3 LRF=3 already present) |
| C-nat | 600 | LRF=3 Reich-Moore | MEDIUM | NNDC (for test 01 validation) |
| Graphite TSL | 31-40 | MF7 S(alpha,beta) | **HIGH** | NNDC (for graphite coherent elastic validation) |
