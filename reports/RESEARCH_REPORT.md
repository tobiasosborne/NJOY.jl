# NJOY2016 Deep Research Report

## Phase 0: Comprehensive Source Code Analysis

**Date:** 2026-03-21
**Source:** NJOY2016 v2016.78 (03Feb25), located at `njoy-reference/src/`
**Tests:** 85 test problems (numbered 01-85, skipping 77), at `njoy-reference/tests/`

---

## 1. Module Dependency Graph

### 1.1 All 39 Source Files

NJOY2016 consists of 39 Fortran 90 files totalling 119 613 lines (100 323 code + 13 900 comment + 5 390 blank, per cloc):

| File | Lines | Module Name | Purpose |
|------|-------|-------------|---------|
| `main.f90` | 278 | (program) | Main driver; reads module names from stdin, dispatches to processing modules |
| `vers.f90` | 8 | `version` | Version string `2016.78` and date `03Feb25` |
| `locale.f90` | 20 | `locale` | Kind parameters: `kr=selected_real_kind(12,300)` (double precision), `k4=4`, `k8=8`; site strings |
| `mainio.f90` | 15 | `mainio` | I/O unit numbers: `nsysi` (stdin), `nsyso` (output file), `nsyse` (terminal) |
| `phys.f90` | 51 | `physics` | CODATA 2014 constants: pi, Boltzmann, neutron/proton/alpha masses, hbar, eV, c, fine-structure |
| `util.f90` | 468 | `util` | Utility: `error`, `mess`, `timer`, `dater`, `wclock`, `repoz`, `skiprz`, `openz`, `closz`, `loada`, `finda`, `scana`, `sigfig`, `a10` |
| `endf.f90` | 2148 | `endf` | ENDF I/O: record types (`contio`, `listio`, `tab1io`, `tab2io`, `moreio`, `tpidio`, `hdatio`, `dictio`, `intgio`); delimiter writers (`asend`..`atend`); skippers (`tosend`..`totend`); navigation (`findf`, `skip6`); interpolation (`terp1`, `terpa`, `intega`, `gral`); data access (`gety1`, `gety2`) |
| `mathm.f90` | 1273 | `mathm` | Math routines: `legndr` (Legendre polynomials), `e1` (exponential integral), `gami` (incomplete gamma), `erfc` (complementary error function) |
| `reconr.f90` | 5747 | `reconm` | Pointwise cross section reconstruction from resonance parameters |
| `broadr.f90` | 1947 | `broadm` | Doppler broadening and thinning of pointwise cross sections |
| `unresr.f90` | 1665 | `unresm` | Unresolved resonance self-shielding (Bondarenko method) |
| `heatr.f90` | 6322 | `heatm` | KERMA (kinetic energy release) and damage energy production |
| `thermr.f90` | 3387 | `thermm` | Thermal scattering cross sections and kernels (free gas + S(alpha,beta)) |
| `groupr.f90` | 12590 | `groupm` | Self-shielded multigroup cross sections, scattering matrices, photon production |
| `gaminr.f90` | 1517 | `gaminm` | Multigroup photon interaction cross sections |
| `errorr.f90` | 11152 | `errorm` | Multigroup covariance matrices |
| `covr.f90` | 2250 | `covm` | Covariance data post-processing from ERRORR |
| `moder.f90` | 1714 | `modem` | ENDF format conversion (ASCII <-> blocked binary) |
| `dtfr.f90` | 1510 | `dtfm` | Discrete ordinates transport format output |
| `ccccr.f90` | 3526 | `ccccm` | CCCC standard interface files (ISOTXS, BRKOXS, DLAYXS) |
| `matxsr.f90` | 2512 | `matxsm` | MATXS cross section interface format |
| `resxsr.f90` | 505 | `resxsm` | Pointwise resonance cross sections for thermal flux calculations |
| `acer.f90` | 546 | `acem` | ACE format for MCNP -- dispatcher module |
| `acefc.f90` | 19669 | `acefcm` | ACE fast/continuous neutron data (largest file) |
| `acedo.f90` | 600 | `acedom` | ACE dosimetry data |
| `aceth.f90` | 2433 | `acethm` | ACE thermal scattering data |
| `acepa.f90` | 1034 | `acepam` | ACE photo-atomic data |
| `acepn.f90` | 3786 | `acepnm` | ACE photo-nuclear data |
| `acecm.f90` | 797 | `acecmm` | ACE charged-particle data |
| `powr.f90` | 4902 | `powm` | EPRI-CELL/EPRI-CPM library format |
| `wimsr.f90` | 2150 | `wimsm` | WIMS-D/WIMS-E library format |
| `plotr.f90` | 3149 | `plotm` | Cross section plotting |
| `viewr.f90` | 1669 | `viewm` | PostScript plot rendering |
| `mixr.f90` | 533 | `mixm` | Cross section mixing (elemental from isotopic) |
| `purr.f90` | 2894 | `purm` | Unresolved resonance probability tables for MCNP |
| `leapr.f90` | 3625 | `leapm` | S(alpha,beta) generation for thermal moderators |
| `gaspr.f90` | 1153 | `gaspm` | Gas production (MT203-207) addition to PENDF |
| `samm.f90` | 7169 | `samm` | SAMMY R-matrix method for resonance cross sections and angular distributions |
| `graph.f90` | 2899 | `graph` | Graphics support routines for plotr/viewr |

### 1.2 Dependency Hierarchy

```
Level 0 (Foundation):
  locale -> [no dependencies]
  version -> [no dependencies]
  mainio -> iso_fortran_env

Level 1 (Infrastructure):
  physics -> locale
  util -> locale, mainio
  mathm -> locale, util

Level 2 (ENDF I/O):
  endf -> locale, util

Level 3 (Algorithms):
  samm -> locale, mainio, endf, util, physics
  graph -> locale, mainio

Level 4 (Processing Modules):
  reconm -> locale, util, endf, samm, mainio, physics
  broadm -> locale, util, endf, mathm, physics, mainio
  unresm -> locale, util, endf, physics, mainio
  heatm  -> locale, util, endf, mathm, physics, mainio
  thermm -> locale, util, endf, mathm, physics, mainio
  groupm -> locale, util, endf, mathm, physics, mainio
  errorm -> locale, util, endf, mathm, physics, samm, mainio
  purm   -> locale, util, endf, physics, mainio
  acem   -> locale, util, endf, physics, mainio + (acefcm, acedom, acethm, acepam, acepnm, acecmm)
  leapm  -> locale, util, endf, physics, mainio, mathm
  gaminm -> locale, util, endf, physics, mainio
  gaspm  -> locale, util, endf, physics, mainio
  modem  -> locale, util, endf, mainio
  ...etc
```

### 1.3 I/O Unit Convention

Modules communicate exclusively through Fortran I/O units (file tapes):
- **Positive unit** = ASCII (coded/formatted) ENDF tape
- **Negative unit** = blocked binary NJOY internal tape
- Units 10-19 are scratch (temporary) files
- Units 20+ are user-specified input/output tapes
- File names follow pattern `tape{N}` (e.g., `tape20`, `tape21`)
- `nsysi` (standard input), `nsyso` (output listing file), `nsyse` (terminal)

### 1.4 Module-Level Globals Shared via ENDF Module

The `endf` module exposes critical global state used by all processing modules:
- `c1h, c2h` -- floating point fields from last CONT record read
- `l1h, l2h, n1h, n2h` -- integer fields from last CONT record
- `math, mfh, mth` -- current MAT, MF, MT numbers
- `nsh, nsp, nsc` -- sequence numbers for output/scratch
- `npage` -- page size for blocked binary (default 306, adjusted to multiple of 102)
- `iverf` -- ENDF version flag (4, 5, or 6)
- `thr6` -- threshold for Coulomb penetrability interpolation

---

## 2. Subroutine Catalogue

### 2.1 RECONR (reconm) -- 5747 lines, ~35 subroutines

The central reconstruction module. Converts ENDF resonance parameters and tabulated cross sections into pointwise (PENDF) format.

| Subroutine | Physics Purpose | Key I/O | Algorithm |
|------------|----------------|---------|-----------|
| `reconr` | Main driver | ENDF in, PENDF out | Orchestrates File 1/2/3 processing |
| `ruina` | Read user input | stdin | Parses mat, err, tempr, errmax, errint; allocates Faddeeva tables if tempr>0 |
| `anlyzd` | Analyze dictionary | MF1/MT451 | Selects redundant reactions (MT1, MT4, MT18, MT27, MT101) |
| `rdfil2` | Read File 2 | MF2 resonance data | Dispatches to formalism-specific readers; builds energy node list |
| `rdf2bw` | BW/RM reader | MF2 LRF=1,2,3 | Reads SLBW, MLBW, Reich-Moore parameters; precomputes shift/penetrability at Er |
| `rdf2aa` | Adler-Adler reader | MF2 LRF=4 | Reads AA coefficients |
| `rdf2hy` | Hybrid R-function | MF2 LRF=6 | Reads hybrid R-matrix parameters |
| `rdf2u0` | Unresolved (E-indep) | MF2 LRU=2, LFW=0 | Energy-independent unresolved parameters |
| `rdf2u1` | Unresolved (fission) | MF2 LRU=2, LFW=1 | Fission-width energy-dependent |
| `rdf2u2` | Unresolved (all) | MF2 LRU=2, LRF=2 | All parameters energy-dependent |
| `order` | Sort energies | energy array | Simple insertion sort of node list |
| `facts` | Penetrability/shift | l, rho | Analytic formulas for l=0..4 hardcoded |
| `facphi` | Phase shift | l, rhoc | Hard-sphere phase shifts for l=0..4 |
| `unfac` | Penetrability (general) | l, rho, rhoc | Combined P, phi for multiple formalisms |
| `genunr` | Unresolved XS table | MF2 -> MT152 | Computes inf-dilute URR cross sections on special grid |
| `sigunr` | URR interpolation | energy -> sig | Interpolates from precomputed URR table |
| `lunion` | Energy grid union | MF3 + MF13 + nodes | Merges all energy grids; linearizes non-linear sections |
| **`resxs`** | **Adaptive reconstruction** | grid + res params -> PENDF | **Stack-based bisection** of energy intervals; tests convergence via fractional + integral criteria |
| **`sigma`** | **Cross section dispatch** | energy -> sig array | Routes to appropriate formalism: SLBW, MLBW, RM, AA, hybrid, SAMMY, URR |
| `csnorp` | No-resonance XS | energy | Returns zero (placeholder for sections without res params) |
| **`csslbw`** | **SLBW cross sections** | e, resonance params | Single-level Breit-Wigner with optional Doppler (psi/chi via quickw) |
| **`csmlbw`** | **MLBW cross sections** | e, resonance params | Multi-level Breit-Wigner (zero temperature only) |
| `csmlbw2` | MLBW+Doppler | e, resonance params | MLBW with Doppler broadening via quickw |
| **`csrmat`** | **Reich-Moore XS** | e, resonance params | R-matrix with 3x3 complex matrix inversion via `frobns` |
| `cshybr` | Hybrid R-function | e, resonance params | Background R-function + resonance terms |
| `csaa` | Adler-Adler XS | e, resonance params | Adler-Adler parameterization |
| `csunr1` | URR type 1 | e, resonance params | Unresolved with energy-independent or fission-dependent widths |
| `csunr2` | URR type 2 | e, resonance params | Unresolved with all widths energy-dependent |
| `gnrl` | Chi-squared sampling | alpha, beta, gamma | Generates chi-squared distributed widths for URR |
| `emerge` | Energy merge | grid + resonance XS | Combines File 3 background with resonance contribution |
| `recout` | Output writer | scratch -> PENDF | Writes final PENDF tape with MF1 dictionary, MF2, MF3 |
| **`quickw`** | **Fast Faddeeva** | x, y -> Re(w), Im(w) | 3-region algorithm: table lookup (x^2+y^2<36), rational approx (<144), asymptotic (>144) |
| **`wtab`** | **Build W table** | -> 62x62 arrays | Precomputes w(x,y) on grid x,y in [-0.1, 5.9] with step 0.1 |
| **`w`** | **Exact Faddeeva** | z -> w(z) | Full complex probability integral; asymptotic or Taylor series depending on region |
| **`frobns`** | **Complex matrix inverse** | A+iB -> C+iD | Frobenius-Schur method for 3x3 complex matrices |
| `thrinv` | Real matrix inverse | D -> D^-1 | Symmetric 3x3 real matrix inversion |
| `abcmat` | Matrix multiply | A, B -> C=AB | 3x3 real matrix multiplication |

### 2.2 BROADR (broadm) -- 1947 lines, ~10 subroutines

Doppler broadens and thins pointwise cross sections.

| Subroutine | Purpose | Algorithm |
|------------|---------|-----------|
| `broadr` | Main driver | Reads input; loops over temperatures; manages 3-page buffering system |
| `stounx` | Store URR data | Copies unresolved resonance data from PENDF |
| `getunx` | Get URR correction | Retrieves unresolved cross sections for energy |
| `bfile3` | Process File 3 | Unionizes low-threshold reactions onto total XS grid |
| **`broadn`** | **Adaptive broadening** | Stack-based bisection (like resxs) ensuring broadened XS converges within tolerance |
| **`bsigma`** | **Doppler kernel** | Exact integral of piecewise-linear XS times Maxwell kernel; uses F-functions (erfc-based) and H-functions (interval integrals) |
| `thinb` | Thinning | Removes redundant points where linear interp suffices |
| **`hunky`** | **H-function** | `h(n,a,b) = f(n,a) - f(n,b)` for n=0..4; falls back to Taylor series when cancellation occurs |
| **`funky`** | **F-function** | `f(n,a) = integral(a to inf) z^n * exp(-z^2) dz / sqrt(pi)` -- analytic in terms of erfc |
| `hnabb` | Taylor H-function | Direct series expansion of h(n,a,b) when b-a is small (avoids cancellation) |

### 2.3 HEATR (heatm) -- 6322 lines, ~20 subroutines

Computes KERMA coefficients and damage energy production.

| Subroutine | Purpose |
|------------|---------|
| `heatr` | Main driver |
| `hinit` | Initialize; read PENDF structure, identify reactions |
| `nheat` | Compute heating numbers for neutron reactions |
| `indx` | Index partial reactions |
| `capdam` | Capture/damage energy for radiative capture |
| `disbar` | Average recoil energy for discrete two-body reactions |
| `df` | Lindhard damage function (partition electronic/nuclear) |
| `conbar` | Average energy for continuum emission |
| `hgtyld` | Get photon/particle yields |
| `anabar` | Analytic average for evaporation spectra |
| `anadam` | Analytic damage energy for evaporation spectra |
| `tabbar` | Tabulated spectrum average energy |
| `tabdam` | Tabulated spectrum damage energy |
| `sixbar` | File 6 format energy deposition |
| `getsix` | Get File 6 distributions |

### 2.4 THERMR (thermm) -- 3387 lines, ~15 subroutines

Thermal neutron scattering: free gas and bound (S(alpha,beta)).

| Subroutine | Purpose |
|------------|---------|
| `thermr` | Main driver |
| `rdelas` | Read elastic thermal data |
| `gateff` | Get effective temperature |
| `coh` | Coherent elastic scattering (Bragg edges); adaptive grid |
| `sigcoh` | Coherent elastic cross section from lattice parameters |
| `form` | Crystal form factor |
| `iel` | Inelastic thermal scattering cross section |
| `terp` | Lagrangian interpolation in S(alpha,beta) tables |
| `calcem` | Compute free-gas or S(alpha,beta) scattering matrix |
| `sig` | Differential scattering cross section from S(alpha,beta) |
| `sigl` | Linear reconstruction of scattering matrix elements |
| `sigu` | Angular distribution reconstruction |

### 2.5 GROUPR (groupm) -- 12590 lines, ~40 subroutines

Self-shielded multigroup cross sections, scattering matrices.

Key subroutines: `groupr` (driver), `ruinb` (input), `nextr` (reaction selector), `namer` (naming), `gengpn`/`gengpg` (group structure generators -- 15+ built-in neutron and 8+ photon group structures), `genwtf` (weight function -- 11 options including constant, 1/E, thermal+1/E+fission, EPRI-cell, CLAW, VITAMIN-E), `getwtf` (evaluate weight), `genflx` (flux calculator), `init` (initialize group arrays).

### 2.6 ACER (acem + 6 sub-modules) -- ~27,000 lines total

Produces ACE format data for MCNP. The `acer.f90` dispatcher calls sub-modules:
- `acefc.f90` (19,669 lines): Fast/continuous neutron data -- the largest single file
- `acedo.f90`: Dosimetry
- `aceth.f90`: Thermal scattering
- `acepa.f90`: Photo-atomic
- `acepn.f90`: Photo-nuclear
- `acecm.f90`: Charged particles

### 2.7 UNRESR (unresm) -- 1665 lines

Bondarenko self-shielding in unresolved range. Key: `unresl` (ladder generation and averaging), `ajku` (J/K integrals), `quikw` (local Faddeeva copy).

### 2.8 PURR (purm) -- 2894 lines

Probability tables for unresolved region (Monte Carlo method). Key: `unresx` (random ladder generation), `ladr2` (ladder construction), `unrest` (probability table from ladder).

### 2.9 ERRORR (errorm) -- 11152 lines, ~50 subroutines

Multigroup covariance processing. Key: `gridd` (energy grid), `covcal` (covariance calculation), `spcint` (spectrum integration).

---

## 3. Test Problem Analysis

### 3.1 Summary Statistics

- **Total tests:** 85 (numbered 01-85, test 77 absent)
- **Most-used modules:** reconr (70 tests), moder (66), broadr (56), acer (37), viewr (34), heatr (13), groupr (17), purr (14), thermr (13), errorr (11), unresr (7), covr (8), gaspr (5), leapr (5), gaminr (2), ccccr (2), matxsr (2), dtfr (1), wimsr (1), plotr (2), mixr (0)

### 3.2 Test-by-Test Catalogue

| Test | Modules Exercised | Description |
|------|-------------------|-------------|
| 01 | moder, reconr, broadr, heatr, thermr, groupr | C-nat: full processing chain with free-gas + bound thermal |
| 02 | moder, reconr, broadr, unresr, groupr, ccccr | Pu-238: 3 temperatures, unresolved, CCCC output |
| 03 | reconr, gaminr, dtfr, matxsr, viewr | Photon interaction: H + U, DTF and MATXS output |
| 04 | moder, reconr, errorr, groupr | U-235: covariance + multigroup nubar |
| 05 | moder, errorr, covr, viewr | C-nat: MF2 covariance processing |
| 06 | plotr, viewr | Cross section plotting only |
| 07 | moder, reconr, broadr, heatr, groupr, acer | Na-23: full chain + ACE |
| 08 | moder, reconr, broadr, heatr, groupr, acer | Be-9: full chain + ACE |
| 09 | moder, reconr, broadr, leapr, thermr | Para-H: thermal scattering with LEAPR |
| 10 | moder, reconr, broadr, unresr, purr, acer | Th-232: unresolved + probability tables + ACE |
| 11 | moder, reconr, broadr, unresr, thermr, groupr, wimsr | U-235: WIMS library production |
| 12 | reconr, gaspr, plotr, viewr | Gas production + plotting |
| 13 | moder, reconr, broadr, gaspr, heatr, acer, viewr | Fe-56: gas production + ACE |
| 14 | acer, viewr | ACE thermal scattering |
| 15 | moder, reconr, broadr, groupr, errorr, covr, viewr | Multi-material covariance |
| 16 | moder, reconr, broadr, groupr, errorr, covr, viewr | Multi-material covariance (variant) |
| 17 | moder, reconr, broadr, groupr, errorr | Covariance with errorr |
| 18 | moder, reconr, broadr, groupr, errorr, covr, viewr | Covariance chain |
| 19 | moder, reconr, broadr, heatr, unresr, purr, acer | Full chain with URR+PURR |
| 20 | reconr, broadr, groupr, errorr, covr, viewr | Covariance chain |
| 21 | reconr, broadr, heatr, purr, acer | Fe-58: ENDF/B-8 full chain |
| 22 | leapr | Para-hydrogen at 20K |
| 23 | leapr | BeO LEAPR |
| 24 | moder, reconr, broadr, heatr, gaspr, unresr, purr, thermr, acer, viewr | Pu-239: comprehensive chain |
| 25 | moder, reconr, broadr, heatr, thermr, acer | Full chain with thermal |
| 26 | moder, reconr, broadr, heatr | Basic reconstruction + broadening + heating |
| 27 | moder, reconr, broadr, groupr, errorr, covr, viewr | Covariance processing |
| 28 | moder, reconr, broadr, acer, viewr | ACE production |
| 29 | moder, reconr, broadr, groupr | Multigroup production |
| 30 | moder, reconr, broadr, gaminr, groupr, matxsr | Photon + neutron multigroup |
| 31 | moder, reconr, broadr, purr | Probability tables |
| 32 | moder, reconr, broadr, thermr | Thermal scattering |
| 33 | leapr | D in D2O |
| 34 | moder, reconr, broadr, groupr, errorr, covr, viewr | Covariance chain |
| 35-42 | moder, reconr, broadr, purr, acer, viewr | PURR test suite (8 tests): various lssf, iinel, iabso options |
| 43-44 | moder, reconr, broadr | BROADR edge cases (thnmax, T=0/0.1K) |
| 45-47 | moder, reconr, broadr, (gaspr/errorr/groupr) | Various processing chains |
| 48 | acer | ACE photo-atomic data |
| 49 | moder, reconr, broadr, heatr, thermr, acer | Full chain |
| 50-54 | moder, acer, viewr | MF6 MT2 LAW=5 tests (angular distributions, identical/different particles) |
| 55-58 | moder, reconr, acer, viewr | ACE production variants |
| 59 | moder, acer | Format conversion + ACE |
| 60 | moder, reconr, broadr, groupr | IRDFF-II processing |
| 61 | acer, viewr | Thermal scattering kernel check |
| 62 | moder, acer, viewr | ACE production |
| 63 | moder, reconr, broadr, purr, acer, viewr | PURR with nunx=2 |
| 64-66 | moder, reconr, acer, viewr | ACE variants |
| 67-70 | moder, reconr, broadr, thermr, acer, viewr | Thermal scattering mode tests (mixed, inelastic only, incoherent, coherent) |
| 71 | moder, acer, viewr | MT11 deuteron production |
| 72 | moder, reconr, broadr, heatr, gaspr, thermr, acer, viewr | Full chain |
| 73 | moder, reconr, broadr, groupr | Multigroup |
| 74 | moder, reconr, broadr, thermr, acer, viewr | Thermal + ACE |
| 75 | reconr, broadr, acer | Ag-109 ACE |
| 76 | moder, reconr, broadr, unresr, groupr, ccccr | CCCC output |
| 78 | moder, reconr, acer, viewr | ACE |
| 79 | moder, reconr, broadr, heatr | Heating |
| 80 | leapr | H in HF |
| 81 | moder, reconr | Basic reconstruction |
| 82 | reconr, broadr, acer | ACE |
| 83 | moder, reconr | Format conversion + reconstruction |
| 84 | reconr | Reconstruction only |
| 85 | moder, reconr | Reconstruction |

### 3.3 Test Classification by Julia Module

For the Julia port, tests map to validation targets:

- **ENDF I/O / MODER:** Tests 81, 83, 85 (pure moder+reconr); all tests use ENDF reading
- **RECONR:** Tests 84 (reconr only), 81, 83, 85 (moder+reconr), plus the reconr component of all multi-module tests
- **BROADR:** Tests 43, 44 (dedicated broadr edge cases), plus broadr component of 56 tests
- **HEATR:** Tests 26, 79 (broadr+heatr), plus heatr component of 13 tests
- **THERMR:** Tests 32, 67-70 (thermal mode variants), plus thermr in 13 tests
- **GROUPR:** 17 tests, often paired with errorr/covr
- **ACER:** 37 tests including dedicated photo-atomic (48), angular distribution (50-54), thermal kernel (14, 61)
- **PURR:** Tests 35-42 (systematic PURR test suite), 10, 31, 63
- **UNRESR:** Tests 02, 10, 11, 19, 24, 47, 76
- **ERRORR/COVR:** Tests 04, 05, 15-18, 20, 27, 34, 46-47, 65
- **LEAPR:** Tests 22, 23, 33, 80 (standalone LEAPR)
- **GASPR:** Tests 12, 13, 24, 45, 72

---

## 4. ENDF Format Mapping

### 4.1 Record Types

The ENDF format uses 80-column card images with specific record structures:

| Record Type | Routine | Structure | Description |
|-------------|---------|-----------|-------------|
| **CONT** | `contio` | `[C1, C2, L1, L2, N1, N2]` (6 fields) | Control record: 2 reals + 4 integers, plus MAT/MF/MT/NS in cols 67-80 |
| **LIST** | `listio` | `CONT header + N1 data values` | List of floating-point values, paged in blocks of `npage` |
| **TAB1** | `tab1io` | `CONT + NR interpolation pairs + NP (x,y) pairs` | Tabulated 1-D function with interpolation law |
| **TAB2** | `tab2io` | `CONT + NR interpolation pairs` | Header for tabulated 2-D function (no data, just metadata) |
| **INTG** | `intgio` | `2 integers + up to 18 small integers` | Compact integer format (for covariance, ndigit controls precision) |
| **HEAD** | same as CONT | `[ZA, AWR, L1, L2, N1, N2]` | Section header (first record of each section) |
| **TPID** | `tpidio` | 66 chars of Hollerith text + MAT/MF/MT | Tape identification |
| **HDAT** | `hdatio` | CONT header + Hollerith description cards | Descriptive text data |
| **DICT** | `dictio` | Series of CONT records | Material directory |
| **SEND** | `asend` | MAT/MF/0 | Section end marker |
| **FEND** | `afend` | MAT/0/0 | File end marker |
| **MEND** | `amend` | 0/0/0 | Material end marker |
| **TEND** | `atend` | -1/0/0 | Tape end marker |

### 4.2 ENDF 11-Column Formatting (a11)

The `a11` subroutine converts floating-point values to the ENDF 11-column string format:
- Normal form: `+1.234567+6` (7 sig figs + 1-digit exponent)
- Reduced: `+1.23456-38` (6 sig figs + 2-digit exponent)
- Extended: `+123456.789` (9 sig figs, F-format, for resonance energies in range 0.1-10^7)

### 4.3 Interpolation Laws (ENDF INT codes)

Used in `terp1`:
- **1** = histogram (y constant)
- **2** = lin-lin (y linear in x)
- **3** = lin-log (y linear in ln(x))
- **4** = log-lin (ln(y) linear in x)
- **5** = log-log (ln(y) linear in ln(x))
- **6** = Coulomb penetrability (charged particles, uses threshold `thr6`)

### 4.4 Key MF/MT Numbers Processed

| MF | Content | Key MTs |
|----|---------|---------|
| 1 | General info | 451 (descriptive data + dictionary), 452 (total nubar), 455 (delayed nubar), 456 (prompt nubar), 458 (fission energy components) |
| 2 | Resonance params | 151 (res params), 152 (URR inf-dilute XS -- NJOY-specific) |
| 3 | Cross sections | 1 (total), 2 (elastic), 4 (inelastic), 5 (other), 16 (n,2n), 17 (n,3n), 18 (fission), 19 (n,f first chance), 20 (n,nf), 21 (n,2nf), 38 (n,3nf), 51-91 (inelastic levels), 102 (n,gamma), 103-107 (n,p/d/t/He3/alpha), 203-207 (gas production), 251 (mubar), 252 (xi), 253 (gamma), 600-849 (charged particle production) |
| 4 | Angular distributions | Same MTs as MF3 |
| 5 | Energy distributions | Same MTs as MF3 |
| 6 | Product energy-angle | Same MTs as MF3 |
| 7 | Thermal scattering | 2 (coherent elastic), 4 (inelastic) |
| 8 | Radioactive products | Fission product yields |
| 10 | Cross sections to isomeric states | |
| 12 | Photon multiplicities | |
| 13 | Photon production XS | |
| 14 | Photon angular dist | |
| 15 | Photon energy dist | |
| 23 | Photo-atomic XS | |
| 27 | Atomic form factors | |
| 33 | Covariance (XS) | |
| 34 | Covariance (angular) | |
| 40 | Covariance (production) | |

### 4.5 Resonance Formalisms (MF2)

| LRU | LRF | Mode | NJOY Handler |
|-----|-----|------|--------------|
| 0 | - | 0 | `csnorp` (no parameters, scattering radius only) |
| 1 | 1 | 1 | `csslbw` (Single-Level Breit-Wigner) |
| 1 | 2 | 2 | `csmlbw`/`csmlbw2` (Multi-Level Breit-Wigner) |
| 1 | 3 | 3 | `csrmat` or `cssammy` (Reich-Moore / R-Matrix Limited) |
| 1 | 4 | 4 | `csaa` (Adler-Adler) |
| 1 | 6 | 6 | `cshybr` (Hybrid R-Function) |
| 1 | 7 | 7 | `cssammy` (R-Matrix Limited via SAMMY) |
| 2 | 1 | 11 | `csunr1` (Unresolved, energy-independent or fission-dependent) |
| 2 | 2 | 12 | `csunr2` (Unresolved, all energy-dependent) |

---

## 5. Global State Audit: RECONR Module

### 5.1 Complete Catalogue of Module-Level Variables (~80 variables)

#### Configuration Parameters (set from user input, constant during processing)

| Variable | Type | Source | Description |
|----------|------|--------|-------------|
| `nin` | integer | user input | ENDF input unit number |
| `nout` | integer | user input | PENDF output unit number |
| `mata` | integer | card 3 | MAT number to process |
| `err` | real(kr) | card 4 | Fractional reconstruction tolerance |
| `tempr` | real(kr) | card 4 | Reconstruction temperature (K) |
| `errmax` | real(kr) | card 4 | Max tolerance when integral criterion met (default 10*err) |
| `errint` | real(kr) | card 4 | Max resonance-integral error per grid point (default err/20000) |
| `ncards` | integer | card 3 | Number of descriptive comment cards |
| `isammy` | integer | hardcoded=1 | Flag for SAMMY method (1=enabled) |
| `nodmax` | parameter=800000 | | Max energy nodes |
| `maxunr` | parameter=500 | | Max unresolved energy points |
| `maxres` | integer | | Max resonance parameter storage |
| `nbufg/nbufr/nbuf/nbufl` | parameter=2000 | | Buffer sizes for loada/finda I/O |

#### Material/Evaluation Metadata (set from ENDF File 1)

| Variable | Type | Description |
|----------|------|-------------|
| `za` | real(kr) | ZA = 1000*Z + A |
| `awr` | real(kr) | Atomic weight ratio to neutron |
| `zain, awin` | real(kr) | Incident particle Z and A (usually 1, 1 for neutrons) |
| `lrp` | integer | Resonance parameter flag (0=none, 1=resolved, 2=urr only, 3=res+urr) |
| `lfi` | integer | Fission flag |
| `lssf` | integer | Self-shielding-only flag for URR |
| `iverf` | integer | ENDF version (4, 5, or 6) -- set in endf module |
| `elis` | real(kr) | Excitation energy of target |
| `sta` | real(kr) | Target stability flag |
| `lis, lis0` | integer | State number, ground state |
| `nfor, lrel, nver` | integer | Format, release, version numbers |
| `efmax` | real(kr) | Upper energy limit |
| `tempi` | real(kr) | Temperature from evaluation |

#### Resonance Parameters and Energy Ranges (set from File 2)

| Variable | Type | Description |
|----------|------|-------------|
| `eresl` | real(kr) | Lowest resonance range lower bound |
| `eresr` | real(kr) | Highest resolved range upper bound |
| `eresh` | real(kr) | Highest resonance range upper bound |
| `eresu` | real(kr) | Lowest unresolved range lower bound |
| `eresm` | real(kr) | Lowest unresolved range upper bound |
| `spin` | real(kr) | Target spin |
| `ascat` | real(kr) | Scattering radius |
| `nsect` | integer | Number of isotope-energy-range sections |
| `isot(20)` | integer array | Isotope index for each section |
| `modet(20)` | integer array | Resonance formalism mode for each section |
| `ibaset(20)` | integer array | Starting index in `res` array for each section |
| `abnt(20)` | real array | Abundance for each section |
| `elt(20), eht(20)` | real arrays | Energy bounds for each section |
| `isect` | integer | Current section being processed |
| `lfw` | integer | Fission width flag |
| `lrx` | integer | Competitive reaction flag |
| `ncoef` | integer | Number of Legendre coefficients (from SAMMY) |
| `nsig` | integer | Number of cross section components (4 or 5+nmtres) |
| `nmtres` | integer | Number of SAMMY resonance reactions |
| `nresp` | integer | Number of SAMMY response functions |
| `mmtres(10)` | integer array | MT numbers for SAMMY reactions |
| `mcards(10)` | integer array | Card counts for SAMMY data |

#### Intermediate Computation State (varies during processing)

| Variable | Type | Description |
|----------|------|-------------|
| `nodes` | integer | Current count of energy nodes |
| `nunr` | integer | Number of unresolved energy points |
| `nxc` | integer | Dictionary entry count |
| `ngo` | integer | Processing flag |
| `mtr4, mtr18` | integer | MT4/MT18 processing flags |
| `mtr(nmtmax), mtrt(nmtmax)` | integer arrays | MT number tracking |
| `nmtr` | integer | Count of tracked MTs |
| `mt103..mt107` | integer | Charged particle reaction flags |
| `mpmin..m4max` | integer | MF/MT range markers |
| `mtr522` | integer | MT522 flag |
| `q18` | real(kr) | Q-value for fission |
| `thr6` | real(kr) | Threshold for Coulomb interpolation |
| `elst(20), slst(20,4)` | real arrays | Last energy/sigma for each section (overlap handling) |
| `enxt(20), snxt(20,4)` | real arrays | Next energy/sigma for each section |
| `ilast(20)` | integer array | Last point flags |

#### Allocated Arrays (I/O buffers and data storage)

| Variable | Type | Description |
|----------|------|-------------|
| `card(:)` | character(4), allocatable | Descriptive comment cards |
| `enode(:)` | real(kr), allocatable | Energy node list (up to nodmax=800,000) |
| `eunr(:)` | real(kr), allocatable | Unresolved energy grid (up to maxunr=500) |
| `res(:)` | real(kr), allocatable | Resonance parameter storage (up to maxres) |
| `tr(62,62), ti(62,62)` | real(kr), allocatable | Faddeeva function lookup tables (real, imaginary) |
| `dict(:)` | real(kr), allocatable | Dictionary data |
| `mfs(:), mts(:), ncs(:)` | integer, allocatable | MF/MT/card-count arrays |
| `sunr(:)` | real(kr), allocatable | Unresolved cross section storage |

### 5.2 Julia Redesign Classification

For the Julia port, these variables should be reorganized into:

1. **Configuration structs** (immutable after setup): `err`, `tempr`, `errmax`, `errint`, `ncards`, `isammy`
2. **Material metadata structs**: `za`, `awr`, `lrp`, `iverf`, `spin`, `ascat`, energy bounds
3. **Resonance data containers** (typed by formalism): `res` array -> proper Julia types per LRF
4. **Section descriptors**: `isot`, `modet`, `ibaset`, `abnt`, `elt`, `eht` -> `Vector{ResonanceSection}`
5. **Computation workspace** (temporary, pass through functions): `nodes`, `enode`, buffers
6. **Output data**: energy grids, cross section arrays

---

## 6. Key Algorithms

### 6.1 Adaptive Grid Construction (Stack-Based Bisection)

Used in both `resxs` (RECONR) and `broadn` (BROADR).

**Algorithm:**
1. Start with an initial coarse energy grid (from File 3 grid points and resonance energies)
2. Maintain a stack of energy points (LIFO, max depth `ndim=30` in reconr, `nstack=12` in broadr)
3. For each pair of adjacent points on the stack:
   a. Compute the midpoint energy `xm = (x_low + x_high) / 2`
   b. Compute the true cross section at `xm` via `sigma()` or `bsigma()`
   c. Compute the linearly interpolated value at `xm`
   d. Test convergence using a 3-tier criterion:
      - **Primary:** `|sigma_true - sigma_interp| <= err * sigma_true` for all reactions
      - **Relaxed:** If primary fails but `|diff| <= errmax * sigma_true` AND `|diff| * dx / (2*xm) <= errint`, accept (resonance integral check)
      - **Step limit:** Reject if energy step is > 4.1x the previous step (guards against missed peaks)
   e. If converged: pop and write the top point; continue with remaining stack
   f. If not converged: push the midpoint onto the stack (bisect the interval)
4. Continue until stack is exhausted, then move to next coarse grid interval
5. Output the adaptively refined grid with cross sections

**Key parameters:**
- `err` = fractional tolerance (default typically 0.001-0.01)
- `errmax` = relaxed tolerance (default 10*err)
- `errint` = integral error tolerance (default err/20000)
- `ndim=30` = maximum stack depth (limits refinement to ~2^30 subdivision)
- `sigfig(x,7,0)` = round to 7 significant figures (prevents infinite subdivision)

**Convergence near thermal:** Tolerance is tightened by factor of 5 when energy < 0.5 eV (`trange=0.4999`).

### 6.2 Faddeeva Function (Complex Probability Integral)

Three-level implementation for `w(z) = exp(-z^2) * erfc(-iz)`:

**Level 1: `w(rez, aim1)` -- Exact computation**
- Region selection based on |Re(z)|, |Im(z)| using 5 break conditions
- Two algorithms:
  - **Asymptotic series** (Im(z) >= 0, far from origin): continued-fraction-like recursion with overflow/underflow scaling (factors of 10^15)
  - **Taylor series** (Im(z) < 0 or intermediate region): expansion with exponential prefactor, convergence tested to eps=10^-7

**Level 2: `wtab(rw, aimw)` -- Lookup table construction**
- Generates 62x62 grid of w(x,y) values
- Grid: x in [-0.1, 5.9], y in [-0.1, 5.9], step = 0.1
- Uses `w()` to fill each grid point

**Level 3: `quickw(ax, y, rew, aimw)` -- Fast evaluation**
- For `x^2 + y^2 < 36`: Bilinear interpolation from 62x62 table (6 coefficients: a1-a5 + pq term)
- For `36 <= x^2 + y^2 < 144`: 2-pole rational approximation (constants c1-c4)
- For `144 <= x^2 + y^2 < 10000`: 1-pole rational approximation (constant c5=2/sqrt(pi))
- For `x^2 + y^2 >= 10000`: Simple `1/(sqrt(pi) * (x^2+y^2))` asymptotic

**Physics connection:** The Faddeeva function gives the Doppler-broadened line shape:
- `psi = sqrt(pi) * theta * Re(w) / 2` (symmetric/Voigt profile)
- `chi = sqrt(pi) * theta * Im(w) / 2` (asymmetric profile)
where `theta = Gamma_total / delta`, `delta = sqrt(4*kT*E/A)` is the Doppler width.

### 6.3 Reich-Moore Matrix Inversion (frobns)

**Purpose:** Invert the complex matrix `(I - R*L)` in the R-matrix collision equation where R is the R-matrix and L is the boundary-condition-dependent matrix.

**Method:** Frobenius-Schur decomposition for 3x3 complex matrices `M = A + iB`:
1. Compute `A^{-1}` via `thrinv` (symmetric real matrix inverse via Gauss elimination)
2. Compute `Q = A^{-1} * B`
3. Compute `D = B * Q` (adds to real part)
4. Compute `C = (A + D)^{-1}` via `thrinv` again
5. Compute `D_inv = -Q * C` (imaginary part of inverse)

The `thrinv` subroutine inverts a symmetric 3x3 real matrix using an in-place method:
1. Transform `D -> I - D`
2. Loop over pivots with `D[lr,lr] = 1/(1 - D[lr,lr])`
3. Update off-diagonal elements

**Why 3x3:** Reich-Moore has at most 3 channels per J-value: neutron + 2 fission channels (partial fission widths gamma_fa and gamma_fb).

### 6.4 Doppler Broadening Kernel (BROADR)

**Method:** Exact analytic Doppler broadening of piecewise-linear cross sections (SIGMA1 method by D.E. Cullen).

The broadened cross section is:
```
sigma_b(v) = (1/sqrt(pi)) * integral sigma(v') * |v'| * exp(-(v-v')^2) dv'
```

where velocities are in units of the Doppler width `v_D = sqrt(2kT/(m_n * A))`.

**Implementation in `bsigma`:**
1. Convert energies to velocity-like variable: `v = sqrt(alpha * E)` where `alpha = AWR / (bk * T_eff)`
2. For each segment of the piecewise-linear cross section:
   - Express `sigma(v')` as `sigma_low + slope * (v'^2 - v_low^2)`
   - Compute contributions using **F-functions**: `f_n(a) = integral(a to inf) z^n * exp(-z^2) dz / sqrt(pi)`
   - And **H-functions**: `h_n(a,b) = f_n(a) - f_n(b)` (difference between F-functions at interval endpoints)
3. The broadened cross section involves two integral terms:
   - `s1 = h_2/y^2 + 2*h_1/y + h_0` (constant term contribution)
   - `s2 = (h_4 + (6y^2-x^2)*h_2)/y^2 + (4*h_3 + (4y^2-2x^2)*h_1)/y + (y^2-x^2)*h_0` (slope contribution)
4. Sum contributions from all intervals below and above the evaluation point
5. Extend as 1/v to E=0 and as constant to E=infinity

**F-functions** (in `funky`):
```
f_0(a) = erfc(a) / 2
f_1(a) = exp(-a^2) / (2*sqrt(pi))
f_2(a) = (f_0 + a * exp(-a^2)/sqrt(pi)) / 2
f_3(a) = (2*f_1 + a^2 * exp(-a^2)/sqrt(pi)) / 2
f_4(a) = (3*f_2 + a^3 * exp(-a^2)/sqrt(pi)) / 2
```

**H-function cancellation handling** (in `hnabb`):
When `b - a` is small, `h_n(a,b) = f_n(a) - f_n(b)` suffers catastrophic cancellation.
The `hnabb` function uses a direct Taylor series expansion of the defining integral with
coefficient recursion, converging to relative precision 10^-8.

### 6.5 Penetrability and Shift Factor Calculations

**`facts(l, rho, se, pe)`** -- Computes P_l(rho) and S_l(rho) analytically:

| l | P_l(rho) | S_l(rho) |
|---|----------|----------|
| 0 | rho | 0 |
| 1 | rho^3 / (1 + rho^2) | -1 / (1 + rho^2) |
| 2 | rho^5 / (9 + 3*rho^2 + rho^4) | -(18 + 3*rho^2) / (9 + 3*rho^2 + rho^4) |
| 3 | rho^7 / (225 + 45*rho^2 + 6*rho^4 + rho^6) | -(675 + 90*rho^2 + 6*rho^4) / D |
| 4 | rho^9 / (11025 + 1575*rho^2 + 135*rho^4 + 10*rho^6 + rho^8) | -(44100 + 4725*rho^2 + ...) / D |

where `rho = k * a` (wavenumber times channel radius).

**`facphi(l, rho, phi)`** -- Hard-sphere phase shifts:
- l=0: phi = rho
- l=1: phi = rho - atan(rho)
- l=2: phi = rho - atan(3*rho / (3 - rho^2))
- l=3,4: similar with higher-order rational functions

**Channel radius:** `a = 0.123 * (A_neutron * AWR)^{1/3} + 0.08` fm (default), or `a = AP` if NAPS=1.

**Wavenumber:** `k = cwaven * AWR/(AWR+1) * sqrt(|E|)` where `cwaven = sqrt(2 * m_n * amu * eV) * 10^{-12} / hbar`

### 6.6 Cross Section Formulas

**Single-Level Breit-Wigner (csslbw):**
- Elastic: `sigma_el = sum_r { 4*pi*g_J*Gamma_n / (Gamma_tot)^2 * [(cos(2*phi)*Gamma_tot - Gamma_x)*psi + sin(2*phi)*chi*Gamma_tot] } + sigma_pot`
- Fission: `sigma_f = sum_r { 4*pi*g_J*Gamma_n*Gamma_f / (Gamma_tot)^2 * psi }`
- Capture: `sigma_gamma = sum_r { 4*pi*g_J*Gamma_n*Gamma_gamma / (Gamma_tot)^2 * psi }`
- `psi` and `chi` are Doppler-broadened line shapes (Voigt profiles) when T>0, or standard Lorentzians when T=0

**Multi-Level Breit-Wigner (csmlbw):**
- Includes interference between resonances of the same J,l
- Total constructed as sum of parts

**Reich-Moore (csrmat):**
- Constructs 3x3 R-matrix per J-value
- Inverts `(I - R*L)` using Frobenius-Schur method
- Collision matrix: `U = exp(2i*phi) * (I - R*L)^{-1} * (I - R*L^*)`
- Cross sections from `U` via optical theorem

---

## 7. Architecture Recommendations for Julia Port

### 7.1 Recommended Julia Type Hierarchy

```julia
# Foundation types
abstract type ENDFRecord end
struct ContRecord <: ENDFRecord ... end
struct Tab1Record <: ENDFRecord ... end
struct Tab2Record <: ENDFRecord ... end
struct ListRecord <: ENDFRecord ... end

# Interpolation
@enum InterpolationLaw Histogram=1 LinLin=2 LinLog=3 LogLin=4 LogLog=5 CoulombPen=6

struct InterpolationTable
    breakpoints::Vector{Int}
    laws::Vector{InterpolationLaw}
end

struct TabulatedFunction
    interp::InterpolationTable
    x::Vector{Float64}
    y::Vector{Float64}
end

# Physics constants (module)
module PhysicsConstants
    const NEUTRON_MASS_AMU = 1.00866491595  # etc.
end

# Material identification
struct MaterialId
    za::Float64          # 1000*Z + A
    awr::Float64         # Atomic weight ratio
    mat::Int             # MAT number
end

# Resonance parameter types (replacing the flat `res` array)
abstract type ResonanceFormalism end

struct BreitWignerParams <: ResonanceFormalism
    er::Float64          # Resonance energy
    aj::Float64          # Spin J
    gn::Float64          # Neutron width
    gg::Float64          # Gamma width
    gf::Float64          # Fission width
    gx::Float64          # Competitive width
end

struct ReichMooreParams <: ResonanceFormalism
    er::Float64
    aj::Float64
    gn::Float64
    gg::Float64
    gfa::Float64         # Fission width A
    gfb::Float64         # Fission width B
end

struct ResonanceRange
    el::Float64          # Lower energy bound
    eh::Float64          # Upper energy bound
    lru::Int             # Resolved(1)/Unresolved(2) flag
    lrf::Int             # Formalism flag
    naps::Int            # Scattering radius flag
    spin::Float64        # Target spin
    ap::Float64          # Scattering radius
    abundance::Float64   # Isotopic abundance
    resonances::Vector{<:ResonanceFormalism}
end

# Cross section result
struct CrossSections
    total::Float64
    elastic::Float64
    fission::Float64
    capture::Float64
    competitive::Float64
end

# Adaptive grid configuration
struct ReconstructionConfig
    err::Float64
    errmax::Float64
    errint::Float64
    tempr::Float64       # Reconstruction temperature
end

# Doppler broadening configuration
struct BroadeningConfig
    errthn::Float64
    errmax::Float64
    errint::Float64
    thnmax::Float64
    alpha::Float64       # AWR / (bk * T_eff)
end
```

### 7.2 Recommended Module Structure

```
NJOY.jl/
  src/
    NJOY.jl                    # Main module
    constants.jl               # Physics constants (from phys.f90)
    endf/
      types.jl                 # Record types, interpolation
      io.jl                    # Reading/writing ENDF format
      interpolation.jl         # terp1, terpa, intega, gral
      navigation.jl            # findf, skip6, tosend etc.
    resonances/
      types.jl                 # ResonanceRange, parameter types
      penetrability.jl         # facts, facphi, unfac
      breit_wigner.jl          # SLBW, MLBW cross sections
      reich_moore.jl           # R-matrix, frobns
      adler_adler.jl           # AA formalism
      hybrid.jl                # Hybrid R-function
      unresolved.jl            # URR cross sections
      faddeeva.jl              # w(z), quickw, wtab
      sammy.jl                 # SAMMY RML method
    processing/
      reconr.jl                # Pointwise reconstruction
      broadr.jl                # Doppler broadening (SIGMA1 kernel)
      heatr.jl                 # KERMA/damage
      thermr.jl                # Thermal scattering
      groupr.jl                # Multigroup
      acer.jl                  # ACE format output
      unresr.jl                # Bondarenko self-shielding
      purr.jl                  # Probability tables
      errorr.jl                # Covariances
    adaptive/
      grid.jl                  # Stack-based bisection algorithm
    utils/
      formatting.jl            # a11, a10, sigfig
      buffered_io.jl           # loada/finda replacement
```

### 7.3 Key Design Principles for Julia Port

1. **Eliminate global state:** Replace module-level variables with explicit function arguments and configuration structs. The ~80 globals in RECONR become fields of 4-5 structs passed explicitly.

2. **Type-dispatch for formalisms:** Use Julia's multiple dispatch instead of mode integers. Each resonance formalism gets its own type, and `sigma(e, range::BreitWignerSLBW)` dispatches automatically.

3. **Replace flat `res` array with typed structs:** The Fortran `res(:)` array stores all resonance data as a flat vector with manual indexing (`res(inow)`, `res(inow+1)`, etc.). This should become `Vector{ResonanceRange}` containing `Vector{BreitWignerParams}` or similar.

4. **Replace I/O-unit communication with in-memory data:** Fortran modules communicate through temporary files (units 10-19). In Julia, pass data structures directly between processing stages.

5. **Preserve the adaptive grid algorithm exactly:** The stack-based bisection with the 3-tier convergence criterion is well-tested and numerically important. Port it faithfully using a `Vector` as a stack.

6. **Use Julia's `SpecialFunctions.erfc` and `erfi`:** For the Faddeeva function, Julia's ecosystem provides optimized implementations. However, the `quickw` table-lookup approach may still be faster for the specific use case (many evaluations at small x,y). Benchmark both.

7. **Matrix operations:** Replace `frobns` 3x3 fixed-size inversion with `StaticArrays.SMatrix{3,3}` and Julia's built-in `inv()`. This should be faster and more readable.

8. **Maintain ENDF compatibility:** The I/O format (80-column cards, 11-character fields) must be preserved exactly for interoperability with other nuclear data codes.

9. **Testing strategy:** The 85 test problems provide a comprehensive regression test suite. Each test has a reference `output` file and reference tapes. Port tests progressively:
   - Wave 1: ENDF I/O (test against moder functionality)
   - Wave 2: RECONR (tests 84, 81, 83, 85)
   - Wave 3: BROADR (tests 43, 44)
   - Wave 4: Full pipeline tests (test 01, 07, etc.)

10. **Buffered I/O replacement:** The `loada`/`finda` system implements a buffered sequential access pattern over binary files. In Julia, this should be replaced with `Vector{NTuple{N,Float64}}` or similar in-memory structures, with optional serialization for very large datasets.

---

## Appendix A: ENDF Module Public Interface Summary

```fortran
! Record I/O
contio(nin, nout, nscr, a, nb, nw)
listio(nin, nout, nscr, a, nb, nw)
tab1io(nin, nout, nscr, a, nb, nw)
tab2io(nin, nout, nscr, a, nb, nw)
moreio(nin, nout, nscr, a, nb, nw)    ! continuation pages
tpidio(nin, nout, nscr, a, nb, nw)    ! tape identification
hdatio(nin, nout, nscr, a, nb, nw)    ! Hollerith data
dictio(nin, nout, nscr, a, nb, nw)    ! dictionary
intgio(nin, nout, nscr, a, nb, nw)    ! compact integer

! Section/file delimiters
tosend, tofend, tomend, totend         ! skip/copy through SEND/FEND/MEND/TEND
asend, afend, amend, atend             ! write delimiter cards

! Navigation
findf(mat, mf, mt, ntape)             ! find section on tape
skip6(nin, nout, nscr, a, law)        ! skip File 6 subsection

! Interpolation
terp1(x1,y1,x2,y2,x,y,i)             ! single point interpolation
terpa(y,x,xnext,idis,a,ip,ir)        ! interpolate in packed TAB1
intega(f,x1,x2,a,ip,ir)              ! integrate packed TAB1
gral(xl,yl,xh,yh,x1,x2,int)         ! panel integral (analytic)

! Sequential access
gety1(x,xnext,idis,y1,itape,a)       ! read TAB1 sequentially (stream 1)
gety2(x,xnext,idis,y2,itape,a)       ! read TAB1 sequentially (stream 2)

! Formatting
a11(x, hx)                            ! real -> 11-char ENDF string
```

## Appendix B: Resonance Cross Section Computation Flow

```
reconr()
  |-> ruina()         -- read user input, build Faddeeva table if T>0
  |-> anlyzd()        -- analyze MF1 dictionary
  |-> rdfil2()        -- read MF2, dispatch to formalism readers
  |     |-> rdf2bw()  -- SLBW/MLBW/RM parameter reading + node generation
  |     |-> rdsammy() -- SAMMY RML parameter reading (if LRF=3,7 + isammy)
  |     |-> rdf2aa()  -- Adler-Adler
  |     |-> rdf2hy()  -- Hybrid R-function
  |     |-> rdf2u*()  -- Unresolved formats
  |
  |-> genunr()        -- compute URR infinite-dilute table (MT152)
  |-> lunion()        -- merge energy grids from MF3 + MF13 + nodes
  |
  |-> resxs()         -- ADAPTIVE GRID RECONSTRUCTION
  |     |             -- for each coarse grid interval:
  |     |->  sigma(e) -- compute cross sections
  |     |      |-> csslbw(e) / csmlbw(e) / csrmat(e) / cssammy(e) / ...
  |     |      |       |-> facts(l,rho) -- penetrability
  |     |      |       |-> facphi(l,rho) -- phase shift
  |     |      |       |-> quickw(x,y) -- Faddeeva (if T>0)
  |     |      |       |-> frobns(A,B) -- complex inverse (RM only)
  |     |      |
  |     |      |-> [sum over isotope sections with abundances]
  |     |
  |     |-- check convergence (err, errmax, errint criteria)
  |     |-- if not converged: push midpoint onto stack
  |     |-- if converged: pop and write point
  |
  |-> emerge()        -- merge resonance XS with File 3 background
  |-> recout()        -- write PENDF output
```
