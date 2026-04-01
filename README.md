# NJOY.jl

A Julia port of [NJOY2016](https://github.com/njoy/NJOY2016), the standard
nuclear data processing system used worldwide for reactor physics, criticality
safety, and radiation transport. The original is ~120,000 lines of Fortran 90;
this port is ~15,000 lines of idiomatic Julia.

**Status**: Active development toward bit-identical output on all 84 NJOY
reference tests. RECONR is the most mature module (19 tests bit-identical).
The full T01 pipeline (reconr+broadr+heatr+thermr) passes at 1e-5 tolerance.

## Test Results

### RECONR (resonance reconstruction) -- 19/84 BIT-IDENTICAL

All 84 tests run without crashes. 19 produce byte-identical MF3 output:

| Status | Tests | Notes |
|--------|-------|-------|
| **BIT-IDENTICAL** (19) | T01-03, T08-13, T18-19, T25-27, T30, T45, T47, T55, T84 | All formalisms: LRU=0, SLBW, MLBW, Reich-Moore, SAMMY/RML, URR modes 11+12 |
| Near-perfect (9) | T04, T07, T15-17, T20, T34, T46, T49 | 89-99% MTs perfect; remaining diffs are +/-1 ULP FP precision |
| Partial (2) | T21, T65 | Grid density or URR boundary diffs |
| Feature gap (5) | T56-58, T60, T64 | Photonuclear MF23 partially implemented; dosimetry MF10-only not yet |
| RAN_OK (31) | T24, T28-29, T31-32, T35-44, T63, T66-70, T72-75, T78-79, T81-83, T85 | Run without oracle comparison |
| Non-RECONR (18) | T05-06, T14, T22-23, T33, T48, T50-54, T59, T61-62, T71, T76, T80 | Exercised through their respective modules |

### Full Pipeline -- T01 passes at 1e-5

T01 (C-nat, reconr+broadr+heatr+thermr x2) produces 32,962 lines matching the
Fortran reference. All 41 sections structurally match. 355 value diffs at 1e-9,
0 at 1e-5. Run via:

```bash
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/01/input"; work_dir="/tmp/t01")'
```

### T02 Broadr -- 3-temperature sequential broadening

T02 (Pu-238, reconr+broadr at 300/900/2100K) grids match oracle exactly
(2925/2592/2418 points). 0 failures at 1e-5 on MF3 data. Run via:

```bash
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/02/input"; work_dir="/tmp/t02")'
```

### T03 Photoatomic -- MF=23 support

T03 (photoatomic H + U, multi-material reconr) produces MF=23 output with
99.5% byte-identical data (2102/2112 lines). 2 failures at 1e-5 from grid
construction path difference. Run via:

```bash
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/03/input"; work_dir="/tmp/t03")'
```

## Architecture

### Orchestration layer (`src/orchestration/`)

`run_njoy(input_path)` parses any Fortran NJOY input deck and dispatches
modules sequentially. Modules communicate exclusively via tape files (no
shared mutable state), matching the Fortran architecture exactly.

```julia
using NJOY
tapes = run_njoy("path/to/input"; work_dir="/tmp/work")
```

### Module status

| Module | Status | Julia files |
|--------|--------|-------------|
| **reconr** | Production -- 19 tests bit-identical | `src/processing/reconr.jl`, `reconr_grid.jl`, `reconr_evaluator.jl` |
| **broadr** | Production -- multi-temp, sigma1 kernel correct | `src/processing/broadr.jl`, `sigma1.jl` |
| **heatr** | Production -- KERMA + damage with MF4/MF12/MF13 | `src/processing/heatr.jl` |
| **thermr** | Production -- free gas + S(a,b) + Bragg edges | `src/processing/thermr.jl` |
| **moder** | Working -- tape copy/register | `src/orchestration/modules/moder.jl` |
| **unresr** | Algorithm implemented, not wired into pipeline | `src/processing/unresr.jl` |
| **groupr** | Algorithm implemented, not wired into pipeline | `src/processing/groupr.jl` |
| **acer** | ACE format reader/writer implemented | `src/formats/ace.jl` |
| **errorr** | Covariance processing implemented | `src/processing/errorr.jl` |
| **leapr** | S(a,b) generation implemented | `src/processing/leapr.jl` |
| **ccccr** | CCCC binary format writer implemented | `src/formats/ccccr.jl` |
| **gaspr, purr, dtfr, matxsr, plotr, viewr, mixr, powr, wimsr, covr, resxsr** | Stub/not yet implemented | -- |

### Resonance formalisms

| Formalism | Status | Tests |
|-----------|--------|-------|
| LRU=0 (no resonances) | Bit-identical | T01, T84, T25, T45 |
| SLBW (LRF=1) | Bit-identical | T02, T10-11, T18-19 |
| MLBW (LRF=2) | Implemented, not separately tested | -- |
| Reich-Moore (LRF=3) | Bit-identical | T08, T12-13, T27, T47, T55 |
| SAMMY/RML (LRF=7) | 98% (4 MTs with +/-1 FP) | T20 |
| URR mode 11 (LRU=2, LRF<=1) | Bit-identical | T02, T10-11 |
| URR mode 12 (LRU=2, LRF=2) | 89% (3 MTs +/-1 at URR boundary) | T04, T07 |

## Installation

Requires Julia 1.10+. Not registered; install from source:

```julia
using Pkg
Pkg.develop(path="/path/to/NJOY.jl")
```

The Fortran NJOY2016 reference binary (for oracle generation and gdb
diagnostics) is in `njoy-reference/build/njoy`. Rebuild with:

```bash
cd njoy-reference/build && cmake --build . --target njoy
```

## Key files

| File | Purpose |
|------|---------|
| `HANDOFF.md` | Detailed session history, traps, and grind state for agents |
| `src/orchestration/pipeline.jl` | `run_njoy()` entry point |
| `src/processing/reconr.jl` | Top-level RECONR pipeline |
| `src/processing/broadr.jl` | Doppler broadening (broadn_grid + sigma1) |
| `src/processing/heatr.jl` | KERMA and damage energy |
| `src/processing/thermr.jl` | Thermal scattering (calcem + sigl + Bragg) |
| `src/processing/pendf_writer.jl` | PENDF output (write_full_pendf, write_broadr_pendf) |
| `test/validation/t01_pipeline.jl` | Hand-wired T01 full pipeline test |
| `test/validation/t02_pipeline.jl` | Hand-wired T02 broadr test |
| `test/validation/sweep_all.jl` | All 84 tests in one run |

## Running tests

```bash
# Full 84-test sweep (reconr + module stubs)
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/sweep_all.jl

# Single test via orchestration
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/01/input"; work_dir="/tmp/t01")'

# Hand-wired pipeline tests
julia --project=. test/validation/t01_pipeline.jl
julia --project=. test/validation/t02_pipeline.jl
```

**Important**: Always clear the precompilation cache before running tests:
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
```

## License

GPL-3.0. See LICENSE file.

## References

1. R.E. MacFarlane et al., *The NJOY Nuclear Data Processing System, Version
   2016*, LA-UR-17-20093, Los Alamos National Laboratory (2016).
2. A. Trkov et al. (Eds.), *ENDF-6 Formats Manual*, BNL-90365-2009 Rev. 2,
   Brookhaven National Laboratory (2018).
