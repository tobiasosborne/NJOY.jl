# NJOY.jl

A Julia port of [NJOY2016](https://github.com/njoy/NJOY2016) — the standard
nuclear data processing system used worldwide for reactor physics, criticality
safety, and radiation transport. The Fortran reference is ~100 k lines of
Fortran 90 in `njoy-reference/src/`; this port is ~24 k lines of idiomatic
Julia in `src/` (cloc, code-only, as of Phase 51).

**Goal**: bit-identical output on all 84 of NJOY's own reference tests. Same
ENDF input, same `tapeNN` output, byte-for-byte (or within published
cross-compiler tolerance; see `reports/ACCEPTANCE_CRITERIA.md`).

## Status at a glance

- **RECONR** is the most mature module — **19 reference tests are
  bit-identical** on the RECONR-produced PENDF, covering all resonance
  formalisms (LRU=0, SLBW, MLBW, Reich-Moore, SAMMY/RML, URR modes 11+12).
  See HANDOFF.md §"Current State" for the test list.
- **Full pipelines** run end-to-end on most reference decks. T01 (C-nat,
  reconr+broadr+heatr+thermr) passes at 1e-5 tolerance (32 962 lines,
  41/41 sections structurally match). T04 (U-235) errorr MF31 tape23 is
  at NUMERIC_PASS 81/82 lines.
- **Active front**: `errorr` covariance. T15/T17 (U-238 JENDL) MF33
  matrix values now match Fortran covcal exactly for LB=5 blocks after
  Phase 51's σ·flx-weighted union-grid collapse fix.
- **Unwired algorithms**: `leapr` (S(α, β) phonon expansion), `purr`
  (probability tables via Monte-Carlo ladders), `covr` (correlation
  matrix + boxer output), `gaspr` (MT203-207 gas-production), and the
  CCCC ISOTXS/BRKOXS/DLAYXS writers each have a working implementation
  in `src/processing/` or `src/formats/` (~200-280 LOC of real physics
  per module), but the orchestration wrappers in
  `src/orchestration/modules/` just copy/touch tapes — the connecting
  plumbing (input-deck parsing, tape splicing, ENDF/CCCC emit) has not
  been wired up.
- **Other gaps**: `plotr` (plot-command tape emission) has no backing
  implementation yet. `thermr` supports free-gas and S(α, β) but
  hard-codes the Debye-Waller integral to graphite. `acer` covers
  `iopt=1` fast-neutron ACE; thermal/photo/dosimetry ACE paths are
  stubs.

See [HANDOFF.md](HANDOFF.md) for the living project state (current phase,
open work, traps, per-test status) and `worklog/T*.md` for per-session
debug journals.

## Installation

Requires Julia **1.10+** (`Project.toml` line 17). Not registered; install
from source:

```julia
using Pkg
Pkg.develop(path="/path/to/NJOY.jl")
```

The Fortran NJOY2016 reference binary is bundled for oracle generation
and gdb diagnostics. Rebuild with:

```bash
cd njoy-reference/build && cmake --build . --target njoy
```

## Quick start

Drop-in replacement for the Fortran CLI — parses an NJOY input deck and
runs the whole pipeline:

```julia
using NJOY
run_njoy("njoy-reference/tests/01/input"; work_dir="/tmp/t01")
```

Modules communicate exclusively via tape files on disk (no shared mutable
state), matching Fortran's `tape{N}` plumbing exactly. Tape paths are
managed by a `TapeManager`.

## Running tests

Before every run (**non-negotiable** — concurrent Julia corrupts the
precompile cache):

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
```

### Single reference test

Drives a specific NJOY reference input deck through `run_njoy` and diffs
every produced `tapeNN` against the Fortran `referenceTapeNN` at four
tolerances (1e-9, 1e-7, 1e-5, 1e-3):

```bash
julia --project=. test/validation/reference_test.jl 7       # T07
julia --project=. test/validation/reference_test.jl 15 17   # T15 and T17
```

### Full 84-test sweep

Runs everything and writes a Markdown report to
`reports/REFERENCE_SWEEP.md`. ~90–150 minutes depending on broadr
performance on heavy isotopes:

```bash
julia --project=. test/validation/sweep_reference_tests.jl
```

### Unit tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Architecture

### Tape-based orchestration

`run_njoy(input_path)` tokenises a Fortran-style input deck, parses each
module's card block, and dispatches sequentially. All state flows through
`tape{N}` files — no shared data structures between modules. The dispatch
table lives in `src/orchestration/pipeline.jl`.

### Module status

Classification refers to the **orchestration wrapper** in
`src/orchestration/modules/`, which is what `run_njoy` actually invokes:
- **Real** — wrapper dispatches to substantial processing code; output
  matches Fortran to a documented tolerance.
- **Partial** — wrapper dispatches real code, but key features or
  reaction classes are missing; may degrade gracefully on unsupported
  inputs.
- **Unwired** — algorithm is implemented in `src/processing/` or
  `src/formats/`, but the orchestration wrapper is a tape copy / empty
  touch / marker file; the connecting plumbing (input-deck parsing,
  tape splicing, ENDF emit) has not been wired up.

| Module | Class | Backing files | Notes |
|--------|-------|---------------|-------|
| `reconr` | Real | `src/processing/reconr*.jl` | 19 RECONR-level bit-identical tests |
| `broadr` | Real | `src/processing/broadr.jl`, `sigma1.jl` | Multi-temperature sequential broadening |
| `unresr` | Real | `src/processing/unresr.jl` | LFW=0/1, LRF=2 Bondarenko |
| `heatr` | Real | `src/processing/heatr.jl` | KERMA + damage; plot-tape stubbed |
| `errorr` | Real | `src/processing/errorr.jl`, `src/orchestration/modules/errorr.jl` | MF31/33 NI + NC LTY=0/v2 + GENDF readback; LTY=1/2/3 standards/ratio return `nothing` |
| `moder` | Partial | `src/processing/moder.jl` | Tape copy; multi-material merge not wired |
| `thermr` | Partial | `src/processing/thermr.jl` | Free-gas + S(α, β); Debye-Waller hardcoded to graphite 2.1997 |
| `groupr` | Partial | `src/processing/groupr.jl` | MF=3 and MF=1 nubar only; transfer matrices (MF=6/8/10/16-36) not ported |
| `acer` | Partial | `src/formats/ace_*.jl` | `iopt=1` fast neutron real; thermal/photoatomic/dosimetry stubs |
| `gaminr` | Partial | `src/processing/gaminr.jl` | MF=23 Gauss-Lobatto group averaging |
| `dtfr` | Partial | `src/formats/dtfr.jl` | DTF table + 3D photon-scatter plot tape; accuracy unverified |
| `matxsr` | Partial | `src/formats/matxsr*.jl` | MATXS record layout; scatter matrices P0-only |
| `viewr` | Partial | `src/viewr/*.jl` | Renders real tapes when upstream plot module populated |
| `purr` | Unwired | `src/processing/purr.jl` (255 LOC) | `generate_ptable` (MC ladders + chi²/Wigner sampling + Doppler + Bondarenko self-shielding) ported; wrapper copies PENDF, MT152 emission not wired |
| `leapr` | Unwired | `src/processing/leapr.jl` (264 LOC) | `generate_sab` (phonon expansion, convolution, discrete oscillators, SCT fallback) ported; wrapper touches empty file, MF7 emission not wired |
| `covr` | Unwired | `src/processing/covr.jl` (283 LOC) | Correlation conversion + relative std dev + boxer format ported; wrapper touches empty file, plotr-format emission not wired |
| `gaspr` | Unwired | `src/processing/gaspr.jl` (252 LOC) | `gas_production` + `gas_multiplicity` (MT203-207) ported; wrapper copies PENDF, MT20x splicing not wired |
| `ccccr` | Unwired | `src/formats/ccccr*.jl` (399 LOC) | `write_isotxs` / `write_brkoxs` / `write_dlayxs` (CCCC-IV binary record format) ported; wrapper writes marker files, GENDF parsing not wired |
| `plotr` | Stub | — | No backing implementation yet; wrapper touches empty file |

### Resonance formalisms (RECONR)

| Formalism | Reader | Evaluator | Status |
|-----------|--------|-----------|--------|
| LRU=0 (no resonances) | — | — | Bit-identical (T01, T84) |
| SLBW (LRF=1) | `_read_bw_params` | `cross_section_slbw` | Bit-identical (T02) |
| MLBW (LRF=2) | `_read_bw_params` | `cross_section_mlbw` | Runs; T49 near-perfect (44/46 MTs) |
| Reich-Moore (LRF=3) | `_read_rm_params` | `cross_section_rm` (Frobenius-Schur) | Bit-identical (T08, T27, T47) |
| SAMMY/RML (LRF=7) | `_read_sammy_params` | `build_rml_evaluator` | Runs |
| URR mode 11 (LFW=1) | `_read_urr_lfw1` | `_csunr1`, `_gnrl` | Bit-identical (T02, T18) |
| URR mode 12 (LRF=2) | `_read_urr_lrf2` | `_csunr2`, `_gnrl` | Runs; ±1 ULP at URR boundary (T07) |

See HANDOFF §"Current State" for detailed per-test results.

## Key files

| File | Purpose |
|------|---------|
| [CLAUDE.md](CLAUDE.md) | Non-negotiable rules: oracle-driven TDD, Fortran-before-Julia, grind method |
| [HANDOFF.md](HANDOFF.md) | Living state: current phase, open work, traps, per-test status |
| `worklog/T*.md` | Per-session debug journals (most recent: Phase 51 LB=5 covcal fix) |
| `reports/ACCEPTANCE_CRITERIA.md` | Tolerance hierarchy (1e-9 stretch, 1e-7 first-round, structural match non-negotiable) |
| `src/orchestration/pipeline.jl` | `run_njoy` entry point + module dispatch |
| `src/processing/reconr.jl` | RECONR top-level pipeline |
| `src/processing/broadr.jl` | Doppler broadening (broadn_grid + sigma1) |
| `src/processing/heatr.jl` | KERMA and damage energy |
| `src/processing/thermr.jl` | Thermal scattering (calcem + sigl + Bragg) |
| `src/processing/errorr.jl`, `src/orchestration/modules/errorr.jl` | Covariance processing (MF31/33, covcal-style σ·flx-weighted collapse) |
| `src/processing/pendf_writer.jl` | PENDF output |
| `test/validation/reference_test.jl` | Single-test Fortran-faithful runner |
| `test/validation/sweep_reference_tests.jl` | Full 84-test sweep writer |
| `njoy-reference/src/*.f90` | **Ground truth** — the canonical Fortran source |

## Development philosophy

1. **Fortran is ground truth.** Every formula, constant, rounding step,
   and IEEE-754 accumulation order matters. Read the specific subroutine
   in `njoy-reference/src/*.f90` before implementing or fixing Julia
   code. Cite line ranges in docstrings.
2. **Oracle-driven TDD.** Every non-trivial change starts with a failing
   reference-tape diff, then a fix, then a passing test. "Obvious" fixes
   are where oracles find surprises — every FP precision floor claimed in
   this project has turned out to be a real bug.
3. **Idiomatic Julia, Fortran semantic.** Multiple dispatch, parametric
   types, broadcasting — yes. Transliterated `do k=1,n` loops over
   workspace arrays — no. But sigfig bias, scratch-tape data paths,
   sentinel values, and non-associative accumulation order all ship
   verbatim.
4. **Tape-based architecture.** Modules communicate via tapes, never
   shared state. If module B needs data from module A, it reads A's
   output tape.

Full ruleset in [CLAUDE.md](CLAUDE.md). If you're a human or agent
picking this up, start there and then HANDOFF.md.

## License

GPL-3.0. See [LICENSE](LICENSE).

## References

1. R.E. MacFarlane et al., *The NJOY Nuclear Data Processing System,
   Version 2016*, LA-UR-17-20093, Los Alamos National Laboratory (2016).
2. A. Trkov et al. (Eds.), *ENDF-6 Formats Manual*, BNL-90365-2009 Rev. 2,
   Brookhaven National Laboratory (2018).
