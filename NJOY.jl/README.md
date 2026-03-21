# NJOY.jl

A Julia reimplementation of the [NJOY2016](https://github.com/njoy/NJOY2016)
nuclear data processing system.  NJOY.jl reads evaluated nuclear data in
ENDF-6 format, reconstructs resonance cross sections, Doppler-broadens them,
and produces pointwise (PENDF) and ACE-format libraries suitable for Monte
Carlo transport codes such as MCNP.

## Key features

- **Complete ENDF-6 I/O** -- reads and writes CONT, LIST, TAB1, TAB2 records
  with full interpolation law support (histogram through Coulomb penetrability).
- **Resonance reconstruction** -- Single-Level Breit-Wigner (SLBW),
  Multi-Level Breit-Wigner (MLBW), and Reich-Moore (RM) formalisms with
  energy-dependent scattering radii.
- **Doppler broadening** -- exact SIGMA1 kernel with adaptive grid refinement,
  matching NJOY2016's BROADR output.
- **Heating and damage** -- KERMA coefficients, Lindhard damage energy
  partition, and sum-rule verification (HEATR).
- **Thermal scattering** -- free-gas and S(alpha,beta) models, coherent
  elastic Bragg edges, incoherent elastic (THERMR).
- **Unresolved resonance self-shielding** -- Bondarenko method (UNRESR) and
  probability-table generation via Monte Carlo ladder sampling (PURR).
- **Multigroup averaging** -- exact piecewise panel integration with pluggable
  weight functions and Bondarenko self-shielding (GROUPR).
- **Covariance processing** -- MF33 reading, NI-type block expansion,
  multigroup collapse, sandwich rule, sensitivity Jacobians (ERRORR).
- **ACE output** -- Type 1 (ASCII) ACE files for MCNP with full NXS/JXS/XSS
  layout, matching NJOY2016 `aceout()` format exactly (ACER).
- **Tape management** -- directory scanning, material extraction, tape merging,
  and structured round-trip I/O (MODER).
- **Differentiable** -- all core functions are pure (no mutation, no global
  state), parameterized on element type, and compatible with ForwardDiff.jl
  automatic differentiation.
- **Composable** -- the processing chain is built from ordinary Julia
  functions and closures; no card-image input decks, no tape numbers.
- **Compact** -- approximately 9,100 lines of Julia replace roughly 120,000
  lines of Fortran (about 13:1 reduction).

## Installation

NJOY.jl is not yet registered. Install it directly from the repository:

```julia
using Pkg
Pkg.develop(path="/path/to/NJOY.jl/NJOY.jl")
```

Requires Julia 1.10 or later. Dependencies (installed automatically):

| Package            | Purpose                              |
|--------------------|--------------------------------------|
| SpecialFunctions   | `erfc`, `erfcx` for Faddeeva/Doppler |
| StaticArrays       | Fixed-size matrices in Reich-Moore    |
| LinearAlgebra      | Eigenvalue checks in ERRORR          |

## Quick start

```julia
using NJOY

# 1. Reconstruct pointwise cross sections from an ENDF file
pendf = reconstruct("n-092_U_238.endf"; err=0.001)

# 2. Doppler-broaden to 300 K
pendf_300K = doppler_broaden(pendf, 300.0; awr=236.0, tol=0.001)

# 3. Write a PENDF tape
open("u238.pendf", "w") do io
    write_pendf(io, pendf_300K)
end

# 4. Produce an ACE file for MCNP
ace = build_ace(pendf_300K; suffix="80c", temperature=300.0)
open("92238.80c", "w") do io
    write_ace(io, ace)
end
```

## Module overview

| Module   | Description                                                  |
|----------|--------------------------------------------------------------|
| ENDF I/O | Parse and write ENDF-6 CONT, LIST, TAB1, TAB2 records       |
| Interp   | ENDF interpolation laws (lin-lin through Coulomb) and integration |
| MF2      | Read File 2 resonance parameters (SLBW, MLBW, RM, URR)      |
| RECONR   | Adaptive resonance cross section reconstruction              |
| BROADR   | Doppler broadening via SIGMA1 kernel                         |
| HEATR    | KERMA coefficients, Lindhard damage energy                   |
| THERMR   | Free-gas and S(alpha,beta) thermal scattering                |
| UNRESR   | Bondarenko self-shielding for unresolved resonances          |
| PURR     | Probability table generation by Monte Carlo ladder sampling  |
| GROUPR   | Multigroup cross section averaging with pluggable weights    |
| MODER    | ENDF tape management (directory, extract, merge)             |
| ERRORR   | Covariance processing (MF33 expand, collapse, sandwich rule) |
| ACER     | ACE format output for MCNP (Type 1 ASCII)                   |

## Comparison with NJOY2016

| Aspect               | NJOY2016 (Fortran)           | NJOY.jl (Julia)                  |
|----------------------|------------------------------|----------------------------------|
| Lines of code        | ~120,000                     | ~9,100                           |
| Input format         | Card-image input deck        | Function keyword arguments       |
| Global state         | Extensive (common blocks)    | None (pure functions + closures) |
| Tape I/O             | Sequential-access unit numbers | Ordinary `IO` streams           |
| Type safety          | Implicit typing              | Parametric types + dispatch      |
| AD compatibility     | Not possible                 | ForwardDiff-compatible           |
| Interpolation        | Procedural (`terp1`)         | Dispatch on `InterpolationLaw`   |
| Resonance formalisms | `if/else` on LRF             | Multiple dispatch on formalism type |
| Adaptive grid        | Tightly coupled to RECONR    | Generic higher-order function    |

## Testing

Run the test suite from the package directory:

```julia
using Pkg
Pkg.test("NJOY")
```

Or from the command line:

```bash
cd NJOY.jl
julia --project -e 'using Pkg; Pkg.test()'
```

The test suite covers:

- **Physics constants** -- exact match against NJOY2016 CODATA 2014 values
- **ENDF I/O** -- float parsing/formatting round-trips, record read/write
- **Interpolation** -- all six ENDF interpolation laws, panel integrals
- **Penetrability** -- P_l, S_l, phi_l for l = 0..4 against analytic formulas
- **Faddeeva function** -- exact evaluation vs. SpecialFunctions.jl reference
- **Cross sections** -- SLBW, MLBW, Reich-Moore at selected energies
- **RECONR pipeline** -- full reconstruction of test ENDF files
- **BROADR** -- Doppler broadening with tolerance checks
- **HEATR** -- elastic/capture/fission heating, Lindhard damage, sum rule
- **THERMR** -- free-gas cross sections, Bragg edges, S(alpha,beta)
- **UNRESR/PURR** -- Bondarenko integrals, probability tables
- **GROUPR** -- group integration, weighted averaging, self-shielding
- **MODER** -- tape directory, material extraction, merge
- **ERRORR** -- covariance block expansion, multigroup collapse, symmetry/PSD
- **ACER** -- ACE table construction, ZAID formatting, Type 1 output

## Contributing

1. Fork the repository and create a feature branch.
2. Write tests for any new functionality in `test/runtests.jl`.
3. Ensure `Pkg.test("NJOY")` passes before submitting a pull request.
4. Follow the existing code style: one function per concept, docstrings on
   all exported functions, no global mutable state.

## License

MIT License. See the LICENSE file for details.

## References

1. R.E. MacFarlane, D.W. Muir, R.M. Boicourt, A.C. Kahler III, *The NJOY
   Nuclear Data Processing System, Version 2016*, LA-UR-17-20093, Los Alamos
   National Laboratory (2016).
2. A. Trkov, M. Herman, D.A. Brown (Eds.), *ENDF-6 Formats Manual*,
   BNL-90365-2009 Rev. 2, Brookhaven National Laboratory (2018).
3. M.B. Chadwick et al., *ENDF/B-VIII.0: The 8th Major Release of the Nuclear
   Reaction Data Library*, Nuclear Data Sheets **148**, 189-240 (2018).
