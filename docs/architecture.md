# Architecture Overview

This document describes the design of NJOY.jl: how its types are organized,
how data flows through the processing chain, and how the architecture differs
from NJOY2016.

---

## Type hierarchy

### Core ENDF types

```
InterpolationLaw (enum)
  Histogram, LinLin, LinLog, LogLin, LogLog, CoulombPen

MaterialId
InterpolationTable
ContRecord
ListRecord
Tab1Record
Tab2Record
TabulatedFunction
```

These types map directly to the ENDF-6 format specification. Every ENDF
record is represented as an immutable struct with named fields. The
`TabulatedFunction` wraps a `Tab1Record` and provides evaluation and
integration methods via `interpolate` and `integrate`.

### Resonance parameter types

```
AbstractResonanceFormalism (abstract)
  |-- SLBWParameters         (LRF=1)
  |-- MLBWParameters         (LRF=2)
  |-- ReichMooreParameters   (LRF=3)
  |-- AdlerAdlerParameters   (LRF=4)
  |-- UnresolvedParameters   (LRU=2)

ResonanceRange{P<:AbstractResonanceFormalism}

CrossSections{T<:Real}

IsotopeData
MF2Data
```

The key design decision is parameterizing `ResonanceRange` on the formalism
type `P`. This enables Julia's multiple dispatch to select the correct cross
section evaluation function at compile time, producing type-stable code paths
without runtime branching on formalism flags.

`CrossSections` is parameterized on `T<:Real` so that ForwardDiff dual numbers
propagate through without truncation. Arithmetic operators (`+`, `*`) are
defined on `CrossSections` for accumulation.

### Processing types

```
MF3Section
ENDFMaterial
PointwiseMaterial

AdaptiveConfig

KERMAResult, LindharParams, FissionQComponents
SABData, BraggData, ThermalResult
URRSpinSequence, URRStatModel
ProbabilityTable
MultiGroupXS

CovarianceBlock, CovarianceMatrix, CovarianceData

TapeEntry, TapeDirectory, ENDFTapeSection, ENDFTapeMaterial
```

### ACE output types

```
ACEHeader
ACETable                   (flat NXS/JXS/XSS)
ACENeutronTable            (structured)
  ReactionXS
  EquiprobableBins
  TabulatedAngular
  AngularBlock
```

Two representations exist for ACE data: `ACENeutronTable` is a type-safe
structured form with named fields; `ACETable` is the flat serialization form
holding raw NXS/JXS/XSS arrays. Both can be written to the same Type 1
ASCII format via `write_ace`.

---

## Data flow

The standard processing pipeline:

```
ENDF file (text)
    |
    v
read_mf2(io) ---------> MF2Data
read_mf3_sections(io) -> Vector{MF3Section}
    |
    v
build_evaluator(mf2) --> closure: f(E) -> (sig_t, sig_e, sig_f, sig_c)
build_grid(mf2, mf3) --> initial energy grid
    |
    v
adaptive_reconstruct(f, grid, config) --> (energies, values)
    |
    v
merge_background!(energies, values, mf3, mf2)
    |
    v
PointwiseMaterial  <-- RECONR output
    |
    +---> doppler_broaden(pendf, T; awr, tol) --> PointwiseMaterial  [BROADR]
    |
    +---> compute_kerma(pendf; awr, Z, Q) ------> KERMAResult       [HEATR]
    |
    +---> compute_thermal(pendf, T, A; model) --> PointwiseMaterial  [THERMR]
    |
    +---> group_average(e, xs, mt, bounds) -----> MultiGroupXS       [GROUPR]
    |
    +---> build_ace(pendf; suffix, T) ----------> ACETable           [ACER]
    |       |
    |       v
    |     write_ace(io, table) -----> ACE file (text)
    |
    +---> write_pendf(io, pendf) ---> PENDF file (text)
```

Each processing step is a function that takes a `PointwiseMaterial` (or raw
arrays) and returns a new result. There is no shared mutable state between
steps. This makes it straightforward to:

1. Run steps in any order or skip steps entirely.
2. Process multiple materials in parallel.
3. Substitute custom implementations for any step.

### RECONR pipeline detail

The RECONR pipeline is decomposed into four composable functions:

```
build_evaluator : MF2Data -> (E -> NTuple{4,Float64})
build_grid      : (MF2Data, [MF3Section]) -> Vector{Float64}
adaptive_reconstruct : (f, grid, config) -> (Vector, Matrix)
merge_background! : (energies, values, mf3, mf2) -> nothing
```

The `build_evaluator` function returns a closure that captures the MF2
data and dispatches to `cross_section_slbw`, `cross_section_mlbw`, or
`cross_section_rm` based on the formalism type. The closure has signature
`f(E::Float64) -> NTuple{4, Float64}`.

The `adaptive_reconstruct` function is a generic higher-order function
that knows nothing about nuclear physics. It accepts any callable returning
an NTuple and refines the grid until the piecewise-linear representation
satisfies the tolerance. This same function is reused by BROADR.

### BROADR pipeline detail

```
sigma1_at      : (E, seg_e, seg_xs, alpha) -> Float64
doppler_broaden : (e, xs, T, awr; tol) -> (new_e, new_xs)
    internally calls: adaptive_reconstruct(broadened_eval, grid, config)
thin_xs        : (e, xs; tol) -> (thinned_e, thinned_xs)
```

The SIGMA1 kernel computes the exact Doppler-broadened cross section at a
single energy by integrating over the piecewise-linear input. The adaptive
grid algorithm (shared with RECONR) refines the output until tolerance is met.

### GROUPR pipeline detail

```
group_integrate : (energies, values, bounds; law) -> Vector{Float64}
group_average   : (e, xs, mt, bounds; weight_fn) -> MultiGroupXS
```

Group integration uses exact analytical panel integrals (the same
`panel_integral` function from the interpolation module, translating NJOY's
`gral`). This avoids quadrature error entirely for piecewise-linear data.

---

## How NJOY.jl differs from NJOY2016

### No global state

NJOY2016 relies heavily on Fortran COMMON blocks and module-level variables
for communication between subroutines. Cross section evaluation in RECONR
reads from global arrays (`res`, `a`, `b` in `reconr.f90`); the adaptive
grid writes results to pre-allocated global buffers.

NJOY.jl eliminates all global mutable state. The evaluator is a closure that
captures its data at construction time. The adaptive grid returns fresh arrays.
Processing parameters are passed as function arguments or in configuration
structs.

### No tape I/O

NJOY2016 uses Fortran sequential-access I/O with numbered "tapes" (logical
units). Data flows through temporary tapes: RECONR writes tape 21, BROADR
reads tape 21 and writes tape 22, etc.

NJOY.jl operates on in-memory data structures. `PointwiseMaterial` holds
the full pointwise representation. File I/O is limited to reading the input
ENDF file and writing the final output (PENDF or ACE). The `write_pendf` and
`write_ace` functions accept any `IO` stream.

### Composable functions instead of modules

In NJOY2016, each processing module (RECONR, BROADR, HEATR, etc.) is a
monolithic subroutine that reads its own card-image input deck, opens its own
tapes, and writes its own output. The modules are linked by tape numbers and
must run in a fixed sequence.

In NJOY.jl, each processing step is an ordinary function:

```julia
pendf = reconstruct("input.endf"; err=0.001)
pendf_T = doppler_broaden(pendf, 300.0; awr=236.0)
kerma = compute_kerma(pendf_T; awr=236.0, Z=92)
ace = build_ace(pendf_T; suffix="80c")
```

Functions can be called in any order, composed with `|>`, or wrapped in
higher-order functions. There is no card-image input syntax.

### Type dispatch instead of if/else chains

NJOY2016 selects the resonance formalism with integer flags and branches:

```fortran
if (lrf.eq.1) call csslbw(...)
if (lrf.eq.2) call csmlbw(...)
if (lrf.eq.3) call csrmat(...)
```

NJOY.jl uses Julia's parametric types and multiple dispatch:

```julia
cross_section(E, range::ResonanceRange{SLBWParameters}; ...) = ...
cross_section(E, range::ResonanceRange{MLBWParameters}; ...) = ...
cross_section(E, range::ResonanceRange{ReichMooreParameters}; ...) = ...
```

The compiler resolves the dispatch at compile time when the type is known,
eliminating runtime overhead. New formalisms can be added by defining a new
parameter type and a new method -- no modification to existing code.

### Generic adaptive grid

NJOY2016's adaptive grid logic (`resxs` in `reconr.f90`) is tightly coupled
to cross section evaluation. It directly calls `sigma` and knows about
resonance parameters, MF3 sections, and energy-range boundaries.

NJOY.jl's `adaptive_reconstruct` is a generic function parameterized by any
callable:

```julia
adaptive_reconstruct(f, initial_grid, config) -> (energies, values)
```

where `f` can be a resonance evaluator (RECONR), a Doppler-broadening kernel
(BROADR), a thermal scattering function (THERMR), or any user-defined
function. The algorithm operates solely on the function's return values and
the tolerance configuration, with no knowledge of the underlying physics.

---

## AD compatibility

NJOY.jl is designed for compatibility with Julia's automatic differentiation
ecosystem, particularly ForwardDiff.jl. The design constraints are:

1. **No mutation in hot paths.** Cross section evaluation functions return
   new values rather than writing to pre-allocated arrays. The
   `CrossSections{T}` struct propagates the element type `T` so that
   `ForwardDiff.Dual` numbers flow through unchanged.

2. **No try-catch in AD paths.** The Faddeeva function and all cross section
   kernels use conditional logic instead of exception handling.

3. **Type-parameterized containers.** `CrossSections{T}` and the SIGMA1
   kernel functions accept `Real` arguments and return type-stable results
   for any `T <: Real`.

4. **Pure functions.** Each processing function takes explicit inputs and
   returns explicit outputs. There are no side effects on global state that
   would confuse AD tracing.

This enables computing sensitivities of processed cross sections to resonance
parameters, temperatures, or other inputs using forward-mode AD:

```julia
using ForwardDiff

# Sensitivity of total XS to energy
dσ_dE = ForwardDiff.derivative(
    E -> cross_section(E, range).total,
    1.0e3  # energy in eV
)
```

The ERRORR module also provides a finite-difference `sensitivity_jacobian`
for cases where AD is not applicable.

---

## Source file organization

```
src/
  NJOY.jl                    Main module: includes + exports
  constants.jl               Physics constants (CODATA 2014)

  endf/
    types.jl                 ENDF record types and TabulatedFunction
    io.jl                    ENDF format readers and writers
    interpolation.jl         terp1, interpolate, integrate, panel_integral

  resonances/
    types.jl                 Formalism types and CrossSections
    penetrability.jl         P_l, S_l, phi_l
    faddeeva_exact.jl        Exact Faddeeva w(z), psi_chi
    faddeeva_table.jl        Lookup table + quickw
    reader.jl                MF2 reader (read_mf2)
    slbw.jl                  SLBW evaluation + cross_section dispatch
    mlbw.jl                  MLBW evaluation
    reich_moore.jl           Reich-Moore evaluation

  processing/
    adaptive_grid.jl         Generic adaptive linearization
    reconr_types.jl          MF3Section, ENDFMaterial, PointwiseMaterial
    reconr_evaluator.jl      build_evaluator, sigma_mf2, merge_background
    reconr_grid.jl           build_grid
    reconr.jl                reconstruct, reconr (top-level pipeline)
    pendf_writer.jl          write_pendf, write_pendf_file
    sigma1.jl                SIGMA1 kernel (f_func, h_func, sigma1_at)
    broadr.jl                doppler_broaden, thin_xs
    heatr.jl                 KERMA, Lindhard, heating functions
    thermr.jl                Free-gas, S(a,b), Bragg edges
    unresr.jl                Bondarenko self-shielding
    purr.jl                  Probability tables
    group_structures.jl      Built-in neutron group structures
    weight_functions.jl      Weight functions for GROUPR
    groupr.jl                group_integrate, group_average
    moder.jl                 Tape management
    errorr.jl                Covariance processing

  formats/
    ace_types.jl             ACEHeader, ACETable, NXS/JXS constants
    ace_neutron.jl           ACENeutronTable, ReactionXS, AngularBlock
    ace_builder.jl           build_ace, build_ace_from_pendf, build_xss
    ace_writer.jl            write_ace, write_ace_table, write_ace_directory
```

Each file has a single responsibility. The `include` order in `NJOY.jl`
reflects the dependency graph: ENDF types first, then resonance types, then
processing modules, then output formats.
