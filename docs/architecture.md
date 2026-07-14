# Architecture Overview

NJOY.jl is a faithful, tape-driven Julia replacement for NJOY2016. Its
canonical integration surface is an NJOY input deck: `run_njoy` parses the
deck, dispatches module calls in order, and connects those modules through
numbered tape files. This mirrors the architecture of the Fortran program.

The repository also exposes typed, in-memory processing functions such as
`reconstruct`, `doppler_broaden`, and `group_average`. They are useful for
testing and direct library use, but they are a secondary API. They do not
replace tape boundaries between NJOY modules and are not, by themselves, the
drop-in NJOY execution path.

The Fortran source in `njoy-reference/src/` defines the required processing
semantics. Architectural convenience never overrides its record ordering,
rounding, accumulation order, sentinels, or scratch-tape behavior.

---

## Canonical deck and tape flow

```text
NJOY input deck
      |
      v
parse_njoy_input
      |
      v
NJOYInputDeck(Vector{ModuleCall})
      |
      v
run_njoy dispatch loop
      |
      +-- parse_reconr(call) --> reconr_module(TapeManager, params)
      |                              tape20 --> tape21
      |
      +-- parse_broadr(call) --> broadr_module(TapeManager, params)
      |                              tape21 --> tape22
      |
      +-- parse_heatr(call)  --> heatr_module(TapeManager, params)
      |                              tape22 --> tape23
      |
      +-- ... one wrapper per NJOY module ...
      |
      v
numbered output tapes (PENDF / ACE / GENDF / covariance / CCCC / plots)
```

`ModuleCall` preserves both tokenized cards and the original deck lines.
Modules with significant blank-card or delimiter semantics can therefore
parse the verbatim representation rather than guessing from normalized
tokens.

`TapeManager` maps Fortran logical-unit numbers to paths. The absolute unit
number selects a tape; the sign convention that distinguishes formatted and
unformatted Fortran units does not create a second Julia path. Unregistered
units resolve to `work_dir/tapeNN`.

The dispatcher is intentionally simple. Processing logic belongs in the
module wrapper and its processing helpers, not in a second orchestration
implementation. A normal wrapper:

1. resolves its input and output unit numbers through `TapeManager`;
2. reads all required information from its input tape or tapes;
3. invokes typed processing helpers;
4. writes a complete output tape; and
5. registers that output path when necessary.

The next module reads the written tape. Module order is therefore constrained
by the input deck and tape dependencies; NJOY-compatible module calls cannot
generally be reordered arbitrarily.

### In-memory work inside a tape boundary

File-backed module boundaries do not require every record to be processed one
line at a time. `PENDFTape`, `PENDFMaterial`, and `PENDFSection` provide a
structured in-memory representation used by several wrappers. A module may
read a complete input tape into these types, transform it, and write a new
tape. What matters architecturally is that downstream modules receive the
result through the output tape, not through a shared Julia object.

---

## Current transitional exceptions

The target architecture is strictly tape-based, but the current dispatcher
still has transitional assembly plumbing. These exceptions are implementation
debt, not patterns for new modules:

- `RunContext` accumulates RECONR, BROADR, HEATR, and THERMR data for a final
  PENDF assembly step.
- `run_njoy` retains the in-memory value returned by `broadr_module` and passes
  it, together with a separately computed RECONR result, to `heatr_module`.
  This bypasses the normal input-tape boundary to retain pre-formatting
  floating-point values.
- For chains needing combined thermal/heating output, `final_assembly!` can
  rebuild the last MODER output from accumulated context after the deck has
  finished.
- `thermr_module` writes an informational `.thermr` sidecar, while the current
  dispatcher recomputes MF6 data for final assembly instead of reading usable
  MF6 records from that sidecar.

New work must not expand these paths. The intended resolution is for each
module to write every required section to its own output tape and for the next
module to consume that tape, matching NJOY2016.

Some module options remain partial or stubbed. Their current maturity is
tracked in `HANDOFF.md` and the reference sweep, rather than encoded as an
architectural exception. Unsupported behavior should fail with enough context
to diagnose it unless the Fortran oracle demonstrates a different outcome.

---

## Core type hierarchy

### ENDF records

```text
InterpolationLaw
  Histogram, LinLin, LinLog, LogLin, LogLog, CoulombPen

MaterialId
InterpolationTable
ContRecord
ListRecord
Tab1Record
Tab2Record
TabulatedFunction
```

These types model ENDF-6 records and interpolation metadata. A
`TabulatedFunction` combines an interpolation table with its `x` and `y`
values and supports evaluation and integration.

### Resonance data

```text
AbstractResonanceFormalism
  SLBWParameters         (LRF=1)
  MLBWParameters         (LRF=2)
  ReichMooreParameters   (LRF=3)
  AdlerAdlerParameters   (LRF=4)
  SAMMYParameters        (LRF=7)
  UnresolvedParameters   (LRU=2)

ResonanceRange{P<:AbstractResonanceFormalism}
CrossSections{T<:Real}
IsotopeData
MF2Data
```

`ResonanceRange` is parameterized by formalism, allowing multiple dispatch to
select SLBW, MLBW, Reich–Moore, or SAMMY evaluation without a public integer
mode switch. The Julia type organization is idiomatic; each evaluator still
matches the corresponding Fortran algorithm and floating-point order.

### Processing and format data

```text
MF3Section, ENDFMaterial, PointwiseMaterial
AdaptiveConfig
KERMAResult, LindharParams, FissionQComponents
SABData, BraggData, ThermalResult
URRSpinSequence, URRStatModel, ProbabilityTable
MultiGroupXS
CovarianceBlock, CovarianceMatrix, CovarianceData

PENDFTape, PENDFMaterial, PENDFSection
TapeEntry, TapeDirectory, ENDFTapeSection, ENDFTapeMaterial

ACEHeader, ACETable, ACENeutronTable
ReactionXS, EquiprobableBins, TabulatedAngular, AngularBlock
```

`PointwiseMaterial` is the common low-level representation for pointwise
energies and reaction columns. It is not the object passed between module
wrappers during `run_njoy`; tapes are.

---

## Secondary processing-kernel API

The processing layer decomposes large NJOY algorithms into typed Julia
functions. These functions support focused tests and direct programmatic use.
They are also called by orchestration wrappers after the wrapper has read its
input tapes.

### RECONR kernels

```text
read_mf2 + read_mf3_sections
        |
        +--> build_evaluator(mf2)
        +--> build_grid(mf2, mf3_sections)
                    |
                    v
adaptive_reconstruct(evaluator, grid, config)
                    |
                    v
merge_background!(energies, values, mf3_sections, mf2)
                    |
                    v
PointwiseMaterial
```

`build_evaluator` returns a closure over parsed MF2 data. Formalism-specific
work is selected through dispatch. `adaptive_reconstruct` is a generic
linearization helper, while the legacy-compatible `reconr` path includes the
additional union-grid, redundant-reaction, sigfig, and writer behavior needed
for reference tapes.

### BROADR kernels

```text
sigma1_at
broadn_grid
doppler_broaden / doppler_broaden_multi
thin_xs
```

`sigma1_at` evaluates the Doppler-broadened cross section at one energy.
`broadn_grid` performs BROADR-specific adaptive refinement for one or more
reaction columns on a shared grid. BROADR does not call the RECONR
`adaptive_reconstruct` helper; the two algorithms have distinct Fortran
semantics.

### GROUPR kernels

```text
group_integrate(energies, values, bounds; law)
group_average(energies, cross_sections, mt_list, bounds; weight_fn)
```

The low-level integration API uses analytical ENDF panel integrals. The full
`groupr_module` wrapper additionally implements deck parsing, group/weight
selection, GENDF layout, reaction derivation, and other Fortran-specific
behavior.

### Direct API versus module equivalence

A direct call such as `doppler_broaden(material, 300.0)` is useful, but it does
not claim to reproduce every card, section-copying rule, header update, or
output-format quirk of a complete BROADR module invocation. Reference-test
equivalence is established through `run_njoy` and its tape outputs.

---

## Julia organization, Fortran semantics

NJOY2016 makes extensive use of module variables, workspace arrays, and
logical units. NJOY.jl replaces subroutine-local global state with explicit
arguments, closures, and typed structs where doing so does not change results.
Modules still communicate through tapes.

Multiple dispatch and parametric types organize the implementation, but they
do not make formalisms extensible in the product sense: supported behavior is
limited to NJOY2016 behavior exercised and verified against the oracle.

Mutation is used where it is natural or necessary to match the Fortran
algorithm and floating-point accumulation order. For example, adaptive grids,
workspace matrices, tape directories, and output buffers are mutated in
place. Automatic-differentiation compatibility is not a current architectural
guarantee, and `ForwardDiff` is not a package dependency.

---

## Source organization

```text
src/
  NJOY.jl          include order and public exports
  constants.jl     constants matching NJOY2016 phys.f90

  endf/            ENDF record types, readers, writers, interpolation
  resonances/      resonance readers and evaluated formalisms
  processing/      reusable processing algorithms and tape writers
  formats/         ACE, CCCC, MATXS, WIMS, DTF, and POWR formats
  orchestration/   input-deck parser, TapeManager, PENDF I/O, dispatcher
    modules/       one orchestration wrapper per NJOY module
  viewr/           Fortran-compatible plot-tape/PostScript engine
  visualization/   Julia-facing plot specifications and backends
```

The include order in `src/NJOY.jl` is the concrete dependency order. In broad
terms it loads ENDF and resonance foundations first, processing and format
implementations next, and orchestration wrappers after their parameter and
tape types are available.

For live module maturity, open work, and known traps, see `HANDOFF.md`. For
numeric and structural acceptance rules, see
`reports/ACCEPTANCE_CRITERIA.md`.
