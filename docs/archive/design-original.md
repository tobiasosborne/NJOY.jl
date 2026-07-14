# NJOY.jl — Agent-Led Port of NJOY Nuclear Data Processing to Julia

## Mission

Port the NJOY2016 nuclear data processing system (~100k lines of Fortran 90 code, `njoy-reference/src/`) to idiomatic, composable, differentiable Julia. The result — `NJOY.jl` — must pass NJOY's own 87 reference test problems and produce bit-for-bit compatible PENDF/ACE output for the full ENDF/B-VIII.1 library. The port must not be a transliteration; it must be a clean reimagining that enables differentiable processing, on-the-fly reconstruction, uncertainty quantification, and compositional module development.

---

## Ground Truth: Local Repositories (MANDATORY)

**NJOY2016 Fortran source and tests MUST be cloned locally before any work begins.**
No agent may hallucinate, paraphrase, or guess at NJOY's implementation. Every design decision must be traced to the actual Fortran source.

```bash
# Clone NJOY2016 — the canonical reference
git clone https://github.com/njoy/NJOY2016.git ./njoy-reference

# Clone ENDFtk for format reference
git clone https://github.com/njoy/ENDFtk.git ./endftk-reference

# Download ENDF/B-VIII.0 test evaluations (at minimum: H-1, C-12, O-16, Fe-56, U-235, U-238, Pu-239)
# Use NNDC or IAEA mirror
mkdir -p ./endf-data
# Fetch key evaluations for testing
```

**Directory layout:**
```
project/
├── njoy-reference/          # NJOY2016 Fortran (READ-ONLY reference)
│   ├── src/                 # 39 Fortran source files
│   └── tests/               # 87 test problems with reference outputs
├── endftk-reference/        # Format reference
├── endf-data/               # Test ENDF evaluations
├── NJOY.jl/                 # Our Julia port (the deliverable)
│   ├── src/
│   ├── test/
│   ├── docs/
│   └── Project.toml
├── reports/                 # Agent decision reports
└── CLAUDE.md                # This file
```

**Compile and verify NJOY2016 builds and its tests pass locally before starting any porting work.** This is the ground truth oracle. If our Julia code disagrees with NJOY Fortran output, our code is wrong until proven otherwise.

---

## Agent Workflow: 3+1 Architecture

Every significant decision in this project uses a structured multi-agent workflow. This is non-negotiable.

### Phase 0: Deep Research Agent (RESEARCHER)

Before any porting begins, spawn a subagent with maximum thinking budget tasked with:

1. **Read every Fortran source file** in `njoy-reference/src/`. Produce a complete module dependency graph showing which modules call which subroutines, which module-level globals are shared, and which I/O units connect modules.

2. **Map every subroutine** in RECONR, BROADR, HEATR, THERMR, GROUPR, ACER, UNRESR/PURR, ERRORR to its physics purpose, its inputs/outputs, and its algorithmic approach. Produce a structured catalogue.

3. **Analyse all 87 test problems.** For each test, document: which modules are exercised, which ENDF materials are used, what the reference output checks, and what formalisms are tested (SLBW, MLBW, Reich-Moore, R-Matrix Limited, Adler-Adler, unresolved formats). Classify tests by which Julia module they validate.

4. **Study the ENDF-6 format specification.** Read the formats manual (available at https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf). Document every record type (CONT, LIST, TAB1, TAB2, INTG), every file number (MF), and every reaction type (MT) that NJOY processes.

5. **Study OpenMC's `openmc.data` Python module** for architectural inspiration. Note its type hierarchy (`IncidentNeutron`, `ResonanceRange`, `ReichMoore`, etc.). Identify what OpenMC got right and what it lacks (no BROADR, no THERMR, no full processing chain).

6. **Produce a RESEARCH_REPORT.md** with complete findings. This report is the shared knowledge base for all subsequent agents.

### Phase 1: Dual-Proposal Architecture (for every major decision)

For every module port and every architectural decision, spawn **two independent subagents** (PROPOSER-A and PROPOSER-B), each with maximum thinking budget. They must:

- Work independently without seeing each other's output
- Each produce a complete, working implementation or design
- Each write a rationale document explaining their choices
- Each identify risks and limitations of their approach

The orchestrator then:
1. Compares both proposals against the Fortran reference
2. Runs both through available tests
3. Selects the better approach OR synthesises a hybrid
4. Documents the decision with explicit reasoning

**This dual-proposal process applies to:**
- The overall Julia type hierarchy for ENDF data
- Each NJOY module port (RECONR, BROADR, HEATR, etc.)
- The ENDF-6 I/O layer design
- The adaptive grid refinement algorithm
- The Faddeeva/Doppler broadening implementation
- Cross-cutting concerns (error handling, logging, configuration)
- The test harness architecture

### Phase 2: Skeptical Reviewer Agent (REVIEWER)

After the orchestrator selects an approach, spawn a **rigorous skeptical reviewer subagent** with maximum thinking budget. The reviewer must:

1. **Diff against Fortran line-by-line.** For every physics formula, verify that the Julia implementation computes the same quantity. Watch for: sign conventions, factor-of-2 errors, energy unit conventions (eV vs MeV), mass unit conventions (amu vs neutron masses), channel radius formulas, penetrability/shift factor definitions.

2. **Check edge cases.** The Fortran code is full of special-case handling accumulated over decades (malformed ENDF files, zero-width resonances, negative cross sections, energy-dependent scattering radii, resolved-unresolved overlap regions). Verify the Julia port handles every one.

3. **Run numerical validation.** For the module under review, compute cross sections for at least 3 test nuclides (one light: C-12, one medium: Fe-56, one heavy: U-238) and compare pointwise against NJOY Fortran output. Differences must be < reconstruction tolerance (default 0.1%).

4. **Review for Julia idiom violations.** Flag: type instability, unnecessary allocations in hot loops, global mutable state, non-composable design, missing docstrings, un-exported public API surface.

5. **Produce a REVIEW_REPORT.md** with: PASS/FAIL verdict, list of issues found, list of issues fixed, residual risks, and sign-off.

**No code may be merged into the main `NJOY.jl` package without a PASS from the reviewer.**

---

## Testing Strategy

### Preferred: Run Julia Implementation Against Original NJOY Test Infrastructure

The strongly preferred approach is to make the Julia code produce output in ENDF/PENDF/ACE format that can be **byte-compared or numerically compared against NJOY's reference tapes.** This means:

1. For each of NJOY's 87 test problems, the Julia code should accept equivalent input parameters and produce output files in the same format as the `referenceTape*` files.

2. Write a Julia test harness that:
   - Reads the NJOY test `input` file to extract processing parameters
   - Runs the corresponding Julia processing chain
   - Writes output in PENDF/ACE/GENDF format
   - Compares against `referenceTape*` using numerical tolerance (not exact byte match, since floating-point representation may differ)
   - Reports per-energy-point maximum deviation for cross sections

3. For cross section comparisons, the tolerance is: `|σ_julia(E) - σ_njoy(E)| / σ_njoy(E) < err` where `err` is the reconstruction tolerance specified in the test input (typically 0.001 to 0.005).

### Fallback: Julia-Native Tests

Where format-level comparison is impractical (e.g., binary PENDF), write Julia-native tests that:

1. Read an ENDF evaluation directly
2. Reconstruct cross sections at specific energies
3. Compare against known values (from NJOY output, from EXFOR measurements, or from analytical formulas for simple cases)

### Unit Tests for Physics Functions

Every cross section formalism (SLBW, MLBW, Reich-Moore, Adler-Adler, hybrid R-function, unresolved) must have standalone unit tests verifying:
- Single-resonance analytical limits
- Multi-resonance interference patterns
- Threshold behaviour
- Doppler-broadened line shapes against analytical Voigt profiles
- Penetrability and shift factors against tabulated values
- The Faddeeva function against known high-precision implementations

### Continuous Integration

Structure tests into tiers:
- **Tier 1 (seconds):** Unit tests for physics functions, format I/O, type constructors
- **Tier 2 (minutes):** Single-material processing chains (H-1, Fe-56, U-238)
- **Tier 3 (hours):** Full 87-test-problem validation suite
- **Tier 4 (overnight):** Full ENDF/B-VIII.1 library processing

---

## Architecture Requirements

### Type System

Replace NJOY's flat `res(:)` array with a proper Julia type hierarchy. The RESEARCHER agent must study the Fortran pointer arithmetic patterns (`jnow`, `inow`, `ibaset`) and design typed alternatives. Minimum type hierarchy:

```
AbstractResonanceFormalism
├── SLBW                    # Single-Level Breit-Wigner (LRF=1)
├── MLBW                    # Multi-Level Breit-Wigner (LRF=2)
├── ReichMoore              # Reich-Moore (LRF=3)
├── AdlerAdler              # Adler-Adler (LRF=4)
├── RMatrixLimited          # R-Matrix Limited (LRF=7, via SAMMY)
├── HybridRFunction         # Hybrid R-function (LRF=6)
├── UnresolvedEnergyIndep   # Unresolved, energy-independent (LRU=2, LFW=0)
├── UnresolvedFissionWidth  # Unresolved, energy-dep fission (LRU=2, LFW=1)
└── UnresolvedEnergyDep     # Unresolved, energy-dependent (LRU=2, LRF=2)
```

Cross section evaluation must use **Julia multiple dispatch**: `cross_section(E, params::ReichMoore)` not a mode switch.

### No Global Mutable State

NJOY's ~80 module-level variables in RECONR alone must be eliminated. All state flows through function arguments. Processing functions must be pure where possible: `σ = cross_section(E, material, params)` with no side effects.

### No Scratch Tape I/O

NJOY's `loada`/`finda` buffered I/O via Fortran unit numbers is eliminated. All intermediate data lives in memory as Julia arrays and structs. The only I/O is reading the input ENDF file and writing the output PENDF/ACE file.

### StaticArrays for Hot Loops

The 3×3 complex matrix inversion in Reich-Moore (`frobns` subroutine) must use `SMatrix{3,3,ComplexF64}` from StaticArrays.jl. The `\` operator replaces 43 lines of manual Frobenius-Schur inversion. Penetrability/shift factor vectors should use `SVector`.

### Modern Faddeeva Function

Replace NJOY's 300 lines of hand-coded complex probability integral (`w`, `quickw`, `wtab` with 62×62 lookup table) with a proper implementation. Options: `SpecialFunctions.erfcx`, the `Faddeeva` package, or a clean standalone implementation. The 1970s lookup table optimisation is unnecessary on modern hardware.

### Composable Processing Chain

Each NJOY module becomes a Julia function that transforms data:

```julia
# The full processing chain is function composition
pendf = reconstruct(endf_material; tol=0.001)        # RECONR
pendf = doppler_broaden(pendf, T=600.0)               # BROADR  
pendf = compute_kerma(pendf, endf_material)            # HEATR
pendf = add_thermal_scattering(pendf, tsl)             # THERMR
gendf = group_average(pendf, flux, group_structure)    # GROUPR
ace  = format_ace(pendf)                               # ACER
```

Each function is independently usable. Unlike NJOY, you can call `reconstruct` without `doppler_broaden`, or call `cross_section(E, ...)` at a single energy without running the full chain.

### Differentiability (Design For, Don't Require Yet)

Structure the code so that `ForwardDiff.jl` or `Enzyme.jl` can differentiate through the cross section evaluation chain. This means: no mutation of arrays in the cross section inner loop, no try-catch in hot paths, no FFI calls in the physics core. The adaptive grid construction can be non-differentiable, but `cross_section(E, params) → σ` must be AD-compatible.

---

## Module Porting Order

Port in dependency order. Each module goes through the full 3+1 agent workflow.

### Wave 1: Foundation
1. **ENDF I/O** — Read/write ENDF-6 format (CONT, LIST, TAB1, TAB2 records). Parse the 80-column format. Handle the ENDF float format (11-char fields, sometimes no `E` in exponent like `1.234567+8`). Write valid PENDF output.
2. **Core Types** — Material, Isotope, ResonanceParameters, CrossSections, EnergyGrid, InterpolationTable.
3. **Physics Constants** — neutron mass, ℏ, Boltzmann constant, matching NJOY's values in `physics.f90` (use CODATA values but document any differences).

### Wave 2: RECONR (the critical module)
4. **Resonance Readers** — Parse File 2 for all formalisms (SLBW, MLBW, RM, AA, hybrid, unresolved ×3).
5. **Cross Section Evaluation** — Implement all formalisms. This is the physics core.
6. **Adaptive Grid Construction** — Stack-based bisection with dual tolerance criterion.
7. **Energy Grid Merging** — Union of resonance grid with File 3 background cross sections.
8. **PENDF Output** — Write reconstructed pointwise cross sections.

### Wave 3: Core Processing Chain
9. **BROADR** — Doppler broadening via Solbrig kernel or σ₁ method. Thinning.
10. **HEATR** — KERMA coefficients, energy-momentum conservation.
11. **THERMR** — Thermal scattering kernels, S(α,β).
12. **UNRESR/PURR** — Unresolved self-shielding (Bondarenko method / probability tables).

### Wave 4: Output Formatters
13. **GROUPR** — Multigroup averaging with self-shielding.
14. **ACER** — ACE format output for Monte Carlo codes.
15. **MODER** — Format conversion utility.
16. **ERRORR** — Covariance processing (lower priority but high value for differentiability story).

### Wave 5: Integration and Validation
17. **Full pipeline integration tests**
18. **87-test-problem validation**
19. **Performance benchmarking vs NJOY Fortran**
20. **Documentation and API reference**

---

## Coding Standards

- **File size limit:** No Julia source file exceeds 300 lines. Split by physics concept, not by legacy module boundaries.
- **Docstrings:** Every exported function has a docstring with: purpose, arguments, return value, physics reference (equation number in NJOY manual or ENDF formats manual), and at least one example.
- **Type stability:** Every function in the cross section evaluation hot path must be type-stable. Verify with `@code_warntype`.
- **Testing:** Every function has at least one test. Physics functions have property-based tests where applicable (e.g., σ_total = σ_elastic + σ_capture + σ_fission).
- **No magic numbers:** Every physical constant, numerical tolerance, and format parameter is a named constant with a comment citing its source.
- **Unicode where it helps:** Use `σ`, `Γ`, `ψ`, `χ`, `ℓ` in variable names for physics quantities. Use ASCII for all public API names.

---

## Decision Log

Every decision made through the 3+1 process must be logged in `reports/decisions/NNN-topic.md` with:
- Decision number and date
- The two proposals (summarised)
- The orchestrator's selection rationale
- The reviewer's verdict
- Any dissenting notes

---

## Success Criteria

The port is complete when:

1. `NJOY.jl` processes all 87 NJOY test problems and produces numerically equivalent output (within reconstruction tolerance) for every test that exercises a ported module.
2. `NJOY.jl` can process the full ENDF/B-VIII.0 neutron sublibrary (~560 materials) through the RECONR→BROADR→ACER chain without failure.
3. `cross_section(E, material)` works as a standalone function call without running a full processing chain.
4. `ForwardDiff.gradient(E -> cross_section(E, material).capture, E₀)` returns a finite, correct derivative.
5. Total Julia source is kept under ~30 k lines (vs ~100 k lines of Fortran code; current: ~24 k).
6. All tests pass in CI.
7. Documentation includes a tutorial processing U-238 from ENDF to ACE.

---

## Key References

- NJOY2016 source: `./njoy-reference/src/`
- NJOY2016 tests: `./njoy-reference/tests/`
- NJOY2016 manual: LA-UR-17-20093 (https://t2.lanl.gov/nis/publications/NJOY2012.pdf — the 2012 manual covers the same algorithms)
- ENDF-6 formats manual: BNL-90365-2009 Rev.2 (https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf)
- OpenMC nuclear data module: https://docs.openmc.org/en/stable/pythonapi/data.html
- ENDFtk: https://github.com/njoy/ENDFtk
- FRENDY manual: JAEA-Data/Code 2018-014

---

## Anti-Patterns (DO NOT)

- **DO NOT** transliterate Fortran to Julia line-by-line. Redesign from the physics up.
- **DO NOT** preserve NJOY's scratch tape I/O, COMMON block patterns, or flat-array-with-pointer-arithmetic data structures.
- **DO NOT** skip the dual-proposal process for "simple" modules. Every module gets two proposals.
- **DO NOT** hallucinate NJOY behaviour. When uncertain, read the Fortran. When still uncertain, run NJOY and examine the output.
- **DO NOT** implement formalisms you cannot test. If no test problem exercises Adler-Adler (LRF=4), note this and defer.
- **DO NOT** optimise prematurely. Correctness first, performance second. A correct Julia implementation will already be fast.
- **DO NOT** use global mutable state. Period.
- **DO NOT** merge any code without a reviewer PASS.
