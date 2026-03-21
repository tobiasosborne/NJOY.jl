# API Reference

Complete reference for all exported functions and types in NJOY.jl, organized
by module.

---

## Physics Constants

```julia
NJOY.PhysicsConstants
NJOY.CODATA2014
```

Module containing fundamental physics constants matching NJOY2016's `phys.f90`
(CODATA 2014 / ENDF-102 Appendix H).

Key constants:

| Constant  | Value                | Description                    |
|-----------|----------------------|--------------------------------|
| `bk`      | 8.617333262e-5       | Boltzmann constant [eV/K]      |
| `ev`      | 1.602176634e-12      | electron-volt [erg/eV]         |
| `clight`  | 2.99792458e10        | speed of light [cm/s]          |
| `amassn`  | 1.00866491595        | neutron mass [amu]             |
| `amassp`  | 1.007276466621       | proton mass [amu]              |
| `amassa`  | 4.001506179127       | alpha particle mass [amu]      |
| `amasse`  | 5.48579909065e-4     | electron mass [amu]            |
| `pi`      | 3.141592653589793238 | pi                             |
| `finstri` | ~137.036             | inverse fine structure constant |

---

## ENDF I/O

### Types

#### `MaterialId`
```julia
struct MaterialId
    mat::Int32   # material number (MAT)
    mf::Int32    # file number (MF)
    mt::Int32    # section number (MT)
end
```
ENDF material/file/section identifier triple.

#### `InterpolationLaw`
```julia
@enum InterpolationLaw::Int32 begin
    Histogram   = 1
    LinLin      = 2
    LinLog      = 3
    LogLin      = 4
    LogLog      = 5
    CoulombPen  = 6
end
```
ENDF interpolation law codes.

#### `InterpolationTable`
```julia
InterpolationTable(nbt::Vector, law::Vector)
```
Breakpoint indices and interpolation laws for a tabulated function.

#### `ContRecord`
```julia
struct ContRecord
    C1::Float64; C2::Float64
    L1::Int32; L2::Int32; N1::Int32; N2::Int32
    id::MaterialId
end
```
ENDF CONT (control) record: 2 floats + 4 integers.

#### `ListRecord`
```julia
struct ListRecord
    C1::Float64; C2::Float64
    L1::Int32; L2::Int32; N1::Int32; N2::Int32
    data::Vector{Float64}
    id::MaterialId
end
```
ENDF LIST record: CONT header followed by N1 data values.

#### `Tab1Record`
```julia
struct Tab1Record
    C1::Float64; C2::Float64; L1::Int32; L2::Int32
    interp::InterpolationTable
    x::Vector{Float64}; y::Vector{Float64}
    id::MaterialId
end
```
ENDF TAB1 record: CONT header + interpolation table + (x,y) pairs.

#### `Tab2Record`
```julia
struct Tab2Record
    C1::Float64; C2::Float64; L1::Int32; L2::Int32
    interp::InterpolationTable; NZ::Int32
    id::MaterialId
end
```
ENDF TAB2 record: CONT header + interpolation table (header for TAB1 sequence).

#### `TabulatedFunction`
```julia
TabulatedFunction(interp, x, y)
TabulatedFunction(rec::Tab1Record)
```
A tabulated function with interpolation metadata, suitable for evaluation and
integration.

### Parsing functions

#### `parse_endf_float`
```julia
parse_endf_float(s::AbstractString) -> Float64
```
Parse an 11-character ENDF float field. Handles compact notation without `E`
(e.g., `" 1.234567+8"` -> `1.234567e8`).

#### `format_endf_float`
```julia
format_endf_float(x::Real) -> String
```
Format a number as an 11-character ENDF float field matching NJOY's `a11`.

### Record readers

#### `read_tpid`
```julia
read_tpid(io::IO) -> NamedTuple{(:text, :mat, :mf, :mt)}
```
Read the TPID (tape identification) record -- the first line of an ENDF tape.

#### `find_section`
```julia
find_section(io::IO, target_mf::Integer, target_mt::Integer) -> Bool
```
Scan from the beginning of `io` until a line with the given MF/MT is found.
Positions the stream at the start of that line. Returns `true` if found.

#### `read_cont`
```julia
read_cont(io::IO) -> ContRecord
```
Read one CONT record (single 80-column line).

#### `read_head`
```julia
read_head(io::IO) -> ContRecord
```
Read a HEAD record (same format as CONT).

#### `read_list`
```julia
read_list(io::IO) -> ListRecord
```
Read a LIST record: CONT header + ceil(NPL/6) data lines.

#### `read_tab1`
```julia
read_tab1(io::IO) -> Tab1Record
```
Read a TAB1 record: CONT header + interpolation table + (x,y) data.

#### `read_tab2`
```julia
read_tab2(io::IO) -> Tab2Record
```
Read a TAB2 record: CONT header + interpolation table.

### Record writers

#### `write_cont`
```julia
write_cont(io::IO, rec::ContRecord; ns::Int=0) -> Int
```
Write a CONT record. Returns the next sequence number.

#### `write_list`
```julia
write_list(io::IO, rec::ListRecord; ns::Int=0) -> Int
```
Write a LIST record.

#### `write_tab1`
```julia
write_tab1(io::IO, rec::Tab1Record; ns::Int=0) -> Int
```
Write a TAB1 record.

#### `write_tab2`
```julia
write_tab2(io::IO, rec::Tab2Record; ns::Int=0) -> Int
```
Write a TAB2 record.

---

## Interpolation

#### `terp1`
```julia
terp1(x1, y1, x2, y2, x, law; coulomb_threshold=0.0) -> Float64
```
Interpolate one point between (x1,y1) and (x2,y2) using ENDF interpolation
`law`. Translation of NJOY2016's `terp1`.

#### `interpolate`
```julia
interpolate(tab::TabulatedFunction, x; coulomb_threshold=0.0) -> Float64
```
Interpolate a tabulated function at `x`. Returns 0.0 outside the data range.

#### `integrate`
```julia
integrate(tab::TabulatedFunction, x1, x2; coulomb_threshold=0.0) -> Float64
```
Integrate the tabulated function from `x1` to `x2` using analytical panel
integrals matching NJOY's `intega`/`gral`.

---

## Resonance Evaluation

### Types

#### `AbstractResonanceFormalism`
Abstract supertype for all resonance parameter formalisms.

#### `CrossSections{T<:Real}`
```julia
struct CrossSections{T<:Real}
    total::T; elastic::T; fission::T; capture::T
end
```
Result container for cross sections at a single energy (barns). Parameterized
on `T` for AD compatibility.

#### `SLBWParameters <: AbstractResonanceFormalism`
Single-Level Breit-Wigner resonance parameters (LRF=1). Contains per-l-group
arrays of resonance energies, spins, and partial widths (Gn, Gg, Gf, Gx).

#### `MLBWParameters <: AbstractResonanceFormalism`
Multi-Level Breit-Wigner parameters (LRF=2). Same storage as SLBW; uses
multi-level interference formula for elastic scattering.

#### `ReichMooreParameters <: AbstractResonanceFormalism`
Reich-Moore parameters (LRF=3). 3-channel R-matrix with two fission widths
(Gfa, Gfb).

#### `AdlerAdlerParameters <: AbstractResonanceFormalism`
Adler-Adler parameters (LRF=4). Background polynomial + resonance-term
parameterization.

#### `UnresolvedParameters <: AbstractResonanceFormalism`
Unresolved resonance region parameters (LRU=2). Three variants:
`:energy_independent`, `:fission_width`, `:energy_dependent`.

#### `ResonanceRange{P}`
```julia
struct ResonanceRange{P<:AbstractResonanceFormalism}
    EL::Float64; EH::Float64
    LRU::Int32; LRF::Int32; LFW::Int32; NRO::Int32; NAPS::Int32
    parameters::P
    ap_tab::Union{Nothing,TabulatedFunction}
end
```
Container for one resonance energy range, parameterized on formalism type.

#### `MF2Data`
```julia
struct MF2Data
    ZA::Float64; AWR::Float64
    isotopes::Vector{IsotopeData}
end
```
All File 2 resonance data for a material.

#### `IsotopeData`
```julia
struct IsotopeData
    ZAI::Float64; ABN::Float64; LFW::Int32
    ranges::Vector{ResonanceRange}
end
```
Resonance data for one isotope.

### Functions

#### `read_mf2`
```julia
read_mf2(io::IO) -> MF2Data
```
Read the entire MF2/MT151 section. The stream must be positioned at the HEAD
record (after `find_section(io, 2, 151)`).

#### `cross_section`
```julia
cross_section(E, range::ResonanceRange; temperature=0.0, table=nothing) -> CrossSections
```
Compute resonance cross sections at energy E [eV]. Dispatches on the
parameter type (SLBW, MLBW, or RM) for type stability.

#### `cross_section_slbw`
```julia
cross_section_slbw(E, params::SLBWParameters, range; temperature=0.0, table=nothing) -> CrossSections
```
Single-Level Breit-Wigner evaluation. Supports T=0 (analytic) and T>0
(Doppler-broadened via Faddeeva function).

#### `cross_section_mlbw`
```julia
cross_section_mlbw(E, params::MLBWParameters, range) -> CrossSections
```
Multi-Level Breit-Wigner evaluation (T=0). Includes inter-resonance
interference for the elastic channel.

#### `cross_section_rm`
```julia
cross_section_rm(E, params::ReichMooreParameters, range) -> CrossSections
```
Reich-Moore evaluation (T=0). Uses 3x3 R-matrix inversion for fissile
materials, scalar R-function for non-fissile.

#### `cwaven_constant`
```julia
cwaven_constant() -> Float64
```
Wavenumber constant: `sqrt(2 * m_n * amu * eV) * 1e-12 / hbar`.

#### `channel_radius`
```julia
channel_radius(awri::Float64) -> Float64
```
Hard-sphere channel radius: `a = 0.123 * (m_n * AWRI)^(1/3) + 0.08` [fm].

### Penetrability functions

#### `penetrability`
```julia
penetrability(l::Integer, rho::Real) -> Float64
```
Penetrability factor P_l(rho) for orbital angular momentum l. Exact formulas
for l = 0..4; upward recursion for l >= 5.

#### `shift_factor`
```julia
shift_factor(l::Integer, rho::Real) -> Float64
```
Shift factor S_l(rho). Exact for l = 0..4.

#### `phase_shift`
```julia
phase_shift(l::Integer, rho::Real) -> Float64
```
Hard-sphere phase shift phi_l(rho). Exact for l = 0..4.

### Faddeeva function

#### `faddeeva_w`
```julia
faddeeva_w(rez, aim1) -> (rew, aimw)
```
Complex probability integral w(z) = exp(-z^2) * erfc(-iz). Direct translation
of NJOY2016's `w` subroutine. AD-compatible.

#### `faddeeva_w_julia`
```julia
faddeeva_w_julia(x, y) -> (rew, aimw)
```
Reference implementation using `SpecialFunctions.erfcx`.

#### `FaddeevaTable`
Precomputed 62x62 lookup table for fast evaluation.

#### `build_faddeeva_table`
```julia
build_faddeeva_table() -> FaddeevaTable
```
Construct the lookup table (translation of NJOY's `wtab`).

#### `quickw`
```julia
quickw(ax, y, table::FaddeevaTable; compute_imaginary=true) -> (rew, aimw)
```
Fast evaluation using table interpolation + rational approximations.
AD-compatible.

#### `psi_chi`
```julia
psi_chi(x, y) -> (psi, chi)
psi_chi(x, y, table::FaddeevaTable) -> (psi, chi)
```
Voigt line profiles via the Faddeeva function.

---

## RECONR

### Types

#### `MF3Section`
```julia
struct MF3Section
    mt::Int32; QM::Float64; QI::Float64
    tab::TabulatedFunction
end
```
A single MF3 cross section section from ENDF.

#### `ENDFMaterial`
```julia
struct ENDFMaterial
    mat::Int32; awr::Float64
    mf2::MF2Data; mf3_sections::Vector{MF3Section}
    redundant_mts::Set{Int}
end
```
All ENDF data needed for RECONR processing.

#### `PointwiseMaterial`
```julia
struct PointwiseMaterial
    mat::Int32
    energies::Vector{Float64}
    cross_sections::Matrix{Float64}   # (n_energies, n_reactions)
    mt_list::Vector{Int}
end
```
Result of RECONR: pointwise cross sections on a linearized grid.

#### `AdaptiveConfig`
```julia
AdaptiveConfig(err; errmax=10*err, errint=err/20000, max_depth=30,
               thermal_threshold=0.4999, thermal_factor=5.0,
               step_guard_limit=Inf)
```
Tolerance parameters for adaptive grid refinement.

### Functions

#### `reconstruct`
```julia
reconstruct(endf_file; mat=0, err=0.001, errmax=nothing,
            errint=nothing, temperature=0.0) -> PointwiseMaterial
```
Full RECONR pipeline: read ENDF, reconstruct resonances, merge backgrounds.

- `endf_file`: path to ENDF-format file
- `mat`: MAT number (0 = auto-detect)
- `err`: fractional reconstruction tolerance
- `temperature`: [K] (0.0 = zero temperature)

#### `reconr`
```julia
reconr(endf_file; mat=0, err=0.001, ...) -> NamedTuple
```
Legacy interface returning a NamedTuple with fields `energies`, `total`,
`elastic`, `fission`, `capture`, `mf2`, `mf3_sections`.

#### `build_evaluator`
```julia
build_evaluator(mf2::MF2Data; temperature=0.0, table=nothing) -> Function
```
Build a closure `f(E) -> NTuple{4, Float64}` that evaluates resonance cross
sections at energy E. Sums over isotopes weighted by abundance.

#### `build_grid`
```julia
build_grid(mf2::MF2Data, mf3_sections::Vector{MF3Section}) -> Vector{Float64}
```
Build the initial energy grid from MF2 resonance nodes + MF3 breakpoints.

#### `adaptive_reconstruct`
```julia
adaptive_reconstruct(f, grid, config::AdaptiveConfig) -> (energies, values)
```
Generic adaptive linearization. `f` is any callable returning an NTuple.

#### `merge_background!`
```julia
merge_background!(energies, values, mf3_sections, mf2)
```
Add MF3 background cross sections to reconstructed resonance data in-place.

#### `read_mf3_sections`
```julia
read_mf3_sections(io::IO, mat::Integer) -> Vector{MF3Section}
```
Read all MF3 sections for a material.

#### `sigma_mf2`
```julia
sigma_mf2(E, mf2::MF2Data) -> CrossSections
```
Evaluate resonance cross sections at E using all ranges in MF2.

#### `round_sigfig`
```julia
round_sigfig(x, n) -> Float64
```
Round to `n` significant figures (translation of NJOY's `sigfig`).

### PENDF writer

#### `write_pendf`
```julia
write_pendf(io::IO, material::PointwiseMaterial;
            endf_source=nothing, temperature=0.0, err=0.001,
            description=String[])
```
Write a PENDF tape in standard ENDF-6 format.

#### `write_pendf_file`
```julia
write_pendf_file(filename, material; ...)
```
File-path convenience wrapper for `write_pendf`.

---

## BROADR

#### `doppler_broaden` (vector form)
```julia
doppler_broaden(energies, xs, T, awr; tol=0.001, do_thin=true, ...) -> (new_e, new_xs)
```
Doppler-broaden a single piecewise-linear cross section from 0 K to T.

#### `doppler_broaden_multi`
```julia
doppler_broaden_multi(energies, xs_matrix, T, awr; ...) -> (new_e, new_xs)
```
Broaden multiple cross sections simultaneously on a shared adaptive grid.

#### `doppler_broaden` (PointwiseMaterial form)
```julia
doppler_broaden(pendf::PointwiseMaterial, T_new; T_old=0.0, awr=1.0,
                tol=0.001, thnmax=Inf) -> PointwiseMaterial
```
Broaden a `PointwiseMaterial` from `T_old` to `T_new`.

#### `thin_xs`
```julia
thin_xs(energies, xs; tol=0.001, step_max=1.24) -> (thinned_e, thinned_xs)
```
Remove redundant points where linear interpolation suffices within `tol`.
Matches NJOY2016's `thinb`.

### SIGMA1 kernel functions

#### `sigma1_at`
```julia
sigma1_at(E, seg_e, seg_xs, alpha) -> Float64
```
Evaluate the Doppler-broadened cross section at energy E using the exact
SIGMA1 integration kernel.

#### `f_func`
```julia
f_func(n::Int, a) -> Real
```
F-function: `f_n(a) = integral(a..inf) z^n exp(-z^2) dz / sqrt(pi)` for
n = 0..4.

#### `f_all`
```julia
f_all(a) -> NTuple{5}
```
Return `(f_0, f_1, f_2, f_3, f_4)` in a single pass (one `exp` + one `erfc`).

#### `h_func`
```julia
h_func(n::Int, a, b) -> Real
```
H-function: `h_n(a,b) = f_n(a) - f_n(b)`. Uses Taylor series when
cancellation is detected.

#### `h_all`
```julia
h_all(a, b) -> NTuple{5}
```
Return `(h_0..h_4)` in a single pass.

#### `h_taylor`
```julia
h_taylor(n::Int, a, b) -> Real
```
Taylor-series expansion of h_n for small `|b-a|`.

---

## HEATR

### Types

#### `KERMAResult`
```julia
struct KERMAResult
    energies::Vector{Float64}
    total_kerma::Vector{Float64}       # MT301
    elastic_kerma::Vector{Float64}     # MT302
    capture_kerma::Vector{Float64}     # MT402
    fission_kerma::Vector{Float64}     # MT318
    inelastic_kerma::Vector{Float64}   # MT304
    damage_energy::Vector{Float64}     # MT444
end
```

#### `LindharParams`
```julia
struct LindharParams
    el_inv::Float64   # 1/E_L [1/eV]
    fl::Float64       # electronic stopping factor
end
```

#### `FissionQComponents`
```julia
struct FissionQComponents
    E_fragments::Float64; E_prompt_n::Float64; E_prompt_gamma::Float64
    E_delayed_beta::Float64; E_delayed_gamma::Float64
    E_delayed_n::Float64; E_neutrino::Float64
    E_pseudo_Q::Float64; E_total::Float64
end
```

### Functions

#### `compute_kerma`
```julia
compute_kerma(pendf; awr=1.0, Z=0, Q_values=Dict(), E_d=nothing,
              local_gamma=false) -> KERMAResult
```
Compute KERMA factors for all reactions in a `PointwiseMaterial`.

#### `verify_kerma_sum_rule`
```julia
verify_kerma_sum_rule(result::KERMAResult; rtol=0.01) -> Bool
```
Check that total KERMA equals the sum of partials within `rtol`.

#### `elastic_heating`
```julia
elastic_heating(E, A) -> Float64
```
Average energy deposited per elastic scatter: `h = 2EA/(A+1)^2`.

#### `elastic_heating_aniso`
```julia
elastic_heating_aniso(E, A, mu_bar) -> Float64
```
Elastic heating with CM anisotropy correction.

#### `elastic_damage`
```julia
elastic_damage(E, A, Z; E_d=25.0) -> Float64
```
Damage energy from elastic scattering via Lindhard partition.

#### `capture_heating`
```julia
capture_heating(E, Q, A; E_gamma=0.0) -> Float64
```
Capture heating: `h = E*A/(A+1) + Q - E_gamma`.

#### `capture_recoil`
```julia
capture_recoil(E, Q, A) -> Float64
```
Compound nucleus recoil energy after capture.

#### `capture_damage`
```julia
capture_damage(E, Q, A, Z; E_d=25.0) -> Float64
```
Damage energy from radiative capture.

#### `fission_heating`
```julia
fission_heating(E, Q_total; E_neutrino=0.0, E_delayed=0.0) -> Float64
fission_heating(E, qc::FissionQComponents) -> Float64
```
Fission heating from energy balance or MT458 components.

#### `inelastic_heating`
```julia
inelastic_heating(E, Q, A; E_secondary=0.0) -> Float64
```

#### `nxn_heating`
```julia
nxn_heating(E, Q, A, n_out; E_secondary=0.0) -> Float64
```

#### `displacement_energy`
```julia
displacement_energy(Z::Integer) -> Float64
```
Displacement threshold energy [eV] for element Z (NJOY2016 defaults).

#### `lindhard_params`
```julia
lindhard_params(Zr, Ar, Zl, Al) -> LindharParams
```
Precompute Lindhard constants for recoil atom (Zr, Ar) in lattice (Zl, Al).

#### `lindhard_damage`
```julia
lindhard_damage(E_recoil, params::LindharParams; E_d=0.0) -> Float64
lindhard_damage(E_recoil, Zr, Ar, Zl, Al; E_d=0.0) -> Float64
```
Lindhard damage energy for a given recoil energy.

---

## THERMR

### Types

#### `SABData`
Tabulated S(alpha,beta) thermal scattering law (ENDF MF7/MT4).

#### `BraggData`
Precomputed Bragg edge data for coherent elastic scattering.

#### `ThermalResult`
Output of `compute_thermal`: thermal cross sections on an energy grid.

### Functions

#### `compute_thermal`
```julia
compute_thermal(pendf, T, A; model=:free_gas, emax=10.0,
                sigma_b=nothing, sab_data=nothing, bragg=nothing,
                debye_waller_prime=nothing) -> PointwiseMaterial
```
Compute thermal cross sections and update the PENDF material.

#### `compute_thermal_xs`
```julia
compute_thermal_xs(energies, sigma_elastic, A, T; sigma_b=nothing,
                   emax=10.0, model=:free_gas, sab_data=nothing) -> (e, xs)
```
Compute thermal elastic cross sections on a merged grid.

#### `free_gas_xs`
```julia
free_gas_xs(E, A, T; sigma_b=1.0) -> Float64
```
Analytical free-gas total thermal cross section [barn].

#### `free_gas_kernel`
```julia
free_gas_kernel(E, E_prime, mu, A, T; sigma_b=1.0) -> Float64
```
Free-gas double-differential kernel d^2sigma/(dE' dmu) [barn/eV].

#### `read_thermal_data`
```julia
read_thermal_data(alpha, beta, sab_values; sigma_b, awr, ...) -> SABData
```
Construct SABData from raw alpha/beta/S(a,b) arrays.

#### `sab_kernel`
```julia
sab_kernel(E, E_prime, mu, data::SABData, T) -> Float64
```
Differential kernel from tabulated S(alpha,beta) [barn/eV].

#### `sab_xs`
```julia
sab_xs(E, data::SABData, T; n_mu=20, n_ep=40) -> Float64
```
Integrated incoherent inelastic cross section from S(alpha,beta) [barn].

#### `bragg_edges`
```julia
bragg_edges(E, bragg::BraggData) -> Float64
```
Coherent elastic cross section from Bragg edges [barn].

#### `bragg_edge_energies`
```julia
bragg_edge_energies(b::BraggData) -> Vector{Float64}
```
Sorted Bragg edge energies [eV].

#### `build_bragg_data`
```julia
build_bragg_data(; a, c, sigma_coh, A_mass, natom=1,
                   debye_waller, emax=5.0, lat=1) -> BraggData
```
Build Bragg data for lattice type `lat`. a, c in [cm], sigma_coh [barn].
Lattice types: 1=graphite, 2=beryllium, 3=beryllium oxide.

#### `structure_factor`
```julia
structure_factor(lat, l1, l2, l3) -> Float64
```
Crystal structure factor for lattice type `lat`.

#### `incoh_elastic_xs`
```julia
incoh_elastic_xs(E, sigma_b, dwp) -> Float64
```
Incoherent elastic cross section: `sigma_b/2 * (1 - exp(-4*E*W'))`.

#### `THERMR_EGRID`
Standard thermal energy grid (116 points from 1e-5 to 10 eV).

---

## UNRESR

### Types

#### `URRSpinSequence`
Average resonance parameters for one (l, J) spin sequence.

#### `URRStatModel`
All spin sequences plus global parameters for one nuclide.

### Functions

#### `bondarenko_xs`
```julia
bondarenko_xs(model, E, T, sigma0_values; bkg=(0,0,0,0)) -> Matrix
```
Self-shielded Bondarenko cross sections. Returns a (5, nsigz) matrix with
rows: total, elastic, fission, capture, transport.

#### `infinite_dilution_xs`
```julia
infinite_dilution_xs(model, E, T) -> CrossSections
```
Infinitely-dilute average cross sections in the URR.

#### `urr_penetrability`
```julia
urr_penetrability(l, rho, rhoc) -> (Vl, phi)
```
Penetrability factor and phase shift for l = 0, 1, 2 (NJOY `uunfac`).

#### `ajku`
```julia
ajku(beta, sti) -> (xj, xk)
```
J and K integral contributions for the Bondarenko kernel.

#### `HWANG_QW`, `HWANG_QP`
10-point Gauss-Laguerre-type quadrature tables (weights and abscissae) for
chi-squared-distributed widths with 1-4 degrees of freedom.

---

## PURR

### Types

#### `ProbabilityTable`
Per-energy probability table with bins of (probability, total, elastic,
fission, capture) cross sections.

### Functions

#### `generate_ptable`
```julia
generate_ptable(model, energies; nladders=64, nbins=20, seed=12345,
                T=300.0, bkg=(0,0,0,0)) -> ProbabilityTable
```
Generate probability tables by Monte Carlo ladder sampling.

#### `generate_ladder`
```julia
generate_ladder(seq, elow, ehigh, rng) -> (er, gnr, gfr, ggr, gxr)
```
Generate one resonance ladder for a spin sequence.

#### `bondarenko_from_ptable`
```julia
bondarenko_from_ptable(ptable, ie, sigma0) -> NTuple{5,Float64}
```
Self-shielded cross sections from a probability table at energy index `ie`.
Returns (total, elastic, fission, capture, transport).

#### `chi2_sample`
```julia
chi2_sample(df, rng) -> Float64
```
Sample from chi-squared(df) using the NJOY 20-point quantile table.

#### `wigner_spacing`
```julia
wigner_spacing(D_mean, rng) -> Float64
```
Sample a level spacing from the Wigner distribution.

#### `CHI2_QUANTILES`
20x4 quantile table for chi-squared distributions with df = 1..4.

---

## GROUPR

### Types

#### `MultiGroupXS`
```julia
struct MultiGroupXS
    group_bounds::Vector{Float64}
    mt_list::Vector{Int}
    xs::Matrix{Float64}        # (n_groups, n_reactions)
    flux::Vector{Float64}      # group fluxes
end
```

### Group structures

Built-in standard neutron group structures (energy bounds in eV, ascending):

| Constant       | Identifier      | Groups | Description                 |
|----------------|-----------------|--------|-----------------------------|
| `LANL_30`      | `IGN_LANL30`    | 30     | LANL 30-group               |
| `WIMS_69`      | `IGN_WIMS69`    | 69     | EPRI-CPM / WIMS 69-group    |
| `VITAMINJ_175` | `IGN_VITAMINJ`  | 175    | VITAMIN-J 175-group (ORNL)  |
| `SANDII_620`   | `IGN_SANDII620` | 620    | SAND-II 620-group           |
| `XMAS_172`     | `IGN_XMAS172`   | 172    | XMAS 172-group              |
| `ECCO_33`      | `IGN_ECCO33`    | 33     | ECCO 33-group               |

#### `get_group_structure`
```julia
get_group_structure(id::GroupStructureId) -> NTuple
```
Return the group bounds for a built-in structure.

#### `num_groups`, `validate_group_bounds`, `find_group`
Utility functions for group structure manipulation.

### Weight functions

#### `constant_weight`
```julia
constant_weight(E) -> 1.0
```

#### `inv_e_weight`
```julia
inv_e_weight(E) -> 1/E
```

#### `maxwell_inv_e_fission`
```julia
maxwell_inv_e_fission(E; Eb=0.0253, Tb=0.0253, Ec=820.3e3, Tc=1.4e6) -> Float64
```
Three-region weight: Maxwellian + 1/E + fission spectrum (NJOY iwt=4).

#### `vitamin_e_weight`, `thermal_fission_fusion`, `tabulated_weight`
Additional weight functions matching NJOY conventions.

#### `get_weight_function`
```julia
get_weight_function(iwt; ...) -> Function
```
Return a weight function by NJOY `iwt` index.

### Integration and averaging

#### `group_integrate`
```julia
group_integrate(energies, values, group_bounds; law=LinLin) -> Vector{Float64}
```
Integrate a piecewise function over each group using exact analytical panel
integrals.

#### `group_average`
```julia
group_average(energies, xs_cols, mt_list, group_bounds;
              weight_fn=inv_e_weight, law=LinLin) -> MultiGroupXS
```
Flux-weighted group-averaged cross sections.

#### `group_average_shielded`
```julia
group_average_shielded(energies, total_xs, reaction_xs, mt_list,
                       group_bounds, sigma0;
                       weight_fn=inv_e_weight, law=LinLin) -> MultiGroupXS
```
Bondarenko self-shielded group averaging with effective weight
`w(E) * sigma_0 / (sigma_t(E) + sigma_0)`.

---

## MODER

### Types

#### `TapeEntry`
```julia
struct TapeEntry
    mat::Int; mf::Int; mt::Int
end
```

#### `TapeDirectory`
```julia
struct TapeDirectory
    tpid::String
    entries::Vector{TapeEntry}
end
```

#### `ENDFTapeSection`
One section stored as raw lines for lossless round-tripping.

#### `ENDFTapeMaterial`
One complete material from a tape.

### Functions

#### `read_tape_directory`
```julia
read_tape_directory(io::IO) -> TapeDirectory
```
Scan an ENDF tape and return a directory of every MAT/MF/MT section.

#### `materials`
```julia
materials(dir::TapeDirectory) -> Vector{Int}
```
Sorted unique MAT numbers on the tape.

#### `sections`
```julia
sections(dir::TapeDirectory, mat::Integer) -> Vector{Tuple{Int,Int}}
```
(MF, MT) pairs for a given material.

#### `extract_material`
```julia
extract_material(io::IO, mat::Integer) -> String
extract_material(io::IO, mat::Integer, out::IO) -> Int
```
Extract all lines for material `mat`.

#### `moder_copy`
```julia
moder_copy(input::IO, output::IO; mat=0) -> Int
```
Copy ENDF material(s) with optional MAT filter.

#### `merge_tapes`
```julia
merge_tapes(out::IO, inputs::IO...; tpid="NJOY.jl merged tape")
merge_tapes(output_path, input_paths...; tpid="NJOY.jl merged tape")
```
Concatenate materials from multiple tapes into one.

#### `write_tpid`
```julia
write_tpid(out::IO, text="")
```

#### `write_tend`
```julia
write_tend(out::IO)
```

#### `validate_tape`
```julia
validate_tape(io::IO) -> NamedTuple{(:valid, :errors)}
```

#### `read_endf_tape`
```julia
read_endf_tape(filename) -> Vector{ENDFTapeMaterial}
```
Structured round-trip read preserving raw lines.

#### `write_endf_tape`
```julia
write_endf_tape(filename, materials::Vector{ENDFTapeMaterial})
```
Write materials with FEND/MEND/TEND delimiters.

---

## ERRORR

### Types

#### `CovarianceBlock`
```julia
struct CovarianceBlock
    mt1::Int; mt2::Int; lb::Int; lt::Int
    energies::Vector{Float64}; data::Vector{Float64}
end
```
One NI-type covariance sub-subsection from ENDF MF33.

#### `CovarianceMatrix`
```julia
struct CovarianceMatrix
    mt1::Int; mt2::Int
    energy_groups::Vector{Float64}
    matrix::Matrix{Float64}
    is_relative::Bool
end
```
Symmetric multigroup covariance matrix.

#### `CovarianceData`
Parsed MF33 section containing NC and NI sub-subsections.

### Functions

#### `read_mf33`
```julia
read_mf33(io::IO, mat::Integer, mt::Integer) -> CovarianceData
```
Read one MF33 covariance section.

#### `expand_covariance_block`
```julia
expand_covariance_block(block::CovarianceBlock, egrid) -> Matrix{Float64}
```
Expand an NI-type block into a matrix on the given energy grid. Supports
LB = 0, 1, 2 (diagonal), LB = 5 (symmetric/full), LB = 6 (asymmetric).

#### `multigroup_covariance`
```julia
multigroup_covariance(blocks, group_bounds; is_relative=nothing, flux=nothing) -> CovarianceMatrix
```
Collapse NI-type blocks into a multigroup covariance matrix.

#### `process_covariance`
```julia
process_covariance(endf_file, group_structure; mts=:all) -> Vector{CovarianceMatrix}
```
Process all MF33 covariance data into multigroup covariance matrices.

#### `sandwich_covariance`
```julia
sandwich_covariance(J, cov_params) -> Matrix{Float64}
```
Compute `J * cov * J'` (sandwich rule). Result is symmetric.

#### `sensitivity_jacobian`
```julia
sensitivity_jacobian(sigma_func, params, group_bounds; dp=1e-6) -> Matrix{Float64}
```
Finite-difference Jacobian `d(sigma_g)/d(params)`.

#### `ni_covariance`
```julia
ni_covariance(ni_data, energy_grid) -> Matrix{Float64}
```

#### `nc_covariance`
```julia
nc_covariance(nc_data, energy_grid; xs_funcs=nothing) -> Matrix{Float64}
```

#### `is_symmetric`
```julia
is_symmetric(cm::CovarianceMatrix; atol=1e-12) -> Bool
```

#### `is_psd`
```julia
is_psd(cm::CovarianceMatrix; atol=1e-10) -> Bool
```
Check positive semi-definiteness.

---

## ACER

### Types

#### `ACEHeader`
```julia
ACEHeader(; zaid, awr, temp_mev, date="", comment="", mat_string="")
```
Metadata for an ACE table.

#### `ACETable`
Flat NXS/JXS/XSS representation for serialization.

#### `ACENeutronTable`
```julia
ACENeutronTable(; header, energy_grid, total_xs, absorption_xs,
                  elastic_xs, heating_numbers, reactions=ReactionXS[],
                  angular_elastic=nothing, angular=Dict(), iz=..., aw=...)
```
Type-safe structured representation of a continuous-energy neutron ACE table.

#### `ReactionXS`
```julia
struct ReactionXS
    mt::Int32; q_value::Float64; ty::Int32; ie_start::Int32
    xs::Vector{Float64}
end
```

#### `EquiprobableBins`
32 equally-probable cosine bins (33 bin edges).

#### `TabulatedAngular`
Tabulated angular distribution: cosines, PDF, CDF.

#### `AngularBlock`
Angular distributions for one reaction across incident energies.

### Builder functions

#### `build_ace`
```julia
build_ace(pendf::PointwiseMaterial; suffix="80c", awr=0.0,
          temperature=300.0, comment="", date="") -> ACETable
```
Build a flat ACETable from pointwise data.

#### `build_ace_from_pendf`
```julia
build_ace_from_pendf(pendf; suffix="80c", temp_kelvin=300.0,
                     comment="", mat_id=0) -> ACENeutronTable
```
Build a structured ACENeutronTable from a PointwiseMaterial.

#### `build_xss`
```julia
build_xss(table::ACENeutronTable) -> (nxs, jxs, xss, is_int)
```
Serialize into flat NXS(16), JXS(32), XSS arrays.

### Writer functions

#### `write_ace`
```julia
write_ace(io::IO, table::ACETable)
write_ace(io::IO, table::ACENeutronTable)
```
Write a Type 1 (ASCII) ACE file matching NJOY2016 `aceout()` format.

#### `write_ace_table`
```julia
write_ace_table(io::IO, table::ACENeutronTable; format=:ascii)
```
Write an ACENeutronTable as Type 1 ASCII.

#### `write_ace_directory`
```julia
write_ace_directory(io::IO, table, nxs; itype=1, filepath="filename",
                    route="route")
```
Write the MCNP xsdir directory line.

### Utility functions

#### `format_zaid`
```julia
format_zaid(za::Integer, suffix) -> String
format_zaid(z::Integer, a::Integer, suffix) -> String
```
Format a ZAID string (e.g., `"92235.80c"`).

#### `parse_zaid`
```julia
parse_zaid(s) -> (za::Int, suffix::String)
```

#### `temp_to_mev`
```julia
temp_to_mev(temp_kelvin) -> Float64
```
Convert Kelvin to MeV (kT).

#### `mev_to_temp`
```julia
mev_to_temp(kT_mev) -> Float64
```
Convert MeV (kT) to Kelvin.

### NXS/JXS constants

Index constants for the NXS and JXS arrays: `NXS_LEN2`, `NXS_IZAID`,
`NXS_NES`, `NXS_NTR`, `NXS_NR`, `NXS_NTRP`, `NXS_NTYPE`, `NXS_NDNF`,
`NXS_IS`, `NXS_IZ`, `NXS_IA`, `JXS_ESZ`, `JXS_NU`, `JXS_MTR`, `JXS_LQR`,
`JXS_TYR`, `JXS_LSIG`, `JXS_SIG`, `JXS_LAND`, `JXS_AND`, `JXS_LDLW`,
`JXS_DLW`, `JXS_GPD`, `JXS_FIS`, `JXS_END`, `ESZ_ENERGY`, `ESZ_TOTAL`,
`ESZ_DISAP`, `ESZ_ELASTIC`, `ESZ_HEATING`.

### Accessor functions

#### `nxs_length`, `nxs_nes`, `nxs_ntr`
```julia
nxs_length(table::ACETable) -> Int
nxs_nes(table::ACETable) -> Int
nxs_ntr(table::ACETable) -> Int
```

#### `esz_energies`, `esz_total`, `esz_elastic`
```julia
esz_energies(table::ACETable) -> view
esz_total(table::ACETable) -> view
esz_elastic(table::ACETable) -> view
```
Views into the ESZ block of a flat ACETable.

#### `ace_nes`, `ace_ntr`
```julia
ace_nes(table::ACENeutronTable) -> Int
ace_ntr(table::ACENeutronTable) -> Int
```
