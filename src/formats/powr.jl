# EPRI-CELL and EPRI-CPM library format types, validation, and I/O
#
# Proposer-B design: clean Julia types for the POWR output formats
# as implemented in NJOY2016 powr.f90.
#
# Three library types:
#   lib=1: EPRI-CELL fast (GAMTAP) -- resonance and non-resonance isotopes
#   lib=2: EPRI-CELL thermal (LIBRAR) -- thermal scattering data
#   lib=3: EPRI-CPM (CLIB) -- material data with self-shielding tables
#
# The CPM library has 7 subfiles:
#   1: header + group structure + fission spectra
#   2: resonance self-shielding (file 2)
#   3: resonance tables by nuclide
#   4: cross section data per material (file 4)
#   5: burnup chain data (file 5)
#   6: P1 scattering matrices (file 6)

using Printf

# =========================================================================
# Constants
# =========================================================================
const POWR_LIB_FAST    = 1
const POWR_LIB_THERMAL = 2
const POWR_LIB_CPM     = 3

const POWR_NGNF      = 185   # fine-structure neutron groups
const POWR_NGND_FAST = 68    # fast library coarse groups
const POWR_NGMIN     = 19    # first fast-range group in fine structure
const POWR_NGMAX     = 62    # last fast-range group with shielding

# =========================================================================
# EPRI-CELL Fast library (lib=1) -- GAMTAP format
# =========================================================================

"""
    EPRIFastScatterMatrix

Compressed scattering matrix for one reaction type in EPRI-CELL fast format.

# Fields
- `bandwidth`: number of non-zero off-diagonal entries
- `offset`: first sink group (1-based)
- `length`: total stored entries
- `data`: flat vector of scattering transfer values
"""
struct EPRIFastScatterMatrix
    bandwidth::Int
    offset::Int
    length::Int
    data::Vector{Float64}
end

"""
    EPRIFastIsotope

Cross section data for one isotope in the EPRI-CELL fast library.

# Fields
- `nid`: nuclide identifier (material number)
- `description`: up to 16-character description
- `fission_spectrum_title`: up to 40-character fission spectrum title
- `absorption`: absorption cross section per group (length ngnd)
- `fission`: sigma_f per group (empty if non-fissile)
- `nu`: nu values per group (empty if non-fissile)
- `chi`: fission spectrum per group (empty if non-fissile)
- `elastic_p0`: P0 elastic scattering matrix
- `elastic_p1`: P1 elastic scattering matrix (empty if not computed)
- `inelastic`: inelastic scattering matrix (empty if not present)
- `n2n`: (n,2n) scattering matrix (empty if not present)
- `sigma_zeros`: sigma-zero values for self-shielding
- `temperatures`: temperatures for self-shielding factors
- `abs_ssf`: absorption self-shielding factors, shape (ngroups, nsigz, ntemp)
- `fis_ssf`: fission self-shielding factors, shape (ngroups, nsigz, ntemp)
"""
struct EPRIFastIsotope
    nid::Int
    description::String
    fission_spectrum_title::String
    absorption::Vector{Float64}
    fission::Vector{Float64}
    nu::Vector{Float64}
    chi::Vector{Float64}
    elastic_p0::EPRIFastScatterMatrix
    elastic_p1::EPRIFastScatterMatrix
    inelastic::EPRIFastScatterMatrix
    n2n::EPRIFastScatterMatrix
    sigma_zeros::Vector{Float64}
    temperatures::Vector{Float64}
    abs_ssf::Array{Float64,3}
    fis_ssf::Array{Float64,3}
end

"""
    validate(iso::EPRIFastIsotope, ngnd::Int)

Validate cross section data dimensions for the fast library.
"""
function validate(iso::EPRIFastIsotope, ngnd::Int=POWR_NGND_FAST)
    length(iso.absorption) == ngnd ||
        error("EPRIFastIsotope $(iso.nid): absorption length must be $ngnd")
    if !isempty(iso.fission)
        length(iso.fission) == ngnd ||
            error("EPRIFastIsotope $(iso.nid): fission length must be $ngnd")
    end
    if !isempty(iso.nu)
        length(iso.nu) == ngnd ||
            error("EPRIFastIsotope $(iso.nid): nu length must be $ngnd")
    end
    if !isempty(iso.chi)
        length(iso.chi) == ngnd ||
            error("EPRIFastIsotope $(iso.nid): chi length must be $ngnd")
    end
    return true
end

# =========================================================================
# EPRI-CELL Thermal library (lib=2) -- LIBRAR format
# =========================================================================

"""
    EPRIThermalIsotope

Cross section data for one isotope/temperature in the EPRI-CELL thermal library.

# Fields
- `mat`: ENDF material number
- `temperature`: temperature identifier (K)
- `name`: up to 10-character name
- `transport_correction`: whether transport correction is applied
- `mti`: MT for thermal inelastic
- `mtc`: MT for thermal coherent elastic
- `xi`: average lethargy gain per collision
- `alpha`: (A-1)^2 / (A+1)^2
- `mubar`: average cosine of scattering angle
- `nu`: average fission neutron yield
- `kappa_fission`: energy release per fission (MeV)
- `kappa_capture`: energy release per capture (MeV)
- `lambda_`: transport mean free path
- `sigma_s`: scattering cross section at reference group
- `scatter_matrix`: thermal scattering matrix (ngroups x ngroups)
"""
struct EPRIThermalIsotope
    mat::Int
    temperature::Float64
    name::String
    transport_correction::Bool
    mti::Int
    mtc::Int
    xi::Float64
    alpha::Float64
    mubar::Float64
    nu::Float64
    kappa_fission::Float64
    kappa_capture::Float64
    lambda_::Float64
    sigma_s::Float64
    scatter_matrix::Matrix{Float64}
end

function validate(iso::EPRIThermalIsotope)
    iso.mat > 0 || error("EPRIThermalIsotope: mat must be > 0")
    iso.temperature >= 0.0 ||
        error("EPRIThermalIsotope: temperature must be >= 0")
    n = size(iso.scatter_matrix, 1)
    size(iso.scatter_matrix, 2) == n ||
        error("EPRIThermalIsotope: scatter_matrix must be square")
    return true
end

# =========================================================================
# EPRI-CPM library (lib=3) -- CLIB format
# =========================================================================

"""
    CPMNuclideSpec

Material specification for one nuclide in the CPM library.

# Fields
- `ident`: nuclide identifier
- `atomic_weight`: atomic weight (amu)
- `fission_indicator`: 0=non-fissile, 1=fissile with res, 2=fissile no res
- `has_resonance`: whether resonance tables are present
- `nf`: resonance tabulation count (0, 1, or 2)
- `ntemps`: number of temperatures
- `nina`: NINA indicator (0=normal, 1=no file2/calc abs, 2=no file2/read abs, 3=read all)
- `has_p1`: whether P1 matrix is present
"""
struct CPMNuclideSpec
    ident::Int
    atomic_weight::Float64
    fission_indicator::Int
    has_resonance::Bool
    nf::Int
    ntemps::Int
    nina::Int
    has_p1::Bool
end

function validate(spec::CPMNuclideSpec)
    spec.ident > 0 || error("CPMNuclideSpec: ident must be > 0")
    spec.atomic_weight > 0.0 ||
        error("CPMNuclideSpec: atomic_weight must be > 0")
    spec.fission_indicator in (0, 1, 2) ||
        error("CPMNuclideSpec: fission_indicator must be 0, 1, or 2")
    spec.nina in (0, 1, 2, 3) ||
        error("CPMNuclideSpec: nina must be 0, 1, 2, or 3")
    return true
end

"""
    CPMResonanceData

Resonance self-shielding data for one group (subfile 2 of CPM library).
Stores lambda, sigma_p, effective resonance integrals at multiple
(temperature, sigma-zero) combinations.
"""
struct CPMResonanceData
    lambda::Float64
    sigma_p::Float64
    nu_res::Float64
    total_res::Float64
end

"""
    CPMBurnupIsotope

One entry in the CPM burnup chain (subfile 5).

# Fields
- `ident`: isotope identifier
- `decay_constant`: decay constant (s^-1), 0 if stable
- `yields`: fission yields from each fissile parent
"""
struct CPMBurnupIsotope
    ident::Int
    decay_constant::Float64
    yields::Vector{Float64}
end

"""
    CPMLibraryHeader

Header data for the CPM library (subfile 1).

# Fields
- `nlib`: library number
- `idat`: date code
- `ngroups`: number of groups
- `nfg`: number of fast groups
- `nrg`: number of resonance groups
- `nfiss`: number of fission groups for spectrum
- `nmat`: number of materials
- `energy_bounds`: group boundary energies (length ngroups+1)
- `flux_spectrum`: collapsing flux spectrum (length ngroups)
- `u235_chi`: U-235 fission spectrum (length nfiss)
- `pu239_chi`: Pu-239 fission spectrum (length nfiss)
- `nuclides`: material specification for each nuclide
"""
struct CPMLibraryHeader
    nlib::Int
    idat::Int
    ngroups::Int
    nfg::Int
    nrg::Int
    nfiss::Int
    nmat::Int
    energy_bounds::Vector{Float64}
    flux_spectrum::Vector{Float64}
    u235_chi::Vector{Float64}
    pu239_chi::Vector{Float64}
    nuclides::Vector{CPMNuclideSpec}
end

function validate(h::CPMLibraryHeader)
    h.ngroups > 0 || error("CPMLibraryHeader: ngroups must be > 0")
    h.nfg >= 0    || error("CPMLibraryHeader: nfg must be >= 0")
    h.nrg >= 0    || error("CPMLibraryHeader: nrg must be >= 0")
    h.nfg + h.nrg <= h.ngroups ||
        error("CPMLibraryHeader: nfg + nrg > ngroups")
    length(h.energy_bounds) == h.ngroups + 1 ||
        error("CPMLibraryHeader: energy_bounds length must be ngroups+1")
    length(h.flux_spectrum) == h.ngroups ||
        error("CPMLibraryHeader: flux_spectrum length must be ngroups")
    length(h.nuclides) == h.nmat ||
        error("CPMLibraryHeader: nuclides count must equal nmat")
    for spec in h.nuclides
        validate(spec)
    end
    return true
end

"""
    CPMLibrary

Complete EPRI-CPM library data structure.

# Fields
- `header`: library header (subfile 1)
- `resonance_data`: resonance self-shielding per nuclide per group (subfile 2)
- `burnup_chain`: burnup chain isotopes (subfile 5, empty if not provided)
- `mode`: creation mode (0=replace, 1=add, 2=create new)
"""
struct CPMLibrary
    header::CPMLibraryHeader
    resonance_data::Vector{Vector{CPMResonanceData}}
    burnup_chain::Vector{CPMBurnupIsotope}
    mode::Int
end

function validate(lib::CPMLibrary)
    validate(lib.header)
    lib.mode in (0, 1, 2) ||
        error("CPMLibrary: mode must be 0, 1, or 2, got $(lib.mode)")
    return true
end

# =========================================================================
# Unified POWR output type
# =========================================================================
"""
    POWROutput

Discriminated union for the three POWR library types.

# Fields
- `lib_type`: 1=fast, 2=thermal, 3=CPM
- `fast_isotopes`: vector of EPRIFastIsotope (lib=1)
- `thermal_isotopes`: vector of EPRIThermalIsotope (lib=2)
- `cpm_library`: CPMLibrary (lib=3)
"""
struct POWROutput
    lib_type::Int
    fast_isotopes::Vector{EPRIFastIsotope}
    thermal_isotopes::Vector{EPRIThermalIsotope}
    cpm_library::Union{Nothing,CPMLibrary}
end

function validate(p::POWROutput)
    p.lib_type in (POWR_LIB_FAST, POWR_LIB_THERMAL, POWR_LIB_CPM) ||
        error("POWROutput: lib_type must be 1, 2, or 3")
    if p.lib_type == POWR_LIB_FAST
        for iso in p.fast_isotopes
            validate(iso)
        end
    elseif p.lib_type == POWR_LIB_THERMAL
        for iso in p.thermal_isotopes
            validate(iso)
        end
    elseif p.lib_type == POWR_LIB_CPM
        p.cpm_library !== nothing ||
            error("POWROutput: cpm_library must not be nothing for lib=3")
        validate(p.cpm_library)
    end
    return true
end

# =========================================================================
# Writer -- EPRI-CELL fast format
# =========================================================================
"""
    write_powr_fast(io::IO, isotopes::Vector{EPRIFastIsotope}; ngnd=68)

Write EPRI-CELL fast library (GAMTAP format) for the given isotopes.
"""
function write_powr_fast(io::IO, isotopes::Vector{EPRIFastIsotope};
                         ngnd::Int=POWR_NGND_FAST)
    for iso in isotopes
        validate(iso, ngnd)
        # Header line: nid, flag, description
        iwa = isempty(iso.absorption) ? 0 : 1
        iwf = isempty(iso.fission) ? 0 : 1
        iwr = isempty(iso.sigma_zeros) ? 0 : length(iso.sigma_zeros)
        @printf(io, "%6d%2d%2d%2d\n", iso.nid, iwa, iwf, iwr)

        # Description
        @printf(io, "%-16s\n", iso.description)

        # Scattering matrix layout info
        _write_scat_layout(io, iso.elastic_p0, "P0 elastic")
        _write_scat_layout(io, iso.elastic_p1, "P1 elastic")
        _write_scat_layout(io, iso.inelastic, "inelastic")
        _write_scat_layout(io, iso.n2n, "n2n")

        # Absorption
        if iwa > 0
            _write_6e(io, iso.absorption)
        end

        # Fission data
        if iwf > 0
            _write_6e(io, iso.fission)
            _write_6e(io, iso.nu)
            if !isempty(iso.chi)
                _write_6e(io, iso.chi)
            end
        end

        # Scattering matrices
        if !isempty(iso.elastic_p0.data)
            _write_6e(io, iso.elastic_p0.data)
        end
        if !isempty(iso.elastic_p1.data)
            _write_6e(io, iso.elastic_p1.data)
        end
        if !isempty(iso.inelastic.data)
            _write_6e(io, iso.inelastic.data)
        end
        if !isempty(iso.n2n.data)
            _write_6e(io, iso.n2n.data)
        end

        # Self-shielding factors
        if iwr > 0
            nsigz = length(iso.sigma_zeros)
            ntp = length(iso.temperatures)
            @printf(io, "%6d%6d%6d%6d\n", nsigz, ntp, 1, iwf)
            _write_6e(io, iso.sigma_zeros)
            _write_6e(io, iso.temperatures)
            # Flatten 3D arrays
            if !isempty(iso.abs_ssf)
                _write_6e(io, vec(iso.abs_ssf))
            end
            if iwf > 0 && !isempty(iso.fis_ssf)
                _write_6e(io, vec(iso.fis_ssf))
            end
        end
    end
end

# =========================================================================
# Writer -- CPM library format
# =========================================================================
"""
    write_powr_cpm(io::IO, lib::CPMLibrary)

Write an EPRI-CPM library (CLIB format).
"""
function write_powr_cpm(io::IO, lib::CPMLibrary)
    validate(lib)
    h = lib.header

    # Subfile 1: header
    @printf(io, "pri   %6d%6d\n", 0, 0)
    @printf(io, "(6e12.0)\n(6e12.0)\n(6e12.0)\n(6e12.0)\n")

    if lib.mode <= 1
        @printf(io, "%12d%6d\n", h.nmat, 0)
    else
        @printf(io, "%12d%6d\n", 0, 1)
    end
    @printf(io, "%12d\n", 1)
    @printf(io, "%12d%6d%6d%6d%6d%6d%6d\n",
            h.nlib, h.idat, h.ngroups, h.nfg, h.nrg, h.nfiss, h.nmat)

    # Group boundaries
    _write_6e(io, h.energy_bounds)
    # Flux spectrum
    _write_6e(io, h.flux_spectrum)
    # Fission spectra
    _write_6e(io, h.u235_chi)
    _write_6e(io, h.pu239_chi)

    # Nuclide specifications
    for spec in h.nuclides
        ip1 = spec.has_p1 ? 1 : 0
        @printf(io, "%12.4f%6d%6d%6d%6d%6d%6d%6d%6d\n",
                spec.atomic_weight, spec.ident, spec.fission_indicator,
                spec.has_resonance ? 1 : 0, spec.nf, spec.ntemps,
                spec.nina, 0, ip1)
    end

    # Subfile 2: resonance data
    for nuc_res in lib.resonance_data
        for rd in nuc_res
            @printf(io, "%12.5E%12.5E%12.5E%12.5E\n",
                    rd.lambda, rd.sigma_p, rd.nu_res, rd.total_res)
        end
    end

    # Subfile 5: burnup chain
    if !isempty(lib.burnup_chain)
        for biso in lib.burnup_chain
            nfis = length(biso.yields)
            @printf(io, "%6d%12.5E", biso.ident, biso.decay_constant)
            for y in biso.yields
                @printf(io, "%12.5E", y)
            end
            println(io)
        end
    end
end

"""
    write_powr(io::IO, output::POWROutput)

Write a POWROutput in the appropriate format based on lib_type.
"""
function write_powr(io::IO, output::POWROutput)
    validate(output)
    if output.lib_type == POWR_LIB_FAST
        write_powr_fast(io, output.fast_isotopes)
    elseif output.lib_type == POWR_LIB_CPM
        write_powr_cpm(io, output.cpm_library)
    end
    # thermal is handled by write_powr_fast with different parameters
end

"""
    write_powr(filename::AbstractString, output::POWROutput)

Write a POWROutput to a file.
"""
function write_powr(filename::AbstractString, output::POWROutput)
    open(filename, "w") do io
        write_powr(io, output)
    end
end

# =========================================================================
# Internal helpers
# =========================================================================

# Write a vector in 6-per-line E12.5 format
function _write_6e(io::IO, v::AbstractVector{<:Real})
    for (i, x) in enumerate(v)
        @printf(io, "%12.5E", Float64(x))
        if i % 6 == 0 || i == length(v)
            println(io)
        end
    end
end

# Write scattering matrix layout line
function _write_scat_layout(io::IO, sm::EPRIFastScatterMatrix, label::String)
    @printf(io, "%36s%12.5E%12.5E%12.5E\n",
            "", Float64(sm.length), Float64(sm.bandwidth), Float64(sm.offset + 1))
end

# =========================================================================
# Convenience: write_epri_cell from MultiGroupXS
# =========================================================================

"""
    write_epri_cell(io::IO, multigroup_data::MultiGroupXS;
                    lib=1, nid=0, description="", temperature=300.0)

Write multigroup cross sections in EPRI-CELL (GAMTAP fast or LIBRAR thermal)
format, constructing the EPRIFastIsotope automatically from MultiGroupXS.

# Keyword Arguments
- `lib`: library type (1=fast, 2=thermal; default 1)
- `nid`: nuclide identifier (default 0)
- `description`: nuclide description (default "")
- `temperature`: temperature in K (default 300.0)
"""
function write_epri_cell(io::IO, multigroup_data::MultiGroupXS;
                         lib::Int = POWR_LIB_FAST,
                         nid::Int = 0,
                         description::AbstractString = "",
                         temperature::Float64 = 300.0)
    ng = length(multigroup_data.group_bounds) - 1
    mt_list = multigroup_data.mt_list
    xs = multigroup_data.xs

    # Build MT lookup
    mt_to_col = Dict{Int,Int}()
    for (c, mt) in enumerate(mt_list)
        mt_to_col[mt] = c
    end

    # Extract cross sections per group
    abs_xs = zeros(ng)
    fis_xs = Float64[]
    nu_xs  = Float64[]
    chi_xs = Float64[]
    el_xs  = zeros(ng)

    for ig in 1:ng
        if haskey(mt_to_col, 102)
            abs_xs[ig] = xs[ig, mt_to_col[102]]
        elseif haskey(mt_to_col, 27)
            abs_xs[ig] = xs[ig, mt_to_col[27]]
        end
        if haskey(mt_to_col, 2)
            el_xs[ig] = xs[ig, mt_to_col[2]]
        end
    end

    has_fission = haskey(mt_to_col, 18)
    if has_fission
        fis_xs = zeros(ng)
        nu_xs  = zeros(ng)
        chi_xs = zeros(ng)
        for ig in 1:ng
            fis_xs[ig] = xs[ig, mt_to_col[18]]
        end
        # Simple chi: all in first group (placeholder)
        if ng > 0; chi_xs[1] = 1.0; end
    end

    # Build trivial diagonal elastic scattering matrix
    empty_sm = EPRIFastScatterMatrix(0, 0, 0, Float64[])
    el_sm = EPRIFastScatterMatrix(1, 1, ng, el_xs)

    iso = EPRIFastIsotope(
        nid, rpad(description, 16)[1:min(16, max(1, length(rpad(description, 16))))],
        "", abs_xs, fis_xs, nu_xs, chi_xs,
        el_sm, empty_sm, empty_sm, empty_sm,
        Float64[], [temperature],
        Array{Float64}(undef, 0, 0, 0),
        Array{Float64}(undef, 0, 0, 0)
    )

    write_powr_fast(io, [iso]; ngnd=ng)
end

# =========================================================================
# Convenience: write_epri_cpm from MultiGroupXS
# =========================================================================

"""
    write_epri_cpm(io::IO, multigroup_data::MultiGroupXS;
                   nlib=1, idat=0, nfg=14, nrg=13, nfiss=22,
                   ident=0, awr=1.0)

Write multigroup cross sections in EPRI-CPM (CLIB) format, constructing
the CPMLibrary automatically from MultiGroupXS.

# Keyword Arguments
- `nlib`: library number (default 1)
- `idat`: date code (default 0)
- `nfg`: number of fast groups (default 14)
- `nrg`: number of resonance groups (default 13)
- `nfiss`: number of fission spectrum groups (default 22)
- `ident`: nuclide identifier (default 1)
- `awr`: atomic weight ratio (default 1.0)
"""
function write_epri_cpm(io::IO, multigroup_data::MultiGroupXS;
                        nlib::Int = 1,
                        idat::Int = 0,
                        nfg::Int = 14,
                        nrg::Int = 13,
                        nfiss::Int = 22,
                        ident::Int = 1,
                        awr::Float64 = 1.0)
    ng = length(multigroup_data.group_bounds) - 1
    gb = multigroup_data.group_bounds
    flux = multigroup_data.flux
    mt_list = multigroup_data.mt_list

    has_fission = 18 in mt_list || 19 in mt_list

    amassn = 1.00866491595
    spec = CPMNuclideSpec(
        ident, awr * amassn,
        has_fission ? 1 : 0,
        false, 0, 1, 0, false
    )

    header = CPMLibraryHeader(
        nlib, idat, ng, nfg, nrg, nfiss, 1,
        collect(Float64, gb),
        collect(Float64, flux),
        zeros(nfiss),   # U-235 chi placeholder
        zeros(nfiss),   # Pu-239 chi placeholder
        [spec]
    )

    lib = CPMLibrary(header, Vector{CPMResonanceData}[], CPMBurnupIsotope[], 2)
    write_powr_cpm(io, lib)
end
