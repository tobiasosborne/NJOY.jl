# Resonance formalism types
#
# Type hierarchy for ENDF File 2 (MF2) resonance parameters.
# Covers all formalisms supported by NJOY2016's reconr module:
#   LRF=1: Single-Level Breit-Wigner (SLBW)
#   LRF=2: Multi-Level Breit-Wigner (MLBW)
#   LRF=3: Reich-Moore (RM)
#   LRF=4: Adler-Adler (AA)
#   LRU=2: Unresolved resonance region (3 variants)

"""
    AbstractResonanceFormalism

Abstract supertype for all ENDF resonance parameter formalisms.
"""
abstract type AbstractResonanceFormalism end

"""
    CrossSections

Result container for computed cross sections at a single energy.
All values in barns.
"""
struct CrossSections
    total::Float64       # sigma_total
    elastic::Float64     # sigma_elastic
    fission::Float64     # sigma_fission
    capture::Float64     # sigma_capture (radiative)
end

CrossSections() = CrossSections(0.0, 0.0, 0.0, 0.0)

function Base.:+(a::CrossSections, b::CrossSections)
    CrossSections(a.total + b.total, a.elastic + b.elastic,
                  a.fission + b.fission, a.capture + b.capture)
end

function Base.:*(s::Real, xs::CrossSections)
    CrossSections(s * xs.total, s * xs.elastic,
                  s * xs.fission, s * xs.capture)
end
Base.:*(xs::CrossSections, s::Real) = s * xs

# --------------------------------------------------------------------------
# Resolved resonance formalisms
# --------------------------------------------------------------------------

"""
    SLBWParameters <: AbstractResonanceFormalism

Single-Level Breit-Wigner resonance parameters (LRF=1).

Each resonance l-group contains:
- `Er`: resonance energies [eV]
- `J`: spin values
- `Gn`: neutron widths [eV]  (Gamma_n)
- `Gg`: gamma (capture) widths [eV]  (Gamma_gamma)
- `Gf`: fission widths [eV]  (Gamma_f)
- `Gx`: competitive widths [eV]  (Gamma_x)
"""
struct SLBWParameters <: AbstractResonanceFormalism
    NLS::Int32                   # number of l-values
    SPI::Float64                 # target spin I
    AP::Float64                  # scattering radius [fm]
    l_values::Vector{Int32}      # orbital angular momentum for each group
    AWRI::Vector{Float64}        # atomic weight ratio per l-group
    QX::Vector{Float64}          # competitive Q-value per l-group
    LRX::Vector{Int32}           # competitive width flag per l-group
    Er::Vector{Vector{Float64}}  # resonance energies per l-group
    AJ::Vector{Vector{Float64}} # spin J per l-group
    Gn::Vector{Vector{Float64}} # neutron widths (Gamma_n)
    Gg::Vector{Vector{Float64}} # gamma widths (Gamma_gamma)
    Gf::Vector{Vector{Float64}} # fission widths (Gamma_f)
    Gx::Vector{Vector{Float64}} # competitive widths (Gamma_x)
end

"""
    MLBWParameters <: AbstractResonanceFormalism

Multi-Level Breit-Wigner resonance parameters (LRF=2).
Same storage layout as SLBW but uses multi-level interference formula.
"""
struct MLBWParameters <: AbstractResonanceFormalism
    NLS::Int32
    SPI::Float64
    AP::Float64
    l_values::Vector{Int32}
    AWRI::Vector{Float64}        # atomic weight ratio per l-group
    QX::Vector{Float64}          # competitive Q-value per l-group
    LRX::Vector{Int32}           # competitive width flag per l-group
    Er::Vector{Vector{Float64}}
    AJ::Vector{Vector{Float64}}
    Gn::Vector{Vector{Float64}}
    Gg::Vector{Vector{Float64}}
    Gf::Vector{Vector{Float64}}
    Gx::Vector{Vector{Float64}}
end

"""
    ReichMooreParameters <: AbstractResonanceFormalism

Reich-Moore resonance parameters (LRF=3).

Uses a 3-channel R-matrix formulation with explicit fission channels.
- `Gfa`: first fission width (Gamma_{fA})
- `Gfb`: second fission width (Gamma_{fB})
"""
struct ReichMooreParameters <: AbstractResonanceFormalism
    NLS::Int32
    SPI::Float64
    AP::Float64
    LAD::Int32                   # angular distributions flag
    l_values::Vector{Int32}
    AWRI::Vector{Float64}        # atomic weight ratio per l-group
    APL::Vector{Float64}         # l-dependent scattering radius per l-group
    Er::Vector{Vector{Float64}}
    AJ::Vector{Vector{Float64}}
    Gn::Vector{Vector{Float64}}
    Gg::Vector{Vector{Float64}}
    Gfa::Vector{Vector{Float64}} # first fission width
    Gfb::Vector{Vector{Float64}} # second fission width
end

"""
    AdlerAdlerParameters <: AbstractResonanceFormalism

Adler-Adler resonance parameters (LRF=4).

Uses background polynomial + resonance-term parameterization with
coefficients AT, BT (total), AC, BC (capture), AF, BF (fission).
"""
struct AdlerAdlerParameters <: AbstractResonanceFormalism
    NLS::Int32
    SPI::Float64
    AP::Float64
    AT::Vector{Float64}   # total background coefficients
    BT::Vector{Float64}   # total resonance coefficients
    AC::Vector{Float64}   # capture background coefficients
    BC::Vector{Float64}   # capture resonance coefficients
    AF::Vector{Float64}   # fission background coefficients
    BF::Vector{Float64}   # fission resonance coefficients
    DET::Vector{Float64}  # resonance determinant energies
    DWT::Vector{Float64}  # resonance determinant widths
end

# --------------------------------------------------------------------------
# Unresolved resonance region
# --------------------------------------------------------------------------

"""
    UnresolvedParameters <: AbstractResonanceFormalism

Unresolved resonance region parameters (LRU=2).

Three variants controlled by LFW and LRF:
- `variant = :energy_independent` (LFW=0, LRF=1): widths do not depend on energy
- `variant = :fission_width`      (LFW=1, LRF=1): only fission widths vary with energy
- `variant = :energy_dependent`   (LFW=0/1, LRF=2): all parameters energy-dependent
"""
struct UnresolvedParameters <: AbstractResonanceFormalism
    variant::Symbol              # :energy_independent, :fission_width, :energy_dependent
    SPI::Float64                 # target spin
    AP::Float64                  # scattering radius
    LSSF::Int32                  # self-shielding flag
    NLS::Int32                   # number of l-values
    l_values::Vector{Int32}      # orbital angular momenta
    D::Vector{Vector{Float64}}   # mean level spacings
    AJ::Vector{Vector{Float64}}  # spin values
    AMUN::Vector{Vector{Float64}}  # degrees of freedom for neutron width
    GNO::Vector{Vector{Float64}}   # average neutron widths
    GG::Vector{Vector{Float64}}    # average gamma widths
    GF::Vector{Vector{Float64}}    # average fission widths
    GX::Vector{Vector{Float64}}    # average competitive widths
end

# --------------------------------------------------------------------------
# Resonance range container
# --------------------------------------------------------------------------

"""
    ResonanceRange

Container for one resonance energy range from ENDF File 2 (MF2/MT151).
Holds the formalism-specific parameters plus metadata about the range.
"""
struct ResonanceRange
    EL::Float64                      # lower energy bound [eV]
    EH::Float64                      # upper energy bound [eV]
    LRU::Int32                       # 0=scat-radius only, 1=resolved, 2=unresolved
    LRF::Int32                       # formalism flag
    LFW::Int32                       # fission width flag (URR only)
    NRO::Int32                       # energy-dependent scatt. radius flag
    NAPS::Int32                      # scattering radius control
    parameters::AbstractResonanceFormalism
end
