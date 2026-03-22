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
    CrossSections{T<:Real}

Result container for computed cross sections at a single energy.
All values in barns.

Parameterized on `T` so that ForwardDiff dual numbers (or BigFloat, etc.)
propagate through without being truncated to Float64.
"""
struct CrossSections{T<:Real}
    total::T       # sigma_total
    elastic::T     # sigma_elastic
    fission::T     # sigma_fission
    capture::T     # sigma_capture (radiative)
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
    GT::Vector{Vector{Float64}} # total widths (Gamma_total) — raw from ENDF
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
    GT::Vector{Vector{Float64}}  # total widths (Gamma_total) — raw from ENDF
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
# SAMMY R-Matrix Limited (LRF=7)
# --------------------------------------------------------------------------

"""
    SAMMYParticlePair

One particle-pair description from LRF=7 data. Stores masses, charges,
spins, Q-value, penetrability/shift flags, MT, and boundary radii.
"""
struct SAMMYParticlePair
    ema::Float64       # mass of particle a [amu ratio to neutron]
    emb::Float64       # mass of particle b [amu ratio to neutron]
    kza::Int32         # charge of particle a
    kzb::Int32         # charge of particle b
    spina::Float64     # spin of particle a
    spinb::Float64     # spin of particle b
    qqq::Float64       # Q-value [eV]
    lpent::Int32       # penetrability flag (0=not calculated, 1=calculated)
    ishift::Int32      # shift factor flag
    mt::Int32          # ENDF reaction MT number
    pa::Float64        # boundary radius a
    pb::Float64        # boundary radius b
end

"""
    SAMMYSpinGroup

One spin-parity group from LRF=7 data. Contains channel definitions
and resonance parameters for all resonances in this group.
"""
struct SAMMYSpinGroup
    AJ::Float64                     # total angular momentum J
    parity::Float64                 # parity (+1 or -1)
    nchan::Int32                    # number of channels (excluding eliminated capture)
    # Per-channel arrays (length = nchan):
    ipp::Vector{Int32}              # particle-pair index for each channel
    lspin::Vector{Int32}            # orbital angular momentum for each channel
    chspin::Vector{Float64}         # channel spin for each channel
    bound::Vector{Float64}          # boundary condition for each channel
    rdeff::Vector{Float64}          # effective radius [fm] for each channel
    rdtru::Vector{Float64}          # true radius [fm] for each channel
    # Background R-matrix (LBK) per channel:
    backgr_type::Vector{Int32}      # LBK flag per channel (0=none, 2=Sammy, 3=Frohner)
    backgr_data::Vector{Vector{Float64}}  # background parameters per channel
    # Resonances in this group:
    nres::Int32                     # number of resonances
    eres::Vector{Float64}           # resonance energies [eV]
    gamgam::Vector{Float64}         # capture (gamma) widths [eV]
    gamma::Vector{Vector{Float64}}  # channel widths gamma[ch][res] [eV]
end

"""
    SAMMYParameters <: AbstractResonanceFormalism

SAMMY R-Matrix Limited resonance parameters (LRF=7, KRM=3 Reich-Moore).

This formalism uses the full R-matrix theory with explicit particle-pair
and spin-group structure, supporting:
- Multiple channels per J-value
- Proper channel-spin coupling
- Eliminated channels for capture
- Background R-matrix contributions (LBK)
"""
struct SAMMYParameters <: AbstractResonanceFormalism
    SPI::Float64                     # target spin
    AP::Float64                      # scattering radius [fm]
    KRM::Int32                       # R-matrix approximation (3=Reich-Moore)
    IFG::Int32                       # reduced width flag (0=widths, 1=reduced widths)
    particle_pairs::Vector{SAMMYParticlePair}
    spin_groups::Vector{SAMMYSpinGroup}
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
    ResonanceRange{P<:AbstractResonanceFormalism}

Container for one resonance energy range from ENDF File 2 (MF2/MT151).
Holds the formalism-specific parameters plus metadata about the range.

Parameterized on `P` (the concrete formalism type) to enable type-stable
dispatch in `cross_section(E, range)`.
"""
struct ResonanceRange{P<:AbstractResonanceFormalism}
    EL::Float64                      # lower energy bound [eV]
    EH::Float64                      # upper energy bound [eV]
    LRU::Int32                       # 0=scat-radius only, 1=resolved, 2=unresolved
    LRF::Int32                       # formalism flag
    LFW::Int32                       # fission width flag (URR only)
    NRO::Int32                       # energy-dependent scatt. radius flag
    NAPS::Int32                      # scattering radius control
    parameters::P
    ap_tab::Union{Nothing,TabulatedFunction}  # energy-dependent AP(E) when NRO!=0
end

# Convenience constructor without ap_tab (backward compatible)
function ResonanceRange(EL, EH, LRU, LRF, LFW, NRO, NAPS, parameters::P) where {P<:AbstractResonanceFormalism}
    ResonanceRange{P}(Float64(EL), Float64(EH), Int32(LRU), Int32(LRF),
                      Int32(LFW), Int32(NRO), Int32(NAPS), parameters, nothing)
end

function ResonanceRange(EL, EH, LRU, LRF, LFW, NRO, NAPS, parameters::P,
                        ap_tab::Union{Nothing,TabulatedFunction}) where {P<:AbstractResonanceFormalism}
    ResonanceRange{P}(Float64(EL), Float64(EH), Int32(LRU), Int32(LRF),
                      Int32(LFW), Int32(NRO), Int32(NAPS), parameters, ap_tab)
end
