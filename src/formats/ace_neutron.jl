# ACE neutron table: structured types for continuous-energy neutron data
#
# ReactionXS     -- per-reaction cross section data
# AngularDist    -- angular distribution for one energy point
# AngularBlock   -- angular distributions for one reaction
# ACENeutronTable -- high-level structured representation
#
# The writer (ace_writer.jl) serializes ACENeutronTable into the flat
# NXS/JXS/XSS arrays and writes Type 1 (ASCII) format matching NJOY2016
# aceout() exactly.

# =========================================================================
# Reaction cross section
# =========================================================================
"""
    ReactionXS

Cross section for a single reaction (MT != 1, 2) on a shared energy grid.
`ie_start` is the 1-based index into the energy grid where this reaction
threshold lies. `xs` has length `NES - ie_start + 1`.
"""
struct ReactionXS
    mt::Int32            # ENDF MT number
    q_value::Float64     # Q-value in MeV (ACE convention)
    ty::Int32            # reaction type: signed neutron multiplicity
    ie_start::Int32      # 1-based threshold index into energy grid
    xs::Vector{Float64}  # cross section values (barns)
end

# =========================================================================
# Angular distribution types
# =========================================================================
"""
    EquiprobableBins

32 equally-probable cosine bins (33 bin edges from -1 to +1).
This is the standard MCNP representation for angular data.
"""
struct EquiprobableBins
    cosines::Vector{Float64}  # exactly 33 values

    function EquiprobableBins(c::AbstractVector{<:Real})
        length(c) == 33 ||
            throw(ArgumentError("need 33 cosine values, got $(length(c))"))
        new(Float64.(c))
    end
end

"""
    TabulatedAngular

Tabulated angular distribution at one incident energy: cosines, PDF, CDF.
Used for Law 61 (new format) angular distributions.
"""
struct TabulatedAngular
    iflag::Int32               # interpolation flag (1=histogram, 2=lin-lin)
    cosine::Vector{Float64}
    pdf::Vector{Float64}
    cdf::Vector{Float64}
end

"""
    AngularBlock

Angular distributions for one reaction across incident energies.
Each entry in `distributions` can be:
  - `EquiprobableBins` (32 equally-probable bins)
  - `TabulatedAngular` (tabulated PDF/CDF)
  - `nothing` (isotropic at that energy)
"""
struct AngularBlock
    energies::Vector{Float64}
    distributions::Vector{Union{EquiprobableBins, TabulatedAngular, Nothing}}
end

# =========================================================================
# Main ACE neutron table (high-level, structured)
# =========================================================================
"""
    ACENeutronTable

Type-safe representation of a continuous-energy neutron ACE table.

The ESZ block holds 5 parallel arrays of length NES:
  energy_grid, total_xs, absorption_xs, elastic_xs, heating_numbers.

Reactions (MT != 1, 2) are stored in `reactions::Vector{ReactionXS}`.
Angular data for elastic (and optionally other reactions) is in
`angular_elastic` and `angular`.
"""
struct ACENeutronTable
    header::ACEHeader
    energy_grid::Vector{Float64}
    total_xs::Vector{Float64}
    absorption_xs::Vector{Float64}
    elastic_xs::Vector{Float64}
    heating_numbers::Vector{Float64}
    reactions::Vector{ReactionXS}
    angular_elastic::Union{AngularBlock, Nothing}
    angular::Dict{Int, AngularBlock}   # MT -> angular block
    iz::Vector{Int32}    # (IZ,AW) pairs: iz values (16 entries)
    aw::Vector{Float64}  # (IZ,AW) pairs: aw values (16 entries)
end

"""
    ACENeutronTable(; header, energy_grid, total_xs, absorption_xs,
                      elastic_xs, heating_numbers, reactions=ReactionXS[],
                      angular_elastic=nothing, angular=Dict{Int,AngularBlock}(),
                      iz=zeros(Int32,16), aw=zeros(Float64,16))

Construct with validation that ESZ arrays all have length NES.
"""
function ACENeutronTable(; header::ACEHeader,
                           energy_grid::AbstractVector{<:Real},
                           total_xs::AbstractVector{<:Real},
                           absorption_xs::AbstractVector{<:Real},
                           elastic_xs::AbstractVector{<:Real},
                           heating_numbers::AbstractVector{<:Real},
                           reactions::Vector{ReactionXS} = ReactionXS[],
                           angular_elastic::Union{AngularBlock, Nothing} = nothing,
                           angular::Dict{Int, AngularBlock} = Dict{Int, AngularBlock}(),
                           iz::AbstractVector{<:Integer} = zeros(Int32, 16),
                           aw::AbstractVector{<:Real} = zeros(Float64, 16))
    nes = length(energy_grid)
    for (nm, a) in [("total_xs", total_xs), ("absorption_xs", absorption_xs),
                    ("elastic_xs", elastic_xs), ("heating_numbers", heating_numbers)]
        length(a) == nes || throw(ArgumentError("$nm length $(length(a)) != NES $nes"))
    end
    iz16 = zeros(Int32, 16); aw16 = zeros(Float64, 16)
    for i in 1:min(length(iz), 16); iz16[i] = Int32(iz[i]); end
    for i in 1:min(length(aw), 16); aw16[i] = Float64(aw[i]); end

    ACENeutronTable(header,
                    Float64.(energy_grid), Float64.(total_xs),
                    Float64.(absorption_xs), Float64.(elastic_xs),
                    Float64.(heating_numbers), reactions,
                    angular_elastic, angular, iz16, aw16)
end

# =========================================================================
# Accessor helpers for an ACENeutronTable
# =========================================================================
"""
    nes(table::ACENeutronTable) -> Int

Number of energy grid points.
"""
nes(table::ACENeutronTable) = length(table.energy_grid)

"""
    ntr(table::ACENeutronTable) -> Int

Number of reactions (excluding elastic).
"""
ntr(table::ACENeutronTable) = length(table.reactions)

# Prefixed aliases to avoid name collision
ace_nes(table::ACENeutronTable) = nes(table)
ace_ntr(table::ACENeutronTable) = ntr(table)
