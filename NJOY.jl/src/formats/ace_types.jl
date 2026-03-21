# ACE format type-safe data structures for MCNP continuous-energy tables
#
# Proposer-B design: composable type hierarchy with strict format correctness.
#
# Architecture:
#   ACEHeader      -- header metadata (ZAID, AWR, temperature, etc.)
#   ReactionXS     -- per-reaction cross section data
#   AngularDist    -- angular distribution for one energy point
#   AngularBlock   -- angular distributions for one reaction
#   ACENeutronTable -- high-level structured representation
#
# The writer (ace_writer.jl) serializes ACENeutronTable into the flat
# NXS/JXS/XSS arrays and writes Type 1 (ASCII) format matching NJOY2016
# aceout() exactly: format strings '(a10,f12.6,1x,1pe11.4,1x,a10)',
# '(a70,a10)', '4(i7,f11.0)', '8i9', '4a20' (i20 or 1pe20.11).

using Printf

# =========================================================================
# NXS index constants (1-based, matching acefc.f90 variable names)
# =========================================================================
const NXS_LEN2  = 1   # length of XSS array
const NXS_IZAID = 2   # 1000*Z + A
const NXS_NES   = 3   # number of energy grid points
const NXS_NTR   = 4   # number of reactions excluding elastic
const NXS_NR    = 5   # number of reactions with secondary neutrons
const NXS_NTRP  = 6   # number of photon production reactions
const NXS_NTYPE = 7   # particle type count (0 for neutron)
const NXS_NDNF  = 8   # number of delayed neutron families
const NXS_IS    = 9   # isomeric state number
const NXS_IZ    = 10  # atomic number Z
const NXS_IA    = 11  # mass number A

# =========================================================================
# JXS index constants (1-based, pointers into XSS)
# =========================================================================
const JXS_ESZ   = 1   # energy/sigma/heating block
const JXS_NU    = 2   # fission nu data
const JXS_MTR   = 3   # MT numbers for reactions
const JXS_LQR   = 4   # Q-values
const JXS_TYR   = 5   # reaction type / neutron yield
const JXS_LSIG  = 6   # cross section locators
const JXS_SIG   = 7   # cross section data
const JXS_LAND  = 8   # angular distribution locators
const JXS_AND   = 9   # angular distribution data
const JXS_LDLW  = 10  # energy distribution locators
const JXS_DLW   = 11  # energy distribution data
const JXS_GPD   = 12  # photon production data
const JXS_MTRP  = 13  # photon production MT array
const JXS_LSIGP = 14  # photon production XS locators
const JXS_SIGP  = 15  # photon production XS data
const JXS_LANDP = 16  # photon angular dist locators
const JXS_ANDP  = 17  # photon angular dist data
const JXS_LDLWP = 18  # photon energy dist locators
const JXS_DLWP  = 19  # photon energy dist data
const JXS_YP    = 20  # photon yield data
const JXS_FIS   = 21  # total fission cross section
const JXS_END   = 22  # end of basic data
const JXS_IURPT = 23  # unresolved probability tables
const JXS_NUD   = 24  # delayed nu data
const JXS_DNDAT = 25  # delayed neutron data
const JXS_LDND  = 26  # delayed neutron dist locators
const JXS_DND   = 27  # delayed neutron dist data
const JXS_PTYPE = 30  # particle type for production
const JXS_NTRO  = 31  # number of particle production reactions
const JXS_PLOCT = 32  # particle production locator

# ESZ sub-block offsets (multiplied by NES, 0-based)
const ESZ_ENERGY  = 0  # energies (MeV)
const ESZ_TOTAL   = 1  # total cross section (barns)
const ESZ_DISAP   = 2  # disappearance/absorption cross section
const ESZ_ELASTIC = 3  # elastic scattering cross section
const ESZ_HEATING = 4  # average heating number (MeV/collision)

# =========================================================================
# ACE Header
# =========================================================================
"""
    ACEHeader

Metadata for an ACE table: ZAID string, atomic weight ratio,
temperature (MeV), date, comment, and MAT string.
"""
struct ACEHeader
    hz::String      # ZAID, e.g. "92235.80c" (padded to 10 chars)
    aw0::Float64    # atomic weight ratio to neutron
    tz::Float64     # temperature in MeV (kT)
    hd::String      # date string (10 chars)
    hk::String      # descriptive comment (up to 70 chars)
    hm::String      # MAT id string (10 chars), e.g. "   mat9228"
end

"""
    ACEHeader(; zaid, awr, temp_mev, date="", comment="", mat_string="")

Construct an ACEHeader with automatic padding to MCNP field widths.
"""
function ACEHeader(; zaid::AbstractString,
                     awr::Float64,
                     temp_mev::Float64,
                     date::AbstractString = "",
                     comment::AbstractString = "",
                     mat_string::AbstractString = "")
    hz = rpad(strip(zaid), 10)[1:10]
    hd = rpad(strip(date), 10)[1:10]
    hk = rpad(length(comment) > 70 ? comment[1:70] : comment, 70)[1:70]
    hm = rpad(strip(mat_string), 10)[1:10]
    ACEHeader(hz, awr, temp_mev, hd, hk, hm)
end

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
# Flat ACE table (NXS/JXS/XSS representation for writer/reader)
# =========================================================================
"""
    ACETable

Flat ACE format data table for MCNP, holding the raw NXS/JXS/XSS arrays.

This is the serialization-level representation: the writer emits this directly
and the reader parses files into this form. Use `ACENeutronTable` for a
type-safe, structured representation.

# Fields
- `zaid::String`     -- ZAID identifier, e.g. "92235.80c"
- `awr::Float64`     -- atomic weight ratio (to neutron mass)
- `temp::Float64`    -- temperature in MeV (= T_kelvin * k_B / 1e6)
- `date::String`     -- processing date string (10 chars)
- `comment::String`  -- descriptive comment (up to 70 chars)
- `mat_id::String`   -- MAT identifier string (10 chars)
- `pairs_iz::NTuple{16,Int32}`  -- IZ values for (IZ,AW) pairs
- `pairs_aw::NTuple{16,Float64}` -- AW values for (IZ,AW) pairs
- `nxs::NTuple{16,Int32}` -- NXS control array
- `jxs::NTuple{32,Int32}` -- JXS locator array
- `xss::Vector{Float64}`  -- XSS data array
"""
struct ACETable
    zaid::String
    awr::Float64
    temp::Float64      # temperature in MeV
    date::String       # 10-char date string
    comment::String    # up to 70 chars
    mat_id::String     # 10-char mat id, e.g. "   mat 125"

    # (IZ, AW) pairs -- 16 pairs for component nuclide information
    pairs_iz::NTuple{16, Int32}
    pairs_aw::NTuple{16, Float64}

    # Control arrays
    nxs::NTuple{16, Int32}
    jxs::NTuple{32, Int32}

    # Data array
    xss::Vector{Float64}
end

# =========================================================================
# Accessor helpers for ACETable (flat representation)
# =========================================================================
"""
    nxs_length(table::ACETable) -> Int

Return the number of XSS entries (NXS(1) = LEN2).
"""
nxs_length(table::ACETable) = Int(table.nxs[NXS_LEN2])

"""
    nxs_nes(table::ACETable) -> Int

Return the number of energy grid points (NXS(3) = NES).
"""
nxs_nes(table::ACETable) = Int(table.nxs[NXS_NES])

"""
    nxs_ntr(table::ACETable) -> Int

Return the number of reactions excluding elastic (NXS(4) = NTR).
"""
nxs_ntr(table::ACETable) = Int(table.nxs[NXS_NTR])

"""
    esz_energies(table::ACETable) -> view

Return a view into the energy grid from the ESZ block.
"""
function esz_energies(table::ACETable)
    i0 = Int(table.jxs[JXS_ESZ])
    n = nxs_nes(table)
    return @view table.xss[i0 : i0 + n - 1]
end

"""
    esz_total(table::ACETable) -> view

Return a view into the total cross section from the ESZ block.
"""
function esz_total(table::ACETable)
    i0 = Int(table.jxs[JXS_ESZ])
    n = nxs_nes(table)
    return @view table.xss[i0 + n : i0 + 2*n - 1]
end

"""
    esz_elastic(table::ACETable) -> view

Return a view into the elastic cross section from the ESZ block.
"""
function esz_elastic(table::ACETable)
    i0 = Int(table.jxs[JXS_ESZ])
    n = nxs_nes(table)
    return @view table.xss[i0 + 3*n : i0 + 4*n - 1]
end

# =========================================================================
# ZAID utility functions
# =========================================================================
"""
    format_zaid(za::Integer, suffix::AbstractString) -> String

Format a ZAID string like "92235.80c" from ZA and suffix.
"""
format_zaid(za::Integer, suffix::AbstractString) = @sprintf("%d.%s", za, suffix)

"""
    format_zaid(z::Integer, a::Integer, suffix::AbstractString) -> String

Format a ZAID from Z, A, and suffix.
"""
format_zaid(z::Integer, a::Integer, suffix::AbstractString) =
    format_zaid(1000 * z + a, suffix)

"""
    parse_zaid(s::AbstractString) -> (za::Int, suffix::String)

Parse "92235.80c" into (92235, "80c").
"""
function parse_zaid(s::AbstractString)
    s = strip(s)
    dot = findfirst('.', s)
    dot === nothing && throw(ArgumentError("no '.' in ZAID: '$s'"))
    (parse(Int, s[1:dot-1]), String(s[dot+1:end]))
end

"""
    temp_to_mev(temp_kelvin::Real) -> Float64

Convert temperature in Kelvin to MeV (kT), matching NJOY: tz = tempd*bk/1e6.
"""
temp_to_mev(temp_kelvin::Real) =
    Float64(temp_kelvin * PhysicsConstants.bk * 1e-6)

"""
    mev_to_temp(kT_mev::Real) -> Float64

Convert MeV (kT) back to Kelvin.
"""
mev_to_temp(kT_mev::Real) =
    Float64(kT_mev / (PhysicsConstants.bk * 1e-6))

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
