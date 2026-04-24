# ACE format constants, header, and flat table data structures
#
# Proposer-B design: composable type hierarchy with strict format correctness.
#
# Architecture:
#   ACEHeader      -- header metadata (ZAID, AWR, temperature, etc.)
#   ACETable       -- flat NXS/JXS/XSS representation for writer/reader
#
# The structured types (ReactionXS, AngularBlock, ACENeutronTable) are
# defined in ace_neutron.jl.

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
    # Fortran NJOY acefc aceout uses `a10` for hz and hd — which stores the
    # string right-justified in the 10-char buffer (leading-space padded).
    # hk and hm are declared character(len=70/10) and stored as-assigned
    # (left-aligned with trailing spaces). Ref: acefc.f90:12946-12957
    #
    # Why: mismatching justification in hz breaks every downstream comparison
    # at byte level; T50 referenceTape34 line 1 starts `  2004.10a` (right-
    # justified, 2 leading spaces) not `2004.10a  ` (left-justified).
    _right10 = s -> (t = strip(s); length(t) >= 10 ? String(t[1:10]) : lpad(t, 10))
    hz = _right10(zaid)
    hd = _right10(date)
    hk = rpad(length(comment) > 70 ? comment[1:70] : comment, 70)[1:70]
    # hm: do NOT strip — caller is responsible for correct alignment. Fortran
    # NJOY assigns `"   mat%4d"` directly to a character(len=10) variable; the
    # three leading spaces ARE the alignment and must survive.
    hm_src = mat_string
    hm = length(hm_src) >= 10 ? String(hm_src[1:10]) : rpad(hm_src, 10)
    ACEHeader(hz, awr, temp_mev, hd, hk, hm)
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
    dot === nothing && throw(ArgumentError("no '.' in ZAID: '\$s'"))
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
