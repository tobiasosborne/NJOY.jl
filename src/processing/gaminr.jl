# GAMINR -- Multigroup photon interaction cross sections from ENDF data.
#
# Correspondence to NJOY2016 gaminr.f90:
#   genggp     -> PHOTON_GROUP_STRUCTURES / get_photon_group_structure
#   gnwtf      -> photon_weight_function (iwt=1..3)
#   gtsig      -> (MF23 data passed directly or read via group_photon_xs)
#   gtff       -> photon_feed_function (MF27 form factors)
#   gpanel     -> (reuses group_integrate from groupr.jl)
#   gaminr     -> gaminr / group_photon_xs (top-level entry points)
#
# Reuses: group_integrate from groupr.jl; infrastructure from group_structures.jl.

# ============================================================================
# Photon group structures (igg=2..10)
# Energy bounds in eV (ascending), matching Fortran genggp (gaminr.f90:585-776)
# ============================================================================

"CSEWG 94-group photon energy structure [eV], 95 bounds."
const CSEWG_94 = NTuple{95,Float64}((
    5e3, 1e4, 1.5e4, 2e4, 3e4, 3.5e4, 4e4, 4.5e4,
    5.5e4, 6e4, 6.5e4, 7.5e4, 8e4, 9e4, 1e5, 1.2e5,
    1.4e5, 1.5e5, 1.6e5, 1.9e5, 2.2e5, 2.6e5, 3e5, 3.25e5,
    3.5e5, 3.75e5, 4e5, 4.25e5, 4.5e5, 5e5, 5.25e5, 5.5e5,
    5.75e5, 6e5, 6.25e5, 6.5e5, 6.75e5, 7e5, 7.5e5, 8e5,
    8.25e5, 8.65e5, 9e5, 1e6, 1.125e6, 1.2e6, 1.25e6, 1.33e6,
    1.42e6, 1.5e6, 1.6e6, 1.66e6, 1.75e6, 1.875e6, 2e6, 2.166e6,
    2.333e6, 2.5e6, 2.666e6, 2.833e6, 3e6, 3.166e6, 3.333e6, 3.5e6,
    3.65e6, 3.8e6, 3.9e6, 4e6, 4.2e6, 4.4e6, 4.5e6, 4.7e6,
    5e6, 5.2e6, 5.4e6, 5.5e6, 5.75e6, 6e6, 6.25e6, 6.5e6,
    6.75e6, 7e6, 7.25e6, 7.5e6, 7.75e6, 8e6, 8.5e6, 9e6,
    9.5e6, 1e7, 1.06e7, 1.1e7, 1.2e7, 1.4e7, 2e7,
))

"LANL 12-group photon energy structure [eV], 13 bounds."
const LANL_12_GAMMA = NTuple{13,Float64}((
    1e4, 1e5, 5e5, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 2e7,
))

"Steiner 21-group photon energy structure [eV] (ORNL-TM-2564), 22 bounds."
const STEINER_21 = NTuple{22,Float64}((
    1e4, 1e5, 2e5, 4e5, 1e6, 1.5e6, 2e6, 2.5e6, 3e6, 3.5e6, 4e6,
    4.5e6, 5e6, 5.5e6, 6e6, 6.5e6, 7e6, 7.5e6, 8e6, 1e7, 1.2e7, 1.4e7,
))

"Straker 22-group photon energy structure [eV], 23 bounds."
const STRAKER_22 = NTuple{23,Float64}((
    1e4, 3e4, 6e4, 1e5, 1.5e5, 3e5, 4.5e5, 6e5, 8e5, 1e6, 1.33e6,
    1.66e6, 2e6, 2.5e6, 3e6, 3.5e6, 4e6, 5e6, 6e6, 7e6, 8e6, 1e7, 1.4e7,
))

"LANL 48-group photon energy structure [eV], 49 bounds."
const LANL_48_GAMMA = NTuple{49,Float64}((
    1e3, 1e4, 2e4, 3e4, 4.5e4, 6e4, 8e4, 1e5, 1.5e5, 2e5, 3e5, 4e5,
    4.5e5, 5e5, 5.25e5, 6e5, 7e5, 8e5, 9e5, 1e6, 1.125e6, 1.2e6,
    1.33e6, 1.5e6, 1.66e6, 1.875e6, 2e6, 2.333e6, 2.5e6, 2.666e6,
    3e6, 3.5e6, 4e6, 4.5e6, 5e6, 5.5e6, 6e6, 6.5e6, 7e6, 7.5e6,
    8e6, 9e6, 1e7, 1.2e7, 1.4e7, 1.7e7, 2e7, 3e7, 5e7,
))

"LANL 24-group photon energy structure [eV], 25 bounds."
const LANL_24_GAMMA = NTuple{25,Float64}((
    1e4, 3e4, 6e4, 1e5, 2e5, 3e5, 5e5, 5.25e5, 7.5e5, 1e6, 1.33e6,
    1.66e6, 2e6, 2.5e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 1e7,
    1.2e7, 1.7e7, 3e7,
))

"VITAMIN-C 36-group photon energy structure [eV], 37 bounds."
const VITAMINC_36 = NTuple{37,Float64}((
    1e4, 2e4, 3e4, 4.5e4, 6e4, 7e4, 1e5, 1.5e5, 2e5, 3e5, 4e5,
    4.5e5, 5.1e5, 5.12e5, 6e5, 7e5, 8e5, 1e6, 1.33e6, 1.5e6, 1.66e6,
    2e6, 2.5e6, 3e6, 3.5e6, 4e6, 4.5e6, 5e6, 5.5e6, 6e6, 6.5e6,
    7e6, 7.5e6, 8e6, 1e7, 1.2e7, 1.4e7,
))

"VITAMIN-E 38-group photon energy structure [eV], 39 bounds."
const VITAMINE_38 = NTuple{39,Float64}((
    1e4, 2e4, 3e4, 4.5e4, 6e4, 7e4, 7.5e4, 1e5, 1.5e5, 2e5, 3e5,
    4e5, 4.5e5, 5.1e5, 5.12e5, 6e5, 7e5, 8e5, 1e6, 1.33e6, 1.5e6,
    1.66e6, 2e6, 2.5e6, 3e6, 3.5e6, 4e6, 4.5e6, 5e6, 5.5e6, 6e6,
    6.5e6, 7e6, 7.5e6, 8e6, 1e7, 1.2e7, 1.4e7, 2e7,
))

"VITAMIN-J 42-group photon energy structure [eV], 43 bounds."
const VITAMINJ_42_GAMMA = NTuple{43,Float64}((
    1e3, 1e4, 2e4, 3e4, 4.5e4, 6e4, 7e4, 7.5e4, 1e5, 1.5e5, 2e5,
    3e5, 4e5, 4.5e5, 5.1e5, 5.12e5, 6e5, 7e5, 8e5, 1e6, 1.33e6,
    1.34e6, 1.5e6, 1.66e6, 2e6, 2.5e6, 3e6, 3.5e6, 4e6, 4.5e6,
    5e6, 5.5e6, 6e6, 6.5e6, 7e6, 7.5e6, 8e6, 1e7, 1.2e7, 1.4e7,
    2e7, 3e7, 5e7,
))

"Named photon group structure identifier."
@enum PhotonGroupId begin
    IGG_CSEWG94=2; IGG_LANL12=3; IGG_STEINER21=4; IGG_STRAKER22=5
    IGG_LANL48=6; IGG_LANL24=7; IGG_VITAMINC36=8; IGG_VITAMINE38=9; IGG_VITAMINJ42=10
end

"""Return photon group boundary tuple for built-in structure `igg`."""
function get_photon_group_structure(igg::PhotonGroupId)
    igg == IGG_CSEWG94    && return CSEWG_94
    igg == IGG_LANL12     && return LANL_12_GAMMA
    igg == IGG_STEINER21  && return STEINER_21
    igg == IGG_STRAKER22  && return STRAKER_22
    igg == IGG_LANL48     && return LANL_48_GAMMA
    igg == IGG_LANL24     && return LANL_24_GAMMA
    igg == IGG_VITAMINC36 && return VITAMINC_36
    igg == IGG_VITAMINE38 && return VITAMINE_38
    igg == IGG_VITAMINJ42 && return VITAMINJ_42_GAMMA
    error("Unknown photon group structure: $igg")
end
get_photon_group_structure(igg::Integer) = get_photon_group_structure(PhotonGroupId(igg))

"""Standard photon interaction MT numbers (ENDF-6 format)."""
const PHOTON_REACTIONS = (
    (mf=23, mt=501, name="total"),
    (mf=23, mt=502, name="coherent"),
    (mf=23, mt=504, name="incoherent"),
    (mf=23, mt=516, name="pair_prod"),
    (mf=23, mt=522, name="photoelectric"),
)

# ============================================================================
# Weight functions (iwt=2..3), Fortran: gnwtf (gaminr.f90:778-823)
# ============================================================================

"Constant weight for photon group averaging (iwt=2)."
photon_weight_constant(E::Real) = 1.0

"1/E weight for photon group averaging (iwt=3)."
photon_weight_inv_e(E::Real) = Float64(E) > 0.0 ? 1.0 / Float64(E) : 0.0

"""Return a photon weight function by NJOY iwt number (2=constant, 3=1/E)."""
function get_photon_weight_function(iwt::Int)
    iwt == 2 && return photon_weight_constant
    iwt == 3 && return photon_weight_inv_e
    throw(ArgumentError("unsupported photon weight function iwt=$iwt (use 2 or 3)"))
end

# ============================================================================
# Result type
# ============================================================================

"""Multigroup photon cross sections collapsed to a photon group structure."""
struct PhotonMultiGroupXS
    group_bounds::Vector{Float64}
    mt_list::Vector{Int}
    mt_names::Vector{String}
    xs::Matrix{Float64}       # (n_groups, n_reactions)
    flux::Vector{Float64}
    heating::Vector{Float64}
end

Base.show(io::IO, p::PhotonMultiGroupXS) = print(io,
    "PhotonMultiGroupXS(", length(p.mt_list), " reactions, ", size(p.xs, 1), " groups)")

# ============================================================================
# Core: photon group averaging -- reuses group_integrate from groupr.jl
# Fortran: gpanel + main loop (gaminr.f90:300-430)
# ============================================================================

"""Average a photon cross section onto a group structure. Returns (avg, flux)."""
function photon_group_average(energies::AbstractVector{<:Real},
                               xs_values::AbstractVector{<:Real}, bounds;
                               weight_fn=photon_weight_constant)
    ne = length(energies)
    @assert ne == length(xs_values) && ne >= 2
    gb = collect(Float64, bounds); ng = length(gb) - 1
    wvals = [Float64(weight_fn(Float64(energies[i]))) for i in 1:ne]
    swvals = [Float64(xs_values[i]) * wvals[i] for i in 1:ne]
    num = group_integrate(energies, swvals, gb)
    den = group_integrate(energies, wvals, gb)
    avg = [den[g] > 0.0 ? num[g] / den[g] : 0.0 for g in 1:ng]
    avg, den
end

"""Compute heating KERMA: integral(sigma*E*w dE) / integral(w dE)."""
function photon_heating_kerma(energies::AbstractVector{<:Real},
                               xs_values::AbstractVector{<:Real}, bounds;
                               weight_fn=photon_weight_constant)
    ne = length(energies); gb = collect(Float64, bounds); ng = length(gb) - 1
    wvals = [Float64(weight_fn(Float64(energies[i]))) for i in 1:ne]
    sewvals = [Float64(xs_values[i]) * Float64(energies[i]) * wvals[i] for i in 1:ne]
    num = group_integrate(energies, sewvals, gb)
    den = group_integrate(energies, wvals, gb)
    [den[g] > 0.0 ? num[g] / den[g] : 0.0 for g in 1:ng]
end

# ============================================================================
# Top-level: gaminr (Fortran: subroutine gaminr, gaminr.f90:25-536)
# ============================================================================

"""
    gaminr(reactions, group_bounds; iwt=2, compute_heating=false) -> PhotonMultiGroupXS

Process photon cross sections into multigroup form from prepared reaction data.
`reactions` is a vector of named tuples `(mt, name, energies, xs)`.
"""
function gaminr(reactions::AbstractVector, group_bounds;
                iwt::Int=2, compute_heating::Bool=false)
    gb = _resolve_photon_bounds(group_bounds)
    ng = length(gb) - 1; weight_fn = get_photon_weight_function(iwt)
    mt_list = Int[]; mt_names = String[]
    xs_avg = Matrix{Float64}(undef, ng, length(reactions))
    flux = zeros(Float64, ng); heating = zeros(Float64, ng)
    for (r, rxn) in enumerate(reactions)
        push!(mt_list, rxn.mt); push!(mt_names, String(rxn.name))
        avg, fl = photon_group_average(rxn.energies, rxn.xs, gb; weight_fn)
        xs_avg[:, r] = avg
        r == 1 && (flux .= fl)
        if compute_heating && rxn.mt in (522, 602)
            heating .= photon_heating_kerma(rxn.energies, rxn.xs, gb; weight_fn)
        end
    end
    PhotonMultiGroupXS(collect(Float64, gb), mt_list, mt_names, xs_avg, flux, heating)
end

"""Dict interface: pass energies + Dict(MT => xs_vector)."""
function gaminr(energies::AbstractVector{<:Real},
                xs_dict::AbstractDict{<:Integer,<:AbstractVector{<:Real}},
                igg; iwt::Int=2, compute_heating::Bool=false)
    reactions = [(mt=Int(mt), name=_photon_mt_name(Int(mt)),
                  energies=collect(Float64, energies),
                  xs=collect(Float64, xs)) for (mt, xs) in sort(collect(xs_dict))]
    gaminr(reactions, igg; iwt=iwt, compute_heating=compute_heating)
end

# ============================================================================
# File-reading interface: group_photon_xs (reads MF23 from ENDF file)
# ============================================================================

"""
    group_photon_xs(endf_file, group_structure; weight_fn=inv_e_weight,
                    mts=[501,502,504,516]) -> PhotonMultiGroupXS

Read MF23 (photo-atomic cross sections) from an ENDF file and average onto
the specified photon group structure. This is the file-reading convenience
entry point corresponding to the full gaminr pipeline.
"""
function group_photon_xs(endf_file::AbstractString, group_structure;
                         weight_fn=photon_weight_inv_e,
                         mts::Vector{Int}=Int[501, 502, 504, 516])
    gb = _resolve_photon_bounds(group_structure)
    ng = length(gb) - 1; nmt = length(mts)
    mt_names = [_photon_mt_name(mt) for mt in mts]
    xs_out = zeros(Float64, ng, nmt); flux = zeros(Float64, ng)
    open(endf_file, "r") do io
        for (col, mt) in enumerate(mts)
            energies, values = _read_mf23_section(io, mt)
            isempty(energies) && continue
            avg, fl = photon_group_average(energies, values, gb; weight_fn)
            xs_out[:, col] = avg
            col == 1 && (flux .= fl)
        end
    end
    PhotonMultiGroupXS(collect(Float64, gb), copy(mts), mt_names,
                        xs_out, flux, zeros(Float64, ng))
end

# ============================================================================
# Helpers
# ============================================================================

function _resolve_photon_bounds(gb)
    gb isa Integer && return collect(Float64, get_photon_group_structure(gb))
    gb isa PhotonGroupId && return collect(Float64, get_photon_group_structure(gb))
    collect(Float64, gb)
end

"Look up standard photon reaction name by MT number."
function _photon_mt_name(mt::Int)
    for r in PHOTON_REACTIONS; r.mt == mt && return r.name; end
    "MT$mt"
end

"Read MF23 tabulated cross section for a given MT. Returns (energies, values)."
function _read_mf23_section(io::IO, mt::Integer)
    seekstart(io)
    find_section(io, 23, mt) || return (Float64[], Float64[])
    read_cont(io)
    tab = read_tab1(io)
    (collect(Float64, tab.x), collect(Float64, tab.y))
end
