# MIXR -- Abundance-weighted mixing of isotopic cross sections
#
# Produces elemental (or mixture) cross sections from isotopic PENDF data
# by linear combination on a union energy grid with lin-lin interpolation.
#
# Correspondence to NJOY2016 mixr.f90:
#   union grid walk (lines 293-313)  -> union_energy_grid (sorted merge)
#   gety + sig accumulation (300-303) -> interpolate_column + mix loop
#   mix_materials driver              -> full pipeline
#
# Design: all core functions are pure. The interpolation uses only arithmetic
# operations (addition, multiplication, division), making it AD-compatible.
# No I/O or global state.

# ==========================================================================
# Types
# ==========================================================================

"""
    MixComponent

One isotopic component for cross-section mixing: a PointwiseMaterial
paired with its atom fraction (abundance weight).
"""
struct MixComponent
    material::PointwiseMaterial
    fraction::Float64

    function MixComponent(material::PointwiseMaterial, fraction::Real)
        fraction >= 0.0 || throw(ArgumentError(
            "fraction must be non-negative, got $fraction"))
        new(material, Float64(fraction))
    end
end

# ==========================================================================
# Union energy grid  (mixr.f90 lines 293-313)
# ==========================================================================

"""
    union_energy_grid(grids::AbstractVector{<:AbstractVector{<:Real}}) -> Vector{Float64}

Merge multiple sorted energy grids into a single sorted, unique grid.
Pure function.
"""
function union_energy_grid(grids::AbstractVector{<:AbstractVector{<:Real}})
    total = sum(length, grids; init=0)
    merged = Vector{Float64}(undef, total)
    k = 0
    for g in grids
        for e in g
            k += 1
            merged[k] = Float64(e)
        end
    end
    sort!(merged)
    unique!(merged)
    merged
end

function union_energy_grid(components::Vector{MixComponent})
    union_energy_grid([c.material.energies for c in components])
end

# ==========================================================================
# Linear interpolation on a sorted grid  (AD-compatible)
# ==========================================================================

"""
    interpolate_column(e_grid, xs, e) -> Float64

Linear-linear interpolation of `xs` defined on sorted `e_grid` at point `e`.
Returns 0.0 outside the data range. Pure, uses only basic arithmetic.
"""
function interpolate_column(e_grid::AbstractVector{<:Real},
                            xs::AbstractVector{<:Real},
                            e::Real)
    n = length(e_grid)
    n == 0 && return zero(Float64)
    e < e_grid[1] && return zero(Float64)
    e > e_grid[n] && return zero(Float64)
    e == e_grid[1] && return Float64(xs[1])
    e == e_grid[n] && return Float64(xs[n])
    # Binary search
    lo, hi = 1, n
    @inbounds while hi - lo > 1
        mid = (lo + hi) >> 1
        e_grid[mid] <= e ? (lo = mid) : (hi = mid)
    end
    # Lin-lin interpolation
    @inbounds begin
        el, eh = Float64(e_grid[lo]), Float64(e_grid[hi])
        f = (Float64(e) - el) / (eh - el)
        Float64(xs[lo]) + f * (Float64(xs[hi]) - Float64(xs[lo]))
    end
end

# ==========================================================================
# Core mixing: Dict-based interface (pure, AD-compatible)
# ==========================================================================

"""
    MixInput

Lightweight description of one material to mix: energy grid, Dict of
MT -> cross section vectors, and abundance fraction.
"""
struct MixInput
    energies::Vector{Float64}
    reactions::Dict{Int, Vector{Float64}}
    fraction::Float64
end

"""
    mix_reactions(inputs::Vector{MixInput}; mt_filter=nothing) ->
        (energies, reactions::Dict{Int, Vector{Float64}})

Mix isotopic reactions onto a union grid by abundance-weighted summation.
Pure function. Returns the union energy grid and mixed cross sections.

If `mt_filter` is provided, only those MTs appear in the output.
"""
function mix_reactions(inputs::Vector{MixInput}; mt_filter=nothing)
    isempty(inputs) && error("mix_reactions: need at least one input")

    # Collect all MTs
    all_mts = Set{Int}()
    for inp in inputs
        union!(all_mts, keys(inp.reactions))
    end
    mts = sort!(collect(all_mts))
    if mt_filter !== nothing
        mt_set = Set{Int}(mt_filter)
        filter!(mt -> mt in mt_set, mts)
    end
    isempty(mts) && error("mix_reactions: no matching MTs")

    # Union grid
    ugrid = union_energy_grid([inp.energies for inp in inputs])
    ne = length(ugrid)

    # Accumulate
    result = Dict{Int, Vector{Float64}}()
    for mt in mts
        mixed = zeros(Float64, ne)
        for inp in inputs
            xs_vec = get(inp.reactions, mt, nothing)
            xs_vec === nothing && continue
            w = inp.fraction
            @inbounds for i in 1:ne
                mixed[i] += w * interpolate_column(inp.energies, xs_vec, ugrid[i])
            end
        end
        result[mt] = mixed
    end

    (ugrid, result)
end

# ==========================================================================
# PointwiseMaterial interface
# ==========================================================================

"""
    mix_materials(components::Vector{MixComponent};
                  mt_filter=nothing, mat_out::Integer=0) -> PointwiseMaterial

Construct a mixed PointwiseMaterial by abundance-weighted combination of
isotopic cross sections on a union energy grid with lin-lin interpolation.

Corresponds to the full NJOY2016 mixr module.
"""
function mix_materials(components::Vector{MixComponent};
                       mt_filter=nothing, mat_out::Integer=0)
    isempty(components) && error("mix_materials: need at least one component")

    # Convert to MixInput
    inputs = MixInput[]
    for c in components
        rxn = Dict{Int, Vector{Float64}}()
        for (j, mt) in enumerate(c.material.mt_list)
            rxn[mt] = c.material.cross_sections[:, j]
        end
        push!(inputs, MixInput(c.material.energies, rxn, c.fraction))
    end

    ugrid, mixed = mix_reactions(inputs; mt_filter=mt_filter)

    # Build PointwiseMaterial
    mt_list = sort!(collect(keys(mixed)))
    ne = length(ugrid)
    nmt = length(mt_list)
    xs_mat = zeros(ne, nmt)
    for (j, mt) in enumerate(mt_list)
        xs_mat[:, j] .= mixed[mt]
    end

    PointwiseMaterial(Int32(mat_out), ugrid, xs_mat, mt_list)
end
