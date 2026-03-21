# GASPR -- Gas production cross sections (MT203-MT207)
#
# Pure-function implementation: computes gas production cross sections by
# summing partial reaction cross sections weighted by charged-particle
# multiplicities.
#
# Correspondence to NJOY2016 gaspr.f90:
#   The multiplicity if/elseif chain (lines 500-826) -> gas_multiplicity()
#   sigma accumulation loop -> accumulate_gas() (pure reduction)
#   compute_gas_production  -> driver wrapping PointwiseMaterial
#
# Design: every function is pure (no I/O, no mutation of inputs).
# All numeric work uses basic arithmetic, making it AD-compatible.

# ==========================================================================
# Result type
# ==========================================================================

"""
    GasProductionResult

Gas production cross sections MT203-MT207 on a shared energy grid.
"""
struct GasProductionResult
    energies::Vector{Float64}
    mt203::Vector{Float64}   # proton production
    mt204::Vector{Float64}   # deuteron production
    mt205::Vector{Float64}   # triton production
    mt206::Vector{Float64}   # He-3 production
    mt207::Vector{Float64}   # alpha production
end

# ==========================================================================
# Multiplicity table  (gaspr.f90 lines 500-826)
# ==========================================================================

"""
    gas_multiplicity(mt::Integer) -> NTuple{5,Int}

Return (proton, deuteron, triton, He-3, alpha) multiplicities for ENDF
reaction `mt`.  Pure function; integer arithmetic only.
"""
function gas_multiplicity(mt::Integer)
    # Two-body charged-particle reactions MT103-107
    mt == 103 && return (1,0,0,0,0)  # (n,p)
    mt == 104 && return (0,1,0,0,0)  # (n,d)
    mt == 105 && return (0,0,1,0,0)  # (n,t)
    mt == 106 && return (0,0,0,1,0)  # (n,3He)
    mt == 107 && return (0,0,0,0,1)  # (n,alpha)
    # Multi-particle MT108-117
    mt == 108 && return (0,0,0,0,2)  # (n,2alpha)
    mt == 109 && return (0,0,0,0,3)  # (n,3alpha)
    mt == 111 && return (2,0,0,0,0)  # (n,2p)
    mt == 112 && return (1,0,0,0,1)  # (n,p+alpha)
    mt == 113 && return (0,0,1,0,2)  # (n,t+2alpha)
    mt == 114 && return (0,1,0,0,2)  # (n,d+2alpha)
    mt == 115 && return (1,1,0,0,0)  # (n,p+d)
    mt == 116 && return (1,0,1,0,0)  # (n,p+t)
    mt == 117 && return (0,1,0,0,1)  # (n,d+alpha)
    # (n,n'+particle) reactions
    mt == 11  && return (0,1,0,0,0)  # (n,2nd)
    mt == 22  && return (0,0,0,0,1)  # (n,n+alpha)
    mt == 23  && return (0,0,0,0,3)  # (n,n+3alpha)
    mt == 24  && return (0,0,0,0,1)  # (n,2n+alpha)
    mt == 25  && return (0,0,0,0,1)  # (n,3n+alpha)
    mt == 28  && return (1,0,0,0,0)  # (n,n+p)
    mt == 29  && return (0,0,0,0,2)  # (n,n+2alpha)
    mt == 30  && return (0,0,0,0,2)  # (n,2n+2alpha)
    mt == 32  && return (0,1,0,0,0)  # (n,n+d)
    mt == 33  && return (0,0,1,0,0)  # (n,n+t)
    mt == 34  && return (0,0,0,1,0)  # (n,n+3He)
    mt == 35  && return (0,1,0,0,2)  # (n,n+d+2alpha)
    mt == 36  && return (0,0,1,0,2)  # (n,n+t+2alpha)
    mt == 41  && return (1,0,0,0,0)  # (n,2n+p)
    mt == 42  && return (1,0,0,0,0)  # (n,3n+p)
    mt == 44  && return (2,0,0,0,0)  # (n,n+2p)
    mt == 45  && return (1,0,0,0,1)  # (n,n+p+alpha)
    # Higher multi-particle MTs 154-200
    mt == 154 && return (0,0,1,0,0)
    mt == 155 && return (0,0,1,0,1)
    mt == 156 && return (1,0,0,0,0)
    mt == 157 && return (0,1,0,0,0)
    mt == 158 && return (0,1,0,0,1)
    mt == 159 && return (1,0,0,0,1)
    mt == 162 && return (1,0,0,0,0)
    mt == 163 && return (1,0,0,0,0)
    mt == 164 && return (1,0,0,0,0)
    mt == 165 && return (0,0,0,0,1)
    mt == 166 && return (0,0,0,0,1)
    mt == 167 && return (0,0,0,0,1)
    mt == 168 && return (0,0,0,0,1)
    mt == 169 && return (0,1,0,0,0)
    mt == 170 && return (0,1,0,0,0)
    mt == 171 && return (0,1,0,0,0)
    mt == 172 && return (0,0,1,0,0)
    mt == 173 && return (0,0,1,0,0)
    mt == 174 && return (0,0,1,0,0)
    mt == 175 && return (0,0,1,0,0)
    mt == 176 && return (0,0,0,1,0)
    mt == 177 && return (0,0,0,1,0)
    mt == 178 && return (0,0,0,1,0)
    mt == 179 && return (2,0,0,0,0)
    mt == 180 && return (0,0,0,0,2)
    mt == 181 && return (1,0,0,0,1)
    mt == 182 && return (0,1,1,0,0)
    mt == 183 && return (1,1,0,0,0)
    mt == 184 && return (1,0,1,0,0)
    mt == 185 && return (0,1,1,0,0)
    mt == 186 && return (1,0,0,1,0)
    mt == 187 && return (0,1,0,1,0)
    mt == 188 && return (0,0,1,1,0)
    mt == 189 && return (0,0,1,0,1)
    mt == 190 && return (2,0,0,0,0)
    mt == 191 && return (1,0,0,1,0)
    mt == 192 && return (0,1,0,1,0)
    mt == 193 && return (0,0,0,1,1)
    mt == 194 && return (2,0,0,0,0)
    mt == 195 && return (0,0,0,0,2)
    mt == 196 && return (1,0,0,0,1)
    mt == 197 && return (3,0,0,0,0)
    mt == 198 && return (3,0,0,0,0)
    mt == 199 && return (2,0,0,0,1)
    mt == 200 && return (2,0,0,0,0)
    return (0,0,0,0,0)
end

# MTs that never produce gas
const _GASPR_SKIP_MTS = Set{Int}([
    1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
    26, 27, 31, 37, 38, 39, 40, 43, 46, 47, 48, 49, 50,
    92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102,
    152, 153, 160, 161, 203, 204, 205, 206, 207,
])

"""
    gas_yield(mt::Integer) -> NTuple{5,Int}

Return gas particle yields for reaction `mt`, with skip-set short-circuit.
"""
function gas_yield(mt::Integer)
    mt in _GASPR_SKIP_MTS && return (0,0,0,0,0)
    gas_multiplicity(mt)
end

# ==========================================================================
# Core pure function: Dict-based interface
# ==========================================================================

"""
    accumulate_gas(energies, reactions) -> GasProductionResult

Compute gas production cross sections MT203-MT207 by summing partial reactions
weighted by particle multiplicities.

# Arguments
- `energies`: sorted energy grid (eV).
- `reactions::Dict{Int,<:AbstractVector}`: MT -> cross section on `energies`.

Pure function. AD-compatible (only uses addition and integer multiplication).
"""
function accumulate_gas(energies::AbstractVector{<:Real},
                        reactions::Dict{Int,<:AbstractVector{<:Real}})
    ne = length(energies)
    T = promote_type(Float64, eltype(valtype(reactions)))
    mt203 = zeros(T, ne)
    mt204 = zeros(T, ne)
    mt205 = zeros(T, ne)
    mt206 = zeros(T, ne)
    mt207 = zeros(T, ne)

    for (mt, xs) in reactions
        p, d, t, h3, a = gas_yield(mt)
        (p == 0 && d == 0 && t == 0 && h3 == 0 && a == 0) && continue
        length(xs) == ne || throw(DimensionMismatch(
            "MT$mt: $(length(xs)) points vs $ne grid points"))
        @inbounds for i in 1:ne
            v = xs[i]
            p != 0 && (mt203[i] += p * v)
            d != 0 && (mt204[i] += d * v)
            t != 0 && (mt205[i] += t * v)
            h3 != 0 && (mt206[i] += h3 * v)
            a != 0 && (mt207[i] += a * v)
        end
    end

    GasProductionResult(Float64.(energies), mt203, mt204, mt205, mt206, mt207)
end

"""
    gas_production(energies, reactions) -> GasProductionResult

High-level driver. Alias for `accumulate_gas`.
"""
gas_production(energies, reactions) = accumulate_gas(energies, reactions)

"""
    gas_production_dict(result::GasProductionResult) -> Dict{Int, Vector{Float64}}

Convert result to Dict{MT => xs}, omitting zero-valued channels.
Useful for composing with downstream modules expecting Dict input.
"""
function gas_production_dict(result::GasProductionResult)
    d = Dict{Int, Vector{Float64}}()
    any(!iszero, result.mt203) && (d[203] = result.mt203)
    any(!iszero, result.mt204) && (d[204] = result.mt204)
    any(!iszero, result.mt205) && (d[205] = result.mt205)
    any(!iszero, result.mt206) && (d[206] = result.mt206)
    any(!iszero, result.mt207) && (d[207] = result.mt207)
    d
end

# ==========================================================================
# PointwiseMaterial interface
# ==========================================================================

"""
    compute_gas_production(pendf::PointwiseMaterial) -> PointwiseMaterial

Add gas production MT203-MT207 to a PointwiseMaterial, replacing any existing.
"""
function compute_gas_production(pendf::PointwiseMaterial)
    ne = length(pendf.energies)
    nmt = length(pendf.mt_list)

    # Build Dict for accumulate_gas
    rxn = Dict{Int, Vector{Float64}}()
    for (j, mt) in enumerate(pendf.mt_list)
        rxn[mt] = pendf.cross_sections[:, j]
    end

    result = accumulate_gas(pendf.energies, rxn)
    gd = gas_production_dict(result)

    # Rebuild: keep non-gas columns, append gas columns
    gas_mts_set = Set([203, 204, 205, 206, 207])
    keep = [j for (j, mt) in enumerate(pendf.mt_list) if !(mt in gas_mts_set)]

    new_mts = [pendf.mt_list[j] for j in keep]
    gas_keys = sort!(collect(keys(gd)))
    append!(new_mts, gas_keys)

    new_xs = zeros(ne, length(new_mts))
    for (jn, jo) in enumerate(keep)
        new_xs[:, jn] .= @view pendf.cross_sections[:, jo]
    end
    offset = length(keep)
    for (jn, mt) in enumerate(gas_keys)
        new_xs[:, offset + jn] .= gd[mt]
    end

    PointwiseMaterial(pendf.mat, copy(pendf.energies), new_xs, new_mts)
end
