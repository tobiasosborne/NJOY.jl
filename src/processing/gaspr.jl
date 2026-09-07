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
# Reaction channel table  (gaspr.f90:495-826)
# ==========================================================================
#
# Fortran does not carry a plain MT -> multiplicity table.  For every MF3
# reaction it tracks TWO things simultaneously (gaspr.f90:495, 507-818):
#
#   izr         the residual-nucleus ZA, initialised to `za + zain` and
#               decremented by the ZA of every particle the reaction emits
#   y203..y207  the count of p / d / t / 3He / alpha *explicitly* emitted
#
# and only afterwards converts a light residual nucleus into gas
# (gaspr.f90:820-825).  Keeping the two halves in one table is what makes the
# port faithful: a lone multiplicity table silently drops every reaction whose
# gas comes from the residual rather than from the named ejectile.

# mt => (ΔZA removed from the residual, (p, d, t, 3He, alpha) emitted)
# Transcribed one-for-one from the `if (mth.eq.N)` chain at gaspr.f90:507-818.
const _GASPR_MT_CHANNEL = Dict{Int,Tuple{Int,NTuple{5,Int}}}(
     11 => (1004, (0,1,0,0,0)),   16 => (   2, (0,0,0,0,0)),
     17 => (   3, (0,0,0,0,0)),   22 => (2005, (0,0,0,0,1)),
     23 => (6013, (0,0,0,0,3)),   24 => (2006, (0,0,0,0,1)),
     25 => (2007, (0,0,0,0,1)),   28 => (1002, (1,0,0,0,0)),
     29 => (4009, (0,0,0,0,2)),   30 => (4010, (0,0,0,0,2)),
     32 => (1003, (0,1,0,0,0)),   33 => (1004, (0,0,1,0,0)),
     34 => (2004, (0,0,0,1,0)),   35 => (5011, (0,1,0,0,2)),
     36 => (5012, (0,0,1,0,2)),   37 => (   4, (0,0,0,0,0)),
     41 => (1003, (1,0,0,0,0)),   42 => (1004, (1,0,0,0,0)),
     44 => (2003, (2,0,0,0,0)),   45 => (3006, (1,0,0,0,1)),
    103 => (1001, (1,0,0,0,0)),  104 => (1002, (0,1,0,0,0)),
    105 => (1003, (0,0,1,0,0)),  106 => (2003, (0,0,0,1,0)),
    107 => (2004, (0,0,0,0,1)),  108 => (4008, (0,0,0,0,2)),
    109 => (6012, (0,0,0,0,3)),  111 => (2002, (2,0,0,0,0)),
    112 => (3005, (1,0,0,0,1)),  113 => (5011, (0,0,1,0,2)),
    114 => (5010, (0,1,0,0,2)),  115 => (2003, (1,1,0,0,0)),
    116 => (2004, (1,0,1,0,0)),  117 => (3006, (0,1,0,0,1)),
    152 => (   5, (0,0,0,0,0)),  153 => (   6, (0,0,0,0,0)),
    154 => (1005, (0,0,1,0,0)),  155 => (3007, (0,0,1,0,1)),
    156 => (1005, (1,0,0,0,0)),  157 => (1005, (0,1,0,0,0)),
    158 => (3007, (0,1,0,0,1)),  159 => (3007, (1,0,0,0,1)),
    160 => (   7, (0,0,0,0,0)),  161 => (   8, (0,0,0,0,0)),
    162 => (1006, (1,0,0,0,0)),  163 => (1007, (1,0,0,0,0)),
    164 => (1008, (1,0,0,0,0)),  165 => (2008, (0,0,0,0,1)),
    166 => (2009, (0,0,0,0,1)),  167 => (2010, (0,0,0,0,1)),
    168 => (2011, (0,0,0,0,1)),  169 => (1006, (0,1,0,0,0)),
    170 => (1007, (0,1,0,0,0)),  171 => (1008, (0,1,0,0,0)),
    172 => (1006, (0,0,1,0,0)),  173 => (1007, (0,0,1,0,0)),
    174 => (1008, (0,0,1,0,0)),  175 => (1009, (0,0,1,0,0)),
    176 => (2005, (0,0,0,1,0)),  177 => (2006, (0,0,0,1,0)),
    178 => (2007, (0,0,0,1,0)),  179 => (2005, (2,0,0,0,0)),
    180 => (4011, (0,0,0,0,2)),  181 => (3008, (1,0,0,0,1)),
    182 => (2005, (0,1,1,0,0)),  183 => (2004, (1,1,0,0,0)),
    184 => (2005, (1,0,1,0,0)),  185 => (2006, (0,1,1,0,0)),
    186 => (3005, (1,0,0,1,0)),  187 => (3006, (0,1,0,1,0)),
    188 => (3007, (0,0,1,1,0)),  189 => (3008, (0,0,1,0,1)),
    190 => (2004, (2,0,0,0,0)),  191 => (3004, (1,0,0,1,0)),
    192 => (3005, (0,1,0,1,0)),  193 => (4007, (0,0,0,1,1)),
    194 => (2006, (2,0,0,0,0)),  195 => (4012, (0,0,0,0,2)),
    196 => (3009, (1,0,0,0,1)),  197 => (3003, (3,0,0,0,0)),
    198 => (3004, (3,0,0,0,0)),  199 => (4009, (2,0,0,0,1)),
    200 => (2007, (2,0,0,0,0)),
)

# lr => (ΔZA removed from the residual, (p, d, t, 3He, alpha) emitted)
# The breakup flag carried in the MF3 TAB1 L2 field of a discrete-level
# section.  Applied only for MT51-91, and only after the emitted neutron has
# already been removed from the residual.  Ref: gaspr.f90:565-611.
# LR=39/40 (`izr=izr` in the Fortran) deliberately change nothing.
const _GASPR_LR_CHANNEL = Dict{Int,Tuple{Int,NTuple{5,Int}}}(
    22 => (2004, (0,0,0,0,1)),   23 => (6012, (0,0,0,0,3)),
    24 => (2005, (0,0,0,0,1)),   25 => (2006, (0,0,0,0,1)),
    28 => (1001, (1,0,0,0,0)),   29 => (4008, (0,0,0,0,2)),
    30 => (4009, (0,0,0,0,2)),   32 => (1002, (0,1,0,0,0)),
    33 => (1003, (0,0,1,0,0)),   34 => (2003, (0,0,0,1,0)),
    35 => (5010, (0,1,0,0,2)),   36 => (5011, (0,0,1,0,2)),
    39 => (   0, (0,0,0,0,0)),   40 => (   0, (0,0,0,0,0)),
)

# Residual-nucleus ZA => extra gas it contributes (gaspr.f90:820-825):
#     if (izr.eq.1001) y203=y203+1   ...   if (izr.eq.4008) y207=y207+2
# The Fortran writes six independent `if`s, but the guards are mutually
# exclusive (izr is a single integer), so a lookup is equivalent.
# The Be-8 entry is the physically interesting one: Be-8 is unbound and
# breaks into two alphas, which is why B-10(n,t) — whose residual is Be-8 —
# produces two alphas on top of its explicit triton.
const _GASPR_RESIDUAL_GAS = Dict{Int,NTuple{5,Int}}(
    1001 => (1,0,0,0,0),   # proton
    1002 => (0,1,0,0,0),   # deuteron
    1003 => (0,0,1,0,0),   # triton
    2003 => (0,0,0,1,0),   # He-3
    2004 => (0,0,0,0,1),   # alpha
    4008 => (0,0,0,0,2),   # Be-8 -> 2 alpha
)

"""
    gas_multiplicity(mt::Integer) -> NTuple{5,Int}

(proton, deuteron, triton, He-3, alpha) particles *explicitly emitted* by
reaction `mt`, i.e. Fortran's `y203..y207` before the residual-nucleus rules
of gaspr.f90:820-825 are applied.

This is only half of a reaction's gas production. Use [`gas_channel`](@ref)
for the value GASPR actually accumulates.
"""
gas_multiplicity(mt::Integer) =
    get(_GASPR_MT_CHANNEL, Int(mt), (0, (0,0,0,0,0)))[2]

"""
    gaspr_skips_mt(mt::Integer; iverf::Integer=6) -> Bool

Whether GASPR's accumulation pass skips reaction `mt` outright.

Ref: njoy-reference/src/gaspr.f90:475-491. The discrete-level bands
(MT600-849 for `iverf>=6`, MT700-798 for `iverf<6`) are skipped so that the
MT103-107 totals are not double counted — for B-10 that is what keeps MT700
`(n,t0)` from being added on top of MT105.
"""
function gaspr_skips_mt(mt::Integer; iverf::Integer=6)
    mt = Int(mt)
    lvmin, lvmax = iverf >= 6 ? (600, 849) : (700, 798)
    mt > 200 && mt < lvmin && return true          # gaspr.f90:475
    (mt > lvmax || mt == 0) && return true         # gaspr.f90:476
    mt <= 4 && return true                         # gaspr.f90:477
    6 <= mt <= 10 && return true
    12 <= mt <= 15 && return true
    18 <= mt <= 21 && return true
    38 <= mt <= 40 && return true
    mt == 43 && return true
    46 <= mt <= 50 && return true
    92 <= mt <= 101 && return true
    lvmin <= mt <= lvmax && return true            # gaspr.f90:486-490
    mt in (152, 153, 160, 161) && return true      # gaspr.f90:491
    false
end

"""
    gas_residual_za(mt, lr, za, zain; iverf=6) -> Int

Residual-nucleus ZA left by reaction `mt`, i.e. Fortran's `izr`.

Starts at `za + zain` (gaspr.f90:495) and loses the ZA of every emitted
particle. For MT51-91 the emitted neutron is removed first (`izr=izr-1`,
gaspr.f90:565) and the breakup flag `lr` removes the rest.
"""
function gas_residual_za(mt::Integer, lr::Integer, za::Real, zain::Real;
                         iverf::Integer=6)
    mt = Int(mt)
    izr = round(Int, za + zain)
    gaspr_skips_mt(mt; iverf) && return izr
    if 51 <= mt <= 91
        izr -= 1
        izr -= get(_GASPR_LR_CHANNEL, Int(lr), (0, (0,0,0,0,0)))[1]
    else
        izr -= get(_GASPR_MT_CHANNEL, mt, (0, (0,0,0,0,0)))[1]
    end
    izr
end

"""
    gas_channel(mt, lr, za, zain; iverf=6) -> NTuple{5,Int}

Total (proton, deuteron, triton, He-3, alpha) yield GASPR accumulates for
reaction `mt`: the explicitly emitted particles plus whatever the residual
nucleus itself contributes.

`za` is the material ZA, `zain` the incident-particle ZA (`int(nsub/10)`,
so 1 for neutrons — gaspr.f90:97-102), and `lr` the MF3 TAB1 breakup flag.

Ref: njoy-reference/src/gaspr.f90:495-826.

Two behaviours here have no counterpart in a plain multiplicity table and
were the T45 defect:

  * MT51-91 gas comes entirely from `lr`, e.g. B-10's LR=28 levels emit a
    proton and its LR=35 levels a deuteron and two alphas;
  * a light residual nucleus is itself gas, e.g. B-10(n,t) leaves Be-8
    (`izr = 5010+1-1003 = 4008`), adding two alphas to MT207.

The MF6/MT5 energy-dependent-yield path (Fortran's `111` sentinel,
gaspr.f90:501-506) is not ported; MT5 contributes nothing here, as before.
"""
function gas_channel(mt::Integer, lr::Integer, za::Real, zain::Real;
                     iverf::Integer=6)
    mt = Int(mt)
    gaspr_skips_mt(mt; iverf) && return (0,0,0,0,0)
    izr = round(Int, za + zain)
    y = (0,0,0,0,0)
    if 51 <= mt <= 91
        izr -= 1                                    # gaspr.f90:565
        dz, y = get(_GASPR_LR_CHANNEL, Int(lr), (0, (0,0,0,0,0)))
        izr -= dz
    else
        dz, y = get(_GASPR_MT_CHANNEL, mt, (0, (0,0,0,0,0)))
        izr -= dz
    end
    r = get(_GASPR_RESIDUAL_GAS, izr, nothing)      # gaspr.f90:820-825
    r === nothing ? y : y .+ r
end

"""
    gas_threshold_candidate(mt, lr, za, zain; iverf=6) -> Bool

Whether reaction `mt` takes part in choosing the gas-production threshold
`thrg` (the low end of the output grid).

Ref: njoy-reference/src/gaspr.f90:421-424. Deliberately *not* the same test
as "does this reaction yield gas": the Fortran keeps a reaction whose
residual is merely light (`0 < izr <= 2004`) even when it produces nothing,
and force-includes the Be-8 residual. Reproduced verbatim rather than
simplified, because `thrg` is a minimum and a wrong predicate silently
shifts the first output energy.
"""
function gas_threshold_candidate(mt::Integer, lr::Integer, za::Real, zain::Real;
                                 iverf::Integer=6)
    gaspr_skips_mt(mt; iverf) && return false
    izr = gas_residual_za(mt, lr, za, zain; iverf)
    izg = any(!iszero, gas_channel(mt, lr, za, zain; iverf))
    izr == 4008 && (izg = true)                     # gaspr.f90:421
    izg || !(izr > 2004 || izr <= 0)                # gaspr.f90:422-423
end

"""
    gas_yield(mt::Integer) -> NTuple{5,Int}

Explicitly emitted gas particles for `mt`, zero for reactions GASPR skips.

Residual-nucleus contributions are *not* included — they need the material
ZA. Prefer [`gas_channel`](@ref); this remains for callers that only have an
MT to hand.
"""
gas_yield(mt::Integer) = gaspr_skips_mt(mt) ? (0,0,0,0,0) : gas_multiplicity(mt)

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
