# Weight/flux functions for GROUPR multigroup averaging
#
# Correspondence to NJOY2016 groupr.f90:
#   genwtf (4765-5013) -> get_weight_function (iwt dispatch)
#   getwtf (5015-5207) -> individual weight functions below
#
# All weight functions have signature w(E::Real) -> Float64.
# The combined Maxwellian+1/E+fission form (iwt=4) matches NJOY convention.

# ==========================================================================
# Basic analytical weight functions
# ==========================================================================

"""
    constant_weight(E) -> 1.0

Constant weight function (iwt=2). Produces unweighted averages.
"""
constant_weight(E::Real) = 1.0

"""
    inv_e_weight(E) -> 1/E

1/E slowing-down spectrum weight (iwt=3). Standard for epithermal averaging.
"""
inv_e_weight(E::Real) = 1.0 / E

# ==========================================================================
# Combined Maxwellian + 1/E + fission spectrum (iwt=4)
# ==========================================================================

"""
    maxwell_inv_e_fission(E; Eb=0.0253, Tb=0.0253, Ec=820.3e3, Tc=1.4e6)

Three-region weight function matching NJOY iwt=4 (getwtf).

Regions:
- E <= Eb: Maxwellian, w(E) = Ab * E * exp(-E/Tb)
- Eb < E < Ec: 1/E slowing-down
- E >= Ec: Fission spectrum, w(E) = Ac * sqrt(E) * exp(-E/Tc)

The normalization constants Ab, Ac are chosen so that the function
is continuous at the breakpoints Eb, Ec.

Default parameters use kT = 0.0253 eV (thermal), fission temperature 1.4 MeV.
"""
function maxwell_inv_e_fission(E::Real;
                               Eb::Real=0.0253, Tb::Real=0.0253,
                               Ec::Real=820.3e3, Tc::Real=1.4e6)
    Eb_f = Float64(Eb)
    Tb_f = Float64(Tb)
    Ec_f = Float64(Ec)
    Tc_f = Float64(Tc)
    Ef = Float64(E)

    if Ef <= Eb_f
        # Maxwellian: continuous with 1/Eb at E=Eb -> Ab*Eb*exp(-Eb/Tb) = 1/Eb
        Ab = 1.0 / (exp(-Eb_f / Tb_f) * Eb_f^2)
        return Ab * Ef * exp(-Ef / Tb_f)
    elseif Ef < Ec_f
        return 1.0 / Ef
    else
        # Fission: continuous with 1/Ec at E=Ec -> Ac*sqrt(Ec)*exp(-Ec/Tc) = 1/Ec
        Ac = 1.0 / (exp(-Ec_f / Tc_f) * Ec_f^1.5)
        return Ac * sqrt(Ef) * exp(-Ef / Tc_f)
    end
end

# ==========================================================================
# VITAMIN-E weight function (iwt=11)
# ==========================================================================

"""
    vitamin_e_weight(E; kT=0.0253)

VITAMIN-E weight function (ORNL-5505, iwt=11).

Five regions:
- E < 0.414 eV:    Maxwellian, C1*E*exp(-E/kT)
- 0.414 <= E < 2.12e6:  1/E
- 2.12e6 <= E < 1e7:    fission, C3*sqrt(E)*exp(-E/theta)
- 1e7 <= E < 1.252e7:   1/E plateau (C4/E)
- 1.252e7 <= E < 1.568e7: fusion peak Gaussian
- E >= 1.568e7:          1/E tail (C6/E)
"""
function vitamin_e_weight(E::Real; kT::Real=0.0253)
    Ef = Float64(E)
    kTf = Float64(kT)
    en1, en2, en3 = 0.414, 2.12e6, 1.0e7
    en4, en5 = 1.252e7, 1.568e7
    theta = 1.415e6
    ep = 1.407e7
    fusion = 2.5e4

    if Ef < en1
        cc = exp(en1 / kTf) / en1^2
        return cc * Ef * exp(-Ef / kTf)
    elseif Ef < en2
        return 1.0 / Ef
    elseif Ef < en3
        con3 = 1.44934e-9
        return con3 * sqrt(Ef) * exp(-Ef / theta)
    elseif Ef < en4
        con4 = 3.90797e-2
        return con4 / Ef
    elseif Ef < en5
        con5 = 2.64052e-5
        return con5 * exp(-5.0 * (sqrt(Ef) - sqrt(ep))^2 / fusion)
    else
        con6 = 6.76517e-2
        return con6 / Ef
    end
end

# ==========================================================================
# Thermal + 1/E + fission + fusion (iwt=6)
# ==========================================================================

"""
    thermal_fission_fusion(E; kT=0.054)

Combined thermal-1/E-fission+fusion weight (iwt=6).
Matches NJOY getwtf logic for iwtt=6.
"""
function thermal_fission_fusion(E::Real; kT::Real=0.054)
    Ef = Float64(E)
    kTf = Float64(kT)
    bb = 2.0 * kTf
    wt6b = 1.578551e-3
    wt6c = 2.1e6
    wt6e = 1.4e6
    wt6f = 2.5e4
    wt6g = 1.407e7

    if Ef <= bb
        cc = wt6b * exp(2.0) / bb^2
        return cc * Ef * exp(-Ef / kTf)
    elseif Ef <= wt6c
        return wt6b / Ef
    else
        wt6d = 2.32472e-12
        wt6h = 2.51697e-11
        w = wt6d * sqrt(Ef) * exp(-Ef / wt6e)
        pow = -(sqrt(Ef / wt6f) - sqrt(wt6g / wt6f))^2 / 2.0
        if pow > -89.0
            w += wt6h * exp(pow)
        end
        return w
    end
end

# ==========================================================================
# Tabulated weight function
# ==========================================================================

"""
    tabulated_weight(tab::TabulatedFunction)

Return a weight function that evaluates a user-supplied tabulated function
via interpolation. The returned callable has signature w(E) -> Float64.
"""
function tabulated_weight(tab::TabulatedFunction)
    return function(E::Real)
        Float64(interpolate(tab, Float64(E)))
    end
end

# ==========================================================================
# Weight function registry (iwt dispatch)
# ==========================================================================

"""
    get_weight_function(iwt::Int; kwargs...) -> callable

Return a weight function by NJOY iwt number. Supported iwt values:

| iwt | Description |
|-----|-------------|
|  2  | Constant |
|  3  | 1/E |
|  4  | Maxwellian + 1/E + fission spectrum |
|  6  | Thermal + 1/E + fission + fusion |
| 11  | VITAMIN-E (ORNL-5505) |

Keyword arguments are forwarded to the underlying function (e.g., kT, Eb, Ec).
"""
function get_weight_function(iwt::Int; kwargs...)
    if iwt == 2
        return constant_weight
    elseif iwt == 3
        return inv_e_weight
    elseif iwt == 4
        kw = Dict{Symbol,Any}(kwargs)
        Eb = get(kw, :Eb, 0.0253)
        Tb = get(kw, :Tb, 0.0253)
        Ec = get(kw, :Ec, 820.3e3)
        Tc = get(kw, :Tc, 1.4e6)
        return E -> maxwell_inv_e_fission(E; Eb=Eb, Tb=Tb, Ec=Ec, Tc=Tc)
    elseif iwt == 6
        kw = Dict{Symbol,Any}(kwargs)
        kT = get(kw, :kT, 0.054)
        return E -> thermal_fission_fusion(E; kT=kT)
    elseif iwt == 11
        kw = Dict{Symbol,Any}(kwargs)
        kT = get(kw, :kT, 0.0253)
        return E -> vitamin_e_weight(E; kT=kT)
    else
        throw(ArgumentError("unsupported weight function iwt=$iwt"))
    end
end
