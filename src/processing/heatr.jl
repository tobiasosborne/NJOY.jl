# HEATR -- KERMA (Kinetic Energy Release in MAterial) and damage energy
#
# PROPOSAL B: Pure-function composable heating kernels.
# Each reaction type has its own pure function; Lindhard damage is standalone.
# Correspondence to NJOY2016 heatr.f90:
#   elastic_heating -> disbar (MT2), capture_heating -> capdam (MT>=102),
#   fission_heating -> nheat fission, lindhard_damage -> df() (2015-2052),
#   compute_kerma -> nheat + gheat driver

# ==========================================================================
# Result type
# ==========================================================================

"""
    KERMAResult

Per-energy KERMA factors [eV-barn] from HEATR processing.
"""
struct KERMAResult
    energies::Vector{Float64}
    total_kerma::Vector{Float64}        # MT301
    elastic_kerma::Vector{Float64}      # MT302
    capture_kerma::Vector{Float64}      # MT402
    fission_kerma::Vector{Float64}      # MT318
    inelastic_kerma::Vector{Float64}    # MT304
    damage_energy::Vector{Float64}      # MT444
end

# ==========================================================================
# Default displacement energies [eV]  --  heatr.f90 lines 329-363
# ==========================================================================

const DEFAULT_DISPLACEMENT_ENERGY = Dict{Int,Float64}(
    4=>31.0, 6=>31.0, 12=>25.0, 13=>27.0, 14=>25.0, 20=>40.0,
    22=>40.0, 23=>40.0, 24=>40.0, 25=>40.0, 26=>40.0, 27=>40.0,
    28=>40.0, 29=>40.0, 40=>40.0, 41=>40.0, 42=>60.0, 47=>60.0,
    73=>90.0, 74=>90.0, 79=>30.0, 82=>25.0,
)

"Displacement threshold energy [eV] for element Z (NJOY2016 defaults)."
displacement_energy(Z::Integer) = get(DEFAULT_DISPLACEMENT_ENERGY, Int(Z), 25.0)

# ==========================================================================
# Lindhard damage energy partition  --  heatr.f90 df() lines 2015-2052
# ==========================================================================

"""
    LindharParams

Precomputed Lindhard electronic-screening damage function constants.
`el_inv` = 1/E_L [1/eV], `fl` = electronic stopping factor.
"""
struct LindharParams
    el_inv::Float64
    fl::Float64
end

"""
    lindhard_params(Zr, Ar, Zl, Al) -> LindharParams

Precompute Lindhard constants for recoil atom (Zr, Ar) in lattice (Zl, Al).
Matches NJOY2016 df() initialisation (heatr.f90 lines 2036-2043).
"""
function lindhard_params(Zr::Real, Ar::Real, Zl::Real, Al::Real)
    Zr == 0 && return LindharParams(0.0, 0.0)
    zr23 = Zr^(2/3);  zl23 = Zl^(2/3)
    el = 30.724 * Zr * Zl * sqrt(zr23 + zl23) * (Ar + Al) / Al
    denom = (zr23 + zl23)^0.75 * Ar^1.5 * sqrt(Al)
    fl = 0.0793 * zr23 * sqrt(Zl) * (Ar + Al)^1.5 / denom
    return LindharParams(1.0 / el, fl)
end

"""
    lindhard_damage(E_recoil, params::LindharParams; E_d=0.0) -> Float64

Lindhard damage energy [eV] for recoil energy `E_recoil`.
Robinson form: E_dam = E/(1 + f_L*(3.4008*eps^{1/6} + 0.40244*eps^{3/4} + eps)).
Returns 0 if E_recoil < E_d (displacement threshold).
"""
function lindhard_damage(E_recoil::Real, params::LindharParams; E_d::Real=0.0)
    E_recoil <= 0 && return 0.0
    params.el_inv == 0.0 && return 0.0
    E_recoil < E_d && return 0.0
    eps = E_recoil * params.el_inv
    return E_recoil / (1.0 + params.fl * (3.4008 * eps^(1/6) +
                                           0.40244 * eps^(3/4) + eps))
end

"Convenience: compute Lindhard params on the fly."
function lindhard_damage(E_recoil::Real, Zr::Real, Ar::Real,
                         Zl::Real, Al::Real; E_d::Real=0.0)
    (Zr == 0 || Zl == 0) && return 0.0
    return lindhard_damage(E_recoil, lindhard_params(Zr, Ar, Zl, Al); E_d=E_d)
end

# ==========================================================================
# Elastic scattering  (MT2)
# ==========================================================================

"""
    elastic_heating(E, A) -> Float64

Average energy deposited per elastic scatter: h = 2EA/(A+1)^2.
Exact for isotropic CM (s-wave) scattering.
"""
elastic_heating(E::Real, A::Real) = 2.0 * E * A / (1.0 + A)^2

"""
    elastic_heating_aniso(E, A, mu_bar) -> Float64

Elastic heating with CM anisotropy: h = 2EA(1 - mu_bar)/(A+1)^2.
Recovers `elastic_heating` when mu_bar=0.
"""
elastic_heating_aniso(E::Real, A::Real, mu_bar::Real) =
    2.0 * E * A * (1.0 - mu_bar) / (1.0 + A)^2

"""
    elastic_damage(E, A, Z; E_d=25.0) -> Float64

Damage energy from elastic scattering: Lindhard partition of average recoil.
Isotropic approximation — use `elastic_damage_angular` for full accuracy.
"""
function elastic_damage(E::Real, A::Real, Z::Real; E_d::Real=25.0)
    E_recoil = elastic_heating(E, A)
    lp = lindhard_params(Z, A, Z, A)
    return lindhard_damage(E_recoil, lp; E_d=E_d)
end

# 64-point Gauss-Legendre quadrature nodes and weights (heatr.f90 lines 1845-1890)
const _GL64_NODES = Float64[
    -9.99305042e-01, -9.96340117e-01, -9.91013371e-01, -9.83336254e-01,
    -9.73326828e-01, -9.61008800e-01, -9.46411375e-01, -9.29569172e-01,
    -9.10522137e-01, -8.89315446e-01, -8.65999398e-01, -8.40629296e-01,
    -8.13265315e-01, -7.83972359e-01, -7.52819907e-01, -7.19881850e-01,
    -6.85236313e-01, -6.48965471e-01, -6.11155355e-01, -5.71895646e-01,
    -5.31279464e-01, -4.89403146e-01, -4.46366017e-01, -4.02270158e-01,
    -3.57220158e-01, -3.11322872e-01, -2.64687162e-01, -2.17423644e-01,
    -1.69644420e-01, -1.21462819e-01, -7.29931218e-02, -2.43502927e-02,
     2.43502927e-02,  7.29931218e-02,  1.21462819e-01,  1.69644420e-01,
     2.17423644e-01,  2.64687162e-01,  3.11322872e-01,  3.57220158e-01,
     4.02270158e-01,  4.46366017e-01,  4.89403146e-01,  5.31279464e-01,
     5.71895646e-01,  6.11155355e-01,  6.48965471e-01,  6.85236313e-01,
     7.19881850e-01,  7.52819907e-01,  7.83972359e-01,  8.13265315e-01,
     8.40629296e-01,  8.65999398e-01,  8.89315446e-01,  9.10522137e-01,
     9.29569172e-01,  9.46411375e-01,  9.61008800e-01,  9.73326828e-01,
     9.83336254e-01,  9.91013371e-01,  9.96340117e-01,  9.99305042e-01]
const _GL64_WEIGHTS = Float64[
     1.78328072e-03,  4.14703326e-03,  6.50445797e-03,  8.84675983e-03,
     1.11681395e-02,  1.34630479e-02,  1.57260305e-02,  1.79517158e-02,
     2.01348232e-02,  2.22701738e-02,  2.43527026e-02,  2.63774697e-02,
     2.83396726e-02,  3.02346571e-02,  3.20579284e-02,  3.38051618e-02,
     3.54722133e-02,  3.70551285e-02,  3.85501532e-02,  3.99537411e-02,
     4.12625632e-02,  4.24735151e-02,  4.35837245e-02,  4.45905582e-02,
     4.54916279e-02,  4.62847966e-02,  4.69681828e-02,  4.75401657e-02,
     4.79993886e-02,  4.83447622e-02,  4.85754674e-02,  4.86909570e-02,
     4.86909570e-02,  4.85754674e-02,  4.83447622e-02,  4.79993886e-02,
     4.75401657e-02,  4.69681828e-02,  4.62847966e-02,  4.54916279e-02,
     4.45905582e-02,  4.35837245e-02,  4.24735151e-02,  4.12625632e-02,
     3.99537411e-02,  3.85501532e-02,  3.70551285e-02,  3.54722133e-02,
     3.38051618e-02,  3.20579284e-02,  3.02346571e-02,  2.83396726e-02,
     2.63774697e-02,  2.43527026e-02,  2.22701738e-02,  2.01348232e-02,
     1.79517158e-02,  1.57260305e-02,  1.34630479e-02,  1.11681395e-02,
     8.84675983e-03,  6.50445797e-03,  4.14703326e-03,  1.78328072e-03]

"""
    _legndr(u, nld) -> Vector{Float64}

Evaluate Legendre polynomials P_0(u)..P_{nld-1}(u).
Matches Fortran legndr (heatr.f90).
"""
function _legndr(u::Real, nld::Int)
    p = Vector{Float64}(undef, nld)
    nld >= 1 && (p[1] = 1.0)
    nld >= 2 && (p[2] = u)
    for l in 3:nld
        p[l] = ((2l - 3) * u * p[l-1] - (l - 2) * p[l-2]) / (l - 1)
    end
    return p
end

"""
    elastic_damage_angular(E, A, Z, fl_coeffs; E_d=25.0) -> Float64

Damage energy from elastic scattering using 64-point Gauss-Legendre
angular integration over the Lindhard function.
Matches Fortran disbar (heatr.f90 lines 1989-2001).

`fl_coeffs` = Legendre coefficients f_l from MF4/MT=2 at this energy
(fl[1]=1 normalization, fl[2]=mu_bar, ...).
"""
function elastic_damage_angular(E::Real, A::Real, Z::Real,
                                 fl_coeffs::AbstractVector{<:Real};
                                 E_d::Real=25.0)
    lp = lindhard_params(Z, A, Z, A)
    afact = 1.0 / (A + 1.0)^2
    arat = 1.0 / A
    # For elastic: r=1, b=A, g=b*arat=1
    g = 1.0
    nld = length(fl_coeffs)
    dame = 0.0
    for iq in 1:64
        u = _GL64_NODES[iq]
        p = _legndr(u, nld)
        # Angular distribution: f = Σ (2l-1)*fl(l)*P_l(u)/2
        f = 0.0
        for il in 1:nld
            f += (2*il - 1) * fl_coeffs[il] * p[il] / 2
        end
        # Recoil energy at this angle
        e2 = E * (1.0 - 2.0 * g * u + g * g) * afact / arat
        dame += _GL64_WEIGHTS[iq] * f * lindhard_damage(e2, lp; E_d=E_d)
    end
    return max(dame, 0.0)
end

"""
    _disbar_damage_fl(E, A, Z, r, fl_coeffs; E_d=25.0) -> Float64

Damage energy via 64-point Gauss-Legendre angular integration with
full MF4 Legendre coefficients. Matches Fortran disbar exactly.
"""
function _disbar_damage_fl(E::Real, A::Real, Z::Real, r::Real,
                            fl_coeffs::AbstractVector{<:Real}; E_d::Real=25.0)
    lp = lindhard_params(Z, A, Z, A)
    afact = 1.0 / (A + 1.0)^2
    arat = 1.0 / A
    g = r  # for neutron: g = b*arat = r
    nld = length(fl_coeffs)
    dame = 0.0
    for iq in 1:64
        u = _GL64_NODES[iq]
        p = _legndr(u, nld)
        f = 0.0
        for il in 1:nld
            f += (2*il - 1) * fl_coeffs[il] * p[il] / 2
        end
        e2 = E * (1.0 - 2.0 * g * u + g * g) * afact / arat
        dame += _GL64_WEIGHTS[iq] * f * lindhard_damage(e2, lp; E_d=E_d)
    end
    return max(dame, 0.0)
end

"""
    _disbar_damage(E, A, Z, r; E_d=25.0) -> Float64

Damage energy via 64-point Gauss-Legendre angular integration.
Uses isotropic angular distribution when MF4 data is not available.
"""
function _disbar_damage(E::Real, A::Real, Z::Real, r::Real; E_d::Real=25.0)
    lp = lindhard_params(Z, A, Z, A)
    afact = 1.0 / (A + 1.0)^2
    arat = 1.0 / A
    g = r  # for neutron: g = b*arat = r*A*(1/A) = r
    dame = 0.0
    for iq in 1:64
        u = _GL64_NODES[iq]
        # Isotropic: f = 1/2
        f = 0.5
        e2 = E * (1.0 - 2.0 * g * u + g * g) * afact / arat
        dame += _GL64_WEIGHTS[iq] * f * lindhard_damage(e2, lp; E_d=E_d)
    end
    return max(dame, 0.0)
end

"""
    capdam_particle(E, Q, A, Z, mt; E_d=25.0) -> Float64

Damage energy for capture-followed-by-particle-emission reactions.
Matches Fortran capdam (heatr.f90 lines 1807-1818):
4-point Gauss-Legendre angular integration of Lindhard damage for
the residual nucleus recoil.
MT mapping: 103→(p,Z=1,A=1), 104→(d,Z=1,A=2), 105→(t,Z=1,A=3),
106→(He3,Z=2,A=3), 107→(α,Z=2,A=4).
"""
function capdam_particle(E::Real, Q::Real, A::Real, Z::Real, mt::Integer;
                          E_d::Real=25.0)
    # Emitted particle Z and A (heatr.f90 lines 1769-1780)
    zx = mt > 102 ? 1.0 : 0.0
    if mt == 106 || mt == 107; zx = 2.0; end
    ax = mt > 102 ? 1.0 : 0.0
    if mt == 104; ax = 2.0; end
    if mt == 105 || mt == 106; ax = 3.0; end
    if mt == 107; ax = 4.0; end
    ax == 0 && return 0.0

    third = 0.333333333  # Fortran truncated constant
    denom = 1.0 / (ax^third + A^third)  # screening length
    econ = 1.029e6  # Coulomb energy constant (heatr.f90 line 1758)
    aw1fac = 1.0 / (A + 1.0)

    ec = econ * zx * Z * denom  # Coulomb barrier
    ea = Q + A * E * aw1fac     # available CM energy
    ea < 0 && return 0.0
    et = (A + 1.0 - ax) * E * aw1fac  # target recoil kinetic energy
    if ea > ec * (1.0 + 1e-10)
        ea = ec  # clamp to Coulomb barrier
    end

    # Residual nucleus: Z_res = Z - zx, A_res = A+1-ax
    lp = lindhard_params(Z - zx, A + 1.0 - ax, Z, A)

    # 4-point Gauss-Legendre angular quadrature (heatr.f90 lines 1753-1756)
    qp4 = (-0.86114, -0.33998, 0.33998, 0.86114)
    qw4 = (0.34785, 0.65215, 0.65215, 0.34785)

    dame = 0.0
    for iq in 1:4
        er = (et - 2.0 * sqrt(et * ax * ea) * qp4[iq] + ax * ea) * aw1fac
        dame += qw4[iq] * lindhard_damage(er, lp; E_d=E_d) / 2.0
    end
    return max(dame, 0.0)
end

# ==========================================================================
# Capture (MT102 and other disappearance reactions)
# ==========================================================================

"""
    capture_heating(E, Q, A; E_gamma=0.0) -> Float64

Capture heating (energy balance): h = E*A/(A+1) + Q - E_gamma.
Set E_gamma=0 for local gamma deposition.
"""
capture_heating(E::Real, Q::Real, A::Real; E_gamma::Real=0.0) =
    E * A / (A + 1.0) + Q - E_gamma

"""
    photon_recoil_heating(E, A, gamma_data) -> Float64

Capture heating after nheat+gheat two-step (Fortran heatr approach).
nheat deposits (E+Q)×σ, gheat subtracts escaping gamma energy
(E+Q-E/(A+1))×σ, net = E/(A+1)×σ + photon recoil corrections.
`gamma_data` = Vector{Tuple{Float64,Float64}} of (E_gamma [eV], yield) pairs.
Recoil per gamma: E_γ²/(2·M_res·c²) where M_res = (A+1)·emc2.
See heatr.f90 gheat lines 5278-5296, disgam lines 5671-5687.
"""
function photon_recoil_heating(E::Real, A::Real,
                                gamma_data::Vector{Tuple{Float64,Float64}})
    emc2 = PhysicsConstants.amassn * PhysicsConstants.amu *
           PhysicsConstants.clight^2 / PhysicsConstants.ev
    M_res_c2 = (A + 1.0) * emc2
    h = E / (A + 1.0)  # compound nucleus recoil energy (heatr.f90 gheat line 5292)
    for (eg, yld) in gamma_data
        h += eg^2 / (2.0 * M_res_c2) * yld
    end
    return h
end

"""
    photon_recoil_damage(E, A, Z, gamma_data; E_d=25.0) -> Float64

Damage energy from capture with gamma cascade (Fortran nheat+gheat net).
nheat adds capdam damage, gheat subtracts it and adds photon recoil damage.
Net result: only photon recoil damages survive (capdam cancels).
See heatr.f90 gheat lines 5285-5289: dame = edam*x*y then dame -= damn*x.
"""
function photon_recoil_damage(E::Real, A::Real, Z::Real,
                               gamma_data::Vector{Tuple{Float64,Float64}};
                               E_d::Real=25.0)
    emc2 = PhysicsConstants.amassn * PhysicsConstants.amu *
           PhysicsConstants.clight^2 / PhysicsConstants.ev
    M_res_c2 = (A + 1.0) * emc2
    lp = lindhard_params(Z, A + 1.0, Z, A)
    d = 0.0
    for (eg, yld) in gamma_data
        er = eg^2 / (2.0 * M_res_c2)
        d += lindhard_damage(er, lp; E_d=E_d) * yld
    end
    # No neutron recoil term: nheat's capdam and gheat's subtraction cancel
    # (heatr.f90 gheat line 5289: dame = dame - damn*x)
    return d
end

"""
    capture_recoil(E, Q, A) -> Float64

Compound nucleus recoil after capture: E/(A+1) + (A*E/(A+1)+Q)^2/(2*M*c^2).
Includes neutron momentum transfer and photon recoil (capdam line 1804).
"""
function capture_recoil(E::Real, Q::Real, A::Real)
    aw1 = A + 1.0
    e_cm = E / aw1
    e_avail = A * E / aw1 + Q
    emc2 = PhysicsConstants.amassn * PhysicsConstants.amu *
           PhysicsConstants.clight^2 / PhysicsConstants.ev
    return e_cm + e_avail^2 / (2.0 * emc2 * aw1)
end

"""
    capture_damage(E, Q, A, Z; E_d=25.0) -> Float64

Damage from radiative capture: Lindhard partition of capture recoil.
"""
function capture_damage(E::Real, Q::Real, A::Real, Z::Real; E_d::Real=25.0)
    E_recoil = capture_recoil(E, Q, A)
    lp = lindhard_params(Z, A + 1.0, Z, A)
    return lindhard_damage(E_recoil, lp; E_d=E_d)
end

# ==========================================================================
# Fission (MT18/19/20/21/38)
# ==========================================================================

"""
    FissionQComponents

MT458 fission energy components [eV]: fragments, prompt/delayed neutrons,
prompt/delayed gammas, betas, neutrinos, pseudo-Q, and total.
"""
struct FissionQComponents
    E_fragments::Float64
    E_prompt_n::Float64
    E_prompt_gamma::Float64
    E_delayed_beta::Float64
    E_delayed_gamma::Float64
    E_delayed_n::Float64
    E_neutrino::Float64
    E_pseudo_Q::Float64
    E_total::Float64
end

"Fission heating (energy balance): h = Q_total - E_neutrino - E_delayed."
fission_heating(E::Real, Q_total::Real;
                E_neutrino::Real=0.0, E_delayed::Real=0.0) =
    Q_total - E_neutrino - E_delayed

"Fission heating from MT458 components: h = E_fragments + E_prompt_gamma."
fission_heating(E::Real, qc::FissionQComponents) =
    qc.E_fragments + qc.E_prompt_gamma

# ==========================================================================
# Inelastic scattering (MT51-91) and (n,xn) reactions
# ==========================================================================

"Inelastic heating: h = E + Q - E_secondary."
inelastic_heating(E::Real, Q::Real, A::Real; E_secondary::Real=0.0) =
    E + Q - E_secondary

"""
    discrete_inelastic_ebar(E, Q, A; mu_bar=0.0) -> Float64

Average secondary neutron energy for discrete two-body inelastic scattering.
Matches Fortran disbar (heatr.f90 lines 1975,1983-1985):
  thresh = (A+1)/A * |Q|
  r = sqrt(1 - thresh/E)
  b = r * sqrt(awr/arat) = r * A  (for neutron, awp=1, arat=1/A)
  cn = (1 + 2*b*wbar + b²) * afact  where afact = 1/(A+1)²
  ebar = E * cn
"""
function discrete_inelastic_ebar(E::Real, Q::Real, A::Real; mu_bar::Real=0.0)
    thresh = (A + 1.0) / A * abs(Q)
    E <= thresh && return E  # below threshold, no scatter
    r = sqrt(1.0 - thresh / E)
    b = r * A  # for neutron: b = r * sqrt(awr/arat) = r * A
    afact = 1.0 / (1.0 + A)^2
    return E * (1.0 + 2.0 * b * mu_bar + b * b) * afact
end

"""
    evaporation_ebar(E, u, theta) -> Float64

Average outgoing neutron energy for an evaporation spectrum (MF5, LF=9).
Matches Fortran anabar (heatr.f90 lines 2450-2462):
  b1 = (E-u)/θ
  ebar = θ * (2 - b1²*exp(-b1) / (1 - (b1+1)*exp(-b1)))
"""
function evaporation_ebar(E::Real, u::Real, theta::Real)
    E <= u && return 0.0
    theta == 0.0 && return 0.0
    b1 = (E - u) / theta
    if b1 >= 0.001
        b3 = exp(-b1)
        return theta * (2.0 - b1 * b1 * b3 / (1.0 - (b1 + 1.0) * b3))
    else
        return 4.0 * (E - u) / 3.0
    end
end

"""
    _sed(E, Ep, theta, u) -> Float64

Evaporation spectrum probability density. Matches Fortran sed (heatr.f90:2601-2620).
"""
function _sed(E::Real, Ep::Real, theta::Real, u::Real)
    xeu = (E - u) / theta
    if abs(xeu) >= 1e-4
        return Ep * exp(-Ep / theta) / (theta^2 * (1.0 - exp(-xeu) * (1.0 + xeu)))
    else
        return Ep * exp(-Ep / theta) / (theta^2 * xeu^2 * (1.0 - xeu) / 2.0)
    end
end

"""
    evaporation_damage(E, u, theta, A, Z; E_d=25.0) -> Float64

Damage energy for evaporation spectrum (MF5, LF=9) via adaptive integration.
Matches Fortran anadam (heatr.f90 lines 2508-2598): adaptive convergence stack
with 4-point Gauss-Legendre angular quadrature at each outgoing energy.
"""
function evaporation_damage(E::Real, u::Real, theta::Real, A::Real, Z::Real;
                            E_d::Real=25.0)
    E <= u && return 0.0
    theta == 0.0 && return 0.0
    (E - u) < 1e-7 * E && return 0.0
    (E - u) < 10.0 && return 0.0  # Fortran line 2539

    lp = lindhard_params(Z, A, Z, A)
    awfac = 1.0 / A
    emax_out = E - u

    qp4 = (-0.86114, -0.33998, 0.33998, 0.86114)
    qw4 = (0.34785, 0.65215, 0.65215, 0.34785)

    function _yval(ep)
        d = 0.0
        for iq in 1:4
            er = (E - 2.0 * sqrt(E * ep) * qp4[iq] + ep) * awfac
            d += qw4[iq] * lindhard_damage(er, lp; E_d=E_d) / 2.0
        end
        return d * _sed(E, ep, theta, u)
    end

    # Adaptive convergence stack (Fortran lines 2547-2596)
    imax = 10
    x = zeros(Float64, imax)
    y = zeros(Float64, imax)

    # Initialize 3 seed points (Fortran lines 2547-2561)
    x[3] = 1.0
    x[2] = theta < emax_out ? theta : emax_out / 2.0
    x[1] = emax_out
    for k in 1:3
        y[k] = _yval(x[k])
    end

    # Adaptive trapezoidal (Fortran lines 2562-2596)
    i = 3; j = 0; dame = 0.0
    xlast = 0.0; ylast = 0.0
    while i > 0
        iflag = false
        if i > 1 && i < imax
            xm = (x[i-1] + x[i]) / 2.0
            ym = (y[i-1] + y[i]) / 2.0
            yt = _yval(xm)
            if abs(yt - ym) > 0.05 * abs(yt) && yt != 0.0
                iflag = true
            end
        end
        if iflag
            # Subdivide: push midpoint onto stack
            i += 1
            x[i] = x[i-1]
            x[i-1] = xm
            y[i] = y[i-1]
            y[i-1] = yt
        else
            # Accept: accumulate trapezoidal contribution
            j += 1
            if j > 1
                dame += (x[i] - xlast) * (y[i] + ylast) / 2.0
            end
            xlast = x[i]
            ylast = y[i]
            i -= 1
        end
    end
    return max(dame, 0.0)
end

"(n,xn) heating: h = E + Q - n_out * E_secondary."
nxn_heating(E::Real, Q::Real, A::Real, n_out::Integer;
            E_secondary::Real=0.0) = E + Q - n_out * E_secondary

# ==========================================================================
# MF4 mu_bar interpolation
# ==========================================================================

"""
    read_mf4_mubar(filename, mat) -> (energies, mu_bar_cm)

Read MF4/MT=2 Legendre data from an ENDF tape and extract the first
Legendre coefficient f₁ = ⟨cos θ_CM⟩ = μ̄_CM at each incident energy.
Returns vectors suitable for linear interpolation.
Matches Fortran heatr disbar: wbar=fl(2) at heatr.f90 line 1983.
"""
function read_mf4_mubar(filename::AbstractString, mat::Integer)
    io = open(filename)
    found = find_section(io, 4, 2; target_mat=mat)
    if !found
        close(io)
        return (Float64[0.0, 2e7], Float64[0.0, 0.0])  # isotropic fallback
    end
    head = read_cont(io)
    lvt = Int(head.L1)  # transformation matrix flag
    cont = read_cont(io)
    nk = Int(cont.N1)   # transformation matrix size
    # Skip transformation matrix if present
    if lvt == 1 && nk > 0
        nlines = cld(nk, 6)
        for _ in 1:nlines
            readline(io)
        end
    end
    # Read TAB2 header
    tab2 = read_cont(io)
    ne = Int(tab2.N2)  # number of incident energies
    # Read interpolation line
    readline(io)
    # Read each LIST record: extract E and f_1
    energies = Vector{Float64}(undef, ne)
    mu_bar = Vector{Float64}(undef, ne)
    for ie in 1:ne
        list_head = read_cont(io)
        e_val = Float64(list_head.C2)
        nl = Int(list_head.N1)  # NPL = number of Legendre coefficients
        energies[ie] = e_val
        # Read NL values packed 6 per line
        nlines = cld(nl, 6)
        f1 = 0.0
        if nlines > 0
            first_line = readline(io)
            p = rpad(first_line, 80)
            f1 = parse_endf_float(p[1:11])  # first coefficient = f_1 = mu_bar
            for _ in 2:nlines
                readline(io)
            end
        end
        mu_bar[ie] = f1
    end
    close(io)
    return (energies, mu_bar)
end

"""
    read_mf4_legendre(filename, mat; mt=2) -> (energies, coeffs)

Read ALL Legendre coefficients from MF4/MT. Returns energy grid and
a vector of vectors (fl[1]=1.0 normalization prepended, fl[2]=f_1=mu_bar, ...).
For angular integration in disbar (heatr.f90 lines 1964-1968).
"""
function read_mf4_legendre(filename::AbstractString, mat::Integer; mt::Integer=2)
    io = open(filename)
    found = find_section(io, 4, mt; target_mat=mat)
    if !found
        close(io)
        # isotropic fallback
        return (Float64[0.0, 2e7], [Float64[1.0, 0.0], Float64[1.0, 0.0]])
    end
    head = read_cont(io)
    lvt = Int(head.L1)
    cont = read_cont(io)
    nk = Int(cont.N1)
    if lvt == 1 && nk > 0
        nlines = cld(nk, 6)
        for _ in 1:nlines; readline(io); end
    end
    tab2 = read_cont(io)
    ne = Int(tab2.N2)
    readline(io)  # interpolation line
    energies = Vector{Float64}(undef, ne)
    coeffs = Vector{Vector{Float64}}(undef, ne)
    for ie in 1:ne
        list_head = read_cont(io)
        e_val = Float64(list_head.C2)
        nl = Int(list_head.N1)
        energies[ie] = e_val
        # Read all NL Legendre coefficients
        vals = Float64[]
        nlines = cld(nl, 6)
        for iline in 1:nlines
            ln = readline(io)
            p = rpad(ln, 80)
            nvals_this = min(6, nl - (iline-1)*6)
            for j in 1:nvals_this
                push!(vals, parse_endf_float(p[(j-1)*11+1:j*11]))
            end
        end
        # Prepend fl[1]=1.0 normalization (Fortran fl(1)=1)
        coeffs[ie] = vcat(1.0, vals)
    end
    close(io)
    return (energies, coeffs)
end

"""
    _interp_legendre(data, E) -> Vector{Float64}

Linearly interpolate MF4 Legendre coefficients at energy E.
Matches Fortran hgtfle which interpolates between the bracketing
LIST records. Coefficients with different nld are zero-padded to
the larger size before interpolation.
"""
function _interp_legendre(data::Tuple{Vector{Float64},Vector{Vector{Float64}}}, E::Real)
    es, cs = data
    E <= es[1] && return cs[1]
    E >= es[end] && return cs[end]
    idx = searchsortedlast(es, E)
    idx < 1 && return cs[1]
    idx >= length(es) && return cs[end]
    # Linear interpolation between cs[idx] and cs[idx+1]
    f = (E - es[idx]) / (es[idx+1] - es[idx])
    lo = cs[idx]; hi = cs[idx+1]
    nmax = max(length(lo), length(hi))
    result = Vector{Float64}(undef, nmax)
    for i in 1:nmax
        vlo = i <= length(lo) ? lo[i] : 0.0
        vhi = i <= length(hi) ? hi[i] : 0.0
        result[i] = vlo + f * (vhi - vlo)
    end
    return result
end

"""
    _disbar_dame_eval(e, A, Z, mf4_data, lp; E_d) -> Float64

Evaluate elastic damage at a single energy using 64-pt Gauss-Legendre
with MF4 Legendre coefficients at that energy.
"""
function _disbar_dame_eval(e::Real, A::Real, Z::Real,
                            mf4_data::Tuple{Vector{Float64},Vector{Vector{Float64}}},
                            lp::LindharParams; E_d::Real=25.0)
    fl = _interp_legendre(mf4_data, e)
    nld = length(fl)
    afact = 1.0 / (A + 1.0)^2
    arat = 1.0 / A
    dame = 0.0
    for iq in 1:64
        u = _GL64_NODES[iq]
        p = _legndr(u, nld)
        f = 0.0
        for il in 1:nld
            f += (2*il - 1) * fl[il] * p[il] / 2
        end
        e2 = e * (1.0 - 2.0 * u + 1.0) * afact / arat  # g=1 for elastic
        dame += _GL64_WEIGHTS[iq] * f * lindhard_damage(e2, lp; E_d=E_d)
    end
    return max(dame, 0.0)
end

"""
    build_disbar_damage_vector(energies, A, Z, mf4_data; E_d, step, thresh) -> Vector{Float64}

Compute damage at each energy, matching Fortran disbar stepping state machine
(heatr.f90 lines 1948-2013). For elastic (thresh=0): r=1. For inelastic (thresh>0):
r=sqrt(1-thresh/E) above threshold, r=0 below. Evaluates at 1.1x-stepped bracket
points clamped to MF4 grid boundaries, linearly interpolates between brackets.
"""
function build_disbar_damage_vector(energies::AbstractVector{<:Real},
                                     A::Real, Z::Real,
                                     mf4_data::Tuple{Vector{Float64},Vector{Vector{Float64}}};
                                     E_d::Real=25.0, step::Real=1.1,
                                     thresh::Real=0.0)
    lp = lindhard_params(Z, A, Z, A)
    afact = 1.0 / (A + 1.0)^2
    arat = 1.0 / A
    enx = E_d * arat / (4.0 * afact)  # minimum energy for nonzero damage
    small = 1e-10

    ne = length(energies)
    dame_out = zeros(Float64, ne)

    # State variables matching Fortran SAVE (lines 1895-1896)
    en = 0.0; damn = 0.0
    el = 0.0; daml = 0.0
    # Ebar state: cn = ebar/E fraction, interpolated same as dame (lines 1985,2007)
    cn_val = 0.0; cl_val = 0.0
    ebar_out = zeros(Float64, ne)

    # MF4 energy grid = discontinuity boundaries for hgtfle
    mf4_es = mf4_data[1]
    mf4_cs = mf4_data[2]
    mf4_idx = 1  # current position in MF4 grid

    # Hgtfle persistent state (Fortran SAVE at lines 4241-4242).
    # The SLIDE at label 120 only copies nhi values from fhi→flo,
    # leaving higher-order flo coefficients STALE from previous brackets.
    max_nld = length(mf4_cs) > 0 ? maximum(length(c) for c in mf4_cs) : 1
    hgt_flo = zeros(Float64, max_nld)
    hgt_fhi = zeros(Float64, max_nld)
    hgt_nlo = 0; hgt_nhi = 0
    hgt_elo = 0.0; hgt_ehi = 0.0
    hgt_pos = 0  # 0 = not yet initialized

    # Advance hgtfle state to bracket containing energy e, return fl with stale behavior
    function _hgtfle_get_fl(e)
        # Initialize on first call (matching hgtfle e=0 init, lines 4248-4318)
        if hgt_pos == 0 && length(mf4_es) >= 2
            c1 = mf4_cs[1]; hgt_nlo = length(c1)
            for i in 1:hgt_nlo; hgt_flo[i] = c1[i]; end
            hgt_elo = mf4_es[1]
            c2 = mf4_cs[2]; hgt_nhi = length(c2)
            for i in 1:hgt_nhi; hgt_fhi[i] = c2[i]; end
            for i in hgt_nhi+1:max_nld; hgt_fhi[i] = 0.0; end
            hgt_ehi = mf4_es[2]
            hgt_pos = 2
        end

        # Advance brackets: if e >= ehi, SLIDE forward (label 110→120→130)
        while e >= hgt_ehi * (1.0 - small) && hgt_pos < length(mf4_es)
            # SLIDE (label 120, lines 4333-4337): PARTIAL copy — only nhi values!
            for i in 1:hgt_nhi
                hgt_flo[i] = hgt_fhi[i]
            end
            # flo[nhi+1:end] retains STALE values from previous brackets
            hgt_nlo = hgt_nhi
            hgt_elo = hgt_ehi
            # Read next high bracket (lines 4350-4361)
            hgt_pos += 1
            cn = mf4_cs[hgt_pos]
            hgt_nhi = length(cn)
            for i in 1:hgt_nhi; hgt_fhi[i] = cn[i]; end
            for i in hgt_nhi+1:max_nld; hgt_fhi[i] = 0.0; end
            hgt_ehi = mf4_es[hgt_pos]
        end

        # Interpolate (label 130, lines 4368-4386)
        nlmax = max(hgt_nlo, hgt_nhi)
        fl = Vector{Float64}(undef, nlmax)
        if e < hgt_elo * (1.0 - small)
            # Below first point: isotropic (label 140)
            fl[1] = 1.0
            for i in 2:nlmax; fl[i] = 0.0; end
            return fl
        end
        if hgt_ehi > hgt_elo
            ff = (e - hgt_elo) / (hgt_ehi - hgt_elo)
            for i in 1:nlmax
                fl[i] = hgt_flo[i] + ff * (hgt_fhi[i] - hgt_flo[i])
            end
        else
            for i in 1:nlmax; fl[i] = hgt_flo[i]; end
        end
        return fl
    end

    for ie in 1:ne
        ee = energies[ie]
        # If ee within current bracket, just interpolate (line 1949)
        if ee <= en * (1.0 + small) && en > 0
            if en > el
                f = (ee - el) / (en - el)
                dame_out[ie] = daml + f * (damn - daml)
                # Ebar: terp1(el,cl,en,cn,ee) → ce, then ebar = ee*ce (line 2007-2008)
                ce = cl_val + f * (cn_val - cl_val)
                ebar_out[ie] = ee * ce
            else
                dame_out[ie] = damn
                ebar_out[ie] = ee * cn_val
            end
            continue
        end

        # Need new bracket point (lines 1950-1955)
        el = en
        daml = damn
        cl_val = cn_val
        e = step * el

        # Find next MF4 discontinuity energy (enext from hgtfle)
        while mf4_idx < length(mf4_es) && mf4_es[mf4_idx] < el * (1.0 + small)
            mf4_idx += 1
        end
        enext = mf4_idx <= length(mf4_es) ? mf4_es[mf4_idx] : 1e20

        # Clamp to MF4 boundary (line 1954)
        if enext < e * (1.0 - small)
            e = enext
        end
        # Jump to ee if beyond stepped range (line 1955)
        if ee > e * (1.0 + small)
            e = ee
        end
        if e < 1e-20
            e = ee  # first call: el=0
        end

        # Get fl coefficients with hgtfle stale behavior (needed for both dame and ebar)
        fl = _hgtfle_get_fl(e)

        # Evaluate damage at e (Fortran disbar lines 1956-2003)
        if e < enx * (1.0 - small)
            damn = 0.0
        else
            if thresh > 0
                # Inelastic: r from threshold (Fortran lines 1957-1963)
                if thresh >= e * (1.0 - small)
                    r = 0.0
                else
                    r = sqrt(1.0 - thresh / e)
                end
                damn = _disbar_damage_fl(e, A, Z, r, fl; E_d=E_d)
            else
                # Elastic: g=1, use same GL integration as _disbar_damage_fl with r=1
                damn = _disbar_damage_fl(e, A, Z, 1.0, fl; E_d=E_d)
            end
        end
        # Compute cn = ebar/E fraction (Fortran disbar lines 1975,1983-1985)
        r_cn = thresh > 0 ? (thresh >= e * (1.0 - small) ? 0.0 : sqrt(1.0 - thresh / e)) : 1.0
        b_cn = r_cn * A  # b = r * sqrt(awr/arat) = r * A for neutron
        wbar_cn = length(fl) >= 2 ? fl[2] : 0.0
        cn_val = (1.0 + 2.0 * b_cn * wbar_cn + b_cn * b_cn) * afact
        en = e

        # Interpolate to ee
        if en > el
            f = (ee - el) / (en - el)
            dame_out[ie] = daml + f * (damn - daml)
            ce = cl_val + f * (cn_val - cl_val)
            ebar_out[ie] = ee * ce
        else
            dame_out[ie] = damn
            ebar_out[ie] = ee * cn_val
        end
    end
    return dame_out, ebar_out
end

"""
    build_conbar_damage_vector(energies, u, theta, A, Z; E_d, step) -> (dame_out, ebar_out)

Precompute evaporation damage and ebar vectors using conbar-style bracket stepping.
Matches Fortran conbar (heatr.f90:2273-2290): evaluates anadam/anabar at 1.5x-stepped
bracket endpoints and linearly interpolates to each query energy.
"""
function build_conbar_damage_vector(energies::AbstractVector{<:Real},
                                    u::Real, theta::Real,
                                    A::Real, Z::Real;
                                    E_d::Real=25.0, step::Real=1.5)
    ne = length(energies)
    dame_out = zeros(Float64, ne)
    ebar_out = zeros(Float64, ne)
    small = 1e-10
    etp = 2e7  # etop

    # Conbar SAVE state (heatr.f90:2068-2071)
    e1 = 0.0;  d1 = 0.0;  eb1 = 0.0
    e2 = 0.0;  d2 = 0.0;  eb2 = 0.0

    for ie in 1:ne
        ee = energies[ie]
        ee <= u && continue

        # Check if we need a new bracket (Fortran conbar lines 2273-2289)
        if ee >= e2 * (1.0 - small)
            done = false
            while !done
                d1 = d2;  eb1 = eb2;  e1 = e2
                e2_new = step * e1
                if e2_new > etp * (1.0 + small)
                    e2_new = etp
                end
                if e1 == 0.0
                    e2_new = ee  # first call: jump to requested energy
                end
                e2 = e2_new
                d2 = evaporation_damage(e2, u, theta, A, Z; E_d=E_d)
                eb2 = evaporation_ebar(e2, u, theta)
                if abs(e2 - etp) < small * etp
                    done = true
                elseif e1 != 0.0 && ee <= e2 * (1.0 + small)
                    done = true
                end
            end
        end

        # Linear interpolation (Fortran conbar line 2290: terp1)
        if e2 > e1
            f = (ee - e1) / (e2 - e1)
            dame_out[ie] = d1 + f * (d2 - d1)
            ebar_out[ie] = eb1 + f * (eb2 - eb1)
        else
            dame_out[ie] = d2
            ebar_out[ie] = eb2
        end
    end
    return dame_out, ebar_out
end

"""
    build_capdam_damage_vector(energies, Q, A, Z, mt; E_d, step) -> dame_out

Precompute charged-particle damage vector using capdam-style bracket stepping.
Matches Fortran capdam (heatr.f90:1792-1826): evaluates GL damage at 1.1x-stepped
bracket endpoints and linearly interpolates to each query energy.
"""
function build_capdam_damage_vector(energies::AbstractVector{<:Real},
                                    Q::Real, A::Real, Z::Real, mt::Integer;
                                    E_d::Real=25.0, step::Real=1.1)
    ne = length(energies)
    dame_out = zeros(Float64, ne)
    small = 1e-10

    # Capdam SAVE state (heatr.f90:1763)
    en = 0.0; damn = 0.0
    el = 0.0; daml = 0.0

    for ie in 1:ne
        ee = energies[ie]

        # Check if we need a new bracket (Fortran capdam lines 1795-1821)
        if ee >= en * (1.0 - small)
            el = en; daml = damn
            e = step * el
            if e < ee * (1.0 - small)
                e = (1.0 - small) * ee
            end
            if el == 0.0
                e = ee  # first call
            end
            damn = capdam_particle(e, Q, A, Z, mt; E_d=E_d)
            en = e
        end

        # Linear interpolation (Fortran capdam line 1825: terp1)
        if en > el
            f = (ee - el) / (en - el)
            dame_out[ie] = daml + f * (damn - daml)
        else
            dame_out[ie] = damn
        end
    end
    return dame_out
end

"Linear interpolation of mu_bar from MF4 (energies, values) table."
function _interp_mubar(data::Tuple{Vector{Float64},Vector{Float64}}, E::Real)
    es, ms = data
    E <= es[1] && return ms[1]
    E >= es[end] && return ms[end]
    idx = searchsortedfirst(es, E)
    idx <= 1 && return ms[1]
    f = (E - es[idx-1]) / (es[idx] - es[idx-1])
    return ms[idx-1] + f * (ms[idx] - ms[idx-1])
end

# ==========================================================================
# Top-level KERMA driver
# ==========================================================================

"""
    compute_kerma(pendf; awr, Z, Q_values, E_d, local_gamma, mu_bar_data,
                  mf4_legendre_data, mf13_gamma) -> KERMAResult

Compute KERMA factors for all reactions in a PointwiseMaterial.
`awr`: atomic weight ratio, `Z`: atomic number (0=skip damage),
`Q_values`: Dict{Int,Float64} of MT->Q [eV], `E_d`: displacement energy [eV],
`mu_bar_data`: (energies, mu_bar) for elastic angular correction from MF4/MT=2,
`mf4_legendre_data`: (energies, coeffs) full Legendre data for angular damage integration,
`mf13_gamma`: Vector of (E_gamma, TabulatedFunction) for MF13 photon production XS.
  Subtracts escaping gamma energy from total KERMA (Fortran gheat line 5362: h=-y*ebar).
"""
function compute_kerma(pendf::PointwiseMaterial;
                       awr::Real=1.0, Z::Int=0,
                       Q_values::Dict{Int,Float64}=Dict{Int,Float64}(),
                       qm_values::Dict{Int,Float64}=Dict{Int,Float64}(),
                       lr_values::Dict{Int,Int32}=Dict{Int,Int32}(),
                       E_d::Union{Nothing,Float64}=nothing,
                       local_gamma::Bool=false,
                       gamma_data::Dict{Int,Vector{Tuple{Float64,Float64}}}=
                           Dict{Int,Vector{Tuple{Float64,Float64}}}(),
                       mu_bar_data::Union{Nothing,Tuple{Vector{Float64},Vector{Float64}}}=nothing,
                       mf4_legendre_data::Union{Nothing,Tuple{Vector{Float64},Vector{Vector{Float64}}}}=nothing,
                       conbar_params::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                       mf4_mubar_all::Union{Nothing,Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}}=nothing,
                       mf4_legendre_all::Union{Nothing,Dict{Int,Tuple{Vector{Float64},Vector{Vector{Float64}}}}}=nothing,
                       mf13_gamma::Vector{Tuple{Float64,TabulatedFunction}}=
                           Tuple{Float64,TabulatedFunction}[])
    ne = length(pendf.energies)
    ed = isnothing(E_d) ? displacement_energy(Z) : E_d

    # Precompute elastic damage+ebar vector matching Fortran disbar stepping scheme
    _el_dame_vec, _el_ebar_vec = if Z > 0 && mf4_legendre_data !== nothing
        build_disbar_damage_vector(pendf.energies, awr, Float64(Z), mf4_legendre_data;
                                    E_d=ed)
    else
        nothing, nothing
    end

    # Precompute inelastic damage+ebar vectors — Fortran disbar uses the SAME stepping
    # state machine for ALL discrete MTs (51-90), not direct evaluation.
    # Each MT gets its own disbar instance with per-MT save variables.
    # The ebar uses stepped cn interpolation matching Fortran (lines 1985,2007-2008).
    # For MTs WITHOUT MF4 data: Fortran disbar still steps with isotropic coefficients
    # (fl(1)=1, fl(2)=0, imiss=1, enext=etop). Generate synthetic isotropic MF4 data.
    _inel_dame_vecs = Dict{Int, Vector{Float64}}()
    _inel_ebar_vecs = Dict{Int, Vector{Float64}}()
    # Synthetic isotropic MF4: 2 energies spanning full range, coefficients [1.0, 0.0]
    _iso_mf4 = (Float64[1e-5, 2e7], Vector{Float64}[[1.0, 0.0], [1.0, 0.0]])
    for (icol, mt) in enumerate(pendf.mt_list)
        (mt < 51 || mt > 90) && continue
        Q = get(Q_values, mt, 0.0)
        abs(Q) < 1e-20 && continue  # skip if no Q (not a real reaction)
        thresh_mt = (awr + 1.0) / awr * abs(Q)
        # Use actual MF4 data if available, otherwise isotropic
        mf4_data = if mf4_legendre_all !== nothing && haskey(mf4_legendre_all, mt)
            mf4_legendre_all[mt]
        else
            _iso_mf4
        end
        if Z > 0
            dame_v, ebar_v = build_disbar_damage_vector(
                pendf.energies, awr, Float64(Z), mf4_data;
                E_d=ed, thresh=thresh_mt)
            _inel_dame_vecs[mt] = dame_v
            _inel_ebar_vecs[mt] = ebar_v
        else
            # No damage, but still need ebar from stepping
            _, ebar_v = build_disbar_damage_vector(
                pendf.energies, awr, 0.0, mf4_data;
                E_d=ed, thresh=thresh_mt)
            _inel_ebar_vecs[mt] = ebar_v
        end
    end

    # Precompute conbar damage vector for MT=91 (evaporation)
    # Fortran conbar uses bracket stepping with step=1.5 for damage (anadam),
    # then linearly interpolates to query energies (lines 2273-2290).
    # Ebar is computed directly via anabar at each query energy (line 2267).
    _conbar_dame_vec = if Z > 0 && conbar_params !== nothing
        u_con, theta_con = conbar_params
        dame_v, _ = build_conbar_damage_vector(pendf.energies, u_con, theta_con,
                                               awr, Float64(Z); E_d=ed)
        dame_v
    else
        nothing
    end

    # Precompute capdam damage vectors for MT=103-107 (charged-particle)
    # Fortran capdam uses bracket stepping with step=1.1 (same pattern as disbar).
    _capdam_dame_vecs = Dict{Int, Vector{Float64}}()
    if Z > 0
        for (icol, mt) in enumerate(pendf.mt_list)
            (mt < 103 || mt > 107) && continue
            Q = get(Q_values, mt, 0.0)
            _capdam_dame_vecs[mt] = build_capdam_damage_vector(
                pendf.energies, Q, awr, Float64(Z), mt; E_d=ed)
        end
    end

    total   = zeros(Float64, ne)
    elastic = zeros(Float64, ne)
    cap     = zeros(Float64, ne)
    fiss    = zeros(Float64, ne)
    inelast = zeros(Float64, ne)
    damage  = zeros(Float64, ne)

    # Skip redundant sum MTs (Fortran heatr skips MT=1,3,4,101)
    skip_mts = Set([1, 3, 4, 101])

    for (icol, mt) in enumerate(pendf.mt_list)
        mt in skip_mts && continue
        Q = get(Q_values, mt, 0.0)
        gd = get(gamma_data, mt, Tuple{Float64,Float64}[])
        for ie in 1:ne
            E = pendf.energies[ie]
            sigma = pendf.cross_sections[ie, icol]
            sigma <= 0.0 && continue

            if mt == 2
                mu = if mu_bar_data !== nothing
                    _interp_mubar(mu_bar_data, E)
                else
                    0.0
                end
                h = elastic_heating_aniso(E, awr, mu) * sigma
                elastic[ie] += h
                if Z > 0
                    if _el_dame_vec !== nothing
                        damage[ie] += _el_dame_vec[ie] * sigma
                    elseif mf4_legendre_data !== nothing
                        fl = _interp_legendre(mf4_legendre_data, E)
                        damage[ie] += elastic_damage_angular(E, awr, Float64(Z), fl; E_d=ed) * sigma
                    else
                        damage[ie] += elastic_damage(E, awr, Float64(Z); E_d=ed) * sigma
                    end
                end
            elseif mt == 102
                if !isempty(gd)
                    # Photon recoil approach (matches Fortran gheat)
                    h = photon_recoil_heating(E, awr, gd) * sigma
                    if Z > 0
                        damage[ie] += photon_recoil_damage(E, awr, Float64(Z), gd; E_d=ed) * sigma
                    end
                elseif local_gamma
                    h = capture_heating(E, Q, awr) * sigma
                    if Z > 0
                        damage[ie] += capture_damage(E, Q, awr, Float64(Z); E_d=ed) * sigma
                    end
                else
                    # No gamma data, no local deposition: use single-photon recoil
                    h = capture_recoil(E, Q, awr) * sigma
                    if Z > 0
                        damage[ie] += capture_damage(E, Q, awr, Float64(Z); E_d=ed) * sigma
                    end
                end
                cap[ie] += h
            elseif mt == 18 || (19 <= mt <= 21) || mt == 38
                h = fission_heating(E, Q) * sigma
                fiss[ie] += h
            elseif 51 <= mt <= 90
                # Fortran nheat lines 1179-1191: icon=1, disbar for discrete levels
                # q0=0 default, q0=t (QM from MF3 HEAD) if LR≠0 and LR≠31
                lr = get(lr_values, mt, Int32(0))
                q0 = (lr != 0 && lr != 31) ? get(qm_values, mt, 0.0) : 0.0
                # Use stepped ebar from disbar state machine (matching Fortran cn
                # interpolation at 1.1x bracket endpoints, lines 1985,2007-2008)
                if haskey(_inel_ebar_vecs, mt)
                    ebar = _inel_ebar_vecs[mt][ie]
                else
                    # Fallback: direct evaluation for MTs without MF4 data
                    wbar = if mf4_mubar_all !== nothing
                        _interp_mubar(get(mf4_mubar_all, mt, (Float64[0.0,2e7], Float64[0.0,0.0])), E)
                    else
                        0.0
                    end
                    ebar = discrete_inelastic_ebar(E, Q, awr; mu_bar=wbar)
                end
                h = (E + q0 - ebar) * sigma
                inelast[ie] += h
                if Z > 0
                    # Fortran disbar: use precomputed stepping state machine
                    if haskey(_inel_dame_vecs, mt)
                        damage[ie] += _inel_dame_vecs[mt][ie] * sigma
                    else
                        # Fallback: direct evaluation (no MF4 data for this MT)
                        thresh_mt = (awr + 1.0) / awr * abs(Q)
                        if E > thresh_mt
                            r = sqrt(1.0 - thresh_mt / E)
                            damage[ie] += _disbar_damage(E, awr, Float64(Z), r; E_d=ed) * sigma
                        end
                    end
                end
            elseif mt == 91
                # Fortran nheat lines 1183,1196-1202: icon=2, conbar for continuum
                # conbar line 2267: anabar computes ebar DIRECTLY at query energy
                # conbar lines 2273-2290: anadam damage uses 1.5x bracket stepping
                lr = get(lr_values, mt, Int32(0))
                q0 = (lr != 0 && lr != 31) ? get(qm_values, mt, Q) : Q
                if conbar_params !== nothing
                    u_con, theta_con = conbar_params
                    ebar = evaporation_ebar(E, u_con, theta_con)
                    h = (E + q0 - ebar) * sigma
                    if Z > 0
                        if _conbar_dame_vec !== nothing
                            damage[ie] += _conbar_dame_vec[ie] * sigma
                        else
                            damage[ie] += evaporation_damage(E, u_con, theta_con,
                                                              awr, Float64(Z); E_d=ed) * sigma
                        end
                    end
                else
                    h = (E + Q) * sigma  # fallback without MF5 data
                end
                inelast[ie] += h
            elseif mt == 16 || mt == 17 || mt == 37
                n_out = mt == 16 ? 2 : (mt == 17 ? 3 : 4)
                h = nxn_heating(E, Q, awr, n_out) * sigma
                inelast[ie] += h
            elseif mt > 100
                # Fortran nheat label 170 (lines 1201-1206): q0=q, icon=0
                # h = (E + q0) * sigma = (E + Q) * sigma
                h = (E + Q) * sigma
                cap[ie] += h
                if Z > 0 && mt >= 103 && mt <= 107
                    if haskey(_capdam_dame_vecs, mt)
                        damage[ie] += _capdam_dame_vecs[mt][ie] * sigma
                    else
                        damage[ie] += capdam_particle(E, Q, awr, Float64(Z), mt; E_d=ed) * sigma
                    end
                end
            end
            total[ie] += h
        end
    end
    # Subtract escaping gamma energy from MF13 photon production sections
    # Fortran gheat line 5362: h = -y * ebar, then c(2) += h
    for (e_gamma, tab) in mf13_gamma
        for ie in 1:ne
            E = pendf.energies[ie]
            sigma_gamma = interpolate(tab, E)
            sigma_gamma <= 0.0 && continue
            total[ie] -= sigma_gamma * e_gamma
        end
    end
    return KERMAResult(copy(pendf.energies), total, elastic, cap,
                       fiss, inelast, damage)
end

# ==========================================================================
# Sum rule verification
# ==========================================================================

"Check total KERMA = sum of partials within rtol."
function verify_kerma_sum_rule(result::KERMAResult; rtol::Real=0.01)
    for i in eachindex(result.energies)
        parts = result.elastic_kerma[i] + result.capture_kerma[i] +
                result.fission_kerma[i] + result.inelastic_kerma[i]
        t = result.total_kerma[i]
        abs(t) < 1e-30 && continue
        abs(t - parts) / abs(t) > rtol && return false
    end
    return true
end
