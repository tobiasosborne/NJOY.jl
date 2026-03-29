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
"""
function elastic_damage(E::Real, A::Real, Z::Real; E_d::Real=25.0)
    E_recoil = elastic_heating(E, A)
    lp = lindhard_params(Z, A, Z, A)
    return lindhard_damage(E_recoil, lp; E_d=E_d)
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
Matches Fortran disbar (heatr.f90 lines 1975-1985,2008) with Q0=0:
  thresh = (A+1)/A * |Q|
  r = sqrt(1 - thresh/E)
  b = r * A,  g = r
  ebar = E * (1 + 2*g*mu_bar + g²) / (A+1)²
"""
function discrete_inelastic_ebar(E::Real, Q::Real, A::Real; mu_bar::Real=0.0)
    thresh = (A + 1.0) / A * abs(Q)
    E <= thresh && return E  # below threshold, no scatter
    r = sqrt(1.0 - thresh / E)
    g = r  # for neutron emission: g = b * arat = r*A * (1/A) = r
    return E * (1.0 + 2.0 * g * mu_bar + g * g) / (1.0 + A)^2
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
    compute_kerma(pendf; awr, Z, Q_values, E_d, local_gamma, mu_bar_data, mf13_gamma) -> KERMAResult

Compute KERMA factors for all reactions in a PointwiseMaterial.
`awr`: atomic weight ratio, `Z`: atomic number (0=skip damage),
`Q_values`: Dict{Int,Float64} of MT->Q [eV], `E_d`: displacement energy [eV],
`mu_bar_data`: (energies, mu_bar) for elastic angular correction from MF4/MT=2,
`mf13_gamma`: Vector of (E_gamma, TabulatedFunction) for MF13 photon production XS.
  Subtracts escaping gamma energy from total KERMA (Fortran gheat line 5362: h=-y*ebar).
"""
function compute_kerma(pendf::PointwiseMaterial;
                       awr::Real=1.0, Z::Int=0,
                       Q_values::Dict{Int,Float64}=Dict{Int,Float64}(),
                       E_d::Union{Nothing,Float64}=nothing,
                       local_gamma::Bool=false,
                       gamma_data::Dict{Int,Vector{Tuple{Float64,Float64}}}=
                           Dict{Int,Vector{Tuple{Float64,Float64}}}(),
                       mu_bar_data::Union{Nothing,Tuple{Vector{Float64},Vector{Float64}}}=nothing,
                       mf13_gamma::Vector{Tuple{Float64,TabulatedFunction}}=
                           Tuple{Float64,TabulatedFunction}[])
    ne = length(pendf.energies)
    ed = isnothing(E_d) ? displacement_energy(Z) : E_d

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
                    damage[ie] += elastic_damage(E, awr, Float64(Z); E_d=ed) * sigma
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
            elseif 51 <= mt <= 91
                h = inelastic_heating(E, Q, awr) * sigma
                inelast[ie] += h
            elseif mt == 16 || mt == 17 || mt == 37
                n_out = mt == 16 ? 2 : (mt == 17 ? 3 : 4)
                h = nxn_heating(E, Q, awr, n_out) * sigma
                inelast[ie] += h
            elseif mt > 100
                # Other absorption: use capture_recoil or local deposit
                h = (!isempty(gd) ? photon_recoil_heating(E, awr, gd) :
                     (Q != 0 ? capture_recoil(E, Q, awr) : (E + Q))) * sigma
                cap[ie] += h
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
