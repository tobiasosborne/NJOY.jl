# THERMR -- Thermal scattering cross sections and kernels
# PROPOSAL B: Composable, AD-compatible pure functions.
# Free-gas and S(alpha,beta) share a generic kernel interface (E,E',mu).
# Detailed balance is exactly satisfied by construction.
# Correspondence: free_gas_kernel->sig()@200, sab_kernel->sig()@2514,
#   bragg_edges->sigcoh(), incoh_elastic_xs->iel(), compute_thermal->thermr()

# ---- Data types ----

"Tabulated S(alpha,beta) thermal scattering law (ENDF MF7/MT4). Stores log(S)."
struct SABData
    alpha::Vector{Float64}; beta::Vector{Float64}
    sab::Matrix{Float64}   # log S(alpha,beta), (n_alpha, n_beta)
    sigma_b::Float64; awr::Float64; T_eff::Float64
    lasym::Int; lat::Int
end

"Precomputed Bragg edge data for coherent elastic scattering."
struct BraggData
    tau_sq::Vector{Float64}; form_factor::Vector{Float64}
    econ::Float64; scon::Float64; n_edges::Int
end

"Output of compute_thermal: thermal cross sections on an energy grid."
struct ThermalResult
    energies::Vector{Float64}; inelastic_xs::Vector{Float64}
    elastic_xs::Vector{Float64}; model::Symbol
end

# ---- Free-gas model (sig() lines 2600-2608) ----

"""
    free_gas_kernel(E, E_prime, mu, A, T; sigma_b=1.0) -> Float64

Free-gas d^2sigma/(dE' dmu) [barn/eV]. S(a,b)=exp(-(a+b)^2/(4a))/sqrt(4*pi*a)
exactly satisfies detailed balance. Matches NJOY2016 sig()@200.
"""
function free_gas_kernel(E::Real, E_prime::Real, mu::Real, A::Real, T::Real;
                         sigma_b::Real=1.0)
    kT = PhysicsConstants.bk * T
    Ef, Epf = Float64(E), Float64(E_prime)
    (Ef <= 0 || Epf <= 0 || kT <= 0) && return 0.0
    alpha = max((Ef + Epf - 2.0*Float64(mu)*sqrt(Ef*Epf)) / (A*kT), 1e-10)
    beta = (Epf - Ef) / kT
    arg = (alpha + beta)^2 / (4.0*alpha)
    arg > 225.0 && return 0.0
    sigma_b * sqrt(Epf/Ef) / (2.0*kT) * exp(-arg) / sqrt(4.0*PhysicsConstants.pi*alpha)
end

"""
    free_gas_xs(E, A, T; sigma_b=1.0) -> Float64

Analytical free-gas total thermal XS [barn].
sigma = (sigma_b/A) * [(1+1/(2x^2))*erf(x) + exp(-x^2)/(x*sqrt(pi))],
where x = sqrt(A*E/kT).
Limits: E >> kT -> sigma_b/A; E << kT -> sigma_b*2/(A*x*sqrt(pi)) (1/v).
"""
function free_gas_xs(E::Real, A::Real, T::Real; sigma_b::Real=1.0)
    kT = PhysicsConstants.bk * T; x2 = A * Float64(E) / kT; x = sqrt(x2)
    x < 1.0e-6 && return sigma_b / A * 2.0 / (x * sqrt(PhysicsConstants.pi))
    isp = 1.0 / sqrt(PhysicsConstants.pi); ef = erf(x); ex = exp(-x2)
    sigma_b / A * ((1.0 + 1.0 / (2.0 * x2)) * ef + ex * isp / x)
end

# ---- S(alpha,beta) model (sig() lines 2514-2598) ----

"""
    read_thermal_data(alpha, beta, sab_values; sigma_b, awr, ...) -> SABData

Construct SABData. If `is_log=false`, raw S values are log-converted.
"""
function read_thermal_data(alpha::AbstractVector, beta::AbstractVector,
                           sab_values::AbstractMatrix;
                           sigma_b::Real, awr::Real, T_eff::Real=0.0,
                           lasym::Int=0, lat::Int=0, is_log::Bool=true)
    na, nb = length(alpha), length(beta)
    @assert size(sab_values) == (na, nb)
    sabflg = -225.0
    sl = is_log ? Float64.(sab_values) : [Float64(v) > exp(sabflg) ?
         log(Float64(v)) : sabflg for v in sab_values]
    sl = reshape(sl, na, nb)
    SABData(Float64.(alpha), Float64.(beta), sl,
            Float64(sigma_b), Float64(awr), Float64(T_eff), lasym, lat)
end

function _interp_sab(a::Real, b::Real, data::SABData)
    sabflg = -225.0
    al, be, sab = data.alpha, data.beta, data.sab
    (a <= 0 || a > al[end]) && return sabflg
    bv = data.lasym == 0 ? abs(b) : b
    (bv > be[end] || bv < be[1]) && return sabflg
    ia = clamp(searchsortedlast(al, a), 1, length(al)-1)
    ib = clamp(searchsortedlast(be, bv), 1, length(be)-1)
    fa = clamp((a - al[ia])/(al[ia+1] - al[ia]), 0.0, 1.0)
    fb = clamp((bv - be[ib])/(be[ib+1] - be[ib]), 0.0, 1.0)
    s00, s10, s01, s11 = sab[ia,ib], sab[ia+1,ib], sab[ia,ib+1], sab[ia+1,ib+1]
    any(v -> v <= sabflg, (s00,s10,s01,s11)) && return sabflg
    (1-fa)*(1-fb)*s00 + fa*(1-fb)*s10 + (1-fa)*fb*s01 + fa*fb*s11
end

"""
    sab_kernel(E, E_prime, mu, data::SABData, T) -> Float64

Differential kernel from tabulated S(a,b) [barn/eV]. Falls back to SCT
approximation outside table range. Detailed balance is exact.
"""
function sab_kernel(E::Real, E_prime::Real, mu::Real, data::SABData, T::Real;
                    sigma_b_override::Union{Real,Nothing}=nothing)
    kT = PhysicsConstants.bk * T; sabflg = -225.0
    Ef, Epf = Float64(E), Float64(E_prime)
    (Ef <= 0 || Epf <= 0 || kT <= 0) && return 0.0
    beta_raw = (Epf - Ef) / kT
    alpha_raw = max((Ef + Epf - 2.0*Float64(mu)*sqrt(Ef*Epf))/(data.awr*kT), 1e-6)
    tevz = 0.0253
    bt = data.lat==1 ? abs(beta_raw)*kT/tevz : abs(beta_raw)
    at = data.lat==1 ? alpha_raw*kT/tevz : alpha_raw
    sb = something(sigma_b_override, data.sigma_b) * ((data.awr+1)/data.awr)^2
    pf = sqrt(Epf/Ef) / (2.0*kT) * sb
    ls = _interp_sab(at, bt, data)
    ls > sabflg && return pf * exp(ls - beta_raw/2)
    # SCT fallback
    Te = data.T_eff > 0 ? data.T_eff*PhysicsConstants.bk : kT
    bs, as = abs(beta_raw), alpha_raw
    arg = (as-bs)^2*kT/(4*as*Te) + (bs+beta_raw)/2
    arg < 225.0 ? pf*exp(-arg)/sqrt(4*PhysicsConstants.pi*as*Te/kT) : 0.0
end

"Integrated incoherent inelastic XS from S(a,b) [barn]."
function sab_xs(E::Real, data::SABData, T::Real; n_mu::Int=20, n_ep::Int=40)
    E <= 0 && return 0.0; kT = PhysicsConstants.bk*T
    Epm = max(E + data.beta[end]*(data.lat==1 ? 0.0253 : kT), 10*kT)
    dEp, dmu, total = Epm/n_ep, 2.0/n_mu, 0.0
    for iep in 1:n_ep; Ep=(iep-0.5)*dEp; Ep<=0 && continue
        ms = 0.0
        for imu in 0:n_mu; mu=-1.0+imu*dmu
            w = (imu==0||imu==n_mu) ? 0.5 : 1.0
            ms += w*sab_kernel(E, Ep, mu, data, T)
        end; total += ms*dmu*dEp
    end; total
end

# ---- Coherent elastic / Bragg edges (sigcoh(), thermr.f90:947-1207) ----

"Coherent elastic XS [barn]: sigma(E) = scon/E * sum_{edges below E} f_i."
function bragg_edges(E::Real, bragg::BraggData)
    E <= 0 && return 0.0; tf = 0.0
    for i in 1:bragg.n_edges
        bragg.tau_sq[i] >= E*bragg.econ && break; tf += bragg.form_factor[i]
    end; bragg.scon*tf/E
end

"Sorted Bragg edge energies [eV]."
bragg_edge_energies(b::BraggData) = [b.tau_sq[i]/b.econ for i in 1:b.n_edges]

"""
    build_bragg_data(; a, c, sigma_coh, A_mass, natom=1, debye_waller, emax=5.0)

Build Bragg data for hexagonal lattice. a,c in [cm], sigma_coh [barn].
"""
function build_bragg_data(; a::Real, c::Real, sigma_coh::Real,
                           A_mass::Real, natom::Int=1,
                           debye_waller::Real, emax::Real=5.0)
    pi_v = PhysicsConstants.pi
    amne = PhysicsConstants.amassn * PhysicsConstants.amu
    econ = PhysicsConstants.ev * 8 * amne / PhysicsConstants.hbar^2
    twopis = (2pi_v)^2
    scon = sigma_coh/natom*(4pi_v)^2 / (2*a^2*c*sqrt(3.0)*econ)
    t2 = PhysicsConstants.hbar / (2*PhysicsConstants.amu*A_mass)
    wint = 0.658173e-15 * A_mass * debye_waller
    c1, c2, ulim = 4.0/(3*a^2), 1.0/c^2, econ*emax
    tsl, ffl = Float64[], Float64[]
    phi = ulim/twopis; i1m = floor(Int, a*sqrt(phi))+1
    for i1 in 1:i1m; l1=i1-1
        disc = 3*(a^2*phi - l1^2); disc < 0 && continue
        i2m = floor(Int, 0.5*(l1+sqrt(disc)))+1
        for i2 in i1:i2m; l2=i2-1
            xv = phi - c1*(l1^2+l2^2-l1*l2)
            i3m = xv > 0 ? floor(Int, c*sqrt(xv))+1 : 1
            for i3 in 1:i3m; l3=i3-1
                w1 = l1==l2 ? 1.0 : 2.0
                w2 = (l1==0&&l2==0) ? 0.5 : (l1==0||l2==0) ? 1.0 : 2.0
                w3 = l3==0 ? 1.0 : 2.0
                for sgn in (1,-1); l2s=sgn*l2
                    tsq = (c1*(l1^2+l2s^2+l1*l2s)+l3^2*c2)*twopis
                    (tsq<=0||tsq>ulim) && continue
                    wt = exp(-tsq*t2*wint)*w1*w2*w3/sqrt(tsq)
                    found = false
                    for k in eachindex(tsl)
                        if abs(tsq-tsl[k]) < 0.05*tsl[k]
                            ffl[k]+=wt; found=true; break
                        end
                    end
                    found || (push!(tsl,tsq); push!(ffl,wt))
    end; end; end; end
    p = sortperm(tsl)
    BraggData(tsl[p], ffl[p], econ, scon, length(tsl))
end

# ---- Incoherent elastic (iel(), thermr.f90:1244-1425) ----

"Incoherent elastic XS: sigma(E) = sigma_b/2 * (1 - exp(-4*E*W'))."
function incoh_elastic_xs(E::Real, sigma_b::Real, dwp::Real)
    E <= 0 && return 0.0; sigma_b/2*(1.0-exp(-4.0*E*dwp))
end

# ---- Standard thermal energy grid (calcem, thermr.f90:1587-1607) ----

const THERMR_EGRID = Float64[
    1e-5,1.78e-5,2.5e-5,3.5e-5,5e-5,7e-5,1e-4,1.26e-4,1.6e-4,2e-4,
    2.53e-4,2.97e-4,3.5e-4,4.2e-4,5.06e-4,6.15e-4,7.5e-4,8.7e-4,
    1.012e-3,1.23e-3,1.5e-3,1.8e-3,2.03e-3,2.277e-3,2.6e-3,3e-3,
    3.5e-3,4.048e-3,4.5e-3,5e-3,5.6e-3,6.325e-3,7.2e-3,8.1e-3,
    9.108e-3,1e-2,1.063e-2,1.15e-2,1.2397e-2,1.33e-2,1.417e-2,1.5e-2,
    1.6192e-2,1.82e-2,1.99e-2,2.0493e-2,2.15e-2,2.28e-2,2.53e-2,
    2.8e-2,3.0613e-2,3.38e-2,3.65e-2,3.95e-2,4.2757e-2,4.65e-2,5e-2,
    5.6925e-2,6.25e-2,6.9e-2,7.5e-2,8.1972e-2,9e-2,9.6e-2,1.035e-1,
    1.11573e-1,1.2e-1,1.28e-1,1.355e-1,1.45728e-1,1.6e-1,1.72e-1,
    1.84437e-1,2e-1,2.277e-1,2.510392e-1,2.705304e-1,2.907501e-1,
    3.011332e-1,3.206421e-1,3.576813e-1,3.9e-1,4.170351e-1,4.5e-1,
    5.032575e-1,5.6e-1,6.25e-1,7e-1,7.8e-1,8.6e-1,9.5e-1,1.05,1.16,
    1.28,1.42,1.55,1.7,1.855,2.02,2.18,2.36,2.59,2.855,3.12,3.42,
    3.75,4.07,4.46,4.9,5.35,5.85,6.4,7.0,7.65,8.4,9.15,9.85,10.0]

# ---- Driver ----

"""
    compute_thermal(pendf, T, A; model=:free_gas, emax=10.0, sigma_b=nothing,
                    sab_data=nothing, bragg=nothing, debye_waller_prime=nothing)

Compute thermal XS and update PENDF material. Replaces elastic channel below emax.
"""
function compute_thermal(pendf::PointwiseMaterial, T::Real, A::Real;
                         model::Symbol=:free_gas, emax::Real=10.0,
                         sigma_b::Union{Real,Nothing}=nothing,
                         sab_data::Union{SABData,Nothing}=nothing,
                         bragg::Union{BraggData,Nothing}=nothing,
                         debye_waller_prime::Union{Real,Nothing}=nothing)
    el_col = findfirst(==(2), pendf.mt_list)
    el_col !== nothing || error("No elastic (MT2) column in PointwiseMaterial")
    tot_col = findfirst(==(1), pendf.mt_list)
    new_e, new_el = compute_thermal_xs(
        pendf.energies, pendf.cross_sections[:,el_col], A, T;
        sigma_b=sigma_b, emax=Float64(emax), model=model, sab_data=sab_data)
    n_new, n_reac = length(new_e), length(pendf.mt_list)
    new_xs = zeros(n_new, n_reac); new_xs[:,el_col] .= new_el
    for j in 1:n_reac; j==el_col && continue
        for i in 1:n_new
            new_xs[i,j] = _linterp(pendf.energies, pendf.cross_sections[:,j], new_e[i])
        end
    end
    if tot_col !== nothing
        for i in 1:n_new
            s = sum(new_xs[i,j] for j in 1:n_reac if j!=tot_col)
            new_xs[i,tot_col] = s
        end
    end
    PointwiseMaterial(pendf.mat, new_e, new_xs, copy(pendf.mt_list))
end

function _linterp(xs::AbstractVector, ys::AbstractVector, x::Float64)
    n = length(xs)
    x <= xs[1] && return ys[1]; x >= xs[n] && return ys[n]
    i = searchsortedfirst(xs, x); i == 1 && return ys[1]
    xs[i] == x && return ys[i]
    f = (x-xs[i-1])/(xs[i]-xs[i-1]); ys[i-1]+f*(ys[i]-ys[i-1])
end

"Compute thermal elastic XS on merged grid. Supports :free_gas and :sab."
function compute_thermal_xs(energies::AbstractVector{<:Real},
                            sigma_elastic::AbstractVector{<:Real},
                            A::Real, T::Real;
                            sigma_b::Union{Real,Nothing}=nothing, emax::Real=10.0,
                            model::Symbol=:free_gas,
                            sab_data::Union{SABData,Nothing}=nothing)
    @assert length(energies)==length(sigma_elastic) && issorted(energies)
    @assert T > 0 && A > 0
    if sigma_b === nothing
        idx = clamp(searchsortedfirst(energies, emax), 1, length(energies))
        sigma_b = A * sigma_elastic[idx]
    end
    em = T > 3000.0 ? emax*T/3000.0 : emax
    te = filter(e -> e <= em, THERMR_EGRID)
    txs = if model==:free_gas
        [free_gas_xs(e, A, T; sigma_b=Float64(sigma_b)) for e in te]
    elseif model==:sab
        sab_data===nothing && error("sab_data required for :sab")
        [sab_xs(e, sab_data, T) for e in te]
    else error("Unknown model: $model") end
    me = vcat(te, energies[energies .> em])
    mxs = vcat(txs, sigma_elastic[energies .> em])
    p = sortperm(me); (me[p], mxs[p])
end
