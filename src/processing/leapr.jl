# LEAPR -- S(alpha,beta) via phonon expansion (Proposal B)
# Numerically stable, AD-compatible, composable with THERMR SABData.
# Correspondence: _fsum->fsum, _start->start, _convol!->convol, _terpt->terpt,
#   generate_sab->contin+SCT, _bfact->bfact, _add_discrete_oscillators!->discre

using SpecialFunctions: besseli

"Phonon density of states on an equally spaced energy grid."
struct PhononDOS; energies::Vector{Float64}; density::Vector{Float64}; end

"Discrete oscillator: energy [eV] and fractional weight."
struct DiscreteOscillator; energy::Float64; weight::Float64; end

"S(alpha,beta) table. sab[i,j] = S_s(alpha_i, -beta_j) (asymmetric)."
struct SABTable
    alpha::Vector{Float64}; beta::Vector{Float64}; sab::Matrix{Float64}
    temperature::Float64; debye_waller::Float64
end

"Trapezoidal moment integral of phonon spectrum. Matches leapr.f90 fsum."
function _fsum(n::Int, p::AbstractVector{T}, np::Int, tau::T, delta::T) where T<:Real
    edsq = exp(delta * tau / 2); an = one(T) - 2 * mod(n, 2)
    be, fs, v = zero(T), zero(T), one(T)
    for ij in 1:np
        w = n > 0 ? be^n : one(T)
        ff = (p[ij] * v * v + p[ij] * an / (v * v)) * w
        (ij == 1 || ij == np) && (ff /= 2)
        fs += ff; be += delta; v *= edsq
    end
    fs * delta
end

"Linear interpolation in T_n(beta). Matches leapr.f90 terpt."
function _terpt(tn::AbstractVector{T}, ntn::Int, delta::T, be::T) where T<:Real
    be > ntn * delta && return zero(T)
    i = floor(Int, be / delta); i < ntn - 1 || return zero(T)
    idx = i + 1; bt = T(i) * delta
    tn[idx] + (be - bt) * (tn[idx+1] - tn[idx]) / delta
end

"Convolve t1 with tlast -> tnext. Returns norm check. Matches leapr.f90 convol."
function _convol!(t1::Vector{Float64}, tlast::Vector{Float64},
                  tnext::Vector{Float64}, n1::Int, nl::Int, nn::Int, delta::Float64)
    ckk = 0.0
    for k in 1:nn
        tnext[k] = 0.0
        for j in 1:n1
            t1[j] > 0.0 || continue
            i1 = k+j-2; i2 = k-j; be = (j-1)*delta
            f1 = i1+1 <= nl ? tlast[i1+1]*exp(-be) : 0.0
            f2 = if i2 >= 0 && i2+1 <= nl; tlast[i2+1]
                 elseif i2 < 0 && 1-i2 <= nl; tlast[1-i2]*exp(i2*delta)
                 else 0.0 end
            cc = t1[j]*(f1+f2); (j == 1 || j == n1) && (cc /= 2.0)
            tnext[k] += cc
        end
        tnext[k] *= delta; tnext[k] < 1e-30 && (tnext[k] = 0.0)
        be = (k-1)*delta; cc = tnext[k] + tnext[k]*exp(-be)
        (k == 1 || k == nn) && (cc /= 2.0); ckk += cc
    end
    ckk * delta
end

"Normalise rho -> P(beta), compute DW lambda and T_eff ratio, build T1."
function _start(rho_raw::Vector{Float64}, np::Int, delta_e::Float64,
                T::Float64, tbeta::Float64)
    kT = PhysicsConstants.bk * T; deltab = delta_e / kT
    p = copy(rho_raw)
    p[1] = np > 1 ? p[2] / deltab^2 : 0.0
    v = exp(deltab/2.0); vv = v; u = deltab
    for j in 2:np
        denom = u*(vv - 1.0/vv); p[j] = denom > 0 ? p[j]/denom : 0.0
        vv *= v; u += deltab
    end
    tau = 0.5; an = _fsum(1, p, np, tau, deltab)/tbeta
    an > 0 || error("LEAPR: DOS normalisation failed"); p ./= an
    f0 = _fsum(0, p, np, tau, deltab)
    tbar = _fsum(2, p, np, tau, deltab)/(2.0*tbeta)
    for i in 1:np; be = deltab*(i-1); p[i] = p[i]*exp(be/2.0)/f0; end
    (p, deltab, f0, tbar)
end

"Bessel function weights for one discrete oscillator. Matches leapr.f90 bfact."
function _bfact(x::Float64, dwc::Float64, betai::Float64)
    bplus = zeros(50); bminus = zeros(50)
    bzero = besseli(0, x) * exp(-dwc)
    for i in 1:50
        bn = besseli(i, x); bn > 1e-30 || continue
        bplus[i] = exp(-dwc - i*betai/2)*bn; bplus[i] < 1e-30 && (bplus[i] = 0.0)
        bminus[i] = exp(-dwc + i*betai/2)*bn; bminus[i] < 1e-30 && (bminus[i] = 0.0)
    end
    (bzero, bplus, bminus)
end

"Log-linear interpolation in ssm[nal,:]. Handles detailed balance for negative beta."
function _interp_sab_log(ssm::Matrix{Float64}, nal::Int,
                         bg::Vector{Float64}, be::Float64)
    slim = -225.0; nb = length(bg)
    if be < 0.0
        s = _interp_sab_log(ssm, nal, bg, -be)
        return s > 0 ? s*exp(be) : 0.0
    end
    be > bg[end] && return 0.0
    idx = clamp(searchsortedlast(bg, be), 1, nb-1)
    s1, s2 = ssm[nal, idx], ssm[nal, idx+1]
    db = bg[idx+1] - bg[idx]; db <= 0 && return max(s1, 0.0)
    ls1 = s1 > 0 ? log(s1) : slim; ls2 = s2 > 0 ? log(s2) : slim
    frac = clamp((be - bg[idx])/db, 0.0, 1.0)
    lval = ls1 + frac*(ls2 - ls1); lval > slim ? exp(lval) : 0.0
end

"Convolve discrete oscillator delta functions with continuous S(a,b)."
function _add_discrete_oscillators!(ssm::Matrix{Float64}, ag::Vector{Float64},
                                    bg::Vector{Float64},
                                    oscs::Vector{DiscreteOscillator},
                                    dw0::Float64, kT::Float64)
    na, nb = size(ssm); nd = length(oscs); maxdd = 500
    bdeln = [o.energy/kT for o in oscs]; adel = [o.weight for o in oscs]
    eb = exp.(bdeln./2); sn = @. (eb-1/eb)/2; cn = @. (eb+1/eb)/2
    ar = @. adel/(sn*bdeln); dbw = @. ar*cn; dw_total = dw0 + sum(dbw)
    sexpb = zeros(nb); bes = zeros(maxdd); wts = zeros(maxdd)
    ben = zeros(maxdd); wtn = zeros(maxdd)
    for nal in 1:na
        al = ag[nal]; fill!(sexpb, 0.0)
        ben[1] = 0.0; wtn[1] = 1.0; nn = 1; n = 0
        for i in 1:nd
            bz, bp, bm = _bfact(al*ar[i], al*dbw[i], bdeln[i])
            for m in 1:nn
                (ben[m] <= 0 || wtn[m]*bz >= 1e-8) && n < maxdd &&
                    (n += 1; bes[n] = ben[m]; wts[n] = wtn[m]*bz)
            end
            for k in 1:50
                bm[k] <= 0 && break
                for m in 1:nn
                    w = wtn[m]*bm[k]; w >= 1e-8 && n < maxdd &&
                        (n += 1; bes[n] = ben[m]-k*bdeln[i]; wts[n] = w)
                end
            end
            for k in 1:50
                bp[k] <= 0 && break
                for m in 1:nn
                    w = wtn[m]*bp[k]; w >= 1e-8 && n < maxdd &&
                        (n += 1; bes[n] = ben[m]+k*bdeln[i]; wts[n] = w)
                end
            end
            nn = n; ben[1:nn] .= bes[1:nn]; wtn[1:nn] .= wts[1:nn]; n = 0
        end
        for m in 1:nn, j in 1:nb
            st = _interp_sab_log(ssm, nal, bg, bg[j]+bes[m])
            add = wts[m]*st; add >= 1e-20 && (sexpb[j] += add)
        end
        ssm[nal, :] .= sexpb
    end
    dw_total
end

# ---- Public API ----

"Compute Debye-Waller factor lambda_s from phonon DOS at temperature T [K]."
function debye_waller_factor(dos::PhononDOS, T::Float64; tbeta::Float64=1.0)
    np = length(dos.density)
    de = np > 1 ? dos.energies[2]-dos.energies[1] : error("Need >= 2 points")
    _, _, f0, _ = _start(dos.density, np, de, T, tbeta); f0
end

"Compute first n_terms phonon expansion T_l(beta). Returns (terms, deltab, f0)."
function phonon_expansion(dos::PhononDOS, T::Float64, n_terms::Int; tbeta::Float64=1.0)
    np = length(dos.density); de = dos.energies[2]-dos.energies[1]
    p, deltab, f0, _ = _start(dos.density, np, de, T, tbeta)
    maxl = n_terms*np; tlast = zeros(maxl); tnow = zeros(maxl)
    tlast[1:np] .= p; terms = Vector{Vector{Float64}}(undef, n_terms); terms[1] = copy(p)
    npl = np
    for n in 2:n_terms
        npn = np+npl-1; _convol!(p, tlast, tnow, np, npl, npn, deltab)
        terms[n] = tnow[1:npn]; tlast[1:npn] .= tnow[1:npn]; npl = npn
    end
    (terms, deltab, f0)
end

default_alpha_grid(; n::Int=20) = 10.0 .^ range(log10(0.01), log10(10.0), length=n)
default_beta_grid(; n::Int=40) = collect(range(0.0, 20.0, length=n))

"""
    generate_sab(dos, T; ...) -> SABTable

Phonon expansion S(alpha,beta) with SCT fallback. Matches NJOY contin()+start().
"""
function generate_sab(dos::PhononDOS, T::Float64;
                      alpha_grid::Vector{Float64}=default_alpha_grid(),
                      beta_grid::Vector{Float64}=default_beta_grid(),
                      n_phonon_terms::Int=100, tbeta::Float64=1.0,
                      oscillators::Vector{DiscreteOscillator}=DiscreteOscillator[])
    np = length(dos.density); de = dos.energies[2]-dos.energies[1]
    p, deltab, f0, tbar = _start(dos.density, np, de, T, tbeta)
    kT = PhysicsConstants.bk * T; na, nb = length(alpha_grid), length(beta_grid)
    ssm = zeros(na, nb); explim = -250.0; tiny = 1.0e-30
    maxl = n_phonon_terms*np+np; tlast = zeros(maxl); tnow = zeros(maxl); xa = zeros(na)
    # l=1 term
    tlast[1:np] .= p
    for j in 1:na
        al = alpha_grid[j]; xa[j] = log(al*f0)
        ex = -f0*al+xa[j]; exx = ex > explim ? exp(ex) : 0.0
        for k in 1:nb
            st = _terpt(p, np, deltab, beta_grid[k])
            ssm[j,k] = st*exx < tiny ? 0.0 : st*exx
        end
    end
    npl = np; maxt = fill(na+1, nb)
    # Phonon expansion l=2..n
    for n in 2:n_phonon_terms
        npn = np+npl-1; _convol!(p, tlast, tnow, np, npl, npn, deltab)
        for j in 1:na
            al = alpha_grid[j]; xa[j] += log(al*f0/n)
            ex = -f0*al+xa[j]; exx = ex > explim ? exp(ex) : 0.0
            for k in 1:nb
                st = _terpt(tnow, npn, deltab, beta_grid[k]); add = st*exx
                add >= tiny && (ssm[j,k] += add)
                n >= n_phonon_terms && ssm[j,k] > 0 && add > ssm[j,k]/1000 &&
                    j < maxt[k] && (maxt[k] = j)
            end
        end
        tlast[1:npn] .= tnow[1:npn]; npl = npn
    end
    # SCT fallback for unconverged alpha
    for i in 2:nb; maxt[i] > maxt[i-1] && (maxt[i] = maxt[i-1]); end
    pi_v = PhysicsConstants.pi
    for k in 1:nb, j in 1:na
        if j >= maxt[k]
            al = alpha_grid[j]; be = beta_grid[k]
            alw = al*tbeta; alp = alw*tbar; ex = -(alw-be)^2/(4*alp)
            ssm[j,k] = ex > explim ? exp(ex)/sqrt(4*pi_v*alp) : 0.0
        end
    end
    dw = f0
    !isempty(oscillators) && (dw = _add_discrete_oscillators!(ssm, alpha_grid,
                                    beta_grid, oscillators, f0, kT))
    SABTable(copy(alpha_grid), copy(beta_grid), ssm, T, dw)
end

"Convolve oscillators into existing SABTable."
function add_discrete_oscillators!(sab::SABTable, oscs::Vector{DiscreteOscillator}, T::Float64)
    ssm = copy(sab.sab); kT = PhysicsConstants.bk*T
    dw = _add_discrete_oscillators!(ssm, sab.alpha, sab.beta, oscs, sab.debye_waller, kT)
    SABTable(sab.alpha, sab.beta, ssm, T, dw)
end

"Convert SABTable to SABData for THERMR. Stores log(S_sym) where S_sym = S_s*exp(-b/2)."
function sab_table_to_thermr(sab::SABTable; sigma_b::Float64=1.0, awr::Float64=1.0, lat::Int=0)
    na, nb = size(sab.sab); sabflg = -225.0
    logsab = Matrix{Float64}(undef, na, nb)
    for j in 1:na, k in 1:nb
        s_sym = sab.sab[j,k]*exp(-sab.beta[k]/2.0)
        logsab[j,k] = s_sym > exp(sabflg) ? log(s_sym) : sabflg
    end
    SABData(sab.alpha, sab.beta, logsab, sigma_b, awr, sab.temperature, 0, lat)
end

"Debye DOS: rho(omega) = 3*omega^2/omega_D^3 for omega <= omega_D."
function debye_dos(T_debye::Float64, n_points::Int; emax_factor::Float64=1.2)
    omega_D = PhysicsConstants.bk * T_debye; emax = omega_D*emax_factor
    energies = collect(range(0.0, emax, length=n_points))
    density = [e <= omega_D ? 3.0*e^2/omega_D^3 : 0.0 for e in energies]
    PhononDOS(energies, density)
end
