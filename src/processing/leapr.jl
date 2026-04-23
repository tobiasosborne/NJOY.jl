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

# =================================================================
# Fortran-faithful transport corrections for leapr (trans/coldh/skold/coher).
#
# These operate on a β-fastest `ssm[β, α, T]` array to match Fortran's
# column-major `ssm(nbeta,nalpha,ntempr)` layout and IEEE-754 accumulation
# order exactly. Ref: njoy-reference/src/leapr.f90
# =================================================================

"Modified Bessel K_1 via NJOY leapr.f90:1253 rational approximation.
For x>1 returns e^x * K_1(x) (scaled); stable() compensates."
function _besk1(x::Float64)
    if x <= 1.0
        c0=0.125
        c1=0.442850424; c2=0.584115288;  c3=6.070134559
        c4=17.864913364; c5=48.858995315; c6=90.924600045
        c7=113.795967431; c8=85.331474517; c9=32.00008698; c10=3.999998802
        c11=1.304923514; c12=1.47785657; c13=16.402802501
        c14=44.732901977; c15=115.837493464; c16=198.437197312
        c17=222.869709703; c18=142.216613971; c19=40.000262262; c20=1.999996391
        c21=1.0; c22=0.5; c23=0.5772156649
        v = c0 * x
        u = v * v
        bi1 = (((((((((c1*u+c2)*u+c3)*u+c4)*u+c5)*u+c6)*u+c7)*u+c8)*u+c9)*u+c10)*v
        bi3 = (((((((((c11*u+c12)*u+c13)*u+c14)*u+c15)*u+c16)*u+c17)*u+c18)*u+c19)*u+c20)
        return c21/x + bi1*(log(c22*x) + c23) - v*bi3
    else
        c25=0.0108241775; c26=0.0788000118; c27=0.2581303765
        c28=0.5050238576; c29=0.663229543; c30=0.6283380681
        c31=0.4594342117; c32=0.2847618149; c33=0.1736431637
        c34=0.1280426636; c35=0.1468582957; c36=0.4699927013
        c37=1.2533141373
        u = 1.0 / x
        bi3 = ((((((((((((-c25*u+c26)*u-c27)*u+c28)*u-c29)*u+c30)*u-c31)*u+c32)*u-c33)*u+c34)*u-c35)*u+c36)*u+c37)
        return sqrt(u) * bi3
    end
end

"""
    _stable!(ap, sd, al, delta, twt, c_diff, ndmax) -> nsd

Fill `sd[1:nsd]` with the translational transport kernel at spacing `delta`:
diffusion kernel (`c_diff>0`, via `besk1`) or free-gas Gaussian (`c_diff==0`).
`ap[1:nsd]` receives the matching β values. Growth stops when
`sd[j] < 1e-7*sd[1]` or `j>=ndmax`; nsd is forced odd (Simpson).

Ref: leapr.f90:1009-1122 (stable).
"""
function _stable!(ap::Vector{Float64}, sd::Vector{Float64},
                  al::Float64, delta::Float64,
                  twt::Float64, c_diff::Float64, ndmax::Int)
    eps = 1e-7
    PI = Float64(pi)
    if c_diff != 0.0
        # Diffusion branch
        d = twt * c_diff
        c2s = sqrt(c_diff*c_diff + 0.25)
        c3 = 2 * d * al
        c4 = c3 * c3
        c8 = c2s * c3 / PI
        c3 = 2 * d * c_diff * al
        be = 0.0
        j = 1
        while true
            c6 = sqrt(be*be + c4)
            c7 = c6 * c2s
            c5 = c7 <= 1.0 ? c8*exp(c3 + be/2) : c8*exp(c3 - c7 + be/2)
            sd[j] = c5 * _besk1(c7) / c6
            ap[j] = be
            be += delta
            j += 1
            if iseven(j)
                (j >= ndmax) && break
                (eps * sd[1] >= sd[j-1]) && break
            end
        end
        return j - 1
    else
        # Free-gas branch
        be = 0.0; j = 1
        wal = twt * al
        while true
            ex = -(wal - be)^2 / (4*wal)
            sd[j] = exp(ex) / sqrt(4*PI*wal)
            ap[j] = be
            be += delta
            j += 1
            if iseven(j)
                (j >= ndmax) && break
                (eps * sd[1] >= sd[j-1]) && break
            end
        end
        return j - 1
    end
end

"Log-linear interpolation in `sd` (spacing delta) at β = be. Returns 0 past table.
Ref: leapr.f90:1124-1162 (terps)."
function _terps(sd::AbstractVector{Float64}, nsd::Int, delta::Float64, be::Float64)
    slim = -225.0
    be > delta * nsd && return 0.0
    i0 = floor(Int, be / delta)
    if i0 < nsd - 1
        bt = i0 * delta
        btp = bt + delta
        i = i0 + 1
        st  = sd[i]   <= 0 ? slim : log(sd[i])
        stp = sd[i+1] <= 0 ? slim : log(sd[i+1])
        stt = st + (be - bt)*(stp - st)/(btp - bt)
        return stt > slim ? exp(stt) : 0.0
    end
    return 0.0
end

"""
    _sbfill!(sb, nbt, delta, be, s, betan, nbeta, ndmax; delta_ref)

Build S(β) on a symmetric grid of 2·nbt-1 points centered on −be, with
spacing `delta`, by log-linear interpolation in the tabulated `s[1:nbeta]`
over `betan[1:nbeta]`. Applies `sb *= exp(-bet)` for `bet>0` to map
S(α,-|β|) onto the positive-β branch via detailed balance.

Ref: leapr.f90:1164-1251 (sbfill).
"""
function _sbfill!(sb::Vector{Float64}, nbt::Int, delta::Float64, be::Float64,
                  s::Vector{Float64}, betan::Vector{Float64},
                  nbeta::Int, ndmax::Int)
    shade = 1.00001; slim = -225.0
    bmin = -be - (nbt-1)*delta
    bmax = -be + (nbt-1)*delta + delta/100
    need = 1 + floor(Int, (bmax - bmin)/delta)
    need > ndmax && error("sbfill: need $need scratch slots, have $ndmax")
    j = nbeta; i = 0
    bet = bmin
    while bet <= bmax
        i += 1
        b = abs(bet)
        idone = 0
        # Bracket search: walk j so that betan[j-1] < b <= betan[j]
        while idone == 0
            if b > betan[j]
                if j == nbeta && b < shade * betan[j]
                    idone = 1
                elseif j == nbeta
                    idone = 2              # out-of-range
                else
                    j += 1
                end
            else
                if b > betan[j-1]
                    idone = 1
                elseif j == 2
                    idone = 1
                else
                    j -= 1
                end
            end
        end
        if idone == 1
            st  = s[j]   <= 0 ? slim : log(s[j])
            stm = s[j-1] <= 0 ? slim : log(s[j-1])
            sb[i] = st + (b - betan[j])*(stm - st)/(betan[j-1] - betan[j])
            bet > 0 && (sb[i] -= bet)    # detailed balance on +β branch
            arg = sb[i]
            sb[i] = arg > slim ? exp(arg) : 0.0
        else
            sb[i] = 0.0
        end
        # Delta safety (FP-degenerate step)
        while bet == bet + delta
            delta *= 10
        end
        bet += delta
    end
    return nothing
end

"""
    trans!(ssm, itemp, alpha_grid, beta_grid, twt, c_diff, tbeta,
           tev, deltab, f0, lat, arat, tempr, tempf)

Add translational (diffusion or free-gas) contribution to the asymmetric
S(α,-β) stored in `ssm[β, α, T]`. Mutates `ssm[:,:,itemp]` in place and
updates `tempf[itemp]`.

Array layout: `ssm[ibeta, ialpha, itemp]` — β is the fastest axis, matching
Fortran's `ssm(nbeta,nalpha,ntempr)` so Simpson accumulation order is
preserved exactly.

Ref: leapr.f90:844-1007 (trans).
"""
function trans!(ssm::AbstractArray{Float64,3}, itemp::Int,
                alpha_grid::Vector{Float64}, beta_grid::Vector{Float64},
                twt::Float64, c_diff::Float64, tbeta::Float64,
                tev::Float64, deltab::Float64, f0::Float64,
                lat::Int, arat::Float64,
                tempr::Vector{Float64}, tempf::Vector{Float64})
    therm = 0.0253
    c0 = 0.4; c1 = 1.0; c2 = 1.42; c3c = 0.2; c4c = 10.0
    TINY = 1e-30

    nbeta  = length(beta_grid)
    nalpha = length(alpha_grid)
    sc     = lat == 1 ? therm/tev : 1.0

    ndmax = max(nbeta, 1_000_000)
    ap    = Vector{Float64}(undef, ndmax)
    sd    = Vector{Float64}(undef, ndmax)
    sb    = Vector{Float64}(undef, ndmax)
    betan = Vector{Float64}(undef, nbeta)

    for ialpha in 1:nalpha
        al = alpha_grid[ialpha] * sc / arat

        # Step δ — two criteria, take the smaller
        w   = twt * c_diff * al
        ded = c_diff != 0.0 ? c0 * w / sqrt(c1 + c2*w*c_diff) : c3c*sqrt(twt*al)
        ded == 0.0 && (ded = c3c*sqrt(twt*al))
        deb = c4c * al * deltab
        delta = min(ded, deb)

        nsd = _stable!(ap, sd, al, delta, twt, c_diff, ndmax)
        nsd <= 1 && continue

        for i in 1:nbeta
            betan[i] = beta_grid[i] * sc
            ap[i]    = ssm[i, ialpha, itemp]
        end

        nbt = nsd
        for ibeta in 1:nbeta
            be = betan[ibeta]
            _sbfill!(sb, nbt, delta, be, ap, betan, nbeta, ndmax)

            # Simpson convolution (left-to-right order preserved for FP)
            s = 0.0
            for i in 1:nbt
                f = (i == 1 || i == nbt) ? 1.0 : 2.0 * (mod(i-1, 2) + 1)
                s += f * sd[i] * sb[nbt + i - 1]
                bb = (i-1) * delta
                s += f * sd[i] * sb[nbt - i + 1] * exp(-bb)
            end
            s *= delta / 3.0
            s < TINY && (s = 0.0)

            st = _terps(sd, nbt, delta, be)
            st > 0.0 && (s += exp(-al * f0) * st)

            ssm[ibeta, ialpha, itemp] = s
        end
    end

    # Effective temperature update (leapr.f90:998)
    tempf[itemp] = (tbeta * tempf[itemp] + twt * tempr[itemp]) / (tbeta + twt)
    return nothing
end

# =================================================================
# coldh — ortho/para H₂ and D₂ discrete-rotational convolutions.
# Ref: leapr.f90:1936-2183 (coldh) and helpers (bt/sumh/cn/sjbes/
# terpk/bfill/exts/sint).
# =================================================================

"Statistical weight p_j for rotational state j at x = δE/kT.
Ref: leapr.f90:2185-2209 (bt)."
function _bt(j::Int, x::Float64)
    a = (2j+1) * exp(-0.5 * j * (j+1) * x)
    b = 0.0
    for i in 1:10
        k = 2i - 2
        isodd(j) && (k += 1)
        b += (2k+1) * exp(-0.5 * k * (k+1) * x)
    end
    return a / (2b)
end

"Spherical Bessel function j_n(x) via Miller's downward recurrence.
Ref: leapr.f90:2342-2442 (sjbes)."
function _sjbes(n::Int, x::Float64)
    HUGE = 1e25; SMALL = 2e-38
    (n >= 30000 || x > 3e4) && return 0.0
    (x < 0.0 || n < 0) && error("sjbes: invalid args n=$n x=$x")

    if x <= 7e-4
        n == 0 && return 1.0
        n > 10 && return 0.0
        t1 = 3.0; t2 = 1.0; t3 = 1.0
        for _ in 1:n
            t3 = t2 * x / t1
            t1 += 2.0
            t2  = t3
        end
        return t3
    else
        w = x < 0.2 ? 1 - x^2*(1 - x^2/20)/6 : sin(x)/x
        n == 0 && return w
        l = x >= 100.0 ? floor(Int, x/50 + 18) :
            x >= 10.0  ? floor(Int, x/10 + 10) :
            x >  1.0   ? floor(Int, x/2  + 5)  : 5
        iii  = floor(Int, x)
        kmax = max(n, iii)
        nm   = kmax + l
        z    = 1.0 / x
        t3   = 0.0; t2 = SMALL
        sj   = 0.0; t1 = 0.0
        for i in 1:nm
            k  = nm - i
            t1 = (2k + 3) * z * t2 - t3
            k == n && (sj = t1)
            if abs(t1) >= HUGE
                t1 /= HUGE; t2 /= HUGE; sj /= HUGE
            end
            t3 = t2; t2 = t1
        end
        return w * sj / t1
    end
end

"Clebsch-Gordon coefficient via log-factorials.
Ref: leapr.f90:2247-2340 (cn). Returns 0 when jj+ll+nn is odd."
function _cn(jj::Int, ll::Int, nn::Int)
    tot = jj + ll + nn
    isodd(tot) && return 0.0
    ka1 = tot
    ka2 = jj + ll - nn
    ka3 = jj - ll + nn
    ka4 = ll - jj + nn
    kb1 = ka1 ÷ 2; kb2 = ka2 ÷ 2; kb3 = ka3 ÷ 2; kb4 = ka4 ÷ 2

    _logfact(k::Int) = k <= 0 ? 0.0 : sum(log, 1:k)
    a1 = sqrt(exp(_logfact(ka1))); a2 = sqrt(exp(_logfact(ka2)))
    a3 = sqrt(exp(_logfact(ka3))); a4 = sqrt(exp(_logfact(ka4)))
    b1 = exp(_logfact(kb1)); b2 = exp(_logfact(kb2))
    b3 = exp(_logfact(kb3)); b4 = exp(_logfact(kb4))

    rat = (2nn + 1) / (tot + 1)
    iwign = (jj + ll - nn) ÷ 2
    sgn = iseven(iwign) ? 1.0 : -1.0
    return sgn * sqrt(rat) * b1 / a1 * a2 / b2 * a3 / b3 * a4 / b4
end

"Sum over Bessel functions × Clebsch-Gordon² for j→jp rotational transition.
Ref: leapr.f90:2211-2245 (sumh)."
function _sumh(j::Int, jp::Int, y::Float64)
    if j == 0
        return (_sjbes(jp, y) * _cn(j, jp, jp))^2
    elseif jp == 0
        return (_sjbes(j, y) * _cn(j, 0, j))^2
    else
        imk = abs(j - jp) + 1
        ipk1 = j + jp + 1
        mpk = ipk1 - imk
        ipk = mpk <= 9 ? ipk1 : imk + 9
        s = 0.0
        for n in imk:ipk
            n1 = n - 1
            s += (_sjbes(n1, y) * _cn(j, jp, n1))^2
        end
        return s
    end
end

"Linear interpolation in `ska(κ=δ·i)` at wave number `be`. Returns 1 past table.
Ref: leapr.f90:2444-2466 (terpk)."
function _terpk(ska::Vector{Float64}, nka::Int, delta::Float64, be::Float64)
    be > nka * delta && return 1.0
    i0 = floor(Int, be / delta)
    i0 < nka - 1 || return 1.0
    bt = i0 * delta
    i  = i0 + 1
    return ska[i] + (be - bt) * (ska[i+1] - ska[i]) / delta
end

"""
Build extended β grid spanning ±β from the input negative-β grid `betan[1..nbeta]`.
`bex[1..nbx]` receives monotone [-βmax ... 0 ... +βmax]; `rdbex[i] = 1/(bex[i+1]-bex[i])`.
Returns `nbx`. If `betan[1] ≈ 0`, the zero point is not duplicated.

Ref: leapr.f90:1798-1832 (bfill).
"""
function _bfill!(bex::Vector{Float64}, rdbex::Vector{Float64},
                 betan::Vector{Float64}, nbeta::Int, maxbb::Int)
    SMALL = 1e-9
    k = nbeta
    for i in 1:nbeta
        bex[i] = -betan[k]
        k -= 1
    end
    if betan[1] <= SMALL
        bex[nbeta] = 0.0
        k = nbeta + 1
    else
        k = nbeta + 2
        bex[nbeta+1] = betan[1]
    end
    for i in 2:nbeta
        bex[k] = betan[i]
        k += 1
    end
    nbxm = k - 2
    nbx  = k - 1
    for i in 1:nbxm
        rdbex[i] = 1.0 / (bex[i+1] - bex[i])
    end
    return nbx
end

"Extend S(α,-β) slice `sexpb` to full ±β array `sex` via detailed balance
(`sex[k] = sexpb[i]·exp(-β_i)` on the positive branch).
Ref: leapr.f90:1834-1865 (exts)."
function _exts!(sexpb::AbstractVector{Float64}, sex::Vector{Float64},
                exb::Vector{Float64}, betan::Vector{Float64},
                nbeta::Int, maxbb::Int)
    SMALL = 1e-9
    k = nbeta
    for i in 1:nbeta
        sex[i] = sexpb[k]
        k -= 1
    end
    if betan[1] <= SMALL
        sex[nbeta] = sexpb[1]
        k = nbeta + 1
    else
        k = nbeta + 2
        sex[nbeta+1] = sexpb[1]
    end
    for i in 2:nbeta
        sex[k] = sexpb[i] * exb[i] * exb[i]
        k += 1
    end
    return nothing
end

"""
Interpolate S(α,β) at β=x in the extended table, or use SCT approximation
when |x| > βmax. Log-linear in β.

Ref: leapr.f90:1867-1934 (sint).
"""
function _sint(x::Float64, bex::Vector{Float64}, rdbex::Vector{Float64},
               sex::Vector{Float64}, nbx::Int, alph::Float64, wt::Float64,
               tbart::Float64, betan::Vector{Float64}, nbeta::Int, maxbb::Int)
    SLIM = -225.0
    PI = Float64(pi)
    # SCT approximation beyond the β-table
    if abs(x) > betan[nbeta]
        alph <= 0 && return 0.0
        ex = -(wt*alph - abs(x))^2 / (4*wt*alph*tbart)
        x > 0 && (ex -= x)
        return exp(ex) / (4*PI*wt*alph*tbart)
    end
    # Bisection for x in bex
    k1 = 1; k2 = nbeta; k3 = nbx
    while true
        if x == bex[k2]
            return sex[k2]
        elseif x > bex[k2]
            k1 = k2
            k2 = (k3 - k2) ÷ 2 + k2
            k3 - k1 <= 1 && break
        else
            k3 = k2
            k2 = (k2 - k1) ÷ 2 + k1
            k3 - k1 <= 1 && break
        end
    end
    ss1 = sex[k1] <= 0 ? SLIM : log(sex[k1])
    ss3 = sex[k3] <= 0 ? SLIM : log(sex[k3])
    ex = ((bex[k3] - x)*ss1 + (x - bex[k1])*ss3) * rdbex[k1]
    return ex > SLIM ? exp(ex) : 0.0
end

"""
    coldh!(ssm, ssp, itemp, alpha_grid, beta_grid,
           ska, nka, dka, ncold, nsk,
           temp, tev, twt, tbeta, lat, arat, tempr, tempf; ifree=0)

Convolve the current S(α,-β) in `ssm[β,α,T]` with discrete rotational modes
for ortho/para H₂ (`ncold=1,2`) or ortho/para D₂ (`ncold=3,4`). Populates
both `ssm` (for negative β) and `ssp` (for positive β) → final S is
**asymmetric**.

`ssp` must be pre-allocated with the same shape as `ssm`. Caller passes in
the static structure factor table `ska[1..nka]` sampled at `κ = dka·i`
(from Card 17-18 of the leapr deck).

Ref: leapr.f90:1936-2183 (coldh).
"""
function coldh!(ssm::AbstractArray{Float64,3}, ssp::AbstractArray{Float64,3},
                itemp::Int,
                alpha_grid::Vector{Float64}, beta_grid::Vector{Float64},
                ska::Vector{Float64}, nka::Int, dka::Float64,
                ncold::Int, nsk::Int,
                temp::Float64, tev::Float64,
                twt::Float64, tbeta::Float64,
                lat::Int, arat::Float64,
                tempr::Vector{Float64}, tempf::Vector{Float64};
                ifree::Int=0, nokap::Int=0, jterm::Int=3)

    # Hardcoded Fortran constants (leapr.f90:1962-1976)
    pmass  = 1.6726231e-24
    dmass  = 3.343586e-24
    deh    = 0.0147
    ded    = 0.0074
    sampch = 0.356
    sampcd = 0.668
    sampih = 2.526
    sampid = 0.403
    therm  = 0.0253
    angst  = 1e-8
    SMALL  = 1e-6
    PI     = Float64(pi)

    nbeta  = length(beta_grid)
    nalpha = length(alpha_grid)
    maxbb  = 2 * nbeta + 1
    betan  = Vector{Float64}(undef, nbeta)
    exb    = Vector{Float64}(undef, nbeta)
    bex    = Vector{Float64}(undef, maxbb)
    rdbex  = Vector{Float64}(undef, maxbb)
    sex    = Vector{Float64}(undef, maxbb)

    sc  = lat == 1 ? therm/tev : 1.0
    law = ncold + 1
    de  = law > 3 ? ded : deh
    x_rot = de / tev

    if law > 3
        amassm = 6.69e-24
        sampc  = sampcd
        bp     = PhysicsConstants.hbar/2 * sqrt(2/ded/PhysicsConstants.ev/dmass) / angst
        sampi  = sampid
    else
        amassm = 3.3464e-24
        sampc  = sampch
        bp     = PhysicsConstants.hbar/2 * sqrt(2/deh/PhysicsConstants.ev/pmass) / angst
        sampi  = sampih
    end
    wt    = twt + tbeta
    tbart = tempf[itemp] / tempr[itemp]
    snorm = sampi^2 + sampc^2

    nbx = 0
    for nal in 1:nalpha
        al    = alpha_grid[nal] * sc / arat
        alp   = wt * al
        waven = angst * sqrt(amassm * tev * PhysicsConstants.ev * al) / PhysicsConstants.hbar
        y     = bp * waven
        sk    = nokap == 1 ? 1.0 : _terpk(ska, nka, dka, waven)

        # Spin-correlation weights (Fortran :2033-:2043)
        swe, swo = if law == 2         # ortho-H₂
            (sampi^2/3, sk*sampc^2 + 2*sampi^2/3)
        elseif law == 3                # para-H₂
            (sk*sampc^2, sampi^2)
        elseif law == 4                # ortho-D₂
            (sk*sampc^2 + 5*sampi^2/8, 3*sampi^2/8)
        elseif law == 5                # para-D₂
            (3*sampi^2/4, sk*sampc^2 + sampi^2/4)
        else
            error("coldh: unsupported ncold=$ncold (expected 1..4)")
        end
        swe /= snorm
        swo /= snorm

        # One-time β-grid extension at first α
        if nal == 1
            for i in 1:nbeta
                be = beta_grid[i]
                lat == 1 && (be *= therm/tev)
                exb[i]   = exp(-be/2)
                betan[i] = be
            end
            nbx = _bfill!(bex, rdbex, betan, nbeta, maxbb)
        end
        _exts!(view(ssm, 1:nbeta, nal, itemp), sex, exb, betan, nbeta, maxbb)

        # Extended-β loop (jj = 1..2·nbeta-1)
        jjmax = 2*nbeta - 1
        for jj in 1:jjmax
            k = jj < nbeta ? (nbeta - jj + 1) : (jj - nbeta + 1)
            be = betan[k]
            jj < nbeta && (be = -be)
            sn = 0.0

            # Rotational state j loop (ipo=1 for even-j start, ipo=2 for odd)
            ipo = (law == 2 || law == 5) ? 2 : 1
            jt1 = 2*jterm
            ipo == 2 && (jt1 += 1)

            l = ipo
            while l <= jt1
                j = l - 1
                pj = _bt(j, x_rot)

                # Even j' transitions (swe weight)
                snlg = 0.0
                for lp in 1:2:10
                    jp = lp - 1
                    betap = (-j*(j+1) + jp*(jp+1)) * x_rot / 2
                    tmp = (2*jp + 1) * pj * swe * 4 * _sumh(j, jp, y)
                    bn = be + betap
                    add = if ifree == 1
                        ex = -(alp - abs(bn))^2 / (4*alp)
                        bn > 0 && (ex -= bn)
                        exp(ex) / sqrt(4*PI*alp)
                    else
                        _sint(bn, bex, rdbex, sex, nbx, al, wt, tbart, betan, nbeta, maxbb)
                    end
                    snlg += tmp * add
                end

                # Odd j' transitions (swo weight)
                snlk = 0.0
                for lp in 2:2:10
                    jp = lp - 1
                    betap = (-j*(j+1) + jp*(jp+1)) * x_rot / 2
                    tmp = (2*jp + 1) * pj * swo * 4 * _sumh(j, jp, y)
                    bn = be + betap
                    add = if ifree == 1
                        ex = -(alp - abs(bn))^2 / (4*alp)
                        bn > 0 && (ex -= bn)
                        exp(ex) / sqrt(4*PI*alp)
                    else
                        _sint(bn, bex, rdbex, sex, nbx, al, wt, tbart, betan, nbeta, maxbb)
                    end
                    snlk += tmp * add
                end

                sn += snlg + snlk
                l += 2
            end

            # Store: negative β → ssm, positive → ssp
            jj <= nbeta && (ssm[k, nal, itemp] = sn)
            jj >= nbeta && (ssp[k, nal, itemp] = sn)
        end
    end
    return nothing
end

