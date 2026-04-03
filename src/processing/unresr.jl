# UNRESR -- Bondarenko self-shielding in the unresolved resonance range
#
# Translation of NJOY2016 unresr.f90 subroutines:
#   unresl  -> bondarenko_xs (ETOX statistical averages)
#   ajku    -> ajku (J and K integral evaluation)
#   uunfac  -> urr_penetrability (penetrability + phase shift for URR)
#
# The Bondarenko method: sigma_eff = <sigma * Phi> / <Phi>
# where Phi = 1/(sigma_t + sigma_0) is the narrow resonance flux.

# ── Hwang quadrature tables ────────────────────────────────────────
# 10-point Gauss-Laguerre-type quadrature for chi-squared-distributed
# widths with 1,2,3,4 degrees of freedom (columns).
const HWANG_QW = [  # weights
    1.1120413e-1   0.033773418  3.3376214e-4   1.7623788e-3
    2.3546798e-1   0.079932171  0.018506108    0.021517749
    2.8440987e-1   0.12835937   0.12309946     0.080979849
    2.2419127e-1   0.17652616   0.29918923     0.18797998
    0.10967668     0.21347043   0.33431475     0.30156335
    0.030493789    0.21154965   0.17766657     0.29616091
    0.0042930874   0.13365186   0.042695894    0.10775649
    2.5827047e-4   0.022630659  4.0760575e-3   2.5171914e-3
    4.9031965e-6   1.6313638e-5 1.1766115e-4   8.9630388e-10
    1.4079206e-8   2.745383e-31 5.0989546e-7   0.0
]

const HWANG_QP = [  # abscissae (scaled chi-squared variable)
    3.0013465e-3   1.3219203e-2  1.0004488e-3   0.013219203
    7.8592886e-2   7.2349624e-2  0.026197629    0.072349624
    4.3282415e-1   0.19089473    0.14427472     0.19089473
    1.3345267      0.39528842    0.44484223     0.39528842
    3.0481846      0.74083443    1.0160615      0.74083443
    5.8263198      1.3498293     1.9421066      1.3498293
    9.9452656      2.5297983     3.3150885      2.5297983
    15.782128      5.2384894     5.2607092      5.2384894
    23.996824      13.821772     7.9989414      13.821772
    36.216208      75.647525     12.072069      75.647525
]

# ── URR Spin-Sequence Model ───────────────────────────────────────
"""
    URRSpinSequence

Average resonance parameters for one (l,J) spin sequence in the URR.
"""
struct URRSpinSequence
    l::Int
    J::Float64
    D::Float64          # mean level spacing [eV]
    GN0::Float64        # reduced neutron width (average) [eV]
    GG::Float64         # gamma width (constant) [eV]
    GF::Float64         # fission width [eV]
    GX::Float64         # competitive width [eV]
    AMUN::Int           # df for neutron width (1--4)
    AMUF::Int           # df for fission width (0 = not sampled)
    AMUX::Int           # df for competitive width (0 = not sampled)
end

"""
    URRStatModel

All spin sequences plus global parameters for one nuclide in the URR.
"""
struct URRStatModel
    SPI::Float64                      # target spin
    AWRI::Float64                     # mass ratio A/m_n
    AP::Float64                       # scattering radius [fm]
    sequences::Vector{URRSpinSequence}
end

# ── Penetrability helpers ─────────────────────────────────────────
"""
    urr_penetrability(l, rho, rhoc) -> (Vl, phi)

Penetrability factor and hard-sphere phase shift for l = 0, 1, 2.
Matches NJOY `uunfac`.
"""
function urr_penetrability(l::Int, rho::Float64, rhoc::Float64)
    if l == 0
        return 1.0, rhoc
    elseif l == 1
        r2 = rho * rho
        return r2 / (1.0 + r2), rhoc - atan(rhoc)
    else  # l >= 2
        r2 = rho * rho
        r4 = r2 * r2
        return r4 / (9.0 + 3.0 * r2 + r4),
               rhoc - atan(3.0 * rhoc / (3.0 - rhoc * rhoc))
    end
end

# ── Complex probability integral w(z) — Faddeeva function ────────
# Matches Fortran unresr.f90 subroutines: uw, uwtab, quikw, ajku.
# The w function table is precomputed once and cached.

"""Compute complex probability integral w(z) = exp(-z²)erfc(-iz).
Matches Fortran `uw` (unresr.f90:1494-1657)."""
function _uw(rez::Float64, aim1::Float64)
    rpi = sqrt(Float64(pi))
    rew = 0.0; aimw = 0.0
    aimz = abs(aim1); abrez = abs(rez)
    (abrez + aimz == 0) && return (1.0, 0.0)
    r2 = rez * rez; ai2 = aimz * aimz
    eps_uw = 1.0e-7

    # Determine computation path (asymptotic vs Taylor)
    use_taylor = false
    if abrez + 1.25 * aimz - 5.0 > 0; use_taylor = false
    elseif abrez + 1.863636 * aimz - 4.1 > 0; use_taylor = false
    elseif r2 + 1.71 * ai2 - 2.89 < 0; use_taylor = true
    elseif r2 + 1.18 * ai2 - 5.76 >= 0; use_taylor = true
    elseif aimz - 1.5 >= 0; use_taylor = false
    else; use_taylor = true
    end

    if !use_taylor && aim1 >= 0
        # Asymptotic series (kw=1, label 370)
        rv = 2 * (r2 - ai2); ak = 4 * rez * aimz; el = ak
        h = 0.0; b = 0.0; a = 0.0
        tempm = 0.0; temel = 0.0; g = 1.0
        c = -1.1283792 * aimz; d = 1.1283792 * rez
        am = rv - 1; aak = 1.0; k = 0
        while true
            ajtemp = 2 * aak; temp4 = (1 - ajtemp) * ajtemp
            ajp = rv - (4 * aak + 1)
            # label 480: recurrence
            tempc = ajp * c + temp4 * a - ak * d
            tempd = ajp * d + temp4 * b + ak * c
            temel_n = ajp * el + temp4 * h + ak * am
            tempm_n = ajp * am + temp4 * g - ak * el
            a = c; b = d; g = am; h = el
            c = tempc; d = tempd; am = tempm_n; el = temel_n
            tempm = tempm_n; temel = temel_n
            # Overflow/underflow scaling
            if abs(tempm) + abs(temel) > 1e15
                c *= 1e-15; d *= 1e-15; am *= 1e-15; el *= 1e-15
                tempc *= 1e-15; tempd *= 1e-15; tempm *= 1e-15; temel *= 1e-15
            elseif abs(tempm) + abs(temel) < 1e-15
                c *= 1e15; d *= 1e15; am *= 1e15; el *= 1e15
                tempc *= 1e15; tempd *= 1e15; tempm *= 1e15; temel *= 1e15
            end
            # label 390: convergence
            aak += 1; k += 1
            pr = rew; pim = aimw
            amagn = tempm^2 + temel^2
            rew = (tempc * tempm + tempd * temel) / amagn
            aimw = (tempm * tempd - temel * tempc) / amagn
            (abs(rew - pr) < eps_uw && abs(aimw - pim) < eps_uw) && break
            k > 200 && break  # safety
        end
    else
        # Taylor series (kw=2, label 420)
        aimz_t = aim1  # use original sign
        temp1 = r2 + ai2; temp2 = 2 * temp1 * temp1
        aj = -(r2 - ai2) / temp2; ak = 2 * rez * aimz_t / temp2
        c = 0.0; b = 0.0; ajsig = 0.0; d = 0.0
        g = 0.0; h = 0.0; el = 0.0; a = 1.0; am = 1.0
        sigp = 1.5
        expon = exp(temp2 * aj)
        expc = expon * cos(temp2 * ak)
        exps = -expon * sin(temp2 * ak)
        sig2p = 2 * sigp
        while true
            aj4sig = 4 * ajsig; aj4sm1 = aj4sig - 1
            temp3 = 1 / (aj4sm1 * (aj4sig + 3))
            tt4 = sig2p * (2 * ajsig - 1)
            denom_tt = aj4sm1 * (aj4sig + 1) * (aj4sig - 3) * aj4sm1
            temp4 = denom_tt == 0 ? 0.0 : tt4 / denom_tt
            ajp = aj + temp3
            # label 480: recurrence
            tempc = ajp * c + temp4 * a - ak * d
            tempd = ajp * d + temp4 * b + ak * c
            temel_n = ajp * el + temp4 * h + ak * am
            tempm_n = ajp * am + temp4 * g - ak * el
            a = c; b = d; g = am; h = el
            c = tempc; d = tempd; am = tempm_n; el = temel_n
            tempm = tempm_n; temel = temel_n
            if abs(tempm) + abs(temel) > 1e15
                c *= 1e-15; d *= 1e-15; am *= 1e-15; el *= 1e-15
                tempc *= 1e-15; tempd *= 1e-15; tempm *= 1e-15; temel *= 1e-15
            elseif abs(tempm) + abs(temel) < 1e-15
                c *= 1e15; d *= 1e15; am *= 1e15; el *= 1e15
                tempc *= 1e15; tempd *= 1e15; tempm *= 1e15; temel *= 1e15
            end
            # label 440: convergence
            ajsig += 1
            temp7 = rpi * (am^2 + el^2)
            ref = (aimz_t * (c * am + d * el) - rez * (am * d - c * el)) / temp7 / temp1
            aimf = (aimz_t * (am * d - c * el) + rez * (c * am + d * el)) / temp7 / temp1
            pr = rew; pim = aimw
            rew = expc - ref; aimw = exps - aimf
            (abs(rew - pr) < eps_uw && abs(aimw - pim) < eps_uw) && break
            sig2p = 2 * ajsig
            ajsig > 200 && break  # safety
        end
    end
    return (rew, aimw)
end

# Precomputed w-function table (62×62 grid, x,y from -0.1 to 6.0)
const _UNRESR_W_TABLE = let
    nx = 62; ny = 62
    dx = 0.1; dy = 0.1; x0 = -0.1; y0 = -0.1
    tr = zeros(nx, ny); ti = zeros(nx, ny)
    for i in 1:nx, j in 1:ny
        xi = x0 + (i - 1) * dx; yj = y0 + (j - 1) * dy
        rw, iw = _uw(xi, yj)
        tr[i, j] = rw; ti[i, j] = iw
    end
    (tr=tr, ti=ti)
end

"""Fast w function lookup. Matches Fortran `quikw` (unresr.f90:1297-1377)."""
function _quikw(ax::Float64, y::Float64)
    rpi = sqrt(Float64(pi))
    aki = ax < 0 ? -1.0 : 1.0
    x = abs(ax)
    test = x * x + y * y

    if test < 36.0
        # Table interpolation
        ii = Int(floor(x * 10)); jj = Int(floor(y * 10))
        i = ii + 2; j = jj + 2; n = j - 1
        p = 10 * x - ii; q = 10 * y - jj
        p2 = p * p; q2 = q * q; pq = p * q
        hp = p / 2; hq = q / 2; hq2 = q2 / 2; hp2 = p2 / 2
        a1 = hq2 - hq; a2 = hp2 - hp
        a3 = 1 + pq - p2 - q2; a4 = hp2 - pq + hp; a5 = hq2 - pq + hq
        tr = _UNRESR_W_TABLE.tr; ti_t = _UNRESR_W_TABLE.ti
        rew = a1 * tr[i, n] + a2 * tr[i-1, j] + a3 * tr[i, j] +
              a4 * tr[i+1, j] + a5 * tr[i, j+1] + pq * tr[i+1, j+1]
        aimw = a1 * ti_t[i, n] + a2 * ti_t[i-1, j] + a3 * ti_t[i, j] +
               a4 * ti_t[i+1, j] + a5 * ti_t[i, j+1] + pq * ti_t[i+1, j+1]
        return (rew, aimw * aki)
    elseif test < 144.0
        a1 = x^2 - y^2; a2 = 2 * x * y; a3 = a2^2
        a4 = a1 - 0.2752551; a5 = a1 - 2.724745
        d1 = 0.5124242 / (a4^2 + a3); d2 = 0.05176536 / (a5^2 + a3)
        rew = d1 * (a2 * x - a4 * y) + d2 * (a2 * x - a5 * y)
        aimw = d1 * (a4 * x + a2 * y) + d2 * (a5 * x + a2 * y)
        return (rew, aimw * aki)
    elseif test < 10000.0
        a1 = (x^2 - y^2) * 2; a2 = 4 * x * y
        a4 = a1 - 1; d1 = 1.1283792 / (a4^2 + a2^2)
        rew = d1 * (a2 * x - a4 * y)
        aimw = d1 * (a4 * x + a2 * y)
        return (rew, aimw * aki)
    else
        a1 = 1 / (rpi * test)
        return (y * a1, x * a1 * aki)
    end
end

# ── ETOX / Bondarenko integration ─────────────────────────────────
"""
    ajku(beta, sti) -> (xj, xk)

J and K integral contributions using Gaussian quadrature and the
complex probability integral (Faddeeva function).
Matches Fortran `ajku` (unresr.f90:1383-1454).
"""
function ajku(beta::Float64, sti::Float64)
    xg = (0.095012510, 0.28160355, 0.45801678, 0.61787624,
          0.75540441, 0.8656312, 0.94457502, 0.98940093)
    wg = (0.18945061, 0.18260342, 0.16915652, 0.14959600,
          0.12462897, 0.09515851, 0.06225352, 0.02715246)
    eps_aj = 0.0002
    pi2 = Float64(pi) / 2; sqpi = sqrt(Float64(pi))

    term = (-sti / 2) < log(floatmin(Float64)) ? 0.0 : exp(-sti / 2)
    ep = (1 - term) / eps_aj - beta

    if ep <= 0
        aj = pi2 / beta
        ak = 2 * aj
    else
        y = 0.0; z = 0.0; yk = 0.0; zk = 0.0
        y1 = sti / 2
        b1 = beta / (sqpi * y1)
        a = sqrt((1 + beta) / beta)
        a2 = a * a
        c = 200 / (a * sti)
        remj = (pi2 - atan(c)) / (beta * a)
        remk = ((1 + a2) * remj - (1 - a2) * c / ((c * c + 1) * beta * a)) / (2 * a2)
        for n in 1:8
            ax = 5 * xg[n] + 5
            rew, _ = _quikw(ax, y1)
            x1 = 1 / (1 + b1 / rew); z1 = x1 * x1 / rew
            ax = -5 * xg[n] + 5
            rew, _ = _quikw(ax, y1)
            x2 = 1 / (1 + b1 / rew); z2 = x2 * x2 / rew
            ax = 45 * xg[n] + 55
            rew, _ = _quikw(ax, y1)
            x3 = 1 / (1 + b1 / rew); z3 = x3 * x3 / rew
            ax = -45 * xg[n] + 55
            rew, _ = _quikw(ax, y1)
            x4 = 1 / (1 + b1 / rew); z4 = x4 * x4 / rew
            y += wg[n] * (x3 + x4)
            yk += wg[n] * (z3 + z4)
            z += wg[n] * (x1 + x2)
            zk += wg[n] * (z1 + z2)
        end
        aj = 5 * (z + 9 * y) / y1 + remj
        ak = aj + 5 * b1 * (zk + 9 * yk) / y1 + remk
    end
    return (aj, ak)
end

"""
    bondarenko_xs(model, E, T, sigma0_values; bkg=(0,0,0,0)) -> Matrix

Compute self-shielded Bondarenko cross sections at energy `E` [eV]
and temperature `T` [K] for background dilutions `sigma0_values`.

Returns a (5, nsigz) matrix: rows are total, elastic, fission, capture, transport.

This is the Julia equivalent of NJOY's `unresl` subroutine.
"""
function bondarenko_xs(model::URRStatModel, E::Float64, T::Float64,
                       sigma0_values::Vector{Float64};
                       bkg::NTuple{4,Float64}=(0.0, 0.0, 0.0, 0.0))
    C = NJOY.PhysicsConstants
    nsigz = length(sigma0_values)

    # Wavenumber: k = cwaven * (A/(A+1)) * sqrt(E)
    cwaven = sqrt(2.0 * C.amassn * C.amu * C.ev) * 1e-12 / C.hbar
    awri = model.AWRI
    rat = awri / (awri + 1.0)
    sqE = sqrt(E)
    k = cwaven * rat * sqE
    ab = 4.0 * Float64(pi) / k^2

    # Channel radius
    aw = awri * C.amassn
    aa = 0.123 * aw^(1.0 / 3.0) + 0.08
    rho = k * aa
    rhoc = k * model.AP

    # Doppler width
    T_eff = max(T, 1.0)
    delta_doppler = 2.0 * sqE * sqrt(T_eff * C.bk / awri)

    # Pass 1: potential scattering + interference
    spot = 0.0; sint = 0.0
    seen_l = Set{Int}()
    for seq in model.sequences
        Vl, ps = urr_penetrability(seq.l, rho, rhoc)
        gnx = seq.GN0 * Vl * sqE * Float64(seq.AMUN)
        gj = (2.0 * seq.J + 1.0) / (4.0 * model.SPI + 2.0)
        if !(seq.l in seen_l)
            push!(seen_l, seq.l)
            spot += ab * (2 * seq.l + 1) * sin(ps)^2
        end
        sint -= ab * Float64(pi) * gj * gnx * sin(ps)^2 / seq.D
    end

    sigbt = bkg[1] + spot + sint
    sigm = [sigbt + s0 for s0 in sigma0_values]

    # Pass 2: Hwang quadrature — accumulate per-sequence tl, tj, tk
    # Matches Fortran unresl lines 1076-1145
    nqp = 10
    ns = length(model.sequences)
    # tl[kx, is0, ks]: per-reaction partial (kx=1..3: fission, capture, elastic)
    tl = zeros(3 * ns, nsigz)
    # tj[ks, is0]: per-sequence total
    tj = zeros(ns, nsigz)
    # tk[ks, is0]: per-sequence transport term
    tk = zeros(ns, nsigz)
    abns = ones(ns)  # abundance per sequence (1.0 for single isotope)

    for (ks, seq) in enumerate(model.sequences)
        Vl, ps = urr_penetrability(seq.l, rho, rhoc)
        gnx = seq.GN0 * Vl * sqE * Float64(seq.AMUN)
        gj = (2.0 * seq.J + 1.0) / (4.0 * model.SPI + 2.0)
        d = seq.D
        mu = seq.AMUF; nu = seq.AMUN; lu = seq.AMUX

        nqf = (mu > 0 && seq.GF > 0.0) ? nqp : 1
        nqn = (nu > 0) ? nqp : 1
        nqx = (lu > 0 && seq.GX > 0.0) ? nqp : 1

        t_part = zeros(4, nsigz)
        tk_part = zeros(nsigz)

        for kf in 1:nqf
            gf = (mu > 0 && seq.GF > 0.0) ? HWANG_QP[kf, mu] * seq.GF : seq.GF
            for kn in 1:nqn
                gn = nu > 0 ? HWANG_QP[kn, nu] * gnx : gnx
                for kl in 1:nqx
                    gx = (lu > 0 && seq.GX > 0.0) ? HWANG_QP[kl, lu] * seq.GX : seq.GX
                    gt = gf + seq.GG + gn + gx
                    s0u = ab * gj * gn / gt
                    sti = gt / delta_doppler
                    for is0 in 1:nsigz
                        beta = sigm[is0] / s0u
                        xj, xk = ajku(beta, sti)
                        wt = 1.0
                        if mu > 0 && seq.GF > 0.0; wt *= HWANG_QW[kf, mu]; end
                        if nu > 0; wt *= HWANG_QW[kn, nu]; end
                        if lu > 0 && seq.GX > 0.0; wt *= HWANG_QW[kl, lu]; end
                        xj *= wt; xk *= wt
                        widths = (gf, seq.GG, gn, gx)
                        for kx in 1:4
                            t_part[kx, is0] += xj * widths[kx]
                        end
                        tk_part[is0] += xk * gt
                    end
                end
            end
        end

        # Normalize and store per-sequence (Fortran lines 1134-1145)
        for is0 in 1:nsigz
            tk[ks, is0] = tk_part[is0] / d
            ttj = 0.0
            for kx in 1:4
                ttj += t_part[kx, is0]
            end
            for kx in 1:3
                tl[kx + (ks-1)*3, is0] = t_part[kx, is0] / d
            end
            tj[ks, is0] = ttj / d
        end
    end

    # Pass 3: Sum over all sequences with cross-sequence correction
    # Matches Fortran unresl lines 1152-1191
    sigu = zeros(5, nsigz)
    for is0 in 1:nsigz
        yy = zeros(3)
        yj = 0.0; yk = 0.0
        for ks in 1:ns
            # xj, xk = sum of OTHER sequences' tj, tk (Fortran lines 1164-1168)
            xj_other = 0.0; xk_other = 0.0
            for ksp in 1:ns
                if ksp != ks
                    xj_other += tj[ksp, is0] * abns[ksp]
                    xk_other += tk[ksp, is0] * abns[ksp]
                end
            end
            # Accumulate with cross-sequence correction (lines 1170-1177)
            for kx in 1:3
                knm1 = kx + (ks - 1) * 3
                yy[kx] += tl[knm1, is0] * (1.0 - xj_other) * abns[ks]
            end
            ttt = tj[ks, is0] * (1.0 - xj_other) * abns[ks]
            yj += ttt
            yk += (tk[ks, is0] - tj[ks, is0]) * (1.0 - xk_other) * abns[ks] + ttt
        end
        # Final XS (Fortran lines 1179-1189)
        denom = 1.0 - yj
        if abs(denom) < 1e-30; denom = 1e-30; end
        sigf1 = sigm[is0] * yy[1] / denom  # fission (gg(1))
        sigf2 = sigm[is0] * yy[2] / denom  # capture (gg(2))
        sigf3 = sigm[is0] * yy[3] / denom  # elastic (gg(3))
        sigf4 = sigf1 + sigf2 + sigf3       # total (sum of partials)
        # Transport (Fortran lines 1184-1189)
        denom_k = 1.0 - yk
        if abs(denom_k) < 1e-30; denom_k = 1e-30; end
        term = denom / denom_k
        term1 = sigma0_values[is0] * (yk - yj) / denom_k
        sigtr = sigbt * term + term1

        sigu[1, is0] = sigf4      # total (URR part only, sigbt added below)
        sigu[2, is0] = sigf3      # elastic
        sigu[3, is0] = sigf1      # fission
        sigu[4, is0] = sigf2      # capture
        sigu[5, is0] = sigtr      # transport
    end

    # Add backgrounds and potential scattering (Fortran lines 1195-1199)
    sigbt_out = bkg[1] + spot + sint
    for is0 in 1:nsigz
        sigu[1, is0] += sigbt_out             # total = URR + sigbt
        sigu[2, is0] += bkg[2] + spot + sint  # elastic
        sigu[3, is0] += bkg[3]                # fission
        sigu[4, is0] += bkg[4]                # capture
    end

    return sigu
end

"""
    infinite_dilution_xs(model, E, T) -> CrossSections

Infinitely-dilute (unshielded) average cross sections in the URR.
"""
function infinite_dilution_xs(model::URRStatModel, E::Float64, T::Float64)
    result = bondarenko_xs(model, E, T, [1.0e10])
    return CrossSections(result[1, 1], result[2, 1], result[3, 1], result[4, 1])
end
