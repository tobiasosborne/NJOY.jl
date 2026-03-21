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

# ── ETOX / Bondarenko integration ─────────────────────────────────
"""
    ajku(beta, sti) -> (xj, xk)

J and K integral contributions for one set of sampled widths.
Numerical quadrature of the Bondarenko kernel.

- `beta = sigma_m / sigma_peak`
- `sti  = Gamma_total / Delta_Doppler`
"""
function ajku(beta::Float64, sti::Float64)
    xj = 0.0
    xk = 0.0
    nq = 40
    du = 12.0 / nq
    for iq in 1:nq
        u = (iq - 0.5) * du
        phi = sti / (u * u + sti * sti)
        denom = 1.0 + beta * phi
        xj += phi / denom * du
        xk += phi / (denom * denom) * du
    end
    xj /= Float64(pi)
    xk /= Float64(pi)
    return xj, xk
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

    # Pass 2: Hwang quadrature
    nqp = 10
    sigu = zeros(5, nsigz)
    for seq in model.sequences
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

        for is0 in 1:nsigz
            tk_part[is0] /= d
            ttj = 0.0
            for kx in 1:4
                t_part[kx, is0] /= d
                ttj += t_part[kx, is0]
            end
            denom = 1.0 - ttj
            if abs(denom) < 1e-30; denom = 1e-30; end
            sigu[3, is0] += sigm[is0] * t_part[1, is0] / denom  # fission
            sigu[4, is0] += sigm[is0] * t_part[2, is0] / denom  # capture
            sigu[2, is0] += sigm[is0] * t_part[3, is0] / denom  # elastic
            sigu[1, is0] += sigm[is0] * ttj / denom             # total (URR part)
            # transport
            denom_k = 1.0 - tk_part[is0]
            if abs(denom_k) < 1e-30; denom_k = 1e-30; end
            fact = sigma0_values[is0]
            sigu[5, is0] += sigbt * denom / denom_k +
                            fact * (tk_part[is0] - ttj) / denom_k
        end
    end

    # Add backgrounds and potential scattering
    for is0 in 1:nsigz
        sigu[1, is0] += bkg[1]    # total includes all backgrounds
        sigu[2, is0] += bkg[2] + spot + sint  # elastic
        sigu[3, is0] += bkg[3]    # fission
        sigu[4, is0] += bkg[4]    # capture
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
