# Single-Level Breit-Wigner (SLBW) cross section evaluation
#
# Direct translation of NJOY2016 reconr.f90 subroutine:
#   csslbw  -- Single-Level Breit-Wigner (LRF=1)
#
# Physics:
#   k = cwaven * AWRI/(AWRI+1) * sqrt(E)   [wavenumber]
#   rho = k * a                              [channel radius parameter]
#   P_l(rho), S_l(rho) from penetrability.jl
#   phi_l(rho_cap) from penetrability.jl (phase_shift)
#
# Wavenumber constant cwaven = sqrt(2 * m_n * amu * eV) * 1e-12 / hbar
#
# Channel radius: a = rc1 * (A_target * m_n)^{1/3} + rc2  [fm]
#   rc1 = 0.123, rc2 = 0.08 (matching NJOY2016)
#   If NAPS=1, use AP from ENDF instead.
#
# References:
#   ENDF-102 Chapter 2
#   NJOY2016 reconr.f90 subroutine csslbw

# Channel radius constants (matching NJOY2016)
const _RC1 = 0.123
const _RC2 = 0.08
const _THIRD = 1.0 / 3.0

"""
    cwaven_constant() -> Float64

Compute the wavenumber constant cwaven = sqrt(2 * m_n * amu * eV) * 1e-12 / hbar.
This matches NJOY2016's formula exactly. Units work out to give k = cwaven * arat * sqrt(E)
in inverse femtometers when E is in eV.
"""
function cwaven_constant()
    C = PhysicsConstants
    return sqrt(2.0 * C.amassn * C.amu * C.ev) * 1.0e-12 / C.hbar
end

"""
    channel_radius(awri::Float64) -> Float64

Compute the hard-sphere channel radius a = rc1 * (m_n * AWRI)^{1/3} + rc2.
Result in femtometers, matching NJOY2016.
"""
function channel_radius(awri::Float64)
    aw = PhysicsConstants.amassn * awri
    return _RC1 * aw^_THIRD + _RC2
end

# ============================================================================
# SLBW cross section (csslbw)
# ============================================================================

"""
    cross_section_slbw(E, params::SLBWParameters, range::ResonanceRange;
                       temperature=0.0, table=nothing) -> CrossSections

Compute Single-Level Breit-Wigner cross sections at energy E [eV].

Follows NJOY2016 csslbw exactly:
- For T=0: uses analytic line shapes (1/(1+x^2) form)
- For T>0: uses Doppler-broadened psi/chi line shapes via Faddeeva function

Returns a CrossSections struct with total, elastic, fission, capture [barns].
"""
function cross_section_slbw(E::Real, params::SLBWParameters,
                            range::ResonanceRange;
                            temperature::Real=0.0,
                            table::Union{Nothing,FaddeevaTable}=nothing)
    E = Float64(E)
    T = Float64(temperature)
    cwaven = cwaven_constant()
    C = PhysicsConstants
    NAPS = range.NAPS

    sig_elastic = 0.0
    sig_fission = 0.0
    sig_capture = 0.0
    sig_competitive = 0.0

    spi = params.SPI
    ap = params.AP
    spifac = 1.0 / (2.0 * spi + 1.0)

    # Use AWRI from first l-group (matching NJOY which reads from first l-group)
    awri = length(params.AWRI) > 0 ? params.AWRI[1] : 1.0

    # Channel radius
    ra = channel_radius(awri)
    if NAPS == 1
        ra = ap
    end

    # Wavenumber and related quantities at energy E
    arat = awri / (awri + 1.0)
    k = cwaven * arat * sqrt(abs(E))
    pifac = C.pi / (k * k)
    rho = k * ra
    rhoc = k * ap

    # Boltzmann constant * temperature for Doppler broadening
    tbk = T * C.bk

    # Loop over l states
    for il in 1:Int(params.NLS)
        ll = Int(params.l_values[il])
        nrs = length(params.Er[il])
        qx = params.QX[il]
        lrx = Int(params.LRX[il])

        # Shift and penetrability at cross section energy
        se, pe = shift_factor(ll, rho), penetrability(ll, rho)

        # Competitive reaction penetrability
        pec = 0.0
        if lrx != 0
            rhop = cwaven * arat * sqrt(abs(E + qx / arat)) * ra
            lp = ll
            # Set competing l' value for a 2+ inelastic residual (U-238 convention)
            if ll == 0; lp = 2; end
            if ll == 2; lp = 0; end
            _sec, pec = shift_factor(lp, rhop), penetrability(lp, rhop)
            if E + qx / arat < 0.0
                pec = 0.0
            end
        end

        # Phase shift at cross section energy
        phi = phase_shift(ll, rhoc)
        cos2p = cos(2.0 * phi)
        sin2p = sin(2.0 * phi)
        sinsq = sin(phi)^2

        # Potential scattering for this l
        spot = 4.0 * (2 * ll + 1) * pifac * sinsq

        # Precompute shift and penetrability at each resonance energy
        for ir in 1:nrs
            er = params.Er[il][ir]
            aj = params.AJ[il][ir]
            gn = params.Gn[il][ir]
            gg = params.Gg[il][ir]
            gf = params.Gf[il][ir]

            # Penetrability and shift at resonance energy
            rho_r = cwaven * arat * sqrt(abs(er)) * ra
            ser = shift_factor(ll, rho_r)
            per = penetrability(ll, rho_r)
            if per == 0.0
                continue  # Skip resonances with zero penetrability
            end
            rper = 1.0 / per

            # Competitive width at resonance energy
            gc = 0.0
            pex = 0.0
            if lrx != 0
                gx_total = params.Gx[il][ir]
                if gx_total > 0.0
                    rhoc_r = cwaven * arat * sqrt(abs(er + qx / arat)) * ra
                    lp = ll
                    if ll == 0; lp = 2; end
                    if ll == 2; lp = 0; end
                    _sec_r, pex = shift_factor(lp, rhoc_r), penetrability(lp, rhoc_r)
                    gc = gx_total
                    if er < -qx / arat
                        gc = 0.0
                    end
                end
            end

            gx = gg + gf  # sum of non-neutron partial widths
            gj = (2.0 * aj + 1.0) * spifac / 2.0

            # Shifted resonance energy
            erp = er + gn * (ser - se) * rper / 2.0

            # Energy-dependent neutron width
            gne = gn * pe * rper

            # Total width
            gtt = gne + gx
            if gc != 0.0 && pex != 0.0
                gtt = gtt + gc * pec / pex
            end

            if T > 0.0
                # Doppler-broadened line shapes
                ex = 2.0 * (E - erp) / gtt
                delta = sqrt(4.0 * tbk * E / awri)
                theta = gtt / delta
                ax = theta * ex / 2.0
                y = theta / 2.0

                # Use Faddeeva function for psi/chi
                if table !== nothing
                    rew, aimw = quickw(ax, y, table)
                else
                    rew, aimw = faddeeva_w(ax, y)
                end
                rpi = sqrt(C.pi)
                psi = rpi * theta * rew / 2.0
                chi = rpi * theta * aimw / 2.0

                smax = 4.0 * pifac * gj * gne / gtt^2
                sig_elastic += smax * ((cos2p * gtt - gx) * psi + sin2p * chi * gtt)
                sig_fission += smax * gf * psi
                sig_capture += smax * gg * psi
                if lrx != 0 && gc != 0.0
                    sig_competitive += smax * gc * pec * psi / pex
                end
            else
                # Zero-temperature line shapes
                edelt = E - erp
                comfac = pifac * gj * gne / (edelt^2 + gtt^2 / 4.0)
                add = comfac * (gne * cos2p - 2.0 * gx * sinsq + 2.0 * edelt * sin2p)
                sig_elastic += add
                sig_fission += comfac * gf
                sig_capture += comfac * gg
                if lrx != 0 && gc != 0.0
                    sig_competitive += comfac * gc * pec / pex
                end
            end
        end

        # Add potential scattering
        sig_elastic += spot
    end

    # Total = elastic + fission + capture + competitive
    sig_total = sig_elastic + sig_fission + sig_capture + sig_competitive
    return CrossSections(sig_total, sig_elastic, sig_fission, sig_capture)
end
