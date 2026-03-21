# Multi-Level Breit-Wigner (MLBW) cross section evaluation and dispatch
#
# Direct translation of NJOY2016 reconr.f90 subroutine:
#   csmlbw  -- Multi-Level Breit-Wigner (LRF=2), T=0 only
#
# The key difference from SLBW is the interference between resonances of the
# same J and l quantum numbers.

# ============================================================================
# MLBW cross section (csmlbw)
# ============================================================================

"""
    cross_section_mlbw(E, params::MLBWParameters, range::ResonanceRange) -> CrossSections

Compute Multi-Level Breit-Wigner cross sections at energy E [eV].

Follows NJOY2016 csmlbw exactly (T=0 only).
The key difference from SLBW is the interference between resonances of the
same J and l quantum numbers. The elastic cross section is computed as:
  sigma_el = pi/k^2 * sum_J g_J * |1 - U_J|^2
where U_J incorporates coherent sum over resonances with that J value.

Returns a CrossSections struct with total, elastic, fission, capture [barns].
"""
function cross_section_mlbw(E::Real, params::MLBWParameters,
                            range::ResonanceRange)
    E = Float64(E)
    cwaven = cwaven_constant()
    C = PhysicsConstants
    NAPS = range.NAPS

    sig_elastic = 0.0
    sig_fission = 0.0
    sig_capture = 0.0
    sig_competitive = 0.0

    spi = params.SPI
    ap = params.AP
    den = 4.0 * spi + 2.0  # denominator for statistical factor

    # Use AWRI from first l-group
    awri = length(params.AWRI) > 0 ? params.AWRI[1] : 1.0

    # Channel radius
    ra = channel_radius(awri)
    if range.NRO == 0
        if NAPS == 1
            ra = ap
        end
    else
        # Energy-dependent scattering radius (NRO != 0)
        ap, ra = _apply_nro(ap, ra, E, range)
    end

    # Wavenumber at energy E
    arat = awri / (awri + 1.0)
    k = cwaven * arat * sqrt(abs(E))
    pifac = C.pi / (k * k)
    rho = k * ra
    rhoc = k * ap

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
            if ll == 0; lp = 2; end
            if ll == 2; lp = 0; end
            _sec, pec = shift_factor(lp, rhop), penetrability(lp, rhop)
            if E + qx / arat < 0.0
                pec = 0.0
            end
        end

        # Phase shift
        phi = phase_shift(ll, rhoc)
        cos2p = 1.0 - cos(2.0 * phi)  # NOTE: csmlbw uses (1 - cos(2*phi))
        sin2p = sin(2.0 * phi)

        # Determine possible J values for this l
        fl = Float64(ll)
        ajmin = abs(abs(spi - fl) - 0.5)
        ajmax = spi + fl + 0.5
        nj = round(Int, ajmax - ajmin + 1.0)

        # Compute statistical spin factors g_J
        gj_arr = zeros(nj)
        sum_gj = 0.0
        aj = ajmin
        for ij in 1:nj
            gj_arr[ij] = (2.0 * aj + 1.0) / den
            sum_gj += gj_arr[ij]
            aj += 1.0
        end
        diff = 2.0 * fl + 1.0 - sum_gj

        # Initialize per-J interference arrays
        sigj = zeros(nj, 2)

        # Loop over resonances
        for ir in 1:nrs
            er = params.Er[il][ir]
            aj_res = params.AJ[il][ir]
            gn = params.Gn[il][ir]
            gg = params.Gg[il][ir]
            gf = params.Gf[il][ir]

            # Map J value to index
            j = round(Int, aj_res - ajmin + 1.0)
            if j < 1 || j > nj
                continue
            end

            # Penetrability and shift at resonance energy
            rho_r = cwaven * arat * sqrt(abs(er)) * ra
            ser = shift_factor(ll, rho_r)
            per = penetrability(ll, rho_r)
            if per == 0.0
                continue
            end
            rper = 1.0 / per

            # Competitive width
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

            # Shifted resonance energy
            erp = er + gn * (ser - se) * rper / 2.0
            edelt = E - erp

            # Energy-dependent neutron width
            gne = gn * pe * rper

            # Non-neutron partial widths
            gx = gg + gf

            # Total width
            gtt = gne + gx
            if gc != 0.0 && pex != 0.0
                gtt = gtt + gc * pec / pex
            end

            # Reduced quantities
            x = 2.0 * edelt / gtt
            comfac = 2.0 * gne / gtt / (1.0 + x * x)

            # Accumulate interference terms per J
            sigj[j, 1] += comfac
            sigj[j, 2] += comfac * x

            # Fission and capture (proportional, no interference)
            comfac_partial = comfac * gj_arr[j] / gtt
            sig_fission += comfac_partial * gf
            sig_capture += comfac_partial * gg
            if lrx != 0 && gc != 0.0
                sig_competitive += comfac_partial * gc * pec / pex
            end
        end

        # Elastic scattering: coherent sum per J
        for j in 1:nj
            add = gj_arr[j] * ((cos2p - sigj[j, 1])^2 + (sin2p + sigj[j, 2])^2)
            sig_elastic += add
        end
        sig_elastic += 2.0 * diff * cos2p
    end

    # Multiply by pi/k^2 factors
    sig_elastic *= pifac
    sig_fission *= 2.0 * pifac
    sig_capture *= 2.0 * pifac

    # Total = elastic + fission + capture + competitive
    sig_total = sig_elastic + sig_fission + sig_capture + sig_competitive
    return CrossSections(sig_total, sig_elastic, sig_fission, sig_capture)
end

# ============================================================================
# Dispatch wrapper
# ============================================================================

"""
    cross_section(E, range::ResonanceRange; temperature=0.0, table=nothing) -> CrossSections

Compute resonance cross sections at energy E [eV] for the given resonance range.
Uses multiple dispatch on the parameter type for type stability.
"""
function cross_section(E::Real, range::ResonanceRange{SLBWParameters};
                       temperature::Real=0.0,
                       table::Union{Nothing,FaddeevaTable}=nothing)
    return cross_section_slbw(E, range.parameters, range;
                              temperature=temperature, table=table)
end

function cross_section(E::Real, range::ResonanceRange{MLBWParameters};
                       temperature::Real=0.0,
                       table::Union{Nothing,FaddeevaTable}=nothing)
    return cross_section_mlbw(E, range.parameters, range)
end

function cross_section(E::Real, range::ResonanceRange{ReichMooreParameters};
                       temperature::Real=0.0,
                       table::Union{Nothing,FaddeevaTable}=nothing)
    return cross_section_rm(E, range.parameters, range)
end

function cross_section(E::Real, range::ResonanceRange{SAMMYParameters};
                       temperature::Real=0.0,
                       table::Union{Nothing,FaddeevaTable}=nothing)
    return cross_section_sammy(E, range.parameters, range)
end

function cross_section(E::Real, range::ResonanceRange;
                       temperature::Real=0.0,
                       table::Union{Nothing,FaddeevaTable}=nothing)
    error("Unsupported resonance formalism: $(typeof(range.parameters))")
end
