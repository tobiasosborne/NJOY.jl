# Reich-Moore (LRF=3) cross section evaluation
#
# Direct translation of NJOY2016 reconr.f90 subroutine csrmat.
#
# Physics:
#   Build R-matrix for each J-value from resonance parameters.
#   3-channel formulation: neutron channel + 2 fission channels.
#   Compute collision matrix U from R-matrix.
#   Extract cross sections from U via the optical theorem.
#
# For non-fissile materials (all GFA=GFB=0), this simplifies to
# a 1-channel R-function calculation (no matrix inversion needed).
#
# Reference: NJOY2016 reconr.f90 subroutines csrmat, frobns
# Reference: ENDF-102 Section 2.4 (Reich-Moore formalism)

"""
    cross_section_rm(E, params::ReichMooreParameters,
                     range::ResonanceRange) -> CrossSections

Compute Reich-Moore (LRF=3) cross sections at energy E [eV].

Follows NJOY2016 csrmat exactly (T=0 only).

The calculation builds the R-matrix for each (l, J) combination:
- For non-fissile materials: uses 1-channel R-function (scalar)
- For fissile materials: uses 3x3 matrix (neutron + 2 fission channels)
  and inverts via Julia's built-in matrix inverse

Returns a CrossSections struct with total, elastic, fission, capture [barns].
"""
function cross_section_rm(E::Real, params::ReichMooreParameters,
                          range::ResonanceRange)
    cwaven = cwaven_constant()
    C = PhysicsConstants
    NAPS = range.NAPS

    sig_total = zero(float(E))
    sig_elastic = zero(float(E))
    sig_fission = zero(float(E))
    sig_capture = zero(float(E))

    spi = params.SPI
    ap = params.AP
    gjd = 2.0 * (2.0 * spi + 1.0)

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

    # Loop over l states
    for il in 1:Int(params.NLS)
        ll = Int(params.l_values[il])
        nrs = length(params.Er[il])
        if nrs == 0
            continue
        end

        apl = params.APL[il]

        # Determine channel radius and rhoc for this l-value
        # Following csrmat exactly:
        rhoc_l = k * ap
        rho_l = k * ra
        if apl != 0.0
            rhoc_l = k * apl
        end
        if apl != 0.0 && NAPS == 1
            rho_l = k * apl
        end

        # Shift and penetrability at cross section energy
        se = shift_factor(ll, rho_l)
        pe = penetrability(ll, rho_l)
        phi = phase_shift(ll, rhoc_l)

        # Phase shift factors
        p1 = cos(2.0 * phi)  # cos(2*phi)
        p2 = sin(2.0 * phi)  # sin(2*phi)

        # Loop over possible J values
        fl = Float64(ll)
        ajmin = abs(abs(spi - fl) - 0.5)
        ajmax = spi + fl + 0.5
        numj = round(Int, ajmax - ajmin + 1.0)

        # Determine jjl for the extra channel spin logic
        if ll != 0 && (fl > spi - 0.5 && fl <= spi)
            jjl = 0
        else
            jjl = 1
        end

        ajc = ajmin - 1.0
        for jj in 1:numj
            ajc += 1.0
            gj = (2.0 * ajc + 1.0) / gjd

            # Loop over channel spins (up to 2)
            kchanl = 0
            idone = false
            while kchanl < 2 && !idone
                kchanl += 1

                # Count positive and negative AJ resonances for this J
                kpstv = 0
                kngtv = 0
                for ir in 1:nrs
                    aj_val = params.AJ[il][ir]
                    if abs(abs(aj_val) - ajc) <= 0.25
                        if aj_val < 0.0
                            kngtv += 1
                        end
                        if aj_val > 0.0
                            kpstv += 1
                        end
                    end
                end

                # Initialize R-matrix scalar accumulators (no mutation, AD-safe)
                zE = zero(float(E))
                r11 = zE; r12 = zE; r13 = zE
                r22 = zE; r23 = zE; r33 = zE
                s11 = zE; s12 = zE; s13 = zE
                s22 = zE; s23 = zE; s33 = zE

                # Accumulate R-matrix from resonances
                has_fission = false
                for ir in 1:nrs
                    aj_val = params.AJ[il][ir]
                    aj_abs = abs(aj_val)

                    # Select only resonances with current J value
                    if abs(aj_abs - ajc) > 0.25
                        continue
                    end

                    # Channel spin selection: skip based on sign of AJ
                    iskip = false
                    if kchanl == 1 && aj_val < 0.0
                        iskip = true
                    end
                    if kchanl == 2 && aj_val > 0.0
                        iskip = true
                    end
                    if iskip
                        continue
                    end

                    # Retrieve resonance parameters
                    er = params.Er[il][ir]
                    gn = params.Gn[il][ir]
                    gg = params.Gg[il][ir]
                    gfa = params.Gfa[il][ir]
                    gfb = params.Gfb[il][ir]

                    # Penetrability at resonance energy
                    rho_r = cwaven * arat * sqrt(abs(er)) * ra
                    if apl != 0.0 && NAPS == 1
                        rho_r = cwaven * arat * sqrt(abs(er)) * apl
                    end
                    per = penetrability(ll, rho_r)
                    if per == 0.0
                        continue
                    end

                    # Reduced width amplitudes
                    a1 = sqrt(gn * pe / per)
                    a2 = 0.0
                    if gfa != 0.0
                        a2 = sqrt(abs(gfa))
                        if gfa < 0.0
                            a2 = -a2
                        end
                    end
                    a3 = 0.0
                    if gfb != 0.0
                        a3 = sqrt(abs(gfb))
                        if gfb < 0.0
                            a3 = -a3
                        end
                    end

                    # Energy factors
                    diff = er - E
                    den_val = diff * diff + 0.25 * gg * gg
                    if den_val == 0.0
                        continue
                    end
                    de2 = 0.5 * diff / den_val
                    gg4 = 0.25 * gg / den_val

                    # Accumulate R-matrix (upper triangular, scalars)
                    r11 += gg4 * a1 * a1
                    s11 -= de2 * a1 * a1

                    if gfa != 0.0 || gfb != 0.0
                        r12 += gg4 * a1 * a2
                        s12 -= de2 * a1 * a2
                        r13 += gg4 * a1 * a3
                        s13 -= de2 * a1 * a3
                        r22 += gg4 * a2 * a2
                        s22 -= de2 * a2 * a2
                        r33 += gg4 * a3 * a3
                        s33 -= de2 * a3 * a3
                        r23 += gg4 * a2 * a3
                        s23 -= de2 * a2 * a3
                        has_fission = true
                    end
                end

                # Determine kkkkkk: whether to add contribution
                # kkkkkk = 0: skip
                # kkkkkk = 1: add resonance contribution, no extra hard-sphere
                # kkkkkk = 2: add resonance + hard-sphere phase shift contribution
                kkkkkk = _rm_kkkkkk(kchanl, kpstv, kngtv, jj, jjl, numj)

                if kkkkkk != 0
                    termt = 0.0
                    termn = 0.0
                    termf = 0.0

                    if has_fission
                        # R-matrix path: Frobenius-Schur inversion matching
                        # Fortran csrmat lines 3429-3440: add identity, symmetrize,
                        # then invert via frobns/thrinv/abcmat (reconr.f90:3503-3607)
                        rmat = MMatrix{3,3,Float64}(
                            r11 + 1.0, r12, r13,
                            r12, r22 + 1.0, r23,
                            r13, r23, r33 + 1.0
                        )
                        smat = MMatrix{3,3,Float64}(
                            s11, s12, s13,
                            s12, s22, s23,
                            s13, s23, s33
                        )
                        ri, si = _frobns(rmat, smat)

                        ri11 = ri[1, 1]
                        si11 = si[1, 1]

                        # Fission term
                        t1 = ri[1, 2]
                        t2 = si[1, 2]
                        t3 = ri[1, 3]
                        t4 = si[1, 3]
                        termf = 4.0 * gj * (t1^2 + t2^2 + t3^2 + t4^2)

                        # U11 = exp(2i*phi) * (2*inv - I)
                        u11r = p1 * (2.0 * ri11 - 1.0) + 2.0 * p2 * si11
                        u11i = p2 * (1.0 - 2.0 * ri11) + 2.0 * p1 * si11

                        termt = 2.0 * gj * (1.0 - u11r)
                        termn = gj * ((1.0 - u11r)^2 + u11i^2)
                    else
                        # R-function path (scalar, no fission)
                        dd = r11
                        rr = 1.0 + dd
                        ss = s11
                        amag = rr^2 + ss^2

                        rri = rr / amag
                        ssi = -ss / amag

                        uur = p1 * (2.0 * rri - 1.0) + 2.0 * p2 * ssi
                        uui = p2 * (1.0 - 2.0 * rri) + 2.0 * p1 * ssi

                        small = 3.0e-4
                        if abs(dd) < small && abs(phi) < small
                            xx = 2.0 * dd
                            xx += 2.0 * (dd^2 + ss^2 + phi^2 + p2 * ss)
                            xx -= 2.0 * phi^2 * (dd^2 + ss^2)
                            xx /= amag
                            termt = 2.0 * gj * xx
                            termn = gj * (xx^2 + uui^2)
                        else
                            termt = 2.0 * gj * (1.0 - uur)
                            termn = gj * ((1.0 - uur)^2 + uui^2)
                        end
                        termf = 0.0
                    end

                    # Add extra hard-sphere phase shift for kkkkkk=2
                    if kkkkkk == 2
                        termn += 2.0 * gj * (1.0 - p1)
                        termt += 2.0 * gj * (1.0 - p1)
                    end

                    # Cross section contributions
                    termg = termt - termf - termn
                    sig_elastic += termn
                    sig_capture += termg
                    sig_fission += termf
                    sig_total += termt
                end

                # Check if we need the second channel spin
                if kchanl == 1 && kngtv == 0
                    idone = true
                end
                if kchanl == 2
                    idone = true
                end
            end
        end
    end

    # Final cross sections: multiply by pi/k^2
    sig_total *= pifac
    sig_elastic *= pifac
    sig_fission *= pifac
    sig_capture *= pifac

    return CrossSections(sig_total, sig_elastic, sig_fission, sig_capture)
end

"""
Determine the kkkkkk flag for the Reich-Moore channel spin logic.
This is a direct translation of the complex if/else block in csrmat.

Returns:
  0 = do not add anything
  1 = add resonance contribution but not extra hard-sphere
  2 = add resonance + hard-sphere phase shift contribution
"""
function _rm_kkkkkk(kchanl, kpstv, kngtv, jj, jjl, numj)
    if kchanl == 1
        if kpstv > 0
            if kngtv == 0
                if jj > jjl && jj < numj
                    return 2
                else
                    return 1
                end
            else  # kngtv > 0
                return 1
            end
        else  # kpstv == 0
            if kngtv == 0
                if jj > jjl && jj < numj
                    return 2
                else
                    return 1
                end
            else  # kngtv > 0
                return 0
            end
        end
    else  # kchanl == 2
        if kpstv > 0
            if kngtv == 0
                return 0
            else  # kngtv > 0
                return 1
            end
        else  # kpstv == 0
            if kngtv == 0
                return 0
            else  # kngtv > 0
                if jj > jjl && jj < numj
                    return 2
                else
                    return 1
                end
            end
        end
    end
end

# ==========================================================================
# Frobenius-Schur complex matrix inversion matching Fortran frobns/thrinv/abcmat
# (reconr.f90:3503-3607). Uses the exact same sequence of FP operations as
# the Fortran to ensure bit-identical intermediate rounding.
# ==========================================================================

"""
    _abcmat(a, b) -> c

3×3 matrix multiplication c = a*b, matching Fortran abcmat (reconr.f90:3589-3607).
"""
function _abcmat(a::MMatrix{3,3,Float64}, b::MMatrix{3,3,Float64})
    c = MMatrix{3,3,Float64}(undef)
    for i in 1:3
        for j in 1:3
            c[i,j] = 0.0
            for k in 1:3
                c[i,j] = c[i,j] + a[i,k] * b[k,j]
            end
        end
    end
    return c
end

"""
    _thrinv!(d, n) -> ind

Invert symmetric matrix d in-place, matching Fortran thrinv (reconr.f90:3539-3587).
The Fortran thrinv first transforms d → (I - d), then inverts via Gaussian elimination.
Returns ind: 0 = success, 1 = singular.
"""
function _thrinv!(d::MMatrix{3,3,Float64}, n::Int)
    # Negate and symmetrize, add identity to diagonal
    for j in 1:n
        for i in 1:j
            d[i,j] = -d[i,j]
            d[j,i] = d[i,j]
        end
        d[j,j] = 1.0 + d[j,j]
    end
    # Gaussian elimination
    s = MVector{3,Float64}(undef)
    for lr in 1:n
        fooey = 1.0 - d[lr,lr]
        if fooey == 0.0
            return 1  # singular
        end
        d[lr,lr] = 1.0 / fooey
        for j in 1:n
            s[j] = d[lr,j]
            if j != lr
                d[j,lr] = d[j,lr] * d[lr,lr]
                d[lr,j] = d[j,lr]
            end
        end
        for j in 1:n
            if j != lr
                for i in 1:j
                    if i != lr
                        d[i,j] = d[i,j] + d[i,lr] * s[j]
                        d[j,i] = d[i,j]
                    end
                end
            end
        end
    end
    return 0
end

"""
    _frobns(a, b) -> (c, d)

Invert complex matrix (a + i*b) using the Frobenius-Schur method,
matching Fortran frobns (reconr.f90:3503-3537). Returns real part c
and imaginary part d of the inverse.
"""
function _frobns(a::MMatrix{3,3,Float64}, b::MMatrix{3,3,Float64})
    c = copy(a)
    ind = _thrinv!(a, 3)
    d = MMatrix{3,3,Float64}(undef)
    if ind != 1
        q = _abcmat(a, b)
        d = _abcmat(b, q)
        for i in 1:3
            for j in 1:3
                c[i,j] = c[i,j] + d[i,j]
            end
        end
        _thrinv!(c, 3)
        d = _abcmat(q, c)
        for i in 1:3
            for j in 1:3
                d[i,j] = -d[i,j]
            end
        end
    end
    return c, d
end
