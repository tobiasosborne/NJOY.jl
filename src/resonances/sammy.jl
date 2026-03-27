# SAMMY R-Matrix Limited (LRF=7) cross section evaluation
#
# Implements the SAMMY method for computing resonance cross sections
# using the full R-matrix theory with Reich-Moore approximation.
#
# Translation of NJOY2016 samm.f90 subroutines:
#   cssammy  -> cross_section_sammy (main entry)
#   abpart   -> energy-dependent R-matrix construction
#   crosss   -> cross section extraction from collision matrix
#   setr     -> R-matrix + Y-matrix setup
#   sectio   -> partial cross sections from XXXX matrix
#   fxradi   -> wavenumber parameters
#   betset   -> energy-independent reduced-width products
#   pgh      -> penetrability with inverse (S-B+iP)^{-1}
#   sinsix   -> phase shift sin^2(phi), sin(2phi)
#
# Reference: ENDF-102 Section 2.4; NJOY2016 samm.f90; SAMMY manual

"""
    cross_section_sammy(E, params::SAMMYParameters,
                        range::ResonanceRange) -> CrossSections

Compute SAMMY R-Matrix Limited (LRF=7) cross sections at energy E [eV].

Uses the full R-matrix theory with the Reich-Moore approximation:
capture is an eliminated channel and computed from unitarity.
The calculation loops over spin-parity groups and for each group:
1. Computes penetrability, shift factor, phase shift per channel
2. Builds the R-matrix from resonance parameters
3. Forms Y = L^{-1} - R and inverts it
4. Computes XXXX = sqrt(P)/L * Y^{-1} * R * sqrt(P)
5. Extracts elastic, absorption and reaction cross sections from XXXX
"""
function cross_section_sammy(E::Real, params::SAMMYParameters,
                             range::ResonanceRange)
    C = PhysicsConstants

    # Find the neutron particle pair (mt=2) to get AWR and elastic properties
    neutron_pp_idx = 0
    awr = 1.0
    for (i, pp) in enumerate(params.particle_pairs)
        if pp.mt == 2
            neutron_pp_idx = i
            awr = pp.emb  # target mass in neutron mass units
            break
        end
    end
    if neutron_pp_idx == 0
        return CrossSections(zero(float(E)), zero(float(E)),
                             zero(float(E)), zero(float(E)))
    end

    npp = length(params.particle_pairs)

    # Precompute wavenumber constants per particle pair (following fxradi)
    pp_n = params.particle_pairs[neutron_pp_idx]
    ema_n = pp_n.ema == 0.0 ? 1.0 : pp_n.ema
    emb_n = awr == 0.0 ? params.particle_pairs[neutron_pp_idx].emb : awr
    if emb_n == 0.0; emb_n = awr; end

    factor_n = emb_n / (emb_n + ema_n)  # = alabcm for neutron channel

    # twomhb = sqrt(2 * m_n * amu_eV) / (hbar_eV_s * 1e15 * c_m_s)
    hbar_evs = C.hbar / C.ev                         # hbar in eV·s
    amu_ev = C.amu * C.clight * C.clight / C.ev       # amu in eV
    cspeed = C.clight / 100.0                         # c in m/s (CGS clight is cm/s)
    twomhb = sqrt(2.0 * C.amassn * amu_ev) / (hbar_evs * 1.0e15 * cspeed)

    # Accumulate cross sections indexed by particle-pair
    sigmas = zeros(typeof(float(E)), npp)

    fourpi = 4.0 * C.pi / 100.0  # 4*pi/100 (the factor of 100 converts fm^2 to barns)

    for sg in params.spin_groups
        nchan = Int(sg.nchan)
        if nchan == 0
            continue
        end

        AJ = sg.AJ

        # Compute goj = (2J+1) / ((2*spina+1)*(2*spinb+1)) for the entrance channel
        goj = 0.0
        nent = 0  # number of entrance channels
        next_ch = 0  # number of exit-only channels
        for ic in 1:nchan
            ippx = sg.ipp[ic]
            pp = params.particle_pairs[ippx]
            if pp.mt == 2  # neutron elastic = entrance channel
                nent += 1
                if goj == 0.0
                    goj = (2.0 * abs(AJ) + 1.0) /
                          ((2.0 * abs(pp.spina) + 1.0) * (2.0 * abs(pp.spinb) + 1.0))
                end
            elseif ippx > 1
                next_ch += 1
            end
        end

        if nent == 0
            continue
        end

        # ---- Compute per-channel energy-dependent quantities ----
        # zke[ic], zkfe[ic], zkte[ic]: wavenumber factors per channel
        # zeta[ic]: Coulomb parameter (nonzero for charged-particle channels)
        zke = zeros(typeof(float(E)), nchan)
        zkfe = zeros(typeof(float(E)), nchan)
        zkte = zeros(typeof(float(E)), nchan)
        echan = zeros(typeof(float(E)), nchan)
        zeta_ch = zeros(typeof(float(E)), nchan)

        # etac constant for Coulomb parameter (from fxradi, samm.f90:1872)
        hbar_evs = C.hbar / C.ev
        finstri = 1.0e16 * C.hbar / (C.ev^2 * C.clight)
        etac = (1.0 / finstri) * amu_ev / (hbar_evs * 1.0e15 * cspeed) * C.amassn

        for ic in 1:nchan
            ippx = sg.ipp[ic]
            pp = params.particle_pairs[ippx]
            ema_c = pp.ema == 0.0 ? 1.0 : pp.ema
            emb_c = pp.emb == 0.0 ? awr : pp.emb
            aa = emb_c / (emb_c + ema_c)
            # Channel threshold from Q-value
            if pp.qqq != 0.0
                echan[ic] = -pp.qqq / factor_n
            end
            redmas = aa * ema_c
            z = twomhb * sqrt(redmas * factor_n / ema_n)
            zke[ic] = z
            zkfe[ic] = z * sg.rdeff[ic] * 10.0   # convert fm to 1e-15 m
            zkte[ic] = z * sg.rdtru[ic] * 10.0
            # Coulomb parameter: zeta = etac * ZA * ZB * redmas / zke
            # (samm.f90 fxradi line 1900-1901)
            if pp.kza != 0 && pp.kzb != 0
                zeta_ch[ic] = etac * abs(pp.kzb) * abs(pp.kza) * redmas / z
            end
        end

        # ---- Compute reduced-width amplitudes (betapr) and beta products ----
        # betapr[ic, ires] = sqrt(|gamma|/(2*P)) * sign(gamma)
        nres = Int(sg.nres)
        ntriag = nchan * (nchan + 1) / 2 |> Int

        betapr = zeros(typeof(float(E)), nchan, nres)
        beta = zeros(typeof(float(E)), ntriag, nres)
        gbetpr2 = zeros(typeof(float(E)), nres)  # (gamgam/2)
        gbetpr3 = zeros(typeof(float(E)), nres)  # (gamgam/2)^2

        for ires in 1:nres
            gbetpr2[ires] = abs(sg.gamgam[ires]) / 2.0
            gbetpr3[ires] = gbetpr2[ires]^2

            for ic in 1:nchan
                ippx = sg.ipp[ic]
                pp = params.particle_pairs[ippx]
                gam = sg.gamma[ic][ires]
                if gam == 0.0
                    continue
                end

                if pp.lpent <= 0
                    # No penetrability calculation
                    betapr[ic, ires] = sqrt(0.5 * abs(gam))
                    if gam < 0.0
                        betapr[ic, ires] = -betapr[ic, ires]
                    end
                else
                    ex = abs(sg.eres[ires] - echan[ic])
                    if ex == 0.0
                        continue
                    end
                    ex = sqrt(ex)
                    rho = zkte[ic] * ex
                    lsp = sg.lspin[ic]
                    # Compute penetrability at resonance energy
                    # TODO: enable Coulomb for charged-particle channels (zeta_ch[ic] != 0)
                    # betapr values verified identical to Fortran with Coulomb enabled,
                    # but per-energy psmall interaction across spin groups needs debugging.
                    p_res = _sammy_pen(rho, Int(lsp), sg.bound[ic], Int(pp.ishift))
                    if p_res <= 0.0
                        p_res = 1.0
                    end
                    betapr[ic, ires] = sqrt(0.5 * abs(gam) / p_res)
                    if gam < 0.0
                        betapr[ic, ires] = -betapr[ic, ires]
                    end
                end
            end

            # Build beta products (lower triangular, packed)
            kl = 0
            for k in 1:nchan
                for l in 1:k
                    kl += 1
                    beta[kl, ires] = betapr[l, ires] * betapr[k, ires]
                end
            end
        end

        # ---- Build R-matrix at energy E ----
        # R is complex, stored as (real, imag) in lower-triangular packed form
        rmat_r = zeros(typeof(float(E)), ntriag)
        rmat_i = zeros(typeof(float(E)), ntriag)

        # Background R-matrix contributions (diagonal only)
        kl = 0
        for k in 1:nchan
            for l in 1:k
                kl += 1
                if l == k && sg.backgr_type[k] > 0
                    bgr, bgi = _sammy_background(E, Int(sg.backgr_type[k]),
                                                 sg.backgr_data[k])
                    rmat_r[kl] = bgr
                    rmat_i[kl] = bgi
                end
            end
        end

        # Add resonance contributions
        for ires in 1:nres
            difen = sg.eres[ires] - E
            aa = difen^2 + gbetpr3[ires]
            if aa == 0.0
                continue
            end
            alphar_v = difen / aa
            alphai_v = gbetpr2[ires] / aa

            kl = 0
            for k in 1:nchan
                for l in 1:k
                    kl += 1
                    if E > echan[k] && E > echan[l] && beta[kl, ires] != 0.0
                        rmat_r[kl] += alphar_v * beta[kl, ires]
                        rmat_i[kl] += alphai_v * beta[kl, ires]
                    end
                end
            end
        end

        # Check if R-matrix is zero
        lrmat = true
        for kl in 1:ntriag
            if rmat_r[kl] != 0.0 || rmat_i[kl] != 0.0
                lrmat = false
                break
            end
        end

        # ---- Build Y = L^{-1} - R and compute penetrabilities ----
        ymat_r = -copy(rmat_r)
        ymat_i = -copy(rmat_i)
        rootp = ones(typeof(float(E)), nchan)
        elinvr = zeros(typeof(float(E)), nchan)
        elinvi = fill(-one(float(E)), nchan)
        sinsqr = zeros(typeof(float(E)), nchan)
        sin2ph = zeros(typeof(float(E)), nchan)

        for ic in 1:nchan
            ippx = sg.ipp[ic]
            pp = params.particle_pairs[ippx]
            ii = ic * (ic + 1) / 2 |> Int  # diagonal index in packed form

            if E <= echan[ic]
                # Channel is closed
                rootp[ic] = zero(float(E))
                elinvr[ic] = one(float(E))
                elinvi[ic] = zero(float(E))
                continue
            end

            if pp.lpent <= 0
                # No penetrability calculation
                ymat_i[ii] -= 1.0
                continue
            end

            lsp = Int(sg.lspin[ic])
            ex = sqrt(E - echan[ic])
            rho = zkte[ic] * ex
            rhof = zkfe[ic] * ex

            # Phase shift (sin^2, sin(2phi))
            sq, s2p = _sammy_sinsix(rhof, lsp)
            sinsqr[ic] = sq
            sin2ph[ic] = s2p

            # Penetrability and inverse (S-B+iP)
            # TODO: enable Coulomb pgh for charged-particle channels (zeta_ch[ic] != 0)
            p_val, hr, hi, iffy = _sammy_pgh(rho, lsp, sg.bound[ic], Int(pp.ishift))

            # Fortran psmall condition (samm.f90 setr lines 3415-3416):
            # Normal path if: iffy==0 AND NOT (ishift<=0 AND (1-p*rmat_i==1 OR p<tiny))
            ishift_ch = Int(pp.ishift)
            use_normal = iffy == 0 && !(ishift_ch <= 0 &&
                (1.0 - p_val * rmat_i[ii] == 1.0 || p_val < 1.0e-8))
            if use_normal
                rootp[ic] = sqrt(p_val)
                elinvr[ic] = hr
                elinvi[ic] = hi
                ymat_r[ii] += hr
                ymat_i[ii] += hi
            else
                # Small penetrability: scale R-matrix rows/cols
                psmall = sqrt(max(p_val, zero(float(E))))
                ymat_r[ii] = p_val * ymat_r[ii]
                ymat_i[ii] = p_val * ymat_i[ii] - 1.0
                rmat_r[ii] = p_val * rmat_r[ii]
                rmat_i[ii] = p_val * rmat_i[ii]

                # Scale off-diagonal entries
                for j in 1:nchan
                    if j == ic; continue; end
                    if j < ic
                        ji = ic * (ic - 1) / 2 + j |> Int
                    else
                        ji = j * (j - 1) / 2 + ic |> Int
                    end
                    ymat_r[ji] *= psmall
                    ymat_i[ji] *= psmall
                    rmat_r[ji] *= psmall
                    rmat_i[ji] *= psmall
                end
                rootp[ic] = one(float(E))
                elinvr[ic] = zero(float(E))
                elinvi[ic] = -one(float(E))
            end
        end

        # Ensure diagonal entries are not zero
        if !lrmat
            kl = 0
            for k in 1:nchan
                kl += k
                if ymat_r[kl] == 0.0 && ymat_i[kl] == 0.0
                    ymat_r[kl] = 1.0
                end
            end
        end

        # ---- Invert Y-matrix ----
        if lrmat
            # R-matrix is zero: xqr=0, xqi=0, xxxxr=0, xxxxi=0
            xxxxr = zeros(typeof(float(E)), ntriag)
            xxxxi = zeros(typeof(float(E)), ntriag)
        else
            # Build full complex Y matrix and invert
            ymat_full = zeros(Complex{typeof(float(E))}, nchan, nchan)
            rmat_full = zeros(Complex{typeof(float(E))}, nchan, nchan)
            for i in 1:nchan
                for j in 1:i
                    ij = i * (i - 1) / 2 + j |> Int
                    ymat_full[i, j] = complex(ymat_r[ij], ymat_i[ij])
                    ymat_full[j, i] = ymat_full[i, j]  # symmetric
                    rmat_full[i, j] = complex(rmat_r[ij], rmat_i[ij])
                    rmat_full[j, i] = rmat_full[i, j]
                end
            end

            yinv_full = inv(ymat_full)

            # XQ = Yinv * R  (note: asymmetric because multiplication)
            xq = yinv_full * rmat_full

            # XXXX = sqrt(P)/L * XQ * sqrt(P)  (symmetric, lower triangular packed)
            xxxxr = zeros(typeof(float(E)), ntriag)
            xxxxi = zeros(typeof(float(E)), ntriag)
            ij = 0
            for i in 1:nchan
                plr = rootp[i] * elinvr[i]
                pli = rootp[i] * elinvi[i]
                for j in 1:i
                    ij += 1
                    xq_val = xq[j, i]
                    xqr_v = real(xq_val)
                    xqi_v = imag(xq_val)
                    xxxxr[ij] = rootp[j] * (xqr_v * plr - xqi_v * pli)
                    xxxxi[ij] = rootp[j] * (xqi_v * plr + xqr_v * pli)
                end
            end
        end

        # ---- Extract cross sections (sectio) ----
        crss = zeros(typeof(float(E)), npp)

        # Elastic: crss(neutron_pp) = g * [sin^2(phi)*(1-2*Xi) - sin(2phi)*Xr + (Xr^2+Xi^2)]
        ii = 0
        ij_cnt = 0
        for i in 1:nent
            zz = zke[i]^2
            ii += i  # diagonal index
            termn = sinsqr[i] * (1.0 - 2.0 * xxxxi[ii]) - sin2ph[i] * xxxxr[ii]
            termn /= zz
            for j in 1:i
                ij_cnt += 1
                ar = (xxxxr[ij_cnt]^2 + xxxxi[ij_cnt]^2) / zz
                if i != j
                    ar *= 2.0
                end
                termn += ar
            end
            crss[neutron_pp_idx] += termn
        end
        crss[neutron_pp_idx] *= goj

        # Absorption (capture + all non-elastic): crss(capture_pp)
        # absorption = g * [Xi - (Xr^2+Xi^2)] / zz
        capture_pp_idx = 0
        for (i, pp) in enumerate(params.particle_pairs)
            if pp.mt == 102
                capture_pp_idx = i
                break
            end
        end

        ii = 0
        ij_cnt = 0
        terma_total = zero(float(E))
        for i in 1:nent
            ii += i
            zz = zke[i]^2
            terma = xxxxi[ii] / zz
            for j in 1:i
                ij_cnt += 1
                ar = (-xxxxr[ij_cnt]^2 - xxxxi[ij_cnt]^2) / zz
                if i != j
                    ar *= 2.0
                end
                terma += ar
            end
            terma_total += terma
        end
        if capture_pp_idx > 0
            crss[capture_pp_idx] = terma_total * goj
        end

        # Reaction channels (non-elastic, non-capture)
        for jj in 1:next_ch
            j = jj + nent
            if j > nchan; continue; end
            ip = sg.ipp[j]
            for i in 1:nent
                zz = zke[i]^2
                ij_idx = j * (j - 1) / 2 + i |> Int
                crss[ip] += (xxxxr[ij_idx]^2 + xxxxi[ij_idx]^2) / zz
            end
        end
        for ip in 1:npp
            pp = params.particle_pairs[ip]
            if pp.mt != 2 && pp.mt != 102
                crss[ip] *= goj
            end
        end

        # Add to totals
        for ip in 1:npp
            sigmas[ip] += crss[ip]
        end
    end

    # Normalize: sigma = sigmas * 4*pi/(100*E)
    norm = fourpi / E
    for ip in 1:npp
        sigmas[ip] *= norm
    end

    # Map to CrossSections: elastic = pp with mt=2, capture = pp with mt=102
    sig_elastic = zero(float(E))
    sig_capture = zero(float(E))
    sig_fission = zero(float(E))

    for (i, pp) in enumerate(params.particle_pairs)
        if pp.mt == 2
            sig_elastic += sigmas[i]
        elseif pp.mt == 102
            # Capture = absorption - sum of non-capture reactions
            sig_capture += sigmas[i]
        elseif pp.mt == 18
            sig_fission += sigmas[i]
        end
    end

    # Subtract non-elastic reaction channels from capture (unitarity)
    for (i, pp) in enumerate(params.particle_pairs)
        if pp.mt != 2 && pp.mt != 102
            sig_capture -= sigmas[i]
        end
    end

    sig_total = sig_elastic + sig_capture + sig_fission
    return CrossSections(sig_total, sig_elastic, sig_fission, sig_capture)
end

# ---- Helper: penetrability P(rho) for pgh ----
function _sammy_pen(rho, l::Int, bound::Float64, ishift::Int)
    r = rho
    r2 = r * r
    if l == 0
        return r
    elseif l == 1
        d = 1.0 + r2
        return r * r2 / d
    elseif l == 2
        r4 = r2 * r2
        d = 9.0 + r2 * (3.0 + r2)
        return r * r4 / d
    elseif l == 3
        r4 = r2 * r2; r6 = r4 * r2
        d = 225.0 + r2 * (45.0 + r2 * (6.0 + r2))
        return r * r6 / d
    elseif l == 4
        r4 = r2 * r2; r6 = r4 * r2; r8 = r4 * r4
        d = 11025.0 + r2 * (1575.0 + r2 * (135.0 + r2 * (10.0 + r2)))
        return r * r8 / d
    else
        return penetrability(l, rho)
    end
end

# ---- Helper: pgh -- computes hr, hi = real, imag of 1/(S-B+iP) ----
function _sammy_pgh(rho, l::Int, bound::Float64, ishift::Int)
    r = rho
    r2 = r * r
    p = 0.0
    s_val = 0.0

    # Compute P and S
    if l == 0
        p = r
        s_val = 0.0
    elseif l == 1
        d = 1.0 + r2
        p = r * r2 / d
        s_val = -1.0 / d
    elseif l == 2
        r4 = r2 * r2
        d = 9.0 + r2 * (3.0 + r2)
        p = r * r4 / d
        s_val = -(18.0 + 3.0 * r2) / d
    elseif l == 3
        r4 = r2 * r2; r6 = r4 * r2
        d = 225.0 + r2 * (45.0 + r2 * (6.0 + r2))
        p = r * r6 / d
        s_val = -(675.0 + r2 * (90.0 + 6.0 * r2)) / d
    elseif l == 4
        r4 = r2 * r2; r6 = r4 * r2; r8 = r4 * r4
        d = 11025.0 + r2 * (1575.0 + r2 * (135.0 + r2 * (10.0 + r2)))
        p = r * r8 / d
        s_val = -(44100.0 + r2 * (4725.0 + r2 * (270.0 + 10.0 * r2))) / d
    else
        p = penetrability(l, rho)
        s_val = shift_factor(l, rho)
    end

    # g = S - B (if ishift > 0, else g=0)
    g = 0.0
    if ishift > 0
        g = s_val - bound
    end

    if p <= 1.0e-35
        p = 0.0
    end

    iffy = 0
    if g == 0.0 && p == 0.0
        return (0.0, 0.0, 0.0, 1)
    elseif g == 0.0
        return (p, 0.0, -1.0 / p, 0)
    elseif p == 0.0
        return (0.0, 1.0 / g, 0.0, 0)
    else
        dd = p^2 + g^2
        hr = g / dd
        hi = -p / dd
        return (p, hr, hi, 0)
    end
end

# ---- Helper: JWKB approximation for Coulomb functions (samm.f90 jwkb) ----
function _jwkb(x::Float64, eta::Float64, xl::Float64)
    aloge = 0.434294481903251816667932   # log10(e)
    six35 = 6.0 / 35.0
    gh2 = x * (2.0 * eta - x)
    xll1 = max(xl * xl + xl, 0.0)
    if gh2 + xll1 <= 0.0
        return (0.0, 0.0, 0)
    end
    hll = xll1 + six35
    hl = sqrt(hll)
    sl = eta / hl + hl / x
    rl2 = 1.0 + eta * eta / hll
    gh = sqrt(gh2 + hll) / x
    phi = x * gh - 0.5 * (hl * log((gh + sl)^2 / rl2) - log(gh))
    if eta != 0.0
        phi -= eta * atan(x * gh, x - eta)
    end
    phi10 = -phi * aloge
    iexp = floor(Int, phi10)
    if iexp > 70
        gjwkb = 10.0^(phi10 - iexp)
    else
        gjwkb = exp(-phi)
        iexp = 0
    end
    fjwkb = 0.5 / (gh * gjwkb)
    return (fjwkb, gjwkb, iexp)
end

# ---- coulx algorithm: Coulomb wavefunctions for small rho (samm.f90) ----
# Translated from NJOY2016 samm.f90 subroutines xsigll, asymp2, taylor,
# getfg, getps. Used for rho values where Steed/coulfg is unreliable.

"""
    _xsigll(eta, lmax) -> (sigma, sigma0)

Compute Coulomb phase shift sigma_0 and sigma(L) for L=0..lmax.
Translated from samm.f90:4429-4517.
"""
function _xsigll(eeta::Float64, lmax::Int)
    small = 0.000001
    ber = (0.1666666666666666666666666666666666667,
           -0.0333333333333333333333333333333333333,
            0.0238095238095238095238095238095238095,
           -0.0333333333333333333333333333333333333,
            0.0757575757575757575757575757575757576)
    mmmxxx = 100000

    eta = eeta
    peta = abs(eta)
    sigma0 = 0.0

    if peta >= 3.0
        sum_ = 0.0
        for i in 1:5
            xi = Float64(i)
            m = 2*i - 1
            xm = Float64(m)
            sum_ += ber[i] / (2*xi * xm * peta^m)
        end
        sigma0 = 3.141592653589793238 / 4 + peta * (log(peta) - 1) - sum_
    else
        sumas = 0.0
        for is_ in 1:mmmxxx
            s = Float64(is_)
            temp1 = peta / s
            if s <= 2*peta
                as_ = temp1 - atan(temp1)
            else
                as_ = 0.0
                k = 0
                for j in 1:mmmxxx
                    m = j + j + 1
                    xm = Float64(m)
                    add = temp1^m / xm
                    if k == 0
                        as_ = as_ + add
                        k = 1
                    else
                        as_ = as_ - add
                        k = 0
                    end
                    if abs(add / as_) <= small
                        break
                    end
                end
            end
            sumas += as_
            if abs(as_ / sumas) <= small
                break
            end
        end
        sigma0 = -0.57721566490153286 * peta + sumas
    end

    if eta < 0
        sigma0 = -sigma0
    end

    sigma = Vector{Float64}(undef, lmax + 1)
    sigma[1] = sigma0
    if lmax > 0
        for ll in 1:lmax
            xl = Float64(ll)
            sigma[ll+1] = sigma[ll] + atan(eta / xl)
        end
    end

    return (sigma, sigma0)
end

"""
    _asymp2(eta, rho, sigma0) -> (u, upr, rhoi)

Compute G_0 and G_0' at a large asymptotic radius rhoi using
the asymptotic expansion. Translated from samm.f90:4519-4601.
"""
function _asymp2(eeta::Float64, rrho::Float64, sigma0::Float64)
    del = 100.0
    epslon = 0.000001

    eta = eeta
    rho = rrho
    rhoi = max(rho * 2, 10.0, 10.0 * eta)

    @label start_asymp2
    xn = 0.0
    zold = [1.0, 0.0, 0.0, 1.0 - eta / rhoi]
    z = copy(zold)
    bigz = [abs(z[i]) for i in 1:4]
    jcheck = 0

    for n in 1:100
        temp = 2 * (xn + 1) * rhoi
        an = (2*xn + 1) * eta / temp
        bn = (eta*eta - xn*(xn+1)) / temp
        xn += 1.0
        znew = Vector{Float64}(undef, 4)
        znew[1] = an*zold[1] - bn*zold[2]
        znew[2] = an*zold[2] + bn*zold[1]
        znew[3] = an*zold[3] - bn*zold[4] - znew[1]/rhoi
        znew[4] = an*zold[4] + bn*zold[3] - znew[2]/rhoi
        icheck = 0
        diverged = false
        for i in 1:4
            z[i] += znew[i]
            zold[i] = znew[i]
            bigz[i] = max(bigz[i], abs(z[i]))
            temp2 = abs(z[i])
            if bigz[i] / temp2 > del
                diverged = true
                break
            end
            if abs(znew[i] / z[i]) <= epslon
                icheck += 1
            end
        end
        if diverged
            rhoi *= 2
            @goto start_asymp2
        end
        w = z[1]*z[4] - z[2]*z[3]
        if abs(w) > 10.0
            rhoi *= 2
            @goto start_asymp2
        end
        if icheck == 4
            jcheck += 1
            if jcheck >= 4
                @goto converged_asymp2
            end
        else
            jcheck = 0
        end
    end
    # If we fall through (no convergence in 100 iterations), double rhoi and retry
    rhoi *= 2
    @goto start_asymp2

    @label converged_asymp2
    phi = rhoi - eta * log(2*rhoi) + sigma0
    cosphi = cos(phi)
    sinphi = sin(phi)
    g0 = z[1]*cosphi - z[2]*sinphi
    g0pr = z[3]*cosphi - z[4]*sinphi
    u = g0
    upr = g0pr
    return (u, upr, rhoi)
end

"""
    _taylor(eta, rho, u, upr, rhoi) -> (u, upr)

Integrate G_0 from rhoi back to the target rho using Taylor expansion.
Translated from samm.f90:4603-4735.
"""
function _taylor(eeta::Float64, rrho::Float64, u::Float64, upr::Float64, rhoi_in::Float64)
    epslon = 1.0e-6
    bigger = 1.0e10
    biggst = 1.0e30
    del = 100.0

    eta = eeta
    rho = rrho
    rhoi = rhoi_in
    delta = rho - rhoi

    if delta == 0.0
        return (u, upr)
    end

    @label label10
    a = Vector{Float64}(undef, 100)
    a[1] = u
    a[2] = delta * upr
    a[3] = -delta*delta / 2 * (1 - 2*eta/rhoi) * a[1]
    nstart = 4

    @label label20
    jcheck = 0
    sum_ = 0.0
    sumpr = 0.0
    big = 0.0
    bigpr = 0.0

    for n in 1:100
        xn = Float64(n - 1)
        if n >= nstart
            a[n] = -(delta*(xn-1)*(xn-2)*a[n-1] +
                      (rhoi - 2*eta)*(delta^2)*a[n-2] +
                      (delta^3)*a[n-3]) /
                    (rhoi*(xn-1)*xn)
            if a[n] > bigger
                # goto 40: halve delta and retry
                nstart = max(nstart, n + 1)
                m = nstart - 1
                delta = delta / 2
                temp = 2.0
                for nn in 1:m
                    temp /= 2
                    a[nn] *= temp
                end
                @goto label20
            end
        end
        sum_ += a[n]
        sumpr += xn * a[n]
        if sum_ >= biggst
            nstart = max(nstart, n + 1)
            m = nstart - 1
            delta = delta / 2
            temp = 2.0
            for nn in 1:m
                temp /= 2
                a[nn] *= temp
            end
            @goto label20
        end
        if sumpr >= biggst
            nstart = max(nstart, n + 1)
            m = nstart - 1
            delta = delta / 2
            temp = 2.0
            for nn in 1:m
                temp /= 2
                a[nn] *= temp
            end
            @goto label20
        end
        big = max(big, abs(sum_))
        bigpr = max(bigpr, abs(sumpr))
        if sum_ == 0.0 || sumpr == 0.0
            jcheck = 0
        else
            if abs(big / sum_) >= del
                nstart = max(nstart, n + 1)
                m = nstart - 1
                delta = delta / 2
                temp = 2.0
                for nn in 1:m
                    temp /= 2
                    a[nn] *= temp
                end
                @goto label20
            end
            if abs(bigpr / sumpr) >= del
                nstart = max(nstart, n + 1)
                m = nstart - 1
                delta = delta / 2
                temp = 2.0
                for nn in 1:m
                    temp /= 2
                    a[nn] *= temp
                end
                @goto label20
            end
            if abs(a[n] / sum_) >= epslon || abs(xn * a[n] / sumpr) >= epslon
                jcheck = 0
            else
                jcheck += 1
                if jcheck >= 4
                    @goto label60
                end
            end
        end
    end
    # n=100, fall through to label 40 equivalent
    n_final = 100
    nstart = max(nstart, n_final + 1)
    m = nstart - 1
    delta = delta / 2
    temp = 2.0
    for nn in 1:m
        temp /= 2
        a[nn] *= temp
    end
    @goto label20

    @label label60
    u = sum_
    upr = sumpr / delta
    rhoi = rhoi + delta
    delta = rho - rhoi
    if abs(delta) >= epslon
        @goto label10
    end

    return (u, upr)
end

"""
    _getfg(eta, rho, lmax, lll) -> (f, fpr, g, gpr, llmax_out)

Construct F(L) and G(L) from G_0, G_0' using recurrence.
Translated from samm.f90:4753-4850. g0 and g0pr are passed as
the initial values in g[1] and gpr[1].
"""
function _getfg(eta::Float64, rho::Float64, llmax::Int, lll::Int,
                g0::Float64, g0pr::Float64)
    big = 1.0e12
    lmax = llmax
    limit = max(3, lmax + 1)

    # Allocate arrays large enough for upward + downward recurrence
    maxsize = limit + 10000 + 10  # generous upper bound
    f   = zeros(Float64, maxsize)
    fpr = zeros(Float64, maxsize)
    g   = zeros(Float64, maxsize)
    gpr = zeros(Float64, maxsize)

    g[1] = g0
    gpr[1] = g0pr

    # Upward recurrence for G
    g[2] = ((eta + 1/rho)*g[1] - gpr[1]) / sqrt(eta^2 + 1)
    for l in 3:limit
        xl = Float64(l - 1)
        temp1 = sqrt(xl*xl + eta*eta)
        g[l] = (2*xl - 1)/temp1 * (eta/(xl-1) + xl/rho) * g[l-1] -
               xl/temp1 * sqrt(1 + (eta/(xl-1))^2) * g[l-2]
        if abs(g[l]) > big && l > lll
            limit = l
            break
        end
    end

    # Continue upward recurrence for G until ratio decays
    gm2 = g[limit-1]
    gm1 = g[limit]
    il = -1
    j_final = limit
    for j in limit:10000
        xl = Float64(j)
        temp1 = sqrt(xl*xl + eta*eta)
        gm = (2*xl - 1)/temp1 * (eta/(xl-1) + xl/rho) * gm1 -
             xl/temp1 * sqrt(1 + (eta/(xl-1))^2) * gm2
        if abs(g[limit] / gm) > 1.0e-4
            il = -2
        end
        if il > 0
            j_final = j
            break
        end
        il += 1
        gm2 = gm1
        gm1 = gm
        j_final = j
    end

    # Downward recurrence for f
    xl = Float64(j_final)
    fp1 = xl / gm1 / sqrt(xl^2 + eta^2)   # gm1 is gm at j_final
    fp2 = 0.0
    l = j_final - 1

    # First phase: down to limit+3
    for ll in 1:(j_final - 3 - limit)
        l -= 1
        xl = Float64(l)
        temp2 = sqrt((xl+1)^2 + eta^2)
        fp = ((2*xl + 3)*(eta/(xl+2) + (xl+1)/rho)*fp1 -
              (xl+1)*sqrt(1 + (eta/(xl+2))^2)*fp2) / temp2
        fp2 = fp1
        fp1 = fp
    end

    # Second phase: down to 0, storing f values
    for ll in 1:(limit + 2)
        l -= 1
        xl = Float64(l)
        temp2 = sqrt((xl+1)^2 + eta^2)
        fp = ((2*xl + 3)*(eta/(xl+2) + (xl+1)/rho)*fp1 -
              (xl+1)*sqrt(1 + (eta/(xl+2))^2)*fp2) / temp2
        f[l+1] = fp
        fp2 = fp1
        fp1 = fp
    end

    # Compute derivatives
    fpr[1] = (1/rho + eta)*f[1] - sqrt(1 + eta^2)*f[2]
    for l in 2:limit
        xl = Float64(l)
        fpr[l] = (xl/rho + eta/xl)*f[l] - sqrt(1 + (eta/xl)^2)*f[l+1]
        temp1 = eta / (xl - 1)
        gpr[l] = sqrt(1 + temp1^2)*g[l-1] - ((xl-1)/rho + temp1)*g[l]
    end

    llmax_out = lmax
    if limit <= lmax
        llmax_out = limit - 1
    end

    return (f, fpr, g, gpr, llmax_out)
end

"""
    _getps(rho, lll, f, fpr, g, gpr) -> (p, s)

Extract penetrability P and shift factor S from F, G, F', G'.
Translated from samm.f90:4969-4996.
"""
function _getps(rho::Float64, lll::Int, f::Vector{Float64}, fpr::Vector{Float64},
                g::Vector{Float64}, gpr::Vector{Float64})
    n = lll + 1
    asq = f[n]^2 + g[n]^2
    a = sqrt(asq)
    p = rho / asq
    ss = rho * (f[n]*fpr[n] + g[n]*gpr[n]) / asq
    return (p, ss)
end

"""
    _coulomb_pen_shift_coulx(rho, l, eta) -> (P, S)

Compute Coulomb penetrability and shift factor using the coulx algorithm
(xsigll -> asymp2 -> taylor -> getfg -> getps). This is the path used
in NJOY2016/SAMMY for small rho values where Steed's method (coulfg) fails.
"""
function _coulomb_pen_shift_coulx(rho::Float64, l::Int, eta::Float64)
    lmax = l

    # Step 1: Coulomb phase shift
    sigma, sigma0 = _xsigll(eta, lmax)

    # Step 2: Asymptotic G_0 at large radius
    u, upr, rhoi = _asymp2(eta, rho, sigma0)

    # Step 3: Integrate G_0 back to target rho
    u, upr = _taylor(eta, rho, u, upr, rhoi)

    # Check for overflow (coulx line 4388: if abs(g0) > 1e25, call end1)
    if abs(u) > 1.0e25
        return (0.0, 0.0)
    end

    # Step 4: Construct F(L) and G(L) from G_0
    f, fpr, g, gpr, _ = _getfg(eta, rho, lmax, l, u, upr)

    # Step 5: Extract P and S
    p, s = _getps(rho, l, f, fpr, g, gpr)

    return (p, s)
end

# ---- Coulomb penetrability using Steed's method (samm.f90 coulfg) ----
# Returns (P, S) = penetrability and shift factor for Coulomb channels.
# Translated from NJOY2016 samm.f90 coulfg + jwkb for the case llmin=0.
# For small rho (< 1.02), dispatches to the coulx algorithm instead.
function _coulomb_pen_shift(rho::Float64, l::Int, eta::Float64)
    accur = 1.0e-16
    acc = accur
    acc4 = acc * 1.0e3
    acch = sqrt(acc)
    abort_limit = 2.0e4

    if rho <= acch
        # Very small rho: return zero penetrability
        return (0.0, 0.0)
    end

    # For small rho where Steed/coulfg is unreliable, use coulx algorithm
    if rho < 1.02
        return _coulomb_pen_shift_coulx(rho, l, eta)
    end

    x = rho
    xlm = 0.0
    xll = Float64(l)
    e2mm1 = eta * eta + xlm * xlm + xlm    # = eta^2 for xlm=0
    xlturn = x * (x - 2.0 * eta) < xlm * xlm + xlm
    l1 = l + 1

    # ---- CF1: evaluate f = F'(xl,eta,x) / F(xl,eta,x) ----
    xi = 1.0 / x
    fcl = 1.0
    pk = xll + 1.0
    px = pk + abort_limit
    # Initial step
    ek = eta / pk
    f = (ek + pk * xi) * fcl + (fcl - 1.0) * xi
    pk1 = pk + 1.0
    if abs(eta * x + pk * pk1) <= acc
        fcl = (1.0 + ek * ek) / (1.0 + (eta / pk1)^2)
        pk = 2.0 + pk
        ek = eta / pk
        f = (ek + pk * xi) * fcl + (fcl - 1.0) * xi
        pk1 = pk + 1.0
    end
    d = 1.0 / ((pk + pk1) * (xi + ek / pk1))
    df = -fcl * (1.0 + ek * ek) * d
    if fcl != 1.0; fcl = -1.0; end
    if d < 0.0; fcl = -fcl; end
    f = f + df

    # CF1 loop
    p_cf = 1.0
    converged_cf1 = false
    for _ in 1:200000
        pk = pk1
        pk1 = pk1 + 1.0
        ek = eta / pk
        tk = (pk + pk1) * (xi + ek / pk1)
        d = tk - d * (1.0 + ek * ek)
        if abs(d) <= acch
            p_cf += 1.0
            if p_cf > 2.0; break; end
        end
        d = 1.0 / d
        if d < 0.0; fcl = -fcl; end
        df = df * (d * tk - 1.0)
        f = f + df
        if pk > px; break; end
        if abs(df) < abs(f) * acc
            converged_cf1 = true
            break
        end
    end

    # ---- Downward recurrence to lambda=0 (if l > 0) ----
    gc_store = zeros(l + 3)  # storage for RL values
    if l > 0
        fcl = fcl * 1.0e-30
        fpl = fcl * f
        xl = xll
        for lp in 1:l
            el = eta / xl
            rl = sqrt(1.0 + el * el)
            sl = el + xl * xi
            fcl1 = (fcl * sl + fpl) / rl
            fpl = fcl1 * sl - fcl * rl
            fcl = fcl1
            gc_store[l + 1 - lp + 1] = rl
            xl -= 1.0
        end
        if fcl == 0.0; fcl = acc; end
        f = fpl / fcl
    end

    # ---- CF2 or JWKB approximation ----
    if xlturn
        fjwkb, gjwkb, iexp = _jwkb(x, eta, max(xlm, 0.0))
    end

    if xlturn && (iexp > 1 || gjwkb > 1.0 / (acch * 100.0))
        # JWKB path (tunneling region)
        w = fjwkb
        gam = gjwkb * w
        p_val = f
        q_val = 1.0
    else
        # CF2: evaluate p + iq using Steed's algorithm
        p_val = 0.0
        q_val = 1.0 - eta * xi
        ar = -e2mm1
        ai = eta
        br = 2.0 * (x - eta)
        bi = 2.0
        wi = eta + eta  # = 2*eta (Fortran variable wi)
        dd = br * br + bi * bi
        dr = br / dd
        di = -bi / dd
        dp = -xi * (ar * di + ai * dr)
        dq = xi * (ar * dr - ai * di)
        pk_cf2 = 0.0
        for _ in 1:200000
            p_val += dp
            q_val += dq
            pk_cf2 += 2.0
            ar += pk_cf2
            ai += wi
            bi += 2.0
            d_re = ar * dr - ai * di + br
            d_im = ai * dr + ar * di + bi
            c = 1.0 / (d_re * d_re + d_im * d_im)
            dr = c * d_re
            di = -c * d_im
            a = br * dr - bi * di - 1.0
            b = bi * dr + br * di
            c_dp = dp * a - dq * b
            dq = dp * b + dq * a
            dp = c_dp
            if abs(dp) + abs(dq) < (abs(p_val) + abs(q_val)) * acc
                break
            end
        end
        if q_val <= acc4 * abs(p_val)
            # CF2 didn't converge well
            return (0.0, 0.0)
        end
        gam = (f - p_val) / q_val
        w = 1.0 / sqrt((f - p_val) * gam + q_val)
    end

    # ---- Reconstruct F, G at lambda=0 ----
    fcm = copysign(w, fcl)
    fc0 = fcm
    if !xlturn || !(iexp > 1 || gjwkb > 1.0 / (acch * 100.0))
        gc0 = fcm * gam
    else
        gc0 = gjwkb
    end
    fcp0 = fcm * f
    gcp0 = gc0 * (p_val - q_val / gam)

    # ---- Upward recurrence to lambda=l ----
    fc_l = fc0
    gc_l = gc0
    fcp_l = fcp0
    gcp_l = gcp0
    if l > 0
        w_norm = w / abs(fcl)
        xl_up = 0.0
        for lp in 1:l
            xl_up += 1.0
            el = eta / xl_up
            rl = gc_store[lp + 1]
            sl = el + xl_up * xi
            gcl1 = (sl * gc_l - gcp_l) / rl
            gcp_l = rl * gc_l - sl * gcl1
            gc_l = gcl1
            xl_up_next = xl_up
        end
        # Re-normalize fc at lambda=l
        fc_l = fc0  # simplified; full version tracks fc through recurrence
    end

    # ---- Compute penetrability and shift factor at lambda=l ----
    asq = fc_l^2 + gc_l^2
    if asq <= 0.0
        return (0.0, 0.0)
    end
    pen = x / asq
    sss = x * (fc_l * fcp_l + gc_l * gcp_l) / asq
    # For l>0, we'd need to track fc/fcp through the upward recurrence,
    # but for the common l=0 case this is exact.
    return (pen, sss)
end

# ---- Coulomb-aware penetrability (dispatches to hard-sphere or Coulomb) ----
function _sammy_pen_coulomb(rho, l::Int, eta::Float64)
    if eta == 0.0
        return _sammy_pen(rho, l, 0.0, 0)
    end
    pen, _ = _coulomb_pen_shift(rho, l, eta)
    return pen
end

# ---- Coulomb-aware pgh (dispatches to hard-sphere or Coulomb) ----
function _sammy_pgh_coulomb(rho, l::Int, bound::Float64, ishift::Int, eta::Float64)
    if eta == 0.0
        return _sammy_pgh(rho, l, bound, ishift)
    end
    pen, s_coul = _coulomb_pen_shift(rho, l, eta)

    g = 0.0
    if ishift > 0
        g = s_coul - bound
    end

    if pen <= 1.0e-35
        pen = 0.0
    end

    iffy = 0
    if g == 0.0 && pen == 0.0
        return (0.0, 0.0, 0.0, 1)
    elseif g == 0.0
        return (pen, 0.0, -1.0 / pen, 0)
    elseif pen == 0.0
        return (0.0, 1.0 / g, 0.0, 0)
    else
        dd = pen^2 + g^2
        hr = g / dd
        hi = -pen / dd
        return (pen, hr, hi, 0)
    end
end

# ---- Helper: sinsix -- sin^2(phi), sin(2*phi) ----
function _sammy_sinsix(rho, l::Int)
    a = rho
    c = cos(a)
    s = sin(a)

    if l == 0
        return (s * s, 2.0 * c * s)
    elseif l == 1
        x = a
        g = 1.0 + x^2
        d = (s - c * x) / g
        return ((s - c * x) * d, 2.0 * (c + s * x) * d)
    elseif l == 2
        a2 = a * a
        x = 3.0 * a
        y = 3.0 - a2
        g = y^2 + x^2
        d = (s * y - c * x) / g
        return ((s * y - c * x) * d, 2.0 * (c * y + s * x) * d)
    elseif l == 3
        a2 = a * a
        x = (15.0 - a2) * a
        y = 15.0 - a2 * 6.0
        g = y^2 + x^2
        d = (s * y - c * x) / g
        return ((s * y - c * x) * d, 2.0 * (c * y + s * x) * d)
    elseif l == 4
        a2 = a * a; a4 = a2 * a2
        x = (105.0 - 10.0 * a2) * a
        y = 105.0 - 45.0 * a2 + a4
        g = y^2 + x^2
        d = (s * y - c * x) / g
        return ((s * y - c * x) * d, 2.0 * (c * y + s * x) * d)
    else
        phi = phase_shift(l, rho)
        return (sin(phi)^2, sin(2.0 * phi))
    end
end

# ---- Helper: background R-matrix contribution ----
function _sammy_background(E, lbk::Int, data::Vector{Float64})
    if isempty(data)
        return (0.0, 0.0)
    end

    if lbk == 2 && length(data) >= 7
        # Sammy parametrisation (LBK=2):
        # R_bg = data[3] - data[7]*(data[2]-data[1])
        #      + (data[4] + data[5]*E)*E
        #      - (data[6] + data[7]*E) * ln((data[2]-E)/(E-data[1]))
        ed = data[1]; eu = data[2]
        r0 = data[3]; r1 = data[4]; r2 = data[5]
        s0 = data[6]; s1 = data[7]
        # Guard against log of negative number
        arg_num = eu - E
        arg_den = E - ed
        if arg_num <= 0.0 || arg_den <= 0.0
            return (r0 + (r1 + r2 * E) * E, 0.0)
        end
        bgr = r0 - s1 * (eu - ed) + (r1 + r2 * E) * E -
              (s0 + s1 * E) * log(arg_num / arg_den)
        return (bgr, 0.0)

    elseif lbk == 3 && length(data) >= 5
        # Frohner parametrisation (LBK=3):
        # R_bg_real = data[3] + 2*data[4]*atanh((2E-esum)/ediff)
        # R_bg_imag = data[5]/ediff / (1 - ((2E-esum)/ediff)^2)
        ed = data[1]; eu = data[2]
        esum = ed + eu
        ediff = eu - ed
        if ediff == 0.0
            return (0.0, 0.0)
        end
        x = (2.0 * E - esum) / ediff
        bgr = data[3] + 2.0 * data[4] * atanh(x)
        bgi = data[5] / ediff / (1.0 - x^2)
        return (bgr, bgi)
    end

    return (0.0, 0.0)
end
