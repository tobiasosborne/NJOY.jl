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
        zke = zeros(typeof(float(E)), nchan)
        zkfe = zeros(typeof(float(E)), nchan)
        zkte = zeros(typeof(float(E)), nchan)
        echan = zeros(typeof(float(E)), nchan)

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
            p_val, hr, hi, iffy = _sammy_pgh(rho, lsp, sg.bound[ic], Int(pp.ishift))

            if iffy == 0 && p_val > 1.0e-8
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
