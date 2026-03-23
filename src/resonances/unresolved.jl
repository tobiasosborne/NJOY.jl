# Unresolved resonance cross section evaluation for RECONR
#
# Implements the Fortran genunr/csunr1/sigunr/gnrl chain from reconr.f90.
# Builds a table of infinitely-dilute average cross sections at the
# unresolved energy nodes, then interpolates at runtime via log-log.
#
# Corresponds to:
#   genunr  (1628-1735)  -> build_unresolved_table
#   csunr1  (3826-4077)  -> _csunr1
#   sigunr  (1737-1769)  -> eval_unresolved
#   gnrl    (4498-4644)  -> _gnrl
#   unfac   (4473-4496)  -> _unfac

# ─── Data types ──────────────────────────────────────────────────────

"""One (l, J) spin sequence in the unresolved range."""
struct URRSequence
    l::Int
    J::Float64
    D::Float64       # mean level spacing [eV]
    AMUN::Float64    # degrees of freedom for neutron width
    GN0::Float64     # average reduced neutron width [eV]
    GG::Float64      # average gamma (capture) width [eV]
    GX::Float64      # average competitive width [eV]
    MUF::Int         # degrees of freedom for fission width
    GF_values::Vector{Float64}  # fission widths at energy nodes (length NE)
end

"""Complete unresolved resonance data for one isotope range."""
struct URRData <: AbstractResonanceFormalism
    EL::Float64
    EH::Float64
    SPI::Float64     # target spin
    AP::Float64      # scattering radius [fm → 1e-12 cm in CGS]
    AWRI::Float64    # atomic weight ratio
    LSSF::Int        # self-shielding flag
    LFW::Int         # fission width flag
    NAPS::Int        # channel radius option
    NRO::Int         # energy-dependent scattering radius flag
    energies::Vector{Float64}   # energy nodes for fission widths (NE+1 entries)
    sequences::Vector{URRSequence}
end

"""One (l, J) spin sequence for mode=12 (LRF=2): ALL widths energy-dependent."""
struct URR2Sequence
    l::Int
    J::Float64
    INT::Int           # interpolation law for parameters (2=linear)
    AMUX::Float64      # DOF for competitive width (lambda)
    AMUN::Float64      # DOF for neutron width (mu)
    AMUF::Float64      # DOF for fission width (nu)
    energies::Vector{Float64}
    D::Vector{Float64}         # level spacing at each node
    GX::Vector{Float64}        # competitive width
    GNO::Vector{Float64}       # reduced neutron width
    GG::Vector{Float64}        # gamma width
    GF::Vector{Float64}        # fission width
end

"""Unresolved data for LRU=2/LRF=2 (mode=12, all widths energy-dependent)."""
struct URR2Data <: AbstractResonanceFormalism
    EL::Float64
    EH::Float64
    SPI::Float64
    AP::Float64
    AWRI::Float64
    LSSF::Int
    NAPS::Int
    NRO::Int
    sequences::Vector{URR2Sequence}
end

"""Pre-computed unresolved XS table (from genunr)."""
struct URRTable
    energies::Vector{Float64}    # energy nodes
    total::Vector{Float64}       # σ_total at each node
    elastic::Vector{Float64}     # σ_elastic at each node
    fission::Vector{Float64}     # σ_fission at each node
    capture::Vector{Float64}     # σ_capture at each node
    intlaw::Int                  # interpolation law (5=log-log for LRF=1)
end

# ─── Fluctuation integrals (gnrl) ────────────────────────────────────

"""
Fluctuation integral for unresolved resonances using Hwang 10-point
Gauss-Laguerre quadrature. Matches Fortran `gnrl` (reconr.f90:4498-4644).

- `galpha`: neutron width
- `gbeta`:  fission width (0 if no fission)
- `gamma`:  capture width
- `mu,nu,lamda`: DOF indices (1-4) for neutron, fission, competitive
- `diff`:   competitive width (0 if none)
- `id`:     1=elastic, 2=capture, 3=fission
"""
function _gnrl(galpha, gbeta, gamma, mu, nu, lamda, diff, id)
    galpha <= 0 && return 0.0
    gamma <= 0 && return 0.0
    gbeta < 0 && return 0.0
    gbeta > 0 && diff < 0 && return 0.0

    qw = HWANG_QW
    qp = HWANG_QP
    s = 0.0

    if gbeta == 0 && diff == 0
        for j in 1:10
            xj = qp[j, mu]
            denom = galpha * xj + gamma
            s += qw[j, mu] * (id == 1 ? xj * xj : xj) / denom
        end
    elseif gbeta == 0 && diff > 0
        for j in 1:10
            xj = qp[j, mu]
            for k in 1:10
                denom = galpha * xj + gamma + diff * qp[k, lamda]
                s += qw[j, mu] * qw[k, lamda] * (id == 1 ? xj * xj : xj) / denom
            end
        end
    elseif gbeta > 0 && diff == 0
        if id == 1
            for j in 1:10
                xj = qp[j, mu]
                for k in 1:10
                    s += qw[j, mu] * qw[k, nu] * xj * xj /
                         (galpha * xj + gbeta * qp[k, nu] + gamma)
                end
            end
        elseif id == 2
            for j in 1:10
                xj = qp[j, mu]
                for k in 1:10
                    s += qw[j, mu] * qw[k, nu] * xj /
                         (galpha * xj + gbeta * qp[k, nu] + gamma)
                end
            end
        elseif id == 3
            for j in 1:10
                for k in 1:10
                    s += qw[j, mu] * qw[k, nu] * qp[j, mu] * qp[k, nu] /
                         (galpha * qp[j, mu] + gbeta * qp[k, nu] + gamma)
                end
            end
        end
    else  # gbeta > 0, diff > 0
        if id == 1
            for j in 1:10
                xj = qp[j, mu]
                for k in 1:10
                    for l in 1:10
                        s += qw[j, mu] * qw[k, nu] * qw[l, lamda] * xj * xj /
                             (galpha * xj + gbeta * qp[k, nu] + gamma + diff * qp[l, lamda])
                    end
                end
            end
        elseif id == 2
            for j in 1:10
                for k in 1:10
                    for l in 1:10
                        s += qw[j, mu] * qw[k, nu] * qw[l, lamda] * qp[j, mu] /
                             (galpha * qp[j, mu] + gbeta * qp[k, nu] + gamma + diff * qp[l, lamda])
                    end
                end
            end
        elseif id == 3
            for j in 1:10
                for k in 1:10
                    for l in 1:10
                        s += qw[j, mu] * qw[k, nu] * qw[l, lamda] * qp[j, mu] * qp[k, nu] /
                             (galpha * qp[j, mu] + gbeta * qp[k, nu] + gamma + diff * qp[l, lamda])
                    end
                end
            end
        end
    end
    s
end

# ─── Penetrability with AMUN (unfac) ────────────────────────────────

"""Penetrability factor including AMUN, and phase shift. Matches Fortran `unfac`."""
function _unfac(l::Int, rho, rhoc, amun)
    r2 = rho * rho
    if l == 0
        return amun, rhoc
    elseif l == 1
        return amun * r2 / (1 + r2), rhoc - atan(rhoc)
    else
        r4 = r2 * r2
        return amun * r4 / (9 + 3r2 + r4), rhoc - atan(3rhoc / (3 - r2))
    end
end

# ─── Cross section evaluation (csunr1) ──────────────────────────────

const _RC1 = 0.123
const _RC2 = 0.08
const _THIRD_URR = 0.333333333  # Fortran truncated 1/3

"""
Evaluate unresolved average cross sections at energy `E` using SLBW
formalism (mode=11). Matches Fortran `csunr1` (reconr.f90:3826-4077).

Returns (total, elastic, fission, capture).
"""
function _csunr1(E::Float64, urr::URRData)
    C = PhysicsConstants
    cwaven = sqrt(2 * C.amassn * C.amu * C.ev) * 1e-12 / C.hbar

    awri = urr.AWRI
    rat = awri / (awri + 1)
    aw = awri * C.amassn
    aa = urr.NAPS == 0 ? _RC1 * aw^_THIRD_URR + _RC2 : urr.AP
    ay = urr.AP
    spi = urr.SPI
    cnst = (2 * C.pi^2) / (cwaven * rat)^2

    sig_e = 0.0; sig_f = 0.0; sig_c = 0.0; spot = 0.0
    prev_l = -1

    for seq in urr.sequences
        ll = seq.l
        gj = (2 * seq.J + 1) / (4 * spi + 2)
        mu = Int(seq.AMUN)
        nu = seq.MUF
        lamda = 0

        # Fission width: interpolate from energy-dependent table
        gfx = if !isempty(seq.GF_values) && length(urr.energies) >= 2
            # Find bracketing interval in fission width energy grid
            ne_fw = min(length(urr.energies), length(seq.GF_values))
            i1 = clamp(searchsortedlast(urr.energies, E), 1, ne_fw - 1)
            e1, e2 = urr.energies[i1], urr.energies[i1 + 1]
            g1, g2 = seq.GF_values[i1], seq.GF_values[i1 + 1]
            # Linear interpolation (intlaw=2 for fission widths)
            e1 == e2 ? g1 : g1 + (g2 - g1) * (E - e1) / (e2 - e1)
        else
            0.0
        end

        e2v = sqrt(E)
        k = rat * e2v * cwaven
        rho = k * aa
        rhoc = k * ay
        vl, ps = _unfac(ll, rho, rhoc, seq.AMUN)
        vl *= e2v

        # Potential scattering (once per l)
        if ll != prev_l
            spot_l = 4 * C.pi * (2ll + 1) * (sin(ps) / k)^2
            spot += spot_l
            prev_l = ll
        end

        gnx = seq.GN0 * vl
        diff = seq.GX
        den = E * seq.D
        temp = cnst * gj * gnx / den
        terg = temp * seq.GG
        ters = temp * gnx
        terf = temp * gfx

        gs  = _gnrl(gnx, gfx, seq.GG, mu, nu, lamda, diff, 1)
        gc  = _gnrl(gnx, gfx, seq.GG, mu, nu, lamda, diff, 2)
        gff = _gnrl(gnx, gfx, seq.GG, mu, nu, lamda, diff, 3)

        gc  *= terg
        gff *= terf
        gs  *= ters

        # Interference correction
        add = cnst * gj * 2 * gnx * sin(ps)^2 / den
        gs -= add

        sig_e += gs
        sig_f += gff
        sig_c += gc
    end

    sig_e += spot
    total = sig_e + sig_f + sig_c
    (total, sig_e, sig_f, sig_c)
end

# ─── Cross section evaluation (csunr2, mode=12) ─────────────────────

"""
Evaluate unresolved average cross sections at energy `E` for LRF=2
(mode=12, all widths energy-dependent). Matches Fortran `csunr2`
(reconr.f90:4079-4317).

Returns (total, elastic, fission, capture).
"""
function _csunr2(E::Float64, urr::URR2Data)
    C = PhysicsConstants
    cwaven = sqrt(2 * C.amassn * C.amu * C.ev) * 1e-12 / C.hbar

    awri = urr.AWRI
    rat = awri / (awri + 1)
    aw = awri * C.amassn
    aa = urr.NAPS == 0 ? _RC1 * aw^_THIRD_URR + _RC2 : urr.AP
    ay = urr.AP
    spi = urr.SPI
    cnst = (2 * C.pi^2) / (cwaven * rat)^2

    e2 = sqrt(E)
    k = cwaven * rat * e2
    rho = k * aa
    rhoc = k * ay

    sig_e = 0.0; sig_f = 0.0; sig_c = 0.0; spot = 0.0
    prev_l = -1; j_in_l = 0

    for seq in urr.sequences
        ll = seq.l
        if ll != prev_l
            j_in_l = 0
            prev_l = ll
        end
        j_in_l += 1

        ne = length(seq.energies)
        ne < 2 && continue

        # Interpolate ALL parameters at energy E (Fortran lines 4231-4242)
        i1 = 1
        while i1 + 1 < ne && seq.energies[i1 + 1] <= E - 1e-8
            i1 += 1
        end
        i2 = i1 + 1
        e1, e2v = seq.energies[i1], seq.energies[i2]
        dx   = e1 == e2v ? seq.D[i1]   : seq.D[i1]   + (seq.D[i2]   - seq.D[i1])   * (E - e1) / (e2v - e1)
        gxx  = e1 == e2v ? seq.GX[i1]  : seq.GX[i1]  + (seq.GX[i2]  - seq.GX[i1])  * (E - e1) / (e2v - e1)
        gnox = e1 == e2v ? seq.GNO[i1] : seq.GNO[i1] + (seq.GNO[i2] - seq.GNO[i1]) * (E - e1) / (e2v - e1)
        ggx  = e1 == e2v ? seq.GG[i1]  : seq.GG[i1]  + (seq.GG[i2]  - seq.GG[i1])  * (E - e1) / (e2v - e1)
        gfx  = e1 == e2v ? seq.GF[i1]  : seq.GF[i1]  + (seq.GF[i2]  - seq.GF[i1])  * (E - e1) / (e2v - e1)

        # Clamp small values (Fortran lines 4243-4244)
        gxx < 1e-8 && (gxx = 0.0)
        gfx < 1e-8 && (gfx = 0.0)

        gj = (2 * seq.J + 1) / (4 * spi + 2)
        mu = Int(seq.AMUN)
        nu = Int(seq.AMUF)
        lamda = Int(seq.AMUX)

        vl, ps = _unfac(ll, rho, rhoc, seq.AMUN)
        vl *= e2

        # Potential scattering (once per l, j<=1 — Fortran line 4249)
        if j_in_l <= 1
            spot_l = 4 * C.pi * (2ll + 1) * (sin(ps) / k)^2
            spot += spot_l
        end

        gnx = gnox * vl
        diff = gxx
        den = E * dx
        temp = cnst * gj * gnx / den
        terg = temp * ggx
        ters = temp * gnx
        terf = temp * gfx

        gs  = _gnrl(gnx, gfx, ggx, mu, nu, lamda, diff, 1)
        gc  = _gnrl(gnx, gfx, ggx, mu, nu, lamda, diff, 2)
        gff = _gnrl(gnx, gfx, ggx, mu, nu, lamda, diff, 3)

        gc  *= terg
        gff *= terf
        gs  *= ters

        # Interference correction — two-step to match Fortran evaluation
        # order (reconr.f90:4271-4272: add=const*gj*2*gnx*(sin(ps))**2;
        # add=add/(e*dx))
        add = cnst * gj * 2 * gnx * sin(ps)^2
        add /= den
        gs -= add

        sig_e += gs
        sig_f += gff
        sig_c += gc
    end

    sig_e += spot
    total = sig_e + sig_f + sig_c
    (total, sig_e, sig_f, sig_c)
end

# ─── Table builder (genunr) ──────────────────────────────────────────

# Standard energy grid from Fortran rdf2u1 (reconr.f90:1326-1338)
const _EGRIDU = [
    1e1, 1.25e1, 1.5e1, 1.7e1, 2e1, 2.5e1, 3e1, 3.5e1, 4e1, 5e1,
    6e1, 7.2e1, 8.5e1, 1e2, 1.25e2, 1.5e2, 1.7e2, 2e2, 2.5e2, 3e2,
    3.5e2, 4e2, 5e2, 6e2, 7.2e2, 8.5e2, 1e3, 1.25e3, 1.5e3, 1.7e3,
    2e3, 2.5e3, 3e3, 3.5e3, 4e3, 5e3, 6e3, 7.2e3, 8.5e3, 1e4,
    1.25e4, 1.5e4, 1.7e4, 2e4, 2.5e4, 3e4, 3.5e4, 4e4, 5e4, 6e4,
    7.2e4, 8.5e4, 1e5, 1.25e5, 1.5e5, 1.7e5, 2e5, 2.5e5, 3e5, 3.5e5,
    4e5, 5e5, 6e5, 7.2e5, 8.5e5, 1e6, 1.25e6, 1.5e6, 1.7e6, 2e6,
    2.5e6, 3e6, 3.5e6, 4e6, 5e6, 6e6, 7.2e6, 8.5e6]

"""
Build the unresolved XS table by evaluating at energy nodes.
Matches Fortran `genunr` (reconr.f90:1628-1735).

Adds intermediate points from egridu between energy nodes when the
step ratio exceeds 1.26 (matching rdf2u1 lines 1383-1425).
"""
function build_unresolved_table(urr::URRData, abn::Float64,
                                 mf3_sections::Vector{MF3Section};
                                 eresr::Float64 = 0.0)
    # Build dense energy grid matching Fortran rdf2u1 + rdfil2 boundary nodes.
    # rdfil2 (lines 754-776) adds sigfig(EL,7,±1) and sigfig(EH,7,±1) to eunr.
    # Entries below eresr are marked negative (resolved-unresolved overlap,
    # rdfil2 lines 864: bg suppressed in genunr for negative-signed energies).
    energies = Float64[]

    # rdfil2 boundary nodes for EL (if EL is not at elow=1e-5)
    el_down = round_sigfig(urr.EL, 7, -1)
    el_up   = round_sigfig(urr.EL, 7, +1)
    push!(energies, el_down < eresr ? -el_down : el_down)
    push!(energies, el_up   < eresr ? -el_up   : el_up)

    # rdf2u1 fission-width energy nodes + egridu intermediates
    wide = 1.26
    raw_e = urr.energies
    for n in 1:length(raw_e)
        ener = raw_e[n]
        ener < urr.EL && continue
        ener >= urr.EH && continue
        push!(energies, round_sigfig(ener, 7, 0))
        if n < length(raw_e)
            enex = raw_e[n + 1]
            if enex > wide * ener
                enut = ener
                while enut < enex
                    idx = findfirst(e -> e > enut + enut/100, _EGRIDU)
                    idx === nothing && break
                    enut = _EGRIDU[idx]
                    enut >= enex && break
                    push!(energies, enut)
                end
            end
        end
    end

    # rdfil2 boundary nodes for EH
    push!(energies, round_sigfig(urr.EH, 7, -1))
    push!(energies, round_sigfig(urr.EH, 7, +1))

    # Sort by absolute value (Fortran rdfil2 line 856: call order(eunr,nunr))
    sort!(energies, by=abs)

    # Deduplicate matching Fortran rdfil2 lines 857-869: skip entries
    # too close to the previous one (within sigfig(abs(e), 7, +1)).
    eps_tol = 1.0e-10
    eresu = urr.EL
    deduped = Float64[]
    prev_abs = 0.0
    for e in energies
        ae = abs(e)
        ae <= (1 - eps_tol) * eresu && continue
        ae >= (1 + eps_tol) * urr.EH && continue
        ae < prev_abs && continue
        push!(deduped, e)
        prev_abs = round_sigfig(ae, 7, +1)
    end
    energies = deduped
    ne = length(energies)
    total = zeros(ne)
    elastic = zeros(ne)
    fission = zeros(ne)
    capture = zeros(ne)

    for i in 1:ne
        E = abs(energies[i])
        if E >= urr.EL && E < urr.EH
            sig = _csunr1(E, urr)
            total[i]   = round_sigfig(abn * sig[1], 7)
            elastic[i] = round_sigfig(abn * sig[2], 7)
            fission[i] = round_sigfig(abn * sig[3], 7)
            capture[i] = round_sigfig(abn * sig[4], 7)
        end
    end

    # Add MF3 backgrounds only when LSSF=0 (Fortran genunr lines 1695-1731).
    # LSSF=1 means backgrounds are NOT included in the table; emerge adds them.
    if urr.LSSF == 0
        for sec in mf3_sections
            mt = Int(sec.mt)
            (mt != 1 && mt != 2 && mt != 18 && mt != 102) && continue
            for i in 1:ne
                E = energies[i]
                E < 0 && continue
                bkg = interpolate(sec.tab, E)
                if     mt == 1;   total[i]   = round_sigfig(total[i]   + bkg, 7)
                elseif mt == 2;   elastic[i] = round_sigfig(elastic[i] + bkg, 7)
                elseif mt == 18;  fission[i] = round_sigfig(fission[i] + bkg, 7)
                elseif mt == 102; capture[i] = round_sigfig(capture[i] + bkg, 7)
                end
            end
        end
    end

    # Update total after adding backgrounds (matching Fortran)
    for i in 1:ne
        total[i] = round_sigfig(elastic[i] + fission[i] + capture[i], 7)
    end

    # intlaw = 5 (log-log) for LRF=1 (mode=11)
    URRTable(energies, total, elastic, fission, capture, 5)
end

"""
Build unresolved XS table for mode=12 (LRF=2, all widths energy-dependent).
Same structure as mode=11 but uses `_csunr2` and different grid generation.
"""
function build_unresolved_table(urr::URR2Data, abn::Float64,
                                 mf3_sections::Vector{MF3Section};
                                 eresr::Float64 = 0.0)
    energies = Float64[]

    # rdfil2 boundary nodes for EL
    el_down = round_sigfig(urr.EL, 7, -1)
    el_up   = round_sigfig(urr.EL, 7, +1)
    push!(energies, el_down < eresr ? -el_down : el_down)
    push!(energies, el_up   < eresr ? -el_up   : el_up)

    # rdf2u2 energy nodes + egridu intermediates
    # Uses first sequence's energy grid (Fortran: n==1, l==1)
    wide = 1.26
    raw_e = isempty(urr.sequences) ? Float64[] : urr.sequences[1].energies
    for n in 1:length(raw_e)
        ener = raw_e[n]
        ener < urr.EL && continue
        ener >= urr.EH && continue
        push!(energies, round_sigfig(ener, 7, 0))
        if n < length(raw_e)
            enex = raw_e[n + 1]
            if enex > wide * ener
                # Gap filling: rdf2u2 uses ener/1000 threshold (vs rdf2u1's enut/100)
                enut = ener
                while true
                    idx = findfirst(e -> e > enut + enut / 1000, _EGRIDU)
                    idx === nothing && break
                    enut = _EGRIDU[idx]
                    enut >= enex && break
                    push!(energies, enut)
                end
            end
        end
    end

    # rdfil2 boundary nodes for EH
    push!(energies, round_sigfig(urr.EH, 7, -1))
    push!(energies, round_sigfig(urr.EH, 7, +1))

    sort!(energies, by=abs)

    # Deduplicate (matching Fortran rdfil2 lines 857-869)
    eps_tol = 1.0e-10
    deduped = Float64[]
    prev_abs = 0.0
    for e in energies
        ae = abs(e)
        ae <= (1 - eps_tol) * urr.EL && continue
        ae >= (1 + eps_tol) * urr.EH && continue
        ae < prev_abs && continue
        push!(deduped, e)
        prev_abs = round_sigfig(ae, 7, +1)
    end
    energies = deduped

    ne = length(energies)
    total = zeros(ne); elastic = zeros(ne); fission = zeros(ne); capture = zeros(ne)

    # Evaluate at each node using csunr2
    for i in 1:ne
        E = abs(energies[i])
        if E >= urr.EL && E < urr.EH
            sig = _csunr2(E, urr)
            total[i]   = round_sigfig(abn * sig[1], 7)
            elastic[i] = round_sigfig(abn * sig[2], 7)
            fission[i] = round_sigfig(abn * sig[3], 7)
            capture[i] = round_sigfig(abn * sig[4], 7)
        end
    end

    # Add MF3 backgrounds only when LSSF=0 (Fortran genunr lines 1695-1731)
    if urr.LSSF == 0
        for sec in mf3_sections
            mt = Int(sec.mt)
            (mt != 1 && mt != 2 && mt != 18 && mt != 102) && continue
            for i in 1:ne
                E = energies[i]
                E < 0 && continue
                bkg = interpolate(sec.tab, E)
                if     mt == 1;   total[i]   = round_sigfig(total[i]   + bkg, 7)
                elseif mt == 2;   elastic[i] = round_sigfig(elastic[i] + bkg, 7)
                elseif mt == 18;  fission[i] = round_sigfig(fission[i] + bkg, 7)
                elseif mt == 102; capture[i] = round_sigfig(capture[i] + bkg, 7)
                end
            end
        end
    end

    for i in 1:ne
        total[i] = round_sigfig(elastic[i] + fission[i] + capture[i], 7)
    end

    # Determine interpolation law from first sequence's INT
    intlaw = isempty(urr.sequences) ? 2 : urr.sequences[1].INT
    URRTable(energies, total, elastic, fission, capture, intlaw)
end

# ─── Table interpolation (sigunr) ───────────────────────────────────

"""
Interpolate unresolved XS from pre-computed table.
Matches Fortran `sigunr` (reconr.f90:1737-1769).
Uses log-log interpolation (intlaw=5).
"""
function eval_unresolved(table::URRTable, E::Float64)
    n = length(table.energies)
    n < 2 && return (0.0, 0.0, 0.0, 0.0)

    # Find bracketing interval
    i2 = 2
    for i in 2:n
        i2 = i
        abs(table.energies[i]) > E && break
    end
    i1 = i2 - 1

    e1 = abs(table.energies[i1])
    e2 = abs(table.energies[i2])
    E < e1 && return (0.0, 0.0, 0.0, 0.0)
    E > e2 && return (0.0, 0.0, 0.0, 0.0)

    e1 == e2 && return (table.total[i1], table.elastic[i1], table.fission[i1], table.capture[i1])

    if table.intlaw == 2
        # Linear interpolation (mode=12, LRF=2)
        frac = (E - e1) / (e2 - e1)
        _lin(a, b) = a + (b - a) * frac
        return (_lin(table.total[i1], table.total[i2]),
                _lin(table.elastic[i1], table.elastic[i2]),
                _lin(table.fission[i1], table.fission[i2]),
                _lin(table.capture[i1], table.capture[i2]))
    end

    # Log-log interpolation (terp1 with intlaw=5, mode=11)
    if e1 > 0 && e2 > 0
        logr = log(E / e1) / log(e2 / e1)
        t = _loglog_interp(table.total[i1],   table.total[i2],   logr)
        el = _loglog_interp(table.elastic[i1], table.elastic[i2], logr)
        f = _loglog_interp(table.fission[i1], table.fission[i2], logr)
        c = _loglog_interp(table.capture[i1], table.capture[i2], logr)
        return (t, el, f, c)
    end
    (0.0, 0.0, 0.0, 0.0)
end

function _loglog_interp(y1, y2, logr)
    (y1 <= 0 || y2 <= 0) && return y1 + (y2 - y1) * logr  # fallback to lin-log
    exp(log(y1) + logr * log(y2 / y1))
end
