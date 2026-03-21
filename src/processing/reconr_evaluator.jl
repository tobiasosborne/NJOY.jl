# RECONR evaluator -- build cross section evaluation closures from MF2 data
#
# Step 1 of the RECONR pipeline: build_evaluator returns a closure that
# evaluates resonance cross sections at a given energy.
#
# Also contains the legacy sigma_mf2 function and merge_background functions.

# ==========================================================================
# Step 1: Build evaluator (closure over MF2 data)
# ==========================================================================

"""
    build_evaluator(mf2::MF2Data; temperature=0.0, table=nothing)

Build a cross section evaluator closure from MF2 resonance data.
Returns a callable `f(E::Float64) -> NTuple{4, Float64}` that computes
`(total, elastic, fission, capture)` at energy `E`.

The returned function:
- Sums over all isotope sections weighted by abundance
- Dispatches to the correct formalism via `cross_section`
- Clamps negative components to zero
- Is fully type-stable (returns `NTuple{4, Float64}`)

This replaces NJOY2016's `sigma` subroutine (reconr.f90:2571-2667).
"""
function build_evaluator(mf2::MF2Data;
                         temperature::Real = 0.0,
                         table::Union{Nothing, FaddeevaTable} = nothing)
    # Pre-collect the isotope/range/abundance triples to avoid
    # iterating nested structures in the hot path
    sections = _collect_sections(mf2)

    function evaluator(E::Float64)
        sig_t = 0.0
        sig_e = 0.0
        sig_f = 0.0
        sig_c = 0.0

        @inbounds for (rng, abn) in sections
            if E >= rng.EL && E < rng.EH && Int(rng.LRU) > 0
                sigp = cross_section(E, rng; temperature=temperature, table=table)
                sig_t += max(0.0, sigp.total) * abn
                sig_e += max(0.0, sigp.elastic) * abn
                sig_f += max(0.0, sigp.fission) * abn
                sig_c += max(0.0, sigp.capture) * abn
            end
        end
        return (sig_t, sig_e, sig_f, sig_c)
    end

    return evaluator
end

# Flatten isotope/range hierarchy into a flat vector of (range, abundance) pairs
function _collect_sections(mf2::MF2Data)
    result = Tuple{ResonanceRange, Float64}[]
    for iso in mf2.isotopes
        abn = iso.ABN
        for rng in iso.ranges
            push!(result, (rng, abn))
        end
    end
    return result
end

# ==========================================================================
# Merge MF3 background into pointwise result
# ==========================================================================

"""
    merge_background!(energies, values, mf3_sections, mf2)

Merge MF3 background cross sections into the resonance cross section
matrix `values` in-place. This replaces NJOY2016's `emerge`.

For each energy and matching MF3 reaction:
- MT 2 -> add to elastic column
- MT 18, 19 -> add to fission column
- MT 102 -> add to capture column
Then recompute total = elastic + fission + capture.
"""
function merge_background!(energies::Vector{Float64},
                            values::Matrix{Float64},
                            mf3_sections::Vector{MF3Section},
                            mf2::MF2Data)
    n = length(energies)
    small = 1.0e-8

    # Find resonance range boundaries for overlap detection
    eresl, eresh, eresr = _resonance_bounds(mf2)

    # Redundant MTs that should never be merged (they are computed as sums)
    _SKIP_MTS = (1, 3, 4, 101, 27)

    # Collect non-primary background contributions that only affect total.
    # We accumulate them separately and add once after the channel columns
    # are finalised, so they participate in the total but do not pollute
    # the elastic/fission/capture columns.
    for i in 1:n
        e = energies[i]
        other_bg = 0.0   # background from non-primary MTs

        for sec in mf3_sections
            mt = Int(sec.mt)
            # Skip redundant reactions
            mt in _SKIP_MTS && continue

            bg = interpolate(sec.tab, e)
            bg == 0.0 && continue

            # In the unresolved-resolved overlap region, backgrounds for
            # primary channels are assigned to the unresolved component
            # (matching NJOY's emerge logic)
            if e >= eresr && e < eresh && (mt == 2 || mt == 18 || mt == 19 ||
                                           mt == 20 || mt == 21 || mt == 38 ||
                                           mt == 102)
                continue
            end

            if mt == 2
                values[i, _COL_ELASTIC] += bg
            elseif mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38
                values[i, _COL_FISSION] += bg
            elseif mt == 102
                values[i, _COL_CAPTURE] += bg
            else
                # All other non-redundant MTs contribute to total only
                other_bg += bg
            end
        end

        # Clamp elastic to a small positive value
        if values[i, _COL_ELASTIC] <= small
            values[i, _COL_ELASTIC] = small
        end

        # Round to 7 significant figures (matching NJOY's emerge:sigfig call)
        for j in _COL_ELASTIC:_COL_CAPTURE
            values[i, j] = round_sigfig(values[i, j], 7, 0)
        end

        # Recompute total = elastic + fission + capture + other backgrounds
        values[i, _COL_TOTAL] = values[i, _COL_ELASTIC] +
                                 values[i, _COL_FISSION] +
                                 values[i, _COL_CAPTURE] +
                                 other_bg
    end
end

# Extract resonance range boundaries from MF2 data
function _resonance_bounds(mf2::MF2Data)
    eresl = Inf
    eresh = 0.0
    eresr = 0.0
    for iso in mf2.isotopes
        for rng in iso.ranges
            Int(rng.LRU) == 0 && continue
            eresl = min(eresl, rng.EL)
            eresh = max(eresh, rng.EH)
            if Int(rng.LRU) == 1
                eresr = max(eresr, rng.EH)
            end
        end
    end
    eresr = clamp(eresr, eresl, eresh)
    return (eresl, eresh, eresr)
end

# ==========================================================================
# Legacy interface functions
# ==========================================================================

"""
    sigma_mf2(E, mf2) -> CrossSections

Evaluate resonance cross sections at energy E by summing over all isotopes.
Legacy interface returning CrossSections struct.
"""
function sigma_mf2(E::Real, mf2::MF2Data)
    sig = CrossSections()
    E_f = Float64(E)
    for iso in mf2.isotopes
        abn = iso.ABN
        for rng in iso.ranges
            if E_f >= rng.EL && E_f < rng.EH && Int(rng.LRU) > 0
                try
                    sigp = cross_section(E_f, rng)
                    total = max(0.0, sigp.total)
                    elastic = max(0.0, sigp.elastic)
                    fission = max(0.0, sigp.fission)
                    capture = max(0.0, sigp.capture)
                    sig = sig + abn * CrossSections(total, elastic, fission, capture)
                catch e
                    @warn "sigma_mf2: cross_section failed at E=$E_f for range [$(rng.EL), $(rng.EH)]" exception=(e, catch_backtrace())
                end
            end
        end
    end
    return sig
end

"""
    merge_background_legacy(energies, res_xs, mf3_sections) -> Vector{CrossSections}

Legacy merge returning Vector{CrossSections}.
"""
# ==========================================================================
# RML (R-Matrix Limited, LRF=7) evaluator for unsupported formalisms
# ==========================================================================

"""
    RMLChannel

One channel within an RML spin group.
"""
struct RMLChannel
    ipp::Int         # particle-pair index (1-based)
    l::Int           # orbital angular momentum
    sch::Float64     # channel spin
    bnd::Float64     # boundary condition
    ape::Float64     # effective channel radius [fm]
    apt::Float64     # true channel radius [fm]
end

"""
    RMLSammyBG

Sammy background R-matrix parameters (LBK=2).

The Sammy background R-matrix formula (NJOY setr, samm.f90):
  R_bg = p1 - p5*(EU-ED) + (p2 + p3*E)*E - (p4 + p5*E) * log((EU-E)/(E-ED))
"""
struct RMLSammyBG
    LCH::Int      # channel number this background applies to (1-based)
    LBK::Int      # background type (2=Sammy, 3=Frohner)
    ED::Float64   # lower energy boundary (from LIST C1)
    EU::Float64   # upper energy boundary (from LIST C2)
    p1::Float64   # parameter 1 (constant term)
    p2::Float64   # parameter 2 (linear coefficient)
    p3::Float64   # parameter 3 (quadratic coefficient)
    p4::Float64   # parameter 4 (log coefficient)
    p5::Float64   # parameter 5 (E*log coefficient)
end

"""
    RMLSpinGroup

One spin group of the R-Matrix Limited formalism.
"""
struct RMLSpinGroup
    AJ::Float64              # total angular momentum
    PJ::Float64              # parity
    KBK::Int                 # background R-matrix flag
    KPS::Int                 # phase-shift option
    channels::Vector{RMLChannel}
    Er::Vector{Float64}      # resonance energies
    GAM::Matrix{Float64}     # width data (nres x (1+nch_nongamma)):
                             #   col 1 = gamgam (capture width)
                             #   col 2.. = channel widths (non-gamma channels)
    sammy_bg::Vector{RMLSammyBG}  # background R-matrix data (from KBK records)
end

"""
    RMLParticlePair

Particle-pair definition for R-Matrix Limited.
"""
struct RMLParticlePair
    MA::Float64; MB::Float64
    ZA::Float64; ZB::Float64
    IA::Float64; IB::Float64
    Q::Float64
    PNT::Float64; SHF::Float64
    MT::Int
    PA::Float64; PB::Float64
end

"""
    RMLData

Complete R-Matrix Limited (LRF=7) data for one resonance range.
"""
struct RMLData
    EL::Float64
    EH::Float64
    IFG::Int                              # reduced-width amplitude flag
    KRL::Int                              # relativistic kinematics flag
    KRM::Int                              # R-matrix formalism option (0=RM)
    AWR::Float64                          # atomic weight ratio
    pairs::Vector{RMLParticlePair}
    spin_groups::Vector{RMLSpinGroup}
end

"""
    read_rml_data(io::IO) -> Union{Nothing, RMLData}

Read R-Matrix Limited (LRF=7) data from an ENDF file positioned at
the start of MF2/MT151. Returns `nothing` if no LRF=7 range is found.
"""
function read_rml_data(io::IO)
    saved_pos = position(io)
    seekstart(io)

    if !find_section(io, 2, 151)
        seek(io, saved_pos)
        return nothing
    end

    head = read_head(io)
    AWR = head.C2
    NIS = Int(head.N1)

    for _ in 1:NIS
        iso_cont = read_cont(io)
        NER = Int(iso_cont.N1)

        for _ in 1:NER
            rng_cont = read_cont(io)
            EL = rng_cont.C1
            EH = rng_cont.C2
            LRU = Int(rng_cont.L1)
            LRF = Int(rng_cont.L2)

            if LRU == 1 && LRF == 7
                # Read RML data
                rml_cont = read_cont(io)
                IFG = Int(rml_cont.L1)
                KRM = Int(rml_cont.L2)
                NJS = Int(rml_cont.N1)
                KRL = Int(rml_cont.N2)

                # Particle pairs
                pp_list = read_list(io)
                NPP = Int(pp_list.L1)
                pairs = RMLParticlePair[]
                for ipp in 1:NPP
                    base = (ipp - 1) * 12
                    push!(pairs, RMLParticlePair(
                        pp_list.data[base+1], pp_list.data[base+2],
                        pp_list.data[base+3], pp_list.data[base+4],
                        pp_list.data[base+5], pp_list.data[base+6],
                        pp_list.data[base+7], pp_list.data[base+8],
                        pp_list.data[base+9], Int(pp_list.data[base+10]),
                        pp_list.data[base+11], pp_list.data[base+12]))
                end

                # Spin groups
                spin_groups = RMLSpinGroup[]
                for _ in 1:NJS
                    sg_list = read_list(io)
                    AJ = sg_list.C1
                    PJ = sg_list.C2
                    KBK = Int(sg_list.L1)
                    KPS = Int(sg_list.L2)
                    NCH = Int(sg_list.N2)

                    channels = RMLChannel[]
                    for ich in 1:NCH
                        base = (ich - 1) * 6
                        push!(channels, RMLChannel(
                            Int(sg_list.data[base+1]),
                            Int(sg_list.data[base+2]),
                            sg_list.data[base+3],
                            sg_list.data[base+4],
                            sg_list.data[base+5],
                            sg_list.data[base+6]))
                    end

                    # Resonance data
                    # NCH includes the gamma channel -- actual reaction channels
                    # are NCH-1 (the gamma channel has IPP=1 and is excluded from
                    # the channel list by the Fortran code, but the ENDF format
                    # counts it in the channel definitions).
                    # However, the channel list already read above may or may not
                    # include the gamma channel depending on the evaluation.
                    # We keep all channels as-is and let cross_section_rml sort
                    # them by particle-pair MT.
                    res_list = read_list(io)
                    nres = Int(res_list.L2)  # number of resonances from L2

                    # Stride: ENDF pads to 6 or 12 values per resonance
                    # Following Fortran: nx=6; if (2+ichan > 6) nx=12
                    # where ichan = number of non-gamma channels
                    # But we can just derive from data length:
                    nx = nres > 0 ? div(length(res_list.data), nres) : 6

                    # Actual number of width values per resonance:
                    # First value is Er, second is gamgam (capture width),
                    # remaining are channel widths (1 per non-gamma channel)
                    n_non_gamma = count(ch -> ch.ipp >= 1 && ch.ipp <= length(pairs) &&
                                        pairs[ch.ipp].MT != 102, channels)

                    Er_vec = Float64[]
                    GAM_mat = zeros(0, max(NCH, 0))
                    if nres > 0 && NCH > 0
                        Er_vec = zeros(nres)
                        # Store gamgam + all channel widths
                        # Column 1 = gamgam, columns 2.. = channel widths
                        GAM_mat = zeros(nres, 1 + n_non_gamma)
                        for ir in 1:nres
                            base = (ir - 1) * nx
                            Er_vec[ir] = res_list.data[base+1]
                            GAM_mat[ir, 1] = res_list.data[base+2]  # gamgam
                            for jch in 1:n_non_gamma
                                GAM_mat[ir, 1+jch] = res_list.data[base+2+jch]
                            end
                        end
                    end

                    # Read background R-matrix data if KBK > 0
                    sammy_bg = RMLSammyBG[]
                    if KBK > 0
                        sammy_bg = _read_rml_background(io, KBK)
                    end
                    # Read tabulated phase shifts if KPS > 0
                    if KPS > 0
                        _skip_rml_phaseshifts(io, KPS)
                    end

                    push!(spin_groups, RMLSpinGroup(AJ, PJ, KBK, KPS,
                                                    channels, Er_vec, GAM_mat,
                                                    sammy_bg))
                end

                seek(io, saved_pos)
                return RMLData(EL, EH, IFG, KRL, KRM, AWR, pairs, spin_groups)
            else
                # Skip this range
                while !eof(io)
                    pos = position(io)
                    line = readline(io)
                    p = rpad(line, 80)
                    mf = _parse_int(p[71:72])
                    mt = _parse_int(p[73:75])
                    if mf == 0 || mt == 0 || mf != 2 || mt != 151
                        seek(io, pos)
                        break
                    end
                end
            end
        end
    end

    seek(io, saved_pos)
    return nothing
end

# Read background R-matrix records for one spin group.
# KBK is the number of background channels.
# Each background channel has:
#   CONT: ED, EU, 0, 0, LBK, 0
#   Then data depending on LBK:
#     LBK=1: LIST (tabulated R-matrix)
#     LBK=2: LIST (Sammy parameterization, 5 values)
#     LBK=3: TAB1 (logarithmic parameterization)
# Returns a vector of RMLSammyBG for LBK=2 records.
function _read_rml_background(io::IO, KBK::Int)
    result = RMLSammyBG[]
    for _ in 1:KBK
        bk_header = read_cont(io)
        LCH = Int(bk_header.L1)  # channel number
        LBK = Int(bk_header.L2)  # background type
        if LBK == 2
            bk_list = read_list(io)
            d = bk_list.data
            ED = bk_list.C1   # lower energy boundary
            EU = bk_list.C2   # upper energy boundary
            p1 = length(d) >= 1 ? d[1] : 0.0
            p2 = length(d) >= 2 ? d[2] : 0.0
            p3 = length(d) >= 3 ? d[3] : 0.0
            p4 = length(d) >= 4 ? d[4] : 0.0
            p5 = length(d) >= 5 ? d[5] : 0.0
            push!(result, RMLSammyBG(LCH, LBK, ED, EU, p1, p2, p3, p4, p5))
        elseif LBK == 3
            bk_list = read_list(io)
            d = bk_list.data
            ED = bk_list.C1
            EU = bk_list.C2
            p1 = length(d) >= 1 ? d[1] : 0.0
            p2 = length(d) >= 2 ? d[2] : 0.0
            p3 = length(d) >= 3 ? d[3] : 0.0
            push!(result, RMLSammyBG(LCH, LBK, ED, EU, p1, p2, p3, 0.0, 0.0))
        elseif LBK == 1
            read_tab1(io)  # real part table
            read_tab1(io)  # imaginary part table
        end
    end
    return result
end

# Skip phase-shift records for one spin group.
# KPS is the number of phase-shift channels.
# Each phase-shift channel has:
#   CONT: 0, 0, 0, 0, LPS, 0
#   Then data depending on LPS:
#     LPS=1: LIST (hard-sphere)
#     LPS=2: TAB1 (tabulated phase shifts)
function _skip_rml_phaseshifts(io::IO, KPS::Int)
    for _ in 1:KPS
        ps_header = read_cont(io)
        LPS = Int(ps_header.L1)
        if LPS == 1
            read_list(io)
        elseif LPS == 2
            read_tab1(io)
        end
    end
end

"""
    cross_section_rml(E, rml::RMLData) -> CrossSections

Evaluate R-Matrix Limited (LRF=7, KRM=3 Reich-Moore) cross sections
using the SAMMY method.

Follows NJOY2016 samm.f90 subroutines: betset, abpart, setr, yinvrs,
setxqx, sectio, crosss.

For each spin group:
1. Compute reduced width amplitudes betapr from widths and P(Er)
2. Form beta products: beta_kl = betapr_k * betapr_l
3. Compute alphar/alphai from the Cauchy-eliminated capture channel
4. Build R-matrix from beta * alpha
5. Form Y = L^{-1} - R where L^{-1} = 1/(S-B+iP)
6. Invert Y to get Y^{-1}
7. Form XXXX = sqrt(P)/L * Y^{-1} * R * sqrt(P)
8. Extract cross sections from XXXX via the collision matrix formalism
"""
function cross_section_rml(E::Float64, rml::RMLData)
    C = PhysicsConstants

    sig_elastic = 0.0
    sig_capture = 0.0
    sig_fission = 0.0

    # --- Wavenumber setup ---
    # Find elastic particle pair to get masses and target spin
    spi = 0.0
    ema_el = 1.0   # neutron mass ratio (amu)
    emb_el = rml.AWR
    for pp in rml.pairs
        if pp.MT == 2
            spi = pp.IB
            ema_el = pp.MA != 0.0 ? pp.MA : 1.0
            emb_el = pp.MB != 0.0 ? pp.MB : rml.AWR
            break
        end
    end

    # Wavenumber constant following NJOY samm.f90 fxradi
    hbarrr = C.hbar / C.ev                          # hbar in eV*s
    amuevv = C.amu * C.clight * C.clight / C.ev      # amu in eV
    cspeed = C.clight / 100.0                         # c in m/s
    twomhb = sqrt(2.0 * C.amassn * amuevv) / (hbarrr * 1.0e15 * cspeed)

    factor_lab = emb_el / (emb_el + ema_el)
    alabcm = factor_lab
    factor = alabcm / ema_el

    fourpi = 4.0 * C.pi / 100.0  # normalization: 4*pi / 100 (fm^2 -> barns)
    if E <= 0.0
        return CrossSections(0.0, 0.0, 0.0, 0.0)
    end

    # Loop over spin groups
    for sg in rml.spin_groups
        nch_all = length(sg.channels)
        nch_all == 0 && continue
        nres = length(sg.Er)

        ajc = abs(sg.AJ)

        # Separate channels: identify non-gamma channels and their properties
        # In LRF=7, channels in the spin group include gamma (IPP corresponding
        # to MT=102) and non-gamma channels. The gamma channel is "eliminated"
        # via the Cauchy approach (its width goes into gamgam).
        #
        # Our reader already separated: GAM column 1 = gamgam,
        # GAM columns 2+ = non-gamma channel widths.
        # The channels vector contains ALL channels from the ENDF data.
        # We need to identify which are elastic (IPP -> MT=2) etc.

        # Build list of non-gamma channels with their particle-pair info
        nongamma_chs = Int[]
        for (ich, ch) in enumerate(sg.channels)
            if ch.ipp >= 1 && ch.ipp <= length(rml.pairs)
                if rml.pairs[ch.ipp].MT != 102  # not gamma/capture
                    push!(nongamma_chs, ich)
                end
            end
        end
        nchan = length(nongamma_chs)
        nchan == 0 && continue

        # Determine entrance channels (elastic, MT=2) and exit channels
        nent = 0   # number of entrance channels
        for idx in nongamma_chs
            ch = sg.channels[idx]
            if rml.pairs[ch.ipp].MT == 2
                nent += 1
            end
        end
        nent == 0 && continue

        # Statistical weight factor: goj = (2J+1) / ((2*IA+1)*(2*IB+1))
        # Use the first entrance channel's particle pair for spins
        first_ent_pp = 0
        for idx in nongamma_chs
            ch = sg.channels[idx]
            if rml.pairs[ch.ipp].MT == 2
                first_ent_pp = ch.ipp
                break
            end
        end
        ia = abs(rml.pairs[first_ent_pp].IA)
        ib = abs(rml.pairs[first_ent_pp].IB)
        goj = (2.0 * ajc + 1.0) / ((2.0 * ia + 1.0) * (2.0 * ib + 1.0))

        # --- Compute channel-specific wavenumber parameters ---
        # Following fxradi in samm.f90
        zke_ch = zeros(nchan)   # k factor for each channel
        zkfe_ch = zeros(nchan)  # k * ape (effective)
        zkte_ch = zeros(nchan)  # k * apt (true)
        echan_ch = zeros(nchan) # channel threshold energy
        lpent_ch = zeros(Int, nchan)
        ishift_ch = zeros(Int, nchan)
        ll_ch = zeros(Int, nchan)
        ape_ch = zeros(nchan)
        apt_ch = zeros(nchan)
        bnd_ch = zeros(nchan)

        for (ic, idx) in enumerate(nongamma_chs)
            ch = sg.channels[idx]
            pp = rml.pairs[ch.ipp]
            ll_ch[ic] = ch.l
            ape_ch[ic] = ch.ape
            apt_ch[ic] = ch.apt
            bnd_ch[ic] = ch.bnd
            lpent_ch[ic] = Int(pp.PNT)
            ishift_ch[ic] = Int(pp.SHF)

            ema_pp = pp.MA != 0.0 ? pp.MA : 1.0
            emb_pp = pp.MB != 0.0 ? pp.MB : rml.AWR
            aa = emb_pp / (emb_pp + ema_pp)

            # Channel threshold from Q value
            if pp.Q != 0.0
                echan_ch[ic] = -pp.Q / alabcm
            end

            redmas = aa * ema_pp
            z = twomhb * sqrt(redmas * factor)
            zke_ch[ic] = z
            zkfe_ch[ic] = z * ch.ape * 10.0  # convert fm to 10^-15 m
            zkte_ch[ic] = z * ch.apt * 10.0
        end

        # --- Compute betapr (reduced width amplitudes) for each resonance ---
        # Following betset in samm.f90
        betapr = zeros(nchan, nres)
        gbetpr2 = zeros(nres)   # |gamgam|/2
        gbetpr3 = zeros(nres)   # (gamgam/2)^2

        for ires in 1:nres
            gg = sg.GAM[ires, 1]  # gamgam (capture width)
            gbetpr2[ires] = abs(gg) / 2.0
            gbetpr3[ires] = gbetpr2[ires]^2

            for ic in 1:nchan
                gam_val = sg.GAM[ires, 1 + ic]  # channel width
                gam_val == 0.0 && continue

                if lpent_ch[ic] <= 0
                    # No penetrability calculation
                    betapr[ic, ires] = sqrt(0.5 * abs(gam_val))
                    if gam_val < 0.0
                        betapr[ic, ires] = -betapr[ic, ires]
                    end
                else
                    # Penetrability at resonance energy
                    er = sg.Er[ires]
                    ex = abs(er - echan_ch[ic])
                    if ex != 0.0
                        ex = sqrt(ex)
                        rho = zkte_ch[ic] * ex
                        lsp = ll_ch[ic]
                        p_er = penetrability(lsp, rho)
                        if p_er <= 0.0
                            p_er = 1.0
                        end
                        betapr[ic, ires] = sqrt(0.5 * abs(gam_val) / p_er)
                        if gam_val < 0.0
                            betapr[ic, ires] = -betapr[ic, ires]
                        end
                    end
                end
            end
        end

        # --- Compute beta products: beta_kl = betapr_k * betapr_l ---
        # Stored in lower-triangular packed form
        ntriag = div(nchan * (nchan + 1), 2)
        beta = zeros(ntriag, nres)
        for ires in 1:nres
            kl = 0
            for ik in 1:nchan
                for il in 1:ik
                    kl += 1
                    beta[kl, ires] = betapr[il, ires] * betapr[ik, ires]
                end
            end
        end

        # --- Energy-dependent part: alphar, alphai ---
        # Following abpart in samm.f90
        # alphar = (Er - E) / ((Er - E)^2 + (gamgam/2)^2)
        # alphai = (gamgam/2) / ((Er - E)^2 + (gamgam/2)^2)
        alphar_v = zeros(nres)
        alphai_v = zeros(nres)
        for ires in 1:nres
            diff = sg.Er[ires] - E
            den = diff * diff + gbetpr3[ires]
            if den != 0.0
                alphar_v[ires] = diff / den
                alphai_v[ires] = gbetpr2[ires] / den
            end
        end

        # --- Build R-matrix (lower-triangular packed complex) ---
        # Following setr in samm.f90
        rmat_r = zeros(ntriag)
        rmat_i = zeros(ntriag)

        # Add background R-matrix elements on diagonal
        for bg in sg.sammy_bg
            # bg.LCH is the channel number (1-based, includes gamma)
            # In the ENDF data, channel 1 is gamma. So non-gamma channel
            # index in our array is bg.LCH - 1 (since we stripped gamma).
            bg_ch = bg.LCH - 1
            if bg_ch < 1 || bg_ch > nchan
                continue
            end
            # Diagonal index in packed lower-triangular
            diag_idx = div(bg_ch * (bg_ch + 1), 2)

            if bg.LBK == 2
                # Sammy parametrization
                r_bg = bg.p1
                if bg.p5 != 0.0
                    r_bg -= bg.p5 * (bg.EU - bg.ED)
                end
                r_bg += (bg.p2 + bg.p3 * E) * E
                if bg.p4 != 0.0 || bg.p5 != 0.0
                    num_v = bg.EU - E
                    den_v = E - bg.ED
                    if num_v > 0.0 && den_v > 0.0
                        r_bg -= (bg.p4 + bg.p5 * E) * log(num_v / den_v)
                    end
                end
                rmat_r[diag_idx] += r_bg
            elseif bg.LBK == 3
                # Frohner parametrization
                esum = bg.ED + bg.EU
                ediff = bg.EU - bg.ED
                rmat_r[diag_idx] += bg.p1 + 2.0 * bg.p2 * atanh((2.0 * E - esum) / ediff)
                rmat_i[diag_idx] += bg.p3 / ediff / (1.0 - ((2.0 * E - esum) / ediff)^2)
            end
        end

        # Add resonance contributions
        for ires in 1:nres
            kl = 0
            for ik in 1:nchan
                for il in 1:ik
                    kl += 1
                    b = beta[kl, ires]
                    if E > echan_ch[ik] && E > echan_ch[il] && b != 0.0
                        rmat_r[kl] += alphar_v[ires] * b
                        rmat_i[kl] += alphai_v[ires] * b
                    end
                end
            end
        end

        # Check if R-matrix is all zeros
        lrmat = all(rmat_r .== 0.0) && all(rmat_i .== 0.0) ? 1 : 0

        # --- Build Y-matrix = L^{-1} - R, compute penetrability/phase ---
        # Following setr in samm.f90
        ymat_r = copy(-rmat_r)
        ymat_i = copy(-rmat_i)

        rootp_ch = ones(nchan)
        elinvr_ch = zeros(nchan)
        elinvi_ch = fill(-1.0, nchan)
        sinsqr_ch = zeros(nchan)
        sin2ph_ch = zeros(nchan)

        for ic in 1:nchan
            ii = div(ic * (ic + 1), 2)  # diagonal index
            if E > echan_ch[ic]
                if lpent_ch[ic] <= 0
                    # No penetrability calculation
                    ymat_i[ii] -= 1.0
                else
                    ex = sqrt(E - echan_ch[ic])
                    rho = zkte_ch[ic] * ex
                    rhof = zkfe_ch[ic] * ex
                    lsp = ll_ch[ic]

                    # Phase shift (from effective radius)
                    phi_ch = phase_shift(lsp, rhof)
                    sp = sin(phi_ch)
                    cp = cos(phi_ch)
                    sinsqr_ch[ic] = sp * sp
                    sin2ph_ch[ic] = 2.0 * cp * sp

                    # Penetrability and shift (from true radius)
                    p_val = penetrability(lsp, rho)
                    s_val = shift_factor(lsp, rho)

                    rootp_ch[ic] = sqrt(p_val)

                    # L^{-1} = 1/(S - B + iP) where B = boundary condition
                    bnd = bnd_ch[ic]
                    gg_val = ishift_ch[ic] > 0 ? s_val - bnd : -bnd
                    hh_val = p_val

                    # Compute g + ih = 1/(gg + i*hh) = (gg - i*hh)/(gg^2+hh^2)
                    if hh_val <= 1.0e-35
                        hh_val = 0.0
                    end

                    if gg_val == 0.0 && hh_val == 0.0
                        # iffy case: very small, skip
                        elinvr_ch[ic] = 0.0
                        elinvi_ch[ic] = -1.0
                        ymat_i[ii] -= 1.0
                    elseif gg_val == 0.0
                        elinvr_ch[ic] = 0.0
                        elinvi_ch[ic] = -1.0 / hh_val
                        ymat_r[ii] += 0.0
                        ymat_i[ii] += -1.0 / hh_val
                    elseif hh_val == 0.0
                        elinvr_ch[ic] = 1.0 / gg_val
                        elinvi_ch[ic] = 0.0
                        ymat_r[ii] += 1.0 / gg_val
                        ymat_i[ii] += 0.0
                    else
                        d = hh_val^2 + gg_val^2
                        g_inv = gg_val / d
                        h_inv = -hh_val / d
                        elinvr_ch[ic] = g_inv
                        elinvi_ch[ic] = h_inv
                        ymat_r[ii] += g_inv
                        ymat_i[ii] += h_inv
                    end
                end
            else
                # Below threshold: channel is closed
                rootp_ch[ic] = 0.0
                elinvr_ch[ic] = 1.0
                elinvi_ch[ic] = 0.0
            end
        end

        # Ensure diagonal elements are not zero (for inversion)
        if lrmat != 1
            for ic in 1:nchan
                ii = div(ic * (ic + 1), 2)
                if ymat_r[ii] == 0.0 && ymat_i[ii] == 0.0
                    ymat_r[ii] = 1.0
                end
            end
        end

        # --- Invert Y-matrix ---
        # If R-matrix is zero, XXXX is zero (no resonance contribution)
        xxxxr = zeros(ntriag)
        xxxxi = zeros(ntriag)

        if lrmat != 1
            # Build complex Y-matrix (full symmetric) for inversion
            Y = zeros(ComplexF64, nchan, nchan)
            for ik in 1:nchan
                for il in 1:ik
                    idx = div(ik * (ik - 1), 2) + il
                    Y[ik, il] = complex(ymat_r[idx], ymat_i[idx])
                    Y[il, ik] = Y[ik, il]  # symmetric
                end
            end

            # Invert Y
            Yinv = inv(Y)

            # Build complex R-matrix (full symmetric) for multiplication
            R = zeros(ComplexF64, nchan, nchan)
            for ik in 1:nchan
                for il in 1:ik
                    idx = div(ik * (ik - 1), 2) + il
                    R[ik, il] = complex(rmat_r[idx], rmat_i[idx])
                    R[il, ik] = R[ik, il]
                end
            end

            # XQ = Yinv * R
            XQ = Yinv * R

            # XXXX = sqrt(P)/L * XQ * sqrt(P) (symmetric, packed)
            for ii in 1:nchan
                plr = rootp_ch[ii] * elinvr_ch[ii]
                pli = rootp_ch[ii] * elinvi_ch[ii]
                pl = complex(plr, pli)
                for jj in 1:ii
                    idx = div(ii * (ii - 1), 2) + jj
                    val = rootp_ch[jj] * pl * XQ[jj, ii]
                    xxxxr[idx] = real(val)
                    xxxxi[idx] = imag(val)
                end
            end
        end

        # --- Extract cross sections from XXXX (following sectio) ---
        crss_elastic = 0.0
        crss_absorb = 0.0

        # Elastic: crss(1) = g * [sin^2(phi)*(1-2*XXXX_i) - sin(2phi)*XXXX_r + |XXXX|^2] / zke^2
        # Absorption: crss(2) = g * [XXXX_i - |XXXX|^2] / zke^2
        ient = 0
        for ic in 1:nchan
            ch = sg.channels[nongamma_chs[ic]]
            if rml.pairs[ch.ipp].MT != 2
                continue
            end
            ient += 1

            zz = zke_ch[ic]^2
            ii_diag = div(ic * (ic + 1), 2)

            # Elastic term (from entrance channel)
            termn = sinsqr_ch[ic] * (1.0 - 2.0 * xxxxi[ii_diag]) -
                    sin2ph_ch[ic] * xxxxr[ii_diag]
            termn /= zz

            # Sum |XXXX_ij|^2 for all j <= i (entrance channels)
            for jc in 1:ic
                ij = div(ic * (ic - 1), 2) + jc
                ar = (xxxxr[ij]^2 + xxxxi[ij]^2) / zz
                if ic != jc
                    ar += ar
                end
                termn += ar
            end
            crss_elastic += termn

            # Absorption term
            terma = xxxxi[ii_diag] / zz
            for jc in 1:ic
                ij = div(ic * (ic - 1), 2) + jc
                ar = -(xxxxr[ij]^2 + xxxxi[ij]^2) / zz
                if ic != jc
                    ar += ar
                end
                terma += ar
            end
            crss_absorb += terma
        end
        crss_elastic *= goj
        crss_absorb *= goj

        # Other reaction channels (fission etc.)
        crss_other = zeros(length(rml.pairs))
        next_count = 0
        for ic in 1:nchan
            ch = sg.channels[nongamma_chs[ic]]
            mt_ch = rml.pairs[ch.ipp].MT
            if mt_ch == 2 || mt_ch == 102
                continue
            end
            # This is an exit-only channel
            for jc in 1:nchan
                ch_j = sg.channels[nongamma_chs[jc]]
                if rml.pairs[ch_j.ipp].MT != 2
                    continue
                end
                # jc is entrance, ic is exit
                if ic > jc
                    ij = div(ic * (ic - 1), 2) + jc
                else
                    ij = div(jc * (jc - 1), 2) + ic
                end
                zz = zke_ch[jc]^2
                crss_other[ch.ipp] += (xxxxr[ij]^2 + xxxxi[ij]^2) / zz
            end
        end
        for ipp in 1:length(rml.pairs)
            crss_other[ipp] *= goj
        end

        # Map to cross section components
        # sigmas(1) = elastic, sigmas(2) = absorption (capture + other)
        # Normalization: fourpi / E (applied later)
        norm = fourpi / E
        sig_elastic += crss_elastic * norm
        sig_capture += crss_absorb * norm

        # Subtract other channel contributions from capture
        for ipp in 1:length(rml.pairs)
            mt_pp = rml.pairs[ipp].MT
            if mt_pp == 18
                sig_fission += crss_other[ipp] * norm
                sig_capture -= crss_other[ipp] * norm
            elseif mt_pp != 2 && mt_pp != 102
                sig_capture -= crss_other[ipp] * norm
            end
        end
    end

    sig_total = sig_elastic + sig_capture + sig_fission

    return CrossSections(sig_total, sig_elastic, sig_fission, sig_capture)
end

"""
    build_rml_evaluator(rml::RMLData) -> Function

Build an evaluator closure for RML data that returns
`(total, elastic, fission, capture)` at a given energy.
"""
function build_rml_evaluator(rml::RMLData)
    function evaluator(E::Float64)
        if E < rml.EL || E >= rml.EH
            return (0.0, 0.0, 0.0, 0.0)
        end
        xs = cross_section_rml(E, rml)
        return (max(0.0, xs.total), max(0.0, xs.elastic),
                max(0.0, xs.fission), max(0.0, xs.capture))
    end
    return evaluator
end

function merge_background_legacy(energies::Vector{Float64},
                                  res_xs::Vector{<:CrossSections},
                                  mf3_sections::Vector{MF3Section})
    _SKIP_MTS_LEGACY = (1, 3, 4, 101, 27)
    n = length(energies)
    result = Vector{CrossSections{Float64}}(undef, n)
    for i in 1:n
        e = energies[i]
        elastic = res_xs[i].elastic
        fission = res_xs[i].fission
        capture = res_xs[i].capture
        other_bg = 0.0
        for sec in mf3_sections
            mt = Int(sec.mt)
            mt in _SKIP_MTS_LEGACY && continue
            bg = interpolate(sec.tab, e)
            bg == 0.0 && continue
            if mt == 2
                elastic += bg
            elseif mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38
                fission += bg
            elseif mt == 102
                capture += bg
            else
                other_bg += bg
            end
        end
        if elastic <= 1.0e-8
            elastic = 1.0e-8
        end
        total = elastic + fission + capture + other_bg
        result[i] = CrossSections(total, elastic, fission, capture)
    end
    return result
end
