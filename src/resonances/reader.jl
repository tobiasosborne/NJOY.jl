# MF2 (File 2) resonance parameter reader
#
# Reads ENDF File 2 / MT151 resonance parameters and populates
# the typed resonance structures defined in types.jl.
#
# Handles:
#   LRU=0: scattering radius only (no resonances)
#   LRU=1, LRF=1: Single-Level Breit-Wigner (SLBW)
#   LRU=1, LRF=2: Multi-Level Breit-Wigner (MLBW)
#   LRU=1, LRF=3: Reich-Moore (RM)
#   NRO=1: energy-dependent scattering radius (TAB1)
#
# Reference: NJOY2016 reconr.f90 subroutines rdfil2, rdf2bw
# Reference: ENDF-102 Section 2 (File 2 format specification)

"""
    IsotopeData

Container for one isotope's resonance data from MF2/MT151.
"""
struct IsotopeData
    ZAI::Float64       # ZA of isotope (e.g. 26056.0 for Fe-56)
    ABN::Float64       # isotopic abundance (fractional)
    LFW::Int32         # fission width flag (for URR)
    ranges::Vector{ResonanceRange}
end

"""
    MF2Data

Container for all File 2 resonance data for a material.
"""
struct MF2Data
    ZA::Float64        # ZA of material
    AWR::Float64       # atomic weight ratio to neutron
    isotopes::Vector{IsotopeData}
end

"""
    read_mf2(io::IO) -> MF2Data

Read the entire MF2/MT151 section from an ENDF file.
The stream `io` must be positioned at the start of the MF2/MT151 HEAD record
(i.e., after calling `find_section(io, 2, 151)`).

Returns an `MF2Data` struct containing all isotope data and resonance ranges.
"""
function read_mf2(io::IO)
    # HEAD record: ZA, AWR, 0, 0, NIS, 0
    head = read_head(io)
    ZA = head.C1
    AWR = head.C2
    NIS = Int(head.N1)

    isotopes = IsotopeData[]

    for _ in 1:NIS
        # Isotope CONT: ZAI, ABN, 0, LFW, NER, 0
        iso_cont = read_cont(io)
        ZAI = iso_cont.C1
        ABN = iso_cont.C2
        LFW = iso_cont.L2
        NER = Int(iso_cont.N1)

        ranges = ResonanceRange[]

        for _ in 1:NER
            # Range CONT: EL, EH, LRU, LRF, NRO, NAPS
            range_cont = read_cont(io)
            EL = range_cont.C1
            EH = range_cont.C2
            LRU = range_cont.L1
            LRF = range_cont.L2
            NRO = range_cont.N1
            NAPS = range_cont.N2

            if LRU == 0
                # Scattering radius only -- no resonances
                # Read the parameters CONT: SPI, AP, 0, 0, NLS, 0
                params_cont = read_cont(io)
                SPI = params_cont.C1
                AP = params_cont.C2
                # Create a minimal SLBW parameter set with no resonances
                params = SLBWParameters(
                    Int32(0), SPI, AP,
                    Int32[], Float64[], Float64[], Int32[],
                    Vector{Float64}[], Vector{Float64}[],
                    Vector{Float64}[], Vector{Float64}[],
                    Vector{Float64}[], Vector{Float64}[]
                )
                push!(ranges, ResonanceRange(EL, EH, LRU, Int32(0), LFW, NRO, NAPS, params))

            elseif LRU == 1 && LRF <= 2
                # SLBW (LRF=1) or MLBW (LRF=2) -- same reading format
                params, ap_tab = _read_bw_params(io, LRF, NRO, NAPS)
                push!(ranges, ResonanceRange(EL, EH, LRU, LRF, LFW, NRO, NAPS, params, ap_tab))

            elseif LRU == 1 && LRF == 3
                # Reich-Moore
                params, ap_tab = _read_rm_params(io, NRO, NAPS)
                push!(ranges, ResonanceRange(EL, EH, LRU, LRF, LFW, NRO, NAPS, params, ap_tab))

            elseif LRU == 1 && LRF == 7
                # SAMMY R-Matrix Limited
                params = _read_sammy_params(io, NAPS)
                push!(ranges, ResonanceRange(EL, EH, LRU, LRF, LFW, NRO, NAPS, params))

            elseif LRU == 2 && LRF <= 1 && LFW == 1
                # Unresolved resonances, mode=11 (LRF=1, LFW=1)
                # Energy-dependent fission widths
                params = _read_urr_lfw1(io, Int(NRO), Int(NAPS))
                push!(ranges, ResonanceRange(EL, EH, LRU, LRF, LFW, NRO, NAPS, params))

            else
                # Unsupported formalisms (Adler-Adler, other URR variants, etc.)
                _skip_unsupported_range(io)
                # Push a placeholder so _resonance_bounds sees the EL/EH
                params = SLBWParameters(
                    Int32(0), 0.0, 0.0,
                    Int32[], Float64[], Float64[], Int32[],
                    Vector{Float64}[], Vector{Float64}[],
                    Vector{Float64}[], Vector{Float64}[],
                    Vector{Float64}[], Vector{Float64}[]
                )
                push!(ranges, ResonanceRange(EL, EH, LRU, LRF, LFW, NRO, NAPS, params))
            end
        end

        push!(isotopes, IsotopeData(ZAI, ABN, LFW, ranges))
    end

    return MF2Data(ZA, AWR, isotopes)
end

"""
Read SLBW (LRF=1) or MLBW (LRF=2) resonance parameters.
Same ENDF format for both; the difference is in cross section evaluation.

Format (from ENDF-102 and NJOY's rdf2bw):
  [optional TAB1 for energy-dependent scattering radius if NRO=1]
  CONT: SPI, AP, 0, 0, NLS, 0
  For each L value:
    LIST header: AWRI, QX, L, LRX, 6*NRS, NRS
    LIST data: (Er, AJ, GT, GN, GG, GF) for each resonance
"""
function _read_bw_params(io::IO, LRF::Int32, NRO::Int32, NAPS::Int32)
    # Energy-dependent scattering radius (if NRO != 0)
    ap_tab = nothing
    if NRO != 0
        ap_tab = TabulatedFunction(read_tab1(io))
    end

    # Parameters CONT: SPI, AP, 0, 0, NLS, 0
    params_cont = read_cont(io)
    SPI = params_cont.C1
    AP = params_cont.C2
    NLS = Int(params_cont.N1)

    l_values = Int32[]
    AWRI_all = Float64[]
    QX_all = Float64[]
    LRX_all = Int32[]
    Er_all = Vector{Float64}[]
    AJ_all = Vector{Float64}[]
    Gn_all = Vector{Float64}[]
    Gg_all = Vector{Float64}[]
    Gf_all = Vector{Float64}[]
    Gx_all = Vector{Float64}[]

    for _ in 1:NLS
        # LIST record: AWRI, QX, L, LRX, 6*NRS, NRS
        # followed by data lines containing resonance parameters
        lst = read_list(io)
        AWRI = lst.C1
        QX = lst.C2
        L = lst.L1
        LRX = lst.L2
        NRS = Int(lst.N2)

        push!(l_values, L)
        push!(AWRI_all, AWRI)
        push!(QX_all, QX)
        push!(LRX_all, LRX)

        Er = Float64[]
        AJ = Float64[]
        Gn = Float64[]
        Gg = Float64[]
        Gf = Float64[]
        Gx = Float64[]

        for ir in 1:NRS
            base = (ir - 1) * 6
            er = lst.data[base + 1]   # resonance energy
            aj = lst.data[base + 2]   # spin J
            gt = lst.data[base + 3]   # total width (GT for SLBW/MLBW)
            gn = lst.data[base + 4]   # neutron width
            gg = lst.data[base + 5]   # gamma (capture) width
            gf = lst.data[base + 6]   # fission width

            # Competitive width: Gx = GT - GN - GG - GF
            gx = gt - gn - gg - gf
            if gx < 0.0
                gx = 0.0
            end

            push!(Er, er)
            push!(AJ, aj)
            push!(Gn, gn)
            push!(Gg, gg)
            push!(Gf, gf)
            push!(Gx, gx)
        end

        push!(Er_all, Er)
        push!(AJ_all, AJ)
        push!(Gn_all, Gn)
        push!(Gg_all, Gg)
        push!(Gf_all, Gf)
        push!(Gx_all, Gx)
    end

    if LRF == Int32(1)
        params = SLBWParameters(Int32(NLS), SPI, AP,
                              l_values, AWRI_all, QX_all, LRX_all,
                              Er_all, AJ_all,
                              Gn_all, Gg_all, Gf_all, Gx_all)
    else
        params = MLBWParameters(Int32(NLS), SPI, AP,
                              l_values, AWRI_all, QX_all, LRX_all,
                              Er_all, AJ_all,
                              Gn_all, Gg_all, Gf_all, Gx_all)
    end
    return params, ap_tab
end

"""
Read Reich-Moore (LRF=3) resonance parameters.

Format (from ENDF-102 and NJOY's rdf2bw):
  [optional TAB1 for energy-dependent scattering radius if NRO=1]
  CONT: SPI, AP, LAD, 0, NLS, NLSC
  For each L value:
    LIST header: AWRI, APL, L, 0, 6*NRS, NRS
    LIST data: (Er, AJ, GN, GG, GFA, GFB) for each resonance

Note: For Reich-Moore, the 6 parameters per resonance are:
  Er, AJ, GN, GG, GFA, GFB (not GT, GN, GG, GF as for BW)
"""
function _read_rm_params(io::IO, NRO::Int32, NAPS::Int32)
    # Energy-dependent scattering radius (if NRO != 0)
    ap_tab = nothing
    if NRO != 0
        ap_tab = TabulatedFunction(read_tab1(io))
    end

    # Parameters CONT: SPI, AP, LAD, 0, NLS, NLSC
    params_cont = read_cont(io)
    SPI = params_cont.C1
    AP = params_cont.C2
    LAD = params_cont.L1
    NLS = Int(params_cont.N1)

    l_values = Int32[]
    AWRI_all = Float64[]
    APL_all = Float64[]
    Er_all = Vector{Float64}[]
    AJ_all = Vector{Float64}[]
    Gn_all = Vector{Float64}[]
    Gg_all = Vector{Float64}[]
    Gfa_all = Vector{Float64}[]
    Gfb_all = Vector{Float64}[]

    for _ in 1:NLS
        # LIST record: AWRI, APL, L, 0, 6*NRS, NRS
        lst = read_list(io)
        AWRI = lst.C1
        APL_val = lst.C2
        L = lst.L1
        NRS = Int(lst.N2)

        push!(l_values, L)
        push!(AWRI_all, AWRI)
        push!(APL_all, APL_val)

        Er = Float64[]
        AJ = Float64[]
        Gn = Float64[]
        Gg = Float64[]
        Gfa = Float64[]
        Gfb = Float64[]

        for ir in 1:NRS
            base = (ir - 1) * 6
            er = lst.data[base + 1]   # resonance energy
            aj = lst.data[base + 2]   # spin J (sign encodes channel spin)
            gn = lst.data[base + 3]   # neutron width
            gg = lst.data[base + 4]   # gamma (capture) width
            gfa = lst.data[base + 5]  # first fission width
            gfb = lst.data[base + 6]  # second fission width

            push!(Er, er)
            push!(AJ, aj)
            push!(Gn, gn)
            push!(Gg, gg)
            push!(Gfa, gfa)
            push!(Gfb, gfb)
        end

        push!(Er_all, Er)
        push!(AJ_all, AJ)
        push!(Gn_all, Gn)
        push!(Gg_all, Gg)
        push!(Gfa_all, Gfa)
        push!(Gfb_all, Gfb)
    end

    params = ReichMooreParameters(Int32(NLS), SPI, AP, LAD,
                                l_values, AWRI_all, APL_all,
                                Er_all, AJ_all,
                                Gn_all, Gg_all, Gfa_all, Gfb_all)
    return params, ap_tab
end

"""
Read SAMMY R-Matrix Limited (LRF=7) resonance parameters.

Format (from ENDF-102 and NJOY's rdsammy for mode=7):
  CONT: SPI, AP, IFG, KRM, NJS, 0    (NJS = number of spin groups)
  LIST: particle pairs (NPP pairs, each 12 values)
  For each spin group:
    LIST: spin group header + channel definitions
    LIST: resonance parameters (Er, gamgam, gamma_1, ..., gamma_nchan)
    [optional: KBK background R-matrix CONT+LIST records per channel]
"""
function _read_sammy_params(io::IO, NAPS::Int32)
    # Parameters CONT: SPI, AP, IFG, KRM, NJS, 0
    params_cont = read_cont(io)
    SPI = params_cont.C1
    AP = params_cont.C2
    IFG = params_cont.L1
    KRM = params_cont.L2
    NJS = Int(params_cont.N1)  # number of spin groups

    # Particle-pair LIST record
    pp_lst = read_list(io)
    NPP = Int(pp_lst.L1)

    particle_pairs = SAMMYParticlePair[]
    for i in 1:NPP
        base = (i - 1) * 12
        ema   = pp_lst.data[base + 1]
        emb   = pp_lst.data[base + 2]
        kza   = Int32(round(pp_lst.data[base + 3]))
        kzb   = Int32(round(pp_lst.data[base + 4]))
        spina = pp_lst.data[base + 5]
        spinb = pp_lst.data[base + 6]
        qqq   = pp_lst.data[base + 7]
        lpent = Int32(round(pp_lst.data[base + 8]))
        ishift = Int32(round(pp_lst.data[base + 9]))
        mt    = Int32(round(pp_lst.data[base + 10]))
        pa    = pp_lst.data[base + 11]
        pb    = pp_lst.data[base + 12]
        # Apply default masses from particle-pair definitions
        if mt == 2  # neutron elastic
            if ema == 0.0; ema = 1.0; end
            if spina == 0.0; spina = 0.5; end
            lpent = Int32(1)
        elseif mt == 18  # fission
            lpent = Int32(0)
        elseif mt == 102  # capture
            # gamma channel: defaults are zero masses
        end
        push!(particle_pairs, SAMMYParticlePair(ema, emb, kza, kzb, spina, spinb,
                                                qqq, lpent, ishift, mt, pa, pb))
    end

    # Find neutron particle-pair (mt=2) to get target spin from its spinb
    target_spin = SPI
    for pp in particle_pairs
        if pp.mt == 2
            target_spin = pp.spinb
            break
        end
    end

    spin_groups = SAMMYSpinGroup[]

    for ig in 1:NJS
        # Spin group header LIST: AJ, parity, KBK, 0, NCH_total*6, NCH_total
        sg_lst = read_list(io)
        AJ = sg_lst.C1
        parity = sg_lst.C2
        KBK = sg_lst.L1  # number of background R-matrix channels
        nch_plus1 = Int(sg_lst.N2)  # NCH+1 (includes eliminated capture)

        # Read channel definitions: first entry is the eliminated capture channel
        ipp_arr = Int32[]
        lspin_arr = Int32[]
        chspin_arr = Float64[]
        bound_arr = Float64[]
        rdeff_arr = Float64[]
        rdtru_arr = Float64[]
        igamma = 0  # index of eliminated capture channel within the full list

        nchan = 0
        for ich in 1:nch_plus1
            base = (ich - 1) * 6
            ippx = Int32(round(sg_lst.data[base + 1]))
            if ippx == 1  # this is the eliminated capture channel (particle pair 1)
                igamma = ich
            else
                nchan += 1
                push!(ipp_arr, ippx)
                push!(lspin_arr, Int32(round(sg_lst.data[base + 2])))
                push!(chspin_arr, sg_lst.data[base + 3])
                push!(bound_arr, sg_lst.data[base + 4])
                push!(rdeff_arr, sg_lst.data[base + 5])
                push!(rdtru_arr, sg_lst.data[base + 6])
            end
        end

        # Resonance parameter LIST: 0, 0, 0, NRES, NPL, NRES
        res_lst = read_list(io)
        nres = Int(res_lst.L2)

        eres = Float64[]
        gamgam = Float64[]
        gamma_ch = [Float64[] for _ in 1:nchan]

        if nres > 0
            # Determine record width: 6 or 12 depending on nchan
            nx = 6
            if 2 + nchan > 6
                nx = 12
            end

            for ir in 1:nres
                base = (ir - 1) * nx
                er = res_lst.data[base + 1]
                gg = res_lst.data[base + 2]

                # Channel widths are in positions 3..2+nchan
                # But if igamma != 1, we need to rearrange
                raw_widths = Float64[]
                for j in 1:(nx - 2)
                    if j <= nchan
                        push!(raw_widths, res_lst.data[base + 2 + j])
                    end
                end

                # Handle the case where igamma is not the first channel
                # (need to swap gamma channel width with capture width)
                if igamma != 1 && igamma > 0
                    # The data was stored with capture in position igamma
                    # and widths shifted. We need to un-swap.
                    # igamma-1 corresponds to a channel index after excluding capture
                    # The raw data has: [ch1, ch2, ..., gamgam_at_igamma-1, ..., chN]
                    # But the capture width gg is for the non-eliminated channel
                    # and the eliminated channel's width is at position igamma in the original list
                    if igamma - 1 <= length(raw_widths)
                        # Swap: the value at position igamma-1 in raw_widths is actually gamgam
                        # and the listed gg is actually a channel width
                        actual_gamgam = raw_widths[igamma - 1]
                        for j in (igamma - 1):-1:2
                            raw_widths[j] = raw_widths[j - 1]
                        end
                        raw_widths[1] = gg
                        gg = actual_gamgam
                    end
                end

                push!(eres, er)
                push!(gamgam, gg)
                for j in 1:nchan
                    if j <= length(raw_widths)
                        push!(gamma_ch[j], raw_widths[j])
                    else
                        push!(gamma_ch[j], 0.0)
                    end
                end
            end
        end

        # Read background R-matrix records if KBK > 0
        backgr_type = zeros(Int32, nchan)
        backgr_data = [Float64[] for _ in 1:nchan]

        if KBK > 0
            for _ in 1:KBK
                bkg_cont = read_cont(io)
                lch = bkg_cont.L1  # channel index (1-based, includes capture)
                lbk = bkg_cont.L2  # background type (0,1,2,3)

                # Convert from full channel index (including capture) to our index
                ch_idx = lch - 1  # subtract 1 because capture channel is removed
                if lch == 1
                    # This is the capture channel - skip but shouldn't happen
                    ch_idx = 0
                end

                if lbk == 0 || ch_idx <= 0
                    # No background or capture channel
                    continue
                end

                if ch_idx > nchan
                    continue
                end

                backgr_type[ch_idx] = Int32(lbk)

                if lbk == 1
                    # Tabulated background: two TAB1 records
                    tab1_r = read_tab1(io)
                    tab1_i = read_tab1(io)
                    # Store as flat array: [type=1, tab_r_data..., tab_i_data...]
                    backgr_data[ch_idx] = Float64[1.0]
                    # For now, skip tabulated background (very rare)
                elseif lbk == 2 || lbk == 3
                    # Sammy (lbk=2) or Frohner (lbk=3) parametrisation: one LIST record
                    bkg_lst = read_list(io)
                    ed = bkg_lst.C1  # ED (lower bound energy)
                    eu = bkg_lst.C2  # EU (upper bound energy)
                    vals = Float64[ed, eu]
                    for j in 1:Int(bkg_lst.N1)
                        push!(vals, bkg_lst.data[j])
                    end
                    backgr_data[ch_idx] = vals
                end
            end
        end

        push!(spin_groups, SAMMYSpinGroup(
            AJ, parity, Int32(nchan),
            ipp_arr, lspin_arr, chspin_arr, bound_arr, rdeff_arr, rdtru_arr,
            backgr_type, backgr_data,
            Int32(nres), eres, gamgam, gamma_ch
        ))
    end

    return SAMMYParameters(target_spin, AP, Int32(KRM), Int32(IFG),
                           particle_pairs, spin_groups)
end

"""
Skip an unsupported resonance range by reading lines until the next
energy range or end of MF2/MT151 section.

This works by reading lines and tracking the ENDF line format:
lines with MF=2, MT=151 belong to this section. We need to consume
all lines for this range without knowing the exact format.

Strategy: read CONT/LIST/TAB1 records according to the structure.
For unresolved ranges, we read the parameter CONT and then consume
the remaining list records.
"""
function _skip_unsupported_range(io::IO)
    # The range CONT (EL, EH, LRU, LRF, NRO, NAPS) has already been read.
    # We need to consume the rest of this range.
    # Read lines until we reach a line that starts a new range or is a SEND record.
    # The safest approach: peek at lines and consume them.
    # Since we don't know the exact structure, we'll read line by line
    # checking the MF/MT fields.

    # Save position and read ahead to find the end of this range
    while !eof(io)
        pos = position(io)
        line = readline(io)
        p = rpad(line, 80)
        mf = _parse_int(p[71:72])
        mt = _parse_int(p[73:75])

        # If MF=0 or MT=0, we've hit a section/file end -- seek back
        if mf == 0 || mt == 0
            seek(io, pos)
            return
        end

        # If we're still in MF2/MT151, keep reading
        if mf != 2 || mt != 151
            seek(io, pos)
            return
        end

        # Otherwise keep consuming this line (it's part of this range)
    end
end

"""
Read LRU=2, LRF<=1, LFW=1 unresolved data (mode=11 with energy-dependent
fission widths). Matches Fortran rdf2u1 (reconr.f90:1312-1425).

The data is a LIST record: header has SPI, AP, LSSF in C1/C2/L1,
NE (energy count) in N1, NLS (l-value count) in N2. Body has NE energies.
Then per-l CONT + per-J LIST records follow.
"""
function _read_urr_lfw1(io::IO, NRO::Int, NAPS::Int)
    # LIST record: C1=SPI, C2=AP, L1=LSSF, L2=0, N1=NE, N2=NLS
    # Matching Fortran rdf2u1 line 1355: call listio → ne=n1h, nls=n2h
    list_cont = read_cont(io)
    SPI = list_cont.C1
    AP = list_cont.C2
    LSSF = Int(list_cont.L1)
    NE = Int(list_cont.N1)    # number of fission width energy nodes
    NLS = Int(list_cont.N2)   # number of l-values

    # LIST body: NE fission-width energy nodes
    energies = _read_list_data(io, NE)

    # For each l-value: CONT + per-J LIST records
    sequences = URRSequence[]
    awri_val = 0.0
    for _ in 1:NLS
        lcont = read_cont(io)
        awri_val = lcont.C1   # AWRI
        ll = Int(lcont.L1)     # L value
        NJS = Int(lcont.N1)    # number of J values

        for _ in 1:NJS
            # LIST: C1=0, C2=0, L1=L, L2=MUF, N1=NPL, N2=0
            jcont = read_cont(io)
            MUF = Int(jcont.L2)
            NPL = Int(jcont.N1)
            jdata = _read_list_data(io, NPL)
            # jdata layout: D, AJ, AMUN, GN0, GG, 0, GF(1)..GF(NE)
            # For LFW=1: GX=0 (no competitive width), GF values start at index 7
            D    = jdata[1]
            AJ   = jdata[2]
            AMUN = jdata[3]
            GN0  = jdata[4]
            GG   = jdata[5]
            gf_vals = length(jdata) >= 7 ? jdata[7:min(end, 6 + NE)] : Float64[]
            push!(sequences, URRSequence(ll, AJ, D, AMUN, GN0, GG, 0.0, MUF, gf_vals))
        end
    end

    URRData(0.0, 0.0, SPI, AP, awri_val, LSSF, 1, NAPS, NRO, energies, sequences)
end

"""Read N floating-point values from LIST body (6 values per line)."""
function _read_list_data(io::IO, n::Int)
    vals = Float64[]
    while length(vals) < n
        line = readline(io)
        p = rpad(line, 80)
        for k in 0:5
            length(vals) >= n && break
            field = p[k*11+1 : k*11+11]
            push!(vals, parse_endf_float(field))
        end
    end
    vals
end
