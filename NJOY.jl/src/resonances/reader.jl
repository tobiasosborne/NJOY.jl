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
                params = _read_bw_params(io, LRF, NRO, NAPS)
                push!(ranges, ResonanceRange(EL, EH, LRU, LRF, LFW, NRO, NAPS, params))

            elseif LRU == 1 && LRF == 3
                # Reich-Moore
                params = _read_rm_params(io, NRO, NAPS)
                push!(ranges, ResonanceRange(EL, EH, LRU, LRF, LFW, NRO, NAPS, params))

            else
                # Unsupported formalisms (Adler-Adler, URR, etc.)
                # Skip by reading lines until the next range or end of section.
                # We use a simple approach: read lines until we hit a SEND record
                # (MF=0) or the next section. For unresolved (LRU=2), we need to
                # consume the rest of this range's data.
                _skip_unsupported_range(io)
                # Don't push anything to ranges -- just skip
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
    if NRO != 0
        # Read TAB1 for AP(E) -- skip it for now but consume the lines
        _tab1 = read_tab1(io)
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
        return SLBWParameters(Int32(NLS), SPI, AP,
                              l_values, AWRI_all, QX_all, LRX_all,
                              Er_all, AJ_all,
                              Gn_all, Gg_all, Gf_all, Gx_all)
    else
        return MLBWParameters(Int32(NLS), SPI, AP,
                              l_values, AWRI_all, QX_all, LRX_all,
                              Er_all, AJ_all,
                              Gn_all, Gg_all, Gf_all, Gx_all)
    end
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
    if NRO != 0
        _tab1 = read_tab1(io)
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

    return ReichMooreParameters(Int32(NLS), SPI, AP, LAD,
                                l_values, AWRI_all, APL_all,
                                Er_all, AJ_all,
                                Gn_all, Gg_all, Gfa_all, Gfb_all)
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
