# unresr_module — Bondarenko self-shielding in the unresolved resonance range
#
# Matches Fortran unresr.f90: reads ENDF MF2 unresolved parameters,
# reads MF3 backgrounds from input PENDF, computes Bondarenko XS at
# each (energy, temperature, sigma0) via bondarenko_xs(), writes output
# PENDF with MF2/MT152 section inserted after MT151.
#
# The output PENDF is a copy of the input with MT152 added.  MT152
# contains one CONT+LIST block per temperature with the self-shielded
# XS at (nunr energies) × (nsigz dilutions) × (5 reaction types).

# ── URR data for unresr ──────────────────────────────────────────────

struct URRSequenceData
    l::Int
    j::Float64
    d::Float64
    gn0::Float64
    gg::Float64
    gf_values::Vector{Float64}   # energy-dependent fission width (length nunr)
    gx::Float64
    amun::Int
    amuf::Int
    amux::Int
end

struct URRForUnresr
    za::Float64
    awr::Float64
    lssf::Int
    intunr::Int
    spi::Float64
    ap::Float64
    energies::Vector{Float64}           # final subdivided grid (nunr pts)
    endf_energies::Vector{Float64}      # original ENDF GF energy nodes
    sequences::Vector{URRSequenceData}
end

# ── Standard 78-point equi-lethargy reference grid (Fortran egridu) ──
const _UNRESR_EGRIDU = Float64[
    1.0e1, 1.25e1, 1.5e1, 1.7e1, 2.0e1, 2.5e1, 3.0e1,
    3.5e1, 4.0e1, 5.0e1, 6.0e1, 7.2e1, 8.5e1, 1.0e2,
    1.25e2, 1.5e2, 1.7e2, 2.0e2, 2.5e2, 3.0e2, 3.5e2,
    4.0e2, 5.0e2, 6.0e2, 7.2e2, 8.5e2, 1.0e3, 1.25e3,
    1.5e3, 1.7e3, 2.0e3, 2.5e3, 3.0e3, 3.5e3, 4.0e3,
    5.0e3, 6.0e3, 7.2e3, 8.5e3, 1.0e4, 1.25e4, 1.5e4,
    1.7e4, 2.0e4, 2.5e4, 3.0e4, 3.5e4, 4.0e4, 5.0e4,
    6.0e4, 7.2e4, 8.5e4, 1.0e5, 1.25e5, 1.5e5, 1.7e5,
    2.0e5, 2.5e5, 3.0e5, 3.5e5, 4.0e5, 5.0e5, 6.0e5,
    7.2e5, 8.5e5, 1.0e6, 1.25e6, 1.5e6, 1.7e6, 2.0e6,
    2.5e6, 3.0e6, 3.5e6, 4.0e6, 5.0e6, 6.0e6, 7.2e6, 8.5e6
]

"""
    _build_unresr_grid(raw_nodes, el, eh, elr, indep) -> Vector{Float64}

Build the unresr energy grid by subdividing between raw nodes using the
78-point reference grid. Matches Fortran rdunf2 lines 699-741.

- `raw_nodes`: energy nodes from ENDF (fission-width energies for mode 11)
- `el, eh`: unresolved range boundaries
- `elr`: resolved range upper boundary (for overlap flagging)
- `indep`: 0 for energy-dependent (always subdivide), 1 for energy-independent
"""
function _build_unresr_grid(raw_nodes::Vector{Float64}, el::Float64, eh::Float64,
                            elr::Float64, indep::Int)
    wide = 1.26; step_ratio = 1.01; onemev = 1.0e6

    # Build initial sorted node list with boundary shading
    # Fortran uses 1e6 sentinel at END (ilist does sorted insertion)
    nodes = Float64[onemev]  # sentinel (goes at end after sorted insertion)
    _ilist!(round_sigfig(el, 7, -1), nodes)
    _ilist!(round_sigfig(el, 7, +1), nodes)
    _ilist!(round_sigfig(eh, 7, -1), nodes)
    _ilist!(round_sigfig(eh, 7, +1), nodes)
    # Fortran rdunf2 line 591: skip first energy (i.gt.1), add rest
    for (idx, e) in enumerate(raw_nodes)
        idx == 1 && continue  # skip first ENDF energy node
        enow = round_sigfig(e, 7, 0)
        _ilist!(enow, nodes)
    end

    # Subdivide: insert reference grid points between consecutive nodes
    # Matching Fortran lines 699-726
    result = copy(nodes)
    i = 1
    elast = result[2]
    while true
        i += 1
        i + 1 > length(result) && break
        enext = result[i + 1]
        enext >= onemev && break

        if enext < wide * elast && indep == 0
            # Close enough — no subdivision needed
        else
            # Subdivide between elast and enext
            et = elast
            while true
                # Find next reference grid point > step_ratio * et
                enut = enext  # default if no reference point found
                for eg in _UNRESR_EGRIDU
                    if eg > step_ratio * et
                        enut = eg
                        break
                    end
                end
                et = enut
                et >= enext && break
                _ilist!(et, result)
                i += 1  # account for inserted point
            end
        end
        i + 1 > length(result) && break
        elast = result[i + 1]
    end

    # Cleanup: remove first sentinel, deduplicate via sigfig(7,2),
    # then remove last entry. Matches Fortran lines 727-741.
    cleaned = Float64[]
    en = 0.0
    lim = length(result) - 1  # skip last sentinel
    for k in 2:lim             # skip first sentinel
        et = result[k]
        if et >= en
            push!(cleaned, et)
            en = round_sigfig(abs(et), 7, 2)
        end
    end
    # Fortran line 741: nunr=nunr-1 (remove last entry = sigfig(eh,7,+1))
    if !isempty(cleaned)
        pop!(cleaned)
    end
    cleaned
end

"""Sorted insertion into ascending list (matches Fortran ilist)."""
function _ilist!(e::Float64, elist::Vector{Float64})
    pos = searchsortedfirst(elist, e)
    if pos <= length(elist) && elist[pos] == e
        return  # already present
    end
    insert!(elist, pos, e)
end

"""Interpolate fission width GF from ENDF energy nodes to unresr grid."""
function _interp_gf(endf_energies::Vector{Float64}, gf_endf::Vector{Float64},
                    grid_energies::Vector{Float64})
    ne = length(endf_energies)
    result = Float64[]
    for e in grid_energies
        if ne <= 1
            push!(result, isempty(gf_endf) ? 0.0 : gf_endf[1])
        elseif e <= endf_energies[1]
            push!(result, gf_endf[1])
        elseif e >= endf_energies[end]
            push!(result, gf_endf[end])
        else
            # Linear interpolation (INT=2)
            idx = searchsortedlast(endf_energies, e)
            idx = clamp(idx, 1, ne - 1)
            t = (e - endf_energies[idx]) / (endf_energies[idx+1] - endf_energies[idx])
            push!(result, gf_endf[idx] + t * (gf_endf[idx+1] - gf_endf[idx]))
        end
    end
    result
end

# ── MF2 reader for unresr ───────────────────────────────────────────

"""
    _read_urr_for_unresr(endf_path, mat) -> URRForUnresr

Read MF2/MT151 from the ENDF tape and extract unresolved resonance
parameters.  Handles mode 11 (LRF=1, LFW=1: energy-dependent fission
widths) and mode 10 (LRF=1, LFW=0: all energy-independent).
"""
function _read_urr_for_unresr(endf_path::String, mat::Int)
    io = open(endf_path)
    try
        found = find_section(io, 2, 151; target_mat=mat)
        found || error("MF2/MT151 not found for MAT=$mat in $endf_path")

        head = read_cont(io)
        za = head.C1; awr = head.C2; nis = Int(head.N1)

        urr_data = nothing
        elr = 0.0  # resolved range upper energy (for overlap)
        urr_el = 0.0; urr_eh = 0.0
        lfw_global = 0

        for _ in 1:nis
            iso_cont = read_cont(io)
            lfw_global = Int(iso_cont.L2)
            ner = Int(iso_cont.N1)

            for _ in 1:ner
                rng = read_cont(io)
                el = rng.C1; eh = rng.C2
                lru = Int(rng.L1); lrf = Int(rng.L2)
                nro = Int(rng.N1); naps = Int(rng.N2)

                if lru == 2
                    urr_el = el; urr_eh = eh
                    if nro > 0
                        _unresr_skip_tab1(io)
                    end
                    urr_data = _read_urr_subsection(io, za, awr, lfw_global, lrf)
                elseif lru == 1
                    elr = eh  # track resolved range boundary
                    _skip_resolved_range(io, lrf, nro)
                end
            end
        end

        urr_data === nothing && error("No unresolved range (LRU=2) found for MAT=$mat")

        # Build the subdivided energy grid (Fortran rdunf2 lines 699-741)
        indep = (lfw_global == 0) ? 1 : 0  # energy-dep → indep=0, energy-indep → indep=1
        endf_energies = urr_data.endf_energies
        grid = _build_unresr_grid(endf_energies, urr_el, urr_eh, elr, indep)

        # Interpolate fission widths from ENDF nodes to subdivided grid
        new_seqs = URRSequenceData[]
        for s in urr_data.sequences
            gf_interp = _interp_gf(endf_energies, s.gf_values, grid)
            push!(new_seqs, URRSequenceData(
                s.l, s.j, s.d, s.gn0, s.gg, gf_interp, s.gx,
                s.amun, s.amuf, s.amux))
        end

        return URRForUnresr(za, awr, urr_data.lssf, urr_data.intunr,
                           urr_data.spi, urr_data.ap, grid, endf_energies, new_seqs)
    finally
        close(io)
    end
end

"""Read one unresolved subsection (after the range CONT).
Matches the existing reconr reader pattern in reader.jl."""
function _read_urr_subsection(io::IO, za::Float64, awr::Float64,
                              lfw::Int, lrf::Int)
    intunr = 2  # default interpolation scheme
    sequences = URRSequenceData[]
    energies = Float64[]
    spi = 0.0; ap = 0.0; lssf = 0

    if lrf <= 1 && lfw == 1
        # Mode 11 (LFW=1, LRF≤1): energy-dependent fission widths
        # Format matches _read_urr_lfw1 in reader.jl:
        # LIST header: C1=SPI, C2=AP, L1=LSSF, L2=0, N1=NE, N2=NLS
        # LIST body: NE energy grid values
        # Then per-l CONT + per-J LIST records
        list_cont = read_cont(io)
        spi = list_cont.C1; ap = list_cont.C2; lssf = Int(list_cont.L1)
        ne = Int(list_cont.N1)   # number of fission width energies
        nls = Int(list_cont.N2)  # number of l-values

        # Read NE energy values (LIST body)
        energies = _read_list_data(io, ne)

        for _ in 1:nls
            lcont = read_cont(io)
            l_val = Int(lcont.L1)
            njs = Int(lcont.N1)

            for _ in 1:njs
                # LIST header: C1=0, C2=0, L1=L, L2=MUF, N1=NPL, N2=0
                jcont = read_cont(io)
                muf = Int(jcont.L2)
                npl = Int(jcont.N1)
                jdata = _read_list_data(io, npl)
                # jdata: D, AJ, AMUN, GN0, GG, 0, GF(1)..GF(NE)
                d_val = jdata[1]; aj = jdata[2]
                amun = round(Int, jdata[3])
                gn0 = jdata[4]; gg = jdata[5]
                gf_vals = length(jdata) >= 7 ? jdata[7:min(end, 6 + ne)] : Float64[]

                push!(sequences, URRSequenceData(
                    l_val, aj, d_val, gn0, gg, gf_vals, 0.0,
                    amun, muf, 0))
            end
        end

    elseif lrf <= 1 && lfw == 0
        # Mode 10 (LFW=0, LRF≤1): all energy-independent
        # CONT: SPI, AP, LSSF, 0, NLS, 0
        sub = read_cont(io)
        spi = sub.C1; ap = sub.C2; lssf = Int(sub.L1)
        nls = Int(sub.N1)

        for _ in 1:nls
            lcont = read_cont(io)
            l_val = Int(lcont.L1)
            njs = Int(lcont.N2)

            for _ in 1:njs
                lst = read_list(io)
                d_val = lst.data[1]; aj = lst.data[2]
                amun = round(Int, lst.data[3])
                gn0 = lst.data[4]; gg = lst.data[5]
                gf0 = lst.data[6]
                ne = Int(lst.N2)

                e_grid = Float64[]
                for ie in 1:ne
                    push!(e_grid, lst.data[6 + 6*(ie-1) + 1])
                end
                if isempty(energies); energies = e_grid; end

                push!(sequences, URRSequenceData(
                    l_val, aj, d_val, gn0, gg, fill(gf0, max(1, length(energies))),
                    0.0, amun, 0, 0))
            end
        end

    elseif lrf == 2
        # Mode 12 (LRF=2): all parameters energy-dependent
        # CONT: SPI, AP, LSSF, 0, NLS, 0
        sub = read_cont(io)
        spi = sub.C1; ap = sub.C2; lssf = Int(sub.L1)
        nls = Int(sub.N1)

        for _ in 1:nls
            lcont = read_cont(io)
            l_val = Int(lcont.L1); njs = Int(lcont.N1)

            for _ in 1:njs
                # LIST: AJ, 0, INT, 0, NPL, NE
                jcont = read_cont(io)
                aj = jcont.C1
                intunr_local = Int(jcont.L1)
                if intunr_local > 0; intunr = intunr_local; end
                npl = Int(jcont.N1); ne = Int(jcont.N2)
                jdata = _read_list_data(io, npl)
                # [0, 0, AMUX, AMUN, AMUG, AMUF, E1,D1,GX1,GN01,GG1,GF1, ...]
                amux = round(Int, jdata[3])
                amun = round(Int, jdata[4])
                amuf = round(Int, jdata[6])

                e_grid = Float64[]; gf_vals = Float64[]
                for ie in 1:ne
                    base = 6*(ie - 1)
                    push!(e_grid, jdata[7 + base])
                    push!(gf_vals, jdata[12 + base])
                end
                if isempty(energies); energies = e_grid; end

                push!(sequences, URRSequenceData(
                    l_val, aj, jdata[8], jdata[10], jdata[11],
                    gf_vals, jdata[9], amun, amuf, amux))
            end
        end
    else
        error("Unsupported URR format: LRF=$lrf, LFW=$lfw")
    end

    # Return with endf_energies = energies (not yet subdivided — caller builds final grid)
    URRForUnresr(za, awr, lssf, intunr, spi, ap, energies, energies, sequences)
end

"""Skip a resolved resonance range. Uses line-scanning to be safe."""
function _skip_resolved_range(io::IO, lrf::Int, nro::Int)
    # Safest approach: scan lines until we hit the end of this range's data.
    # We know the next record after this range is either another range CONT
    # or the isotope SEND record. We can't easily detect these without
    # knowing the format, so we use the existing _skip_unsupported_range
    # approach: read lines while they belong to MF2/MT151.
    # Actually, just read CONT+LIST records per the format.

    if nro > 0
        _unresr_skip_tab1(io)
    end

    sub = read_cont(io)  # SPI, AP, 0, 0, NLS, 0
    nls = Int(sub.N1)

    if lrf <= 2
        # SLBW/MLBW: each l-value is ONE list record (CONT header + resonance data)
        for _ in 1:nls
            lcont = read_cont(io)  # AWRI, QX, L, LRX, 6*NRS, NRS
            nrs = Int(lcont.N2)
            if nrs > 0
                _read_list_data(io, 6 * nrs)  # skip the resonance data
            end
        end
    elseif lrf == 3
        # Reich-Moore: similar to BW
        for _ in 1:nls
            lcont = read_cont(io)  # APL, 0, L, 0, 6*NRS, NRS
            nrs = Int(lcont.N2)
            if nrs > 0
                _read_list_data(io, 6 * nrs)
            end
        end
    else
        # SAMMY/other: use line-scanning fallback
        _unresr_skip_to_end_of_section(io)
    end
end

"""Skip lines until we leave MF2/MT151 (fallback for complex formats)."""
function _unresr_skip_to_end_of_section(io::IO)
    while !eof(io)
        pos = position(io)
        line = readline(io)
        p = rpad(line, 80)
        mf = _parse_int(p[71:72]); mt = _parse_int(p[73:75])
        if mf != 2 || mt != 151
            seek(io, pos)
            return
        end
    end
end

function _unresr_skip_tab1(io::IO)
    read_tab1(io)  # just read and discard
end

function _unresr_skip_list(io::IO)
    read_list(io)  # just read and discard
end

function _unresr_skip_list_or_cont(io::IO)
    # Read a record — could be LIST or CONT
    line = readline(io)
    p = rpad(line, 80)
    n1 = _parse_int(p[45:55])
    n2 = _parse_int(p[56:66])
    # If N1 > 0 and it looks like a LIST, skip data lines
    npl = n1
    if npl > 0
        nlines = div(npl + 5, 6)
        for _ in 1:nlines
            readline(io)
        end
    end
end

# ── MF3 background reader ───────────────────────────────────────────

"""
Read MF3 backgrounds (MT=1,2,18,102) from the ENDF tape at URR energies.
Returns matrix (nunr, 4): total, elastic, fission, capture.

Matches Fortran rdunf3 which reads from nendf (ENDF), not the PENDF.
Fortran nudges first and last energies by ×1.00001 and ×0.99999.
"""
function _read_unresr_backgrounds(endf_path::String, mat::Int,
                                  energies::Vector{Float64})
    nunr = length(energies)
    bkg = zeros(nunr, 4)  # total, elastic, fission, capture

    # Fortran mtx = [1, 2, 18, 19, 102] → columns [1, 2, 3, skip, 4]
    mt_map = [(1, 1), (2, 2), (18, 3), (102, 4)]

    for (mt, col) in mt_map
        io = open(endf_path)
        try
            found = find_section(io, 3, mt; target_mat=mat)
            found || continue
            read_cont(io)  # skip HEAD record; TAB1 follows
            tab = read_tab1(io)
            length(tab.x) < 1 && continue
            tf = TabulatedFunction(tab)
            for (ie, e) in enumerate(energies)
                # Fortran rdunf3 line 823-824: nudge first/last energies
                eq = abs(e)
                if ie == 1; eq *= 1.00001; end
                if ie == nunr; eq *= 0.99999; end
                if eq >= tf.x[1] && eq <= tf.x[end]
                    bkg[ie, col] = interpolate(tf, eq)
                end
            end
        finally
            close(io)
        end
    end

    bkg
end

# ── Bondarenko computation ───────────────────────────────────────────

"""
Compute Bondarenko XS at all URR energies for one temperature.
Returns matrix (nunr, 5*nsigz): for each energy, 5 reaction types × nsigz.
Reaction order matches Fortran: total, elastic, capture, fission, transport.
"""
function _bondarenko_all_energies(urr::URRForUnresr, temp::Float64,
                                  sigz::Vector{Float64},
                                  bkg::Matrix{Float64})
    nunr = length(urr.energies)
    nsigz = length(sigz)
    # Output: (nunr, 5*nsigz) — 5 reaction types interleaved with nsigz values each
    result = zeros(nunr, 5 * nsigz)

    for ie in 1:nunr
        # Build URRStatModel at this energy (fission widths may be energy-dependent)
        seqs = URRSpinSequence[]
        for s in urr.sequences
            gf = ie <= length(s.gf_values) ? s.gf_values[ie] : s.gf_values[end]
            push!(seqs, URRSpinSequence(
                s.l, s.j, s.d, s.gn0, s.gg, gf, s.gx,
                s.amun, s.amuf, s.amux))
        end
        model = URRStatModel(urr.spi, urr.awr, urr.ap, seqs)

        # Background XS: (total, elastic, fission, capture)
        bkg_tuple = (bkg[ie, 1], bkg[ie, 2], bkg[ie, 3], bkg[ie, 4])

        # Compute Bondarenko XS
        sigu = bondarenko_xs(model, urr.energies[ie], temp, sigz; bkg=bkg_tuple)
        # sigu rows match Fortran: 1=total, 2=elastic, 3=fission, 4=capture, 5=transport
        # No swap needed — both use the same order
        for ix in 1:5
            for is in 1:nsigz
                result[ie, (ix-1)*nsigz + is] = round_sigfig(sigu[ix, is], 7, 0)
            end
        end
    end

    result
end

# ── PENDF writer ─────────────────────────────────────────────────────

"""
Write output PENDF: rewrite MF1 header with MT152 directory entry,
copy MF2/MT151 from input, insert MT152 data, copy MF3+ verbatim.

Matches Fortran unresr output structure (lines 241-321):
- TPID with reconr description
- MF1/MT451 with NWD description lines + directory including MT152
- MF2/MT151 copied from input
- MF2/MT152 with Bondarenko tables
- MF3+ copied verbatim
"""
function _write_unresr_pendf(pendf_in::String, pendf_out::String,
                             params::UnresrParams, urr::URRForUnresr,
                             mt152_blocks::Vector)
    in_lines = readlines(pendf_in)
    mat = params.mat

    # Parse input PENDF structure
    tpid_line = in_lines[1]
    mf1_start = 0; mf1_end = 0
    mf2_start = 0; mf2_fend = 0
    mf3_start = 0
    old_mt152_start = 0; old_mt152_end = 0

    for (i, line) in enumerate(in_lines)
        length(line) < 14 && continue
        # Parse MAT/MF/MT from last 14 chars (handles variable-width lines)
        n = length(line)
        pmat = _parse_int(line[n-13:n-10]); pmf = _parse_int(line[n-9:n-8]); pmt = _parse_int(line[n-7:n-5])
        pmat == mat || continue
        if pmf == 1 && pmt == 451 && mf1_start == 0; mf1_start = i; end
        if pmf == 1 && pmt == 0 && mf1_start > 0 && mf1_end == 0; mf1_end = i; end
        if pmf == 2 && pmt == 151 && mf2_start == 0; mf2_start = i; end
        if pmf == 2 && pmt == 152; old_mt152_start == 0 && (old_mt152_start = i); old_mt152_end = i; end
        if pmf == 0 && pmt == 0 && mf2_start > 0 && mf2_fend == 0; mf2_fend = i; end
        if pmf == 3 && mf3_start == 0; mf3_start = i; end
    end

    # Extract existing descriptions from input MF1
    descriptions = String[]
    za = urr.za; awr = urr.awr
    old_nwd = 0; temp_val = 0.0; tol_val = 0.0
    if mf1_start > 0
        p2 = rpad(in_lines[mf1_start], 80)[1:80]
        za = parse_endf_float(p2[1:11]); awr = parse_endf_float(p2[12:22])
        if mf1_start + 1 <= length(in_lines)
            p3 = rpad(in_lines[mf1_start + 1], 80)[1:80]
            temp_val = parse_endf_float(p3[1:11])
            tol_val = parse_endf_float(p3[12:22])
            old_nwd = Int(_parse_int(p3[45:55]))
        end
        for k in 1:old_nwd
            idx = mf1_start + 1 + k
            idx <= length(in_lines) && push!(descriptions, rpad(in_lines[idx], 80)[1:min(66, length(rpad(in_lines[idx], 80)))])
        end
    end

    # Collect directory entries from input (MF2 + MF3 sections)
    dir_entries = Tuple{Int,Int,Int,Int}[]  # (mf, mt, nc, mod)
    dir_start = mf1_start + 2 + old_nwd  # first directory line in input
    if dir_start <= mf1_end
        for i in dir_start:mf1_end-1  # exclude SEND
            ln = in_lines[i]; nl = length(ln); nl < 14 && continue
            _parse_int(ln[nl-13:nl-10]) == mat || continue
            _parse_int(ln[nl-9:nl-8]) == 1 && _parse_int(ln[nl-7:nl-5]) == 451 || continue
            p = rpad(ln, 80)[1:80]
            # Parse directory: try standard 6×11 ENDF fields then broadr %32s+%11d format
            # Broadr writes %32s%11d%11d%11d%11d → MF at [34:44], MT at [45:55], NC at [56:66]
            # But sometimes shifted: MF at [33:43], etc. Use field parsing for robustness.
            vals = Int[]
            for f in 1:6
                v = Int(_parse_int(p[(f-1)*11+1:f*11]))
                v != 0 && push!(vals, v)
            end
            length(vals) < 3 && continue
            dmf = vals[1]; dmt = vals[2]; dnc = vals[3]
            (dmf == 0 && dmt == 0) && continue
            (dmf == 1 && dmt == 451) && continue  # skip self-reference
            push!(dir_entries, (dmf, dmt, dnc, 0))
        end
    end

    # Add MT152 directory entry if not already present.
    # If MT152 already exists (from broadr), keep its NC value — the Fortran
    # unresr does NOT update the directory NC when new=0 (line 283: j only
    # advances if new>0, so the computed ncds gets overwritten by the old entry).
    mt152_nc = 2 + div(length(params.sigz) + length(urr.energies) * (1 + 5 * length(params.sigz)) - 1, 6)
    has_mt152 = any(d -> d[1] == 2 && d[2] == 152, dir_entries)
    if !has_mt152
        idx152 = findfirst(d -> d[1] == 2 && d[2] == 151, dir_entries)
        if idx152 !== nothing
            insert!(dir_entries, idx152 + 1, (2, 152, mt152_nc, 0))
        else
            push!(dir_entries, (2, 152, mt152_nc, 0))
        end
    else
        # Keep existing NC from broadr (Fortran behavior: new=0 → no update)
        # Don't modify dir_entries — broadr NC is preserved
    end

    # Update temperature to first temperature (Fortran writes MF1 at it=1)
    temp_val = isempty(params.temperatures) ? temp_val : params.temperatures[1]

    nwd = length(descriptions)
    nxc = length(dir_entries) + 1  # +1 for self-reference

    # Strategy: copy input line-by-line, replacing the first block's
    # MF1/MT451 and MT152 while copying everything else verbatim.
    # This handles multi-temperature PENDFs correctly.

    # Find first block's MT152 boundaries (for replacement) and MT151 SEND
    blk1_mt152_start = 0; blk1_mt152_send = 0; blk1_mt151_send = 0
    for (i, line) in enumerate(in_lines)
        length(line) < 14 && continue
        nl = length(line)
        lmat = _parse_int(line[nl-13:nl-10])
        lmf = _parse_int(line[nl-9:nl-8])
        lmt = _parse_int(line[nl-7:nl-5])
        lmat == mat || continue
        if lmf == 2 && lmt == 152 && blk1_mt152_start == 0
            blk1_mt152_start = i
        end
        if lmf == 2 && lmt == 0 && blk1_mt152_start > 0 && blk1_mt152_send == 0
            blk1_mt152_send = i
        end
        # Find first MF2 SEND (end of MT151) — for inserting MT152 when none exists
        if lmf == 2 && lmt == 0 && blk1_mt151_send == 0
            blk1_mt151_send = i
        end
    end

    open(pendf_out, "w") do out
        i = 1
        while i <= length(in_lines)
            # Replace first MF1/MT451 section
            if i == mf1_start
                # Write TPID (line before MF1)
                # (TPID is already written at line 1 = tpid_line)

                # Write rewritten MF1/MT451
                ns = 1
                @printf(out, "%s%s%11d%11d%11d%11d%4d%2d%3d%5d\n",
                    format_endf_float(za), format_endf_float(awr),
                    3, 1, 0, nxc, mat, 1, 451, ns); ns += 1
                @printf(out, "%s%s%11d%11d%11d%11d%4d%2d%3d%5d\n",
                    format_endf_float(temp_val), format_endf_float(tol_val),
                    0, 0, nwd, nxc, mat, 1, 451, ns); ns += 1
                for desc in descriptions
                    @printf(out, "%-66s%4d%2d%3d%5d\n", desc, mat, 1, 451, ns); ns += 1
                end
                self_nc = 2 + nwd + nxc
                @printf(out, "%22s%11d%11d%11d%11d%4d%2d%3d%5d\n",
                    "", 1, 451, self_nc, 0, mat, 1, 451, ns); ns += 1
                for (dmf, dmt, dnc, dmod) in dir_entries
                    @printf(out, "%22s%11d%11d%11d%11d%4d%2d%3d%5d\n",
                        "", dmf, dmt, dnc, dmod, mat, 1, 451, ns); ns += 1
                end
                @printf(out, "%66s%4d%2d%3d%5d\n", "", mat, 1, 0, 99999)
                @printf(out, "%66s%4d%2d%3d%5d\n", "", mat, 0, 0, 0)  # FEND
                # Skip to after MF1 SEND + FEND + zero/gap lines from input
                i = mf1_end + 1
                # Skip FEND and gap lines between MF1 and MF2
                while i <= length(in_lines)
                    ln = in_lines[i]; nl = length(ln)
                    nl < 14 && (i += 1; continue)
                    lmf = _parse_int(ln[nl-9:nl-8]); lmt = _parse_int(ln[nl-7:nl-5])
                    (lmf == 0 && lmt == 0) && (i += 1; continue)  # skip FEND/gap
                    break
                end
                continue

            # Replace first MT152 block with full unresr data
            elseif blk1_mt152_start > 0 && i == blk1_mt152_start
                _write_mt152(out, params, urr, mt152_blocks)
                # Skip old MT152 data; the SEND at blk1_mt152_send terminates
                i = blk1_mt152_send
                continue

            # Insert MT152 after MT151 SEND when no existing MT152
            elseif blk1_mt152_start == 0 && blk1_mt151_send > 0 && i == blk1_mt151_send
                # Write the MT151 SEND line first
                println(out, in_lines[i])
                # Insert new MT152 + SEND (MF2 FEND follows from broadr copy)
                _write_mt152(out, params, urr, mt152_blocks)
                @printf(out, "%66s%4d%2d%3d%5d\n", "", mat, 2, 0, 99999)
                i += 1
                continue

            else
                println(out, in_lines[i])
                i += 1
            end
        end
    end
end

"""Format and write MT152 records to output stream."""
function _write_mt152(io::IO, params::UnresrParams, urr::URRForUnresr,
                      mt152_blocks::Vector)
    nunr = length(urr.energies)
    nsigz = length(params.sigz)
    mat = params.mat
    trailer(seq) = @sprintf("%4d%2d%3d%5d", mat, 2, 152, seq)

    for (bi, block) in enumerate(mt152_blocks)
        temp = block.temp
        xs = block.xs  # (nunr, 5*nsigz)

        # CONT record: ZA, AWR, LSSF, 0, 0, INTUNR
        seq = 1
        cont = format_endf_float(urr.za) * format_endf_float(urr.awr) *
               lpad(string(urr.lssf), 11) * lpad("0", 11) *
               lpad("0", 11) * lpad(string(urr.intunr), 11)
        println(io, cont * trailer(seq))

        # LIST record header: temp, 0, NREACT=5, NSIGZ, NPL, NE=nunr
        npl = nsigz + nunr * (1 + 5 * nsigz)
        seq += 1
        list_head = format_endf_float(temp) * format_endf_float(0.0) *
                    lpad("5", 11) * lpad(string(nsigz), 11) *
                    lpad(string(npl), 11) * lpad(string(nunr), 11)
        println(io, list_head * trailer(seq))

        # LIST data: pack into 6-values-per-line
        data = Float64[]

        # First: sigma0 values
        append!(data, params.sigz)

        # Then: for each energy, energy value + 5×nsigz XS values
        for ie in 1:nunr
            push!(data, urr.energies[ie])
            for k in 1:5*nsigz
                push!(data, xs[ie, k])
            end
        end

        # Write data in 6-values-per-line format
        idx = 1
        while idx <= length(data)
            buf = ""
            for col in 1:6
                if idx <= length(data)
                    buf *= format_endf_float(data[idx])
                    idx += 1
                end
            end
            seq += 1
            println(io, rpad(buf, 66) * trailer(seq))
        end
    end

    # SEND record
    println(io, rpad("", 66) * @sprintf("%4d%2d%3d%5d", mat, 2, 0, 99999))
end

# ── Main entry point ────────────────────────────────────────────────

function unresr_module(tapes::TapeManager, params::UnresrParams)
    endf_path = resolve(tapes, params.nendf)
    pendf_in_path = resolve(tapes, params.npendf_in)
    pendf_out_path = joinpath(tapes.work_dir, "tape$(params.npendf_out)")

    # Read URR data from ENDF MF2/MT151
    urr = _read_urr_for_unresr(endf_path, params.mat)

    # Read MF3 backgrounds from input PENDF
    bkg = _read_unresr_backgrounds(endf_path, params.mat, urr.energies)

    # Compute Bondarenko XS for each temperature
    mt152_blocks = []
    for temp in params.temperatures
        xs = _bondarenko_all_energies(urr, temp, params.sigz, bkg)
        push!(mt152_blocks, (temp=temp, xs=xs))
    end

    # Write output PENDF (copy input + insert MT152)
    _write_unresr_pendf(pendf_in_path, pendf_out_path, params, urr, mt152_blocks)

    register!(tapes, params.npendf_out, pendf_out_path)
    @info "unresr: MAT=$(params.mat), $(length(params.temperatures)) temps, $(length(params.sigz)) sigma0, $(length(urr.energies)) energies"
    return nothing
end
