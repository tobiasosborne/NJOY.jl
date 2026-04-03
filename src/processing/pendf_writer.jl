# PENDF writer -- streaming output of pointwise cross sections
#
# PROPOSAL B DESIGN: The writer is designed for streaming I/O:
#   - Each section is written independently, never buffering the whole file
#   - The `write_pendf` function works with `PointwiseMaterial` directly
#   - MF2 is copied verbatim from the original ENDF source if available
#   - Provides both PointwiseMaterial and NamedTuple interfaces
#
# Correspondence to NJOY2016 reconr.f90:
#   recout (4984-5441) -> write_pendf, write_pendf_file

using Printf

# ==========================================================================
# Main entry points
# ==========================================================================

"""
    write_pendf(io::IO, material::PointwiseMaterial;
                endf_source=nothing, temperature=0.0, err=0.001,
                description=String[])

Write a PENDF tape for `material` to the output stream `io`.

The output is in standard ENDF-6 format with:
- MF1/MT451: descriptive data with dictionary
- MF2/MT151: resonance parameters (copied from source or minimal placeholder)
- MF3/MT{N}: pointwise cross sections with lin-lin interpolation
- Proper SEND/FEND/MEND/TEND delimiters
"""
function write_pendf(io::IO, material::PointwiseMaterial;
                     endf_source::Union{IO, Nothing} = nothing,
                     temperature::Float64 = 0.0,
                     err::Float64 = 0.001,
                     description::Vector{String} = String[],
                     za::Float64 = 0.0,
                     awr::Float64 = 0.0)
    mat = Int(material.mat)
    n_pts = length(material.energies)
    ns = Ref(1)  # line sequence counter

    # TPID record
    _write_tpid_line(io, "PENDF tape produced by NJOY.jl", mat)

    # MF1/MT451
    _write_mf1_section(io, mat, material, temperature, err, description, ns)

    # MF2/MT151
    if endf_source !== nothing
        _copy_mf2_from_source(io, endf_source, mat)
    else
        _write_minimal_mf2(io, mat, ns)
    end

    # MF3 for each reaction -- continuous sequence numbering across sections
    ns[] = 1
    for (col, mt) in enumerate(material.mt_list)
        _write_mf3_from_matrix(io, mat, mt, material.energies,
                               material.cross_sections, col, ns;
                               za=za, awr=awr)
    end
    _write_fend(io, mat; ns=ns)

    # MEND + TEND
    _write_record_sep(io, 0, 0, 0)
    _write_record_sep(io, -1, 0, 0)
end

"""
    write_pendf(io::IO, result::NamedTuple; mat=0, label="",
                err=0.001, tempr=0.0)

Write PENDF from the legacy reconr() NamedTuple result format.
"""
function write_pendf(io::IO, result::NamedTuple;
                     mat::Integer = 0,
                     label::AbstractString = "reconstructed data",
                     err::Float64 = 0.001,
                     tempr::Float64 = 0.0,
                     title = nothing,
                     descriptions::Vector{String} = String[])
    mf2 = result.mf2
    actual_mat = mat > 0 ? Int32(mat) : Int32(max(1, round(Int, mf2.ZA / 10)))
    ns = Ref(1)

    # Determine reactions
    reactions = _collect_reactions(result)

    # TPID: use title if provided, otherwise label
    tpid_text = title !== nothing ? title : label
    _write_tpid_line(io, tpid_text, Int(actual_mat))

    # Build MF lookup from sections
    mt_mf = Dict{Int,Int}()
    for sec in result.mf3_sections
        mt_mf[Int(sec.mt)] = Int(sec.mf)
    end

    # MF1/MT451
    _write_legacy_mf1(io, mf2, actual_mat, err, tempr, length(result.energies), reactions, ns;
                      mt_to_mf=mt_mf, descriptions=descriptions)

    # MF2/MT151
    _write_legacy_mf2(io, mf2, actual_mat, ns)

    # MF3/MF23 sections
    _write_legacy_mf3(io, result, actual_mat, reactions, ns)

    # MEND and TEND (NS=0 matching Fortran)
    _write_fend_zero(io, 0)   # MEND: MAT=0, MF=0, MT=0, NS=0
    blanks = repeat(" ", 66)
    @printf(io, "%s%4d%2d%3d%5d\n", blanks, -1, 0, 0, 0)  # TEND
end

"""
    write_pendf_file(filename, data; kwargs...)

Write PENDF to a file. `data` can be `PointwiseMaterial` or `NamedTuple`.
"""
function write_pendf_file(filename::AbstractString, data; kwargs...)
    open(filename, "w") do io
        write_pendf(io, data; kwargs...)
    end
end

# ==========================================================================
# Full-pipeline PENDF writer: per-MT grids + MF12/MF13/MF6
# ==========================================================================

"""
    write_full_pendf(io, reconr_result; mat=0, label="", err=0.001, tempr=0.0,
                     override_mf3=Dict(), extra_mf3=Dict(),
                     mf6_records=Dict(), mf12_lines=String[], mf13_lines=String[])

Write a complete PENDF file from reconr output + per-module augmentations.
Each MF3 section has its OWN energy grid (not shared).

- `override_mf3`: Dict{Int, Tuple{Vector,Vector}} — MT → (energies, xs) replacing reconr data
- `extra_mf3`: Dict{Int, Tuple{Vector,Vector}} — MT → (energies, xs) for new sections (301,444,221,etc.)
- `mf6_records`: Dict{Int, Vector{MF6ListRecord}} — MT → angular distribution data
- `mf12_lines`: raw ENDF lines for MF12 passthrough (from reconr oracle or ENDF)
- `mf13_lines`: raw ENDF lines for MF13 passthrough
"""
function write_full_pendf(io::IO, result::NamedTuple;
                          mat::Integer = 0,
                          label::AbstractString = "reconstructed data",
                          err::Float64 = 0.001,
                          tempr::Float64 = 0.0,
                          override_mf3::Dict{Int, <:Tuple} = Dict{Int, Tuple{Vector{Float64},Vector{Float64}}}(),
                          extra_mf3::Dict{Int, <:Tuple} = Dict{Int, Tuple{Vector{Float64},Vector{Float64}}}(),
                          mf6_records::Dict{Int, <:Any} = Dict{Int, Any}(),
                          mf6_stubs::Dict{Int, <:NamedTuple} = Dict{Int, NamedTuple}(),
                          mf12_lines::Vector{String} = String[],
                          mf13_lines::Vector{String} = String[],
                          descriptions::Vector{String} = String[],
                          mf6_xsi::Dict{Int, Vector{Float64}} = Dict{Int, Vector{Float64}}(),
                          mf6_emax::Dict{Int, Float64} = Dict{Int, Float64}(),
                          thermr_mts::Set{Int} = Set{Int}(),
                          thermr_coh_ne::Int = -1)
    mf2 = result.mf2
    actual_mat = mat > 0 ? Int32(mat) : Int32(max(1, round(Int, mf2.ZA / 10)))
    ns = Ref(1)

    # Collect all MF3 reactions (reconr + extra)
    reactions = _collect_reactions(result)
    extra_mts = sort(collect(keys(extra_mf3)))

    # Insert extra MTs at their natural position
    for emt in extra_mts
        any(r -> r[1] == Int32(emt), reactions) && continue
        idx = findfirst(r -> r[1] > Int32(emt), reactions)
        entry = (Int32(emt), "MT$emt")
        idx === nothing ? push!(reactions, entry) : insert!(reactions, idx, entry)
    end

    # Count lines for each section (needed for MF1 directory)
    section_lines = Dict{Tuple{Int,Int}, Int}()  # (MF,MT) → line count

    # Pre-compute MF3 line counts
    # Fortran tpend NC formula: nc = 3 + (ne_param+2)/3 where ne_param is the
    # grid count BEFORE tpend adds sentinel points. For free-gas thermr,
    # ne_param = NP. For coh-grid thermr, ne_param = NP - 2 (coh sentinels).
    for (mt, _) in reactions
        imt = Int(mt)
        if haskey(override_mf3, imt)
            e, xs = override_mf3[imt]
            np = length(e)
        elseif haskey(extra_mf3, imt)
            e, xs = extra_mf3[imt]
            np = length(e)
        else
            sec = _get_legacy_section(result, imt)
            sec === nothing && continue
            np = length(sec[1])
        end
        if imt in thermr_mts && thermr_coh_ne >= 0
            # Use Fortran tpend formula: nc = 3 + (ne+2)/3
            # For coh-grid sections (NP > coh_ne), ne = coh_ne
            # For free-gas sections (NP < coh_ne), ne = NP
            ne_param = np > thermr_coh_ne ? thermr_coh_ne : np
            section_lines[(3, imt)] = 3 + div(ne_param + 2, 3)
        else
            section_lines[(3, imt)] = 3 + cld(np, 3)
        end
    end

    # MF6 stub line counts: Fortran ncdse=3 (undercounts TAB1 by 1 vs actual 4)
    for (mt, _) in mf6_stubs
        section_lines[(6, mt)] = 3
    end

    # MF6 line counts: Fortran ncds undercounts TAB1 yield by 1
    # (counts packed interp+data as 1 line, but output has 2 separate lines)
    for (mt, records) in mf6_records
        ne = length(records)
        # HEAD(1) + TAB1(3 lines) + TAB2(2 lines)
        n_lines = 1 + 3 + 2
        for rec in records
            n_ep = length(rec.entries)
            nw = n_ep * 10
            n_lines += 1 + cld(nw, 6)  # LIST header + data lines
        end
        # Subtract 1: Fortran ncds counts TAB1 as 2 lines (packed), actual is 3
        section_lines[(6, mt)] = n_lines - 1
    end

    # MF12/MF13 line counts
    if !isempty(mf12_lines)
        # Parse MT from the lines
        for line in mf12_lines
            length(line) >= 75 || continue
            p = rpad(line, 80)
            mt = _parse_int(p[73:75])
            mt > 0 && (section_lines[(12, mt)] = get(section_lines, (12, mt), 0) + 1)
        end
    end
    if !isempty(mf13_lines)
        for line in mf13_lines
            length(line) >= 75 || continue
            p = rpad(line, 80)
            mt = _parse_int(p[73:75])
            mt > 0 && (section_lines[(13, mt)] = get(section_lines, (13, mt), 0) + 1)
        end
    end

    # TPID
    _write_tpid_line(io, label, Int(actual_mat))

    # MF1/MT451 with full directory
    _write_full_mf1(io, mf2, actual_mat, err, tempr, section_lines, ns;
                    descriptions=descriptions, thermr_mts=thermr_mts)

    # MF2/MT151
    _write_legacy_mf2(io, mf2, actual_mat, ns)

    # MF3 sections with per-MT grids
    za = mf2.ZA; awr = mf2.AWR
    for (mt, _) in reactions
        imt = Int(mt)
        ns[] = 1

        qm = 0.0; qi = 0.0; l2_head = 0; lr_tab = 0
        if haskey(override_mf3, imt)
            sec_e, sec_xs = override_mf3[imt]
            # Get QM/QI from _get_legacy_section (handles redundant MTs correctly)
            sec_orig = _get_legacy_section(result, imt)
            if sec_orig !== nothing
                qm, qi = sec_orig[3], sec_orig[4]
            end
            # Broadened sections: L2 from reconr convention
            # Fortran recout: ALL redundant/computed reactions get L2=99
            is_redundant = imt == 1 || imt == 4
            if !is_redundant && imt in (103,104,105,106,107)
                lo, hi = Dict(103=>(600,649),104=>(650,699),
                        105=>(700,749),106=>(750,799),107=>(800,849))[imt]
                is_redundant = any(s -> lo <= Int(s.mt) <= hi, result.mf3_sections)
            end
            if is_redundant
                l2_head = 99
            else
                for sec in result.mf3_sections
                    if Int(sec.mt) == imt
                        l2_head = Int(sec.L2)
                        lr_tab = Int(sec.LR)
                        break
                    end
                end
            end
        elseif haskey(extra_mf3, imt)
            sec_e, sec_xs = extra_mf3[imt]
            # Extra sections: L2=0
            l2_head = 0
            # Thermr sections get TAB1 C1=temperature (Fortran tpend convention)
            if imt in thermr_mts
                qm = tempr
            end
        else
            sec = _get_legacy_section(result, imt)
            sec === nothing && continue
            sec_e, sec_xs, qm, qi = sec
            # Reconr sections: L2 from original ENDF MF3 HEAD record
            # Redundant MTs (1, 4, 103-107): L2=0, QI=0 (Fortran recout)
            # Exception: MT=1 gets L2=99 (Fortran recout dummy MT1)
            # Fortran recout: ALL redundant/computed reactions get L2=99 (hardcoded scr(4)=99)
            is_redundant = imt == 1 || imt == 4
            if !is_redundant && imt in (103,104,105,106,107)
                lo, hi = Dict(103=>(600,649),104=>(650,699),
                        105=>(700,749),106=>(750,799),107=>(800,849))[imt]
                is_redundant = any(s -> lo <= Int(s.mt) <= hi, result.mf3_sections)
            end
            if is_redundant
                l2_head = 99
            else
                # Non-redundant: L2 from ENDF MF3 HEAD (level number for MT=51-68)
                for sec in result.mf3_sections
                    if Int(sec.mt) == imt
                        l2_head = Int(sec.L2)
                        lr_tab = Int(sec.LR)
                        break
                    end
                end
            end
        end
        _write_cont_line(io, za, awr, 0, l2_head, 0, 0, Int(actual_mat), 3, imt, ns)
        # TAB1: QM, QI, L1=0, LR from MF3 section
        np = length(sec_e)
        _write_cont_line(io, qm, qi, 0, lr_tab, 1, np, Int(actual_mat), 3, imt, ns)
        # Interpolation
        @printf(io, "%11d%11d%44s%4d%2d%3d%5d\n", np, 2, "", Int(actual_mat), 3, imt, ns[])
        ns[] += 1
        # Data pairs
        data = Float64[]
        for i in 1:np
            push!(data, sec_e[i]); push!(data, sec_xs[i])
        end
        _write_data_values(io, data, Int(actual_mat), 3, imt, ns; pair_data=true)
        _write_send(io, Int(actual_mat), 3)
    end
    _write_fend_zero(io, Int(actual_mat))

    # MF6 sections (thermal angular distributions + coherent elastic stubs)
    all_mf6_mts = sort(collect(union(keys(mf6_records), keys(mf6_stubs))))
    for mt in all_mf6_mts
        if haskey(mf6_records, mt)
            _write_mf6_section(io, Int(actual_mat), mt, za, awr, mf6_records[mt], tempr;
                               xsi=get(mf6_xsi, mt, Float64[]),
                               emax_yield=get(mf6_emax, mt, 0.0))
        elseif haskey(mf6_stubs, mt)
            s = mf6_stubs[mt]
            _write_mf6_coherent_stub(io, Int(actual_mat), mt, za, awr,
                                      s.nbragg, s.emin, s.emax)
        end
    end
    if !isempty(all_mf6_mts)
        _write_fend_zero(io, Int(actual_mat))
    end

    # MF12 passthrough (raw lines with corrected sequence numbers)
    if !isempty(mf12_lines)
        for line in mf12_lines
            print(io, line)
            endswith(line, '\n') || println(io)
        end
        _write_send(io, Int(actual_mat), 12)
        _write_fend_zero(io, Int(actual_mat))
    end

    # MF13 passthrough (raw lines with corrected sequence numbers)
    if !isempty(mf13_lines)
        for line in mf13_lines
            print(io, line)
            endswith(line, '\n') || println(io)
        end
        _write_send(io, Int(actual_mat), 13)
        _write_fend_zero(io, Int(actual_mat))
    end

    # MEND and TEND
    _write_fend_zero(io, 0)
    blanks = repeat(" ", 66)
    @printf(io, "%s%4d%2d%3d%5d\n", blanks, -1, 0, 0, 0)
end

"""Write MF1/MT451 with full directory covering MF3+MF6+MF12+MF13."""
function _write_full_mf1(io::IO, mf2, mat::Int32, err, tempr,
                          section_lines::Dict{Tuple{Int,Int}, Int}, ns::Ref{Int};
                          descriptions::Vector{String}=String[],
                          thermr_mts::Set{Int}=Set{Int}())
    za = mf2.ZA; awr = mf2.AWR
    ns[] = 1

    # Count directory entries: MF1/MT451 + MF2/MT151 + all other sections
    nxc = 2 + length(section_lines)
    nwd = length(descriptions)

    # HEAD (Fortran tpend: L1=0, L2=0)
    _write_cont_line(io, za, awr, 0, 0, 0, 0, Int(mat), 1, 451, ns)
    # Blank CONT (Fortran convention)
    _write_cont_line(io, 0.0, 0.0, 0, 0, 0, 0, Int(mat), 1, 451, ns)
    # Control: tempr, err, 0, 0, NWD, NXC
    _write_cont_line(io, tempr, err, 0, 0, nwd, nxc, Int(mat), 1, 451, ns)
    # Description text lines
    for desc in descriptions
        @printf(io, "%-66s%4d%2d%3d%5d\n", desc, Int(mat), 1, 451, ns[])
        ns[] += 1
    end

    # Directory entries (6 integers per line: MF, MT, NC, MOD per entry)
    # Self count: 3 header + NWD text + NXC directory entries + 1 SEND = 4 + NWD + NXC
    entries = Tuple{Int,Int,Int,Int}[]
    push!(entries, (1, 451, nwd + nxc, 0))  # self: NC = NWD + NXC (Fortran convention)
    push!(entries, (2, 151, 4, 0))  # MF2

    for ((mf, mt), nc) in sort(collect(section_lines))
        push!(entries, (mf, mt, nc, 0))
    end

    for entry in entries
        blanks = repeat(" ", 22)
        @printf(io, "%s%11d%11d%11d%11d%4d%2d%3d%5d\n",
                blanks, entry..., Int(mat), 1, 451, ns[])
        ns[] += 1
    end

    _write_send(io, Int(mat), 1)
    _write_fend_zero(io, Int(mat))
end

# ==========================================================================
# MF6 writer: thermal angular distributions (TAB2 + LIST records)
# Matches Fortran thermr.f90 calcem output (ltt=5, LANG=3).
# ==========================================================================

"""
    _write_mf6_section(io, mat, mt, za, awr, records, temperature)

Write one MF6 section (HEAD + product TAB1 + TAB2 + LIST records).
`records` is a Vector of MF6ListRecord from compute_mf6_thermal.
"""
function _write_mf6_section(io::IO, mat::Int, mt::Int, za::Float64, awr::Float64,
                             records::AbstractVector, temperature::Float64;
                             xsi::Vector{Float64}=Float64[], emax_yield::Float64=0.0)
    id = MaterialId(Int32(mat), Int32(6), Int32(mt))
    ne = length(records)  # number of incident energies
    ns = 1

    # HEAD: ZA, AWR, JP=0, LCT=1, NK=1, 0 (thermr.f90:1927-1934)
    ns = write_cont(io, ContRecord(za, awr, Int32(0), Int32(1), Int32(1), Int32(0), id); ns=ns)

    # Product yield TAB1: ZAP=1, AWP=1, LIP=-1, LAW=1 (thermr.f90:1936-1949)
    emin = ne > 0 ? records[1].E_incident : 1e-5
    emax_y = emax_yield > 0 ? emax_yield : (ne > 0 ? records[end].E_incident * 1.2 : 1.2)
    interp_yield = InterpolationTable(Int32[2], [LinLin])
    yield_tab1 = Tab1Record(1.0, 1.0, Int32(-1), Int32(1),
                            interp_yield, [emin, emax_y], [1.0, 1.0], id)
    ns = write_tab1(io, yield_tab1; ns=ns)

    # TAB2 outer envelope: T, 0, LANG=3, LEP=1, NR=1, NE (thermr.f90:1951-1960)
    interp_outer = InterpolationTable(Int32[ne], [LinLin])
    tab2 = Tab2Record(temperature, 0.0, Int32(3), Int32(1),
                      interp_outer, Int32(ne), id)
    ns = write_tab2(io, tab2; ns=ns)

    # Normalize sigma by xsi (Fortran tpend lines 3315,3329: rxsec=1/xsi; yy=sigma*rxsec)
    has_xsi = length(xsi) == ne

    # For each incident energy, a LIST record (thermr.f90:2225-2237)
    for (ie, rec) in enumerate(records)
        n_ep = length(rec.entries)  # number of secondary energies
        nl = 10  # values per secondary energy: E', σ, μ₁...μ₈
        nw = n_ep * nl
        rxsec = has_xsi && xsi[ie] != 0.0 ? 1.0 / xsi[ie] : 1.0
        # Pack data: flatten all entries into a single vector
        data = Float64[]
        sizehint!(data, nw)
        for entry in rec.entries
            for k in 1:nl
                if k == 2 && has_xsi
                    # Normalize sigma: Fortran tpend line 3329
                    yy = entry[k] * rxsec
                    push!(data, abs(yy) >= 1e-9 ? round_sigfig(yy, 7, 0) :
                                                   round_sigfig(yy, 6, 0))
                elseif k >= 3
                    # Cosines: sigfig(9,0) matching Fortran line 2182
                    push!(data, round_sigfig(entry[k], 9, 0))
                else
                    push!(data, entry[k])
                end
            end
        end
        list_rec = ListRecord(0.0, rec.E_incident, Int32(0), Int32(0),
                              Int32(nw), Int32(nl), data, id)
        ns = write_list(io, list_rec; ns=ns)
    end

    # SEND
    _write_send(io, mat, 6)
end

"""
    _write_mf6_coherent_stub(io, mat, mt, za, awr, nbragg, emin, emax)

Write MF6/MT=230 coherent elastic stub (4 lines).
LIP=-nbragg, LAW=0 (no angular distribution). Matches Fortran thermr coh output.
"""
function _write_mf6_coherent_stub(io::IO, mat::Int, mt::Int, za::Float64, awr::Float64,
                                    nbragg::Int, emin::Float64, emax::Float64)
    id = MaterialId(Int32(mat), Int32(6), Int32(mt))
    ns = Ref(1)
    # HEAD: ZA, AWR, 0, LCT=1, NK=1, 0
    ns[] = write_cont(io, ContRecord(za, awr, Int32(0), Int32(1), Int32(1), Int32(0), id); ns=ns[])
    # Product TAB1: ZAP=1, AWP=1, LIP=-nbragg, LAW=0
    interp = InterpolationTable(Int32[2], [LinLin])
    tab1 = Tab1Record(1.0, 1.0, Int32(-nbragg), Int32(0), interp,
                      [emin, emax], [1.0, 1.0], id)
    ns[] = write_tab1(io, tab1; ns=ns[])
    _write_send(io, mat, 6)
end

# ==========================================================================
# Low-level ENDF formatting
# ==========================================================================

function _write_tpid_line(io::IO, text::AbstractString, mat::Int)
    t = rpad(text, 66)[1:66]
    # TPID always uses MAT=1, MF=0, MT=0, NS=0 per NJOY convention
    @printf(io, "%s%4d%2d%3d%5d\n", t, 1, 0, 0, 0)
end

function _write_record_sep(io::IO, mat::Int, mf::Int, mt::Int;
                           ns::Union{Ref{Int}, Nothing}=nothing)
    blanks = repeat(" ", 66)
    seq = ns === nothing ? 99999 : ns[]
    @printf(io, "%s%4d%2d%3d%5d\n", blanks, mat, mf, mt, seq)
    if ns !== nothing
        ns[] += 1
    end
end

function _write_fend(io::IO, mat::Int; ns::Union{Ref{Int}, Nothing}=nothing)
    # FEND only: each MT section writes its own SEND; no extra SEND here
    _write_record_sep(io, mat, 0, 0; ns=ns)
end

"""Write SEND record with NS=99999 (matching Fortran convention)."""
function _write_send(io::IO, mat::Int, mf::Int)
    blanks = repeat(" ", 66)
    @printf(io, "%s%4d%2d%3d%5d\n", blanks, mat, mf, 0, 99999)
end

"""Write FEND record with NS=0 (matching Fortran convention)."""
function _write_fend_zero(io::IO, mat::Int)
    blanks = repeat(" ", 66)
    @printf(io, "%s%4d%2d%3d%5d\n", blanks, mat, 0, 0, 0)
end

function _write_cont_line(io::IO, c1, c2, l1, l2, n1, n2,
                          mat::Int, mf::Int, mt::Int, ns::Ref{Int})
    s1 = format_endf_float(Float64(c1))
    s2 = format_endf_float(Float64(c2))
    @printf(io, "%s%s%11d%11d%11d%11d%4d%2d%3d%5d\n",
            s1, s2, l1, l2, n1, n2, mat, mf, mt, ns[])
    ns[] += 1
end

function _write_data_values(io::IO, vals, mat::Int, mf::Int, mt::Int, ns::Ref{Int};
                            pair_data::Bool=false)
    n = length(vals)
    i = 1
    while i <= n
        for j in 1:6
            if i <= n
                # For (E, XS) pair data: odd positions are energies (may use 9-sigfig
                # extended format), even positions are XS (always 7-sigfig scientific).
                # Matching Fortran a11 behavior where only energies from round_sigfig(x,9)
                # have genuine 9-sigfig precision.
                ext = !pair_data || isodd(i)
                write(io, format_endf_float(Float64(vals[i]); extended=ext))
                i += 1
            else
                write(io, "           ")
            end
        end
        @printf(io, "%4d%2d%3d%5d\n", mat, mf, mt, ns[])
        ns[] += 1
    end
end

# ==========================================================================
# MF3 writer for PointwiseMaterial
# ==========================================================================

function _write_mf3_from_matrix(io::IO, mat::Int, mt::Int,
                                 energies::Vector{Float64},
                                 xs::Matrix{Float64}, col::Int,
                                 ns::Ref{Int};
                                 za::Float64=0.0, awr::Float64=0.0)
    n = length(energies)
    n >= 2 || return

    # HEAD: ZA, AWR, 0, LR=99, 0, 0 (matching reference tape format)
    _write_cont_line(io, za, awr, 0, 99, 0, 0, mat, 3, mt, ns)

    # TAB1 header: QM=0, QI=0, 0, 0, NR=1, NP=n
    _write_cont_line(io, 0.0, 0.0, 0, 0, 1, n, mat, 3, mt, ns)

    # Interpolation table: NBT=n, INT=2
    vals = Float64[Float64(n), 2.0]
    _write_data_values(io, vals, mat, 3, mt, ns)

    # Data pairs
    data = Float64[]
    sizehint!(data, 2 * n)
    for i in 1:n
        push!(data, energies[i])
        push!(data, xs[i, col])
    end
    _write_data_values(io, data, mat, 3, mt, ns; pair_data=true)

    # SEND with continuous sequence numbering
    _write_record_sep(io, mat, 3, 0; ns=ns)
end

# ==========================================================================
# MF1 section writer for PointwiseMaterial
# ==========================================================================

function _write_mf1_section(io::IO, mat::Int, material::PointwiseMaterial,
                            temperature::Float64, err::Float64,
                            description::Vector{String}, ns::Ref{Int})
    ns[] = 1
    nxc = length(material.mt_list) + 2  # MF1 + MF2 + MF3 entries
    nwd = length(description)

    # HEAD
    _write_cont_line(io, 0.0, 0.0, 0, 0, 0, 0, mat, 1, 451, ns)
    # Control
    _write_cont_line(io, temperature, err, 0, 0, nwd, nxc, mat, 1, 451, ns)

    # Description cards
    for desc in description
        padded = rpad(desc, 66)[1:66]
        @printf(io, "%s%4d%2d%3d%5d\n", padded, mat, 1, 451, ns[])
        ns[] += 1
    end

    # Dictionary entries
    n_pts = length(material.energies)
    nc_per_section = 3 + cld(n_pts, 3)
    _write_cont_line(io, 0.0, 0.0, 1, 451, 3, 0, mat, 1, 451, ns)
    _write_cont_line(io, 0.0, 0.0, 2, 151, 4, 0, mat, 1, 451, ns)
    for mt in material.mt_list
        _write_cont_line(io, 0.0, 0.0, 3, mt, nc_per_section, 0, mat, 1, 451, ns)
    end

    # SEND + FEND
    _write_record_sep(io, mat, 1, 0)
    _write_record_sep(io, mat, 0, 0)
end

# ==========================================================================
# MF2 handling
# ==========================================================================

function _write_minimal_mf2(io::IO, mat::Int, ns::Ref{Int})
    ns[] = 1
    _write_cont_line(io, 0.0, 0.0, 0, 0, 1, 0, mat, 2, 151, ns)
    _write_cont_line(io, 0.0, 0.0, 0, 0, 0, 0, mat, 2, 151, ns)
    _write_record_sep(io, mat, 2, 0)
    _write_record_sep(io, mat, 0, 0)
end

function _copy_mf2_from_source(io_out::IO, io_src::IO, mat::Int)
    seekstart(io_src)
    copying = false
    while !eof(io_src)
        line = readline(io_src)
        p = rpad(line, 80)
        mf = _pendf_parse_int(p[71:72])
        mat_line = _pendf_parse_int(p[67:70])

        if mat_line == mat && mf == 2
            copying = true
        end

        if copying
            println(io_out, rstrip(p))
            if mf == 0 && copying
                break
            end
        end
    end

    if !copying
        ns = Ref(1)
        _write_minimal_mf2(io_out, mat, ns)
    end
end

# ==========================================================================
# Legacy interface (NamedTuple result)
# ==========================================================================

function _collect_reactions(result)
    # Collect reactions in ENDF file order (matching Fortran emerge/recout).
    # For MF=3: MT=1 (total) first; for MF=23: MT=501 (photoatomic total) first.
    reactions = Tuple{Int32, String}[]
    is_photoatomic = !isempty(result.mf3_sections) &&
                     all(s -> Int(s.mf) == 23, result.mf3_sections)
    if is_photoatomic
        push!(reactions, (Int32(501), "total"))
    else
        push!(reactions, (Int32(1), "total"))
    end

    has_inelastic = any(s -> Int(s.mt) >= 51 && Int(s.mt) <= 91, result.mf3_sections)
    has_partial_fission = any(s -> Int(s.mt) == 19, result.mf3_sections)

    # Detect redundant charged-particle groups (Fortran anlyzd lines 567-590)
    _has_range(lo, hi) = any(s -> lo <= Int(s.mt) <= hi, result.mf3_sections)
    has_mp = _has_range(600, 649)   # MT=103 redundant
    has_md = _has_range(650, 699)   # MT=104 redundant
    has_mt5 = _has_range(700, 749)  # MT=105 redundant
    has_m3 = _has_range(750, 799)   # MT=106 redundant
    has_m4 = _has_range(800, 849)   # MT=107 redundant

    for sec in result.mf3_sections
        mt = Int(sec.mt)
        # Skip truly redundant MTs (never output by Fortran)
        (mt == 1 || mt == 3 || mt == 101 || mt == 120 ||
         mt == 151 || mt == 27 || (mt >= 251 && mt <= 300 && mt != 261)) && continue
        mt == 501 && continue  # photoatomic total (redundant, handled above)
        mt == 460 && continue  # delayed photon total
        # MT=4: output as redundant sum if inelastic levels present
        if mt == 4
            has_inelastic && push!(reactions, (Int32(4), "inelastic"))
            continue
        end
        # Skip MT=103-107 when partials exist (Fortran lunion lines 1884-1888)
        (mt == 103 && has_mp) && continue
        (mt == 104 && has_md) && continue
        (mt == 105 && has_mt5) && continue
        (mt == 106 && has_m3) && continue
        (mt == 107 && has_m4) && continue
        # MT=19: insert redundant MT=18 before it (matching Fortran recout order)
        if mt == 19 && has_partial_fission
            push!(reactions, (Int32(18), "fission"))
        end
        push!(reactions, (Int32(mt), "MT$mt"))
    end

    # Insert redundant charged-particle MTs at their natural MT position
    # (Fortran recout interleaves redundant MTs by sorted MT number)
    for (flag, redmt) in [(has_mp, 103), (has_md, 104), (has_mt5, 105),
                           (has_m3, 106), (has_m4, 107)]
        flag || continue
        any(r -> r[1] == Int32(redmt), reactions) && continue
        idx = findfirst(r -> r[1] > Int32(redmt), reactions)
        entry = (Int32(redmt), "MT$redmt")
        idx === nothing ? push!(reactions, entry) : insert!(reactions, idx, entry)
    end

    return reactions
end

function _write_legacy_mf1(io::IO, mf2::MF2Data, mat::Int32, err, tempr,
                            n_pts, reactions, ns::Ref{Int};
                            mt_to_mf::Dict{Int,Int}=Dict{Int,Int}(),
                            descriptions::Vector{String}=String[])
    ns[] = 1
    nxc = 2 + length(reactions)
    nwd = length(descriptions)

    _write_cont_line(io, mf2.ZA, mf2.AWR, 2, 0, 0, 0,
                     Int(mat), 1, 451, ns)
    _write_cont_line(io, tempr, err, 0, 0, nwd, nxc,
                     Int(mat), 1, 451, ns)

    # Description text lines (NWD lines)
    for desc in descriptions
        @printf(io, "%-66s%4d%2d%3d%5d\n", rpad(desc, 66)[1:66], Int(mat), 1, 451, ns[])
        ns[] += 1
    end

    nc_per = 3 + cld(n_pts, 3)
    nc_self = 2 + nwd + nxc
    _write_cont_line(io, 0.0, 0.0, 1, 451, nc_self, 0, Int(mat), 1, 451, ns)
    _write_cont_line(io, 0.0, 0.0, 2, 151, 4, 0, Int(mat), 1, 451, ns)
    for (mt, _) in reactions
        dir_mf = get(mt_to_mf, Int(mt), 3)
        _write_cont_line(io, 0.0, 0.0, dir_mf, Int(mt), nc_per, 0,
                         Int(mat), 1, 451, ns)
    end

    _write_send(io, Int(mat), 1)
    _write_fend_zero(io, Int(mat))
end

function _write_legacy_mf2(io::IO, mf2::MF2Data, mat::Int32, ns::Ref{Int})
    ns[] = 1
    nis = length(mf2.isotopes)
    _write_cont_line(io, mf2.ZA, mf2.AWR, 0, 0, nis, 0, Int(mat), 2, 151, ns)

    for iso in mf2.isotopes
        # Fortran recout uses ZA (not ZAI) for the isotope CONT record
        _write_cont_line(io, mf2.ZA, iso.ABN, 0, 0, 1, 0, Int(mat), 2, 151, ns)
        el = length(iso.ranges) > 0 ? iso.ranges[1].EL : 1.0e-5
        # Fortran recout uses eresr for EH: resolved range EH if LRU=1 exists,
        # otherwise ehigh=2e7 for LRU=0-only materials
        eh = 2.0e7  # ehigh: default for no-resonance materials
        if length(iso.ranges) > 0
            resolved_ehs = [r.EH for r in iso.ranges if r.LRU == 1]
            if !isempty(resolved_ehs)
                eh = maximum(resolved_ehs)
            end
            # If only LRU=0 ranges, keep eh = 2.0e7 (Fortran eresr = ehigh)
        end
        _write_cont_line(io, el, eh, 0, 0, 0, 0, Int(mat), 2, 151, ns)
        spi = 1.0; ap = 0.0
        if length(iso.ranges) > 0
            p = iso.ranges[1].parameters
            if hasproperty(p, :SPI); spi = p.SPI; end
            if hasproperty(p, :AP); ap = p.AP; end
        end
        _write_cont_line(io, spi, ap, 0, 0, 0, 0, Int(mat), 2, 151, ns)
    end

    _write_send(io, Int(mat), 2)
    _write_fend_zero(io, Int(mat))
end

function _write_legacy_mf3(io::IO, result, mat::Int32, reactions, ns::Ref{Int})
    energies = result.energies
    awr = result.mf2.AWR
    za = result.mf2.ZA

    # Build MF lookup: MT → actual MF from section data
    mt_to_mf = Dict{Int, Int}()
    for sec in result.mf3_sections
        mt_to_mf[Int(sec.mt)] = Int(sec.mf)
    end

    cur_mf = 0

    for (mt, _) in reactions
        ns[] = 1

        sec_data = _get_legacy_section(result, Int(mt))
        sec_data === nothing && continue

        sec_e, sec_xs, qm, qi = sec_data

        # Determine actual MF (3 for neutron, 23 for photoatomic)
        out_mf = get(mt_to_mf, Int(mt), 3)

        # Handle MF transition: write FEND when MF changes
        if cur_mf != 0 && out_mf != cur_mf
            _write_fend_zero(io, Int(mat))
        end
        cur_mf = out_mf

        # HEAD: ZA, AWR, 0, L2=99, 0, 0
        _write_cont_line(io, za, awr, 0, 99, 0, 0,
                         Int(mat), out_mf, Int(mt), ns)

        # TAB1: QM, QI, 0, LR, NR=1, NP
        qi_out = (mt == 1 || mt == 4 || mt == 501) ? 0.0 : qi
        np = length(sec_e)
        _write_cont_line(io, qm, qi_out, 0, 0, 1, np,
                         Int(mat), out_mf, Int(mt), ns)

        # Interpolation table
        @printf(io, "%11d%11d%44s%4d%2d%3d%5d\n",
                np, 2, "", Int(mat), out_mf, Int(mt), ns[])
        ns[] += 1

        # Data pairs
        data = Float64[]
        sizehint!(data, 2 * np)
        for i in 1:np
            push!(data, sec_e[i])
            push!(data, sec_xs[i])
        end
        _write_data_values(io, data, Int(mat), out_mf, Int(mt), ns; pair_data=true)

        _write_send(io, Int(mat), out_mf)
    end

    if cur_mf != 0
        _write_fend_zero(io, Int(mat))
    end
end

"""Get energies, XS, QM, QI for a given MT, handling thresholds."""
function _get_legacy_section(result, mt::Int)
    energies = result.energies
    awr = result.mf2.AWR

    # Look up QM/QI from the MF3 section data
    qm, qi = 0.0, 0.0
    for sec in result.mf3_sections
        if Int(sec.mt) == mt
            qm, qi = sec.QM, sec.QI
            break
        end
    end

    if mt == 1
        return energies, result.total, 0.0, 0.0
    elseif mt == 501
        # Photoatomic total: MT=501 = sum of all MF=23 except MT=515,517
        return energies, result.total, 0.0, 0.0
    elseif mt == 4
        # Redundant sum: MT=4 = sum(MT=51-91)
        # Use the MINIMUM threshold from any inelastic level (not just
        # the first). Fortran emerge accumulates from ALL levels; mtrt
        # tracking picks up wherever ANY contribution first appears.
        first_inel = findfirst(s -> Int(s.mt) >= 51 && Int(s.mt) <= 91, result.mf3_sections)
        first_inel === nothing && return nothing
        qi_4 = result.mf3_sections[first_inel].QI
        min_thrxx = Inf
        for sec in result.mf3_sections
            smt = Int(sec.mt)
            (smt < 51 || smt > 91) && continue
            s_qi = sec.QI
            if s_qi < 0.0
                s_awr = sec.awr > 0.0 ? sec.awr : awr
                s_thrx = s_awr > 0.0 ? -s_qi * (s_awr + 1) / s_awr : -s_qi
                s_thrxx = round_sigfig(s_thrx, 7, +1)
                min_thrxx = min(min_thrxx, s_thrxx)
            end
        end
        thrxx_4 = isfinite(min_thrxx) ? min_thrxx : 0.0
        # Build grid: all energies at or above lowest threshold
        sec_e = Float64[]
        sec_xs = Float64[]
        for e in energies
            thrxx_4 > 0.0 && thrxx_4 - e > 1.0e-9 * thrxx_4 && continue
            # Sum all MT=51-91 at this energy
            total_inel = 0.0
            for sec in result.mf3_sections
                smt = Int(sec.mt)
                (smt < 51 || smt > 91) && continue
                # Per-level threshold handling (matching Fortran emerge)
                s_qi = sec.QI
                s_thrxx = 0.0
                if s_qi < 0.0
                    s_awr = sec.awr > 0.0 ? sec.awr : awr
                    s_thrx = s_awr > 0.0 ? -s_qi * (s_awr + 1) / s_awr : -s_qi
                    s_thrxx = round_sigfig(s_thrx, 7, +1)
                    # Below threshold → skip (Fortran emerge line 4792)
                    if s_thrxx > 1.0 && (s_thrxx - e) > 1.0e-10 * s_thrxx
                        continue
                    end
                    # At threshold → zero (Fortran emerge line 4795)
                    if s_thrxx > 0.0 && abs(s_thrxx - e) < 1.0e-9 * s_thrxx
                        continue
                    end
                end
                bg = if s_thrxx > 0.0 && e >= s_thrxx && sec.tab.x[1] < s_thrxx
                    _threshold_interp(sec.tab, e, s_thrxx)
                else
                    interpolate(sec.tab, e)
                end
                total_inel += round_sigfig(bg, 7)
            end
            push!(sec_e, e)
            push!(sec_xs, round_sigfig(total_inel, 7))
        end
        isempty(sec_e) && return nothing
        # Pseudo-threshold: skip leading zero values (matching Fortran
        # recout ith tracking — start at the point before first nonzero)
        first_nz = findfirst(>(0.0), sec_xs)
        if first_nz !== nothing && first_nz > 1
            sec_e = sec_e[first_nz-1:end]
            sec_xs = sec_xs[first_nz-1:end]
        end
        isempty(sec_e) && return nothing
        return sec_e, sec_xs, 0.0, 0.0  # Fortran recout writes QI=0 for redundant MT=4
    elseif mt in (103, 104, 105, 106, 107) && begin
            # Only handle as redundant sum when partials exist
            _lo, _hi = Dict(103=>(600,649), 104=>(650,699),
                105=>(700,749), 106=>(750,799), 107=>(800,849))[mt]
            any(s -> _lo <= Int(s.mt) <= _hi, result.mf3_sections)
        end
        # Redundant sum: MT=103-107 from charged-particle partials
        # (Fortran recout lines 5270-5318, anlyzd lines 567-590)
        # MT=103=sum(600-649), 104=sum(650-699), 105=sum(700-749),
        # 106=sum(750-799), 107=sum(800-849)
        part_lo, part_hi = Dict(103=>(600,649), 104=>(650,699),
            105=>(700,749), 106=>(750,799), 107=>(800,849))[mt]
        partials = filter(s -> part_lo <= Int(s.mt) <= part_hi, result.mf3_sections)
        sec_e = Float64[]
        sec_xs = Float64[]
        has_rxs = hasproperty(result, :reaction_xs)
        for (idx, e) in enumerate(energies)
            total_cp = 0.0
            for sec in partials
                smt = Int(sec.mt)
                s_qi = sec.QI
                s_thrxx = 0.0
                if s_qi < 0.0
                    s_awr = sec.awr > 0.0 ? sec.awr : awr
                    s_thrx = s_awr > 0.0 ? -s_qi * (s_awr + 1) / s_awr : -s_qi
                    s_thrxx = round_sigfig(s_thrx, 7, +1)
                    # Below threshold → skip
                    if s_thrxx > 1.0 && (s_thrxx - e) > 1.0e-10 * s_thrxx
                        continue
                    end
                    # At threshold → zero
                    if s_thrxx > 0.0 && abs(s_thrxx - e) < 1.0e-9 * s_thrxx
                        continue
                    end
                end
                bg = if s_thrxx > 0.0 && e >= s_thrxx && sec.tab.x[1] < s_thrxx
                    _threshold_interp(sec.tab, e, s_thrxx)
                else
                    interpolate(sec.tab, e)
                end
                # Add RML reaction channel contribution (matching Fortran
                # emerge: sn = gety1(bg) + res(1+itype) for reaction channels)
                if has_rxs && haskey(result.reaction_xs, smt)
                    bg += result.reaction_xs[smt][idx]
                end
                total_cp += round_sigfig(bg, 7)
            end
            push!(sec_e, e)
            push!(sec_xs, round_sigfig(total_cp, 7))
        end
        isempty(sec_e) && return nothing
        first_nz = findfirst(>(0.0), sec_xs)
        if first_nz !== nothing && first_nz > 1
            sec_e = sec_e[first_nz-1:end]
            sec_xs = sec_xs[first_nz-1:end]
        end
        isempty(sec_e) && return nothing
        return sec_e, sec_xs, 0.0, 0.0
    elseif mt == 2
        return energies, result.elastic, qm, qi
    elseif mt == 18
        # Redundant sum: MT=18 = MT=19 + MT=20 + MT=21 + MT=38
        # (matching Fortran emerge accumulation lines 4886-4893)
        has_partial = any(s -> Int(s.mt) == 19, result.mf3_sections)
        if has_partial
            # Get QM/QI from MT=19 (Fortran saves q18 from MT=19, line 1916)
            for sec in result.mf3_sections
                if Int(sec.mt) == 19
                    qm, qi = sec.QM, sec.QI
                    break
                end
            end
            fission_sum = copy(result.fission)
            for sec in result.mf3_sections
                mt_sec = Int(sec.mt)
                (mt_sec == 20 || mt_sec == 21 || mt_sec == 38) || continue
                s_qi = sec.QI
                s_thrxx = 0.0
                if s_qi < 0.0
                    s_awr = sec.awr > 0.0 ? sec.awr : awr
                    s_thrx = s_awr > 0.0 ? -s_qi * (s_awr + 1) / s_awr : -s_qi
                    s_thrxx = round_sigfig(s_thrx, 7, +1)
                end
                for (i, e) in enumerate(energies)
                    if s_thrxx > 0.0 && abs(s_thrxx - e) < 1.0e-7 * s_thrxx
                        continue  # threshold suppression
                    end
                    bg = interpolate(sec.tab, e)
                    if s_thrxx > 0.0 && e >= s_thrxx && length(sec.tab.x) >= 2 &&
                       e < sec.tab.x[2] && sec.tab.x[1] < s_thrxx
                        bg = _threshold_interp(sec.tab, e, s_thrxx)
                    end
                    if bg != 0.0
                        fission_sum[i] += round_sigfig(bg, 7)
                    end
                end
            end
            # Final sigfig on accumulated sum (Fortran recout line 5308)
            for i in eachindex(fission_sum)
                fission_sum[i] = round_sigfig(fission_sum[i], 7)
            end
            return energies, fission_sum, qm, qi
        else
            return energies, result.fission, qm, qi
        end
    elseif mt == 19
        return energies, result.fission, qm, qi
    elseif mt == 102
        return energies, result.capture, qm, qi
    end

    # Non-primary MT: find the MF3 section
    for sec in result.mf3_sections
        Int(sec.mt) != mt && continue

        qm = sec.QM
        qi = sec.QI

        # Compute threshold for reactions with Q < 0
        thrxx = 0.0
        # thresh: Fortran emerge's threshold = sigfig(tab.x[1], 7, 0)
        # from gety1 initialization (line 4753). Used for BOTH the
        # below-threshold skip (line 4792) AND the at-threshold
        # suppression (line 4794). When tab.x[1] > physical threshold
        # (thrx), Fortran effectively uses the MF3 start energy, not
        # the Q-value threshold.
        thresh = 0.0
        if qi < 0.0
            # Use per-section AWR (from MF3 HEAD) matching Fortran emerge
            # which reads awrx from each section's CONT record (line 4753)
            sec_awr = sec.awr > 0.0 ? sec.awr : awr
            thrx = sec_awr > 0.0 ? -qi * (sec_awr + 1) / sec_awr : -qi
            thrxx = round_sigfig(thrx, 7, +1)
            thresh = round_sigfig(sec.tab.x[1], 7)
            if thresh < thrxx
                thresh = thrxx  # physical threshold is higher
            end
        end

        # Build the section output matching Fortran emerge:
        # 1. Skip points below threshold (line 4792)
        # 2. Set xs=0 at threshold energy (line 4795)
        # 3. Pseudo-threshold skip via ith logic (lines 4834, 4918)
        sec_e = Float64[]
        sec_xs = Float64[]
        for (idx, e) in enumerate(energies)
            if thresh > 0.0 && thresh - e > 1.0e-10 * thresh
                continue  # below threshold (Fortran line 4792)
            end
            # Near threshold, modify first breakpoint to (thrxx, 0) and adjust
            # subsequent breakpoints that fall below thrxx (Fortran lunion
            # lines 1936-1943), then interpolate using the MF3 interpolation law.
            # The condition extends beyond x[2] because thrxx can exceed x[2].
            bg = if thrxx > 0.0 && e >= thrxx && sec.tab.x[1] < thrxx
                _threshold_interp(sec.tab, e, thrxx)
            else
                interpolate(sec.tab, e)
            end
            # For RML reaction channels, add the resonance contribution
            # (matching Fortran emerge: sn = gety1(bg) + res(1+itype))
            if hasproperty(result, :reaction_xs) && haskey(result.reaction_xs, mt)
                bg += result.reaction_xs[mt][idx]
            end
            bg = round_sigfig(bg, 7)

            # Set xs=0 at threshold (line 4795): Fortran uses
            # thresh=sigfig(tab.x[1],7,0), NOT thrxx=sigfig(thrx,7,+1).
            # Since thrxx > thresh, the grid point at thrxx is NOT
            # suppressed, matching Fortran behavior.
            if thresh > 1.0 && abs(thresh - e) < 1.0e-10 * thresh
                bg = 0.0
            end

            push!(sec_e, e)
            push!(sec_xs, bg)
        end

        isempty(sec_e) && return nothing

        # Pseudo-threshold skip: matching Fortran emerge ith logic
        # (lines 4834, 4918). Find first nonzero XS (ith), back up one
        # if possible, output from there.
        ith = findfirst(x -> x > 0.0, sec_xs)
        if ith !== nothing && ith > 1
            ith -= 1  # Fortran line 4918: if (ith.gt.1) ith=ith-1
            sec_e = sec_e[ith:end]
            sec_xs = sec_xs[ith:end]
        elseif ith === nothing
            return nothing  # all zero
        end

        return sec_e, sec_xs, qm, qi
    end
    return nothing
end

function _pendf_parse_int(s::AbstractString)
    t = strip(s)
    isempty(t) && return 0
    return parse(Int, t)
end

"""
    _threshold_interp(tab, e, thrxx)

Threshold-adjusted interpolation: temporarily modify the first breakpoint
to (thrxx, 0.0) and interpolate using the MF3's own interpolation law.
Matches Fortran emerge line 1936 (modify ex(1)) + gety1.
"""
function _threshold_interp(tab::TabulatedFunction, e::Float64, thrxx::Float64)
    # Save original breakpoints that will be modified
    n = length(tab.x)
    saved = Tuple{Int,Float64}[]
    # Replace x[1] with threshold (Fortran lunion line 1936)
    push!(saved, (1, tab.x[1]))
    tab.x[1] = thrxx
    push!(saved, (-1, tab.y[1]))  # use -1 as sentinel for y[1]
    tab.y[1] = 0.0
    # Adjust subsequent breakpoints that are now <= the threshold
    # (Fortran lunion lines 1937-1943: do while scr(l+2) <= scr(l))
    k = 1
    while k < n && tab.x[k+1] <= tab.x[k]
        push!(saved, (k+1, tab.x[k+1]))
        tab.x[k+1] = round_sigfig(tab.x[k], 7, +1)
        k += 1
        k > 21 && break  # safety limit matching Fortran lsave+20
    end
    result = interpolate(tab, e)
    # Restore all modified breakpoints
    for (idx, val) in reverse(saved)
        if idx == -1
            tab.y[1] = val
        else
            tab.x[idx] = val
        end
    end
    return result
end
