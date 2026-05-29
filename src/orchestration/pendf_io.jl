# pendf_io.jl -- PENDF tape read/write for the module orchestration layer
#
# Reads a PENDF tape into PENDFTape struct, writes it back.
# Provides extract_mf3 for parsing cross-section data from raw lines,
# and copy_with_modifications for the Fortran-style "stream input to output,
# replacing/adding sections" pattern.
#
# Adapts the existing read_endf_tape/write_endf_tape (moder.jl) for
# PENDF-specific handling (MF1/MT451 directory separation).

using Printf

# =========================================================================
# Reading
# =========================================================================

"""
    read_pendf(path::AbstractString) -> PENDFTape

Read a PENDF tape file into a structured PENDFTape.
Separates MF1/MT451 header from other sections.
"""
function read_pendf(path::AbstractString)::PENDFTape
    materials = PENDFMaterial[]
    tpid = ""

    open(path, "r") do io
        # First line is TPID
        tpid_line = readline(io)
        tpid = rpad(tpid_line, 80)[1:66]

        cur_mat = 0; cur_mf = 0; cur_mt = 0
        cur_lines = String[]
        mf1_lines = String[]
        sections = PENDFSection[]

        function flush_section!()
            if cur_mt > 0 && !isempty(cur_lines)
                if cur_mf == 1 && cur_mt == 451
                    append!(mf1_lines, cur_lines)
                else
                    push!(sections, PENDFSection(cur_mf, cur_mt, copy(cur_lines)))
                end
            end
            empty!(cur_lines)
        end

        function flush_material!()
            flush_section!()
            if cur_mat > 0 && (!isempty(mf1_lines) || !isempty(sections))
                push!(materials, PENDFMaterial(cur_mat, copy(mf1_lines), copy(sections)))
            end
            empty!(mf1_lines); empty!(sections)
        end

        for raw in eachline(io)
            line = rpad(raw, 80)
            mat = _parse_int(line[67:70])
            mf  = _parse_int(line[71:72])
            mt  = _parse_int(line[73:75])

            # TEND record: end of tape
            if mat < 0
                flush_material!()
                break
            end

            # MEND record: end of material
            if mat == 0 && mf == 0 && mt == 0
                flush_material!()
                cur_mat = 0; cur_mf = 0; cur_mt = 0
                continue
            end

            # FEND record: end of file (MF)
            if mf == 0 && mt == 0 && mat > 0
                flush_section!()
                cur_mf = 0; cur_mt = 0
                continue
            end

            # SEND record: end of section (MT). Treated as a structural
            # marker — we do NOT keep it in section data. The writer
            # re-emits SEND/FEND/MEND/TEND from the structural shape.
            if mt == 0 && mf > 0 && mat > 0
                flush_section!()
                cur_mt = 0
                continue
            end

            # New material
            if mat != cur_mat && mat > 0
                flush_material!()
                cur_mat = mat
            end

            # New MF
            if mf != cur_mf
                flush_section!()
                cur_mf = mf
            end

            # New MT
            if mt != cur_mt
                flush_section!()
                cur_mt = mt
            end

            push!(cur_lines, line)
        end
        flush_material!()
    end

    PENDFTape(tpid, materials)
end

# =========================================================================
# Writing
# =========================================================================

"""
    write_pendf_tape(path::AbstractString, tape::PENDFTape)

Write a PENDFTape to a file with proper TPID, FEND, MEND, TEND records.
Sequence numbers (columns 76-80) are regenerated.
"""
function write_pendf_tape(path::AbstractString, tape::PENDFTape)
    open(path, "w") do io
        # TPID line: ENDF spec says MAT=1, MF=0, MT=0, seq=0.
        tpid_padded = rpad(tape.tpid, 66)[1:66]
        @printf(io, "%s%4d%2d%3d%5d\n", tpid_padded, 1, 0, 0, 0)

        for material in tape.materials
            # MF1/MT451 — seq restarts at 1; close with SEND only. The MF
            # transition (or the trailing FEND below) emits the FEND.
            prev_mf = 0
            if !isempty(material.mf1_lines)
                seq = 0
                for line in material.mf1_lines
                    seq += 1
                    _write_line(io, line, material.mat, seq)
                end
                _write_send(io, material.mat, 1)
                prev_mf = 1
            end

            for sec in material.sections
                if sec.mf != prev_mf && prev_mf > 0
                    _write_fend(io, material.mat)
                end
                prev_mf = sec.mf

                seq = 0
                for line in sec.lines
                    seq += 1
                    _write_line(io, line, material.mat, seq)
                end
                _write_send(io, material.mat, sec.mf)
            end

            # Trailing FEND for the last MF written.
            prev_mf > 0 && _write_fend(io, material.mat)

            # MEND
            @printf(io, "%66s%4d%2d%3d%5d\n", "", 0, 0, 0, 0)
        end

        # TEND
        @printf(io, "%66s%4d%2d%3d%5d\n", "", -1, 0, 0, 0)
    end
end

function _write_line(io::IO, line::AbstractString, mat::Int, seq::Int)
    # Preserve cols 1-75 (data + MAT/MF/MT); regenerate seq in cols 76-80.
    data = rpad(line, 80)
    @printf(io, "%s%5d\n", data[1:75], seq)
end

# SEND/FEND helpers reused from src/processing/pendf_writer.jl:
#   _write_send(io, mat, mf)  — SEND record (MAT/MF/0, seq=99999)
#   _write_fend_zero(io, mat) — FEND record (MAT/0/0, seq=0)
const _write_fend = _write_fend_zero

# =========================================================================
# Data Extraction
# =========================================================================

"""
    extract_mf3(tape::PENDFTape, mat::Int, mt::Int) -> (energies, xs)

Parse MF3 data from raw PENDF lines for a given MAT/MT.
Returns (energies::Vector{Float64}, xs::Vector{Float64}).
Returns (Float64[], Float64[]) if section not found.
"""
function extract_mf3(tape::PENDFTape, mat::Int, mt::Int)
    for material in tape.materials
        material.mat != mat && continue
        for sec in material.sections
            (sec.mf != 3 || sec.mt != mt) && continue
            return _parse_mf3_lines(sec.lines)
        end
    end
    return (Float64[], Float64[])
end

"""
    extract_mf3_all(tape::PENDFTape, mat::Int) -> Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}

Extract all MF3 sections for a material, keyed by MT number.
"""
function extract_mf3_all(tape::PENDFTape, mat::Int)
    result = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
    for material in tape.materials
        material.mat != mat && continue
        for sec in material.sections
            sec.mf != 3 && continue
            e, xs = _parse_mf3_lines(sec.lines)
            isempty(e) && continue
            result[sec.mt] = (e, xs)
        end
    end
    result
end

"""
Extract the temperature stored in a PENDFMaterial's MF1/MT451 block.
Locates the TEMP CONT record by walking lines 2..5 and selecting the last
line where (L2=0, N1>0, 0 ≤ C1 < 1e6) — TEMP CONT is the last numeric CONT
before description Hollerith. Handles both ENDF-IV layout (TEMP at line 3 of
the tape, mf1_lines[2]) and ENDF-V/VI layout (TEMP at mf1_lines[3] or [4]).
Returns 0.0 if no plausible TEMP CONT is found.
"""
function _pendf_material_temperature(material::PENDFMaterial)::Float64
    isempty(material.mf1_lines) && return 0.0
    found_temp = 0.0
    found_any = false
    for li in 2:min(5, length(material.mf1_lines))
        p = rpad(material.mf1_lines[li], 80)
        # Parse integer CONT fields first; description text lines won't parse
        # and are skipped before we attempt the float parse on cols 1-11.
        l1 = tryparse(Int, strip(p[23:33]))
        l2 = tryparse(Int, strip(p[34:44]))
        n1 = tryparse(Int, strip(p[45:55]))
        (l1 === nothing || l2 === nothing || n1 === nothing) && continue
        c1 = try
            parse_endf_float(p[1:11])
        catch
            continue
        end
        if l2 == 0 && n1 > 0 && c1 >= 0.0 && c1 < 1.0e6
            found_temp = c1
            found_any = true
        end
    end
    found_any ? found_temp : 0.0
end

"""
    extract_mf3_at_temperature(tape::PENDFTape, mat::Int, target_temp::Float64;
                                tol::Float64=0.5) -> Dict{Int, Tuple{...}}

Extract all MF3 sections for the material entry whose stored TEMP matches
`target_temp` within `tol` Kelvin. Multi-T PENDFs (broadr output) place each
temperature in a separate MEND-bounded material entry. If no exact match is
found and the target is the only candidate (single-T PENDF), return its
sections regardless. Errors loudly per Rule 6 if no candidate exists.
"""
function extract_mf3_at_temperature(tape::PENDFTape, mat::Int, target_temp::Float64;
                                     tol::Float64=0.5)
    candidates = [m for m in tape.materials if m.mat == mat]
    isempty(candidates) && error(
        "extract_mf3_at_temperature: no PENDF material with mat=$mat (target_temp=$target_temp K)")
    chosen = nothing
    for m in candidates
        t = _pendf_material_temperature(m)
        if abs(t - target_temp) <= max(tol, 1e-3 * target_temp)
            chosen = m
            break
        end
    end
    if chosen === nothing
        if length(candidates) == 1
            chosen = candidates[1]   # single-T PENDF: return what we have
        else
            avail = join([string(_pendf_material_temperature(m), " K") for m in candidates], ", ")
            error("extract_mf3_at_temperature: mat=$mat target_temp=$target_temp K not found "
                  * "among PENDF temperatures [$avail]")
        end
    end
    result = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
    for sec in chosen.sections
        sec.mf != 3 && continue
        e, xs = _parse_mf3_lines(sec.lines)
        isempty(e) && continue
        result[sec.mt] = (e, xs)
    end
    result
end

"""Parse MF3 TAB1 data from raw ENDF lines into (energies, xs) vectors."""
function _parse_mf3_lines(lines::Vector{String})
    energies = Float64[]
    xs = Float64[]
    isempty(lines) && return (energies, xs)

    # Line 1: HEAD (ZA, AWR, 0, 0, 0, 0) — skip
    # Line 2: TAB1 header (QM, QI, 0, LR, NR, NP)
    length(lines) < 3 && return (energies, xs)

    p2 = rpad(lines[2], 80)
    nr = _parse_int(p2[45:55])
    np = _parse_int(p2[56:66])
    np <= 0 && return (energies, xs)

    # Skip interpolation lines: ceil(2*NR / 6) lines
    n_interp_lines = cld(2 * max(nr, 0), 6)
    data_start = 3 + n_interp_lines  # 1-indexed

    # Parse data lines: 6 values per line, alternating E and XS
    sizehint!(energies, np)
    sizehint!(xs, np)
    for li in data_start:length(lines)
        p = rpad(lines[li], 80)
        # Check for SEND record
        mt_check = _parse_int(p[73:75])
        mt_check == 0 && break
        # Parse up to 3 (E, XS) pairs per line
        for col in 0:2
            length(energies) >= np && break
            e_str = p[1 + col*22 : 11 + col*22]
            x_str = p[12 + col*22 : 22 + col*22]
            push!(energies, parse_endf_float(e_str))
            push!(xs, parse_endf_float(x_str))
        end
    end

    (energies, xs)
end

"""
    _parse_mf3_tab1(lines) -> NamedTuple(qm, qi, lr, energies, xs, nbt, law)

Parse an MF3 section's TAB1 **with** its interpolation table and Q-values.
Unlike `_parse_mf3_lines` (which discards the interpolation law and assumes
lin-lin), this preserves the native ENDF interpolation regions (NBT/INT) and
the QM/QI fields. Required for the charged-particle ACE path, where MF3
reaction cross sections are stored on coarse grids with non-linear laws
(INT=5 log-log, INT=6 Coulomb penetrability) and must be sampled onto the
finer ESZ grid via `gety1`/`terp1` rather than linear interpolation.

The reconr/broadr neutron PENDF linearizes everything to INT=2, so the
neutron path can keep using `_parse_mf3_lines`; this is only consumed by
the `izai != neutron` branch in acer_module.

Returns empty vectors on a malformed/empty section.

Ref: ENDF-6 §3.2 (MF3 TAB1 layout); njoy-reference/src/acefc.f90:5494-5505
(acelod's per-reaction gety1 sampling + sigfig-7 store).
"""
function _parse_mf3_tab1(lines::Vector{String})
    empty_nbt = Int32[]; empty_law = Int32[]
    energies = Float64[]; xs = Float64[]
    isempty(lines) && return (; qm=0.0, qi=0.0, lr=0, energies, xs,
                                nbt=empty_nbt, law=empty_law)
    length(lines) < 3 && return (; qm=0.0, qi=0.0, lr=0, energies, xs,
                                   nbt=empty_nbt, law=empty_law)

    # Line 2: TAB1 CONT (QM, QI, L1, LR, NR, NP).
    p2 = rpad(lines[2], 80)
    qm = parse_endf_float(p2[1:11])
    qi = parse_endf_float(p2[12:22])
    lr = _parse_int(p2[34:44])
    nr = _parse_int(p2[45:55])
    np = _parse_int(p2[56:66])
    np <= 0 && return (; qm, qi, lr, energies, xs, nbt=empty_nbt, law=empty_law)

    # Interpolation table: NR (NBT, INT) integer pairs, 3 pairs per line.
    n_interp_lines = cld(2 * max(nr, 0), 6)
    nbt = Int32[]; law = Int32[]
    sizehint!(nbt, nr); sizehint!(law, nr)
    interp_collected = 0
    for li in 3:(2 + n_interp_lines)
        li > length(lines) && break
        p = rpad(lines[li], 80)
        for col in 0:2
            interp_collected >= nr && break
            nbt_v = _parse_int(p[1 + col*22 : 11 + col*22])
            int_v = _parse_int(p[12 + col*22 : 22 + col*22])
            push!(nbt, nbt_v); push!(law, int_v)
            interp_collected += 1
        end
    end

    data_start = 3 + n_interp_lines
    sizehint!(energies, np); sizehint!(xs, np)
    for li in data_start:length(lines)
        p = rpad(lines[li], 80)
        _parse_int(p[73:75]) == 0 && break   # SEND
        for col in 0:2
            length(energies) >= np && break
            push!(energies, parse_endf_float(p[1 + col*22 : 11 + col*22]))
            push!(xs,       parse_endf_float(p[12 + col*22 : 22 + col*22]))
        end
    end

    (; qm, qi, lr, energies, xs, nbt, law)
end

"""
    extract_mf3_tab1_all(tape::PENDFTape, mat::Int)
        -> Dict{Int, NamedTuple(qm, qi, lr, tab::TabulatedFunction)}

Charged-particle counterpart to `extract_mf3_all`: returns every MF3 section
for `mat` as a `TabulatedFunction` carrying its native interpolation table,
plus the QM/QI Q-values. The charged ACE path samples these via
`interpolate` (which honours INT=5/INT=6) to reproduce Fortran `gety1`.
"""
function extract_mf3_tab1_all(tape::PENDFTape, mat::Int)
    result = Dict{Int, NamedTuple{(:qm, :qi, :lr, :tab),
                  Tuple{Float64, Float64, Int, TabulatedFunction}}}()
    for material in tape.materials
        material.mat != mat && continue
        for sec in material.sections
            sec.mf != 3 && continue
            rec = _parse_mf3_tab1(sec.lines)
            isempty(rec.energies) && continue
            # Fall back to a single lin-lin region if NR was unreadable.
            nbt = isempty(rec.nbt) ? Int32[length(rec.energies)] : rec.nbt
            law = isempty(rec.law) ? Int32[2] : rec.law
            tab = TabulatedFunction(InterpolationTable(nbt, law),
                                    rec.energies, rec.xs)
            result[sec.mt] = (; qm=rec.qm, qi=rec.qi, lr=rec.lr, tab)
        end
    end
    result
end

# =========================================================================
# Copy With Modifications
# =========================================================================

"""
    copy_with_modifications(tape::PENDFTape, mat::Int;
        modified_mf3=Dict(), added_mf3=Dict(), added_mf6=Dict(),
        mf6_stubs=Dict(), mf12_lines=String[], mf13_lines=String[],
        temperature=0.0, descriptions=String[],
        thermr_mts=Set{Int}(), thermr_coh_ne=0) -> PENDFTape

Build a new PENDFTape by copying `tape` with specified modifications:
- `modified_mf3`: Dict{MT => (energies, xs)} — replace existing MF3 sections
- `added_mf3`: Dict{MT => (energies, xs)} — add new MF3 sections
- `added_mf6`: Dict{MT => lines} — add new MF6 sections (raw lines)
- Other kwargs for future extension

This is the core of the Fortran "copy input PENDF to output, modifying as you go" pattern.
"""
function copy_with_modifications(tape::PENDFTape, mat::Int;
        modified_mf3::Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}} = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}(),
        added_mf3::Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}} = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}(),
        added_mf6_lines::Dict{Int, Vector{String}} = Dict{Int, Vector{String}}(),
        mf6_stub_lines::Dict{Int, Vector{String}} = Dict{Int, Vector{String}}(),
        added_mf12_lines::Vector{String} = String[],
        added_mf13_lines::Vector{String} = String[],
        temperature::Float64 = 0.0)

    new_materials = PENDFMaterial[]

    for material in tape.materials
        if material.mat != mat
            push!(new_materials, material)
            continue
        end

        new_sections = PENDFSection[]

        for sec in material.sections
            if sec.mf == 3 && haskey(modified_mf3, sec.mt)
                # Replace this MF3 section with broadened/modified data
                e, xs = modified_mf3[sec.mt]
                new_lines = _build_mf3_lines(sec.lines, e, xs, mat)
                push!(new_sections, PENDFSection(3, sec.mt, new_lines))
            else
                # Copy section unchanged
                push!(new_sections, sec)
            end
        end

        # Add new MF3 sections (sorted by MT, inserted after existing MF3)
        for mt in sort(collect(keys(added_mf3)))
            e, xs = added_mf3[mt]
            new_lines = _build_new_mf3_section(e, xs, mat, mt)
            push!(new_sections, PENDFSection(3, mt, new_lines))
        end

        # Add new MF6 sections
        for mt in sort(collect(keys(added_mf6_lines)))
            push!(new_sections, PENDFSection(6, mt, added_mf6_lines[mt]))
        end

        # Add MF6 stubs
        for mt in sort(collect(keys(mf6_stub_lines)))
            push!(new_sections, PENDFSection(6, mt, mf6_stub_lines[mt]))
        end

        # Add MF12/MF13 passthrough lines as sections
        if !isempty(added_mf12_lines)
            push!(new_sections, PENDFSection(12, 102, added_mf12_lines))
        end
        if !isempty(added_mf13_lines)
            push!(new_sections, PENDFSection(13, 51, added_mf13_lines))
        end

        # Sort sections by (MF, MT) to maintain ENDF order
        sort!(new_sections, by = s -> (s.mf, s.mt))

        # TODO: update MF1/MT451 directory to reflect added/modified sections
        # For now, carry the original MF1 header
        push!(new_materials, PENDFMaterial(mat, copy(material.mf1_lines), new_sections))
    end

    PENDFTape(tape.tpid, new_materials)
end

"""
Build MF3 data lines from (energies, xs) arrays, preserving the HEAD
line from the original section.
"""
function _build_mf3_lines(orig_lines::Vector{String}, energies::Vector{Float64},
                          xs::Vector{Float64}, mat::Int)
    lines = String[]
    np = length(energies)
    np == 0 && return lines

    # HEAD line — copy from original (ZA, AWR, L1, L2, N1, N2)
    push!(lines, orig_lines[1])

    # TAB1 header — copy QM, QI, L1, LR from original, update NR=1, NP
    p2 = rpad(orig_lines[2], 80)
    qm = p2[1:11]; qi = p2[12:22]; l1 = p2[23:33]; lr = p2[34:44]
    nr_str = lpad("1", 11); np_str = lpad(string(np), 11)
    mf_mt_str = @sprintf("%4d%2d%3d", mat, 3, _parse_int(p2[73:75]))
    push!(lines, qm * qi * l1 * lr * nr_str * np_str * mf_mt_str * "    0")

    # Interpolation line: NR=1 region, NP points, INT=2 (linear-linear)
    interp_line = lpad(string(np), 11) * lpad("2", 11) * " "^44 *
                  @sprintf("%4d%2d%3d", mat, 3, _parse_int(p2[73:75])) * "    0"
    push!(lines, interp_line)

    # Data lines: 3 (E, XS) pairs per line. SEND is added by writer.
    mt = _parse_int(p2[73:75])
    trailer = @sprintf("%4d%2d%3d", mat, 3, mt)
    idx = 1
    while idx <= np
        buf = ""
        for col in 1:3
            idx > np && break
            buf *= format_endf_float(energies[idx]) * format_endf_float(xs[idx])
            idx += 1
        end
        push!(lines, rpad(buf, 66) * trailer * "    0")
    end

    lines
end

"""
Build a complete new MF3 section with HEAD + TAB1. SEND is emitted by writer.
"""
function _build_new_mf3_section(energies::Vector{Float64}, xs::Vector{Float64},
                                mat::Int, mt::Int)
    np = length(energies)
    lines = String[]
    trailer = @sprintf("%4d%2d%3d", mat, 3, mt)

    # HEAD line: ZA=0, AWR=0, L1=0, L2=0, N1=0, N2=0
    head = format_endf_float(0.0) * format_endf_float(0.0) *
           lpad("0", 11) * lpad("0", 11) * lpad("0", 11) * lpad("0", 11)
    push!(lines, head * trailer * "    0")

    # TAB1 header: QM=0, QI=0, L1=0, LR=0, NR=1, NP
    tab1h = format_endf_float(0.0) * format_endf_float(0.0) *
            lpad("0", 11) * lpad("0", 11) * lpad("1", 11) * lpad(string(np), 11)
    push!(lines, tab1h * trailer * "    0")

    # Interpolation line
    interp = lpad(string(np), 11) * lpad("2", 11) * " "^44
    push!(lines, interp * trailer * "    0")

    # Data lines (SEND is emitted by writer)
    idx = 1
    while idx <= np
        buf = ""
        for col in 1:3
            idx > np && break
            buf *= format_endf_float(energies[idx]) * format_endf_float(xs[idx])
            idx += 1
        end
        push!(lines, rpad(buf, 66) * trailer * "    0")
    end

    lines
end
