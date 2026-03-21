# MODER -- Tape management and format conversion for ENDF-format data
#
# In NJOY Fortran, MODER converts between binary and ASCII formats and
# merges/extracts materials. Since NJOY.jl works exclusively with ASCII ENDF,
# MODER reduces to copy/reformat/material-selection utilities.
#
# Correspondence to NJOY2016 moder.f90:
#   tpidio + loop  -> read_tape_directory (scan for MAT/MF/MT)
#   matd selection -> extract_material / moder_copy (copy with MAT filter)
#   multi-tape cat -> merge_tapes (concatenate)
#   tape read/write -> read_endf_tape / write_endf_tape (structured I/O)

# ==========================================================================
# Types
# ==========================================================================

"One section on an ENDF tape, identified by MAT/MF/MT numbers."
struct TapeEntry
    mat::Int; mf::Int; mt::Int
end

"Directory of an ENDF tape: tape-ID string and section list."
struct TapeDirectory
    tpid::String
    entries::Vector{TapeEntry}
end

Base.length(td::TapeDirectory) = length(td.entries)
Base.iterate(td::TapeDirectory, args...) = iterate(td.entries, args...)
Base.show(io::IO, td::TapeDirectory) =
    print(io, "TapeDirectory(\"", strip(td.tpid), "\", ", length(td.entries), " sections)")

"Return sorted unique MAT numbers on the tape."
materials(dir::TapeDirectory) = sort!(unique(e.mat for e in dir.entries))

"Return (MF,MT) pairs for a given material."
sections(dir::TapeDirectory, mat::Integer) =
    [(e.mf, e.mt) for e in dir.entries if e.mat == mat]

"One ENDF section stored as raw lines for lossless round-tripping."
struct ENDFTapeSection
    mf::Int; mt::Int; lines::Vector{String}
end

"One complete ENDF material from a tape."
struct ENDFTapeMaterial
    mat::Int; sections::Vector{ENDFTapeSection}
end

# ==========================================================================
# read_tape_directory
# ==========================================================================

"""
    read_tape_directory(io::IO) -> TapeDirectory

Scan an ENDF tape and return a directory listing every MAT/MF/MT section.
"""
function read_tape_directory(io::IO)
    tpid_text = ""; entries = TapeEntry[]
    first_line = true; seen = Set{Tuple{Int,Int,Int}}()
    for line in eachline(io)
        p = rpad(line, 80)
        if first_line; tpid_text = p[1:66]; first_line = false; continue; end
        mat = _parse_int(p[67:70]); mf = _parse_int(p[71:72]); mt = _parse_int(p[73:75])
        (mat <= 0 || mf <= 0 || mt <= 0) && continue
        key = (mat, mf, mt)
        key in seen && continue
        push!(seen, key); push!(entries, TapeEntry(mat, mf, mt))
    end
    TapeDirectory(tpid_text, entries)
end

# ==========================================================================
# extract_material
# ==========================================================================

"""
    extract_material(io::IO, mat::Integer) -> String

Return all lines for material `mat` as a string (including FEND/MEND).
"""
function extract_material(io::IO, mat::Integer)
    buf = IOBuffer(); in_mat = false
    for line in eachline(io)
        p = rpad(line, 80); line_mat = _parse_int(p[67:70])
        if line_mat == mat
            in_mat = true; println(buf, line)
        elseif in_mat
            println(buf, line)
            line_mat == 0 && break
        end
    end
    String(take!(buf))
end

"""
    extract_material(io::IO, mat::Integer, out::IO) -> Int

Stream version: write extracted material lines directly to `out`.
"""
function extract_material(io::IO, mat::Integer, out::IO)
    nlines = 0; in_mat = false
    for line in eachline(io)
        p = rpad(line, 80); line_mat = _parse_int(p[67:70])
        if line_mat == mat
            in_mat = true; println(out, line); nlines += 1
        elseif in_mat
            println(out, line); nlines += 1
            line_mat == 0 && break
        end
    end
    nlines
end

# ==========================================================================
# moder_copy: stream copy with optional MAT filter
# ==========================================================================

"""
    moder_copy(input::IO, output::IO; mat=0) -> Int

Copy ENDF material(s) from `input` to `output`, optionally filtering to a
single MAT. When `mat=0`, all materials are copied. Returns lines written.
"""
function moder_copy(input::IO, output::IO; mat::Integer=0)
    nlines = 0; in_mat = false
    for line in eachline(input)
        padded = rpad(line, 80); line_mat = _parse_int(padded[67:70])
        if mat == 0
            println(output, padded); nlines += 1
        elseif line_mat == mat
            in_mat = true; println(output, padded); nlines += 1
        elseif in_mat && line_mat == 0
            mf = _parse_int(padded[71:72])
            println(output, padded); nlines += 1
            mf == 0 && (in_mat = false)
        elseif line_mat == -1
            println(output, padded); nlines += 1
        end
    end
    nlines
end

# ==========================================================================
# merge_tapes
# ==========================================================================

"Write a TPID record."
function write_tpid(out::IO, text::AbstractString="")
    println(out, rpad(text, 66)[1:66], "   0 0  0    0")
end

"Write a TEND record."
function write_tend(out::IO)
    println(out, format_endf_float(0.0) * format_endf_float(0.0) *
           lpad("0", 11) * lpad("0", 11) * lpad("0", 11) * lpad("0", 11) *
           @sprintf("%4d%2d%3d%5d", -1, 0, 0, 0))
end

"""
    merge_tapes(out::IO, inputs::IO...; tpid="NJOY.jl merged tape")

Concatenate all materials from `inputs` onto `out` with TPID/TEND framing.
"""
function merge_tapes(out::IO, inputs::IO...; tpid::AbstractString="NJOY.jl merged tape")
    write_tpid(out, tpid)
    for inp in inputs
        first_line = true
        for line in eachline(inp)
            if first_line; first_line = false; continue; end
            p = rpad(line, 80)
            _parse_int(p[67:70]) == -1 && continue
            println(out, line)
        end
    end
    write_tend(out)
end

"File-path convenience: merge multiple ENDF files into one."
function merge_tapes(output_path::AbstractString, input_paths::AbstractString...;
                     tpid::AbstractString="NJOY.jl merged tape")
    open(output_path, "w") do out
        ios = [open(p, "r") for p in input_paths]
        try; merge_tapes(out, ios...; tpid=tpid)
        finally; for f in ios; close(f); end; end
    end
end

# ==========================================================================
# read_endf_tape / write_endf_tape: structured round-trip I/O
# ==========================================================================

"""
    read_endf_tape(filename::AbstractString) -> Vector{ENDFTapeMaterial}

Read an ENDF file and return all materials as structured data with raw
lines preserved for lossless round-tripping.
"""
function read_endf_tape(filename::AbstractString)
    mats = ENDFTapeMaterial[]
    open(filename, "r") do io
        cur_mat = 0; cur_mf = 0; cur_mt = 0
        sl = String[]; ms = ENDFTapeSection[]
        flush_s!() = begin
            cur_mt > 0 && !isempty(sl) && push!(ms, ENDFTapeSection(cur_mf, cur_mt, copy(sl)))
            empty!(sl)
        end
        flush_m!() = begin
            flush_s!()
            cur_mat > 0 && !isempty(ms) && push!(mats, ENDFTapeMaterial(cur_mat, copy(ms)))
            empty!(ms)
        end
        for raw in eachline(io)
            line = rpad(raw, 80)
            mat = Int(_parse_int(line[67:70]))
            mf  = Int(_parse_int(line[71:72]))
            mt  = Int(_parse_int(line[73:75]))
            if mat < 0; flush_m!(); break; end
            if mat == 0 && mf == 0 && mt == 0
                flush_m!(); cur_mat = 0; cur_mf = 0; cur_mt = 0; continue
            end
            if mf == 0 && mt == 0 && mat > 0; flush_s!(); cur_mf = 0; cur_mt = 0; continue; end
            if mt == 0 && mf > 0 && mat > 0; push!(sl, line); flush_s!(); cur_mt = 0; continue; end
            if mat != cur_mat && mat > 0; flush_m!(); cur_mat = mat; end
            if mf != cur_mf; flush_s!(); cur_mf = mf; end
            if mt != cur_mt; flush_s!(); cur_mt = mt; end
            push!(sl, line)
        end
        flush_m!()
    end
    mats
end

"""
    write_endf_tape(filename::AbstractString, materials::Vector{ENDFTapeMaterial})

Write materials to an ENDF-format file with FEND, MEND, and TEND records.
"""
function write_endf_tape(filename::AbstractString, materials::Vector{ENDFTapeMaterial})
    open(filename, "w") do io
        for md in materials
            prev_mf = 0
            for sec in md.sections
                if sec.mf != prev_mf && prev_mf > 0
                    println(io, " "^66 * @sprintf("%4d%2d%3d%5d", md.mat, 0, 0, 0))
                end
                prev_mf = sec.mf
                for line in sec.lines; println(io, rpad(line, 80)); end
            end
            prev_mf > 0 && println(io, " "^66 * @sprintf("%4d%2d%3d%5d", md.mat, 0, 0, 0))
            println(io, " "^66 * @sprintf("%4d%2d%3d%5d", 0, 0, 0, 0))
        end
        println(io, " "^66 * @sprintf("%4d%2d%3d%5d", -1, 0, 0, 0))
    end
end

"""
    validate_tape(io::IO) -> NamedTuple{(:valid, :errors)}

Quick structural validation of an ENDF tape.
"""
function validate_tape(io::IO)
    errs = String[]; lineno = 0
    for line in eachline(io)
        lineno += 1; p = rpad(line, 80)
        _parse_int(p[67:70]) == -1 && break
    end
    lineno == 0 && push!(errs, "empty tape")
    (valid=isempty(errs), errors=errs)
end
