# ENDF-6 format I/O routines
# Reference: NJOY2016 endf.f90 (contio, listio, tab1io, tab2io, a11)
#
# ENDF lines are 80 columns: 6x11-char data fields, then I4+I2+I3+I5 trailer.
# Float fields omit 'E': "1.234567+8" means 1.234567e8.

"""
    parse_endf_float(s::AbstractString) -> Float64

Parse an 11-character ENDF float field. Handles compact notation without 'E'.
"""
function parse_endf_float(s::AbstractString)
    t = strip(s)
    isempty(t) && return 0.0
    if occursin(r"[eEdD]", t)
        return parse(Float64, replace(t, r"[dD]" => "e"))
    end
    for i in lastindex(t):-1:2
        c = t[i]
        if (c == '+' || c == '-') && (isdigit(t[i-1]) || t[i-1] == '.')
            return parse(Float64, t[1:i-1] * "e" * t[i:end])
        end
    end
    return parse(Float64, t)
end

"""
    format_endf_float(x::Real) -> String

Format a number as an 11-character ENDF float field (matching NJOY's `a11`).
"""
function format_endf_float(x::Real)
    x == 0.0 && return " 0.000000+0"
    ax = abs(Float64(x))
    n = floor(Int, log10(ax))
    f = x / 10.0^n
    esign, en = n >= 0 ? ('+', n) : ('-', -n)
    if abs(f) >= 9.9999995
        f /= 10.0
        if esign == '+'; en += 1
        else en -= 1; if en < 0; esign = '+'; en = -en; end
        end
    end
    if en < 10;      return @sprintf("%9.6f", f) * esign * string(en)
    elseif en < 100; return @sprintf("%8.5f", f) * esign * string(en)
    else;            return @sprintf("%7.4f", f) * esign * string(en)
    end
end

# --- Line parsing ---

function parse_endf_line(line::AbstractString)
    p = rpad(line, 80)
    fields = ntuple(i -> p[(i-1)*11+1 : i*11], 6)
    mat = _parse_int(p[67:70])
    mf  = _parse_int(p[71:72])
    mt  = _parse_int(p[73:75])
    ns  = _parse_int(p[76:80])
    return (fields, mat, mf, mt, ns)
end

_parse_int(s::AbstractString) = (t = strip(s); isempty(t) ? Int32(0) : parse(Int32, t))

# --- Readers ---

"""
    read_tpid(io::IO) -> NamedTuple{(:text, :mat, :mf, :mt)}

Read the TPID (tape identification) record — the first line of an ENDF tape.
Returns the 66-character text and the MAT/MF/MT identifiers.
"""
function read_tpid(io::IO)
    line = readline(io)
    p = rpad(line, 80)
    text = p[1:66]
    mat = _parse_int(p[67:70])
    mf  = _parse_int(p[71:72])
    mt  = _parse_int(p[73:75])
    return (text=text, mat=mat, mf=mf, mt=mt)
end

"""
    find_section(io::IO, target_mf::Integer, target_mt::Integer) -> Bool

Scan from the beginning of `io` until a line with the given MF/MT is found.
Leaves the stream positioned at the start of that line (the HEAD record).
Returns `true` if found, `false` if EOF is reached.
"""
function find_section(io::IO, target_mf::Integer, target_mt::Integer)
    seekstart(io)
    while !eof(io)
        pos = position(io)
        line = readline(io)
        p = rpad(line, 80)
        mf = _parse_int(p[71:72])
        mt = _parse_int(p[73:75])
        if mf == target_mf && mt == target_mt
            seek(io, pos)
            return true
        end
    end
    return false
end

"""
    read_cont(io::IO) -> ContRecord

Read one CONT record (single 80-column line).
"""
function read_cont(io::IO)
    fields, mat, mf, mt, _ = parse_endf_line(readline(io))
    ContRecord(parse_endf_float(fields[1]), parse_endf_float(fields[2]),
               _parse_int(fields[3]), _parse_int(fields[4]),
               _parse_int(fields[5]), _parse_int(fields[6]),
               MaterialId(mat, mf, mt))
end

"Read a HEAD record (same format as CONT)."
read_head(io::IO) = read_cont(io)

"""Read a LIST record: CONT header + ceil(NPL/6) data lines."""
function read_list(io::IO)
    head = read_cont(io)
    npl = Int(head.N1)
    data = Vector{Float64}(undef, npl)
    k = 0
    while k < npl
        fields, _, _, _, _ = parse_endf_line(readline(io))
        for j in 1:min(6, npl - k)
            k += 1; data[k] = parse_endf_float(fields[j])
        end
    end
    ListRecord(head.C1, head.C2, head.L1, head.L2, head.N1, head.N2, data, head.id)
end

# Shared helper: read NR pairs of (NBT, INT) for interpolation tables
function _read_interp_table(io::IO, nr::Int)
    nbt = Vector{Int32}(undef, nr)
    int_arr = Vector{Int32}(undef, nr)
    k = 0
    while k < nr
        fields, _, _, _, _ = parse_endf_line(readline(io))
        for j in 1:min(3, nr - k)
            k += 1
            nbt[k] = _parse_int(fields[2j-1])
            int_arr[k] = _parse_int(fields[2j])
        end
    end
    InterpolationTable(nbt, int_arr)
end

"""Read a TAB1 record: CONT header + interp table + (x,y) data."""
function read_tab1(io::IO)
    head = read_cont(io)
    nr, np = Int(head.N1), Int(head.N2)
    interp = _read_interp_table(io, nr)
    xarr = Vector{Float64}(undef, np)
    yarr = Vector{Float64}(undef, np)
    k = 0
    while k < np
        fields, _, _, _, _ = parse_endf_line(readline(io))
        for j in 1:min(3, np - k)
            k += 1
            xarr[k] = parse_endf_float(fields[2j-1])
            yarr[k] = parse_endf_float(fields[2j])
        end
    end
    Tab1Record(head.C1, head.C2, head.L1, head.L2, interp, xarr, yarr, head.id)
end

"""Read a TAB2 record: CONT header + interp table."""
function read_tab2(io::IO)
    head = read_cont(io)
    interp = _read_interp_table(io, Int(head.N1))
    Tab2Record(head.C1, head.C2, head.L1, head.L2, interp, head.N2, head.id)
end

# --- Writers ---

function _write_line(io::IO, parts::Vector{String}, id::MaterialId, ns::Int)
    while length(parts) < 6; push!(parts, "           "); end
    println(io, join(parts) * @sprintf("%4d%2d%3d%5d", id.mat, id.mf, id.mt, ns))
    ns + 1
end

function _write_data_line(io::IO, vals, id::MaterialId, ns::Int)
    _write_line(io, [format_endf_float(Float64(v)) for v in vals], id, ns)
end

function _write_int_line(io::IO, vals, id::MaterialId, ns::Int)
    _write_line(io, [lpad(string(Int(v)), 11) for v in vals], id, ns)
end

"""Write a CONT record."""
function write_cont(io::IO, rec::ContRecord; ns::Int=0)
    line = format_endf_float(rec.C1) * format_endf_float(rec.C2) *
           lpad(string(Int(rec.L1)), 11) * lpad(string(Int(rec.L2)), 11) *
           lpad(string(Int(rec.N1)), 11) * lpad(string(Int(rec.N2)), 11) *
           @sprintf("%4d%2d%3d%5d", rec.id.mat, rec.id.mf, rec.id.mt, ns)
    println(io, line)
    ns + 1
end

# Shared helper: write interpolation table
function _write_interp_table(io::IO, interp::InterpolationTable, id::MaterialId, ns::Int)
    nr = length(interp)
    k = 0
    while k < nr
        bs = min(3, nr - k)
        vals = Int[]
        for j in 1:bs; push!(vals, Int(interp.nbt[k+j]), Int(interp.law[k+j])); end
        ns = _write_int_line(io, vals, id, ns)
        k += bs
    end
    ns
end

"""Write a LIST record."""
function write_list(io::IO, rec::ListRecord; ns::Int=0)
    ns = write_cont(io, ContRecord(rec.C1, rec.C2, rec.L1, rec.L2,
                                   rec.N1, rec.N2, rec.id); ns=ns)
    npl = Int(rec.N1)
    k = 0
    while k < npl
        ns = _write_data_line(io, rec.data[k+1:min(k+6, npl)], rec.id, ns)
        k += 6
    end
    ns
end

"""Write a TAB1 record."""
function write_tab1(io::IO, rec::Tab1Record; ns::Int=0)
    nr, np = length(rec.interp), length(rec.x)
    ns = write_cont(io, ContRecord(rec.C1, rec.C2, rec.L1, rec.L2,
                                   Int32(nr), Int32(np), rec.id); ns=ns)
    ns = _write_interp_table(io, rec.interp, rec.id, ns)
    k = 0
    while k < np
        bs = min(3, np - k)
        vals = Float64[]
        for j in 1:bs; push!(vals, rec.x[k+j], rec.y[k+j]); end
        ns = _write_data_line(io, vals, rec.id, ns)
        k += bs
    end
    ns
end

"""Write a TAB2 record."""
function write_tab2(io::IO, rec::Tab2Record; ns::Int=0)
    nr = length(rec.interp)
    ns = write_cont(io, ContRecord(rec.C1, rec.C2, rec.L1, rec.L2,
                                   Int32(nr), rec.NZ, rec.id); ns=ns)
    _write_interp_table(io, rec.interp, rec.id, ns)
end
