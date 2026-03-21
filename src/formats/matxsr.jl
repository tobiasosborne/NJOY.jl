# MATXSR -- MATXS (Material Cross Section) interface format output
# Correspondence to NJOY2016 matxsr.f90: cmatxs -> write_matxs
# Supports ASCII and binary (Fortran unformatted) output modes.

using Printf

# =========================================================================
# ASCII format helpers
# =========================================================================

"""Format Hollerith items as MATXS ASCII lines (tag + 8-char words, 9 per continuation)."""
function _matxs_hollerith_line(tag::AbstractString, items::Vector{<:AbstractString})
    lines = String[]
    isempty(items) && (push!(lines, tag); return join(lines, "\n"))
    parts = [rpad(tag, 4)]
    for i in 1:min(8, length(items)); push!(parts, rpad(items[i], 8)[1:8]); end
    push!(lines, join(parts, ""))
    idx = min(8, length(items)) + 1
    while idx <= length(items)
        chunk = String[]
        for _ in 1:9; idx > length(items) && break; push!(chunk, rpad(items[idx], 8)[1:8]); idx += 1; end
        push!(lines, join(chunk, ""))
    end
    join(lines, "\n")
end

"""Format integers in i6 format, 12 per line."""
function _matxs_int_line(values::Vector{<:Integer}; per_line::Int=12)
    lines = String[]
    for i in 1:per_line:length(values)
        push!(lines, join([@sprintf("%6d", v) for v in values[i:min(i+per_line-1, length(values))]], ""))
    end
    join(lines, "\n")
end

"""Format reals in 1pE12.5: tag line with 5 values, continuation with 6."""
function _matxs_real_line(tag::AbstractString, values::Vector{<:Real})
    lines = String[]
    isempty(values) && (push!(lines, tag); return join(lines, "\n"))
    n1 = min(5, length(values))
    push!(lines, rpad(tag, 4) * "        " * join([@sprintf("%12.5E", Float64(v)) for v in values[1:n1]], ""))
    idx = n1 + 1
    while idx <= length(values)
        ce = min(idx + 5, length(values))
        push!(lines, join([@sprintf("%12.5E", Float64(v)) for v in values[idx:ce]], ""))
        idx = ce + 1
    end
    join(lines, "\n")
end

"""Convert ENDF MT number to MATXS-style reaction name."""
function _mt_to_matxs_name(mt::Integer)
    names = Dict{Int,String}(1=>"ntot", 2=>"nelas", 4=>"ninel", 16=>"n2n",
        17=>"n3n", 18=>"nfiss", 102=>"ngamma", 103=>"np", 104=>"nd",
        105=>"nt", 106=>"nhe3", 107=>"nalpha", 251=>"mubar", 452=>"nubar")
    get(names, Int(mt), @sprintf("mt%d", mt))
end

# =========================================================================
# ASCII MATXS writer
# =========================================================================
"""
    write_matxs(io::IO, multigroup_data::MultiGroupXS; title, ivers, particle, data_type, material_name, mat_id, atom_mass)

Write a MATXS ASCII file. Records: file identification, file control, set hollerith ID,
file data, group structure, material control, vector control, vector block.
"""
function write_matxs(io::IO, multigroup_data::MultiGroupXS;
                     title::AbstractString="NJOY.jl", ivers::Integer=0,
                     particle::AbstractString="n", data_type::AbstractString="nscat",
                     material_name::AbstractString="mat1", mat_id::Integer=0,
                     atom_mass::Real=1.0)
    gb = multigroup_data.group_bounds; ngroup = length(gb) - 1
    mt_list = multigroup_data.mt_list; xs = multigroup_data.xs
    nr = length(mt_list); maxw = 5000
    # Record 1: File identification
    println(io, @sprintf(" 0v %s*%s%s*%6d", rpad("matxs",8)[1:8],
            rpad("NJOY.jl",8)[1:8], rpad("",8)[1:8], ivers))
    # Record 2: File control
    println(io, @sprintf(" 1d   %6d%6d%6d%6d%6d%6d", 1, 1, 1, 1, maxw, 8))
    # Record 3: Set hollerith ID
    pt = rpad(title, 72)[1:72]
    println(io, " 2d "); println(io, join([pt[(i-1)*8+1:i*8] for i in 1:9], ""))
    # Record 4: File data
    h = [rpad(particle,8)[1:8], rpad(data_type,8)[1:8], rpad(material_name,8)[1:8]]
    println(io, " 3d     " * join(h, ""))
    println(io, _matxs_int_line(Int[ngroup, 1, 1, 1, 1]))
    # Record 5: Group structure
    gpb = [gb[ngroup-g+2] for g in 1:ngroup]
    println(io, _matxs_real_line(" 4d ", Float64.(vcat(gpb, [gb[1]]))))
    # Record 6: Material control
    println(io, @sprintf(" 5d %s%12.5E", rpad(material_name,8)[1:8], Float64(atom_mass)))
    println(io, @sprintf("%12.5E%12.5E%6d%6d%6d%6d", 0.0, 1.0e10, 1, nr, 0, 1))
    # Record 7: Vector control
    println(io, _matxs_hollerith_line(" 6d ", [_mt_to_matxs_name(mt) for mt in mt_list]))
    println(io, _matxs_int_line(vcat(ones(Int, nr), fill(ngroup, nr))))
    # Record 8: Vector block
    vdata = Float64[xs[g, r] for r in 1:nr for g in 1:ngroup]
    println(io, _matxs_real_line(" 7d ", vdata))
    return nothing
end

# =========================================================================
# Binary MATXS writer (Fortran unformatted)
# =========================================================================
"""
    write_matxs_binary(io::IO, multigroup_data::MultiGroupXS; ...)

Write a MATXS binary (Fortran unformatted) file with the same record structure as write_matxs.
"""
function write_matxs_binary(io::IO, multigroup_data::MultiGroupXS;
                            title::AbstractString="NJOY.jl", ivers::Integer=0,
                            particle::AbstractString="n", data_type::AbstractString="nscat",
                            material_name::AbstractString="mat1", mat_id::Integer=0,
                            atom_mass::Real=1.0)
    gb = multigroup_data.group_bounds; ngroup = length(gb) - 1
    mt_list = multigroup_data.mt_list; xs = multigroup_data.xs; nr = length(mt_list)
    # Record 1: File identification
    _write_record(io, _record_buf("matxs  ", "NJOY.jl ", "        ", Int32(ivers)))
    # Record 2: File control
    _write_record(io, _record_buf(Int32[1, 1, 1, 1, 5000, 8]))
    # Record 3: Set hollerith ID
    pt = rpad(title, 72)[1:72]
    _write_record(io, _record_buf([pt[(i-1)*8+1:i*8] for i in 1:9]...))
    # Record 4: File data
    buf4 = IOBuffer()
    for s in [particle, data_type, material_name]; write(buf4, codeunits(_pad8(s))); end
    for v in Int32[ngroup, 1, 1, 1, 1]; write(buf4, htol(v)); end
    _write_record(io, take!(buf4))
    # Record 5: Group structure
    buf5 = IOBuffer()
    for g in 1:ngroup; write(buf5, htol(Float32(gb[ngroup-g+2]))); end
    write(buf5, htol(Float32(gb[1])))
    _write_record(io, take!(buf5))
    # Record 6: Material control
    buf6 = IOBuffer()
    write(buf6, codeunits(_pad8(material_name)))
    write(buf6, htol(Float32(atom_mass))); write(buf6, htol(Float32(0.0)))
    write(buf6, htol(Float32(1.0e10)))
    for v in Int32[1, nr, 0, 1]; write(buf6, htol(v)); end
    _write_record(io, take!(buf6))
    # Record 7: Vector control
    buf7 = IOBuffer()
    for mt in mt_list; write(buf7, codeunits(_pad8(_mt_to_matxs_name(mt)))); end
    for _ in mt_list; write(buf7, htol(Int32(1))); end
    for _ in mt_list; write(buf7, htol(Int32(ngroup))); end
    _write_record(io, take!(buf7))
    # Record 8: Vector block
    buf8 = IOBuffer()
    for r in 1:nr, g in 1:ngroup; write(buf8, htol(Float32(xs[g, r]))); end
    _write_record(io, take!(buf8))
    return nothing
end
