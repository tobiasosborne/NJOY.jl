# MATXSR -- MATXS (Material Cross Section) interface format output
#
# Produces MATXS files from multigroup cross section data. MATXS is the NJOY
# generalized multigroup interface format based on CCCC conventions.
#
# Correspondence to NJOY2016 matxsr.f90:
#   cmatxs   -> write_matxs       (main MATXS file writer)
#   ruinm    -> (input embedded in arguments)
#   mtxdat   -> (data preparation embedded)
#
# The MATXS format supports both binary (Fortran unformatted) and ASCII
# (formatted text) modes. This implementation produces ASCII output,
# matching the nmatx>0 branch of NJOY2016.
#
# File structure (from matxsr.f90 spec):
#   1. File identification
#   2. File control
#   3. Set hollerith identification
#   4. File data (particle/type/material names, group counts, pointers)
#   5. Group structures (one per particle)
#   6. Material control (one per material)
#   7. Vector control + vector blocks (cross section vectors)
#   8. Matrix control + matrix sub-blocks (scattering matrices)

using Printf

# =========================================================================
# ASCII MATXS format helpers
# =========================================================================

"""
    _matxs_hollerith_line(tag, items::Vector{<:AbstractString})

Format Hollerith items as MATXS ASCII lines: tag line followed by
continuation lines of 9 items per line (a8 each).
"""
function _matxs_hollerith_line(tag::AbstractString,
                               items::Vector{<:AbstractString})
    lines = String[]
    if isempty(items)
        push!(lines, tag)
        return join(lines, "\n")
    end
    # First line: tag + up to 8 items
    first_n = min(8, length(items))
    parts = [rpad(tag, 4)]
    for i in 1:first_n
        push!(parts, rpad(items[i], 8)[1:8])
    end
    push!(lines, join(parts, ""))
    # Continuation lines: 9 items each
    idx = first_n + 1
    while idx <= length(items)
        chunk = String[]
        for _ in 1:9
            idx > length(items) && break
            push!(chunk, rpad(items[idx], 8)[1:8])
            idx += 1
        end
        push!(lines, join(chunk, ""))
    end
    join(lines, "\n")
end

"""
    _matxs_int_line(values::Vector{<:Integer}; per_line=12)

Format integers in MATXS i6 format, 12 per line.
"""
function _matxs_int_line(values::Vector{<:Integer}; per_line::Int=12)
    lines = String[]
    for i in 1:per_line:length(values)
        chunk_end = min(i + per_line - 1, length(values))
        push!(lines, join([@sprintf("%6d", v) for v in values[i:chunk_end]], ""))
    end
    join(lines, "\n")
end

"""
    _matxs_real_line(tag::AbstractString, values::Vector{<:Real})

Format real values in MATXS 1pE12.5 format: tag line with 5 values,
continuation lines with 6 values each.
"""
function _matxs_real_line(tag::AbstractString, values::Vector{<:Real})
    lines = String[]
    if isempty(values)
        push!(lines, tag)
        return join(lines, "\n")
    end
    # First line: tag + 8 spaces + up to 5 values
    first_n = min(5, length(values))
    parts = [@sprintf("%12.5E", Float64(v)) for v in values[1:first_n]]
    push!(lines, rpad(tag, 4) * "        " * join(parts, ""))
    # Continuation lines: 6 values each
    idx = first_n + 1
    while idx <= length(values)
        chunk_end = min(idx + 5, length(values))
        push!(lines, join([@sprintf("%12.5E", Float64(v)) for v in values[idx:chunk_end]], ""))
        idx = chunk_end + 1
    end
    join(lines, "\n")
end

# =========================================================================
# Main MATXS writer
# =========================================================================

"""
    write_matxs(io::IO, multigroup_data::MultiGroupXS;
                title="NJOY.jl", ivers=0, particle="n",
                data_type="nscat", material_name="mat1",
                mat_id=0, atom_mass=1.0)

Write a MATXS-format ASCII file from multigroup cross section data.

Produces a complete MATXS file with one particle, one data type, and one
material. Cross section vectors are written for each MT in the data.

# Arguments
- `io`: output IO stream
- `multigroup_data`: `MultiGroupXS` struct from groupr
- `title`: set hollerith identification (up to 72 chars)
- `ivers`: file version number
- `particle`: particle identifier (e.g. "n" for neutron)
- `data_type`: data type identifier (e.g. "nscat")
- `material_name`: material hollerith identifier
- `mat_id`: integer material identifier (ENDF MAT number)
- `atom_mass`: atomic weight ratio
"""
function write_matxs(io::IO, multigroup_data::MultiGroupXS;
                     title::AbstractString="NJOY.jl",
                     ivers::Integer=0,
                     particle::AbstractString="n",
                     data_type::AbstractString="nscat",
                     material_name::AbstractString="mat1",
                     mat_id::Integer=0,
                     atom_mass::Real=1.0)
    gb = multigroup_data.group_bounds
    ngroup = length(gb) - 1
    mt_list = multigroup_data.mt_list
    xs = multigroup_data.xs
    n_reactions = length(mt_list)

    npart = 1    # one particle type
    ntype = 1    # one data type
    nholl = 1    # one line of hollerith ID
    nmat = 1     # one material
    maxw = 5000  # max record size for sub-blocking
    nsubm = 1    # one submaterial

    # Compute file length: records count
    # file_id + file_control + hollerith + file_data + group_struct
    #   + mat_control + vec_control + vec_block = 8 records
    length_val = 8

    # --- Record 1: File identification ---
    # format: ' 0v ', hfile(a8), '*', huse(2*a8), '*', ivers(i6)
    hfile = rpad("matxs", 8)[1:8]
    huse1 = rpad("NJOY.jl", 8)[1:8]
    huse2 = rpad("", 8)[1:8]
    println(io, @sprintf(" 0v %s*%s%s*%6d", hfile, huse1, huse2, ivers))

    # --- Record 2: File control ---
    # format: ' 1d   ', npart, ntype, nholl, nmat, maxw, length
    println(io, @sprintf(" 1d   %6d%6d%6d%6d%6d%6d",
                         npart, ntype, nholl, nmat, maxw, length_val))

    # --- Record 3: Set hollerith identification ---
    padded_title = rpad(title, 72)[1:72]
    # Split into 9 a8 words
    hwords = [padded_title[(i-1)*8+1 : i*8] for i in 1:9]
    println(io, " 2d ")
    println(io, join(hwords, ""))

    # --- Record 4: File data ---
    # Hollerith part: hprt(npart), htype(ntype), hmatn(nmat)
    hprt = [rpad(particle, 8)[1:8]]
    htype = [rpad(data_type, 8)[1:8]]
    hmatn = [rpad(material_name, 8)[1:8]]
    all_holl = vcat(hprt, htype, hmatn)
    println(io, " 3d     " * join(all_holl, ""))
    # Integer part: ngrp(npart), jinp(ntype), joutp(ntype), nsubm(nmat), locm(nmat)
    ngrp_vals = Int[ngroup]
    jinp_vals = Int[1]       # input particle = particle 1
    joutp_vals = Int[1]      # output particle = particle 1
    nsubm_vals = Int[nsubm]
    locm_vals = Int[1]       # location of material 1
    all_ints = vcat(ngrp_vals, jinp_vals, joutp_vals, nsubm_vals, locm_vals)
    println(io, _matxs_int_line(all_ints))

    # --- Record 5: Group structure ---
    # gpb(ngroup) in descending energy order, then emin
    gpb = [gb[ngroup - g + 2] for g in 1:ngroup]
    emin = gb[1]
    all_bounds = vcat(gpb, [emin])
    println(io, _matxs_real_line(" 4d ", Float64.(all_bounds)))

    # --- Record 6: Material control ---
    # hmat(a8), amass, temp, sigz, itype, n1d, n2d, locs
    temp_val = 0.0
    sigz_val = 1.0e10  # infinite dilution
    n1d = n_reactions   # number of vectors
    n2d = 0             # no matrices (vectors only for simplicity)
    locs_val = 1
    println(io, @sprintf(" 5d %s%12.5E", rpad(material_name, 8)[1:8],
                         Float64(atom_mass)))
    println(io, @sprintf("%12.5E%12.5E%6d%6d%6d%6d",
                         temp_val, sigz_val, 1, n1d, n2d, locs_val))

    # --- Record 7: Vector control ---
    # hvps(n1d), nfg(n1d), nlg(n1d)
    mt_names = [_mt_to_matxs_name(mt) for mt in mt_list]
    println(io, _matxs_hollerith_line(" 6d ", mt_names))
    # nfg: first group in band (all start at 1)
    nfg_vals = ones(Int, n_reactions)
    # nlg: last group in band (all end at ngroup)
    nlg_vals = fill(ngroup, n_reactions)
    println(io, _matxs_int_line(vcat(nfg_vals, nlg_vals)))

    # --- Record 8: Vector block ---
    # All vector data concatenated: for each reaction, ngroup values
    all_vdata = Float64[]
    for r in 1:n_reactions
        for g in 1:ngroup
            push!(all_vdata, xs[g, r])
        end
    end
    println(io, _matxs_real_line(" 7d ", all_vdata))

    return nothing
end

# =========================================================================
# MATXS reaction name mapping
# =========================================================================

"""
    _mt_to_matxs_name(mt::Integer) -> String

Convert ENDF MT number to MATXS-style reaction name (up to 8 characters).
"""
function _mt_to_matxs_name(mt::Integer)
    mt_names = Dict{Int,String}(
        1   => "ntot",
        2   => "nelas",
        4   => "ninel",
        16  => "n2n",
        17  => "n3n",
        18  => "nfiss",
        102 => "ngamma",
        103 => "np",
        104 => "nd",
        105 => "nt",
        106 => "nhe3",
        107 => "nalpha",
        251 => "mubar",
        252 => "xi",
        253 => "gamma",
        452 => "nubar",
    )
    return get(mt_names, Int(mt), @sprintf("mt%d", mt))
end

# =========================================================================
# Binary MATXS writer (Fortran unformatted)
# =========================================================================

"""
    write_matxs_binary(io::IO, multigroup_data::MultiGroupXS;
                       title="NJOY.jl", ivers=0, particle="n",
                       data_type="nscat", material_name="mat1",
                       mat_id=0, atom_mass=1.0)

Write a MATXS-format binary (Fortran unformatted) file.
Uses the same record structure as write_matxs but in binary form.
"""
function write_matxs_binary(io::IO, multigroup_data::MultiGroupXS;
                            title::AbstractString="NJOY.jl",
                            ivers::Integer=0,
                            particle::AbstractString="n",
                            data_type::AbstractString="nscat",
                            material_name::AbstractString="mat1",
                            mat_id::Integer=0,
                            atom_mass::Real=1.0)
    gb = multigroup_data.group_bounds
    ngroup = length(gb) - 1
    mt_list = multigroup_data.mt_list
    xs = multigroup_data.xs
    n_reactions = length(mt_list)

    npart = 1
    ntype = 1
    nholl = 1
    nmat = 1
    maxw = 5000
    nsubm = 1
    length_val = 8

    # Record 1: File identification
    rec1 = _record_buf("matxs  ", "NJOY.jl ", "        ", Int32(ivers))
    _write_record(io, rec1)

    # Record 2: File control
    rec2 = _record_buf(Int32[npart, ntype, nholl, nmat, maxw, length_val])
    _write_record(io, rec2)

    # Record 3: Set hollerith identification
    padded = rpad(title, 72)[1:72]
    hwords = [padded[(i-1)*8+1:i*8] for i in 1:9]
    rec3 = _record_buf(hwords...)
    _write_record(io, rec3)

    # Record 4: File data
    buf4 = IOBuffer()
    write(buf4, codeunits(_pad8(particle)))
    write(buf4, codeunits(_pad8(data_type)))
    write(buf4, codeunits(_pad8(material_name)))
    for v in Int32[ngroup, 1, 1, nsubm, 1]
        write(buf4, htol(v))
    end
    _write_record(io, take!(buf4))

    # Record 5: Group structure
    buf5 = IOBuffer()
    for g in 1:ngroup
        write(buf5, htol(Float32(gb[ngroup - g + 2])))
    end
    write(buf5, htol(Float32(gb[1])))
    _write_record(io, take!(buf5))

    # Record 6: Material control
    buf6 = IOBuffer()
    write(buf6, codeunits(_pad8(material_name)))
    write(buf6, htol(Float32(atom_mass)))
    write(buf6, htol(Float32(0.0)))      # temp
    write(buf6, htol(Float32(1.0e10)))   # sigz
    write(buf6, htol(Int32(1)))          # itype
    write(buf6, htol(Int32(n_reactions)))  # n1d
    write(buf6, htol(Int32(0)))          # n2d
    write(buf6, htol(Int32(1)))          # locs
    _write_record(io, take!(buf6))

    # Record 7: Vector control
    buf7 = IOBuffer()
    for mt in mt_list
        write(buf7, codeunits(_pad8(_mt_to_matxs_name(mt))))
    end
    for _ in mt_list
        write(buf7, htol(Int32(1)))       # nfg
    end
    for _ in mt_list
        write(buf7, htol(Int32(ngroup)))  # nlg
    end
    _write_record(io, take!(buf7))

    # Record 8: Vector block
    buf8 = IOBuffer()
    for r in 1:n_reactions
        for g in 1:ngroup
            write(buf8, htol(Float32(xs[g, r])))
        end
    end
    _write_record(io, take!(buf8))

    return nothing
end
