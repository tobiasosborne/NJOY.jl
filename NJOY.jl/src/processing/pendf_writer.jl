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
                     description::Vector{String} = String[])
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

    # MF3 for each reaction
    for (col, mt) in enumerate(material.mt_list)
        _write_mf3_from_matrix(io, mat, mt, material.energies,
                               material.cross_sections, col, ns)
    end
    _write_fend(io, mat)

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
                     tempr::Float64 = 0.0)
    mf2 = result.mf2
    actual_mat = mat > 0 ? Int32(mat) : Int32(max(1, round(Int, mf2.ZA / 10)))
    ns = Ref(1)

    # Determine reactions
    reactions = _collect_reactions(result)

    # TPID
    _write_tpid_line(io, label, Int(actual_mat))

    # MF1/MT451
    _write_legacy_mf1(io, mf2, actual_mat, err, tempr, length(result.energies), reactions, ns)

    # MF2/MT151
    _write_legacy_mf2(io, mf2, actual_mat, ns)

    # MF3 sections
    _write_legacy_mf3(io, result, actual_mat, reactions, ns)

    # Terminators
    _write_record_sep(io, 0, 0, 0)
    _write_record_sep(io, -1, 0, 0)
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
# Low-level ENDF formatting
# ==========================================================================

function _write_tpid_line(io::IO, text::AbstractString, mat::Int)
    t = rpad(text, 66)[1:66]
    @printf(io, "%s%4d%2d%3d%5d\n", t, 1, 0, 0, 0)
end

function _write_record_sep(io::IO, mat::Int, mf::Int, mt::Int)
    blanks = repeat(" ", 66)
    @printf(io, "%s%4d%2d%3d%5d\n", blanks, mat, mf, mt, 99999)
end

function _write_fend(io::IO, mat::Int)
    _write_record_sep(io, mat, 3, 0)
    _write_record_sep(io, mat, 0, 0)
end

function _write_cont_line(io::IO, c1, c2, l1, l2, n1, n2,
                          mat::Int, mf::Int, mt::Int, ns::Ref{Int})
    s1 = format_endf_float(Float64(c1))
    s2 = format_endf_float(Float64(c2))
    @printf(io, "%s%s%11d%11d%11d%11d%4d%2d%3d%5d\n",
            s1, s2, l1, l2, n1, n2, mat, mf, mt, ns[])
    ns[] += 1
end

function _write_data_values(io::IO, vals, mat::Int, mf::Int, mt::Int, ns::Ref{Int})
    n = length(vals)
    i = 1
    while i <= n
        for j in 1:6
            if i <= n
                write(io, format_endf_float(Float64(vals[i])))
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
                                 ns::Ref{Int})
    n = length(energies)
    n >= 2 || return
    ns[] = 1

    # HEAD: ZA=0, AWR=0, 0, LR=99, 0, 0
    _write_cont_line(io, 0.0, 0.0, 0, 99, 0, 0, mat, 3, mt, ns)

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
    _write_data_values(io, data, mat, 3, mt, ns)

    # SEND
    _write_record_sep(io, mat, 3, 0)
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
    reactions = Tuple{Int32, String}[]
    push!(reactions, (Int32(1), "total"))
    push!(reactions, (Int32(2), "elastic"))
    if any(x -> x > 0.0, result.fission)
        push!(reactions, (Int32(18), "fission"))
    end
    push!(reactions, (Int32(102), "capture"))
    for sec in result.mf3_sections
        mt = Int(sec.mt)
        if mt != 1 && mt != 2 && mt != 18 && mt != 19 && mt != 102
            push!(reactions, (Int32(mt), "MT$mt"))
        end
    end
    return reactions
end

function _write_legacy_mf1(io::IO, mf2::MF2Data, mat::Int32, err, tempr,
                            n_pts, reactions, ns::Ref{Int})
    ns[] = 1
    nxc = 2 + length(reactions)

    _write_cont_line(io, mf2.ZA, mf2.AWR, 2, 0, 0, 0,
                     Int(mat), 1, 451, ns)
    _write_cont_line(io, tempr, err, 0, 0, 0, nxc,
                     Int(mat), 1, 451, ns)

    nc_per = 3 + cld(n_pts, 3)
    _write_cont_line(io, 0.0, 0.0, 1, 451, 3, 0, Int(mat), 1, 451, ns)
    _write_cont_line(io, 0.0, 0.0, 2, 151, 4, 0, Int(mat), 1, 451, ns)
    for (mt, _) in reactions
        _write_cont_line(io, 0.0, 0.0, 3, Int(mt), nc_per, 0,
                         Int(mat), 1, 451, ns)
    end

    _write_record_sep(io, Int(mat), 1, 0)
    _write_record_sep(io, Int(mat), 0, 0)
end

function _write_legacy_mf2(io::IO, mf2::MF2Data, mat::Int32, ns::Ref{Int})
    ns[] = 1
    nis = length(mf2.isotopes)
    _write_cont_line(io, mf2.ZA, mf2.AWR, 0, 0, nis, 0, Int(mat), 2, 151, ns)

    for iso in mf2.isotopes
        _write_cont_line(io, iso.ZAI, iso.ABN, 0, 0, 1, 0, Int(mat), 2, 151, ns)
        el = length(iso.ranges) > 0 ? iso.ranges[1].EL : 1.0e-5
        eh = length(iso.ranges) > 0 ? iso.ranges[end].EH : 2.0e7
        _write_cont_line(io, el, eh, 0, 0, 0, 0, Int(mat), 2, 151, ns)
        spi = 1.0; ap = 0.0
        if length(iso.ranges) > 0
            p = iso.ranges[1].parameters
            if hasproperty(p, :SPI); spi = p.SPI; end
            if hasproperty(p, :AP); ap = p.AP; end
        end
        _write_cont_line(io, spi, ap, 0, 0, 0, 0, Int(mat), 2, 151, ns)
    end

    _write_record_sep(io, Int(mat), 2, 0)
    _write_record_sep(io, Int(mat), 0, 0)
end

function _write_legacy_mf3(io::IO, result, mat::Int32, reactions, ns::Ref{Int})
    energies = result.energies
    n = length(energies)

    for (mt, _) in reactions
        xs_data = _get_legacy_xs(result, Int(mt))
        xs_data === nothing && continue

        ns[] = 1

        _write_cont_line(io, result.mf2.ZA, result.mf2.AWR, 0, 99, 0, 0,
                         Int(mat), 3, Int(mt), ns)
        _write_cont_line(io, 0.0, 0.0, 0, 0, 1, n,
                         Int(mat), 3, Int(mt), ns)

        vals = Float64[Float64(n), 2.0]
        _write_data_values(io, vals, Int(mat), 3, Int(mt), ns)

        data = Float64[]
        sizehint!(data, 2 * n)
        for i in 1:n
            push!(data, energies[i])
            push!(data, xs_data[i])
        end
        _write_data_values(io, data, Int(mat), 3, Int(mt), ns)

        _write_record_sep(io, Int(mat), 3, 0)
    end

    _write_record_sep(io, Int(mat), 0, 0)
end

function _get_legacy_xs(result, mt::Int)
    mt == 1 && return result.total
    mt == 2 && return result.elastic
    (mt == 18 || mt == 19) && return result.fission
    mt == 102 && return result.capture
    for sec in result.mf3_sections
        Int(sec.mt) == mt && return [interpolate(sec.tab, e) for e in result.energies]
    end
    return nothing
end

function _pendf_parse_int(s::AbstractString)
    t = strip(s)
    isempty(t) && return 0
    return parse(Int, t)
end
