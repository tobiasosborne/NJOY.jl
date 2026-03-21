# RECONR data types and MF3 reader
#
# Data structures for the RECONR pipeline:
#   MF3Section    -- a single MF3 cross section section from ENDF
#   ENDFMaterial  -- all ENDF data needed for RECONR processing
#   PointwiseMaterial -- result of RECONR processing

# ==========================================================================
# Data types
# ==========================================================================

"""
    MF3Section

A single MF3/MF13/MF23 cross section section from ENDF.
Stores the reaction type and the tabulated cross section data.
"""
struct MF3Section
    mt::Int32
    QM::Float64
    QI::Float64
    tab::TabulatedFunction
end

"""
    ENDFMaterial

All ENDF data needed for RECONR processing of one material.
"""
struct ENDFMaterial
    mat::Int32                              # MAT number
    awr::Float64                            # atomic weight ratio
    mf2::MF2Data                            # File 2 resonance parameters
    mf3_sections::Vector{MF3Section}        # File 3 cross section sections
    redundant_mts::Set{Int}                 # MTs that are sums (1, 3, 4, 18, 101, 27)
end

"""
    PointwiseMaterial

Result of RECONR processing: pointwise cross sections on a linearized grid.
"""
struct PointwiseMaterial
    mat::Int32
    energies::Vector{Float64}
    cross_sections::Matrix{Float64}  # (n_energies, n_reactions)
    mt_list::Vector{Int}             # which MT numbers are stored in each column
end

# Column indices for the primary 4 resonance channels
const _COL_TOTAL   = 1
const _COL_ELASTIC = 2
const _COL_FISSION = 3
const _COL_CAPTURE = 4

# ==========================================================================
# MF3 reader
# ==========================================================================

"""
    read_mf3_sections(io::IO, mat::Integer) -> Vector{MF3Section}

Read all MF3 cross section sections for material `mat`.
"""
function read_mf3_sections(io::IO, mat::Integer)
    sections = MF3Section[]
    seekstart(io)

    while !eof(io)
        pos = position(io)
        line = readline(io)
        p = rpad(line, 80)
        mf = _parse_int(p[71:72])
        mt = _parse_int(p[73:75])
        mat_line = _parse_int(p[67:70])

        mat_line != mat && continue

        if mf == 3 && mt > 0
            seek(io, pos)
            try
                head = read_cont(io)
                tab1 = read_tab1(io)
                tf = TabulatedFunction(tab1)
                push!(sections, MF3Section(Int32(mt), head.C1, tab1.C2, tf))
                # Skip to SEND
                _skip_to_send(io)
            catch e
                @warn "read_mf3_sections: skipping MF3/MT=$mt due to parse error" exception=(e, catch_backtrace())
                _skip_to_send(io)
            end
        end
    end
    return sections
end

function _skip_to_send(io::IO)
    while !eof(io)
        line = readline(io)
        p = rpad(line, 80)
        mt = _parse_int(p[73:75])
        mt == 0 && break
    end
end

function _detect_mat(io::IO, requested::Integer)
    requested > 0 && return Int32(requested)
    seekstart(io)
    while !eof(io)
        line = readline(io)
        p = rpad(line, 80)
        mf_val = _parse_int(p[71:72])
        mat_val = _parse_int(p[67:70])
        if mat_val > 0 && mf_val > 0
            return Int32(mat_val)
        end
    end
    error("reconstruct: could not detect MAT number from file")
end

# _parse_int is already defined in endf/io.jl -- no need to redefine here
