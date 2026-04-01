# types.jl -- Core types for the NJOY module orchestration layer
#
# TapeManager: maps Fortran logical unit numbers to file paths
# PENDFTape: structured representation of a PENDF tape file
#
# Design: matches Fortran NJOY architecture — modules communicate via tape files,
# no shared mutable state between modules.

# =========================================================================
# Tape Manager
# =========================================================================

"""
    TapeManager

Maps Fortran-style logical unit numbers to file paths.

In Fortran NJOY, modules reference tapes by integer unit numbers
(e.g., 20 for the ENDF input, -21 for binary PENDF). The sign indicates
format (positive=ASCII, negative=binary) which is irrelevant in Julia.

The TapeManager resolves unit numbers to actual file paths, generating
default paths in `work_dir` for tapes not explicitly registered.
"""
struct TapeManager
    unit_to_path::Dict{Int, String}
    work_dir::String
end

"""
    resolve(tm::TapeManager, unit::Int) -> String

Resolve a tape unit number to a file path. Uses the absolute value of
the unit (sign is the Fortran binary/ASCII convention, irrelevant here).
Returns a registered path if available, otherwise `work_dir/tapeNN`.
"""
function resolve(tm::TapeManager, unit::Int)::String
    au = abs(unit)
    au == 0 && error("tape unit 0 is invalid")
    get(tm.unit_to_path, au, joinpath(tm.work_dir, "tape$au"))
end

"""
    register!(tm::TapeManager, unit::Int, path::String)

Register a file path for a tape unit number.
"""
function register!(tm::TapeManager, unit::Int, path::String)
    tm.unit_to_path[abs(unit)] = path
end

# =========================================================================
# PENDF Tape Representation
# =========================================================================

"""
    PENDFSection

One MF/MT section from a PENDF tape, stored as raw ENDF-format lines.
Lines include columns 1-66 (data) + 67-75 (MAT/MF/MT) but NOT the
sequence number in columns 76-80 (regenerated on write).
"""
struct PENDFSection
    mf::Int
    mt::Int
    lines::Vector{String}   # raw lines, each ≥ 75 chars
end

"""
    PENDFMaterial

One material (MAT) from a PENDF tape. The MF1/MT451 header is stored
separately because modules need to update the directory when adding
or modifying sections.
"""
struct PENDFMaterial
    mat::Int
    mf1_lines::Vector{String}       # MF1/MT451 header + directory
    sections::Vector{PENDFSection}   # all other sections in tape order
end

"""
    PENDFTape

Structured representation of a complete PENDF tape. Contains a tape ID
line and one or more materials.

The Fortran NJOY modules stream input PENDF → output PENDF, copying
unmodified sections and modifying/adding their own. In Julia, we read
the whole tape into this struct, modify it, and write it back.
"""
struct PENDFTape
    tpid::String                         # tape identification line
    materials::Vector{PENDFMaterial}
end
