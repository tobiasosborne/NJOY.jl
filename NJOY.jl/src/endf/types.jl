# ENDF-6 record types and core data structures
#
# Based on the ENDF-6 format specification and NJOY2016's endf.f90.
# An ENDF record line is 80 columns:
#   columns 1-66: data (6 fields of 11 characters each)
#   columns 67-70: MAT number (I4)
#   columns 71-72: MF number (I2)
#   columns 73-75: MT number (I3)
#   columns 76-80: line sequence number (I5)

"""
    InterpolationLaw

ENDF interpolation law codes (INT field in TAB1/TAB2 records).
"""
@enum InterpolationLaw::Int32 begin
    Histogram   = 1   # y is constant (histogram)
    LinLin      = 2   # y is linear in x
    LinLog      = 3   # y is linear in ln(x)
    LogLin      = 4   # ln(y) is linear in x
    LogLog      = 5   # ln(y) is linear in ln(x)
    CoulombPen  = 6   # Coulomb penetrability law (charged particles)
end

"""
    MaterialId

ENDF material/file/section identifier triple.
"""
struct MaterialId
    mat::Int32   # material number (MAT)
    mf::Int32    # file number (MF)
    mt::Int32    # section number (MT)
end

"""
    ContRecord

ENDF CONT (control) record: 2 floats + 4 integers.

Fields match NJOY convention: C1, C2, L1, L2, N1, N2.
"""
struct ContRecord
    C1::Float64
    C2::Float64
    L1::Int32
    L2::Int32
    N1::Int32
    N2::Int32
    id::MaterialId
end

"""
    InterpolationTable

Interpolation table: a set of (NBT, INT) pairs specifying which
interpolation law applies in each region of a tabulated function.
- `nbt`: breakpoint indices (last point index for each region)
- `law`: interpolation law for each region
"""
struct InterpolationTable
    nbt::Vector{Int32}
    law::Vector{InterpolationLaw}

    function InterpolationTable(nbt::AbstractVector{<:Integer},
                                law::AbstractVector{<:Integer})
        length(nbt) == length(law) ||
            throw(ArgumentError("nbt and law must have same length"))
        new(Int32.(nbt), InterpolationLaw.(law))
    end

    function InterpolationTable(nbt::AbstractVector{<:Integer},
                                law::AbstractVector{InterpolationLaw})
        length(nbt) == length(law) ||
            throw(ArgumentError("nbt and law must have same length"))
        new(Int32.(nbt), collect(law))
    end
end

Base.length(it::InterpolationTable) = length(it.nbt)

"""
    ListRecord

ENDF LIST record: CONT header followed by N1 data values.
"""
struct ListRecord
    C1::Float64
    C2::Float64
    L1::Int32
    L2::Int32
    N1::Int32     # number of data items (NPL)
    N2::Int32
    data::Vector{Float64}
    id::MaterialId
end

"""
    Tab1Record

ENDF TAB1 record: CONT header + interpolation table + (x,y) pairs.
- N1 = NR (number of interpolation ranges)
- N2 = NP (number of data points)
"""
struct Tab1Record
    C1::Float64
    C2::Float64
    L1::Int32
    L2::Int32
    interp::InterpolationTable
    x::Vector{Float64}
    y::Vector{Float64}
    id::MaterialId

    function Tab1Record(C1, C2, L1, L2, interp, x, y, id)
        length(x) == length(y) ||
            throw(ArgumentError("x and y must have same length"))
        new(Float64(C1), Float64(C2), Int32(L1), Int32(L2),
            interp, Float64.(x), Float64.(y), id)
    end
end

"""
    Tab2Record

ENDF TAB2 record: CONT header + interpolation table.
Used as a header for a sequence of TAB1 records (e.g., angular distributions).
- N1 = NR (number of interpolation ranges)
- N2 = NZ (number of subsequent records)
"""
struct Tab2Record
    C1::Float64
    C2::Float64
    L1::Int32
    L2::Int32
    interp::InterpolationTable
    NZ::Int32
    id::MaterialId
end

"""
    TabulatedFunction

A tabulated function with interpolation metadata, suitable for
evaluation and integration. Wraps the data from a TAB1 record.
"""
struct TabulatedFunction
    interp::InterpolationTable
    x::Vector{Float64}
    y::Vector{Float64}

    function TabulatedFunction(interp::InterpolationTable,
                               x::AbstractVector{<:Real},
                               y::AbstractVector{<:Real})
        n = length(x)
        n == length(y) ||
            throw(ArgumentError("x and y must have same length"))
        n >= 1 || throw(ArgumentError("must have at least one point"))
        new(interp, Float64.(x), Float64.(y))
    end
end

function TabulatedFunction(rec::Tab1Record)
    TabulatedFunction(rec.interp, rec.x, rec.y)
end

Base.length(tf::TabulatedFunction) = length(tf.x)
