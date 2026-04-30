# covr_press.jl -- "boxer" compressed-format library output for covr.
#
# Ports covr.f90:1991-2247 (`press`, `setfor`).
#
# Boxer format compresses a (nrow × ncol) matrix into two arrays:
#   xval  -- distinct data values, formatted via the `nvf` format code.
#   icon  -- control codes:  positive = "carry-down N" repeats from row above;
#                            negative = "repeated-value N" copies of the
#                                       most recent xval entry.
# The (nrow,ncol) field accepts ncol=0 as the "symmetric" sentinel (only
# the upper triangle is stored, filled in by the reader).
#
# Output format (covr.f90:2199-2207):
#   first page :   '(i1,1x,a12,1x,a21,2(i5,i4),2(i4,i3),3i4)'
#   subsequent :   '(i1,1x,34("-"),         2(i5,i4),2(i4,i3),3i4)'
#   xval line  :   per `setfor(ift(nvf))` table — defaults to "(1p8e10.3)"
#   icon line  :   per `setfor(ift(ncf))`       — defaults to "(20i4)" / etc
#
# `itype` field meanings (covr.f90:2017-2020):
#   0 = energy-group boundaries
#   1 = cross sections
#   2 = relative standard deviations
#   3 = relative covariance matrix
#   4 = correlation matrix

using Printf

# Format-table identical to covr.f90:2230-2234. Indices 1..14 → number of
# values per line / per-value width.  setfor errors out if nvf ∉ [7,14] or
# ncf ∉ [1,6] (covr.f90:2236-2241).
const _COVR_PRESS_FMT = NamedTuple[
    (count=80, width=1,  kind=:int),
    (count=40, width=2,  kind=:int),
    (count=26, width=3,  kind=:int),
    (count=20, width=4,  kind=:int),
    (count=16, width=5,  kind=:int),
    (count=13, width=6,  kind=:int),
    (count=11, width=7,  kind=:f4),
    (count=10, width=8,  kind=:f5),
    (count=8,  width=9,  kind=:e2),
    (count=8,  width=10, kind=:e3),
    (count=7,  width=11, kind=:e4),
    (count=6,  width=12, kind=:e5),
    (count=6,  width=13, kind=:e6),
    (count=5,  width=14, kind=:e7),
]

"""
    covr_setfor(nvf::Integer, ncf::Integer) -> (val_fmt, ctl_fmt)

Return the per-line (count, width, kind) tuples for value (xval) and
control (icon) arrays. Mirrors covr.f90:2220-2247. Errors on out-of-range.
"""
function covr_setfor(nvf::Integer, ncf::Integer)
    (7 <= nvf <= 14 && 1 <= ncf <= 6) ||
        error("covr_setfor: nvf=$nvf or ncf=$ncf illegal (covr.f90:2236-2241)")
    (_COVR_PRESS_FMT[nvf], _COVR_PRESS_FMT[ncf])
end

# Format one float/int per its kind+width spec.
function _press_fmt(v::Real, fmt::NamedTuple)
    w = fmt.width
    if fmt.kind === :int
        lpad(string(Int(v)), w)
    elseif fmt.kind === :f4
        @sprintf("%*.*f", w, 4, Float64(v))
    elseif fmt.kind === :f5
        @sprintf("%*.*f", w, 5, Float64(v))
    elseif fmt.kind === :e2
        @sprintf("%*.*E", w, 2, Float64(v))
    elseif fmt.kind === :e3
        @sprintf("%*.*E", w, 3, Float64(v))
    elseif fmt.kind === :e4
        @sprintf("%*.*E", w, 4, Float64(v))
    elseif fmt.kind === :e5
        @sprintf("%*.*E", w, 5, Float64(v))
    elseif fmt.kind === :e6
        @sprintf("%*.*E", w, 6, Float64(v))
    elseif fmt.kind === :e7
        @sprintf("%*.*E", w, 7, Float64(v))
    end
end

"""
    BoxerWriter

Stateful boxer-format writer. One BoxerWriter per (matype, hlibid, hdescr)
library; calling `boxer_press!` for each (mat, mt, mat1, mt1, type, matrix)
emits one page (or paginated set, when xval/icon overflow nvmax/ncmax).
"""
mutable struct BoxerWriter
    io::IO
    matype::Int
    hlibid::String
    hdescr::String
    first::Bool   # first page sets the header; subsequent pages use the dashes form
end

BoxerWriter(io::IO; matype::Integer, hlibid::AbstractString, hdescr::AbstractString) =
    BoxerWriter(io, Int(matype),
                rpad(String(hlibid), 12)[1:12],
                rpad(String(hdescr), 21)[1:21],
                true)

const _NVMAX = 880
const _NCMAX = 900
const _PRESS_EPS = 1.0e-10   # covr.f90:2040

"""
    boxer_press!(w, mat, mt, mat1, mt1, itype, xa; nrow, ncol, nvf, ncf)

Compress a 2-D matrix `xa` (size `nrow × ncol`) into the boxer format and
write to `w.io`. `ncol == 0` activates symmetric mode (only upper triangle
stored). Mirrors covr.f90:1991-2218.

`itype`, `nvf`, `ncf` follow Fortran semantics (see file header).

When the (xval, icon) arrays would overflow `_NVMAX`/`_NCMAX`, the page
flushes and the loop continues on a new page. Paginated headers use the
"dashes" form per covr.f90:2203-2205.
"""
function boxer_press!(w::BoxerWriter,
                      mat::Integer, mt::Integer, mat1::Integer, mt1::Integer,
                      itype::Integer, xa::Matrix{Float64};
                      nrow::Integer, ncol::Integer,
                      nvf::Integer=10, ncf::Integer=3)::Nothing
    val_fmt, ctl_fmt = covr_setfor(nvf, ncf)
    ndig = nvf < 9 ? nvf - 3 : nvf - 6
    lvmax = 10^(ncf - 1) - 1
    lcmax = 10^ncf - 1

    nsym = ncol == 0 ? 1 : 0
    nsym == 1 && (ncol = nrow)

    # Sigfig-round the entire matrix once (covr.f90:2071-2081 for the
    # rectangular branch; symmetric branch sigfigs row-by-row in the original
    # but the result is identical after one bulk pass).
    xa_round = similar(xa, Float64)
    @inbounds for j in 1:size(xa, 2), i in 1:size(xa, 1)
        xa_round[i, j] = round_sigfig(xa[i, j], ndig, 0)
    end

    istart = 1
    while true
        xval = Float64[]
        icon = Int[]
        lv = 0
        lc = 0
        ivopt = true   # repeated-value still active?
        icopt = true   # carry-down still active?
        icd = 0
        xvtest = xa_round[istart, 1]
        nrowm = 0

        i = istart - 1
        while i < nrow
            i += 1
            jlow = (nsym == 1) ? i : 1

            # carry-down test reference (xa[i-1, j], or 0 on first row).
            for j in jlow:ncol
                v = xa_round[i, j]
                # Repeated-value test.
                if ivopt
                    test = _PRESS_EPS * abs(xvtest)
                    if abs(v - xvtest) <= test && lv != lvmax
                        lv += 1; icd = 0
                    else
                        ivopt = false
                    end
                end
                # Carry-down test.
                if icopt
                    xctest = i > istart ? xa_round[i - 1, j] : 0.0
                    test = _PRESS_EPS * abs(xctest)
                    if abs(v - xctest) <= test && lc != lcmax
                        lc += 1; icd = 1
                    else
                        icopt = false
                    end
                end
                if !icopt && !ivopt
                    # End of pattern run — store + reinit.
                    push!(icon, icd != 0 ? lc : -lv)
                    if icd == 0
                        push!(xval, xvtest)
                    end
                    lv = 0; lc = 0; ivopt = true; icopt = true
                    xvtest = v
                    # Re-test this same v under the new pattern start.
                    test = _PRESS_EPS * abs(xvtest)
                    if abs(v - xvtest) <= test && lv != lvmax
                        lv += 1; icd = 0
                    else
                        ivopt = false
                    end
                    if icopt
                        xctest = i > istart ? xa_round[i - 1, j] : 0.0
                        test = _PRESS_EPS * abs(xctest)
                        if abs(v - xctest) <= test && lc != lcmax
                            lc += 1; icd = 1
                        else
                            icopt = false
                        end
                    end
                    if !icopt && !ivopt
                        # double-fail: still need to flush
                        push!(icon, icd != 0 ? lc : -lv)
                        icd == 0 && push!(xval, xvtest)
                        lv = 0; lc = 0; ivopt = true; icopt = true
                        xvtest = v
                    end
                end
            end
            # Page-overflow check (covr.f90:2180-2183).
            leng = nsym == 1 ? (ncol + 1 - i) : (ncol + 1)
            if length(xval) + leng > _NVMAX || length(icon) + leng > _NCMAX
                break
            end
        end

        # Flush trailing pattern (covr.f90:2185-2192).
        push!(icon, icd != 0 ? lc : -lv)
        icd == 0 && push!(xval, xvtest)
        nrowm = nrow - i

        ncol_out = ncol
        nsym == 1 && (ncol_out = 0)

        _emit_boxer_header!(w, mat, mt, mat1, mt1, itype,
                            length(xval), nvf, length(icon), ncf,
                            nrowm, nrow, ncol_out)
        !isempty(xval) && _emit_boxer_data!(w.io, xval, val_fmt)
        !isempty(icon) && _emit_boxer_data!(w.io, icon, ctl_fmt)

        nrowm == 0 && break
        istart = i + 1
    end
    nothing
end

# Header emission; first page includes hlibid/hdescr, subsequent use 34 dashes.
function _emit_boxer_header!(w::BoxerWriter,
                              mat::Integer, mt::Integer,
                              mat1::Integer, mt1::Integer,
                              itype::Integer, nval::Integer, nvf::Integer,
                              ncon::Integer, ncf::Integer,
                              nrowm::Integer, nrow::Integer, ncol::Integer)
    if w.first
        @printf(w.io,
                "%1d %s %s%5d%4d%5d%4d%4d%3d%4d%3d%4d%4d%4d\n",
                itype, w.hlibid, w.hdescr,
                Int(mat), Int(mt), Int(mat1), Int(mt1),
                Int(nval), Int(nvf), Int(ncon), Int(ncf),
                Int(nrowm), Int(nrow), Int(ncol))
        w.first = false
    else
        @printf(w.io,
                "%1d %s%5d%4d%5d%4d%4d%3d%4d%3d%4d%4d%4d\n",
                itype, "----------------------------------",
                Int(mat), Int(mt), Int(mat1), Int(mt1),
                Int(nval), Int(nvf), Int(ncon), Int(ncf),
                Int(nrowm), Int(nrow), Int(ncol))
    end
end

# Data emission: pack `vals` per (count, width, kind) format spec.
function _emit_boxer_data!(io::IO, vals::AbstractVector, fmt::NamedTuple)
    n = length(vals)
    cnt = fmt.count
    k = 0
    while k < n
        bs = min(cnt, n - k)
        for j in 1:bs
            print(io, _press_fmt(vals[k + j], fmt))
        end
        println(io)
        k += bs
    end
end
