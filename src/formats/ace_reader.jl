# Type-1 (ASCII) ACE file reader.
#
# Mirrors Fortran acefix's itype=1 read path (acefc.f90:13992-14013, mcnpx=0):
#
#   read(nin,'(a10,f12.6,1x,1p,e11.4,1x,a10/a70,a10)') hz,aw0,tz,hd,hko,hm
#   read(nin,'(4(i7,f11.0))') (izo(i),awo(i),i=1,16)
#   read(nin,'(8i9)') len2,izaid,nes,...,ploct            ! NXS(16)+JXS(32)
#   n=(len2+3)/4 ; do i=1,n ; read(nin,'(4e20.0)') xss(...) ; enddo
#
# The companion writer is src/formats/ace_writer.jl (write_ace / write_ace_table);
# this reader inverts it, parsing a file back into the flat `ACETable` that
# `_acer_aplots` consumes. Because the ACE body we read is the byte-identical
# tape34 copy, every value here is provably correct — the reader is a pure
# inverse of an already-verified serialization.

using Printf

"""
    read_ace_ascii(lines::Vector{String}) -> ACETable

Parse a Type-1 (ASCII) ACE table from its raw text lines into the flat
`ACETable` (header, 16 (IZ,AW) pairs, NXS(16), JXS(32), XSS).

Field widths are fixed-column (Fortran formatted I/O), not whitespace-delimited:
header line 1 is `a10,f12.6,1x,e11.4,1x,a10`; the 16 (IZ,AW) pairs are
`4(i7,f11.0)` (4 pairs per line, 4 lines); the 48 control integers are `8i9`
(8 per line, 6 lines); the XSS body is `4e20.0` (4 fields of width 20 per line).

Ref: njoy-reference/src/acefc.f90:13992-14013 (acefix, itype=1, mcnpx=0).
"""
function read_ace_ascii(lines::Vector{String})
    length(lines) >= 12 || error("read_ace_ascii: ACE file has only $(length(lines)) lines; need >= 12 for the header block")

    # --- Header line 1: hz(a10) aw0(f12.6) 1x tz(e11.4) 1x hd(a10) ---
    l1   = rpad(lines[1], 80)
    hz   = l1[1:10]
    aw0  = parse(Float64, strip(l1[11:22]))
    tz   = parse(Float64, strip(l1[24:34]))
    hd   = l1[36:45]

    # --- Header line 2: hko(a70) hm(a10) ---
    l2      = rpad(lines[2], 80)
    comment = l2[1:70]
    hm      = l2[71:80]

    # --- Lines 3-6: 16 (IZ,AW) pairs, 4(i7,f11.0) ---
    izs = zeros(Int32, 16)
    aws = zeros(Float64, 16)
    pair = 0
    for row in 1:4
        ln = rpad(lines[2 + row], 72)
        for col in 1:4
            pair += 1
            off = (col - 1) * 18
            izs[pair] = parse(Int32, strip(ln[off+1 : off+7]))
            aws[pair] = parse(Float64, strip(ln[off+8 : off+18]))
        end
    end

    # --- Lines 7-12: NXS(16) then JXS(32) = 48 integers, 8i9 ---
    ints = zeros(Int32, 48)
    idx = 0
    for row in 1:6
        ln = rpad(lines[6 + row], 72)
        for col in 1:8
            idx += 1
            ints[idx] = parse(Int32, strip(ln[(col-1)*9 + 1 : col*9]))
        end
    end
    nxs = ntuple(i -> ints[i], 16)
    jxs = ntuple(i -> ints[16 + i], 32)

    # --- XSS body: len2 = NXS(1), packed 4e20.0 ---
    len2 = Int(nxs[NXS_LEN2])
    xss  = Vector{Float64}(undef, len2)
    nlines = (len2 + 3) ÷ 4
    got = 0
    for row in 1:nlines
        ln = lines[12 + row]
        for col in 1:4
            got >= len2 && break
            off = (col - 1) * 20
            field = ln[off+1 : min(off+20, length(ln))]
            got += 1
            xss[got] = _parse_ace_real(field)
        end
    end
    got == len2 || error("read_ace_ascii: read $got XSS values, expected len2=$len2")

    ACETable(strip(hz), aw0, tz, strip(hd), comment, hm,
             ntuple(i -> izs[i], 16), ntuple(i -> aws[i], 16),
             nxs, jxs, xss)
end

"""
    _parse_ace_real(s::AbstractString) -> Float64

Parse one `e20.0` XSS field. The writer (acecm.f90 typen) drops the `E` for
3-digit exponents (e.g. `1.15510000000-117`); Fortran list-directed read
accepts that, but Julia `parse` does not, so re-insert the `E` before the
final signed exponent when the `E` is missing.
"""
function _parse_ace_real(s::AbstractString)
    t = strip(s)
    isempty(t) && return 0.0
    if occursin('E', t) || occursin('e', t)
        return parse(Float64, t)
    end
    # No exponent letter: a sign in the interior (after the first char) marks
    # the start of a 3-digit exponent, e.g. "1.15510000000-117".
    m = findlast(c -> c == '+' || c == '-', t)
    if m !== nothing && m > 1
        return parse(Float64, t[1:m-1] * "E" * t[m:end])
    end
    parse(Float64, t)
end
