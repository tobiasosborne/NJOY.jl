# ACE Type 1 (ASCII) writer for MCNP continuous-energy neutron tables
# Matches NJOY2016 aceout()/typen()/change() format exactly.
# aceout -> write_ace_table; typen -> _ace_fmt_int/_ace_fmt_real

using Printf

# =========================================================================
# Low-level Type 1 format helpers (matching NJOY2016 typen subroutine)
# =========================================================================

"""Format an integer for XSS data block: right-justified in 20 chars (i20)."""
_ace_fmt_int(val::Integer) = lpad(string(Int(val)), 20)

"""Format a real for XSS data block: scientific in 20 chars (1pe20.11)."""
_ace_fmt_real(val::Real) = @sprintf("%20.11E", Float64(val))

# Also keep Proposer-A naming for compatibility
_format_xss_value(x::Float64) = _ace_fmt_real(x)
_format_xss_integer(x::Int) = _ace_fmt_int(x)

# =========================================================================
# Type 1 ACE file writer (Proposer-B: from ACENeutronTable)
# =========================================================================

"""
    write_ace_table(io::IO, table::ACENeutronTable; format=:ascii)

Write an ACENeutronTable as a Type 1 (ASCII) ACE file.
Output format matches NJOY2016 aceout() exactly for MCNP compatibility.
"""
function write_ace_table(io::IO, table::ACENeutronTable; format::Symbol=:ascii)
    format === :ascii || throw(ArgumentError("only :ascii (Type 1) supported"))
    nxs, jxs, xss, is_int = build_xss(table)
    h = table.header

    # Line 1: hz(a10) aw0(f12.6) ' ' tz(1pe11.4) ' ' hd(a10)
    @printf(io, "%s%12.6f %11.4E %s\n",
            rpad(h.hz, 10)[1:10], h.aw0, h.tz, rpad(h.hd, 10)[1:10])

    # Line 2: hk(a70) hm(a10)
    print(io, rpad(h.hk, 70)[1:70], rpad(h.hm, 10)[1:10], "\n")

    # Lines 3-6: IZ/AW pairs, 4(i7,f11.0) per line
    for row in 1:4
        for col in 1:4
            idx = (row - 1) * 4 + col
            @printf(io, "%7d%11.0f", table.iz[idx], table.aw[idx])
        end
        print(io, "\n")
    end

    # Lines 7-12: NXS(16)+JXS(32) = 48 integers, 8i9 per line
    all_ij = vcat(Int.(nxs), Int.(jxs))
    for row in 1:6
        for col in 1:8
            @printf(io, "%9d", all_ij[(row-1)*8 + col])
        end
        print(io, "\n")
    end

    # Data: XSS values, 4 per line, 20 chars each
    for i in 1:length(xss)
        if is_int[i]
            print(io, _ace_fmt_int(round(Int, xss[i])))
        else
            print(io, _ace_fmt_real(xss[i]))
        end
        if i % 4 == 0 || i == length(xss)
            print(io, "\n")
        end
    end
end

"""
    write_ace_directory(io::IO, table::ACENeutronTable, nxs;
                        itype=1, filepath="filename", route="route")

Write the MCNP xsdir directory line for an ACE table.
"""
function write_ace_directory(io::IO, table::ACENeutronTable,
                              nxs::AbstractVector{<:Integer};
                              itype::Int=1, filepath::String="filename",
                              route::String="route")
    h = table.header
    @printf(io, "%s%12.6f %s %s%2d%4d %8d%6d%6d%10.3E\n",
            rpad(strip(h.hz), 10)[1:10], h.aw0,
            filepath, route, itype, 1, nxs[NXS_LEN2], 0, 0, h.tz)
end

# write_ace dispatches to write_ace_table for ACENeutronTable
write_ace(io::IO, table::ACENeutronTable) = write_ace_table(io, table)

# =========================================================================
# Write flat ACETable (NXS/JXS/XSS representation)
# =========================================================================
"""
    write_ace(io::IO, table::ACETable)

Write a flat ACETable (NXS/JXS/XSS) as a Type 1 (ASCII) ACE file.
Output format matches NJOY2016 aceout() for itype=1.
"""
function write_ace(io::IO, table::ACETable)
    # Line 1: ZAID(a10), AWR(f12.6), ' ', TZ(e11.4), ' ', date(a10)
    @printf(io, "%s%12.6f %11.4E %s\n",
            rpad(table.zaid, 10)[1:10], table.awr, table.temp,
            rpad(table.date, 10)[1:10])

    # Line 2: comment(a70), mat_id(a10)
    print(io, rpad(table.comment, 70)[1:70], rpad(table.mat_id, 10)[1:10], "\n")

    # Lines 3-6: IZ/AW pairs, 4(i7,f11.0) per line
    for row in 1:4
        for col in 1:4
            idx = (row - 1) * 4 + col
            @printf(io, "%7d%11.0f", table.pairs_iz[idx], table.pairs_aw[idx])
        end
        print(io, "\n")
    end

    # Lines 7-12: NXS(16)+JXS(32) = 48 integers, 8i9 per line
    all_ints = Int32[table.nxs..., table.jxs...]
    for row in 1:6
        for col in 1:8
            @printf(io, "%9d", all_ints[(row - 1) * 8 + col])
        end
        print(io, "\n")
    end

    # Data: XSS values, 4 per line, 20 chars each (1pe20.11)
    n = length(table.xss)
    for i in 1:n
        print(io, _ace_fmt_real(table.xss[i]))
        if i % 4 == 0 || i == n
            print(io, "\n")
        end
    end
    return nothing
end
