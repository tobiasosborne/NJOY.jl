#=
    Low-level PostScript output primitives for the NJOY graph engine.
    Faithful port of graph.f90 routines: gplot, gdone, newp, endp,
    moveh, drawh, fillh, gset, gend, cclip, nclip, ncurve.

    All format strings match Fortran exactly for bit-identical PS output.
=#

using Printf

"""
    gplot!(gs, io, lori, xpage, ypage)

Initialize PostScript output. Writes the DSC header, BoundingBox,
orientation, and `2 setlinecap`.
Matches graph.f90 `gplot` (lines 2103-2147).
"""
function gplot!(gs::GraphState, io::IO, lori::Int, xpage::Float64, ypage::Float64)
    gs.ipage = 0
    gs.land = lori

    if lori == 0  # portrait
        i1 = trunc(Int, (XPAPER - xpage) * 72 / 2)
        i2 = trunc(Int, (YPAPER - ypage) * 72 / 2)
        i3 = trunc(Int, i1 + xpage * 72)
        i4 = trunc(Int, i2 + ypage * 72)
        gs.ushift = Float64(i1)
        gs.vshift = Float64(i2)
        gs.uwidth = xpage * 72
    else  # landscape
        i1 = trunc(Int, (YPAPER - xpage) * 72 / 2)
        i2 = trunc(Int, (XPAPER - ypage) * 72 / 2)
        i3 = trunc(Int, i1 + ypage * 72)
        i4 = trunc(Int, i2 + xpage * 72)
        gs.ushift = Float64(i2)
        gs.vshift = Float64(i1)
        gs.uwidth = ypage * 72
    end

    println(io, "%!PS-Adobe-")
    @printf(io, "%%%%BoundingBox: %5d%5d%5d%5d\n", i1, i2, i3, i4)
    println(io, "%%Pages: (atend)")
    if lori == 1
        println(io, "%%Orientation: Landscape")
    else
        println(io, "%%Orientation: Portrait")
    end
    println(io, "2 setlinecap")
end

"""
    gdone!(gs, io)

Finalize PostScript output. Writes trailer with page count and %%EOF.
Matches graph.f90 `gdone` (lines 2149-2166).
"""
function gdone!(gs::GraphState, io::IO)
    println(io, "%gdone")
    println(io, "%%Trailer")
    if gs.ipage < 10
        @printf(io, "%%%%Pages: %1d\n", gs.ipage)
    elseif gs.ipage < 100
        @printf(io, "%%%%Pages: %2d\n", gs.ipage)
    else
        @printf(io, "%%%%Pages: %3d\n", gs.ipage)
    end
    println(io, "%%EOF")
end

"""
    newp!(gs, io)

Start a new PostScript page. Increments page counter.
Matches graph.f90 `newp` (lines 2197-2204).
"""
function newp!(gs::GraphState, io::IO)
    gs.ipage += 1
    @printf(io, "%%%%Page:%4d%4d\n", gs.ipage, gs.ipage)
end

"""
    endp!(gs, io)

End the current PostScript page: stroke, showpage, %endp.
Matches graph.f90 `endp` (lines 2206-2212).
"""
function endp!(gs::GraphState, io::IO)
    println(io, "stroke")
    println(io, "showpage")
    println(io, "%endp")
end

# ---- Coordinate transform helpers (portrait/landscape + clamping) ----

@inline function _to_ps_coords(gs::GraphState, x::Float64, y::Float64)
    if gs.land == 1
        u = gs.uwidth - 72 * y + gs.ushift
        v = 72 * x + gs.vshift
    else
        u = 72 * x + gs.ushift
        v = 72 * y + gs.vshift
    end
    u = clamp(u, -1000.0, 2000.0)
    v = clamp(v, -1000.0, 2000.0)
    return u, v
end

"""
    moveh!(gs, io, x, y)

Low-level PostScript move-to with landscape rotation and color setting.
Strokes previous path, starts new path, sets foreground color, emits moveto.
Matches graph.f90 `moveh` (lines 2214-2247).
"""
function moveh!(gs::GraphState, io::IO, x::Float64, y::Float64)
    u, v = _to_ps_coords(gs, x, y)

    println(io, "stroke")
    println(io, "newpath")

    rgb = 255.0
    r = IFRGB[1, gs.ifg] / rgb
    g = IFRGB[2, gs.ifg] / rgb
    b = IFRGB[3, gs.ifg] / rgb
    @printf(io, "%6.3f%6.3f%6.3f setrgbcolor\n", r, g, b)
    @printf(io, "%9.2f%9.2f moveto\n", u, v)
end

"""
    drawh!(gs, io, x, y, w, idash)

Low-level PostScript line-to with linewidth and dash pattern caching.
If `w < 0`, records state without drawing.
Matches graph.f90 `drawh` (lines 2249-2309).
"""
function drawh!(gs::GraphState, io::IO, x::Float64, y::Float64, w::Float64, idash::Int)
    if w < 0.0
        gs.wlast = w
        gs.ldash = idash
        return
    end

    # Emit setlinewidth only when width changes
    if w != gs.wlast
        @printf(io, "%8.3f setlinewidth\n", 72 * w)
        gs.wlast = w
    end

    # Emit setdash only when pattern changes
    if idash != gs.ldash
        wp = 72 * w  # width in points
        if idash == 0
            println(io, "[] 0 setdash")
        elseif idash == 1
            @printf(io, "[%9.2f%9.2f] 0 setdash\n", wp * 5, wp * 3)
        elseif idash == 2
            @printf(io, "[%9.2f%9.2f%9.2f%9.2f] 0 setdash\n",
                    wp * 6, wp * 3, wp * 3, wp * 3)
        elseif idash == 3
            @printf(io, "[%9.2f%9.2f%9.2f%9.2f] 0 setdash\n",
                    wp * 6, wp * 3, wp * 1, wp * 3)
        elseif idash == 4
            @printf(io, "[%9.2f%9.2f] 0 setdash\n", wp * 1, wp * 3)
        end
        gs.ldash = idash
    end

    u, v = _to_ps_coords(gs, x, y)
    @printf(io, "%9.2f%9.2f lineto\n", u, v)
end

"""
    fillh!(gs, io, color)

Fill current path. `color > 0.99` = background, `color < 0.01` = foreground.
Matches graph.f90 `fillh` (lines 2311-2350).
"""
function fillh!(gs::GraphState, io::IO, color::Float64)
    rgb = 255.0
    r, g, b = 0.0, 0.0, 0.0

    if color > 0.99
        # Background color
        r = IBRGB[1, gs.ibg] / rgb
        g = IBRGB[2, gs.ibg] / rgb
        b = IBRGB[3, gs.ibg] / rgb
    elseif color < 0.01 && gs.ifg <= 10
        # Curve colors
        r = IFRGB[1, gs.ifg] / rgb
        g = IFRGB[2, gs.ifg] / rgb
        b = IFRGB[3, gs.ifg] / rgb
    elseif color < 0.01 && gs.ifg <= 20
        # Shades of gray
        r = (20 - gs.ifg) / 10.0
        g = r
        b = r
    elseif color < 0.01 && gs.ifg > 20
        # Shade colors
        r = ISRGB[1, gs.ifg - 40] / rgb
        g = ISRGB[2, gs.ifg - 40] / rgb
        b = ISRGB[3, gs.ifg - 40] / rgb
    end

    @printf(io, "gsave%6.3f%6.3f%6.3f setrgbcolor fill grestore\n", r, g, b)
end

"""
    gset!(gs, io, ull, vll, uur, vur)

Set up a rectangular PostScript clipping path.
Matches graph.f90 `gset` (lines 2168-2187).
"""
function gset!(gs::GraphState, io::IO, ull::Float64, vll::Float64,
               uur::Float64, vur::Float64)
    println(io, "%gset")
    w = 0.005
    moveh!(gs, io, ull, vll)
    drawh!(gs, io, uur, vll, w, 0)
    drawh!(gs, io, uur, vur, w, 0)
    drawh!(gs, io, ull, vur, w, 0)
    drawh!(gs, io, ull, vll, w, 0)
    println(io, " gsave clip newpath")
end

"""
    gend!(gs, io)

End a clipping path. Matches graph.f90 `gend` (lines 2189-2195).
"""
function gend!(gs::GraphState, io::IO)
    println(io, "stroke grestore newpath")
    println(io, "%gend")
end

"""
    cclip!(gs, io)

Push a clipping context on the current path.
Matches graph.f90 `cclip` (lines 2352-2358).
"""
function cclip!(gs::GraphState, io::IO)
    println(io, "%cclip")
    println(io, "gsave clip newpath")
end

"""
    nclip!(gs, io)

Pop the clipping context. Matches graph.f90 `nclip` (lines 2360-2366).
"""
function nclip!(gs::GraphState, io::IO)
    println(io, "stroke grestore")
    println(io, "%nclip")
end

"""
    ncurve!(gs, io)

Cancel the current path (for shaded fills without borders).
Matches graph.f90 `ncurve` (lines 2368-2374).
"""
function ncurve!(gs::GraphState, io::IO)
    println(io, "newpath")
end
