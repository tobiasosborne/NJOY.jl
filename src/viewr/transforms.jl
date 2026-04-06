#=
    Coordinate transforms and page geometry for the NJOY graph engine.
    Faithful port of graph.f90 subroutines: transw, trans3, initp,
    window, endw, init2, frame2, endfr, vect2, vect3, poly2, poly3.

    All format strings and numerical constants match Fortran exactly
    for bit-identical PostScript output.
=#

using Printf

# ---- Coordinate transforms ----

"""
    transw(gs, x, y) -> (u, v)

Window coordinate transform with optional rotation.
Maps window-local `(x, y)` to page coordinates `(u, v)`.
If no window is active (`www == 0`), returns identity.
Matches graph.f90 `transw` (lines 363-384).
"""
function transw(gs::GraphState, x::Float64, y::Float64)
    if gs.www == 0.0
        return (x, y)
    else
        ct = cos(2π * gs.wwr / 360)
        st = sin(2π * gs.wwr / 360)
        u = gs.xwll + x * ct - y * st
        v = gs.ywll + x * st + y * ct
        return (u, v)
    end
end

"""
    trans3(gs, x, y, z) -> (u, v)

3D perspective projection.
Maps 3D coordinates `(x, y, z)` to 2D window coordinates `(u, v)`
using the current viewing parameters (ro, rs, ct, st, cp, sp, du, dv).
Matches graph.f90 `trans3` (lines 571-583).
"""
function trans3(gs::GraphState, x::Float64, y::Float64, z::Float64)
    xo = gs.ct * (gs.cp * x + gs.sp * y) + gs.st * z
    s = gs.rs / (gs.ro - xo)
    u = s * (-gs.sp * x + gs.cp * y) + gs.du
    v = s * (-gs.st * (gs.cp * x + gs.sp * y) + gs.ct * z) + gs.dv
    return (u, v)
end

# ---- Page initialization ----

"""
    initp!(gs, io, iori, xpage, ypage, istyle, htt, wtt, ibord, ipcol)

Initialize the page. Sets up page size, orientation, font style,
font height, normal line weight, and page background color.
Draws page border if requested.
Matches graph.f90 `initp` (lines 212-283).
"""
function initp!(gs::GraphState, io::IO, iori::Int, xpage::Float64,
                ypage::Float64, istyle::Int, htt::Float64, wtt::Float64,
                ibord::Int, ipcol::Int)
    # width of border line
    wbord = 0.005

    gs.land = iori
    gs.ifont = istyle
    gs.width = wtt
    gs.ht = htt
    gs.ibg = 1 + ipcol
    gs.ifg = 1
    gs.ro = 1.0
    gs.rs = 1.0
    gs.ct = 0.0
    gs.st = 1.0
    gs.cp = 0.0
    gs.sp = -1.0
    gs.du = 0.0
    gs.dv = 0.0
    gs.xwll = 0.0
    gs.ywll = 0.0
    gs.www = 0.0
    gs.wwh = 0.0
    gs.wwr = 0.0

    newp!(gs, io)

    n = 5
    x = Vector{Float64}(undef, 5)
    y = Vector{Float64}(undef, 5)
    x[1] = 0.0
    y[1] = 0.0
    x[2] = xpage
    y[2] = 0.0
    x[3] = xpage
    y[3] = ypage
    x[4] = 0.0
    y[4] = ypage
    x[5] = x[1]
    y[5] = y[1]

    w = wbord
    if gs.ibg == 1
        w = -1.0
    end
    backgr = 1.0

    poly2!(gs, io, x, y, n, w, backgr)
end

# ---- Window setup ----

"""
    window!(gs, io, xll, yll, ww, wh, wr, t1, n1, t2, n2, ibord)

Set up a window on the page and write titles (if any) in the upper
left corner. Draws a window border if requested. Sets up a clipping
path at the window border.
Matches graph.f90 `window` (lines 286-341).
"""
function window!(gs::GraphState, io::IO, xll::Float64, yll::Float64,
                 ww::Float64, wh::Float64, wr::Float64,
                 t1::String, n1::Int, t2::String, n2::Int, ibord::Int)
    # width of window border
    wbord = 0.005

    gs.ifg = 1

    # load window transformation
    gs.xwll = xll
    gs.ywll = yll
    gs.www = ww
    gs.wwh = wh
    gs.wwr = wr

    # draw window border
    if ibord != 0
        w = wbord
        vect2!(gs, io, 0.0, 0.0, 0.0, gs.wwh, w, 0)
        vect2!(gs, io, 0.0, gs.wwh, gs.www, gs.wwh, w, 0)
        vect2!(gs, io, gs.www, gs.wwh, gs.www, 0.0, w, 0)
        vect2!(gs, io, gs.www, 0.0, 0.0, 0.0, w, 0)
    end

    # set clipping to window border
    ull, vll = transw(gs, 0.0, 0.0)
    uur, vur = transw(gs, gs.www, gs.wwh)
    gset!(gs, io, ull, vll, uur, vur)

    # draw titles, if any
    gs.tspace = 0.0
    if n1 != 0
        x = 3 * gs.ht
        y = gs.wwh - 3 * gs.ht / 2
        gs.tspace += 3 * gs.ht / 2
        text2!(gs, io, t1, n1, gs.ifont, gs.ht, x, y, 1.0, 0.0, 0.0, 1.0)
        if n2 != 0
            y = y - (gs.ht + gs.ht / 5)
            gs.tspace += 3 * gs.ht / 2
            text2!(gs, io, t2, n2, gs.ifont, gs.ht, x, y, 1.0, 0.0, 0.0, 1.0)
        end
    end
end

"""
    endw!(gs, io)

End a window by resetting the clipping path and restoring the
default 3D transformation to identity.
Matches graph.f90 `endw` (lines 343-361).
"""
function endw!(gs::GraphState, io::IO)
    gend!(gs, io)
    gs.ro = 1.0
    gs.rs = 1.0
    gs.ct = 0.0
    gs.st = 1.0
    gs.cp = 0.0
    gs.sp = -1.0
    gs.du = 0.0
    gs.dv = 0.0
    gs.www = 0.0
    gs.wwh = 0.0
    gs.wwr = 0.0
end

# ---- 2D frame setup ----

"""
    init2!(gs, io, iright, iwcol) -> (uo, vo, xg, yg)

Initialize a 2D plot by choosing an origin and axis lengths to
nicely position the graph in the window. Returns the plot origin
`(uo, vo)` and axis dimensions `(xg, yg)` in window coordinates.
`iright > 0` reserves space for a right-hand axis.
`iwcol` sets the background color inside the graph frame.
Matches graph.f90 `init2` (lines 386-458).
"""
function init2!(gs::GraphState, io::IO, iright::Int, iwcol::Int)
    # width of frame border
    wframe = 0.005

    # adjust axis lengths to fit graph inside window
    uo = 4 * gs.ht
    vo = 4 * gs.ht
    xg = gs.www - uo - 2 * gs.ht
    if iright > 0
        xg = xg - 2 * gs.ht
    end
    yg = gs.wwh - gs.tspace - vo - gs.ht
    if gs.tspace == 0.0
        yg = yg - gs.ht
    end

    # set up axis parameters
    gs.wt = gs.width
    gs.wg = gs.wt / 2
    gs.wa = 2 * gs.wt
    gs.tic = gs.ht / 2
    gs.gap = 4 * gs.ht / 10
    gs.lfont = gs.ifont
    gs.hf = 3.0 / 4.0
    gs.hl = gs.ht
    gs.hn = 8 * gs.hl / 10

    # set up coordinate transformation (no perspective)
    gs.ro = 1.0
    gs.rs = 1.0
    gs.ct = 0.0
    gs.st = 1.0
    gs.cp = 0.0
    gs.sp = -1.0
    gs.du = uo
    gs.dv = vo

    # color in the background inside the graph frame
    gs.ibg = 1 + iwcol
    gs.ifg = 1
    ull, vll = transw(gs, uo, vo)
    uur, vur = transw(gs, uo + xg, vo + yg)
    n = 5
    x = Vector{Float64}(undef, 5)
    y = Vector{Float64}(undef, 5)
    x[1] = ull
    y[1] = vll
    x[2] = uur
    y[2] = vll
    x[3] = uur
    y[3] = vur
    x[4] = ull
    y[4] = vur
    x[5] = x[1]
    y[5] = y[1]
    w = wframe
    backgr = 1.0
    poly2!(gs, io, x, y, n, w, backgr)

    return (uo, vo, xg, yg)
end

"""
    frame2!(gs, io, xg, yg, grace, iccol)

Set up a clipping path just outside the frame boundary and set the
foreground color for this curve.
Matches graph.f90 `frame2` (lines 460-478).
"""
function frame2!(gs::GraphState, io::IO, xg::Float64, yg::Float64,
                 grace::Float64, iccol::Int)
    gs.ifg = 1 + iccol
    xll, yll = trans3(gs, -grace, -grace, 0.0)
    xur, yur = trans3(gs, xg + grace, yg + grace, 0.0)
    ull, vll = transw(gs, xll, yll)
    uur, vur = transw(gs, xur, yur)
    gset!(gs, io, ull, vll, uur, vur)
end

"""
    endfr!(gs, io)

End the 2D plotting frame by clearing the frame clipping path.
Matches graph.f90 `endfr` (lines 480-485).
"""
function endfr!(gs::GraphState, io::IO)
    gend!(gs, io)
end

# ---- Vector drawing ----

"""
    vect2!(gs, io, x1, y1, x2, y2, w, ivec)

Draw a 2D vector by delegating to `vect3!` with `z = 0`.
Matches graph.f90 `vect2` (lines 585-597).
"""
function vect2!(gs::GraphState, io::IO, x1::Float64, y1::Float64,
                x2::Float64, y2::Float64, w::Float64, ivec::Int)
    vect3!(gs, io, x1, y1, 0.0, x2, y2, 0.0, w, ivec)
end

"""
    vect3!(gs, io, x1, y1, z1, x2, y2, z2, w, ivec)

Draw a 3D vector with optional arrowhead (`ivec != 0`).
Arrow head size is 0.2 inches, drawn at the endpoint.
Matches graph.f90 `vect3` (lines 599-636).
"""
function vect3!(gs::GraphState, io::IO,
                x1::Float64, y1::Float64, z1::Float64,
                x2::Float64, y2::Float64, z2::Float64,
                w::Float64, ivec::Int)
    # size of vector arrowhead
    ahead = 0.2

    wu1, wv1 = trans3(gs, x1, y1, z1)
    wu2, wv2 = trans3(gs, x2, y2, z2)
    u1, v1 = transw(gs, wu1, wv1)
    u2, v2 = transw(gs, wu2, wv2)

    if ivec == 0
        moveh!(gs, io, u1, v1)
        drawh!(gs, io, u2, v2, w, 0)
    else
        head = ahead
        r = sqrt((u2 - u1)^2 + (v2 - v1)^2)
        wx = w * (u2 - u1) / r
        wy = w * (v2 - v1) / r
        dx = head * (u2 - u1) / r
        dy = head * (v2 - v1) / r
        hx = dx / 2
        hy = dy / 2
        moveh!(gs, io, u1, v1)
        drawh!(gs, io, u2 - wx, v2 - wy, w, 0)
        moveh!(gs, io, u2 - dx - hy, v2 - dy + hx)
        drawh!(gs, io, u2, v2, w, 0)
        drawh!(gs, io, u2 - dx + hy, v2 - dy - hx, w, 0)
    end
end

# ---- Polygon drawing ----

"""
    poly2!(gs, io, x, y, n, w, fill)

Draw a filled 2D polygon. If `w < 0`, suppress the outline
(drawh records state without drawing when w is negative).
Matches graph.f90 `poly2` (lines 1545-1562).
"""
function poly2!(gs::GraphState, io::IO,
                x::AbstractVector{Float64}, y::AbstractVector{Float64},
                n::Int, w::Float64, fill::Float64)
    moveh!(gs, io, x[1], y[1])
    for i in 2:n
        drawh!(gs, io, x[i], y[i], w, 0)
    end
    fillh!(gs, io, fill)
end

"""
    poly3!(gs, io, x, y, z, n, w)

Draw a 3D polygon with hiding. Each vertex is transformed through
`trans3` then `transw` before drawing.
Matches graph.f90 `poly3` (lines 1564-1584).
"""
function poly3!(gs::GraphState, io::IO,
                x::AbstractVector{Float64}, y::AbstractVector{Float64},
                z::AbstractVector{Float64}, n::Int, w::Float64)
    wu, wv = trans3(gs, x[1], y[1], z[1])
    u, v = transw(gs, wu, wv)
    moveh!(gs, io, u, v)
    for i in 2:n
        wu, wv = trans3(gs, x[i], y[i], z[i])
        u, v = transw(gs, wu, wv)
        drawh!(gs, io, u, v, w, 0)
    end
end
