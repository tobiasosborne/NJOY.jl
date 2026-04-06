#=
    Curve, grid, and hatching routines for the NJOY graph engine.
    Faithful port of graph.f90: grid2, grid3, curv2, curv3, hatch.

    Scale functions (xscale, yscale, zscale, xinvrs, yinvrs) are
    defined in axes.jl and used here without redefinition.
=#

# ---- Grid routines ----

"""
    grid2!(gs, io, nx, ny)

Draw a 2D grid. Wrapper around grid3! with z=0.
Matches graph.f90 `grid2` (lines 1218-1229).
"""
function grid2!(gs::GraphState, io::IO, nx::Int, ny::Int)
    gs.zmin = 0.0
    gs.zmax = 0.0
    grid3!(gs, io, nx, ny, 0)
end

"""
    grid3!(gs, io, nx, ny, nz)

Draw a 3D grid. One of nx/ny/nz must be 0 or -1 (the fixed plane).
Matches graph.f90 `grid3` (lines 1231-1428).
"""
function grid3!(gs::GraphState, io::IO, nx::Int, ny::Int, nz::Int)
    subs = (0.301, 0.477, 0.602, 0.699, 0.778, 0.845, 0.903, 0.954, 1.000)

    # X grid lines
    if nx > 0
        if nz == 0
            y1, y2 = gs.ymin, gs.ymax
            z1, z2 = gs.zmin, gs.zmin
        elseif nz == -1
            y1, y2 = gs.ymin, gs.ymax
            z1, z2 = gs.zmax, gs.zmax
        else
            if ny == 0
                y1, y2 = gs.ymin, gs.ymin
                z1, z2 = gs.zmin, gs.zmax
            else
                y1, y2 = gs.ymax, gs.ymax
                z1, z2 = gs.zmin, gs.zmax
            end
        end
        x1 = gs.xmin
        if gs.logx == 1
            v = gs.fvx
            i = round(Int, v)
            k = 0
            if abs(v - i - subs[1]) < 0.1
                k = 1
            end
            if abs(v - i - subs[4]) < 0.1
                k = 4
            end
            xt = x1
            if k > 0
                xt = x1 - subs[k] / gs.dvx
            end
            k += 1
        end
        idone = false
        while !idone
            if gs.logx == 0
                x1 += gs.xstp / nx
            else
                x1 = xt + subs[k] / gs.dvx
                k += 1
                if k > 9
                    xt = x1
                    k = 1
                end
            end
            if abs(gs.xmax - x1) < 0.01
                idone = true
            else
                vect3!(gs, io, x1, y1, z1, x1, y2, z2, gs.wg, 0)
            end
        end
    end

    # Y grid lines
    if ny > 0
        if nz == 0
            x1, x2 = gs.xmin, gs.xmax
            z1, z2 = gs.zmin, gs.zmin
        elseif nz == -1
            x1, x2 = gs.xmin, gs.xmax
            z1, z2 = gs.zmax, gs.zmax
        else
            if nx == 0
                x1, x2 = gs.xmin, gs.xmin
                z1, z2 = gs.zmin, gs.zmax
            else
                x1, x2 = gs.xmax, gs.xmax
                z1, z2 = gs.zmin, gs.zmax
            end
        end
        y1 = gs.ymin
        if gs.logy == 1
            v = gs.fvy
            i = round(Int, v)
            k = 0
            if abs(v - i - subs[1]) < 0.1
                k = 1
            end
            if abs(v - i - subs[4]) < 0.1
                k = 4
            end
            yt = y1
            if k > 0
                yt = y1 - subs[k] / gs.dvy
            end
            k += 1
        end
        idone = false
        while !idone
            if gs.logy == 0
                y1 += gs.ystp / ny
            else
                y1 = yt + subs[k] / gs.dvy
                k += 1
                if k > 9
                    yt = y1
                    k = 1
                end
            end
            if abs(gs.ymax - y1) < 0.01
                idone = true
            else
                vect3!(gs, io, x1, y1, z1, x2, y1, z2, gs.wg, 0)
            end
        end
    end

    # Z grid lines
    if nz > 0
        if nx == 0
            y1, y2 = gs.ymin, gs.ymax
            x1, x2 = gs.xmin, gs.xmin
        elseif nx == -1
            y1, y2 = gs.ymin, gs.ymax
            x1, x2 = gs.xmax, gs.xmax
        else
            if ny == 0
                y1, y2 = gs.ymin, gs.ymin
                x1, x2 = gs.xmin, gs.xmax
            else
                y1, y2 = gs.ymax, gs.ymax
                x1, x2 = gs.xmin, gs.xmax
            end
        end
        z1 = gs.zmin
        if gs.logz == 1
            v = gs.fvz
            i = round(Int, v)
            k = 0
            if abs(v - i - subs[1]) < 0.1
                k = 1
            end
            if abs(v - i - subs[4]) < 0.1
                k = 4
            end
            zt = z1
            if k > 0
                zt = z1 - subs[k] / gs.dvz
            end
            k += 1
        end
        idone = false
        while !idone
            if gs.logz == 0
                z1 += gs.zstp / nz
            else
                z1 = zt + subs[k] / gs.dvz
                k += 1
                if k > 9
                    zt = z1
                    k = 1
                end
            end
            if abs(gs.zmax - z1) < 0.01
                idone = true
            else
                vect3!(gs, io, x1, y1, z1, x2, y2, z1, gs.wg, 0)
            end
        end
    end
end

# ---- Curve routines ----

"""
    curv2!(gs, io, x, y, n, idash, width, icon, isym, ssym, ishade)

Draw a 2D curve with optional symbols, shading, and hatching.
Matches graph.f90 `curv2` (lines 1430-1505).
"""
function curv2!(gs::GraphState, io::IO,
                x::Vector{Float64}, y::Vector{Float64}, n::Int,
                idash::Int, width::Float64, icon::Int, isym::Int,
                ssym::Float64, ishade::Int)
    w = width

    # Draw curve if requested
    if icon >= 0
        wu1, wu2, wv1, wv2 = 100.0, -100.0, 100.0, -100.0

        xn = xscale(gs, x[1])
        yn = yscale(gs, y[1])
        wu, wv = trans3(gs, xn, yn, 0.0)
        wu1 = min(wu1, wu); wu2 = max(wu2, wu)
        wv1 = min(wv1, wv); wv2 = max(wv2, wv)
        u, v = transw(gs, wu, wv)
        moveh!(gs, io, u, v)

        for i in 2:n
            xn = xscale(gs, x[i])
            yn = yscale(gs, y[i])
            wu, wv = trans3(gs, xn, yn, 0.0)
            wu1 = min(wu1, wu); wu2 = max(wu2, wu)
            wv1 = min(wv1, wv); wv2 = max(wv2, wv)
            u, v = transw(gs, wu, wv)
            drawh!(gs, io, u, v, w, idash)
        end

        # Fill with shading pattern
        if ishade > 0 && ishade <= 10
            gs.ifg = ishade + 10
            fillh!(gs, io, 0.0)
            gs.ifg = 0
        elseif ishade > 40
            gs.ifg = ishade
            fillh!(gs, io, 0.0)
            gs.ifg = 0
        end

        # Fill with hatching
        if ishade > 10
            hatch!(gs, io, ishade, wu1, wv1, wu2, wv2)
        end

        # Make curve invisible (used with shaded areas)
        if width == 0.0
            ncurve!(gs, io)
        end
    end

    # Draw symbols if requested
    if icon != 0
        for i in 1:n
            if mod(i, abs(icon)) == 0
                xn = xscale(gs, x[i])
                yn = yscale(gs, y[i])
                wu, wv = trans3(gs, xn, yn, 0.0)
                u, v = transw(gs, wu, wv)
                dsym!(gs, io, u, v, isym, ssym)
            end
        end
    end
end

"""
    curv3!(gs, io, x, y, z, n, w)

Draw a 3D curve as a filled polygon (painter's algorithm hidden surface).
Matches graph.f90 `curv3` (lines 1507-1543).
"""
function curv3!(gs::GraphState, io::IO,
                x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64},
                n::Int, w::Float64)
    # Close the polygon
    push!(x, x[1])
    push!(y, y[1])
    push!(z, z[1])
    nmax = n + 1

    xn = xscale(gs, x[1])
    yn = yscale(gs, y[1])
    zn = zscale(gs, z[1])
    wu, wv = trans3(gs, xn, yn, zn)
    u, v = transw(gs, wu, wv)
    moveh!(gs, io, u, v)

    for i in 2:nmax
        xn = xscale(gs, x[i])
        yn = yscale(gs, y[i])
        zn = zscale(gs, z[i])
        wu, wv = trans3(gs, xn, yn, zn)
        u, v = transw(gs, wu, wv)
        drawh!(gs, io, u, v, w, 0)
    end
    fillh!(gs, io, 1.0)  # white fill for hiding
end

"""
    hatch!(gs, io, ihatch, wu1, wv1, wu2, wv2)

Fill current clipped path with diagonal hatching or cross-hatching.
Matches graph.f90 `hatch` (lines 1586-1632).
"""
function hatch!(gs::GraphState, io::IO, ihatch::Int,
                wu1::Float64, wv1::Float64, wu2::Float64, wv2::Float64)
    oh2 = 0.02
    w = 0.003

    iop = div(ihatch - 1, 10)
    dx = oh2 * (10 - ihatch + 10 * iop + 1)

    cclip!(gs, io)

    if iop == 1 || iop == 3
        # Forward diagonal (/)
        u1, v1 = transw(gs, wu1, wv1)
        u2, v2 = transw(gs, wu2, wv2)
        d = max(abs(u2 - u1), abs(v2 - v1))
        nd = round(Int, d / dx)
        for i in 1:nd
            xl = u1 + i * dx
            moveh!(gs, io, xl, v1)
            drawh!(gs, io, xl - d, v1 + d, w, 0)
        end
    end

    if iop == 2 || iop == 3
        # Backward diagonal (\)
        u1, v1 = transw(gs, wu1, wv1)
        u2, v2 = transw(gs, wu2, wv2)
        d = max(abs(u2 - u1), abs(v2 - v1))
        nd = round(Int, d / dx)
        for i in 1:nd
            xl = u1 + i * dx
            moveh!(gs, io, xl, v2)
            drawh!(gs, io, xl - d, v2 - d, w, 0)
        end
    end

    nclip!(gs, io)
end
