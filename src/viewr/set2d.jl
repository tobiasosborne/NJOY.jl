#=
    2D/3D plot drivers, auto-scaling, error bars, legend, and tag routines.
    Faithful port of viewr.f90: ascalv, set2d, erbar, legndb, tagit, set3d.
    Also includes init3! from graph.f90 (lines 487-569).

    All numerical constants and control flow match Fortran exactly.
=#

using Printf

# ---- Auto-scaling (Los Alamos SC4020 library) ----

"""
    ascalv(m, z1, z2) -> (z1_new, z2_new, major, minor)

Automatic scaling routine for a linear axis.
Borrowed from Los Alamos SC4020 library (with modifications).

On input: `m` = minimum number of major divisions desired,
`z1`/`z2` = min/max values of data to be plotted.
Returns new axis limits and division counts.

Matches viewr.f90 `ascalv` (lines 1526-1667).
"""
function ascalv(m::Int, z1::Float64, z2::Float64)
    zmin = z1
    zmax = z2

    # check input -- exit if parameters are unreasonable
    if zmax <= zmin || m <= 0 || m > 20
        return (0.0, 2 * z2, 1, 0)
    end

    # find minimum span per interval
    fm = Float64(m)
    if zmax == 0.0 || zmin == 0.0
        zmax = zmax - abs(zmax) / 1000000
        zmin = zmin + abs(zmin) / 1000000
    else
        zbar = zmax / zmin
        if abs(zbar) / 1000 >= 1.0
            zmin = 0.0
        elseif abs(zbar) * 1000 <= 1.0
            zmax = 0.0
            zmax = zmax - abs(zmax) / 1000000
            zmin = zmin + abs(zmin) / 1000000
        else
            if (abs(zbar - 1) - 5 * fm / 100000) > 0.0
                zmax = zmax - abs(zmax) / 1000000
                zmin = zmin + abs(zmin) / 1000000
            else
                zbar = (zmax + zmin) / 2
                z = 26 * fm * abs(zbar) / 1000000
                zmax = zbar + z
                zmin = zbar - z
            end
        end
    end
    p = (zmax - zmin) / fm

    # determine k such that 10^k <= p < 10^(k+1)
    iflag = 0
    tenk = 1.0
    k = 0
    if p < 1.0
        iflag = 1
        p = 1 / p
    end
    while p >= 10000.0
        p = p / 10000
        tenk = tenk * 10000
        k = k + 4
    end
    while p >= 10.0
        p = p / 10
        tenk = tenk * 10
        k = k + 1
    end
    if iflag != 0
        p = 10 / p
        tenk = 1.0 / 10 / tenk
        k = k - 1
    end

    # determine dz = 1*10^k, 2*10^k, or 5*10^k
    nm = 0
    test = 2 + 1.0 / 100
    if p <= test
        test = 2 - 1.0 / 100
        if p <= test
            p = 1.0
            nm = 5
        else
            p = 2.0
            nm = 4
        end
    else
        test = 5 - 1.0 / 100
        if p < test
            p = 2.0
            nm = 4
        else
            p = 5.0
            nm = 5
        end
    end
    dz = p * tenk

    # find n1 and n2 such that n1*dz <= zmin < (n1+1)*dz
    # and (n2-1)*dz < zmax <= n2*dz
    n1 = trunc(Int, zmin / dz)
    fn = Float64(n1)
    z = fn * dz
    if z > zmin
        z = z - dz
        n1 = n1 - 1
    end
    zmin = z
    z1_new = zmin

    n2 = trunc(Int, zmax / dz)
    fn = Float64(n2)
    z = fn * dz
    if z < zmax
        n2 = n2 + 1
        z = z + dz
    end
    zmax = z
    z2_new = zmax

    major = n2 - n1
    minor = nm * major
    return (z1_new, z2_new, major, minor)
end

# ---- 3D perspective initialization ----

"""
    init3!(gs, io, bx, by, bz, vx, vy, vz, iwcol)

Initialize a 3D plot by setting up the perspective transformation
and auto-centering the projection in the window.
Matches graph.f90 `init3` (lines 487-569).
"""
function init3!(gs::GraphState, io::IO,
                bx::Float64, by::Float64, bz::Float64,
                vx::Float64, vy::Float64, vz::Float64, iwcol::Int)
    # set color of 3d slices to the window color
    gs.ibg = 1 + iwcol
    gs.ifg = 1

    # set up axis parameters
    gs.wt = gs.width
    gs.wg = gs.wt / 2
    gs.wa = gs.wt
    gs.tic = gs.ht / 2
    gs.gap = 4 * gs.ht / 10
    gs.lfont = gs.ifont
    gs.hf = 3.0 / 4.0
    gs.hl = gs.ht
    gs.hn = gs.ht

    # set up 3-d transformation
    gs.ro = sqrt(vx^2 + vy^2 + vz^2)
    gs.rp = sqrt(vx^2 + vy^2)
    gs.ct = gs.rp / gs.ro
    gs.st = vz / gs.ro
    gs.cp = vx / gs.rp
    gs.sp = vy / gs.rp
    gs.rs = gs.ro / sqrt(bx^2 + by^2 + bz^2)
    gs.du = 0.0
    gs.dv = 0.0

    # adjust view scale and position to center the projection in the window
    umin = 1000.0
    umax = -1000.0
    vmin = 1000.0
    vmax = -1000.0

    # corner 1
    yy = -5 * gs.ht
    if bx > 0.0
        u, v = trans3(gs, 0.0, yy, 0.0)
    else
        u, v = trans3(gs, bx, yy, 0.0)
    end
    if u < umin; umin = u; end
    if u > umax; umax = u; end
    if v < vmin; vmin = v; end
    if v > vmax; vmax = v; end

    # corner 2
    if bx > 0.0
        u, v = trans3(gs, 0.0, by, bz + gs.ht)
    else
        u, v = trans3(gs, bx, by, bz + gs.ht)
    end
    if u < umin; umin = u; end
    if u > umax; umax = u; end
    if v < vmin; vmin = v; end
    if v > vmax; vmax = v; end

    # corner 3
    xx = bx + 2 * gs.ht
    yy = -7 * gs.ht / 2
    if bx > 0.0
        u, v = trans3(gs, xx, yy, 0.0)
    else
        u, v = trans3(gs, 0.0, yy, 0.0)
    end
    if u < umin; umin = u; end
    if u > umax; umax = u; end
    if v < vmin; vmin = v; end
    if v > vmax; vmax = v; end

    # corner 4
    dx = 3 * gs.ht
    dy = 2 * gs.ht
    if bx > 0.0
        u, v = trans3(gs, bx + dx, by + dy, 0.0)
    else
        u, v = trans3(gs, dx, by + dy, 0.0)
    end
    if u < umin; umin = u; end
    if u > umax; umax = u; end
    if v < vmin; vmin = v; end
    if v > vmax; vmax = v; end

    su = gs.www / (umax - umin)
    sv = gs.wwh / (vmax - vmin)
    s = su
    if sv < su; s = sv; end
    gs.rs = s * gs.rs
    gs.du = (gs.www - s * (umax - umin)) / 2 - s * umin
    gs.dv = (gs.wwh - s * (vmax - vmin)) / 2 - s * vmin
end

# ---- 2D plot rendering ----

"""
    set2d!(gs, io, iplot, x, y, dxm, dxp, dym, dyp, n, ierrb, vs,
           xll, yll, ww, wh, wr)

Set up 2-dimensional multi-curve plots.
Window geometry (xll, yll, ww, wh, wr) comes from card 2.
Matches viewr.f90 `set2d` (lines 757-1040).
"""
function set2d!(gs::GraphState, io::IO, iplot::Int,
                x::Vector{Float64}, y::Vector{Float64},
                dxm::Vector{Float64}, dxp::Vector{Float64},
                dym::Vector{Float64}, dyp::Vector{Float64},
                n::Int, ierrb::Int, vs::ViewrState,
                xll::Float64, yll::Float64, ww::Float64,
                wh::Float64, wr::Float64)
    wt = vs.wline
    bard = 0.03
    small = 1.0e-15
    twenty = 20.0e6
    shade = 0.99
    ten = 10.0
    grace = 0.0

    # initial settings
    if iplot == 1
        if vs.csize > 0.0
            vs.hlab = vs.csize
        end
        if vs.csize < 0.0
            vs.hlab = -vs.csize * wh
        end
        vs.hleg = 2 * vs.hlab / 3
        initp!(gs, io, vs.lori, vs.xpage, vs.ypage,
               vs.istyle, vs.hlab, vs.wline, 0, vs.ipcol)
    end

    # set up subplot area for current axes
    if abs(iplot) == 1
        window!(gs, io, xll, yll, ww, wh, wr,
                vs.t1, vs.n1, vs.t2, vs.n2, 0)
        if vs.itype == 0
            return
        end
        _, _, vs.xg, vs.yg = init2!(gs, io, vs.jtype, vs.iwcol)
    end

    # set up regular axes or alternate y scale
    if iplot <= 1
        if iplot < -1
            vs.ymin = vs.zmin
            vs.ymax = vs.zmax
            vs.ystp = vs.zstp
        end

        # check for automatic scaling
        xstep = vs.xstp
        ystep = vs.ystp
        if (vs.xstp == 0.0 || vs.ystp == 0.0) && n > 0

            # search for min and max values
            if vs.xstp == 0.0
                vs.xmin = x[1]
                vs.xmax = x[1]
            end
            if vs.ystp == 0.0 || iplot >= -1
                if vs.ystp == 0.0
                    vs.ymax = y[1]
                    vs.ymin = y[1]
                end
                for i in 2:n
                    if vs.xstp == 0.0 && iplot >= -1
                        if x[i] < vs.xmin; vs.xmin = x[i]; end
                        if x[i] > vs.xmax; vs.xmax = x[i]; end
                    end
                    if vs.ystp == 0.0
                        if y[i] < vs.ymin && y[i] != small; vs.ymin = y[i]; end
                        if y[i] > vs.ymax; vs.ymax = y[i]; end
                    end
                end

                # choose limits and step sizes for linear axes
                if vs.itype <= 2 && iplot >= -1
                    if vs.xstp == 0.0
                        vs.xmin, vs.xmax, major, _ = ascalv(4, vs.xmin, vs.xmax)
                        xstep = (vs.xmax - vs.xmin) / major
                    end
                end
                if iplot < -1 || vs.itype == 1 || vs.itype == 3
                    if iplot >= -1 || vs.jtype == 1
                        if vs.ystp == 0.0
                            vs.ymin, vs.ymax, major, _ = ascalv(4, vs.ymin, vs.ymax)
                            ystep = (vs.ymax - vs.ymin) / major
                        end
                    end
                end
            end
        end

        # adjust limits for log scales
        # limit number of x and y decades
        if iplot >= -1 && vs.itype >= 3 && vs.xstp == 0.0
            if vs.xmax <= 0.0; vs.xmax = 1.0; end
            if vs.xmin <= 0.0; vs.xmin = vs.xmax / 1e8; end
            top = log10(vs.xmax)
            if vs.xmax != twenty
                if top >= 0.0
                    vs.xmax = ten^trunc(Int, top + shade)
                end
                if top < 0.0
                    vs.xmax = ten^trunc(Int, top)
                end
            end
            xlft = log10(vs.xmin)
            ixlft = trunc(Int, xlft)
            if xlft > 0.0
                vs.xmin = ten^ixlft
            end
            if xlft < 0.0
                vs.xmin = ten^(ixlft - 1)
            end
            if vs.xmax / vs.xmin > ten^13
                vs.xmin = ten^trunc(Int, top - 13)
            end
        end
        if iplot < -1 || vs.itype == 2 || vs.itype == 4
            if iplot >= -1 || vs.jtype == 2
                if vs.ystp == 0.0
                    if vs.ymax <= 0.0; vs.ymax = 1.0; end
                    if vs.ymin <= 0.0; vs.ymin = vs.ymax / 1e10; end
                    top = log10(vs.ymax)
                    if top >= 0.0
                        vs.ymax = ten^trunc(Int, top + shade)
                    end
                    if top < 0.0
                        vs.ymax = ten^trunc(Int, top)
                    end
                    bot = log10(vs.ymin)
                    ibot = trunc(Int, bot)
                    if bot > 0.0
                        vs.ymin = ten^ibot
                    end
                    if bot < 0.0
                        vs.ymin = ten^(ibot - 1)
                    end
                    if vs.ymax / vs.ymin > ten^8
                        vs.ymin = ten^trunc(Int, top - 8)
                    end
                end
            end
        end

        # set up the axes and grid
        if iplot >= -1
            vs.itic = 0
            if vs.igrid == 2; vs.itic = -1; end
            if vs.igrid == 3; vs.itic = 1; end
            xop = xstep
            yop = ystep
            if vs.itype == 3 || vs.itype == 4; xop = 0.0; end
            if vs.itype == 2 || vs.itype == 4; yop = 0.0; end

            # bottom x axis (labeled)
            axis2!(gs, io, vs.xmin, vs.xmax, xop, vs.xl, -vs.nx,
                   vs.itic, 0, vs.xg,
                   0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0)
            # top x axis (unlabeled)
            axis2!(gs, io, vs.xmin, vs.xmax, xop, " ", 0,
                   -vs.itic, 0, vs.xg,
                   0.0, vs.yg, 1.0, 0.0, 0.0, 1.0, 0)
            # left y axis (labeled)
            axis2!(gs, io, vs.ymin, vs.ymax, yop, vs.yl, vs.ny,
                   -vs.itic, 0, vs.yg,
                   0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 1)
            # right y axis (unlabeled)
            axis2!(gs, io, vs.ymin, vs.ymax, yop, " ", 0,
                   vs.itic, 0, vs.yg,
                   vs.xg, 0.0, 0.0, 1.0, -1.0, 0.0, 1)

            if vs.igrid == 1
                grid2!(gs, io, 1, 1)
            end
        else
            # alternate y axis
            rmin = vs.ymin
            rmax = vs.ymax
            rstep = ystep
            rop = rstep
            if vs.jtype == 2; rop = 0.0; end
            axis2!(gs, io, rmin, rmax, rop, vs.rl, -vs.nr,
                   -vs.itic, 0, vs.yg,
                   vs.xg, 0.0, 0.0, 1.0, -1.0, 0.0, 1)
            axis2!(gs, io, rmin, rmax, rop, " ", 0,
                   vs.itic, 0, vs.yg,
                   0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 1)
        end

        # locate legend block, if any
        if vs.ileg == 1
            if vs.xtag == 0.0 && vs.ytag == 0.0
                vs.xpt = xscale(gs, vs.xmin) + vs.hleg
                vs.ypt = yscale(gs, vs.ymax) - 3 * vs.hleg / 2
            else
                vs.xpt = xscale(gs, vs.xtag)
                vs.ypt = yscale(gs, vs.ytag)
            end
        end
    end

    # draw the curves on the new or existing frame
    # using a clipping path just outside of the frame
    frame2!(gs, io, vs.xg, vs.yg, grace, vs.iccol)
    if n > 0
        bar = bard
        if ierrb > 0
            erbar!(gs, io, x, y, n, dxm, dxp, dym, dyp, ierrb, bar, vs)
        end
        w = vs.ithick * wt
        ssym = 3 * vs.hleg / 10
        curv2!(gs, io, x, y, n, vs.idash, w, vs.icon, vs.isym, ssym, vs.ishade)
    end

    # add curve tags or legend lines
    if vs.ileg == 1
        legndb!(gs, io, iplot, vs)
    end
    if vs.ileg == 2
        tagit!(gs, io, x, y, n, vs.itype, vs)
    end

    # finished with this call to setup
    endfr!(gs, io)
end

# ---- Error bars ----

"""
    erbar!(gs, io, xx, yy, n, dxm, dxp, dym, dyp, ixy, w, vs)

Put error bars on previously drawn plot.
Matches viewr.f90 `erbar` (lines 1042-1103).
"""
function erbar!(gs::GraphState, io::IO,
                xx::Vector{Float64}, yy::Vector{Float64}, n::Int,
                dxm::Vector{Float64}, dxp::Vector{Float64},
                dym::Vector{Float64}, dyp::Vector{Float64},
                ixy::Int, w::Float64, vs::ViewrState)
    wt = vs.wline
    capd = 0.02

    # set up cap size
    wd = w
    if wd == 0.0; wd = capd; end

    # loop over data points
    for np in 1:n
        xn = xscale(gs, xx[np])
        yn = yscale(gs, yy[np])

        # draw x error bars, if any
        if ixy > 1
            if dxm[np] != 0.0 || dxp[np] != 0.0
                x1 = xscale(gs, xx[np] - dxm[np])
                x2 = xscale(gs, xx[np] + dxp[np])
                vect3!(gs, io, x1, yn, 0.0, x2, yn, 0.0, wt, 0)
                vect3!(gs, io, x1, yn + wd, 0.0, x1, yn - wd, 0.0, wt, 0)
                vect3!(gs, io, x2, yn + wd, 0.0, x2, yn - wd, 0.0, wt, 0)
            end
        end

        # draw y error bars, if any
        if ixy != 0 && ixy != 2
            if xx[np] >= vs.xmin && xx[np] <= vs.xmax
                if yy[np] >= vs.ymin && yy[np] <= vs.ymax
                    if dym[np] != 0.0 || dyp[np] != 0.0
                        y1 = yscale(gs, yy[np] - dym[np])
                        y2 = yscale(gs, yy[np] + dyp[np])
                        vect3!(gs, io, xn, y1, 0.0, xn, y2, 0.0, wt, 0)
                        vect3!(gs, io, xn - wd, y1, 0.0, xn + wd, y1, 0.0, wt, 0)
                        vect3!(gs, io, xn - wd, y2, 0.0, xn + wd, y2, 0.0, wt, 0)
                    end
                end
            end
        end
    end
end

# ---- Legend block ----

"""
    legndb!(gs, io, iplot, vs)

Write the legend, a line at a time.
Uses gs.leg_x, gs.leg_y as persistent state (Fortran SAVE x,y).
Matches viewr.f90 `legndb` (lines 1105-1177).
"""
function legndb!(gs::GraphState, io::IO, iplot::Int, vs::ViewrState)
    wt = vs.wline
    wd = 0.005
    backgr = 1.0

    if abs(iplot) == 1
        gs.leg_x = vs.xpt
        gs.leg_y = vs.ypt
        xta, yta = trans3(gs, gs.leg_x, gs.leg_y, 0.0)
        wleg = txtlen(gs, vs.aleg, vs.nleg, vs.istyle, vs.hleg)
        xb = Vector{Float64}(undef, 5)
        yb = Vector{Float64}(undef, 5)
        xb[1] = xta - 2 * vs.hleg / 10
        yb[1] = yta - 3 * vs.hleg / 10
        xb[2] = xta + wleg + 7.0 / 10 + 4 * vs.hleg / 10
        yb[2] = yta - 3 * vs.hleg / 10
        xb[3] = xta + wleg + 7.0 / 10 + 4 * vs.hleg / 10
        yb[3] = yta
        xb[4] = xta - 2 * vs.hleg / 10
        yb[4] = yta
        xb[5] = xta - 2 * vs.hleg / 10
        yb[5] = yta - 3 * vs.hleg / 10
        w = wd
        poly2!(gs, io, xb, yb, 5, -w, backgr)
    end

    xc = Vector{Float64}(undef, 2)
    yc = Vector{Float64}(undef, 2)
    xc[1] = xinvrs(gs, gs.leg_x)
    yc[1] = yinvrs(gs, gs.leg_y - 7 * vs.hleg / 10)
    xc[2] = xinvrs(gs, gs.leg_x + 5.0 / 10)
    yc[2] = yinvrs(gs, gs.leg_y - 7 * vs.hleg / 10)

    xta, yta = trans3(gs, gs.leg_x, gs.leg_y, 0.0)
    wleg = txtlen(gs, vs.aleg, vs.nleg, vs.istyle, vs.hleg)
    xb = Vector{Float64}(undef, 5)
    yb = Vector{Float64}(undef, 5)
    xb[1] = xta - 2 * vs.hleg / 10
    yb[1] = yta - 14 * vs.hleg / 10
    xb[2] = xta + wleg + 7.0 / 10 + 4 * vs.hleg / 10
    yb[2] = yta - 14 * vs.hleg / 10
    xb[3] = xta + wleg + 7.0 / 10 + 4 * vs.hleg / 10
    yb[3] = yta - 3 * vs.hleg / 10
    xb[4] = xta - 2 * vs.hleg / 10
    yb[4] = yta - 3 * vs.hleg / 10
    xb[5] = xta - 2 * vs.hleg / 10
    yb[5] = yta - 14 * vs.hleg / 10
    w = wd
    poly2!(gs, io, xb, yb, 5, -w, backgr)

    w = vs.ithick * wt
    if vs.icon >= 0
        curv2!(gs, io, xc, yc, 2, vs.idash, w, 0, 0, 0.0, 0)
    end
    ssym = 3 * vs.hleg / 10
    u, v = trans3(gs, gs.leg_x + 1.0 / 4, gs.leg_y - 7 * vs.hleg / 10, 0.0)
    if vs.icon != 0
        dsym!(gs, io, u, v, vs.isym, ssym)
    end
    text3!(gs, io, vs.aleg, vs.nleg, vs.istyle, vs.hleg,
           gs.leg_x + 7.0 / 10, gs.leg_y - 11 * vs.hleg / 10, 0.0,
           1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    gs.leg_y = gs.leg_y - 11 * vs.hleg / 10
end

# ---- Curve tag ----

"""
    tagit!(gs, io, x, y, n, it, vs)

Place a tag string at xtag,ytag and draw a vector to the curve at xpoint.
Matches viewr.f90 `tagit` (lines 1179-1247).
"""
function tagit!(gs::GraphState, io::IO,
                x::Vector{Float64}, y::Vector{Float64},
                n::Int, it::Int, vs::ViewrState)
    wt = vs.wline
    backgr = 0.0
    wd = 0.005
    ivec = 1021

    # search curve for xpoint, ypoint
    ypoint = 0.0
    if vs.xpoint > 0.0
        j = 1
        while x[j] < vs.xpoint && j < n
            j = j + 1
        end
        if j == 1; j = 2; end
        jt = it + 1
        ypoint = terp1(x[j-1], y[j-1], x[j], y[j], vs.xpoint, jt)
    end

    # write tag (with a blanked background) and draw vector
    xtp = xscale(gs, vs.xtag)
    ytp = yscale(gs, vs.ytag)
    xta, yta = trans3(gs, xtp, ytp, 0.0)
    wleg = txtlen(gs, vs.aleg, vs.nleg, vs.istyle, vs.hleg)
    xb = Vector{Float64}(undef, 5)
    yb = Vector{Float64}(undef, 5)
    xb[1] = xta - 2 * vs.hleg / 10
    yb[1] = yta - 2 * vs.hleg / 10
    xb[2] = xta + wleg + 4 * vs.hleg / 10
    yb[2] = yta - 2 * vs.hleg / 10
    xb[3] = xta + wleg + 4 * vs.hleg / 10
    yb[3] = yta + 12 * vs.hleg / 10
    xb[4] = xta - 2 * vs.hleg / 10
    yb[4] = yta + 12 * vs.hleg / 10
    xb[5] = xta - 2 * vs.hleg / 10
    yb[5] = yta - 2 * vs.hleg / 10
    w = wd
    poly2!(gs, io, xb, yb, 5, -w, backgr)
    text3!(gs, io, vs.aleg, vs.nleg, vs.istyle, vs.hleg,
           xtp, ytp, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0)

    if vs.xpoint <= 0.0; return; end
    xpp = xscale(gs, vs.xpoint)
    ypp = yscale(gs, ypoint)
    xxx = xtp - vs.hleg / 2
    yyy = ytp + vs.hleg / 2
    vect3!(gs, io, xxx, yyy, 0.0, xpp, ypp, 0.0, wt, ivec)
end

# ---- 3D surface plot ----

"""
    set3d!(gs, io, iplot, xyz, nxyz, vs,
           xll, yll, ww, wh, wr)

Set up 3-dimensional surface plots.
Matches viewr.f90 `set3d` (lines 1249-1524).
"""
function set3d!(gs::GraphState, io::IO, iplot::Int,
                xyz::Vector{Float64}, nxyz::Int, vs::ViewrState,
                xll::Float64, yll::Float64, ww::Float64,
                wh::Float64, wr::Float64)
    wt = vs.wline
    big = 1.0e10
    d0 = 0.001
    d3 = 0.301
    d7 = 0.699
    kmax = 1999
    ten = 10.0

    # set up the window and titles
    if vs.csize > 0.0; vs.hlab = vs.csize; end
    if vs.csize < 0.0; vs.hlab = -vs.csize * wh; end
    if iplot == 1
        initp!(gs, io, vs.lori, vs.xpage, vs.ypage,
               vs.istyle, vs.hlab, vs.wline, 0, vs.ipcol)
    end
    window!(gs, io, xll, yll, ww, wh, wr,
            vs.t1, vs.n1, vs.t2, vs.n2, 0)

    # set up the 3-d perspective projection
    init3!(gs, io, vs.x3, vs.y3, vs.z3, vs.xv, vs.yv, vs.zv, vs.iwcol)

    # locate the different curves in the family
    # find lo and hi limits for each axis
    xlo = big; xhi = -big
    ylo = big; yhi = -big
    zlo = big; zhi = -big
    length_max = 2000
    lll = Vector{Int}(undef, length_max)

    i = 1
    j = 1
    nn = 1
    while nn != 0
        lll[j] = i
        yn = xyz[i]
        nn = round(Int, xyz[i+1])
        if nn != 0
            if yn < ylo; ylo = yn; end
            if yn > yhi; yhi = yn; end
            for k in 1:nn
                xn = xyz[i + 2*k]
                if xn < xlo; xlo = xn; end
                if xn > xhi; xhi = xn; end
                zn = xyz[i + 2*k + 1]
                if zn < zlo; zlo = zn; end
                if zn > zhi; zhi = zn; end
            end
            i = i + 2 + 2 * nn
            j = j + 1
            if j > length_max
                error("set3d: array overflow, increase the length parameter")
            end
        end
    end
    ncurv = j - 1

    # choose z axis limits and step sizes
    if vs.zstp == 0.0
        vs.zmin = zlo
        vs.zmax = zhi
        if vs.jtype == 1
            vs.zmin, vs.zmax, major, _ = ascalv(2, vs.zmin, vs.zmax)
            vs.zstp = (vs.zmax - vs.zmin) / major
        else
            if vs.zmax <= 0.0; vs.zmax = 1.0; end
            if vs.zmin <= 0.0; vs.zmin = vs.zmax / 1e8; end
            top = log10(vs.zmax)
            itop = trunc(Int, top)
            if top < 0.0; itop = itop - 1; end
            vs.zmax = Float64(itop + 1)
            if top <= itop + d7 + d0; vs.zmax = Float64(itop) + d7; end
            if top <= itop + d3 + d0; vs.zmax = Float64(itop) + d3; end
            if top <= itop + d0; vs.zmax = Float64(itop); end
            vs.zmax = ten^vs.zmax
            bot = log10(vs.zmin)
            ibot = trunc(Int, bot)
            if bot < 0.0; ibot = ibot - 1; end
            vs.zmin = Float64(ibot)
            if bot >= ibot + d3 - d0; vs.zmin = Float64(ibot) + d3; end
            if bot >= ibot + d7 - d0; vs.zmin = Float64(ibot) + d7; end
            if bot >= ibot + 1 - d0; vs.zmin = Float64(ibot + 1); end
            vs.zmin = ten^vs.zmin
            if vs.zmax / vs.zmin > ten^13
                vs.zmin = ten^(itop - 13)
            end
        end
    end

    # choose y axis limits and step sizes
    if vs.ystp == 0.0
        vs.ymin = ylo
        vs.ymax = yhi
        if vs.itype == 1 || vs.itype == 3
            vs.ymin, vs.ymax, major, _ = ascalv(2, vs.ymin, vs.ymax)
            vs.ystp = (vs.ymax - vs.ymin) / major
        else
            if vs.ymax <= 0.0; vs.ymax = 1.0; end
            if vs.ymin <= 0.0; vs.ymin = vs.ymax / 1e8; end
            top = log10(vs.ymax)
            itop = trunc(Int, top)
            if top < 0.0; itop = itop - 1; end
            vs.ymax = Float64(itop + 1)
            if top <= itop + d7 + d0; vs.ymax = Float64(itop) + d7; end
            if top <= itop + d3 + d0; vs.ymax = Float64(itop) + d3; end
            if top <= itop + d0; vs.ymax = Float64(itop); end
            vs.ymax = ten^vs.ymax
            bot = log10(vs.ymin)
            ibot = trunc(Int, bot)
            if bot < 0.0; ibot = ibot - 1; end
            vs.ymin = Float64(ibot)
            if bot >= ibot + d3 - d0; vs.ymin = Float64(ibot) + d3; end
            if bot >= ibot + d7 - d0; vs.ymin = Float64(ibot) + d7; end
            if bot >= ibot + 1 - d0; vs.ymin = Float64(ibot + 1); end
            vs.ymin = ten^vs.ymin
            if vs.ymax / vs.ymin > ten^13
                vs.ymin = ten^(itop - 13)
            end
        end
    end

    # choose x axis limits and step sizes
    if vs.xstp == 0.0
        vs.xmin = xlo
        vs.xmax = xhi
        if vs.itype == 1 || vs.itype == 2
            vs.xmin, vs.xmax, major, _ = ascalv(2, vs.xmin, vs.xmax)
            vs.xstp = (vs.xmax - vs.xmin) / major
        else
            if vs.xmax <= 0.0; vs.xmax = 1.0; end
            if vs.xmin <= 0.0; vs.xmin = vs.xmax / 1e8; end
            top = log10(vs.xmax)
            itop = trunc(Int, top)
            if top < 0.0; itop = itop - 1; end
            vs.xmax = Float64(itop + 1)
            if top <= itop + d7 + d0; vs.xmax = Float64(itop) + d7; end
            if top <= itop + d3 + d0; vs.xmax = Float64(itop) + d3; end
            if top <= itop + d0; vs.xmax = Float64(itop); end
            vs.xmax = ten^vs.xmax
            bot = log10(vs.xmin)
            ibot = trunc(Int, bot)
            if bot < 0.0; ibot = ibot - 1; end
            vs.xmin = Float64(ibot)
            if bot >= ibot + d3 - d0; vs.xmin = Float64(ibot) + d3; end
            if bot >= ibot + d7 - d0; vs.xmin = Float64(ibot) + d7; end
            if bot >= ibot + 1 - d0; vs.xmin = Float64(ibot + 1); end
            vs.xmin = ten^vs.xmin
            if vs.xmax / vs.xmin > ten^13
                vs.xmin = ten^(itop - 13)
            end
        end
    end

    # draw 3d z axis (uses 2d right axis parameters)
    zop = vs.zstp
    if vs.jtype == 2; zop = 0.0; end
    if vs.x3 > 0.0
        axis3!(gs, io, vs.zmin, vs.zmax, zop, vs.rl, vs.nr, 1, 1,
               vs.z3, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1)
    else
        axis3!(gs, io, vs.zmin, vs.zmax, zop, vs.rl, vs.nr, 1, 1,
               vs.z3, vs.x3, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1)
    end

    # draw 3d x axis
    xop = vs.xstp
    if vs.itype == 3 || vs.itype == 4; xop = 0.0; end
    axis3!(gs, io, vs.xmin, vs.xmax, xop, vs.xl, -vs.nx, -1, 0,
           vs.x3, 0.0, 0.0, 0.0,
           1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0)

    # draw 3d y axis
    yop = vs.ystp
    if vs.itype == 2 || vs.itype == 4; yop = 0.0; end
    if vs.x3 > 0.0
        axis3!(gs, io, vs.ymin, vs.ymax, yop, vs.yl, -vs.ny, -1, 0,
               vs.y3, vs.x3, 0.0, 0.0,
               0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0)
    else
        axis3!(gs, io, vs.ymin, vs.ymax, yop, vs.yl, -vs.ny, -1, 0,
               vs.y3, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0)
    end

    # draw extra axes to complete the 3-d frame
    if vs.x3 > 0.0
        axis3!(gs, io, vs.ymin, vs.ymax, yop, " ", 0, 0, 0,
               vs.y3, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0)
    else
        axis3!(gs, io, vs.ymin, vs.ymax, yop, " ", 0, 0, 0,
               vs.y3, vs.x3, 0.0, 0.0,
               0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0)
    end
    if vs.x3 > 0.0
        axis3!(gs, io, vs.ymin, vs.ymax, yop, " ", 0, 0, 0,
               vs.y3, 0.0, 0.0, vs.z3,
               0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0)
    else
        axis3!(gs, io, vs.ymin, vs.ymax, yop, " ", 0, 0, 0,
               vs.y3, vs.x3, 0.0, vs.z3,
               0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0)
    end
    if vs.x3 > 0.0
        axis3!(gs, io, vs.zmin, vs.zmax, zop, " ", 0, 0, 0,
               vs.z3, 0.0, vs.y3, 0.0,
               0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0)
    else
        axis3!(gs, io, vs.zmin, vs.zmax, zop, " ", 0, 0, 0,
               vs.z3, vs.x3, vs.y3, 0.0,
               0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0)
    end
    if vs.x3 > 0.0
        axis3!(gs, io, vs.zmin, vs.zmax, zop, " ", 0, 0, 0,
               vs.z3, vs.x3, vs.y3, 0.0,
               0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0)
    else
        axis3!(gs, io, vs.zmin, vs.zmax, zop, " ", 0, 0, 0,
               vs.z3, 0.0, vs.y3, 0.0,
               0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0)
    end
    axis3!(gs, io, vs.xmin, vs.xmax, xop, " ", 0, 0, 0,
           vs.x3, 0.0, vs.y3, 0.0,
           1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0)
    axis3!(gs, io, vs.xmin, vs.xmax, xop, " ", 0, 0, 0,
           vs.x3, 0.0, vs.y3, vs.z3,
           1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0)

    # add the base-plane grids if requested
    if vs.igrid == 1
        grid3!(gs, io, 1, 1, 0)
        if vs.x3 > 0.0
            grid3!(gs, io, 0, 1, 1)
        else
            grid3!(gs, io, -1, 1, 1)
        end
        grid3!(gs, io, 1, -1, 1)
    end

    # draw the curves back-to-front
    x3d = Vector{Float64}(undef, 2000)
    y3d = Vector{Float64}(undef, 2000)
    z3d = Vector{Float64}(undef, 2000)
    for jj in 1:ncurv
        l = lll[ncurv - jj + 1]
        yy = xyz[l]
        if yy >= vs.ymin && yy <= vs.ymax
            nn = round(Int, xyz[l + 1])
            x3d[1] = xyz[l + 2]
            y3d[1] = yy
            z3d[1] = vs.zmin
            k = 2
            for ii in 1:nn
                if k < kmax
                    x3d[k] = xyz[l + 2*ii]
                    if x3d[k] < vs.xmin; x3d[k] = vs.xmin; end
                    if x3d[k] > vs.xmax; x3d[k] = vs.xmax; end
                    y3d[k] = yy
                    z3d[k] = xyz[l + 2*ii + 1]
                    if z3d[k] < vs.zmin; z3d[k] = vs.zmin; end
                    if z3d[k] > vs.zmax; z3d[k] = vs.zmax; end
                    k = k + 1
                end
            end
            if k >= kmax
                @warn "set3d: curve truncated"
            end
            x3d[k] = xyz[l + 2*nn]
            y3d[k] = yy
            z3d[k] = vs.zmin
            curv3!(gs, io, x3d[1:k], y3d[1:k], z3d[1:k], k, wt)
        end
    end
end
