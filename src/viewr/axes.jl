#=
    Axis drawing and coordinate scale functions for the NJOY graph engine.
    Faithful port of graph.f90 subroutines: xscale, yscale, zscale,
    xinvrs, yinvrs, axis2, axis3 (lines 638-1216).

    All numerical constants match Fortran exactly for bit-identical output.
=#

using Printf

# ---- Coordinate scale functions ----

"""
    xscale(gs, x) -> Float64

Convert user x coordinate to axis coordinate (0 to asize range).
Uses fvx/dvx/logx state set by axis3!.
Matches graph.f90 `xscale` (lines 1141-1154).
"""
function xscale(gs::GraphState, x::Float64)
    xn = x
    if gs.logx == 1
        xn = xn > 0.0 ? log10(xn) : -1000.0
    end
    return (xn - gs.fvx) / gs.dvx
end

"""
    yscale(gs, y) -> Float64

Convert user y coordinate to axis coordinate (0 to asize range).
Uses fvy/dvy/logy state set by axis3!.
Matches graph.f90 `yscale` (lines 1156-1169).
"""
function yscale(gs::GraphState, y::Float64)
    yn = y
    if gs.logy == 1
        yn = yn > 0.0 ? log10(yn) : -1000.0
    end
    return (yn - gs.fvy) / gs.dvy
end

"""
    zscale(gs, z) -> Float64

Convert user z coordinate to axis coordinate (0 to asize range).
Uses fvz/dvz/logz state set by axis3!.
Matches graph.f90 `zscale` (lines 1171-1184).
"""
function zscale(gs::GraphState, z::Float64)
    zn = z
    if gs.logz == 1
        zn = zn > 0.0 ? log10(zn) : -1000.0
    end
    return (zn - gs.fvz) / gs.dvz
end

"""
    xinvrs(gs, x) -> Float64

Convert x axis coordinate back to user coordinate.
Inverse of xscale.
Matches graph.f90 `xinvrs` (lines 1186-1200).
"""
function xinvrs(gs::GraphState, x::Float64)
    xn = x * gs.dvx + gs.fvx
    if gs.logx == 1
        return 10.0^xn
    end
    return xn
end

"""
    yinvrs(gs, y) -> Float64

Convert y axis coordinate back to user coordinate.
Inverse of yscale.
Matches graph.f90 `yinvrs` (lines 1202-1216).
"""
function yinvrs(gs::GraphState, y::Float64)
    yn = y * gs.dvy + gs.fvy
    if gs.logy == 1
        return 10.0^yn
    end
    return yn
end

# ---- Number formatting helpers ----
# stripv is defined in text.jl and used here without redefinition.

"""
    format_linear_number(v, scale, ifracs) -> (num::String, lnum::Int)

Format a linear axis number using Fortran-matching format strings.
If ifracs == 1, uses f6.1 format; otherwise uses i5 format.
After formatting, strips leading/trailing blanks via stripv.
"""
function format_linear_number(v::Float64, scale::Float64, ifracs::Int)
    if ifracs == 1
        vv = v / scale
        num = @sprintf("%6.1f", vv)
    else
        iv = round(Int, v / scale)
        num = @sprintf("%5d", iv)
    end
    return stripv(num)
end

"""
    format_scale_suffix(nscale) -> (num::String, lnum::Int)

Format the "*10^n" scale annotation for linear axes.
Uses i1/i2/i3 width depending on nscale magnitude.
Matches graph.f90 axis3 (lines 831-840).
"""
function format_scale_suffix(nscale::Int)
    if (nscale > -10 && nscale < 0) || nscale >= 10
        num = @sprintf("*10#EH.8<%2d#HXEX<", nscale)
        lnum = 17
    elseif nscale <= -10
        num = @sprintf("*10#EH.8<%3d#HXEX<", nscale)
        lnum = 18
    else
        num = @sprintf("*10#EH.8<%1d#HXEX<", nscale)
        lnum = 16
    end
    return (num, lnum)
end

"""
    format_log_number(i) -> (num::String, lnum::Int)

Format a log-scale axis number as "10^n" with superscript markup.
Uses i1/i2/i3 width depending on exponent magnitude.
Matches graph.f90 axis3 log-number formatting.
"""
function format_log_number(i::Int)
    if (i > -10 && i < 0) || i >= 10
        num = @sprintf("10#EH.8<%2d#HXEX<", i)
        lnum = 16
    elseif i <= -10
        num = @sprintf("10#EH.8<%3d#HXEX<", i)
        lnum = 17
    else
        num = @sprintf("10#EH.8<%1d#HXEX<", i)
        lnum = 15
    end
    return (num, lnum)
end

# ---- Axis drawing routines ----

"""
    axis2!(gs, io, amin, amax, astp, label, nlabel, itic, nonum,
           asize, xo, yo, xx, yx, xy, yy, numr)

Draw a 2D axis by delegating to axis3! with z components set to zero.
Matches graph.f90 `axis2` (lines 638-653).
"""
function axis2!(gs::GraphState, io::IO,
                amin::Float64, amax::Float64, astp::Float64,
                label::String, nlabel::Int, itic::Int, nonum::Int,
                asize::Float64, xo::Float64, yo::Float64,
                xx::Float64, yx::Float64, xy::Float64, yy::Float64,
                numr::Int)
    axis3!(gs, io, amin, amax, astp, label, nlabel, itic, nonum,
           asize, xo, yo, 0.0, xx, yx, 0.0, xy, yy, 0.0, numr)
end

"""
    axis3!(gs, io, amin, amax, astp, label, nlabel, itic, nonum,
           asize, xo, yo, zo, xx, yx, zx, xy, yy, zy, numr)

Draw an axis in 3D perspective with tick marks, numbers, and label.

Arguments:
- amin/amax: axis data range
- astp: step between numbered ticks (0 = log scale)
- label: axis label text ("." = no label)
- nlabel: >0 labels above, <0 labels below, 0 no labels/numbers
- itic: >0 ticks above, <0 ticks below, 0 no ticks
- nonum: 0=all numbered, 1=skip first, 2=skip last, 3=skip both ends
- asize: axis length in absolute units
- xo,yo,zo: axis origin in absolute coordinates
- xx,yx,zx: unit vector along axis direction
- xy,yy,zy: unit vector for character vertical direction
- numr: 0=horizontal numbers, 1=rotated perpendicular

Sets gs.fvx/dvx/logx (or fvy/dvy/logy or fvz/dvz/logz) as side effect.
Matches graph.f90 `axis3` (lines 655-1139).
"""
function axis3!(gs::GraphState, io::IO,
                amin::Float64, amax::Float64, astp::Float64,
                label::String, nlabel::Int, itic::Int, nonum::Int,
                asize::Float64,
                xo::Float64, yo::Float64, zo::Float64,
                xx::Float64, yx::Float64, zx::Float64,
                xy::Float64, yy::Float64, zy::Float64,
                numr::Int)
    # Fortran constants
    ds   = 0.01
    sc1  = 0.099
    sc2  = 0.901

    # Log subdivision table (8 entries, no trailing 1.0)
    subs = (0.301, 0.477, 0.602, 0.699, 0.778, 0.845, 0.903, 0.954)

    # Set color to black
    gs.ifg = 1

    # Draw the axis vector
    xe = xo + asize * xx
    ye = yo + asize * yx
    ze = zo + asize * zx
    vect3!(gs, io, xo, yo, zo, xe, ye, ze, gs.wa, 0)

    if astp != 0.0
        # ---- Linear scale ----
        _axis3_linear!(gs, io, amin, amax, astp, label, nlabel, itic, nonum,
                       asize, xo, yo, zo, xx, yx, zx, xy, yy, zy, numr,
                       subs, ds, sc1, sc2)
    else
        # ---- Log scale ----
        _axis3_log!(gs, io, amin, amax, label, nlabel, itic, nonum,
                    asize, xo, yo, zo, xx, yx, zx, xy, yy, zy, numr,
                    subs, ds)
    end
end

"""
Internal: linear scale branch of axis3.
"""
function _axis3_linear!(gs::GraphState, io::IO,
                        amin::Float64, amax::Float64, astp::Float64,
                        label::String, nlabel::Int, itic::Int, nonum::Int,
                        asize::Float64,
                        xo::Float64, yo::Float64, zo::Float64,
                        xx::Float64, yx::Float64, zx::Float64,
                        xy::Float64, yy::Float64, zy::Float64,
                        numr::Int,
                        subs, ds, sc1, sc2)
    n = round(Int, (amax - amin) / astp)
    step = asize / n
    n = n + 1
    test = 1.0 / 1000.0

    # Detect which axis and set coordinate scale state
    if abs(xx - 1) < test
        gs.fvx = amin
        gs.dvx = (amax - amin) / asize
        gs.logx = 0
        gs.xmin = xo
        gs.xmax = xo + asize
        gs.xstp = step
    elseif abs(yx - 1) < test
        gs.fvy = amin
        gs.dvy = (amax - amin) / asize
        gs.logy = 0
        gs.ymin = yo
        gs.ymax = yo + asize
        gs.ystp = step
    elseif abs(zx - 1) < test
        gs.fvz = amin
        gs.dvz = (amax - amin) / asize
        gs.logz = 0
        gs.zmin = zo
        gs.zmax = zo + asize
        gs.zstp = step
    end

    x = xo
    y = yo
    z = zo
    v = amin
    shl = 0.0

    # Compute engineering notation scale
    nscale = 0
    if abs(astp) < sc1 || abs(astp) > sc2
        nscale = round(Int, log10(abs(astp)))
        if nscale < 0
            nscale = nscale - 2
        end
        nscale = 3 * div(nscale, 3)
    end
    scale = 10.0^nscale

    # Determine integer vs fractional formatting
    iv = round(Int, astp / scale)
    vv = astp / scale
    test = 1.0 / 100.0
    ifracs = 0
    if abs(vv - iv) > test
        ifracs = 1
    end

    # Track xc/yc/zc for scale suffix positioning
    xc = x
    yc = y
    zc = z

    # Draw tick marks and numbers for each step
    for i in 1:n
        # Draw tick mark
        if itic < 0
            vect3!(gs, io, x, y, z,
                   x - gs.tic * xy, y - gs.tic * yy, z - gs.tic * zy,
                   gs.wt, 0)
        elseif itic > 0
            vect3!(gs, io, x, y, z,
                   x + gs.tic * xy, y + gs.tic * yy, z + gs.tic * zy,
                   gs.wt, 0)
        end

        # Draw number label
        if nlabel != 0
            if i != 1 || (nonum != 1 && nonum != 3)
                num, lnum = format_linear_number(v, scale, ifracs)
                ww = txtlen(gs, num, lnum, gs.lfont, gs.hn) / 2
                wwt = ww

                # Adjust width at endpoints for alignment
                if i == 1 && itic > 0 && nlabel < 0
                    wwt = ww / 2
                end
                if i == 1 && itic < 0 && nlabel > 0
                    wwt = ww / 2
                end
                if i == n && itic > 0 && nlabel < 0
                    wwt = 3 * ww / 2
                end
                if i == n && itic < 0 && nlabel > 0
                    wwt = 3 * ww / 2
                end

                if nlabel < 0
                    sh = gs.gap + gs.hf * gs.hn
                    if itic < 0
                        sh = sh + gs.tic - gs.gap / 2
                    end
                    if numr == 0
                        xc = x - sh * xy - wwt * xx
                        yc = y - sh * yy - wwt * yx
                        zc = z - sh * zy - wwt * zx
                        if sh + gs.hf * gs.hn > shl
                            shl = sh + gs.hf * gs.hn
                        end
                    else
                        xc = x - sh * xy - gs.hf * gs.hn * xx / 2
                        yc = y - sh * yy - gs.hf * gs.hn * yx / 2
                        zc = z - sh * zy - gs.hf * gs.hn * zx / 2
                        test = sh + 2 * ww + gs.hn / 2
                        if test > shl
                            shl = test
                        end
                    end
                else
                    sh = gs.gap
                    if itic > 0
                        sh = sh + gs.tic - gs.gap / 2
                    end
                    if numr == 0
                        xc = x + sh * xy - wwt * xx
                        yc = y + sh * yy - wwt * yx
                        zc = z + sh * zy - wwt * zx
                        test = sh + 2 * gs.hf * gs.hn
                        if test > shl
                            shl = test
                        end
                    else
                        xc = x + (sh + 2 * ww) * xy - gs.hf * gs.hn * xx / 2
                        yc = y + (sh + 2 * ww) * yy - gs.hf * gs.hn * yx / 2
                        zc = z + (sh + 2 * ww) * zy - gs.hf * gs.hn * zx / 2
                        test = sh + 2 * ww + gs.hn / 2
                        if test > shl
                            shl = test
                        end
                    end
                end

                if numr == 0
                    text3!(gs, io, num, lnum, gs.lfont, gs.hn,
                           xc, yc, zc, xx, yx, zx, xy, yy, zy)
                else
                    text3!(gs, io, num, lnum, gs.lfont, gs.hn,
                           xc, yc, zc, -xy, -yy, -zy, xx, yx, zx)
                end
            end
        end

        x = x + step * xx
        y = y + step * yx
        z = z + step * zx
        v = v + astp
    end

    # Draw scale annotation "*10^n" if needed
    if nscale != 0 && nlabel != 0
        if nonum != 2 && nonum != 3
            num, lnum = format_scale_suffix(nscale)
            www = txtlen(gs, num, lnum, gs.lfont, gs.hn)

            if nlabel < 0
                if numr == 0
                    xc = xc - 6 * gs.hn * xy / 5
                    yc = yc - 6 * gs.hn * yy / 5
                    zc = zc - 6 * gs.hn * zy / 5
                else
                    xc = xc - 6 * gs.hn * xx / 5
                    yc = yc - 6 * gs.hn * yx / 5
                    zc = zc - 6 * gs.hn * zx / 5
                end
            else
                if numr == 0
                    xc = xc + 6 * gs.hn * xy / 5
                    yc = yc + 6 * gs.hn * yy / 5
                    zc = zc + 6 * gs.hn * zy / 5
                else
                    dd = gs.hf * gs.hn / 2 + 6 * gs.hn / 5
                    xc = x + (www + gs.gap) * xy - dd * xx - step * xx
                    yc = y + (www + gs.gap) * yy - dd * yx - step * yx
                    zc = z + (www + gs.gap) * zy - dd * zx - step * zx
                end
            end

            if numr == 0
                text3!(gs, io, num, lnum, gs.lfont, gs.hn,
                       xc, yc, zc, xx, yx, zx, xy, yy, zy)
            else
                text3!(gs, io, num, lnum, gs.lfont, gs.hn,
                       xc, yc, zc, -xy, -yy, -zy, xx, yx, zx)
            end
        end
    end

    # Add centered axis label
    _axis3_label!(gs, io, label, nlabel, asize, xo, yo, zo,
                  xx, yx, zx, xy, yy, zy, shl)
end

"""
Internal: log scale branch of axis3.
"""
function _axis3_log!(gs::GraphState, io::IO,
                     amin::Float64, amax::Float64,
                     label::String, nlabel::Int, itic::Int, nonum::Int,
                     asize::Float64,
                     xo::Float64, yo::Float64, zo::Float64,
                     xx::Float64, yx::Float64, zx::Float64,
                     xy::Float64, yy::Float64, zy::Float64,
                     numr::Int,
                     subs, ds)
    shl = 0.0

    # Guard: log axes require positive limits
    if amin <= 0.0; amin = 1e-30; end
    if amax <= 0.0; amax = 1.0; end
    if amax <= amin; amax = amin * 10.0; end

    # Compute origen (start of log range, snapped to decades/sub-decades)
    origen = log10(amin)
    i = trunc(Int, origen)
    if origen < 0.0
        i = i - 1
    end
    v = origen - i
    if v < ds
        v = 0.0
    elseif v < subs[1] + ds
        v = subs[1]
    elseif v < subs[4] + ds
        v = subs[4]
    else
        i = i + 1
        v = 0.0
    end
    origen = i + v
    imin = i
    if v > 0.0
        imin = i + 1
    end

    # Measure width of first decade label
    num, lnum = format_log_number(imin)
    ww = txtlen(gs, num, lnum, gs.lfont, gs.hn)
    www = ww

    # Compute aend (end of log range, snapped)
    aend = log10(amax)
    j = trunc(Int, aend)
    if aend < 0.0
        j = j - 1
    end
    v = aend - j
    if v < ds
        v = 0.0
    elseif v < subs[1] + ds
        v = subs[1]
    elseif v < subs[4] + ds
        v = subs[4]
    else
        j = j + 1
        v = 0.0
    end
    aend = j + v
    imax = j
    if v > 0.0
        imax = imax - 1
    end

    # Measure width of last decade label
    num, lnum = format_log_number(imax)
    ww = txtlen(gs, num, lnum, gs.lfont, gs.hn)
    if ww > www
        www = ww
    end

    # Compute cycles and label skip
    cycles = (aend - origen) / asize
    test = 0.1
    if abs(xx * xy + yx * yy + zx * zy - 1) < test
        room = 6 * gs.hn / 5
    else
        room = 6 * www / 5
    end
    iskip = trunc(Int, room * abs(cycles))
    iskip = 1 + iskip

    # Set coordinate scale state
    test = 0.001
    if abs(xx - 1) < test
        gs.fvx = origen
        gs.dvx = cycles
        gs.logx = 1
        gs.xmin = xo
        gs.xmax = xo + asize
    elseif abs(yx - 1) < test
        gs.fvy = origen
        gs.dvy = cycles
        gs.logy = 1
        gs.ymin = yo
        gs.ymax = yo + asize
    elseif abs(zx - 1) < test
        gs.fvz = origen
        gs.dvz = cycles
        gs.logz = 1
        gs.zmin = zo
        gs.zmax = zo + asize
    end

    # Compute tick direction vector
    if itic < 0
        xt = -gs.tic * xy
        yt = -gs.tic * yy
        zt = -gs.tic * zy
    else
        xt = gs.tic * xy
        yt = gs.tic * yy
        zt = gs.tic * zy
    end

    # Start position for decade loop
    i = trunc(Int, origen)
    test = 0.1
    if origen - i > test
        i = i + 1
    end
    v_val = 10.0^i
    x = xo + ((i - origen) / cycles) * xx
    y = yo + ((i - origen) / cycles) * yx
    z = zo + ((i - origen) / cycles) * zx

    # Draw sub-decade ticks before first decade (for < 1 decade span)
    if itic != 0
        test_cyc = 1.0
        if cycles <= test_cyc
            k1 = 0
            test = 0.1
            if abs(origen - i + 1 - subs[1]) < test
                k1 = 1
            end
            if abs(origen - i + 1 - subs[4]) < test
                k1 = 4
            end
            if k1 != 0
                for k in k1:8
                    sh = (subs[k] - 1) / cycles
                    xs = x + sh * xx
                    ys = y + sh * yx
                    zs = z + sh * zx
                    dx = xt / 2
                    dy = yt / 2
                    dz = zt / 2
                    vect3!(gs, io, xs, ys, zs,
                           xs + dx, ys + dy, zs + dz, gs.wt, 0)
                end
            end
        end
    end

    # Main decade loop
    idone = 0
    while idone == 0
        # Draw major tick and sub-decade ticks
        if itic != 0
            vect3!(gs, io, x, y, z, x + xt, y + yt, z + zt, gs.wt, 0)

            test = 0.1
            if abs(aend - i) >= test
                test_cyc = 1.0
                if cycles <= test_cyc
                    k2 = 8
                    test = 0.1
                    if abs(aend - i - subs[1]) < test
                        k2 = 1
                    end
                    if abs(aend - i - subs[4]) < test
                        k2 = 4
                    end
                    for k in 1:k2
                        sh = subs[k] / cycles
                        xs = x + sh * xx
                        ys = y + sh * yx
                        zs = z + sh * zx
                        dx = xt / 2
                        dy = yt / 2
                        dz = zt / 2
                        vect3!(gs, io, xs, ys, zs,
                               xs + dx, ys + dy, zs + dz, gs.wt, 0)
                    end
                end
            end
        end

        # Draw number label
        if nlabel != 0
            test = 0.001
            if abs(i - origen) >= test || (nonum != 1 && nonum != 3)
                if log10(v_val) < aend || (nonum != 2 && nonum != 3)
                    iii = trunc(Int, i - origen)
                    if mod(iii, iskip) == 0
                        num, lnum = format_log_number(i)
                        ww = txtlen(gs, num, lnum, gs.lfont, gs.hn) / 2

                        # Adjust ww at first tick position
                        if i == origen
                            if itic == 0
                                ww = ww / 2
                            else
                                if itic < 0 && nlabel > 0
                                    ww = ww / 2
                                elseif itic > 0 && nlabel < 0
                                    ww = ww / 2
                                end
                            end
                        end

                        if nlabel < 0
                            sh = gs.gap + gs.hf * gs.hn
                            if itic < 0
                                sh = sh + gs.tic - gs.gap / 2
                            end
                            if numr == 0
                                xc = x - sh * xy - ww * xx
                                yc = y - sh * yy - ww * yx
                                zc = z - sh * zy - ww * zx
                                test = sh + gs.hf * gs.hn
                                if test > shl
                                    shl = test
                                end
                            else
                                xc = x - sh * xy - gs.hf * gs.hn * xx / 2
                                yc = y - sh * yy - gs.hf * gs.hn * yx / 2
                                zc = z - sh * zy - gs.hf * gs.hn * zx / 2
                                test = sh + www + gs.hn / 2
                                if test > shl
                                    shl = test
                                end
                            end
                        else
                            sh = gs.gap
                            if itic > 0
                                sh = sh + gs.tic - gs.gap / 2
                            end
                            if numr == 0
                                xc = x + sh * xy - ww * xx
                                yc = y + sh * yy - ww * yx
                                zc = z + sh * zy - ww * zx
                                test = sh + 2 * gs.hf * gs.hn
                                if test > shl
                                    shl = test
                                end
                            else
                                xc = x + (sh + www) * xy - gs.hf * gs.hn * xx / 2
                                yc = y + (sh + www) * yy - gs.hf * gs.hn * yx / 2
                                zc = z + (sh + www) * zy - gs.hf * gs.hn * zx / 2
                                test = sh + www + gs.hn / 2
                                if test > shl
                                    shl = test
                                end
                            end
                        end

                        if numr == 0
                            text3!(gs, io, num, lnum, gs.lfont, gs.hn,
                                   xc, yc, zc, xx, yx, zx, xy, yy, zy)
                        else
                            text3!(gs, io, num, lnum, gs.lfont, gs.hn,
                                   xc, yc, zc, -xy, -yy, -zy, xx, yx, zx)
                        end
                    end
                end
            end
        end

        # Advance to next decade or finish
        if i >= j
            idone = 1
        else
            x = x + (1 / cycles) * xx
            y = y + (1 / cycles) * yx
            z = z + (1 / cycles) * zx
            i = i + 1
            v_val = v_val * 10
        end
    end

    # Add centered axis label
    _axis3_label!(gs, io, label, nlabel, asize, xo, yo, zo,
                  xx, yx, zx, xy, yy, zy, shl)
end

"""
Internal: draw centered axis label (shared by linear and log branches).
Matches graph.f90 axis3 (lines 1118-1137).
"""
function _axis3_label!(gs::GraphState, io::IO,
                       label::String, nlabel::Int,
                       asize::Float64,
                       xo::Float64, yo::Float64, zo::Float64,
                       xx::Float64, yx::Float64, zx::Float64,
                       xy::Float64, yy::Float64, zy::Float64,
                       shl::Float64)
    if nlabel != 0
        if label != "."
            ww = txtlen(gs, label, abs(nlabel), gs.lfont, gs.hl)
            ww = (asize - ww) / 2
            if nlabel < 0
                sh = shl + gs.hf * gs.hl
                xl = xo - sh * xy + ww * xx
                yl = yo - sh * yy + ww * yx
                zl = zo - sh * zy + ww * zx
            else
                sh = shl
                xl = xo + sh * xy + ww * xx
                yl = yo + sh * yy + ww * yx
                zl = zo + sh * zy + ww * zx
            end
            text3!(gs, io, label, abs(nlabel), gs.lfont, gs.hl,
                   xl, yl, zl, xx, yx, zx, xy, yy, zy)
        end
    end
end
