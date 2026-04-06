#=
    Main viewr rendering loop — reads an NJOY plot tape and renders PostScript.
    Faithful port of viewr.f90 `viewr` subroutine (lines 42-755).

    All control flow and data reading patterns match the Fortran exactly.
=#

using Printf

"""
    viewr_render!(io, plot_tape_path)

Read an NJOY plot tape file and render PostScript output to `io`.
Matches the main `viewr` subroutine in viewr.f90 (lines 42-755).
"""
function viewr_render!(io::IO, plot_tape_path::String)
    lines = readlines(plot_tape_path)
    isempty(lines) && return

    gs = GraphState()
    vs = ViewrState()
    mmax = 20000
    maxaa = 500000
    mxc = 80
    flag = -99.0

    # default viewpoint and workbox
    dxv = 15.0; dyv = -15.0; dzv = 15.0
    dx3 = 2.5; dy3 = 6.5; dz3 = 2.5

    # ---- Card 1: page setup ----
    li = 1
    c1, li = _read_card_reals(lines, li, 4,
                              [1.0, 2.0, Float64(DSIZE), 0.0])
    vs.lori = round(Int, c1[1])
    vs.istyle = round(Int, c1[2])
    vs.csize = c1[3]
    vs.ipcol = round(Int, c1[4])

    # default page size
    if vs.lori == 1
        vs.xpage = YPAPER - YMARG
        vs.ypage = XPAPER - XMARG
    else
        vs.xpage = XPAPER - XMARG
        vs.ypage = YPAPER - YMARG
    end
    vs.wline = DLINE

    # initialize plotting
    gplot!(gs, io, vs.lori, vs.xpage, vs.ypage)

    # ---- Loop over plots ----
    nplot = 0
    i3d = 0
    iplot = 1

    while iplot != 99 && li <= length(lines)
        # defaults for card 2
        iplot = 1
        iwcol = 0
        factx = 1.0
        facty = 1.0
        xll = 0.0
        yll = 0.0
        ww = vs.xpage
        wh = vs.ypage
        wr = 0.0

        # read card 2
        c2, li = _read_card_reals(lines, li, 9,
                                  [1.0, 0.0, 1.0, 1.0, 0.0, 0.0,
                                   vs.xpage, vs.ypage, 0.0])
        iplot = round(Int, c2[1])
        iwcol = round(Int, c2[2])
        factx = c2[3]
        facty = c2[4]
        xll = c2[5]
        yll = c2[6]
        ww = c2[7]
        wh = c2[8]
        wr = c2[9]
        vs.iwcol = iwcol
        vs.factx = factx
        vs.facty = facty

        if iplot == 99; break; end

        ierrb = 0

        # ---- Cards 3-7a: if this is first plot on these axes ----
        if abs(iplot) == 1
            # card 3: title line 1
            text = " "
            if li <= length(lines)
                text, _, li = _read_card_string(lines, li)
            end
            vs.n1 = 0
            for i in 1:min(length(text), mxc)
                if i <= length(text) && text[i] != ' '
                    vs.n1 = i
                end
            end
            vs.t1 = vs.n1 > 0 ? text[1:vs.n1] : " "

            # card 3a: title line 2
            text = " "
            if li <= length(lines)
                text, _, li = _read_card_string(lines, li)
            end
            vs.n2 = 0
            for i in 1:min(length(text), mxc)
                if i <= length(text) && text[i] != ' '
                    vs.n2 = i
                end
            end
            vs.t2 = vs.n2 > 0 ? text[1:vs.n2] : " "

            # card 4: plot type and grids
            vs.itype = 4
            vs.jtype = 0
            vs.igrid = 2
            vs.ileg = 0
            vs.xtag = 0.0
            vs.ytag = 0.0
            c4, li = _read_card_reals(lines, li, 6,
                                      [4.0, 0.0, 2.0, 0.0, 0.0, 0.0])
            vs.itype = round(Int, c4[1])
            vs.jtype = round(Int, c4[2])
            vs.igrid = round(Int, c4[3])
            vs.ileg = round(Int, c4[4])
            vs.xtag = c4[5]
            vs.ytag = c4[6]

            # no more input needed if itype==0
            if vs.itype != 0
                i3d = 0
                if vs.itype <= 0
                    i3d = 1
                    vs.itype = -vs.itype
                end

                # card 5: x-axis limits
                vs.xmin = 0.0; vs.xmax = 0.0; vs.xstp = 0.0
                c5, li = _read_card_reals(lines, li, 3, [0.0, 0.0, 0.0])
                vs.xmin = c5[1]; vs.xmax = c5[2]; vs.xstp = c5[3]
                if (vs.itype == 3 || vs.itype == 4) &&
                   (vs.xmin != 0.0 || vs.xmax != 0.0)
                    vs.xstp = 1.0
                end

                # card 5a: x label
                text = " "
                if li <= length(lines)
                    text, _, li = _read_card_string(lines, li)
                end
                vs.nx = 0
                for i in 1:min(length(text), mxc)
                    if i <= length(text) && text[i] != ' '
                        vs.nx = i
                    end
                end
                vs.xl = vs.nx > 0 ? text[1:vs.nx] : " "

                # card 6: y-axis limits
                vs.ymin = 0.0; vs.ymax = 0.0; vs.ystp = 0.0
                c6, li = _read_card_reals(lines, li, 3, [0.0, 0.0, 0.0])
                vs.ymin = c6[1]; vs.ymax = c6[2]; vs.ystp = c6[3]
                if (vs.itype == 2 || vs.itype == 4) &&
                   (vs.ymin != 0.0 || vs.ymax != 0.0)
                    vs.ystp = 1.0
                end

                # card 6a: y label
                text = " "
                if li <= length(lines)
                    text, _, li = _read_card_string(lines, li)
                end
                vs.ny = 0
                for i in 1:min(length(text), mxc)
                    if i <= length(text) && text[i] != ' '
                        vs.ny = i
                    end
                end
                vs.yl = vs.ny > 0 ? text[1:vs.ny] : " "

                # card 7/7a: alternate y axis or z axis
                if vs.jtype != 0
                    vs.zmin = 0.0; vs.zmax = 0.0; vs.zstp = 0.0
                    c7, li = _read_card_reals(lines, li, 3, [0.0, 0.0, 0.0])
                    vs.zmin = c7[1]; vs.zmax = c7[2]; vs.zstp = c7[3]
                    if vs.jtype == 2 && (vs.zmin != 0.0 || vs.zmax != 0.0)
                        vs.zstp = 1.0
                    end

                    # card 7a: right/z label
                    text = " "
                    if li <= length(lines)
                        text, _, li = _read_card_string(lines, li)
                    end
                    vs.nr = 0
                    for i in 1:min(length(text), mxc)
                        if i <= length(text) && text[i] != ' '
                            vs.nr = i
                        end
                    end
                    vs.rl = vs.nr > 0 ? text[1:vs.nr] : " "
                end
            end
        end

        # ---- Additional input for real plots ----
        if vs.itype != 0

            # card 8: dummy card (always 0/)
            if li <= length(lines)
                li += 1
            end

            # read plotting parameters for next 2d curve
            if i3d == 0
                vs.icon = 0; vs.isym = 0; vs.idash = 0
                vs.iccol = 0; vs.ithick = 1; vs.ishade = 0
                c9, li = _read_card_reals(lines, li, 6,
                                          [0.0, 0.0, 0.0, 0.0, 1.0, 0.0])
                vs.icon = round(Int, c9[1])
                vs.isym = round(Int, c9[2])
                vs.idash = round(Int, c9[3])
                vs.iccol = round(Int, c9[4])
                vs.ithick = round(Int, c9[5])
                vs.ishade = round(Int, c9[6])

                # read legend or tag title lines
                if vs.ileg != 0
                    text = " "
                    if li <= length(lines)
                        text, _, li = _read_card_string(lines, li)
                    end
                    vs.nleg = 0
                    for i in 1:min(length(text), mxc)
                        if i <= length(text) && text[i] != ' '
                            vs.nleg = i
                        end
                    end
                    vs.aleg = vs.nleg > 0 ? text[1:vs.nleg] : " "

                    if vs.ileg == 2
                        vs.xtag = 0.0; vs.ytag = 0.0; vs.xpoint = 0.0
                        c10a, li = _read_card_reals(lines, li, 3,
                                                    [0.0, 0.0, 0.0])
                        vs.xtag = c10a[1]
                        vs.ytag = c10a[2]
                        vs.xpoint = c10a[3]
                    end
                end

            else
                # read parameters for plotting 3d surface
                vs.xv = dxv; vs.yv = dyv; vs.zv = dzv
                vs.x3 = dx3; vs.y3 = dy3; vs.z3 = dz3
                c11, li = _read_card_reals(lines, li, 6,
                                           [dxv, dyv, dzv, dx3, dy3, dz3])
                vs.xv = c11[1]; vs.yv = c11[2]; vs.zv = c11[3]
                vs.x3 = c11[4]; vs.y3 = c11[5]; vs.z3 = c11[6]
            end

            # card 12: determine input type
            nform = 0
            c12, li = _read_card_reals(lines, li, 1, [0.0])
            nform = round(Int, c12[1])

            # Pre-declare arrays for both 2D and 3D paths
            x = Vector{Float64}(undef, mmax)
            y = Vector{Float64}(undef, mmax)
            dxm = Vector{Float64}(undef, mmax)
            dxp = Vector{Float64}(undef, mmax)
            dym = Vector{Float64}(undef, mmax)
            dyp = Vector{Float64}(undef, mmax)
            aa = Vector{Float64}(undef, maxaa)
            n = 0; naa = 0; ierrb = 0

            # ---- format 0: 2d data ----
            if nform == 0
                i = 0
                idone = 0
                while idone == 0 && li <= length(lines)
                    i += 1
                    z = [flag, flag, 0.0, 0.0, 0.0, 0.0]
                    parsed = read_plot_reals(lines[li])
                    li += 1
                    for j in 1:min(length(parsed), 6)
                        z[j] = parsed[j]
                    end
                    zz1 = z[1]; zz2 = z[2]
                    if zz1 == flag && zz2 == flag
                        idone = 1
                    else
                        if i > mmax
                            @warn "viewr: too many data points, truncating"
                            idone = 1
                            i -= 1
                        else
                            x[i] = z[1] * factx
                            y[i] = z[2] * facty
                            dym[i] = z[3] * facty
                            dyp[i] = z[4] * facty
                            if dyp[i] == 0.0; dyp[i] = dym[i]; end
                            dxm[i] = z[5] * factx
                            dxp[i] = z[6] * factx
                            if dxp[i] == 0.0; dxp[i] = dxm[i]; end
                            if dym[i] != 0.0; ierrb = 1; end
                            if dxm[i] != 0.0; ierrb = 2; end
                            if dym[i] != 0.0 && dxm[i] != 0.0; ierrb = 3; end
                        end
                    end
                end
                n = i - 1

            # ---- format 1: 3d data ----
            elseif nform == 1
                l = 0
                iskip = 0
                idone = 0
                while idone == 0 && li <= length(lines)
                    parsed = read_plot_reals(lines[li])
                    li += 1
                    yy = isempty(parsed) ? flag : parsed[1]
                    if yy == flag
                        idone = 1
                    else
                        if iskip == 0
                            if l + 5000 >= maxaa
                                @warn "viewr: too much 3d data, truncating"
                                iskip = 1
                            end
                            l += 1; aa[l] = yy
                            l += 1; aa[l] = 0.0
                            ll = l
                        end
                        nn_pts = 0
                        inside = 0
                        while inside == 0 && li <= length(lines)
                            parsed2 = read_plot_reals(lines[li])
                            li += 1
                            xx = isempty(parsed2) ? flag : parsed2[1]
                            zz = length(parsed2) >= 2 ? parsed2[2] : flag
                            if xx == flag && zz == flag
                                inside = 1
                            else
                                if iskip == 0
                                    nn_pts += 1
                                    l += 1; aa[l] = xx
                                    l += 1; aa[l] = zz
                                end
                            end
                        end
                        if iskip == 0; aa[ll] = Float64(nn_pts); end
                    end
                end
                l += 1; aa[l] = 0.0
                l += 1; aa[l] = 0.0
                naa = l
            end
        end

        # ---- do plot ----
        if nplot != 0 && iplot == 1
            endp!(gs, io)
        end
        if iplot == -1
            endw!(gs, io)
        end
        if vs.itype != 0
            if nform == 0
                set2d!(gs, io, iplot,
                       x[1:n], y[1:n], dxm[1:n], dxp[1:n],
                       dym[1:n], dyp[1:n],
                       n, ierrb, vs, xll, yll, ww, wh, wr)
            elseif naa > 0
                set3d!(gs, io, iplot, aa[1:naa], naa, vs,
                       xll, yll, ww, wh, wr)
            end
        else
            # itype == 0: titles only, still need window setup
            if abs(iplot) == 1
                if iplot == 1
                    if vs.csize > 0.0; vs.hlab = vs.csize; end
                    if vs.csize < 0.0; vs.hlab = -vs.csize * wh; end
                    vs.hleg = 2 * vs.hlab / 3
                    initp!(gs, io, vs.lori, vs.xpage, vs.ypage,
                           vs.istyle, vs.hlab, vs.wline, 0, vs.ipcol)
                end
                window!(gs, io, xll, yll, ww, wh, wr,
                        vs.t1, vs.n1, vs.t2, vs.n2, 0)
            end
        end
        nplot = 1
    end

    # viewr is finished
    endw!(gs, io)
    endp!(gs, io)
    gdone!(gs, io)
end
