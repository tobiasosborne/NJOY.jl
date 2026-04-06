#=
    GraphState — mutable state for the NJOY PostScript rendering engine.
    Encapsulates all module-level variables from Fortran graph.f90 + viewr.f90.
    Passed as first argument to all rendering functions.
=#

"""
    GraphState

All mutable state for the NJOY PostScript graph engine.
Replaces ~60 module-level variables from Fortran graph.f90.
"""
mutable struct GraphState
    # ---- PS file state ----
    ipage::Int          # current page counter
    land::Int           # orientation: 0=portrait, 1=landscape
    ushift::Float64     # PS x-offset (points)
    vshift::Float64     # PS y-offset (points)
    uwidth::Float64     # plotting width (points, for landscape rotation)

    # ---- Page/font state ----
    ifont::Int          # global font: 1=Times-Roman, 2=Helvetica, 3=Symbol
    width::Float64      # base line width (inches)
    ht::Float64         # character height (inches)
    ifg::Int            # foreground color index (into IFRGB, 1-based)
    ibg::Int            # background color index (into IBRGB, 1-based)

    # ---- Tick/axis geometry (set by init2/init3) ----
    wt::Float64         # tick line weight
    wg::Float64         # grid line weight
    wa::Float64         # axis line weight
    tic::Float64        # tick mark length
    gap::Float64        # gap between axis and number text
    hf::Float64         # height fraction for vertical centering (3/4)
    hl::Float64         # label height
    hn::Float64         # number height
    lfont::Int          # font for axis labels/numbers
    tspace::Float64     # vertical space consumed by window title

    # ---- Coordinate scale state (set by axis3, used by xscale/yscale/zscale) ----
    fvx::Float64;  dvx::Float64;  logx::Int
    fvy::Float64;  dvy::Float64;  logy::Int
    fvz::Float64;  dvz::Float64;  logz::Int
    xmin::Float64; xmax::Float64; xstp::Float64
    ymin::Float64; ymax::Float64; ystp::Float64
    zmin::Float64; zmax::Float64; zstp::Float64

    # ---- 3D projection state ----
    ro::Float64         # view distance from origin
    rs::Float64         # composite scale factor
    rp::Float64         # horizontal projection length
    ct::Float64         # cos(elevation)
    st::Float64         # sin(elevation)
    cp::Float64         # cos(azimuth)
    sp::Float64         # sin(azimuth)
    du::Float64         # 2D translation x after projection
    dv::Float64         # 2D translation y after projection

    # ---- Window state ----
    xwll::Float64       # window lower-left x (page coords)
    ywll::Float64       # window lower-left y
    www::Float64        # window width
    wwh::Float64        # window height
    wwr::Float64        # window rotation (degrees)

    # ---- drawh persistent state (Fortran SAVE variables) ----
    wlast::Float64      # last linewidth written to PS
    ldash::Int          # last dash pattern written to PS

    # ---- legndb persistent state (Fortran SAVE variables) ----
    leg_x::Float64      # current legend block x position
    leg_y::Float64      # current legend block y position
end

function GraphState()
    GraphState(
        # PS file state
        0,          # ipage
        0,          # land
        0.0, 0.0,   # ushift, vshift
        0.0,        # uwidth
        # Page/font state
        2,          # ifont (Helvetica default)
        DLINE,      # width
        DSIZE,      # ht
        1,          # ifg (black)
        1,          # ibg (white)
        # Tick/axis geometry
        0.0, 0.0, 0.0,  # wt, wg, wa
        0.0, 0.0,       # tic, gap
        0.75,            # hf = 3/4
        0.0, 0.0,       # hl, hn
        2,               # lfont
        0.0,             # tspace
        # Coordinate scale
        0.0, 1.0, 0,   # fvx, dvx, logx
        0.0, 1.0, 0,   # fvy, dvy, logy
        0.0, 1.0, 0,   # fvz, dvz, logz
        0.0, 0.0, 0.0, # xmin, xmax, xstp
        0.0, 0.0, 0.0, # ymin, ymax, ystp
        0.0, 0.0, 0.0, # zmin, zmax, zstp
        # 3D projection (identity defaults matching endw reset)
        1.0,        # ro
        1.0,        # rs
        0.0,        # rp
        0.0,        # ct
        1.0,        # st
        0.0,        # cp
        -1.0,       # sp
        0.0, 0.0,   # du, dv
        # Window state
        0.0, 0.0,   # xwll, ywll
        0.0, 0.0,   # www, wwh
        0.0,        # wwr
        # drawh persistent state
        0.0,        # wlast
        -1,         # ldash
        # legndb persistent state
        0.0, 0.0,   # leg_x, leg_y
    )
end

"""
    ViewrState

Holds viewr.f90 module-level variables for the plot tape parser.
Separate from GraphState to keep the graph engine reusable.
"""
mutable struct ViewrState
    # ---- Page setup (card 1) ----
    xpage::Float64
    ypage::Float64
    wline::Float64
    csize::Float64
    lori::Int
    istyle::Int
    ipcol::Int

    # ---- Current window (card 2) ----
    iwcol::Int
    factx::Float64
    facty::Float64

    # ---- Labels and titles (cards 3-7) ----
    t1::String; n1::Int
    t2::String; n2::Int
    xl::String; nx::Int
    yl::String; ny::Int
    rl::String; nr::Int

    # ---- Axis config (card 4) ----
    itype::Int
    jtype::Int
    igrid::Int
    ileg::Int
    xtag::Float64
    ytag::Float64

    # ---- Axis limits (cards 5-7) ----
    xmin::Float64; xmax::Float64; xstp::Float64
    ymin::Float64; ymax::Float64; ystp::Float64
    zmin::Float64; zmax::Float64; zstp::Float64

    # ---- Curve style (card 9) ----
    icon::Int
    isym::Int
    idash::Int
    iccol::Int
    ithick::Int
    ishade::Int

    # ---- Legend (card 10/10a) ----
    aleg::String; nleg::Int
    xpoint::Float64

    # ---- 3D viewpoint (card 11) ----
    x3::Float64; y3::Float64; z3::Float64
    xv::Float64; yv::Float64; zv::Float64

    # ---- Computed quantities ----
    hlab::Float64       # computed label character height
    hleg::Float64       # computed legend character height
    itic::Int           # tick direction (Fortran SAVE in set2d)
    xpt::Float64        # legend block x (absolute coords)
    ypt::Float64        # legend block y (absolute coords)
    xg::Float64         # axis width (set by init2)
    yg::Float64         # axis height (set by init2)
end

function ViewrState()
    ViewrState(
        0.0, 0.0, DLINE, DSIZE, 1, 2, 0,  # page setup
        0, 1.0, 1.0,                        # window
        "", 0, "", 0, "", 0, "", 0, "", 0,  # labels
        0, 0, 0, 0, 0.0, 0.0,              # axis config
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,     # axis limits x/y
        0.0, 0.0, 0.0,                      # axis limits z
        0, 0, 0, 0, 1, 0,                  # curve style
        "", 0, 0.0,                          # legend
        2.5, 6.5, 2.5, 15.0, -15.0, 15.0,  # 3D viewpoint
        0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0,  # computed
    )
end
