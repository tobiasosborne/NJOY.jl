# covr_plot.jl -- Emit viewr plot-tape commands for one covariance frame group.
#
# Ports covr.f90:1014-1647 (`plotit`, `matshd`, `level`, `patlev`).
#
# Each call to `write_covr_frames!` emits 3 (or more, when matshd splits the
# correlation matrix into many shaded regions) frames of viewr commands:
#   1. semi-log std-dev (or cross-section, when mtflg=1) of MT, horizontal
#   2. same for MT1, rotated 90° (matches the right-hand axis of frame 3)
#   3. shaded contour of the correlation matrix
#   + axis text overlays  (covr.f90:1264-1300)
#   + correlation legend (covr.f90:1524-1593) — drawn once per group.
#
# All Fortran format strings are reproduced exactly:
#   '5f8.3'    -> @sprintf("%8.3f%8.3f%8.3f%8.3f%8.3f", ...)
#   '1p,2e12.4' -> @sprintf("%12.4E%12.4E", ...)
#   '1p,2e13.4' -> @sprintf("%13.4E%13.4E", ...)
#   '1p,3e12.4' -> @sprintf("%12.4E%12.4E%12.4E", ...)
#   'f6.3'     -> @sprintf("%6.3f", ...)
# Trailing slash terminator is appended verbatim per Fortran "write … /".

using Printf

# ---------- pattern look-up (covr.f90:1621-1647) ----------

"""
    covr_patlev(ilevel::Integer, icolor::Integer, nlev::Integer, ndiv::Integer) -> Int

Convert a correlation level index to the viewr shading-pattern code
(`jpat` in Fortran). Mirrors covr.f90:1621-1647.

Conventions:
- ilevel=±1 → no shading (jpat=0; the lowest band is hidden, matching
  Fortran's `if ilevel.gt.1` / `if ilevel.lt.-1` guards).
- Monochrome (icolor=0): hatching codes 20-39, calculated as
  `20 - (nlev-ilevel)*2/ndiv` for positives (with integer division ÷).
- Color (icolor≥1): codes 50-59 with `scale = 10/(nlev-1)`.
"""
function covr_patlev(ilevel::Integer, icolor::Integer, nlev::Integer, ndiv::Integer)::Int
    jpat = 0
    if icolor == 0
        if ilevel > 1
            jpat = 20 - (nlev - ilevel) * 2 ÷ ndiv
        elseif ilevel < -1
            jpat = 40 - (nlev + ilevel) * 2 ÷ ndiv
        end
    else
        scale = 10.0 / (nlev - 1)
        if ilevel > 1
            t = Float64(nlev - ilevel)
            jpat = 50 - round(Int, scale * t)
        elseif ilevel < -1
            t = Float64(nlev + ilevel)
            jpat = 60 - round(Int, scale * t)
        end
    end
    jpat
end

# ---------- correlation-level lookup (covr.f90:1601-1619) ----------
"""
    covr_level(c::Real, xlev::AbstractVector{<:Real}, nlev::Integer) -> Int

Return the (signed) level index of correlation coefficient `c` on the
boundary scale `xlev[1..nlev]`. The first `i` such that `|c| < xlev[i]`
gives `ilev = i`; if `c < 0` the result is negated. If no boundary is
exceeded, `ilev = nlev` (note: Fortran sets `ilev = i` inside the loop and
keeps the loop's last index; Julia mirrors that).
"""
function covr_level(c::Real, xlev::AbstractVector{<:Real}, nlev::Integer)::Int
    ilev = 1
    for i in 1:nlev
        ilev = i
        abs(c) < xlev[i] && break
    end
    c < 0 && (ilev = -ilev)
    ilev
end

# ---------- expanded shade levels (covr.f90:286-306) ----------
"""
    covr_build_xlev(tlev::Vector{Float64}, ndiv::Integer) -> Vector{Float64}

Apply the `ndiv` subdivision of each tlev band. After this, the returned
`xlev` has `length(tlev)*ndiv` entries, each `tlow + (tlev[k]-tlow)*j/ndiv`.
Final entry is bumped by `da = 1/1000` (covr.f90:302-304) so that values
exactly equal to the top tlev still fall inside the highest band.
"""
function covr_build_xlev(tlev::Vector{Float64}, ndiv::Integer)::Vector{Float64}
    nlev = length(tlev)
    xlev = Float64[]
    tlow = 0.0
    for i in 1:nlev
        tx = (tlev[i] - tlow) / ndiv
        for j in 1:ndiv
            push!(xlev, tlow + tx * j)
        end
        tlow = tlev[i]
    end
    xlev[end] += 0.001     # da = 1.0/1000
    xlev
end

# ---------- axis log-round helpers (covr.f90:1106-1122 / 1206-1222) ----------

# Round xmin DOWN to a nice decade or 2× / 5× decade boundary.
function _axis_round_low(v::Float64)::Float64
    t = log10(v)
    t1 = t < 0.0 ? t - 0.99999 : t + 0.00001
    i = trunc(Int, t1)
    out = 10.0^i
    if t - i > log10(2.0); out = 2.0 * 10.0^i; end
    if t - i > log10(5.0); out = 5.0 * 10.0^i; end
    out
end

# Round xmax UP to a nice decade / 2× / 5×.
function _axis_round_high(v::Float64)::Float64
    t = log10(v)
    t = t - t / 200.0
    i = trunc(Int, t)
    out = 10.0^(i + 1)
    if t - i < log10(5.0); out = 5.0 * 10.0^i; end
    if t - i < log10(2.0); out = 2.0 * 10.0^i; end
    out
end

# Clamp std-dev percentage to plot bounds and flag warning. Mirrors
# covr.f90:1152-1175 for mtflg=0 (rsd in %) and 1167-1175 for mtflg=1
# (raw cross section / nubar).
function _clamp_yyy(y::Float64, mtflg::Int, yrange::Float64, yymin::Float64,
                    iwarn_in::Int)::Tuple{Float64,Int}
    yyym  = 59.999
    ylogmn = 0.1001
    ylogmx = 99.9
    yrtest = 10.0
    iwarn = iwarn_in
    if mtflg == 0
        yyy = 100 * y
        if yrange > yrtest && yyy < ylogmn
            yyy != 0.0 && (iwarn = 1)
            yyy = ylogmn
        end
        if yrange > yrtest && yyy > ylogmx
            iwarn = 1; yyy = ylogmx
        elseif yrange <= yrtest && yyy > yyym
            iwarn = 1; yyy = yyym
        end
        return (yyy, iwarn)
    else
        yyy = y
        if yrange > yrtest && yyy <= yymin
            yyy != 0.0 && (iwarn = 1)
            yyy = 1.0001 * yymin
        end
        return (yyy, iwarn)
    end
end

# ---------- plot-tape line helpers ----------

# Write Fortran "5f8.3" — five field-width-8 floats with 3 decimal places
# (covr.f90:1126 etc.).
_fmt_f5(a, b, c, d, e) = @sprintf("%8.3f%8.3f%8.3f%8.3f%8.3f", a, b, c, d, e)

# Write Fortran "1p,2e12.4" — two width-12 scientific values (covr.f90:1136).
_fmt_2e12(a, b) = @sprintf("%12.4E%12.4E", a, b)

# Write Fortran "1p,2e13.4" — two width-13 scientific values (covr.f90:1176).
_fmt_2e13(a, b) = @sprintf("%13.4E%13.4E", a, b)

# Write Fortran "1p,3e12.4" — three width-12 scientific (covr.f90:1541).
_fmt_3e12(a, b, c) = @sprintf("%12.4E%12.4E%12.4E", a, b, c)

# Single quote char (Fortran `qu = '''`).  We use ASCII apostrophe.
const _COVR_Q = "'"

# ---------- writing one frame group ----------

"""
    PlotState

Per-(case,sub-case) running state that the plot writer needs to keep
across calls. Mirrors variables `nfig` and `ipat` in covr.f90.
"""
mutable struct PlotState
    icolor::Int
    nlev::Int
    ndiv::Int
    xlev::Vector{Float64}
    nstart::Int
    noleg::Int
    nfig::Int
    nfigmx::Int
    iverf::Int
end

PlotState(p::CovrParams; iverf::Integer=6, nfigmx::Integer=100) = begin
    xlev = covr_build_xlev(p.tlev, p.ndiv)
    PlotState(p.icolor, length(xlev), p.ndiv, xlev,
              p.nstart, p.noleg, p.nstart - 1, nfigmx, iverf)
end

"""
    write_covr_header!(io, icolor)

Write the very first record of a covr plot tape (covr.f90:335):
    '1 2 .22'<i3>'/'   — viewr "begin plot, color mode" command.
"""
function write_covr_header!(io::IO, icolor::Integer)
    @printf(io, "1 2 .22%3d/\n", icolor)
end

"""
    write_covr_footer!(io)

Write the final record of a covr plot tape (covr.f90:482):
    '99/'   — viewr "end of document" command.
"""
write_covr_footer!(io::IO) = println(io, "99/")

"""
    write_covr_frames!(io, state, cr, mat, mt, mat1, mt1, mfflg, mf35, izap, ne)
        -> ishade::Int

Emit the three viewr-format frames for one (mat,mt,mat1,mt1) reaction
pair (after `corr` has produced the CorrResult `cr`). Returns the number
of distinct shading regions written (`ishade`). When `ishade == 0` the
caller treats the matrix as plot-empty (covr.f90:411-417): mark izero=2
and skip the matrix legend. Even with `ishade == 0` the two semi-log
frames are *already* on the tape — Fortran does likewise (the matrix
frame is the only one suppressed).
"""
function write_covr_frames!(io::IO, state::PlotState, cr::CorrResult,
                            mat::Integer, mt::Integer,
                            mat1::Integer, mt1::Integer,
                            mfflg::Integer, mf35::Integer,
                            izap::Integer, ne::Integer,
                            iza::Integer, iza1::Integer,
                            ixmin::Integer)::Int
    # Geometry constants (covr.f90:1030-1042).
    dx1, dy1 = 1.0, 0.75
    dx2, dy2 = -0.25, 0.5
    dx3, dy3 = -0.75, 1.9
    dy4      = 0.75
    xsize    = 5.00
    ysize    = 3.38

    ixmax = cr.ixmax
    ixn   = ixmax - ixmin + 1

    # ---- truncated views (covr.f90:414 — `call plotit(x(ixmin:),...)`) ----
    # plotit is called with sliced x/y/rsdx/rsdy so its local index 1
    # corresponds to global ixmin.  We mirror with truncated copies so the
    # local i-th element holds the global (ixmin + i - 1)-th element.
    xs_x = view(cr.x, ixmin:(ixmax + 1))
    ys_y = view(cr.y, ixmin:(ixmax + 1))
    xx = cr.xx
    xy = cr.xy
    rsdx = collect(view(cr.rsdx, ixmin:ixmax))   # length ixn
    rsdy = collect(view(cr.rsdy, ixmin:ixmax))

    # mtflg=1 self-cov: overwrite rsdx with the per-group cross section
    # (or, for mf35=5, the cross section / dE) — covr.f90:1048-1059.
    mtflg = (mat == mat1 && mt == mt1) ? 1 : 0
    if mtflg == 1
        if mf35 != 5
            for i in 1:ixn
                rsdx[i] = xx[i + (ixmax - ixn)]
            end
        else
            for i in 1:ixn
                de = ys_y[i + 1] - ys_y[i]
                rsdx[i] = xx[i + (ixmax - ixn)] / de
            end
        end
    end

    # ----- index of first non-zero rsdx, rsdy, then ii := max(ii, jj) -----
    ii = 1
    while ii <= ixn && rsdx[ii] <= 0.0; ii += 1; end
    jj = 1
    while jj <= ixn && rsdy[jj] <= 0.0; jj += 1; end
    jj > ii && (ii = jj)

    # =====================================================================
    # FRAME 1 — MT std-dev / cross section, horizontal axis = energy
    # (covr.f90:1080-1141)
    # =====================================================================
    xig, yig, npts, ymin1, ymax1, yrange1 =
        _build_step_curve(xs_x, rsdx, ii, ixn)
    xmin = _axis_round_low(xig[1])
    xmax = _axis_round_high(xig[npts])
    ydec = log10(ymin1)
    ydec < 0.0 && (ydec -= 1.0)
    yymin1 = 10.0^trunc(Int, ydec)

    xpos1 = ysize - dy1
    ypos1 = xsize - dx1
    @printf(io, "1 0 1. 1.%s/\n", _fmt_f5(xpos1, ypos1, xsize, ysize, 0.0))
    title1 = covr_smilab(iza, mat, mt, izap, mtflg, mfflg, mf35, cr.einc;
                        iverf=state.iverf)
    # Fortran plotit declares `character(80)::strng` (covr.f90:1025); the
    # write at covr.f90:1129 emits the full 80 chars (incl. trailing spaces).
    println(io, _COVR_Q, rpad(title1, 80), _COVR_Q, "/")
    println(io, "/")
    println(io, yrange1 > 10.0 ? "4/" : "3/")
    println(io, _fmt_2e12(xmin, xmax), "/")
    println(io, _COVR_Q, ".", _COVR_Q, "/")
    println(io, "/")
    println(io, _COVR_Q, ".", _COVR_Q, "/")
    println(io, "/")
    println(io, "0 0 0 3 2/")
    println(io, "0")
    iwarn = 0
    for i in 1:npts
        yyy, iwarn = _clamp_yyy(yig[i], mtflg, yrange1, yymin1, iwarn)
        println(io, _fmt_2e13(xig[i], yyy), "/")
    end
    println(io, "/")

    # =====================================================================
    # FRAME 2 — MT1 std-dev, vertical axis (rotated 90°)
    # (covr.f90:1180-1261)
    # =====================================================================
    xig2, yig2, npts2, ymin2_lo, ymax2, yrange2 =
        _build_step_curve(ys_y, rsdy, ii, ixn)
    # In Fortran:  ymin = xig(1) (used as plot lower bound), and if mtflg=1
    # ymin = xmin (the FRAME-1 xmin, so the matrix renders square).
    ymin_axis = xig2[1]
    if mtflg == 1
        ymin_axis = xig[1]   # frame-1 xmin1 (pre-axis-rounding)
    end
    ymin_round = _axis_round_low_e(ymin_axis)
    ymax_round = _axis_round_high(xig2[npts2])

    xpos2 = ysize + dy2
    ypos2 = dx2
    @printf(io, "-1 0 1. 1.%s/\n", _fmt_f5(xpos2, ypos2, xsize, ysize, 90.0))
    title2 = covr_smilab(iza1, mat1, mt1, izap, 0, mfflg, mf35, cr.einc;
                        iverf=state.iverf)
    # Same character(80) padding as frame 1 (covr.f90:1025/1229).
    println(io, _COVR_Q, rpad(title2, 80), _COVR_Q, "/")
    println(io, "/")
    println(io, yrange2 > 10.0 ? "4/" : "3/")
    println(io, _fmt_2e12(ymin_round, ymax_round), "/")
    println(io, _COVR_Q, ".", _COVR_Q, "/")
    println(io, "/")
    println(io, _COVR_Q, ".", _COVR_Q, "/")
    println(io, "/")
    println(io, "0 0 0 3 2/")
    println(io, "0")
    for i in 1:npts2
        yyy, iwarn = _clamp_yyy(yig2[i], 0, yrange2, 0.0, iwarn)
        println(io, _fmt_2e13(xig2[i], yyy), "/")
    end
    println(io, "/")

    # =====================================================================
    # AXIS-DESCRIPTION TEXT BLOCKS (covr.f90:1263-1300)
    # =====================================================================
    xpos3 = dy3
    ypos3 = xsize + dx3
    @printf(io, "-1 0 1. 1.%s/\n", _fmt_f5(xpos3, ypos3, ysize, 1.0, 90.0))
    if mtflg == 0
        println(io, _COVR_Q, "#H.75<o>rdinate scale is %", _COVR_Q, "/")
        println(io, _COVR_Q, "#H.75<>relative standard deviation.", _COVR_Q, "/")
    elseif mt == 251
        println(io, _COVR_Q, "#H.75<o>rdinate scales are % relative", _COVR_Q, "/")
        println(io, _COVR_Q, "#H.75<>standard deviation and mu-bar.", _COVR_Q, "/")
    elseif 452 <= mt <= 456
        println(io, _COVR_Q, "#H.75<o>rdinate scales are % relative", _COVR_Q, "/")
        println(io, _COVR_Q, "#H.75<>standard deviation and nu-bar.", _COVR_Q, "/")
    elseif mf35 == 5
        println(io, _COVR_Q, "#H.75<o>rdinate scales are % standard", _COVR_Q, "/")
        println(io, _COVR_Q, "#H.75<>deviation and spectrum/eV.", _COVR_Q, "/")
    else
        println(io, _COVR_Q, "#H.75<o>rdinate scales are % relative", _COVR_Q, "/")
        println(io, _COVR_Q, "#H.75<>standard deviation and barns.", _COVR_Q, "/")
    end
    println(io, "0/")
    xpos3 += dy4
    @printf(io, "-1 0 1. 1.%s/\n", _fmt_f5(xpos3, ypos3, ysize, 1.0, 90.0))
    println(io, _COVR_Q, "#H.75<a>bscissa scales are energy (e<v>).", _COVR_Q, "/")
    println(io, "/")

    if iwarn != 0
        xpos3 += 0.45
        println(io, "0/")
        @printf(io, "-1 0 1. 1.%s/\n", _fmt_f5(xpos3, ypos3, ysize, 1.0, 90.0))
        println(io, _COVR_Q, "#H.75<w>arning:  some uncertainty", _COVR_Q, "/")
        println(io, _COVR_Q, "#H.75<>data were suppressed.", _COVR_Q, "/")
    end
    println(io, "0/")

    # =====================================================================
    # FRAME 3 — correlation matrix (shaded contour) + legend
    # (covr.f90:1302 → matshd)
    # =====================================================================
    ishade = _matshd!(io, state, cr, ixmin, ixn, mat, mt, mat1, mt1,
                     xmin, xmax, ymin_round, ymax_round)
    ishade
end

# Build a step-curve representation of (energy_bound, value) used for the
# semi-log frames. Mirrors the inner loops at covr.f90:1078-1102.
function _build_step_curve(grid::AbstractVector{Float64}, vals::AbstractVector{Float64},
                            ii::Int, ixn::Int)
    j = 0
    xig = Float64[grid[ii]]
    yig = Float64[0.0]
    ymin =  1.0e20
    ymax = -1.0e20
    for i in ii:ixn
        j += 1
        push!(xig, grid[i]);     push!(yig, vals[i])
        if vals[i] > 0.0 && vals[i] < ymin; ymin = vals[i]; end
        vals[i] > ymax && (ymax = vals[i])
        push!(xig, grid[i + 1]); push!(yig, vals[i])
    end
    npts = 2 * j + 2
    push!(xig, grid[ixn + 1]); push!(yig, 0.0)
    yrange = ymin != 0.0 ? ymax / ymin : 0.0
    (xig, yig, npts, ymin, ymax, yrange)
end

# Axis-low rounding for the ROTATED frame: same logic as
# _axis_round_low but emits the t-i>0 short-circuit (covr.f90:1213).
function _axis_round_low_e(v::Float64)::Float64
    t = log10(v)
    t1 = t < 0.0 ? t - 0.99999 : t + 0.00001
    i = trunc(Int, t1)
    out = v
    if t - i > 0.0; out = 10.0^i; end
    if t - i > log10(2.0); out = 2.0 * 10.0^i; end
    if t - i > log10(5.0); out = 5.0 * 10.0^i; end
    out
end

# ---------- matshd ----------------------------------------------------------
# Connected-region shading of the correlation matrix.  Direct port of
# covr.f90:1310-1599.

function _matshd!(io::IO, state::PlotState, cr::CorrResult,
                  ixmin::Int, ixnow::Int,
                  mat::Integer, mt::Integer,
                  mat1::Integer, mt1::Integer,
                  xmin::Float64, xmax::Float64,
                  ymin::Float64, ymax::Float64)::Int
    ixmax = cr.ixmax
    nlev  = state.nlev
    xlev  = state.xlev
    icolor = state.icolor
    ndiv   = state.ndiv

    xsize = 5.0
    ysize = 3.38
    dx1m, dy1m = -0.25, -0.75
    dx2m, dy2m =  0.75,  0.60
    dx3m, dy3m =  0.625, 1.75
    dy4m       = 1.375
    eps = 1.0e-6

    # ---- compute level grid (ilcf), clip to [-1,1], track extremes ----
    ilcf = zeros(Int, ixnow * ixnow)
    none_count = 0
    ntwo_count = 0
    cofm = 0.0
    ixmx = 0; jxmx = 0
    ii = 0
    @inbounds for i in ixmin:ixmax
        ii += 1
        jj = 0
        for j in ixmin:ixmax
            jj += 1
            cof = cr.cf[ixmax * (i - 1) + j]
            cofa = abs(cof) - eps
            if cofa > abs(cofm)
                cofm = cof
                ixmx = i
                jxmx = j
            end
            if cofa > 2.0
                ntwo_count += 1
            elseif cofa > 1.0
                none_count += 1
            end
            cof = cof >  1.0 ?  1.0 : cof
            cof = cof < -1.0 ? -1.0 : cof
            ilcf[ixnow * (ii - 1) + jj] = covr_level(cof, xlev, nlev)
        end
    end
    if none_count > 0
        @info @sprintf("matshd: mat/mt %d/%d vs mat1/mt1 %d/%d  largest=%.5e at idx %d %d  (%d > 1, reset)",
                       mat, mt, mat1, mt1, cofm, ixmx, jxmx, none_count)
    end
    ntwo_count > 0 && @info @sprintf("matshd: %d coefficients > 2 (reset and continue)", ntwo_count)

    # ---- pattern-search loop (covr.f90:1392-1522) ----
    ipat = 0
    inext = 0
    jnext = 1
    ishade = 0
    ixmip = zeros(Int, ixnow)
    ixmap = zeros(Int, ixnow)

    while true
        ipat += 1
        ipat > 99999 && error("_matshd!: ipat > 99999")

        # ---- find next un-allocated square (covr.f90:1396-1402) ----
        found = false
        while true
            inext += 1
            if inext > ixnow
                jnext += 1
                jnext > ixnow && @goto draw_legend
                inext = 1
            end
            if ilcf[ixnow * (inext - 1) + jnext] != 0
                found = true; break
            end
        end
        @assert found
        il = inext
        jl = jnext
        ilevel = ilcf[ixnow * (il - 1) + jl]
        ilcf[ixnow * (il - 1) + jl] = 0

        # ---- right-extent of pattern in starter row (covr.f90:1408-1420) ----
        ixmip[jl] = il
        ixmap[jl] = il
        jfirst = jl
        jlast  = jl
        while true
            il += 1
            il > ixnow && break
            ilcf[ixnow * (il - 1) + jl] != ilevel && break
            ilcf[ixnow * (il - 1) + jl] = 0
            ixmap[jl] = il
            inext = il
        end

        # ---- subsequent rows (covr.f90:1422-1461) ----
        while true
            jl += 1
            jl > ixnow && break
            i1 = ixmip[jl - 1]
            i2 = ixmap[jl - 1]
            istart = 0
            for i in i1:i2
                if ilcf[ixnow * (i - 1) + jl] == ilevel
                    istart = i; break
                end
            end
            istart == 0 && break  # end of pattern
            jlast = jl
            ixmip[jl] = istart
            ixmap[jl] = istart
            ilcf[ixnow * (istart - 1) + jl] = 0
            # extend left
            if istart != 1
                ilim = istart - 1
                for i in 1:ilim
                    ilcf[ixnow * (istart - i - 1) + jl] != ilevel && break
                    ixmip[jl] = istart - i
                    ilcf[ixnow * (istart - i - 1) + jl] = 0
                end
            end
            # extend right
            if istart != ixnow
                ilim = ixnow - istart
                for i in 1:ilim
                    ilcf[ixnow * (istart + i - 1) + jl] != ilevel && break
                    ixmap[jl] = istart + i
                    ilcf[ixnow * (istart + i - 1) + jl] = 0
                end
            end
        end

        # ---- build polygon contour (covr.f90:1463-1493) ----
        xig = Float64[]
        yig = Float64[]
        # forward sweep (right-edge of pattern, top → bottom)
        jl = jfirst
        while true
            ilast = ixmap[jl]
            push!(xig, cr.x[ilast + ixmin])    # x[ilast+1] in Fortran (1-based)
            push!(yig, cr.y[jl + ixmin - 1])   # y[jl] in Fortran
            push!(xig, cr.x[ilast + ixmin])
            push!(yig, cr.y[jl + ixmin])       # y[jl+1] in Fortran
            jl >= jlast && break
            jl += 1
        end
        # backward sweep (left-edge, bottom → top)
        jl = jlast + 1
        while true
            jl -= 1
            ifirst = ixmip[jl]
            push!(xig, cr.x[ifirst + ixmin - 1])  # x[ifirst] in Fortran
            push!(yig, cr.y[jl + ixmin])          # y[jl+1] in Fortran
            push!(xig, cr.x[ifirst + ixmin - 1])
            push!(yig, cr.y[jl + ixmin - 1])      # y[jl] in Fortran
            jl <= jfirst && break
        end

        # ---- emit polygon (covr.f90:1497-1521) ----
        if ipat == 1
            xpos = ysize + dy1m
            ypos = dx1m
            @printf(io, "-1 0 1. 1.%s/\n", _fmt_f5(xpos, ypos, xsize, xsize, 0.0))
            println(io, "/")
            println(io, "/")
            println(io, "4 0 0/")
            println(io, _fmt_2e12(xmin, xmax), "/")
            println(io, "/")
            println(io, _fmt_2e12(ymin, ymax), "/")
            println(io, "/")
        else
            println(io, "2/")
        end
        println(io, "/")
        jpat = covr_patlev(ilevel, icolor, nlev, ndiv)
        jpat != 0 && (ishade += 1)
        @printf(io, "0 0 0 0 0%3d/\n", jpat)
        println(io, "0")
        for k in 1:length(xig)
            println(io, _fmt_2e13(xig[k], yig[k]), "/")
        end
        println(io, "/")
    end
    @label draw_legend

    # ---- legend (covr.f90:1525-1593) ----
    _matshd_legend!(io, icolor, nlev, ndiv, xlev, xsize, ysize, dx2m, dy2m, dx3m, dy3m, dy4m)
    ishade
end

function _matshd_legend!(io::IO, icolor::Int, nlev::Int, ndiv::Int,
                          xlev::Vector{Float64},
                          xsize::Float64, ysize::Float64,
                          dx2::Float64, dy2::Float64, dx3::Float64, dy3::Float64,
                          dy4::Float64)
    xpos = ysize - dx2 + xsize + dy2
    ypos = dx2
    @printf(io, "-1 0 1. 1.%s/\n", _fmt_f5(xpos, ypos, xpos, 1.0, 90.0))
    println(io, _COVR_Q, "<c>orrelation <m>atrix", _COVR_Q, "/")
    println(io, "/")
    println(io, "0/")

    xpos += dy3
    ypos = dx3
    @printf(io, "-1 0 1. 1.%s/\n", _fmt_f5(xpos, ypos, 1.75, 2.50, 90.0))
    println(io, "/")
    println(io, "/")
    println(io, "1 0 0/")
    println(io, _fmt_3e12(0.0, 1.0, 1.0), "/")
    println(io, "/")
    println(io, _fmt_3e12(0.0, 1.0, 0.2), "/")
    println(io, _COVR_Q, ".", _COVR_Q, "/")
    println(io, "/")

    # positive part (covr.f90:1546-1564)
    tlow = 0.0
    for ilevel in 1:nlev
        jpat = covr_patlev(ilevel, icolor, nlev, ndiv)
        if icolor == 0 && jpat != 0
            jpat += 10
        end
        @printf(io, "0 0 0 0 0%3d/\n", jpat)
        println(io, "0")
        thi = xlev[ilevel]
        @printf(io, "0. %6.3f/\n", tlow)
        @printf(io, "1. %6.3f/\n", tlow)
        @printf(io, "1. %6.3f/\n", thi)
        @printf(io, "0. %6.3f/\n", thi)
        @printf(io, "0. %6.3f/\n", tlow)
        println(io, "/")
        if ilevel < nlev
            println(io, "2/")
            println(io, "/")
        end
        tlow = thi
    end

    ypos += dy4
    @printf(io, "-1 0 1. 1.%s/\n", _fmt_f5(xpos, ypos, 1.75, 2.50, 90.0))
    println(io, "/")
    println(io, "/")
    println(io, "1 0 0/")
    println(io, _fmt_3e12(0.0, 1.0, 1.0), "/")
    println(io, "/")
    println(io, _fmt_3e12(0.0, -1.0, -0.2), "/")
    println(io, _COVR_Q, ".", _COVR_Q, "/")
    println(io, "/")

    # negative part (covr.f90:1576-1592)
    tlow = 0.0
    for iloop in 1:nlev
        ilevel = -iloop
        jpat = covr_patlev(ilevel, icolor, nlev, ndiv)
        @printf(io, "0 0 0 0 0%3d/\n", jpat)
        println(io, "0")
        thi = -xlev[iloop]
        @printf(io, "0. %6.3f/\n", tlow)
        @printf(io, "1. %6.3f/\n", tlow)
        @printf(io, "1. %6.3f/\n", thi)
        @printf(io, "0. %6.3f/\n", thi)
        @printf(io, "0. %6.3f/\n", tlow)
        println(io, "/")
        if iloop < nlev
            println(io, "2/")
            println(io, "/")
        end
        tlow = thi
    end
    nothing
end
