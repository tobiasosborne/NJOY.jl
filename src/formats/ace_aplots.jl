# Port of Fortran `aplots` (acefc.f90:15212-16830) and its delegates
# `aplotr` / `aplopf` / `aplof4` (16832-17428): emit a viewr-format plot tape
# describing an ACE file. Invoked from acer iopt=7 (acefix, acefc.f90:14205).
#
# The plot tape is a sequence of viewr "cards". aplots walks the ACE arrays and
# writes, in order: a global header, then a page per plottable quantity
# (principal cross sections log-log and lin-lin, optional resonance/URR/heating/
# damage/threshold pages, then angular-distribution 3D pages), terminated by
# `99/` and `stop`.
#
# Float formatting is matched byte-for-byte to Fortran:
#   1p,2e14.6 / 1p,e14.6  -> %14.6E             (_e14_6)
#   1p,3e12.3             -> %12.3E per field   (_e12_3)
#   2e12.4   (NO 1p)      -> 0.ddddE±dd form    (_e12_4_no1p)
#
# Ref: njoy-reference/src/acefc.f90:15212-17428, 19541-19724.

using Printf

"""
    AplotsNotPortedError(msg)

Raised by `_acer_aplots` when it reaches an aplots block that is not yet ported
(heating/damage curves, threshold/inelastic-level pages, particle-production
pages, etc. — the data exists in the ACE arrays, but the corresponding Fortran
sub-plot has not been transcribed yet). This is distinct from a generic
`error`: the caller (`_acer_iopt7`) catches *this* type specifically, falls back
to the documented empty plot-tape stub, and keeps writing the remaining output
tapes — so a deferred plot block degrades the plot tape to a stub rather than
crashing the whole acer module (which would lose tape34/tape35). Any other
exception still propagates and aborts loudly (Rule 6).
"""
struct AplotsNotPortedError <: Exception
    msg::String
end
Base.showerror(io::IO, e::AplotsNotPortedError) = print(io, "aplots: ", e.msg)
_aplots_unported(msg) = throw(AplotsNotPortedError(msg))

# =========================================================================
# Fortran float-format helpers
# =========================================================================

"""
`1p,e14.6` data field. For 2-digit exponents C `%14.6E` reproduces Fortran
exactly (e.g. `  2.950000E-01`, `  8.927518E+08`). For 3-digit exponents Fortran
runs out of field width and *drops the `E`* — `1p,e14.6` of `1.155100e-117`
prints `  1.155100-117` (mantissa, sign of exponent, three digits, no `E`), still
right-justified in 14 columns. C `%14.6E` keeps the `E` (` 1.155100E-117`, only
one leading space), so we strip the `E` and re-pad to width 14.

Ref: gfortran `1p,e14.6` edit descriptor (Ew.d with the exponent overflowing the
default 2-digit `Eee` field: the `E` is omitted, leaving `±eee`).
"""
function _e14_6(v::Real)
    s = @sprintf("%14.6E", Float64(v))
    # 2-digit exponent fits with the `E`; only 3+ digit exponents drop it.
    ei = findlast('E', s)
    ei === nothing && return s
    # exponent digits after the sign character (s[ei+1] is '+' or '-')
    ndig = length(s) - (ei + 1)
    ndig <= 2 && return s
    return lpad(string(s[1:ei-1], s[ei+1:end]), 14)
end

"""`1p,e12.3` axis field: C `%12.3E` reproduces it exactly (e.g. `   2.000E-01`)."""
_e12_3(v::Real) = @sprintf("%12.3E", Float64(v))

"""
    _e12_4_no1p(v::Real) -> String

Fortran `e12.4` WITHOUT the `1p` scale factor: the mantissa is normalised to
`0.dddd` (i.e. in `[0.1, 1.0)`), 4 significant digits, then `E`, sign, 2-digit
exponent — total width 12, right-justified. This is the legend-tag (xtag,ytag)
format on the scale card. C `%12.4E` gives the `1p` form ` 9.9990E-01`, which is
wrong; we build the no-`1p` form `  0.9999E+00` directly.

Ref: acefc.f90:15312 `write(nout,'(''4 0 2 1'',2e12.4,''/'')') xtag,ytag`.
"""
function _e12_4_no1p(v::Real)
    x = Float64(v)
    if x == 0.0
        return lpad("0.0000E+00", 12)
    end
    neg = x < 0
    a = abs(x)
    # Normalise to mantissa in [0.1, 1.0): a = m * 10^e, 0.1 <= m < 1.
    e = floor(Int, log10(a)) + 1
    m = a / 10.0^e
    # Guard against FP edge cases pushing m out of [0.1,1.0).
    if m >= 1.0
        m /= 10.0; e += 1
    elseif m < 0.1
        m *= 10.0; e -= 1
    end
    # Round mantissa to 4 digits; a carry can bump it to 1.0000 -> renormalise.
    digits4 = round(Int, m * 10000)
    if digits4 >= 10000
        digits4 = div(digits4, 10); e += 1
    end
    mant = @sprintf("%04d", digits4)
    es = e < 0 ? "-" : "+"
    body = string(neg ? "-" : "", "0.", mant, "E", es, @sprintf("%02d", abs(e)))
    return lpad(body, 12)
end

# =========================================================================
# Axis scaling: ascll (log) and ascle (linear)
# =========================================================================

"""
    _ascll(amin, amax) -> (amin', amax')

Adjust axis limits for a log scale. Snaps log10(limits) to decade boundaries at
0, .301, .699, 1. Ref: acefc.f90:19684-19724.
"""
function _ascll(amin::Float64, amax::Float64)
    b1 = 0.0; b2 = 0.301; b3 = 0.699; b4 = 1.0
    ohoh1 = 0.001; ten = 10.0

    top = log10(amax)
    itop = trunc(Int, top)
    top < 0 && (itop -= 1)
    amaxl = itop + 1.0
    top <= itop + b3 + ohoh1 && (amaxl = itop + b3)
    top <= itop + b2 + ohoh1 && (amaxl = itop + b2)
    top <= itop + b1 + ohoh1 && (amaxl = itop + b1)
    amax_out = ten^amaxl

    bot = log10(amin)
    ibot = trunc(Int, bot)
    bot < 0 && (ibot -= 1)
    aminl = Float64(ibot)
    bot >= ibot + b2 - ohoh1 && (aminl = ibot + b2)
    bot >= ibot + b3 - ohoh1 && (aminl = ibot + b3)
    bot >= ibot + b4 - ohoh1 && (aminl = ibot + b4)
    amin_out = ten^aminl

    if amin_out == amax_out
        amin_out = 0.1 * amin_out
        amax_out = 10.0 * amax_out
    end
    return amin_out, amax_out
end

"""
    _ascle(m, z1, z2) -> (z1', z2', major, minor)

Automatic linear-axis scaling (Los Alamos SC4020). Returns adjusted limits and
the major/minor division counts. Ref: acefc.f90:19541-19682.
"""
function _ascle(m::Int, z1::Float64, z2::Float64)
    one = 1.0; ten = 10.0; tentho = 10000.0
    zmin = z1; zmax = z2

    if zmax <= zmin || m <= 0 || m > 20
        return 0.0, 2*z2, 1, 0
    end

    fm = Float64(m)
    if zmax == 0.0 || zmin == 0.0
        zmax = zmax - abs(zmax)/1_000_000
        zmin = zmin + abs(zmin)/1_000_000
    else
        zbar = zmax/zmin
        if abs(zbar)/1000 >= one
            zmin = 0.0
        end
        if abs(zbar)*1000 <= one
            zmax = 0.0
            zmax = zmax - abs(zmax)/1_000_000
            zmin = zmin + abs(zmin)/1_000_000
        else
            if (abs(zbar - 1) - 5*fm/100000) > 0.0
                zmax = zmax - abs(zmax)/1_000_000
                zmin = zmin + abs(zmin)/1_000_000
            else
                zbar = (zmax + zmin)/2
                z = 26*fm*abs(zbar)/1_000_000
                zmax = zbar + z
                zmin = zbar - z
            end
        end
    end
    p = (zmax - zmin)/fm

    iflag = 0; tenk = 1.0; k = 0
    if p < one
        iflag = 1
        p = 1/p
    end
    while p >= tentho
        p /= 10000; tenk *= 10000; k += 4
    end
    while p >= ten
        p /= 10; tenk *= 10; k += 1
    end
    if iflag != 0
        p = 10/p
        tenk = one/10/tenk
        k -= 1
    end

    nm = 5
    test = 2 + one/100
    if p <= test
        test = 2 - one/100
        if p <= test
            p = 1.0; nm = 5
        else
            p = 2.0; nm = 4
        end
    else
        test = 5 - one/100
        if p < test
            p = 2.0; nm = 4
        else
            p = 5.0; nm = 5
        end
    end
    dz = p*tenk

    n1 = trunc(Int, zmin/dz)
    z = n1*dz
    if z > zmin
        z -= dz; n1 -= 1
    end
    zmin = z
    z1o = zmin
    n2 = trunc(Int, zmax/dz)
    z = n2*dz
    if z < zmax
        n2 += 1; z += dz
    end
    zmax = z
    z2o = zmax

    major = n2 - n1
    minor = nm*major
    return z1o, z2o, major, minor
end

# =========================================================================
# izai (incident-particle id) from the ZAID class letter (mcnpx=0).
# Ref: acefc.f90:14042-14075.
# =========================================================================
function _aplots_izai(ht::Char)
    ht == 'p' && return 0
    ht == 'u' && return 0
    ht == 'c' && return 1
    ht == 't' && return 1
    ht == 'y' && return 1
    ht == 'h' && return 1001
    ht == 'o' && return 1002
    ht == 'r' && return 1003
    ht == 's' && return 2003
    ht == 'a' && return 2004
    error("aplots: problem with particle id in zaid (class letter '$ht')")
end

# =========================================================================
# Main entry
# =========================================================================

"""
    _acer_aplots(io::IO, table::ACETable, hk::AbstractString, ht::Char)

Write the viewr plot tape for an ACE table. `hk` is the comment string (used as
the page titles, wrapped in `<...>`); `ht` is the ZAID class letter (column 10),
which selects `izai`.

Ref: njoy-reference/src/acefc.f90:15212-16830 (aplots).
"""
function _acer_aplots(io::IO, table::ACETable, hk::AbstractString, ht::Char)
    xss = table.xss
    nxs = table.nxs
    jxs = table.jxs
    nes = Int(nxs[NXS_NES])
    ntr = Int(nxs[NXS_NTR])
    nr  = Int(nxs[NXS_NR])
    esz  = Int(jxs[JXS_ESZ])
    gpd  = Int(jxs[JXS_GPD])
    iurpt = Int(jxs[JXS_IURPT])
    mtr  = Int(jxs[JXS_MTR])
    lsig = Int(jxs[JXS_LSIG])
    sig  = Int(jxs[JXS_SIG])
    land = Int(jxs[JXS_LAND])
    andb = Int(jxs[JXS_AND])
    izai = _aplots_izai(ht)

    big = 1.0e10; scale = 1.0e6; hmin = 1.0e-10
    ten = 10.0; small = 1.0e-12
    nden = 4000
    iwcol = 3            # colored pages (acefc.f90:15252-15253)
    ipcol = 2

    # hk truncated at the last non-blank (it), as Fortran's `hk(1:it)`.
    _it(s) = (j = 1; for i in 1:min(70, length(s)); s[i] != ' ' && (j = i); end; j)
    hk70 = rpad(String(hk), 70)
    it = _it(hk70)
    title = hk70[1:it]
    qu = '\''

    pr(s) = print(io, s, '\n')

    # ----------------------------------------------------------------
    # Global header card.  Ref: acefc.f90:15273.
    # ----------------------------------------------------------------
    pr(@sprintf("1 2 .30%3d/", ipcol))

    # esz sub-block accessors (1-based, matching xss(esz-1+i) etc.)
    e_i(i)    = xss[esz - 1 + i]
    tot_i(i)  = xss[esz + nes - 1 + i]
    abso_i(i) = xss[esz + 2*nes - 1 + i]
    elas_i(i) = xss[esz + 3*nes - 1 + i]
    heat_i(i) = xss[esz + 4*nes - 1 + i]
    gprod_i(i) = xss[gpd - 1 + i]

    # ----------------------------------------------------------------
    # PAGE: log-log principal cross sections.  Ref: acefc.f90:15275-15402.
    # ----------------------------------------------------------------
    xmin = big; xmax = 0.0; ymin = big; ymax = -big
    for i in 1:nes
        e = e_i(i); tot = tot_i(i); abso = abso_i(i); elas = elas_i(i)
        e < xmin && (xmin = e);   e > xmax && (xmax = e)
        tot < ymin && (ymin = tot); tot > ymax && (ymax = tot)
        abso < ymin && (ymin = abso); abso > ymax && (ymax = abso)
        elas < ymin && (ymin = elas); elas > ymax && (ymax = elas)
        if gpd != 0
            gp = gprod_i(i)
            gp < ymin && (ymin = gp); gp > ymax && (ymax = gp)
        end
    end
    xmin, xmax = _ascll(xmin, xmax)
    ymin < ymax/scale && (ymin = ymax/scale)
    ymin, ymax = _ascll(ymin, ymax)
    pr(@sprintf("1%3d/", iwcol))
    pr("$(qu)<$(title)>$(qu)/")
    pr("$(qu)<p>rincipal cross sections$(qu)/")
    xtag = 5*xmin
    ytag = 7*log10(ymin)/10 + 3*log10(ymax)/10
    ytag = ten^ytag
    pr("4 0 2 1$(_e12_4_no1p(xtag))$(_e12_4_no1p(ytag))/")
    pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(1.0))/")
    pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
    pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(1.0))/")
    pr("$(qu)<c>ross section (barns)$(qu)/")
    pr("/"); pr("/")
    thin = ten^(log10(xmax/xmin)/nden)

    # total
    pr("$(qu)total$(qu)/"); pr("0/")
    xlast = small
    for i in 1:nes
        e = e_i(i); tot = tot_i(i)
        tot < ymin && (tot = ymin)
        if nes <= nden || e >= thin*xlast
            pr("$(_e14_6(e))$(_e14_6(tot))/")
            xlast = e
        end
    end
    pr("/")
    # absorption
    pr("2/"); pr("/"); pr("0 0 0 3/")
    pr("$(qu)absorption$(qu)/"); pr("0/")
    xlast = small
    for i in 1:nes
        e = e_i(i); abss = abso_i(i)
        abss < ymin && (abss = ymin)
        if nes <= nden || e >= thin*xlast
            pr("$(_e14_6(e))$(_e14_6(abss))/")
            xlast = e
        end
    end
    pr("/")
    # elastic
    pr("3/"); pr("/"); pr("0 0 0 2/")
    pr("$(qu)elastic$(qu)/"); pr("0/")
    xlast = small
    for i in 1:nes
        e = e_i(i); elas = elas_i(i)
        elas < ymin && (elas = ymin)
        if nes <= nden || e >= thin*xlast
            pr("$(_e14_6(e))$(_e14_6(elas))/")
            xlast = e
        end
    end
    pr("/")
    # gamma production (only if gpd != 0)
    if gpd != 0
        pr("4/"); pr("/"); pr("0 0 0 1/")
        pr("$(qu)gamma production$(qu)/"); pr("0/")
        xlast = small
        for i in 1:nes
            e = e_i(i); gp = gprod_i(i)
            gp < ymin && (gp = ymin)
            if nes <= nden || e >= thin*xlast
                pr("$(_e14_6(e))$(_e14_6(gp))/")
                xlast = e
            end
        end
        pr("/")
    end

    # ----------------------------------------------------------------
    # Resonance / URR / log heating / log damage / log non-threshold pages.
    # All guarded; for non-resonance light targets (nes<1500, iurpt=0, heating
    # all-zero, no MT=444, ntr=0) they emit nothing.  Ref: 15404-16288.
    # ----------------------------------------------------------------
    if nes >= 1500
        _aplots_unported("nes=$nes >= 1500 resonance-expansion plotting not yet ported")
    end
    if iurpt != 0
        _aplots_unported("iurpt=$iurpt URR plotting not yet ported")
    end

    # log-log heating per reaction.  Ref: 16177-16226.
    xmin = big; xmax = 0.0; ymin = big; ymax = -big
    for i in 1:nes
        e = e_i(i); heat = heat_i(i)
        heat < hmin && (heat = hmin)
        e < xmin && (xmin = e); e > xmax && (xmax = e)
        heat < ymin && (ymin = heat); heat > ymax && (ymax = heat)
    end
    if ymax > ymin
        # Ref: acefc.f90:16191-16226.  Both axes log (ascll); `4 0 2 1/` (no xtag);
        # log thinning thin=10**(log10(xmax/xmin)/nden); data 1p,2e14.6.
        xmin, xmax = _ascll(xmin, xmax)
        ymin < ymax/scale && (ymin = ymax/scale)
        ymin, ymax = _ascll(ymin, ymax)
        pr(@sprintf("1%3d/", iwcol))
        pr("$(qu)<$(title)>$(qu)/")
        pr("$(qu)<h>eating$(qu)/")
        pr("4 0 2 1/")
        pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(1.0))/")
        pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
        pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(1.0))/")
        pr("$(qu)<h>eating (<m>e<v>/reaction)$(qu)/")
        pr("/"); pr("/")
        pr("$(qu)heating$(qu)/"); pr("0/")
        thin = ten^(log10(xmax/xmin)/nden)
        xlast = small
        for i in 1:nes
            e = e_i(i)
            if nes <= nden || e >= thin*xlast
                heat = heat_i(i)
                heat < hmin && (heat = hmin)
                heat < ymin && (heat = ymin)
                pr("$(_e14_6(e))$(_e14_6(heat))/")
                xlast = e
            end
        end
        pr("/")
    end

    # log-log damage (MT=444).  Ref: 16228-16288.
    _find_mt444 = () -> begin
        for i in 1:ntr
            if Int(round(xss[mtr + i - 1])) == 444
                return Int(round(xss[lsig + i - 1] + sig - 1))
            end
        end
        return 0
    end
    _find_mt444() != 0 && _aplots_unported("MT=444 log-log damage plot not yet ported")

    # log-log non-threshold reactions.  Ref: 16290-16398.  For ntr=0 there are
    # no reactions, so nlev stays 0 and nothing is emitted.
    _aplots_nonthreshold!(io, xss, nes, ntr, mtr, lsig, sig, izai, qu, title,
                          iwcol, nden, scale)

    # ----------------------------------------------------------------
    # PAGE: lin-lin principal cross sections.  Ref: acefc.f90:16400-16515.
    # ----------------------------------------------------------------
    xmin = big; xmax = 0.0; ymin = big; ymax = -big
    for i in 1:nes
        e = e_i(i)
        if e > 2
            tot = tot_i(i); abso = abso_i(i); elas = elas_i(i)
            e < xmin && (xmin = e); e > xmax && (xmax = e)
            tot < ymin && (ymin = tot); tot > ymax && (ymax = tot)
            abso < ymin && (ymin = abso); abso > ymax && (ymax = abso)
            elas < ymin && (ymin = elas); elas > ymax && (ymax = elas)
            if gpd != 0
                gp = gprod_i(i)
                gp < ymin && (ymin = gp); gp > ymax && (ymax = gp)
            end
        end
    end
    ymin = 0.0
    xmin, xmax, major, _ = _ascle(4, xmin, xmax)
    xstep = (xmax - xmin)/major
    ymin, ymax, major, _ = _ascle(4, ymin, ymax)
    ystep = (ymax - ymin)/major
    pr(@sprintf("1%3d/", iwcol))
    pr("$(qu)<$(title)>$(qu)/")
    pr("$(qu)<p>rincipal cross sections$(qu)/")
    xtag = 35*xmin/100 + 65*xmax/100
    ytag = ymin/10 + 9*ymax/10
    pr("1 0 2 1$(_e12_4_no1p(xtag))$(_e12_4_no1p(ytag))/")
    pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(xstep))/")
    pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
    pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(ystep))/")
    pr("$(qu)<c>ross section (barns)$(qu)/")
    pr("/"); pr("/")
    thin = (xmax - xmin)/nden
    # total (test = 1/5)
    pr("$(qu)total$(qu)/"); pr("0/")
    xlast = small
    for i in 1:nes
        e = e_i(i)
        test = 1.0/5
        if e >= test
            if nes <= nden || i == nes || e >= xlast + thin
                pr("$(_e14_6(e))$(_e14_6(tot_i(i)))/")
                xlast = e
            end
        end
    end
    pr("/")
    # absorption (test = 1/10)
    pr("2/"); pr("/"); pr("0 0 0 3/")
    pr("$(qu)absorption$(qu)/"); pr("0/")
    xlast = small
    for i in 1:nes
        e = e_i(i)
        test = 1.0/10
        if e >= test
            if nes <= nden || e >= xlast + thin || i == nes
                pr("$(_e14_6(e))$(_e14_6(abso_i(i)))/")
                xlast = e
            end
        end
    end
    pr("/")
    # elastic (test = 1/10)
    pr("3/"); pr("/"); pr("0 0 0 2/")
    pr("$(qu)elastic$(qu)/"); pr("0/")
    xlast = small
    for i in 1:nes
        e = e_i(i)
        test = 1.0/10
        if e >= test
            if nes <= nden || e >= xlast + thin || i == nes
                pr("$(_e14_6(e))$(_e14_6(elas_i(i)))/")
                xlast = e
            end
        end
    end
    pr("/")
    if gpd != 0
        pr("4/"); pr("/"); pr("0 0 0 1/")
        pr("$(qu)gamma production$(qu)/"); pr("0/")
        xlast = small
        for i in 1:nes
            e = e_i(i)
            test = 1.0/10
            if e >= test
                if nes <= nden || e >= xlast + thin
                    pr("$(_e14_6(e))$(_e14_6(gprod_i(i)))/")
                    xlast = e
                end
            end
        end
        pr("/")
    end

    # ----------------------------------------------------------------
    # lin-lin heating / damage / non-threshold pages.  Guarded; skip for T50.
    # Ref: 16544-16793.
    # ----------------------------------------------------------------
    xmin = big; xmax = 0.0; ymin = big; ymax = -big
    for i in 1:nes
        e = e_i(i); heat = heat_i(i)
        e < xmin && (xmin = e); e > xmax && (xmax = e)
        heat < ymin && (ymin = heat)
        (heat > ymax && e > 1.0) && (ymax = heat)
    end
    if ymin != 0.0 || ymax != 0.0
        # Ref: acefc.f90:16557-16595.  ascle(4,...) both axes; `1 0 2 1/` (no xtag);
        # energy threshold e>=0.2 (test=1/5); linear thinning; always write last
        # point (i==nes); NO hmin floor.
        xmin, xmax, major, _ = _ascle(4, xmin, xmax)
        xstep = (xmax - xmin)/major
        ymin, ymax, major, _ = _ascle(4, ymin, ymax)
        ystep = (ymax - ymin)/major
        pr(@sprintf("1%3d/", iwcol))
        pr("$(qu)<$(title)>$(qu)/")
        pr("$(qu)<h>eating$(qu)/")
        pr("1 0 2 1/")
        pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(xstep))/")
        pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
        pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(ystep))/")
        pr("$(qu)<h>eating (<m>e<v>/reaction)$(qu)/")
        pr("/"); pr("/")
        pr("$(qu)heating$(qu)/"); pr("0/")
        thin = (xmax - xmin)/nden
        xlast = small
        for i in 1:nes
            e = e_i(i)
            test = 1.0/5
            if e >= test
                if nes <= nden || e >= xlast + thin || i == nes
                    pr("$(_e14_6(e))$(_e14_6(heat_i(i)))/")
                    xlast = e
                end
            end
        end
        pr("/")
    end
    _find_mt444() != 0 && _aplots_unported("MT=444 lin-lin damage plot not yet ported")
    # lin-log non-threshold reactions (ntr=0 -> nothing).  Ref: 16670-16793.
    _aplots_nonthreshold_linlog!(io, xss, nes, ntr, mtr, lsig, sig, izai, qu,
                                 title, iwcol, nden)

    # ----------------------------------------------------------------
    # Delegated pages.  Ref: 16797-16815.
    # ----------------------------------------------------------------
    _aplotr(io, xss, ntr, mtr, lsig, sig, izai, qu, hk70, iwcol)   # 16797
    _aplopf(io, xss, ntr, mtr, lsig, sig, izai, qu, hk70, iwcol)   # 16800
    _aplof4(io, xss, nr, mtr, land, andb, izai, qu, hk70, iwcol)   # 16803
    # aplonu/aplodd/aplodn/aplopp/aploxp: guarded on nu>0, nr!=0, ndnf>0,
    # photon data, ntype>0 — none present for a bare charged-particle elastic
    # file.  Implemented below only where a test exercises them.
    Int(jxs[JXS_NU]) > 0 && _aplots_unported("aplonu (nubar) not yet ported")
    # aploxp — particle-production plots (heating, recoil, production XS, 3D
    # angular).  Ref: acefc.f90:16818 `if (ntype.gt.0) call aploxp(...)`.
    ntype = Int(nxs[NXS_NTYPE])
    if ntype > 0
        ptype = Int(jxs[JXS_PTYPE]); ntro = Int(jxs[JXS_NTRO]); ploct = Int(jxs[JXS_PLOCT])
        _aploxp(io, xss, nes, esz, ntype, ptype, ntro, ploct, izai, qu, hk70, iwcol)
    end
    # aplodd needs nr!=0 (reactions with secondary neutrons); aplof4 already
    # handled the angular pages.  aplopp (detailed photon) emits nothing when
    # there are no photon-production reactions (ntrp=0).
    Int(nxs[NXS_NDNF]) > 0 && _aplots_unported("aplodn (delayed neutron) not yet ported")

    # ----------------------------------------------------------------
    # Terminator.  Ref: acefc.f90:16827-16828.
    # ----------------------------------------------------------------
    pr("99/")
    pr("stop")
    nothing
end

# =========================================================================
# Non-threshold reaction pages (log-log and lin-log).  For ntr=0 the reaction
# loop never iterates, so nlev stays 0 and these emit nothing.  Ported in full
# so they work for tests with ntr>0.  Ref: 16290-16398 (log-log), 16670-16793.
# =========================================================================
function _aplots_nonthreshold!(io, xss, nes, ntr, mtr, lsig, sig, izai, qu,
                               title, iwcol, nden, scale)
    big = 1.0e10; ten = 10.0; small = 1.0e-12
    mtlast = 0
    pr(s) = print(io, s, '\n')
    nlev = 1
    while nlev > 0
        xmin = 1000.0; xmax = 0.0; ymin = 1000.0; ymax = 0.0
        nlev = 0
        for i in 1:ntr
            mt = Int(round(xss[mtr + i - 1]))
            k = Int(round(xss[lsig + i - 1] + sig - 1))
            n = Int(round(xss[k + 1])); iaa = Int(round(xss[k]))
            iflag = 0
            nlev == 5 && (iflag = 1)
            mt <= mtlast && (iflag = 1)
            mt > 207 && (iflag = 1)
            xss[iaa] > 1.0e-6 && (iflag = 1)
            if iflag == 0
                nlev += 1
                for j in 1:n
                    x = xss[iaa + j - 1]; y = xss[k + 2 + j - 1]
                    if y != 0.0 || j <= 1
                        x < xmin && (xmin = x); x > xmax && (xmax = x)
                        y < ymin && (ymin = y); y > ymax && (ymax = y)
                    end
                end
            end
        end
        if nlev != 0
            xmin, xmax = _ascll(xmin, xmax)
            ymin < ymax/scale && (ymin = ymax/scale)
            ymin, ymax = _ascll(ymin, ymax)
            thin = ten^(log10(xmax/xmin)/nden)
            pr(@sprintf("1%3d/", iwcol))
            pr("$(qu)<$(title)>$(qu)/")
            pr("$(qu)<n>on-threshold reactions$(qu)/")
            pr("4 0 2 1/")
            pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(1.0))/")
            pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
            pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(1.0))/")
            pr("$(qu)<c>ross section (barns)$(qu)/")
            pr("/")
            nlev = 0; mtl = 0
            for i in 1:ntr
                mt = Int(round(xss[mtr + i - 1]))
                k = Int(round(xss[lsig + i - 1] + sig - 1))
                n = Int(round(xss[k + 1])); iaa = Int(round(xss[k]))
                iflag = 0
                nlev == 5 && (iflag = 1)
                mt <= mtlast && (iflag = 1)
                mt > 207 && (iflag = 1)
                xss[iaa] > 1.0e-6 && (iflag = 1)
                if iflag == 0
                    mtl = mt; nlev += 1
                    nlev > 1 && pr(@sprintf("%2d/", nlev))
                    nlev > 1 && pr("/")
                    icurv = (nlev - 1) % 5
                    pr(@sprintf("0 0 0%2d/", icurv))
                    pr("$(qu)$(_aplots_mtname(mt, izai))$(qu)/")
                    pr("0/")
                    xlast = small
                    for j in 1:n
                        x = xss[iaa + j - 1]
                        if n <= nden || x >= xlast + thin || j == n
                            y = xss[k + 2 + j - 1]
                            y < ymin && (y = ymin)
                            pr("$(_e14_6(x))$(_e14_6(y))/")
                            xlast = x
                        end
                    end
                    pr("/")
                end
            end
            mtlast = mtl
        end
    end
end

function _aplots_nonthreshold_linlog!(io, xss, nes, ntr, mtr, lsig, sig, izai,
                                      qu, title, iwcol, nden)
    small = 1.0e-12
    mtlast = 0
    pr(s) = print(io, s, '\n')
    idone = 0
    while idone == 0
        xmin = 1000.0; xmax = 0.0; ymin = 1000.0; ymax = 0.0; nlev = 0
        for i in 1:ntr
            mt = Int(round(xss[mtr + i - 1]))
            k = Int(round(xss[lsig + i - 1] + sig - 1))
            n = Int(round(xss[k + 1])); iaa = Int(round(xss[k]))
            iflag = 0
            nlev == 5 && (iflag = 1)
            mt <= mtlast && (iflag = 1)
            mt > 207 && (iflag = 1)
            xss[iaa] > 1.0e-6 && (iflag = 1)
            if iflag == 0
                nlev += 1
                for j in 1:n
                    x = xss[iaa + j - 1]; y = xss[k + 2 + j - 1]
                    test = 1.0/5; f2 = 0
                    x < test && (f2 = 1)
                    (y == 0.0 && j > 1) && (f2 = 1)
                    if f2 == 0
                        x < xmin && (xmin = x); x > xmax && (xmax = x)
                        y < ymin && (ymin = y)
                        (y > ymax && x > 1.0) && (ymax = y)
                    end
                end
            end
        end
        if nlev == 0 || ymax == 0.0
            idone = 1
        else
            xmin, xmax, major, _ = _ascle(4, xmin, xmax)
            xstep = (xmax - xmin)/major
            ymin, ymax = _ascll(ymin, ymax)
            if ymax == ymin
                idone = 1
            else
                thin = (xmax - xmin)/nden
                pr(@sprintf("1%3d/", iwcol))
                pr("$(qu)<$(title)>$(qu)/")
                pr("$(qu)<n>on-threshold reactions$(qu)/")
                pr("2 0 2 1/")
                pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(xstep))/")
                pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
                pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(1.0))/")
                pr("$(qu)<c>ross section (barns)$(qu)/")
                pr("/")
                nlev = 0; mtl = 0
                for i in 1:ntr
                    mt = Int(round(xss[mtr + i - 1]))
                    k = Int(round(xss[lsig + i - 1] + sig - 1))
                    n = Int(round(xss[k + 1])); iaa = Int(round(xss[k]))
                    iflag = 0
                    nlev == 5 && (iflag = 1)
                    mt <= mtlast && (iflag = 1)
                    mt > 207 && (iflag = 1)
                    xss[iaa] > 1.0e-6 && (iflag = 1)
                    if iflag == 0
                        mtl = mt; nlev += 1
                        nlev > 1 && pr(@sprintf("%2d/", nlev))
                        nlev > 1 && pr("/")
                        icurv = (nlev - 1) % 5
                        pr(@sprintf("0 0 0%2d/", icurv))
                        pr("$(qu)$(_aplots_mtname(mt, izai))$(qu)/")
                        pr("0/")
                        xlast = small
                        for j in 1:n
                            x = xss[iaa + j - 1]; test = 1.0/5
                            if x >= test
                                if n <= nden || x >= xlast + thin || j == n
                                    y = xss[k + 2 + j - 1]
                                    y < ymin && (y = ymin)
                                    pr("$(_e14_6(x))$(_e14_6(y))/")
                                    xlast = x
                                end
                            end
                        end
                        pr("/")
                    end
                end
                mtlast = mtl
            end
        end
    end
end

# =========================================================================
# Caller-side incident-particle name swap.  After `mtname`, the aplots/aplotr
# call sites overwrite the second character of a `(...)` reaction name with the
# incident-particle letter (e.g. `(n,p*0)` -> `(d,p*0)` for a deuteron file).
# This is distinct from mtname's own izai>1 substitution (acecm.f90:209-221).
# Ref: acefc.f90:17056-17071 (identical block at 16365-16378, 16754-16767,
# 16933-16946).
# =========================================================================
function _aplots_name_izai(name::String, izai::Int)
    (isempty(name) || name[1] != '(') && return name
    c = izai == 1    ? 'n' : izai == 1001 ? 'p' : izai == 1002 ? 'd' :
        izai == 1003 ? 't' : izai == 2003 ? 's' : izai == 2004 ? 'a' : name[2]
    return string(name[1], c, name[3:end])
end

# =========================================================================
# aplotr — inelastic levels + threshold reactions.  Ref: 16832-17094.
# For ntr=0 both loops never iterate -> emits nothing.
# =========================================================================
function _aplotr(io, xss, ntr, mtr, lsig, sig, izai, qu, hk70, iwcol)
    small = 1.0e-12; xsmin = 1.0e-6; nden = 4000
    pr(s) = print(io, s, '\n')
    _it(s) = (j = 1; for i in 1:min(70, length(s)); s[i] != ' ' && (j = i); end; j)
    title = hk70[1:_it(hk70)]
    ilev = 0; ilast = 0

    # inelastic levels (only for izai==1 neutrons)
    while ilast == 0
        xmin = 1000.0; xmax = 0.0; ymin = 1000.0; ymax = 0.0; nlev = 0
        for i in 1:ntr
            mt = Int(round(xss[mtr + i - 1]))
            if mt >= 51 + ilev && mt <= 90 && izai == 1
                if nlev < 5
                    nlev += 1
                    k = Int(round(xss[lsig + i - 1] + sig - 1))
                    n = Int(round(xss[k + 1])); iaa = Int(round(xss[k]))
                    ylast = -1.0
                    for j in 1:n
                        x = xss[iaa + j - 1]; y = xss[k + 2 + j - 1]
                        if y != ylast
                            x < xmin && (xmin = x); x > xmax && (xmax = x)
                            y < ymin && (ymin = y); y > ymax && (ymax = y)
                        end
                        ylast = y
                    end
                end
            end
        end
        if nlev == 0
            ilast = 1
        else
            _aplots_unported("aplotr inelastic-level plotting (izai=1) not yet ported")
        end
        ilev += nlev
    end

    # threshold reactions
    idone = 0; mtlast = 0
    while idone == 0
        xmin = 1000.0; xmax = 0.0; ymin = 1000.0; ymax = 0.0; nlev = 0
        for i in 1:ntr
            mt = Int(round(xss[mtr + i - 1]))
            k = Int(round(xss[lsig + i - 1] + sig - 1))
            n = Int(round(xss[k + 1])); iaa = Int(round(xss[k]))
            iflag = 0
            nlev == 5 && (iflag = 1)
            mt < 5 && (iflag = 1)
            (mt >= 18 && mt <= 21) && (iflag = 1)
            mt == 38 && (iflag = 1)
            (mt >= 50 && mt < 91 && izai == 1) && (iflag = 1)
            (mt > 207 && mt < 600) && (iflag = 1)
            mt <= mtlast && (iflag = 1)
            xss[iaa] < xsmin && (iflag = 1)
            if iflag == 0
                nlev += 1
                for j in 1:n
                    x = xss[iaa + j - 1]; y = xss[k + 2 + j - 1]
                    if y != 0.0 || j <= 1
                        x < xmin && (xmin = x); x > xmax && (xmax = x)
                        y < ymin && (ymin = y); y > ymax && (ymax = y)
                    end
                end
                if mt < 203 && Int(round(xss[mtr + i])) >= 203 &&
                   Int(round(xss[mtr + i])) <= 207
                    nlev = 5
                end
            end
        end
        if nlev == 0 || ymax == 0.0
            idone = 1
        else
            # --- emit the threshold-reaction page.  Ref: acefc.f90:17012-17090.
            xmin, xmax, major, _ = _ascle(4, xmin, xmax)   # 17012
            xstep = (xmax - xmin) / major
            ymin, ymax, major, _ = _ascle(4, ymin, ymax)   # 17014
            ystep = (ymax - ymin) / major
            thin = (xmax - xmin) / nden                    # 17016
            # page header records (17017-17029)
            pr(@sprintf("1%3d/", iwcol))
            pr("$(qu)<$(title)>$(qu)/")
            pr("$(qu)<t>hreshold reactions$(qu)/")
            pr("1 0 2 1/")
            pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(xstep))/")
            pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
            pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(ystep))/")
            pr("$(qu)<c>ross section (barns)$(qu)/")
            pr("/")
            # per-curve loop, same iflag filter as the scan (17030-17089)
            nlev = 0; mtl = 0
            for i in 1:ntr
                mt = Int(round(xss[mtr + i - 1]))
                k = Int(round(xss[lsig + i - 1] + sig - 1))
                n = Int(round(xss[k + 1])); iaa = Int(round(xss[k]))
                iflag = 0
                nlev == 5 && (iflag = 1)
                mt < 5 && (iflag = 1)
                (mt >= 18 && mt <= 21) && (iflag = 1)
                mt == 38 && (iflag = 1)
                (mt >= 50 && mt < 91 && izai == 1) && (iflag = 1)
                (mt > 207 && mt < 600) && (iflag = 1)
                mt <= mtlast && (iflag = 1)
                xss[iaa] < xsmin && (iflag = 1)
                if iflag == 0
                    mtl = mt; nlev += 1
                    nlev > 1 && pr(@sprintf("%2d/", nlev))   # 17048-17049
                    nlev > 1 && pr("/")
                    icurv = (nlev - 1) % 5                    # 17050
                    pr(@sprintf("0 0 0%2d/", icurv))         # iwcol=3 color mode (17054)
                    name = _aplots_name_izai(_aplots_mtname(mt, izai), izai)  # 17056-17071
                    pr("$(qu)$(name)$(qu)/")                  # 17072
                    pr("0/")                                  # 17073
                    xlast = small
                    for j in 1:n
                        x = xss[iaa + j - 1]
                        if n <= nden || x >= xlast + thin || j == 1 || j == n  # 17078-17079
                            y = xss[k + 2 + j - 1]
                            pr("$(_e14_6(x))$(_e14_6(y))/")  # 17081
                            xlast = x
                        end
                    end
                    pr("/")                                   # 17085 end-curve
                    if mt < 203 && Int(round(xss[mtr + i])) >= 203 &&
                       Int(round(xss[mtr + i])) <= 207        # 17086-17087
                        nlev = 5
                    end
                end
            end
            mtlast = mtl                                      # 17090
        end
    end
end

# =========================================================================
# aplopf — higher fission reactions.  Ref: 17096-17203.  No fission (mt 20/21/
# 38 absent) -> nlev=0 -> emits nothing.
# =========================================================================
function _aplopf(io, xss, ntr, mtr, lsig, sig, izai, qu, hk70, iwcol)
    nlev = 0; ymax = 0.0
    for i in 1:ntr
        mt = Int(round(xss[mtr + i - 1]))
        iflag = 0
        mt < 20 && (iflag = 1)
        (mt > 21 && mt < 38) && (iflag = 1)
        mt > 38 && (iflag = 1)
        if iflag == 0
            nlev += 1
            k = Int(round(xss[lsig + i - 1] + sig - 1))
            n = Int(round(xss[k + 1])); iaa = Int(round(xss[k]))
            for j in 1:n
                y = xss[k + 2 + j - 1]
                (y != 0.0 || j <= 1) && y > ymax && (ymax = y)
            end
        end
    end
    if nlev != 0 && ymax != 0.0
        _aplots_unported("aplopf higher-fission plotting not yet ported")
    end
end

# =========================================================================
# aplof4 — angular distributions in 3D form.  Ref: 17205-17428.
# Ports the elastic (n=0) tabulated perspective path (k<0) used by T50.
# =========================================================================
function _aplof4(io, xss, nr, mtr, land, andb, izai, qu, hk70, iwcol)
    eps = 1.0e-5; big = 1.0e10; dn = 0.99; up = 1.01; one = 1.0
    pr(s) = print(io, s, '\n')
    _it(s) = (j = 1; for i in 1:min(70, length(s)); s[i] != ' ' && (j = i); end; j)
    it = 0

    nr1 = nr + 1
    for nn in 1:nr1
        n = nn - 1
        na0 = Int(round(xss[land + n]))
        na0 <= 0 && continue
        if n == 0
            mt = 2
            name = "elastic   "   # Fortran character(10) name='elastic' (acefc.f90:17242)
        else
            mt = abs(Int(round(xss[mtr + n - 1])))
            name = _aplots_mtname(mt, izai)
        end
        na = na0 + andb - 1
        ne = Int(round(xss[na]))
        nb = na + ne
        k = Int(round(xss[nb + 1]))

        if k >= 0
            # equiprobable-contours path — not exercised by current tests
            _aplots_unported("aplof4 equiprobable-contour angular plot (k>=0) not yet ported")
        else
            # perspective view of tabulated distribution.  The aplof4 elastic
            # subtitle keeps the 10-char name verbatim (trailing blanks).
            _aplof4_perspective!(io, xss, na, ne, nb, andb,
                                 "angular distribution for $(name)", hk70, iwcol, qu)
        end
    end
end

# =========================================================================
# _aplof4_perspective! — the perspective view of a tabulated angular
# distribution (the `k<0` body of aplof4, acefc.f90:17268-17428).  Factored
# out so aploxp's particle-production angular pages (acefc.f90:19347-19455,
# byte-identical apart from the subtitle and the andb→andh locator) reuse it.
# `subtitle` is the full text between the quotes on the third card.
# =========================================================================
function _aplof4_perspective!(io, xss, na, ne, nb, andb, subtitle::AbstractString,
                              hk70, iwcol, qu)
    eps = 1.0e-5; big = 1.0e10; dn = 0.99; up = 1.01; one = 1.0
    pr(s) = print(io, s, '\n')
    _it(s) = (j = 1; for i in 1:min(70, length(s)); s[i] != ' ' && (j = i); end; j)

    xmin = -1.0; xmax = 1.0; xstep = 1.0/2
    ymin = xss[na + 1]; ymax = xss[na + ne]
    itwo = 0
    if ne > 4
        test = 3.0
        ymax > test*xss[na + ne - 1] && (ymax = xss[na + ne - 1])
        test = 50.0
        if ymax >= test
            break_ = 20.0
            (ymin < dn*break_ && ymax > up*break_) && (itwo = 1)
        end
    end
    break_ = 20.0
    while itwo >= 0 && itwo <= 2
        if itwo == 1
            ymax = break_
        elseif itwo == 2
            ymin = break_; ymax = xss[na + ne]
        end
        ymin, ymax, major, _ = _ascle(4, ymin, ymax)
        ystep = (ymax - ymin)/major
        zmin = big; zmax = 0.0
        for i in 1:ne
            e = xss[na + i]
            if e <= (1 + eps)*ymax && e >= ymin
                kk = Int(round(abs(xss[nb + i]))) + andb - 1
                np = Int(round(xss[kk + 1]))
                kk += 1
                for j in 1:np
                    pp = xss[kk + np + j]
                    pp < zmin && (zmin = pp)
                    pp > zmax && (zmax = pp)
                end
            end
        end
        rat = 100000.0
        zmin <= 0.0 && (zmin = zmax/rat)
        zmax/zmin > rat && (zmin = zmax/rat)
        zmin, zmax = _ascll(zmin, zmax)
        if zmin == zmax
            zmax = 2*zmax; zmin = zmax/10
        end
        it = _it(hk70)
        title = hk70[1:it]
        pr(@sprintf("1%3d/", iwcol))
        pr("$(qu)<$(title)>$(qu)/")
        pr("$(qu)$(subtitle)$(qu)/")
        pr("-1 2/")
        pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(xstep))/")
        pr("$(qu)<c>osine$(qu)/")
        pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(ystep))/")
        pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
        pr("$(_e12_3(zmin))$(_e12_3(zmax))$(_e12_3(one))/")
        pr("$(qu)<p>rob/<C>os$(qu)/")
        pr("/")
        pr(" 15. -15. 15. -2.5 6.5 2.5/")
        pr("1/")
        stepm = (ymax - ymin)/150
        elast = 0.0
        for i in 1:ne
            e = xss[na + i]
            iflag = 0
            e > (1 + eps)*ymax && (iflag = 1)
            e < ymin && (iflag = 1)
            (i > 1 && i < ne && ne > 150 && e - elast < stepm) && (iflag = 1)
            if iflag == 0
                pr("$(_e14_6(e))/")
                kk = Int(round(abs(xss[nb + i]))) + andb - 1
                intt = Int(round(xss[kk]))
                np = Int(round(xss[kk + 1]))
                kk += 1
                ylast = zmin
                for j in 1:np
                    cc = xss[kk + j]
                    pp = xss[kk + np + j]
                    pp < zmin && (pp = zmin)
                    intt == 1 && pr("$(_e14_6(cc))$(_e14_6(ylast))/")
                    pr("$(_e14_6(cc))$(_e14_6(pp))/")
                    ylast = pp
                end
                pr("/")
                elast = e
            end
        end
        pr("/")
        itwo == 0 && (itwo = 3)
        itwo == 2 && (itwo = 3)
        itwo == 1 && (itwo = 2)
    end
end

# =========================================================================
# aploxp — particle-production plots.  Ref: acefc.f90:18787-19460.
# Emits, in order: (C1) lin-lin particle heating contributions, one curve per
# production type; (C2) lin-lin recoil heating = total heating − Σ particle
# heating; (C3) lin-lin particle production cross sections; then per-type/per-MT
# 3D angular-distribution pages (C5).  Locators come from the production arrays
# ptype/ntro/ploct (jxs(30..32)); per type i the 10-word locator block lives at
# ploct+10*(i-1)+k.  Energy column is the ESZ grid, indexed via iaa=xss(hpd),
# naa=xss(hpd+1): e = xss(esz+iaa-2+i).
# =========================================================================
function _aploxp(io, xss, nes, esz, ntype, ptype, ntro, ploct, izai, qu, hk70, iwcol)
    big = 1.0e10; small = 1.0e-12; nden = 4000
    pr(s) = print(io, s, '\n')
    _it(s) = (j = 1; for i in 1:min(70, length(s)); s[i] != ' ' && (j = i); end; j)
    title = hk70[1:_it(hk70)]

    # ipt -> single-word label (C1 heating / C3 production).  Ref: 18866-18872.
    plabel(ipt) = ipt == 1 ? "neutrons" : ipt == 2 ? "photons" : ipt == 9 ? "protons" :
        ipt == 31 ? "deuterons" : ipt == 32 ? "tritons" : ipt == 33 ? "he-3" :
        ipt == 34 ? "alphas" : error("aploxp: unknown ipt=$ipt")
    # ipt -> production-XS label (photons divided by 5).  Ref: 19058-19064.
    plabel_xs(ipt) = ipt == 2 ? "photons/5" : plabel(ipt)
    # ipt -> singular particle word for the angular subtitle.  Ref: 19387-19410.
    pword(ipt) = ipt == 1 ? "neutron" : ipt == 9 ? "proton" : ipt == 31 ? "deuteron" :
        ipt == 32 ? "triton" : ipt == 33 ? "3he" : ipt == 34 ? "alpha" :
        error("aploxp: unknown ipt=$ipt for angular subtitle")

    # locator accessors for production type i (1-based).  Ref: 19124-19133.
    hpd_i(i)   = Int(round(xss[ploct + 10*(i-1)]))
    mtrh_i(i)  = Int(round(xss[ploct + 10*(i-1) + 1]))
    nmtr_i(i)  = Int(round(xss[ntro + i - 1]))
    lsigh_i(i) = Int(round(xss[ploct + 10*(i-1) + 3]))
    sigh_i(i)  = Int(round(xss[ploct + 10*(i-1) + 4]))
    landh_i(i) = Int(round(xss[ploct + 10*(i-1) + 5]))
    andh_i(i)  = Int(round(xss[ploct + 10*(i-1) + 6]))
    ldlwh_i(i) = Int(round(xss[ploct + 10*(i-1) + 7]))
    dlwh_i(i)  = Int(round(xss[ploct + 10*(i-1) + 8]))
    ipt_i(i)   = Int(round(xss[ptype + i - 1]))

    # ===== C1: lin-lin particle heating contributions.  Ref: 18816-18934. =====
    xmin = big; xmax = 0.0; ymin = big; ymax = -big
    for i in 1:ntype
        hpd = hpd_i(i); iaa = Int(round(xss[hpd])); naa = Int(round(xss[hpd+1]))
        test = 1.0/5
        for ie in 1:naa
            e = xss[esz + iaa + ie - 2]
            if e > test
                xs = xss[hpd + 1 + naa + ie]
                e < xmin && (xmin = e); e > xmax && (xmax = e)
                xs < ymin && (ymin = xs); xs > ymax && (ymax = xs)
            end
        end
    end
    if ymax != 0.0
        xmin, xmax, major, _ = _ascle(4, xmin, xmax); xstep = (xmax-xmin)/major
        ymin, ymax, major, _ = _ascle(4, ymin, ymax); ystep = (ymax-ymin)/major
        pr(@sprintf("1%3d/", iwcol))
        pr("$(qu)<$(title)>$(qu)/")
        pr("$(qu)<p>article heating contributions$(qu)/")
        xtag = 95*xmin/100 + 5*xmax/100
        ytag = ymin/10 + 9*ymax/10
        pr("1 0 2 1$(_e12_4_no1p(xtag))$(_e12_4_no1p(ytag))/")
        pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(xstep))/")
        pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
        pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(ystep))/")
        pr("$(qu)<m>e<v>/collision$(qu)/")
        pr("/")
        thin = (xmax - xmin)/nden
        for ii in 1:ntype
            ipt = ipt_i(ii); hpd = hpd_i(ii)
            iaa = Int(round(xss[hpd])); naa = Int(round(xss[hpd+1]))
            if ii > 1
                pr("2/"); pr("/")
            end
            pr(@sprintf("0 0 0%2d/", ii-1))
            pr("$(qu)$(plabel(ipt))$(qu)/")
            pr("0/")
            xlast = small; test = 1.0/5
            for i in 1:naa
                e = xss[esz + iaa - 2 + i]
                if e > test
                    xs = xss[hpd + 1 + naa + i]
                    if naa <= nden || e >= xlast + thin || i == naa
                        pr("$(_e14_6(e))$(_e14_6(xs))/")
                        xlast = e
                    end
                end
            end
            pr("/")
        end
    end

    # ===== C2: lin-lin recoil heating = esz heating − Σ particle heating. =====
    # Ref: 18936-19004.
    xmin = big; xmax = 0.0; ymin = big; ymax = -big
    for i in 1:nes
        e = xss[esz - 1 + i]
        heat = xss[esz + 4*nes - 1 + i]
        for j in 1:ntype
            hpd = hpd_i(j); iaa = Int(round(xss[hpd])); naa = Int(round(xss[hpd+1]))
            ie = i - iaa - 1
            (ie >= 1 && ie <= naa) && (heat -= xss[hpd + 1 + naa + ie])
        end
        e < xmin && (xmin = e); e > xmax && (xmax = e)
        heat < ymin && (ymin = heat); heat > ymax && (ymax = heat)
    end
    if ymin != 0.0 || ymax != 0.0
        (ymin < 0.0 && ymax < -ymin/2) && (ymax = -ymin/2)
        xmin, xmax, major, _ = _ascle(4, xmin, xmax); xstep = (xmax-xmin)/major
        ymin, ymax, major, _ = _ascle(4, ymin, ymax); ystep = (ymax-ymin)/major
        pr(@sprintf("1%3d/", iwcol))
        pr("$(qu)<$(title)>$(qu)/")
        pr("$(qu)<r>ecoil <h>eating$(qu)/")
        pr("1 0 2 1/")
        pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(xstep))/")
        pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
        pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(ystep))/")
        pr("$(qu)<h>eating (<m>e<v>/reaction)$(qu)/")
        pr("/"); pr("/")
        pr("$(qu)recoil heating$(qu)/"); pr("0/")
        thin = (xmax - xmin)/nden
        xlast = small
        for i in 1:nes
            e = xss[esz - 1 + i]; test = 1.0/5
            if e >= test
                if nes <= nden || e >= xlast + thin || i == nes
                    heat = xss[esz + 4*nes - 1 + i]
                    for k in 1:ntype
                        hpd = hpd_i(k); iaa = Int(round(xss[hpd])); naa = Int(round(xss[hpd+1]))
                        ie = i - iaa - 1
                        (ie >= 1 && ie <= naa) && (heat -= xss[hpd + 1 + naa + ie])
                    end
                    pr("$(_e14_6(e))$(_e14_6(heat))/")
                    xlast = e
                end
            end
        end
        pr("/")
    end

    # ===== C3: lin-lin particle production cross sections (always fires). =====
    # Ref: 19006-19119.  xs = xss(hpd+1+i); photons (ipt==2) divided by 5.
    xmin = big; xmax = 0.0; ymin = big; ymax = -big
    for i in 1:ntype
        ipt = ipt_i(i); hpd = hpd_i(i)
        iaa = Int(round(xss[hpd])); naa = Int(round(xss[hpd+1]))
        test = 1.0/5
        for ie in 1:naa
            e = xss[esz + iaa + ie - 2]
            if e > test
                xs = xss[hpd + 1 + ie]
                ipt == 2 && (xs /= 5)
                e < xmin && (xmin = e); e > xmax && (xmax = e)
                xs < ymin && (ymin = xs); xs > ymax && (ymax = xs)
            end
        end
    end
    xmin, xmax, major, _ = _ascle(4, xmin, xmax); xstep = (xmax-xmin)/major
    ymin, ymax, major, _ = _ascle(4, ymin, ymax); ystep = (ymax-ymin)/major
    pr(@sprintf("1%3d/", iwcol))
    pr("$(qu)<$(title)>$(qu)/")
    pr("$(qu)<p>article production cross sections$(qu)/")
    xtag = 95*xmin/100 + 5*xmax/100
    xmax > 30.0 && (xtag = 35*xmin/100 + 65*xmax/100)
    ytag = ymin/10 + 9*ymax/10
    pr("1 0 2 1$(_e12_4_no1p(xtag))$(_e12_4_no1p(ytag))/")
    pr("$(_e12_3(xmin))$(_e12_3(xmax))$(_e12_3(xstep))/")
    pr("$(qu)<e>nergy (<m>e<v>)$(qu)/")
    pr("$(_e12_3(ymin))$(_e12_3(ymax))$(_e12_3(ystep))/")
    pr("$(qu)<c>ross section (barns)$(qu)/")
    pr("/")
    thin = (xmax - xmin)/nden
    for ii in 1:ntype
        ipt = ipt_i(ii); hpd = hpd_i(ii)
        iaa = Int(round(xss[hpd])); naa = Int(round(xss[hpd+1]))
        if ii > 1
            pr("2/"); pr("/")
        end
        pr(@sprintf("0 0 0%2d/", ii-1))
        pr("$(qu)$(plabel_xs(ipt))$(qu)/")
        pr("0/")
        xlast = small; test = 1.0/5
        for i in 1:naa
            e = xss[esz + iaa - 2 + i]
            if e > test
                xs = xss[hpd + 1 + i]
                ipt == 2 && (xs /= 5)
                if naa <= nden || e >= xlast + thin || i == naa
                    pr("$(_e14_6(e))$(_e14_6(xs))/")
                    xlast = e
                end
            end
        end
        pr("/")
    end

    # ===== C5: per-type per-MT 3D plots.  Ref: 19121-19458. =====
    for i in 1:ntype
        ipt = ipt_i(i)
        nmtr = nmtr_i(i)
        landh = landh_i(i); andh = andh_i(i)
        ldlwh = ldlwh_i(i); dlwh = dlwh_i(i)
        for imt in 1:nmtr
            # law of this reaction.  Ref: 19141-19143.
            l1 = Int(round(xss[ldlwh + imt - 1]))
            l2 = dlwh + l1 - 1
            law = Int(round(xss[l2 + 1]))
            if law == 4 || law == 44 || law == 61
                # 3D emission-spectrum page — no T53 reaction has these laws.
                _aplots_unported("aploxp law=$law emission-spectrum 3D plot not yet ported")
            end
            # angular distribution.  Ref: 19307-19456 (== aplof4 perspective body).
            na0 = Int(round(xss[landh + imt - 1]))
            if na0 > 0
                mt = abs(Int(round(xss[mtrh_i(i) + imt - 1])))
                name = _aplots_name_izai(_aplots_mtname(mt, izai), izai)
                na = na0 + andh - 1
                ne = Int(round(xss[na]))
                nb = na + ne
                # len_trim'd name + ' ' + particle word.  Ref: 19386-19410.
                nm = rstrip(name)
                subtitle = "angular distribution for $(nm) $(pword(ipt))"
                _aplof4_perspective!(io, xss, na, ne, nb, andh, subtitle,
                                     hk70, iwcol, qu)
            end
        end
    end
    nothing
end

# =========================================================================
# mtname — MT -> reaction name (acecm.f90:21-224), with charged-particle
# alternate names.  Returns the 10-char Fortran name verbatim (trailing blanks
# preserved, as the Fortran write emits `a` = the full character(10)).
# =========================================================================
# Full 500-entry `hndf` table, transcribed verbatim from acecm.f90:31-138.
const _APLOTS_MT_HNDF = String[
     "total     ","elastic   ","nonelastic","inelastic ","(n,x)     ",
     "(n,1/2*1) ","(n,1/2*2) ","(n,1/2*3) ","(n,1/2*4) ","(n,x)     ",
     "(n,2nd)   ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,2n)    ","(n,3n)    ","fission   ","(n,f)     ","(n,n*f)   ",
     "(n,2nf)   ","(n,n*)a   ","(n,n*)3a  ","(n,2n)a   ","(n,3n)a   ",
     "(n,2n)iso ","(n,abs)   ","(n,n*)p   ","(n,n*)2a  ","(n,2n)2a  ",
     "(n,x)     ","(n,n*)d   ","(n,n*)t   ","(n,n*)he3 ","(n,n*)d2a ",
     "(n,n*)t2a ","(n,4n)    ","(n,3nf)   ","(n,x)     ","(n,x)     ",
     "(n,2np)   ","(n,3np)   ","(n,x)     ","(n,n2p)   ","(n,npa)   ",
     "(n,2/2*1) ","(n,2/2*2) ","(n,2/2*3) ","(n,2/2*4) ","(n,n*0)   ",   # 50
     "(n,n*1)   ","(n,n*2)   ","(n,n*3)   ","(n,n*4)   ","(n,n*5)   ",
     "(n,n*6)   ","(n,n*7)   ","(n,n*8)   ","(n,n*9)   ","(n,n*10)  ",
     "(n,n*11)  ","(n,n*12)  ","(n,n*13)  ","(n,n*14)  ","(n,n*15)  ",
     "(n,n*16)  ","(n,n*17)  ","(n,n*18)  ","(n,n*19)  ","(n,n*20)  ",
     "(n,n*21)  ","(n,n*22)  ","(n,n*23)  ","(n,n*24)  ","(n,n*25)  ",
     "(n,n*26)  ","(n,n*27)  ","(n,n*28)  ","(n,n*29)  ","(n,n*30)  ",
     "(n,n*31)  ","(n,n*32)  ","(n,n*33)  ","(n,n*34)  ","(n,n*35)  ",
     "(n,n*36)  ","(n,n*37)  ","(n,n*38)  ","(n,n*39)  ","(n,n*40)  ",
     "(n,n*c)   ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,n*)gma ","(n,x)     ",   # 100
     "(n,parab) ","(n,gma)   ","(n,p)     ","(n,d)     ","(n,t)     ",
     "(n,he3)   ","(n,a)     ","(n,2a)    ","(n,3a)    ","(n,x)     ",
     "(n,2p)    ","(n,pa)    ","(n,t2a)   ","(n,d2a)   ","(n,pd)    ",
     "(n,pt)    ","(n,da)    ","(n,x)     ","(n,x)     ","(n,dest)  ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",   # 150
     "(n,x)     ","(n,5n)    ","(n,6n)    ","(n,2nt)   ","(n,ta)    ",
     "(n,4np)   ","(n,3nd)   ","(n,nda)   ","(n,2npa)  ","(n,7n)    ",
     "(n,8n)    ","(n,5np)   ","(n,6np)   ","(n,7np)   ","(n,4na)   ",
     "(n,5na)   ","(n,6na)   ","(n,7na)   ","(n,4nd)   ","(n,5nd)   ",
     "(n,6nd)   ","(n,3nt)   ","(n,4nt)   ","(n,5nt)   ","(n,6nt)   ",
     "(n,2nhe3) ","(n,3nhe3) ","(n,4nhe3) ","(n,3n2p)  ","(n,3n2a)  ",
     "(n,3npa)  ","(n,dt)    ","(n,npd)   ","(n,npt)   ","(n,ndt)   ",
     "(n,nphe3) ","(n,ndhe3) ","(n,nthe3) ","(n,nta)   ","(n,2n2p)  ",
     "(n,phe3)  ","(n,dhe3)  ","(n,he3a)  ","(n,4n2p)  ","(n,4n2a)  ",
     "(n,4npa)  ","(n,3p)    ","(n,n3p)   ","(n,3n2pa) ","(n,5n2p)  ",   # 200
     "(n,p*0)   ",
     "(n,p*1)   ","(n,p*2)   ","(n,p*3)   ","(n,p*4)   ","(n,p*5)   ",
     "(n,p*6)   ","(n,p*7)   ","(n,p*8)   ","(n,p*9)   ","(n,p*10)  ",
     "(n,p*11)  ","(n,p*12)  ","(n,p*13)  ","(n,p*14)  ","(n,p*15)  ",
     "(n,p*16)  ","(n,p*17)  ","(n,p*18)  ","(n,p*19)  ","(n,p*20)  ",
     "(n,p*21)  ","(n,p*22)  ","(n,p*23)  ","(n,p*24)  ","(n,p*25)  ",
     "(n,p*26)  ","(n,p*27)  ","(n,p*28)  ","(n,p*29)  ","(n,p*30)  ",
     "(n,p*31)  ","(n,p*32)  ","(n,p*33)  ","(n,p*34)  ","(n,p*35)  ",
     "(n,p*36)  ","(n,p*37)  ","(n,p*38)  ","(n,p*39)  ","(n,p*40)  ",
     "(n,p*41)  ","(n,p*42)  ","(n,p*43)  ","(n,p*44)  ","(n,p*45)  ",
     "(n,p*46)  ","(n,p*47)  ","(n,p*48)  ","(n,p*c)   ",                # 250
     "(n,d*0)   ",
     "(n,d*1)   ","(n,d*2)   ","(n,d*3)   ","(n,d*4)   ","(n,d*5)   ",
     "(n,d*6)   ","(n,d*7)   ","(n,d*8)   ","(n,d*9)   ","(n,d*10)  ",
     "(n,d*11)  ","(n,d*12)  ","(n,d*13)  ","(n,d*14)  ","(n,d*15)  ",
     "(n,d*16)  ","(n,d*17)  ","(n,d*18)  ","(n,d*19)  ","(n,d*20)  ",
     "(n,d*21)  ","(n,d*22)  ","(n,d*23)  ","(n,d*24)  ","(n,d*25)  ",
     "(n,d*26)  ","(n,d*27)  ","(n,d*28)  ","(n,d*29)  ","(n,d*30)  ",
     "(n,d*31)  ","(n,d*32)  ","(n,d*33)  ","(n,d*34)  ","(n,d*35)  ",
     "(n,d*36)  ","(n,d*37)  ","(n,d*38)  ","(n,d*39)  ","(n,d*40)  ",
     "(n,d*41)  ","(n,d*42)  ","(n,d*43)  ","(n,d*44)  ","(n,d*45)  ",
     "(n,d*46)  ","(n,d*47)  ","(n,d*48)  ","(n,d*c)   ",                # 300
     "(n,t*0)   ",
     "(n,t*1)   ","(n,t*2)   ","(n,t*3)   ","(n,t*4)   ","(n,t*5)   ",
     "(n,t*6)   ","(n,t*7)   ","(n,t*8)   ","(n,t*9)   ","(n,t*10)  ",
     "(n,t*11)  ","(n,t*12)  ","(n,t*13)  ","(n,t*14)  ","(n,t*15)  ",
     "(n,t*16)  ","(n,t*17)  ","(n,t*18)  ","(n,t*19)  ","(n,t*20)  ",
     "(n,t*21)  ","(n,t*22)  ","(n,t*23)  ","(n,t*24)  ","(n,t*25)  ",
     "(n,t*26)  ","(n,t*27)  ","(n,t*28)  ","(n,t*29)  ","(n,t*30)  ",
     "(n,t*31)  ","(n,t*32)  ","(n,t*33)  ","(n,t*34)  ","(n,t*35)  ",
     "(n,t*36)  ","(n,t*37)  ","(n,t*38)  ","(n,t*39)  ","(n,t*40)  ",
     "(n,t*41)  ","(n,t*42)  ","(n,t*43)  ","(n,t*44)  ","(n,t*45)  ",
     "(n,t*46)  ","(n,t*47)  ","(n,t*48)  ","(n,t*c)   ",                # 350
     "(n,he3*0) ",
     "(n,he3*1) ","(n,he3*2) ","(n,he3*3) ","(n,he3*4) ","(n,he3*5) ",
     "(n,he3*6) ","(n,he3*7) ","(n,he3*8) ","(n,he3*9) ","(n,he3*10)",
     "(n,he3*11)","(n,he3*12)","(n,he3*13)","(n,he3*14)","(n,he3*15)",
     "(n,he3*16)","(n,he3*17)","(n,he3*18)","(n,he3*19)","(n,he3*20)",
     "(n,he3*21)","(n,he3*22)","(n,he3*23)","(n,he3*24)","(n,he3*25)",
     "(n,he3*26)","(n,he3*27)","(n,he3*28)","(n,he3*29)","(n,he3*30)",
     "(n,he3*31)","(n,he3*32)","(n,he3*33)","(n,he3*34)","(n,he3*35)",
     "(n,he3*36)","(n,he3*37)","(n,he3*38)","(n,he3*39)","(n,he3*40)",
     "(n,he3*41)","(n,he3*42)","(n,he3*43)","(n,he3*44)","(n,he3*45)",
     "(n,he3*46)","(n,he3*47)","(n,he3*48)","(n,he3*c) ",                # 400
     "(n,a*0)   ",
     "(n,a*1)   ","(n,a*2)   ","(n,a*3)   ","(n,a*4)   ","(n,a*5)   ",
     "(n,a*6)   ","(n,a*7)   ","(n,a*8)   ","(n,a*9)   ","(n,a*10)  ",
     "(n,a*11)  ","(n,a*12)  ","(n,a*13)  ","(n,a*14)  ","(n,a*15)  ",
     "(n,a*16)  ","(n,a*17)  ","(n,a*18)  ","(n,a*19)  ","(n,a*20)  ",
     "(n,a*21)  ","(n,a*22)  ","(n,a*23)  ","(n,a*24)  ","(n,a*25)  ",
     "(n,a*26)  ","(n,a*27)  ","(n,a*28)  ","(n,a*29)  ","(n,a*30)  ",
     "(n,a*31)  ","(n,a*32)  ","(n,a*33)  ","(n,a*34)  ","(n,a*35)  ",
     "(n,a*36)  ","(n,a*37)  ","(n,a*38)  ","(n,a*39)  ","(n,a*40)  ",
     "(n,a*41)  ","(n,a*42)  ","(n,a*43)  ","(n,a*44)  ","(n,a*45)  ",
     "(n,a*46)  ","(n,a*47)  ","(n,a*48)  ","(n,a*c)   ",                # 450
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",   # 475
     "(n,2n*0)  ",
     "(n,2n*1)  ","(n,2n*2)  ","(n,2n*3)  ","(n,2n*4)  ","(n,2n*5)  ",
     "(n,2n*6)  ","(n,2n*7)  ","(n,2n*8)  ","(n,2n*9)  ","(n,2n*10) ",
     "(n,2n*11) ","(n,2n*12) ","(n,2n*13) ","(n,2n*14) ","(n,2n*15) ",
     "(n,2n*c)  ",
     "(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ","(n,x)     ",
     "(n,x)     ","(n,x)     ","(n,x)     "]                             # 500

"""
    _aplots_mtname(mt::Int, izai::Int) -> String

Return the 10-char reaction name for MT, mirroring acecm.f90 mtname for the
iverf>=6 path (the common case), plus the charged-particle alternates. The
returned string keeps the Fortran trailing-blank padding so the viewr write
`'(a,a,a,''/'')' qu,name,qu` is reproduced byte-for-byte.

Only the entries that current tests exercise are tabulated below; an
un-tabulated MT raises (Rule 6) so a missing case surfaces loudly instead of
silently emitting a wrong label.
"""
function _aplots_mtname(mt::Int, izai::Int)
    # iverf>=6 path (acecm.f90:157-167)
    name = ""
    if mt >= 201 && mt <= 207
        name = ("(n,xn)    ", "(n,xgma)  ", "(n,xp)    ", "(n,xd)    ",
                "(n,xt)    ", "(n,xhe3)  ", "(n,xa)    ")[mt - 200]
    elseif mt == 444
        name = "damage    "
    else
        i = mt
        i > 999 && (i -= 1000*(i ÷ 1000))
        i >= 600 && (i -= 399)
        if 1 <= i <= length(_APLOTS_MT_HNDF)
            name = _APLOTS_MT_HNDF[i]
        else
            error("aplots/mtname: MT=$mt (mapped index $i) not tabulated; add it from acecm.f90 before plotting this reaction")
        end
    end
    # charged-particle alternate names (acecm.f90:209-221)
    if izai > 1
        if mt == 4
            name = "(z,n)     "
        elseif mt == 5
            name = "(z,x)     "
        elseif (izai == 1001 && mt == 103) || (izai == 1002 && mt == 104) ||
               (izai == 1003 && mt == 105) || (izai == 2003 && mt == 106) ||
               (izai == 2004 && mt == 107)
            name = _APLOTS_MT_HNDF[4]  # "inelastic "
        end
    end
    return name
end
