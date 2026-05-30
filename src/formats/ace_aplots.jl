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

"""`1p,e14.6` data field: C `%14.6E` reproduces it exactly (e.g. `  2.950000E-01`)."""
_e14_6(v::Real) = @sprintf("%14.6E", Float64(v))

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
        _aplots_unported("log-log heating plot not yet ported (heating non-trivial)")
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
        _aplots_unported("lin-lin heating plot not yet ported (heating non-trivial)")
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
    Int(nxs[NXS_NTYPE]) > 0 && _aplots_unported("aploxp (particle production) not yet ported")
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
            _aplots_unported("aplotr threshold-reaction plotting not yet ported")
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
            # perspective view of tabulated distribution
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
                pr("$(qu)angular distribution for $(name)$(qu)/")
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
    end
end

# =========================================================================
# mtname — MT -> reaction name (acecm.f90:21-224), with charged-particle
# alternate names.  Returns the 10-char Fortran name verbatim (trailing blanks
# preserved, as the Fortran write emits `a` = the full character(10)).
# =========================================================================
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
     "(n,2/2*1) ","(n,2/2*2) ","(n,2/2*3) ","(n,2/2*4) ","(n,n*0)   "]

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
