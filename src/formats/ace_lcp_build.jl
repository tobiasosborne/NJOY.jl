# acelcp per-type block builder (flat XSS).
#
# Ref: njoy-reference/src/acefc.f90:9121-11056.
#
# Appends PTYPE/NTRO/PLOCT + per-type (HPD/MTRH/TYRH/LSIGH/SIGH/LANDH/ANDH/
# LDLWH/DLWH/YH) blocks to the XSS vector built by build_xss, mirroring the
# Fortran index arithmetic exactly. Operates on a mutable XSS wrapper so the
# `xss(i)=v` and `next` semantics match the Fortran scratch model.

"""Mutable XSS accessor mirroring Fortran's `xss(i)` (1-based) + `next`.
Auto-grows on write past the end (Fortran pre-allocates a huge array)."""
mutable struct _Xss
    v::Vector{Float64}
    is_int::Vector{Bool}
end

@inline function _ensure!(x::_Xss, i::Int)
    while length(x.v) < i
        push!(x.v, 0.0); push!(x.is_int, false)
    end
end
@inline _get(x::_Xss, i::Int) = x.v[i]
@inline function _set!(x::_Xss, i::Int, val::Real; isint::Bool=false)
    _ensure!(x, i)
    x.v[i] = Float64(val); x.is_int[i] = isint
end

"""
    append_particle_blocks!(xss_v, is_int, nxs, jxs, prod, nes, mf3_tab) -> nothing

Build the acelcp particle blocks and append them to `xss_v`/`is_int`, setting
NXS[NTYPE] and JXS[PTYPE/NTRO/PLOCT]. `jxs` holds the already-built basic-data
pointers (ESZ/MTR/LQR/TYR/LSIG/SIG/LAND/AND/END). `prod` is the production
table. `nes` is the energy-grid length.

This mutates `xss_v`, `is_int`, `nxs`, `jxs` in place; build_xss then sets
END/LEN2 afterward (acelcp updates `len2=next-1` and the basic-data `end` is
left at the pre-particle value — see acefc.f90:5893 vs 9192).
"""
function append_particle_blocks!(xss_v::Vector{Float64}, is_int::Vector{Bool},
                                  nxs::Vector{Int32}, jxs::Vector{Int32},
                                  prod::ParticleProduction, nes::Int,
                                  mf3_tab::Dict)
    # ---- pre-pass: graceful-skip unported-but-valid LAW forms -----------
    # The acelcp ANDH/DLWH builders only port a subset of the MF6 LAW forms
    # that Fortran acelcp handles (acefc.f90:9466-10836): MT=2 elastic
    # recoil, LAW=2 (two-body angular, lawnow=33), and LAW=4 (two-body
    # recoil, back-reads its LAW=2 partner). The other valid LAW forms
    # (LAW 1 Kalbach/tabulated, LAW 3, LAW 5 phase-space, LAW 6 n-body,
    # LAW 7 laboratory-angle) hit a loud `error(...)` in `_build_andh!`
    # (line ~504) / `_build_dlwh!` (line ~731), which would abort the whole
    # acer run (T51 p+H2: MT=28 LAW=6 → CRASH, was STRUCTURAL_FAIL). Mirror
    # the Phase 77 read_mf32 graceful-skip: if ANY contributing production
    # uses an unported LAW, @warn(maxlog=1) with MAT/MT/LAW context and bail
    # out BEFORE touching xss/nxs/jxs, leaving NTYPE=0 (no particle blocks).
    # The basic ACE tape is then exactly what a non-production evaluation
    # would produce — restoring T51's basic tape34 instead of crashing.
    # Genuine structural violations still error loudly downstream (Rule 6).
    # Proper fix = port LAW 1/3/5/6/7 in acelcp (bead NJOY_jl-cnh.2).
    if !_acelcp_laws_supported(prod)
        jxs[JXS_PTYPE] = 0; jxs[JXS_NTRO] = 0; jxs[JXS_PLOCT] = 0
        nxs[NXS_NTYPE] = 0
        return nothing
    end

    x = _Xss(copy(xss_v), copy(is_int))

    esz  = Int(jxs[JXS_ESZ])
    mtr  = Int(jxs[JXS_MTR])
    lqr  = Int(jxs[JXS_LQR])
    lsig = Int(jxs[JXS_LSIG])
    sig  = Int(jxs[JXS_SIG])
    land = Int(jxs[JXS_LAND])
    andb = Int(jxs[JXS_AND])
    endb = Int(jxs[JXS_END])   # = length(basic xss) (acefc.f90:5893)

    izai = prod.izai
    awi  = prod.awi
    awr  = prod.awr
    za   = prod.za
    emev = _LCP_EMEV
    delt = _LCP_DELT
    small = _LCP_SMALL
    emc2 = PhysicsConstants.amassn * PhysicsConstants.amu *
           PhysicsConstants.clight^2 / PhysicsConstants.ev / emev

    nprod = length(prod.mprod)

    # ---- store particle production types (acefc.f90:9191-9217) --------
    ptype = endb + 1
    ntype = 0
    pt_codes = Int[]
    if prod.t201 < _LCP_ETOP && izai != 1
        ntype += 1; _set!(x, ptype + ntype - 1, 1; isint=true); push!(pt_codes, 1)
    end
    if prod.t203 < _LCP_ETOP && izai != 1001
        ntype += 1; _set!(x, ptype + ntype - 1, 9; isint=true); push!(pt_codes, 9)
    end
    if prod.t204 < _LCP_ETOP && izai != 1002
        ntype += 1; _set!(x, ptype + ntype - 1, 31; isint=true); push!(pt_codes, 31)
    end
    if prod.t205 < _LCP_ETOP && izai != 1003
        ntype += 1; _set!(x, ptype + ntype - 1, 32; isint=true); push!(pt_codes, 32)
    end
    if prod.t206 < _LCP_ETOP && izai != 2003
        ntype += 1; _set!(x, ptype + ntype - 1, 33; isint=true); push!(pt_codes, 33)
    end
    if prod.t207 < _LCP_ETOP && izai != 2004
        ntype += 1; _set!(x, ptype + ntype - 1, 34; isint=true); push!(pt_codes, 34)
    end

    if ntype == 0
        jxs[JXS_PTYPE] = 0; jxs[JXS_NTRO] = 0; jxs[JXS_PLOCT] = 0
        nxs[NXS_NTYPE] = 0
        return nothing
    end

    # map a PTYPE code to the production-ZAP `ip`
    code_to_ip(ipt) = ipt == 1 ? 1 : ipt == 9 ? 1001 : ipt == 31 ? 1002 :
                      ipt == 32 ? 1003 : ipt == 33 ? 2003 : 2004

    # ---- count productions per type (acefc.f90:9228-9246) -------------
    ntro = ptype + ntype
    ploct = ntro + ntype
    for i in 1:ntype
        ip = code_to_ip(round(Int, _get(x, ptype + i - 1)))
        ntrh = 0
        for j in 1:nprod
            prod.iprod[j] == ip && prod.lprod[j] > 0 && (ntrh += 1)
        end
        _set!(x, ntro + i - 1, ntrh; isint=true)
    end
    next = ploct + 10 * ntype

    # ---- loop over each production type --------------------------------
    for itype in 1:ntype
        ipt  = round(Int, _get(x, ptype + itype - 1))
        ntrh = round(Int, _get(x, ntro + itype - 1))
        ip = code_to_ip(ipt)
        ipp = ip

        thresh = ip == 1 ? prod.t201 / emev : ip == 1001 ? prod.t203 / emev :
                 ip == 1002 ? prod.t204 / emev : ip == 1003 ? prod.t205 / emev :
                 ip == 2003 ? prod.t206 / emev : prod.t207 / emev

        # threshold grid index it (acefc.f90:9268-9271)
        it = 1
        while _get(x, esz + it - 1) < thresh * (1 - delt)
            it += 1
        end

        hpd = next
        ploct_base = ploct + 10 * (itype - 1)
        _set!(x, ploct_base, hpd; isint=true)
        _set!(x, hpd, it; isint=true)
        _set!(x, hpd + 1, nes - it + 1; isint=true)
        mtrh = hpd + 2 + 2 * (nes - it + 1)
        _set!(x, ploct_base + 1, mtrh; isint=true)
        tyrh = mtrh + ntrh
        _set!(x, ploct_base + 2, tyrh; isint=true)
        lsigh = tyrh + ntrh
        _set!(x, ploct_base + 3, lsigh; isint=true)
        sigh = lsigh + ntrh
        _set!(x, ploct_base + 4, sigh; isint=true)
        next = sigh
        # zero xs + heating columns
        for ie in it:nes
            _set!(x, hpd + 2 + ie - it, 0.0)
            _set!(x, hpd + 2 + (nes - it + 1) + ie - it, 0.0)
        end

        # ---- SIG/MTRH/TYRH/LSIGH/YIELD loop (acefc.f90:9292-9464) -----
        jp = 0
        for j in 1:nprod
            ipn = prod.iprod[j]; mt = prod.mprod[j]; mf = prod.kprod[j]
            (ipn == ip && mf != 0) || continue
            jp += 1
            _set!(x, mtrh + jp - 1, mt; isint=true)
            _set!(x, lsigh + jp - 1, next - sigh + 1; isint=true)
            # locate the basic-data reaction (mtr/lsig/sig) for this MT
            iaa, kk, _n = _find_reaction(x, mtr, lsig, sig, mt)
            mt == 2 && (iaa = 1)

            # MF=6 normal branch (LAW≠4-direct handled via the subsection scan)
            _build_yield_and_sig!(x, prod, j, mt, mf, ip, it, nes, esz,
                                   hpd, sigh, mtrh, tyrh, lsigh, jp,
                                   iaa, kk, () -> next, n -> (next = n))
        end

        # ---- LANDH/ANDH loop (acefc.f90:9466-10056) -------------------
        landh = next
        _set!(x, ploct_base + 5, landh; isint=true)
        andh = landh + ntrh
        _set!(x, ploct_base + 6, andh; isint=true)
        next = andh
        jp = 0
        for j in 1:nprod
            ipj = prod.iprod[j]; mt = prod.mprod[j]; mf = prod.kprod[j]
            (ipj == ip && mf != 0) || continue
            iaa, kk, _n = _find_reaction(x, mtr, lsig, sig, mt)
            q = _get(x, lqr + _react_ir(x, mtr, mt) - 1)
            jp += 1
            mt == 2 && (iaa = 1)
            next = _build_andh!(x, prod, j, mt, mf, ip, it, nes, esz,
                                hpd, land, andb, landh, andh, jp, iaa, kk, q,
                                awr, awi, emc2, next)
        end

        # ---- LDLWH/DLWH loop (acefc.f90:10058-10836) ------------------
        ldlwh = next
        _set!(x, ploct_base + 7, ldlwh; isint=true)
        dlwh = ldlwh + ntrh
        _set!(x, ploct_base + 8, dlwh; isint=true)
        next = dlwh
        jp = 0
        for j in 1:nprod
            ipj = prod.iprod[j]; mt = prod.mprod[j]; mf = prod.kprod[j]
            iskip = (ipj != ip || mf == 0)
            if !iskip
                jp += 1
                if mt == 2
                    _set!(x, ldlwh + jp - 1, 0; isint=true)
                    iskip = true
                end
            end
            iskip && continue
            iaa, kk, _n = _find_reaction(x, mtr, lsig, sig, mt)
            q = _get(x, lqr + _react_ir(x, mtr, mt) - 1)
            next = _build_dlwh!(x, prod, j, mt, mf, ip, it, nes, esz,
                                hpd, ldlwh, dlwh, landh, jp, iaa, kk, q,
                                awr, awi, emc2, next)
        end

        # ---- divide heating by total xs (acefc.f90:10838-10853) -------
        naa = round(Int, _get(x, hpd + 1))
        for ie in it:nes
            tot = _get(x, esz + nes + ie - 1)
            hv = _get(x, hpd + 2 + naa + ie - it)
            tot != 0.0 && (hv = hv / tot)
            hv < delt && (hv = 0.0)
            ip == 1 && (hv = 0.0)
            if izai > 1
                _set!(x, esz + 4*nes + ie - 1,
                      _get(x, esz + 4*nes + ie - 1) + hv)
            end
            _set!(x, hpd + 2 + naa + ie - it, round_sigfig(hv, 7, 0))
            _set!(x, esz + 4*nes + ie - 1,
                  round_sigfig(_get(x, esz + 4*nes + ie - 1), 7, 0))
        end

        # ---- YH block (acefc.f90:10855-10866) -------------------------
        yh = next
        _set!(x, ploct_base + 9, yh; isint=true)
        _set!(x, yh, 0; isint=true)
        next += 1
        for j in 1:nprod
            if prod.iprod[j] == ipp
                _set!(x, next, prod.mprod[j]; isint=true)
                next += 1
            end
        end
        _set!(x, yh, next - yh - 1; isint=true)
    end

    len2 = next - 1
    # write JXS/NXS for particle blocks
    jxs[JXS_PTYPE] = ptype
    jxs[JXS_NTRO]  = ntro
    jxs[JXS_PLOCT] = ploct
    nxs[NXS_NTYPE] = ntype

    # copy back the grown XSS
    resize!(xss_v, length(x.v)); copyto!(xss_v, x.v)
    resize!(is_int, length(x.is_int)); copyto!(is_int, x.is_int)
    nxs[NXS_LEN2] = len2
    nothing
end

# ---- helpers --------------------------------------------------------------

"""
    _acelcp_laws_supported(prod) -> Bool

Pre-pass guard for `append_particle_blocks!`. Return `true` iff every
contributing MF6 production in `prod` uses a LAW form this port can
faithfully render through `_build_andh!`/`_build_dlwh!`; otherwise emit a
loud `@warn(maxlog=1)` (with MAT/MT/LAW context) and return `false` so the
caller skips ALL particle blocks (NTYPE=0) without aborting the acer run.

Supported (matching the ported ANDH/DLWH branches):
  - MT=2          → elastic recoil (special-cased; no MF6 LAW lookup).
  - LAW=2         → two-body angular (ptleg2/pttab2) + law-33 energy.
  - LAW=4         → two-body recoil, *provided* its MT has a LAW=2 partner
                    subsection to back-read (`_find_law2_partner`).
Unported-but-valid (→ graceful skip): LAW 1/3/5/6/7 (Kalbach, phase-space,
n-body, laboratory-angle, …). Fortran acelcp handles these (acefc.f90:9990
LAW=7, :10258 LAW=1, :10541 LAW=6, …); the port does not yet — bead
NJOY_jl-cnh.2. We read each contributing MT's MF6 subsections once here
(same walk `_read_mf6_subs` does) purely to inspect LAWs — nothing is
mutated, so a `false` result leaves the tape untouched.
"""
function _acelcp_laws_supported(prod::ParticleProduction)
    nprod = length(prod.mprod)
    # cache subsections per MT so we read each MF6 section only once
    subs_cache = Dict{Int, Vector{_MF6Sub}}()
    get_subs(mt) = get!(subs_cache, mt) do
        _, s = _read_mf6_subs(prod.endf_path, prod.mat, mt)
        s
    end
    for j in 1:nprod
        prod.lprod[j] > 0 || continue
        mt = prod.mprod[j]
        mt == 2 && continue                     # elastic recoil: always OK
        subs = get_subs(mt)
        lprodj = prod.lprod[j]
        (1 <= lprodj <= length(subs)) || continue  # malformed index: let the
                                                    # real builder raise it
        law = subs[lprodj].law
        ok = if law == 2
            true
        elseif law == 4
            _find_law2_partner(subs) !== nothing   # needs a LAW=2 partner
        else
            false                                  # LAW 1/3/5/6/7: unported
        end
        if !ok
            @warn "acelcp: unported MF6 LAW=$law for MAT=$(prod.mat) MT=$mt \
                   (only MT=2 recoil + LAW=2/4 ported) — skipping particle \
                   production blocks (NTYPE=0)" maxlog=1
            return false
        end
    end
    true
end

"""Find the basic-data reaction index `ir` (1-based) for MT in the MTR block."""
function _react_ir(x::_Xss, mtr::Int, mt::Int)
    ir = 0; mtt = 0
    while mtt != mt
        ir += 1
        mtt = round(Int, _get(x, mtr + ir - 1))
        ir > 10000 && error("acelcp: MT=$mt not found in MTR block")
    end
    ir
end

"""Return (iaa, k, n) for MT: iaa = SIG block ie_start, k = sig-block base
index (pointing at the `ie` word), n = ne. Mirrors acefc.f90:9303-9308."""
function _find_reaction(x::_Xss, mtr::Int, lsig::Int, sig::Int, mt::Int)
    ir = _react_ir(x, mtr, mt)
    k = round(Int, _get(x, lsig + ir - 1)) + sig - 1
    n = round(Int, _get(x, k + 1))
    iaa = round(Int, _get(x, k))
    (iaa, k, n)
end

# =========================================================================
# SIG/yield builder (acefc.f90:9336-9460, MF=6 normal branch)
# =========================================================================
"""Accumulate the production cross section into HPD and store the yield in
the SIGH block. Reads the MF6/MT subsection ik=lprod(j) for the yield TAB1.
Sets TYRH = -lct (+1 if lct=1). Handles the LAW=4 case by using the paired
distribution's yield (the yield TAB1 is the same record regardless of LAW)."""
function _build_yield_and_sig!(x::_Xss, prod::ParticleProduction, j::Int,
                               mt::Int, mf::Int, ip::Int, it::Int, nes::Int,
                               esz::Int, hpd::Int, sigh::Int, mtrh::Int,
                               tyrh::Int, lsigh::Int, jp::Int,
                               iaa::Int, kk::Int, getnext, setnext)
    next = getnext()
    emev = _LCP_EMEV; delt = _LCP_DELT
    lprodj = prod.lprod[j]

    lct, subs = _read_mf6_subs(prod.endf_path, prod.mat, mt)
    # TYRH = -1 for lct=2, +1 for lct=1 (acefc.f90:9389-9394)
    _set!(x, tyrh + jp - 1, lct == 1 ? 1 : -1; isint=true)

    # nk==1 special: izap=za (mt=2). For T54 mt=2 lprod handled separately.
    sub = subs[lprodj]
    y = sub.yield   # the yield TAB1 (Tab1Record)

    # accumulate production xs (acefc.f90:9416-9428)
    for ie in iaa:nes
        e = _get(x, esz + ie - 1) * emev
        yv = _interp_tab1(y, e)
        yv < delt && (yv = 0.0)
        ss = mt == 2 ? _get(x, esz + 3*nes + ie - 1) : _get(x, 2 + kk + ie - iaa)
        tt = _get(x, hpd + 2 + ie - it) + yv * ss
        _set!(x, hpd + 2 + ie - it, round_sigfig(tt, 7, 0))
    end

    # store the yield (acefc.f90:9430-9454)
    _set!(x, next, 12; isint=true)
    _set!(x, next + 1, mt; isint=true)
    next += 2
    nrint = length(y.interp.nbt)
    int1  = isempty(y.interp.law) ? 2 : Int(y.interp.law[1])
    if nrint == 1 && int1 == 2
        _set!(x, next, 0; isint=true)
    else
        _set!(x, next, nrint; isint=true)
        for i in 1:nrint
            _set!(x, next + i, Int(y.interp.nbt[i]); isint=true)
            _set!(x, next + nrint + i, Int(y.interp.law[i]); isint=true)
        end
        next += 2 * nrint
    end
    next += 1
    ne = length(y.x)
    _set!(x, next, ne; isint=true)
    for i in 1:ne
        _set!(x, next + i,      round_sigfig(y.x[i] / emev, 7, 0))
        _set!(x, next + i + ne, round_sigfig(y.y[i], 7, 0))
    end
    next += 1 + 2 * ne

    setnext(next)
    nothing
end

"""Lin-lin/terp1 evaluate a Tab1Record at energy e (eV) honouring its
interpolation table (mirrors NJOY terpa). Used for the yield."""
function _interp_tab1(y::Tab1Record, e::Float64)
    np = length(y.x)
    np == 0 && return 0.0
    e <= y.x[1]   && return y.y[1]
    e >= y.x[np]  && return y.y[np]
    j = 1
    while j < np && y.x[j+1] < e
        j += 1
    end
    x1 = y.x[j]; y1 = y.y[j]; x2 = y.x[j+1]; y2 = y.y[j+1]
    # interpolation law for this panel
    law = _interp_law_at(y.interp, j)
    return _terp1(x1, y1, x2, y2, e, law)
end

function _interp_law_at(interp::InterpolationTable, point::Int)
    for r in 1:length(interp.nbt)
        if point < Int(interp.nbt[r])
            return Int(interp.law[r])
        end
    end
    isempty(interp.law) ? 2 : Int(interp.law[end])
end

"""terp1 (endf.f90:1604-1632) for INT=1..5."""
function _terp1(x1, y1, x2, y2, x, int)
    int == 1 && return y1
    (x2 == x1) && return y1
    if int == 2
        f = (x2 - x) / (x2 - x1)
        return y1 * f + y2 * (1 - f)
    elseif int == 3
        return y1 + (y2 - y1) * log(x / x1) / log(x2 / x1)
    elseif int == 4
        (y2 == y1) && return y1
        return y1 * exp((x - x1) / (x2 - x1) * log(y2 / y1))
    elseif int == 5
        (y2 == y1 || x == x1 || y1 == 0.0) && return y1
        return y1 * exp(log(x / x1) * log(y2 / y1) / log(x2 / x1))
    end
    # default lin-lin
    f = (x2 - x) / (x2 - x1)
    y1 * f + y2 * (1 - f)
end

# =========================================================================
# ANDH builder (acefc.f90:9466-10056)
# =========================================================================
"""Build one reaction's angular block in ANDH. Handles MT=2 elastic recoil
(reverse the basic elastic distribution), and LAW=2/LAW=4 (Legendre/tabulated
→ tabulated via ptleg2/pttab2). Returns updated `next`."""
function _build_andh!(x::_Xss, prod::ParticleProduction, j::Int, mt::Int,
                      mf::Int, ip::Int, it::Int, nes::Int, esz::Int, hpd::Int,
                      land::Int, andb::Int, landh::Int, andh::Int, jp::Int,
                      iaa::Int, kk::Int, q::Float64,
                      awr::Float64, awi::Float64, emc2::Float64, next::Int)
    small = _LCP_SMALL

    if mt == 2
        # --- elastic recoil: reverse the scattered-particle distribution
        #     (acefc.f90:9495-9574) ---
        _set!(x, landh + jp - 1, next - andh + 1; isint=true)
        loce = round(Int, _get(x, land)) + andb - 1
        ne = round(Int, _get(x, loce))
        _set!(x, next, ne; isint=true)
        leee = next
        for ie in 1:ne
            _set!(x, next + ie, _get(x, loce + ie))   # copy energies
        end
        nb = next + ne
        next = next + 1 + 2 * ne
        loce = loce + 1 + 2 * ne
        # heating scratch table
        scr = zeros(Float64, 8 + 2*ne + 16)
        scr[5] = 1; scr[6] = ne; scr[7] = ne; scr[8] = 5
        amass = awr / awi
        ubar = 0.0
        for ie in 1:ne
            intc = round(Int, _get(x, loce))
            n = round(Int, _get(x, loce + 1))
            # AND-block LC locator: integer (cf. _andh_law2! line ~675; Fortran
            # acefc.f90:9520 sets xss(nb+ie), written via write_integer_list /
            # typen iflag=1 at acefc.f90:13724 — i20 integer format, negative =
            # tabulated angular distribution).
            _set!(x, nb + ie, -(next - andh + 1); isint=true)
            _set!(x, next, intc; isint=true)
            _set!(x, next + 1, n; isint=true)
            for i in 1:n
                _set!(x, next + 1 + i, -_get(x, loce + 2 + n - i))
                v = _get(x, loce + 2 + 2*n - i)
                _set!(x, next + 1 + n + i, v)
                _get(x, next + 1 + n + i) < small && _set!(x, next + 1 + n + i, 0.0)
                if i == 1
                    _set!(x, next + 1 + 2*n + 1, 0.0)
                    ubar = 0.0
                end
                if i > 1 && intc == 1
                    s = _get(x, next + 1 + 2*n + i - 1) +
                        _get(x, next + 1 + n + i - 1) *
                        (_get(x, next + 1 + i) - _get(x, next + 1 + i - 1))
                    _set!(x, next + 1 + 2*n + i, round_sigfig(s, 7, 0))
                    ubar += _get(x, next + 1 + n + i - 1) *
                            (_get(x, next + 1 + i) - _get(x, next + 1 + i - 1)) +
                            (_get(x, next + 1 + i) + _get(x, next + 1 + i - 1)) / 2
                end
                if i > 1 && intc == 2
                    s = _get(x, next + 1 + 2*n + i - 1) +
                        (_get(x, next + 1 + n + i) + _get(x, next + 1 + n + i - 1)) *
                        (_get(x, next + 1 + i) - _get(x, next + 1 + i - 1)) / 2
                    _set!(x, next + 1 + 2*n + i, round_sigfig(s, 7, 0))
                    ubar += (_get(x, next + 1 + n + i) + _get(x, next + 1 + n + i - 1)) *
                            (_get(x, next + 1 + i) - _get(x, next + 1 + i - 1)) *
                            (_get(x, next + 1 + i) + _get(x, next + 1 + i - 1)) / 4
                end
            end
            renorm = 1 / _get(x, next + 1 + 3*n)
            for i in 1:n
                _set!(x, next + 1 + n + i, round_sigfig(renorm * _get(x, next + 1 + n + i), 7, 0))
                _set!(x, next + 1 + 2*n + i, round_sigfig(renorm * _get(x, next + 1 + 2*n + i), 9, 0))
            end
            next = next + 2 + 3*n
            loce = loce + 2 + 3*n
            e = _get(x, leee + ie)
            scr[8 + 2*ie - 1] = e
            scr[8 + 2*ie]     = 2 * amass * e * (1 + ubar) / (1 + amass)^2
        end
        # heating contribution (acefc.f90:9566-9574)
        naa = round(Int, _get(x, hpd + 1))
        for ie in it:nes
            e = _get(x, esz + ie - 1)
            h = _terpa_scr2(scr, ne, e)
            _set!(x, hpd + 2 + naa + ie - it,
                  _get(x, hpd + 2 + naa + ie - it) + h * _get(x, hpd + 2 + ie - it))
        end
        return next
    end

    # --- LAW=2 / LAW=4: re-read the MF6 subsection and build the table ---
    lct, subs = _read_mf6_subs(prod.endf_path, prod.mat, mt)
    lprodj = prod.lprod[j]
    sub = subs[lprodj]

    # default isotropic locator = -1 (acefc.f90:9837); set when LAW=2 found
    _set!(x, landh + jp - 1, -1; isint=true)

    izarec = 0; awprec = 0.0
    law = sub.law; awp = sub.awp
    if law == 4
        # back up to the paired LAW=2 subsection (acefc.f90:9844-9860).
        # The recoil's LAW=2 partner is the OTHER subsection of the same MT
        # whose IZAP differs; in T54 the LAW=2 sub is ik=1 (deuteron MT650) or
        # MT50 sub1 (neutron) — but for the recoil we need the SAME-MT LAW=2.
        izarec = sub.izap; awprec = awp
        partner = _find_law2_partner(subs)
        partner === nothing && error("acelcp ANDH: no LAW=2 partner for LAW=4 recoil MT=$mt")
        sub = partner
        law = sub.law; awp = sub.awp
    end

    if (law == 2)
        return _andh_law2!(x, sub, ip, izarec, awprec, awp, awr, awi, emc2,
                           q, it, nes, esz, hpd, andh, landh, jp, iaa, kk, next)
    else
        error("acelcp ANDH: unsupported LAW=$law for MT=$mt (only 2/4 + MT2 recoil ported)")
    end
end

"""Find the LAW=2 subsection among a MT's subsections (the recoil partner)."""
function _find_law2_partner(subs::Vector{_MF6Sub})
    for s in subs
        s.law == 2 && return s
    end
    nothing
end

"""terpa on a scratch heating table laid out as Fortran scr(llht): header at
scr[5]=NR=1, scr[6]=NP=ne, scr[7]=NBT, scr[8]=INT, then (E,V) pairs at
scr[8+2i-1],scr[8+2i]. Honours the INT field (lin-lin=2, log-log=5, …) via
terp1 — the MT=2 elastic recoil heating table uses INT=5 (acefc.f90:9515),
the LAW=2/4 tables use INT=2."""
function _terpa_scr2(scr::Vector{Float64}, ne::Int, e::Float64)
    int = round(Int, scr[8])
    Ei(i) = scr[8 + 2*i - 1]; Vi(i) = scr[8 + 2*i]
    e < Ei(1) && return 0.0   # Fortran terpa label 170 (endf.f90:1812-1817): below first point → 0; e == Ei(1) falls through to interp → Vi(1) (label 140)
    e >= Ei(ne) && return Vi(ne)
    j = 1
    while j < ne && Ei(j+1) < e
        j += 1
    end
    _terp1(Ei(j), Vi(j), Ei(j+1), Vi(j+1), e, int)
end

# =========================================================================
# LAW=2 angular reconstruction → ANDH (acefc.f90:9864-9987)
# =========================================================================
"""Build the LAW=2 angular block for one reaction into ANDH. `sub` is the
LAW=2 _MF6Sub; `izarec`>0 means this is a LAW=4 recoil (negate odd Legendre
coefficients). Uses ptleg2 (LANG=0 Legendre) or pttab2 (LANG>0 tabulated)."""
function _andh_law2!(x::_Xss, sub::_MF6Sub, ip::Int, izarec::Int,
                     awprec::Float64, awp::Float64, awr::Float64, awi::Float64,
                     emc2::Float64, q::Float64, it::Int, nes::Int, esz::Int,
                     hpd::Int, andh::Int, landh::Int, jp::Int,
                     iaa::Int, kk::Int, next::Int)
    small = _LCP_SMALL
    ne = sub.ne
    _set!(x, landh + jp - 1, next - andh + 1; isint=true)
    _set!(x, next, ne; isint=true)
    iebase = next
    ilbase = iebase + ne
    next = ilbase + ne + 1
    amass = awr / awi
    aprime = izarec == 0 ? awp / awi : awprec / awi

    # heating scratch table (scr(llht)): header + ne (E,V) pairs
    scr = zeros(Float64, 8 + 2*ne + 16)
    scr[5] = 1; scr[6] = ne; scr[7] = ne; scr[8] = 2

    for iie in 1:ne
        rec = sub.lists[iie]
        lang = Int(rec.L1)
        # Build the ptleg2/pttab2 scratch buffer `a` (1-based; a[1]=C1).
        nw = Int(rec.N1); nl = Int(rec.N2)
        # `a` must hold ptleg2's reconstructed output (8 + 2*ii words, ii up
        # to ptleg2's maxang); size generously.
        a = zeros(Float64, 8 + 2*7000 + 16)
        a[1] = rec.C1; a[2] = rec.C2; a[3] = Float64(rec.L1); a[4] = Float64(rec.L2)
        a[5] = Float64(rec.N1); a[6] = Float64(rec.N2)
        for d in 1:length(rec.data)
            a[6 + d] = rec.data[d]
        end
        if lang == 0
            if izarec != 0
                # negate odd-order Legendre coefficients (recoil cosine flip).
                # acefc.f90:9897-9901: scr(lld+5+iil) → a[6+iil].
                nl2 = Int(rec.N2)
                for iil in 1:nl2
                    isodd(iil) && (a[6 + iil] = -a[6 + iil])
                end
            end
            ptleg2!(a)
        else
            # LANG>0: already-tabulated (μ,p) data → normalize via pttab2.
            # acefc.f90:9904-9914. The LIST data is NL pairs (μ_k,p_k) at
            # rec.data[2k-1],rec.data[2k]; lay them out as pttab2 expects
            # (nr=1, np, pairs at a[7+2i],a[8+2i]) with interpolation iint.
            iint = lang - 10
            np = Int(rec.N2)
            a[5] = 1; a[6] = np; a[7] = np; a[8] = iint
            for i2 in 1:np
                a[7 + 2*i2] = rec.data[2*i2 - 1]
                a[8 + 2*i2] = rec.data[2*i2]
            end
            pttab2!(a)
        end
        _set!(x, iebase + iie, round_sigfig(a[2] / _LCP_EMEV, 7, 0))
        m = round(Int, a[5]); n = round(Int, a[6])
        loc = next - andh + 1
        _set!(x, ilbase + iie, -loc; isint=true)
        intc = round(Int, a[8])
        _set!(x, next, intc; isint=true)
        _set!(x, next + 1, n; isint=true)
        ubar = 0.0
        # cosine = scr(lld+4+2m+2i) → a[5+2m+2i]; pdf = scr(lld+5+2m+2i) →
        # a[6+2m+2i] (a[k]=scr[lld+k-1]). acefc.f90:9928-9931.
        for i in 1:n
            _set!(x, next + 1 + i,     round_sigfig(a[5 + 2*m + 2*i], 7, 0))
            _set!(x, next + 1 + n + i, round_sigfig(a[6 + 2*m + 2*i], 7, 0))
            _get(x, next + 1 + n + i) < small && _set!(x, next + 1 + n + i, 0.0)
            if i == 1
                _set!(x, next + 1 + 2*n + i, 0.0); ubar = 0.0
            end
            if i > 1 && intc == 1
                s = _get(x, next + 1 + 2*n + i - 1) +
                    _get(x, next + 1 + n + i - 1) *
                    (_get(x, next + 1 + i) - _get(x, next + 1 + i - 1))
                _set!(x, next + 1 + 2*n + i, round_sigfig(s, 7, 0))
                ubar += _get(x, next + 1 + n + i - 1) *
                        (_get(x, next + 1 + i) - _get(x, next + 1 + i - 1)) +
                        (_get(x, next + 1 + i) + _get(x, next + 1 + i - 1)) / 2
            end
            if i > 1 && intc == 2
                s = _get(x, next + 1 + 2*n + i - 1) +
                    (_get(x, next + 1 + n + i) + _get(x, next + 1 + n + i - 1)) *
                    (_get(x, next + 1 + i) - _get(x, next + 1 + i - 1)) / 2
                _set!(x, next + 1 + 2*n + i, round_sigfig(s, 7, 0))
                ubar += (_get(x, next + 1 + n + i) + _get(x, next + 1 + n + i - 1)) *
                        (_get(x, next + 1 + i) - _get(x, next + 1 + i - 1)) *
                        (_get(x, next + 1 + i) + _get(x, next + 1 + i - 1)) / 4
            end
        end
        renorm = 1 / _get(x, next + 1 + 3*n)
        for i in 1:n
            _set!(x, next + 1 + n + i,   round_sigfig(renorm * _get(x, next + 1 + n + i), 7, 0))
            _set!(x, next + 1 + 2*n + i, round_sigfig(renorm * _get(x, next + 1 + 2*n + i), 9, 0))
        end
        next = next + 2 + 3*n
        e = _get(x, iebase + iie)
        th = (1 + amass) * q / amass
        r1 = amass * (amass + 1 - aprime) / aprime
        r2 = aprime / (1 + amass)^2
        betasq = r1 * (1 + th / e)
        betasq < 0 && (betasq = 0.0)
        scr[8 + 2*iie - 1] = e
        scr[8 + 2*iie]     = e * r2 * (betasq + 1 + 2 * ubar * sqrt(betasq))
    end
    # heating contribution (acefc.f90:9977-9987)
    naa = round(Int, _get(x, hpd + 1))
    for ie in it:nes
        e = _get(x, esz + ie - 1)
        h = _terpa_scr2(scr, ne, e)
        ss = ie >= iaa ? _get(x, 2 + kk + ie - iaa) : 0.0
        _set!(x, hpd + 2 + naa + ie - it, _get(x, hpd + 2 + naa + ie - it) + h * ss)
    end
    next
end

# =========================================================================
# DLWH builder (acefc.f90:10058-10536)
# =========================================================================
"""Build one reaction's energy-distribution block in DLWH. Handles LAW=2
(two-body, lawnow=33, E1/E2 from the LAW=2 incident-energy span) and LAW=4
(two-body recoil, lawnow=33, E1=E2=0). Returns updated `next`.
MT=2 is skipped before this is called (ldlwh=0)."""
function _build_dlwh!(x::_Xss, prod::ParticleProduction, j::Int, mt::Int,
                      mf::Int, ip::Int, it::Int, nes::Int, esz::Int, hpd::Int,
                      ldlwh::Int, dlwh::Int, landh::Int, jp::Int, iaa::Int, kk::Int,
                      q::Float64, awr::Float64, awi::Float64, emc2::Float64,
                      next::Int)
    emev = _LCP_EMEV
    _set!(x, ldlwh + jp - 1, next - dlwh + 1; isint=true)
    last = next
    _set!(x, next, 0; isint=true)       # lnw
    _set!(x, next + 1, 0; isint=true)   # lawnow (set below)
    next += 3                            # skip lnw, lawnow, idat

    lct, subs = _read_mf6_subs(prod.endf_path, prod.mat, mt)
    lprodj = prod.lprod[j]
    sub = subs[lprodj]
    law = sub.law; awp = sub.awp

    if law == 2
        # lawnow=33, build TAB1 threshold + law-33 data (acefc.f90:10258-10475)
        lawnow = 33
        _set!(x, last + 1, lawnow; isint=true)
        lee = next
        _set!(x, next, 0; isint=true)       # nr
        _set!(x, next + 1, 2; isint=true)   # ne
        next = next + 2 + 2*2
        _set!(x, last + 2, next - dlwh + 1; isint=true)
        # E1/E2/P1/P2 from the LAW=2 incident-energy span (lines 10304-10310)
        ne = sub.ne
        e1 = sub.e_inc[1]; e2 = sub.e_inc[ne]
        _set!(x, lee + 2, round_sigfig(e1 / emev, 7, 0))
        _set!(x, lee + 4, 1.0)
        _set!(x, lee + 3, round_sigfig(e2 / emev, 7, 0))
        _set!(x, lee + 5, 1.0)
        # law-33 data
        amass = awr / awi; aprime = awp / awi
        _set!(x, next,     round_sigfig((1 + amass) * (-q) / amass, 7, 0))
        _set!(x, next + 1, round_sigfig(amass * (amass + 1 - aprime) / (1 + amass)^2, 7, 0))
        next += 2
        # no DLWH-stage heating for law=2 (handled in ANDH)
        return next

    elseif law == 4 || law == 3
        # lawnow=33 recoil (acefc.f90:10512-10536). E1=E2=0 (lee not filled).
        lawnow = 33
        _set!(x, last + 1, lawnow; isint=true)
        law == 3 && _set!(x, landh + jp - 1, 0; isint=true)
        _set!(x, next, 0; isint=true)
        _set!(x, next + 1, 2; isint=true)
        next = next + 2 + 2*2
        _set!(x, last + 2, next - dlwh + 1; isint=true)
        # law=3/4 uses the current subsection's own awp (no back-read here;
        # only ANDH backs up to the LAW=2 partner). acefc.f90:10519-10520.
        amass = awr / awi; aprime = awp / awi
        _set!(x, next,     round_sigfig((1 + amass) * (-q) / amass, 7, 0))
        _set!(x, next + 1, round_sigfig(amass * (amass + 1 - aprime) / (1 + amass)^2, 7, 0))
        if law != 4
            naa = round(Int, _get(x, hpd + 1))
            for ie in it:nes
                e = _get(x, esz + ie - 1)
                ss = ie >= iaa ? _get(x, 2 + kk + ie - iaa) : 0.0
                tt = _get(x, next + 1) * (e - _get(x, next)) * ss
                _set!(x, hpd + 2 + naa + ie - it, _get(x, hpd + 2 + naa + ie - it) + tt)
            end
        end
        next += 2
        return next
    else
        error("acelcp DLWH: unsupported LAW=$law for MT=$mt (only 2/3/4 ported)")
    end
end
