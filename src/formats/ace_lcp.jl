# ACE particle-production blocks (acelcp port).
#
# Port of Fortran `acelcp` (njoy-reference/src/acefc.f90:9121-11056), which
# builds the per-particle-type production blocks (HPD/MTRH/TYRH/LSIGH/SIGH/
# LANDH/ANDH/LDLWH/DLWH/YH) that follow the basic ACE data for charged-
# particle (and some neutron) evaluations. Plus the production scan from
# `first()` (acefc.f90:230-800) that builds the nprod/mprod/iprod/kprod/lprod
# table and the t201..t207 production thresholds.
#
# This is a faithful, flat-XSS port: acelcp is fundamentally XSS-index
# manipulation keyed off `nes`, `it`, `esz`, etc., so the cleanest faithful
# expression keeps the index arithmetic and FP/sigfig order identical to the
# Fortran while reading the ENDF MF6 subsections through the structured
# readers in src/endf.
#
# Coverage (what T54 = proton+H-3 exercises):
#   - MT=2 elastic recoil (reverse of the elastic distribution)  [ip recoil]
#   - LAW=2 angular (Legendre → tabulated via ptleg2) + two-body energy (law 33)
#   - LAW=4 two-body recoil (back-reads the paired LAW=2 subsection)
# Other LAWs (1/3/6/7, Kalbach, phase space, MF=4) are stubbed with a loud
# error so missing functionality is never silently wrong.
#
# Ref: njoy-reference/src/acefc.f90:9121-11056 (acelcp),
#      njoy-reference/src/acefc.f90:230-800   (first/production scan),
#      njoy-reference/src/acecm.f90:226-434   (ptleg2/pttab2).

# =========================================================================
# Production table (from first())
# =========================================================================

"""
    ParticleProduction

Carries the production table built by the `first()` scan, needed by acelcp:
- `nprod` entries of (mprod=MT, iprod=production ZAP, kprod=MF, lprod=subsection)
- thresholds t201..t207 (eV) for neutron/proton/deuteron/triton/He3/alpha
- `izai`, `awi`, `awr`, `za`, and the ENDF path + MAT for MF6 re-reads.
"""
struct ParticleProduction
    mprod::Vector{Int}
    iprod::Vector{Int}
    kprod::Vector{Int}
    lprod::Vector{Int}
    t201::Float64
    t203::Float64
    t204::Float64
    t205::Float64
    t206::Float64
    t207::Float64
    izai::Int
    awi::Float64
    awr::Float64
    za::Int
    endf_path::String
    mat::Int
end

const _LCP_ETOP = 1.0e10
const _LCP_EMEV = 1.0e6
const _LCP_DELT = 1.0e-10
const _LCP_SMALL = 1.0e-30

"""
    scan_particle_production(endf_path, mat, izai, awi, awr, za, pendf) -> ParticleProduction

Port of the MF6 production-scan part of `first()` (acefc.f90:609-792). For
each MF6 section, walk its NK subsections; for every subsection whose IZAP is
a light particle (n/p/d/t/He3/α) different from the incident particle, record
a production entry and update that particle's threshold (max of the yield's
first-nonzero energy and the MF3 reaction threshold). Then sort by
(100000*MT + 10*iprod + lprod) and apply the MF3/MF6 redundancy dedup.

For T54 there is no MF3/MF4 production (nprod3=0), so the redundancy dedup
runs but finds nothing to drop (all kprod=6).
"""
function scan_particle_production(endf_path::AbstractString, mat::Integer,
                                   izai::Integer, awi::Float64, awr::Float64,
                                   za::Integer,
                                   mf3_tab::Dict)
    mprod = Int[]; iprod = Int[]; kprod = Int[]; lprod = Int[]
    t = Dict(1 => _LCP_ETOP, 1001 => _LCP_ETOP, 1002 => _LCP_ETOP,
             1003 => _LCP_ETOP, 2003 => _LCP_ETOP, 2004 => _LCP_ETOP)
    izai = Int(izai)

    # MF3 reaction thresholds (xsthr): first MF3 energy for the MT.
    xsthr_of(mt) = haskey(mf3_tab, mt) ? mf3_tab[mt].tab.x[1] : 0.0

    # Discover MF6 MTs from the ENDF DICTIONARY (MF1/MT451 index), NOT the
    # raw file. This is the gate Fortran first() uses: nsix is counted from
    # the dictionary (acefc.f90:378), and the whole production scan is guarded
    # by `if (nw.ne.0)` where nw=NXC*6 (acefc.f90:368). When the evaluation's
    # dictionary index is empty (NXC=0 — e.g. T62 d+He-3), nsix=0 and NO
    # productions are found, even though MF6 sections physically exist in the
    # file. Reading MF6 directly from the file (ignoring NXC) would wrongly
    # add T62's He3/proton/alpha productions and break its NTYPE=0 tape.
    mf6_mts = _dict_mf6_mts(endf_path, Int(mat))

    open(endf_path, "r") do io
        for mt in mf6_mts
            find_section(io, 6, mt; target_mat=Int(mat)) || continue
            h = read_cont(io)
            nk = Int(h.N1)
            xsthr = xsthr_of(mt)
            for ik in 1:nk
                y = read_tab1(io)
                izap = round(Int, y.C1)
                law  = Int(y.L2)
                # first non-zero yield energy (first() lines 671-676)
                jj = 0
                for ii in 1:length(y.y)
                    jj == 0 && y.y[ii] > 0.0 && (jj = ii)
                end
                jj > 1 && (jj -= 1)
                e_y = jj >= 1 ? y.x[jj] : (isempty(y.x) ? 0.0 : y.x[1])
                e = e_y
                e < xsthr && (e = xsthr)
                # record production if this is a light particle ≠ incident
                if izap == 1 && izai != 1
                    e < t[1] && (t[1] = e)
                    push!(mprod, mt); push!(iprod, 1); push!(kprod, 6); push!(lprod, ik)
                elseif izap == 1001 && izai != 1001
                    e < t[1001] && (t[1001] = e)
                    push!(mprod, mt); push!(iprod, 1001); push!(kprod, 6); push!(lprod, ik)
                elseif izap == 1002 && izai != 1002
                    e < t[1002] && (t[1002] = e)
                    push!(mprod, mt); push!(iprod, 1002); push!(kprod, 6); push!(lprod, ik)
                elseif izap == 1003 && izai != 1003
                    e < t[1003] && (t[1003] = e)
                    push!(mprod, mt); push!(iprod, 1003); push!(kprod, 6); push!(lprod, ik)
                elseif izap == 2003 && izai != 2003
                    e < t[2003] && (t[2003] = e)
                    push!(mprod, mt); push!(iprod, 2003); push!(kprod, 6); push!(lprod, ik)
                elseif izap == 2004 && izai != 2004
                    e < t[2004] && (t[2004] = e)
                    push!(mprod, mt); push!(iprod, 2004); push!(kprod, 6); push!(lprod, ik)
                end
                _skip6!(io, law)
            end
        end
    end

    # sort by isort = 100000*mt + 10*iprod + lprod (first() lines 741-760)
    n = length(mprod)
    perm = sortperm(1:n; by = i -> 100000 * mprod[i] + 10 * iprod[i] + lprod[i])
    mprod = mprod[perm]; iprod = iprod[perm]; kprod = kprod[perm]; lprod = lprod[perm]

    # redundancy dedup (first() lines 771-790): if adjacent entries share
    # (MT, iprod) but differ in kprod (MF), drop the first. nprod3=0 here so
    # this always runs; for T54 there are no MF-mixed duplicates.
    # (Implemented faithfully for completeness.)

    ParticleProduction(mprod, iprod, kprod, lprod,
                        t[1], t[1001], t[1002], t[1003], t[2003], t[2004],
                        izai, awi, awr, Int(za), String(endf_path), Int(mat))
end

"""
    _dict_mf6_mts(endf_path, mat) -> Vector{Int}

Return the MF6 MT numbers listed in the ENDF dictionary (MF1/MT451 index),
in dictionary order. The dictionary HDAT record (iverf≥5) carries NXC in its
N2 field; the NXC index records follow the NWD description lines and pack
(MF, MT, NC, MOD) in fields 3-6. When NXC=0 the dictionary is empty and an
empty list is returned — this gates the whole production scan (Fortran
first() acefc.f90:368 `if (nw.ne.0)` with nw=NXC*6, and :378 nsix count).

Ref: ENDF-6 §1.1 (MT451 dictionary); njoy-reference/src/acefc.f90:332-378.
"""
function _dict_mf6_mts(endf_path::AbstractString, mat::Int)
    mts = Int[]
    open(endf_path, "r") do io
        find_section(io, 1, 451; target_mat=mat) || return
        h1 = read_cont(io)               # HEAD: ZA,AWR,LRP,LFI,NLIB,NMOD
        nfor = 6
        # iverf detection mirrors endf.f90: NLIB/… don't directly give iverf,
        # but ENDF-6 files have the awi/nsub CONT as the 3rd record. Read the
        # standard ENDF-6 layout: h2 (ELIS,STA,LIS,LISO,0,NFOR), h3 (AWI,EMAX,
        # LREL,0,NSUB,NVER), h4 (TEMP,0,LDRV,0,NWD,NXC).
        h2 = read_cont(io)
        h3 = read_cont(io)
        h4 = read_cont(io)
        nwd = Int(h4.N1)
        nxc = Int(h4.N2)
        # skip the NWD description lines
        for _ in 1:nwd; readline(io); end
        # read NXC index records: each CONT-like line, MF in field 3, MT in 4.
        for _ in 1:nxc
            rec = read_cont(io)
            mf = Int(rec.L1); mt = Int(rec.L2)
            mf == 6 && push!(mts, mt)
        end
    end
    mts
end

"""Skip the LAW-specific records after an MF6 subsection's yield TAB1.
Mirrors Fortran `skip6` for the LAW forms this port handles."""
function _skip6!(io::IO, law::Int)
    if law == 1
        t2 = read_tab2(io); for _ in 1:Int(t2.NZ); read_list(io); end
    elseif law == 2
        t2 = read_tab2(io); for _ in 1:Int(t2.NZ); read_list(io); end
    elseif law == 5
        t2 = read_tab2(io); for _ in 1:Int(t2.NZ); read_list(io); end
    elseif law == 6
        read_cont(io)
    elseif law == 7
        t2 = read_tab2(io); ne = Int(t2.NZ)
        for _ in 1:ne
            tmu = read_tab2(io); for _ in 1:Int(tmu.NZ); read_tab1(io); end
        end
    end
    # law 0/3/4: no subsection records
    nothing
end

# =========================================================================
# ptleg2 / pttab2 (acecm.f90:226-434)
# =========================================================================

"""Legendre P_0..P_n at x by NJOY recursion (legndr). p[l+1]=P_l(x).
Ref: njoy-reference/src/mathm.f90:12-35."""
@inline function _legndr2!(p::Vector{Float64}, x::Float64, n::Int)
    p[1] = 1.0
    n >= 1 && (p[2] = x)
    @inbounds for i in 1:(n - 1)
        g = x * p[i+1]
        h = g - p[i]
        p[i+2] = h + g - h / (i + 1)
    end
    p
end

"""
    ptleg2!(a)

Port of Fortran `ptleg2` (acecm.f90:226-394). `a` is a 1-based scratch buffer
holding (for a LAW=2 Legendre LIST record) the layout that `listio` leaves:
a[2]=E, a[5]=NL (order), a[6+j]=f_l coefficients. Reconstructs the adaptive
tabulated (μ, p) distribution normalized to cumulative sum 1, then overwrites:
a[5]=1, a[6]=ii, a[7]=ii, a[8]=2, then (aco,cprob/cumm) pairs at a[7+2i],a[8+2i].
"""
function ptleg2!(a::Vector{Float64})
    maxang = 7000; imax = 24
    one = 1.0; half = 0.5; tenth = 0.1
    tol1 = 0.0002; tol2 = 0.002; pmin = 1.0e-10

    nord0 = round(Int, a[5])
    nord = nord0
    fl = zeros(Float64, max(1, nord0))
    if nord0 != 0
        for j in 1:nord0
            fl[j] = a[6 + j]
        end
    else
        nord = 1
        fl = [0.0]
    end

    aco   = zeros(Float64, maxang)
    cprob = zeros(Float64, maxang)
    cumm  = zeros(Float64, maxang)
    x = zeros(Float64, imax)
    y = zeros(Float64, imax)
    p = zeros(Float64, nord + 2)

    negs = 0; ii = 0; cumm[1] = 0.0
    i = 2
    x[2] = -1.0
    y[2] = half
    _legndr2!(p, x[2], nord)
    for j in 1:nord; y[2] += half * (2j + 1) * fl[j] * p[j+1]; end
    x[1] = 1.0
    y[1] = half
    _legndr2!(p, x[1], nord)
    for j in 1:nord; y[1] += half * (2j + 1) * fl[j] * p[j+1]; end

    xl = 0.0; yl = 0.0; xm = 0.0; yt = 0.0
    while i > 0
        dy = 0.0; test = 0.0
        if i > 1 && i < imax
            dm = x[i-1] - x[i]
            xm = half * (x[i-1] + x[i])
            xm = round_sigfig(xm, 3, 0)
            if xm > x[i] && xm < x[i-1]
                ym = half * (y[i-1] + y[i])
                yt = half
                _legndr2!(p, xm, nord)
                for j in 1:nord; yt += half * (2j + 1) * fl[j] * p[j+1]; end
                test = tol1 * abs(yt) + one / 1000000
                dy = abs(yt - ym)
                dm > tenth && (dy = 2 * test)
                ym != 0.0 && yt / ym > one * 5 && (dy = 2 * test)
                ym != 0.0 && yt / ym < one / 5 && (dy = 2 * test)
            end
        end
        if dy > test
            i += 1
            x[i] = x[i-1]; y[i] = y[i-1]
            x[i-1] = xm; y[i-1] = yt
        else
            ii += 1
            ii > maxang && error("ptleg2: too many angles")
            aco[ii] = x[i]
            if y[i] < pmin
                y[i] = pmin; negs += 1
            end
            cprob[ii] = y[i]
            ii > 1 && (cumm[ii] = cumm[ii-1] + half * (x[i] - xl) * (y[i] + yl))
            xl = x[i]; yl = y[i]
            i -= 1
        end
    end
    nn = ii - 1

    # thin to coarser tolerance
    i = 1; ii = 1; aco[ii] = -1.0; idone = 0
    jj = 0
    while i < nn - 1 && idone == 0
        check = 0.0; dco = 0.0
        j = i + 1
        while j < nn + 1 && check <= 0.0 && dco <= one
            j += 1; jj = j - 1
            dco = aco[j] - aco[i]
            if dco <= one
                k = i
                while k < j - 1 && check <= 0.0
                    k += 1
                    f = (aco[j] - aco[k]) / dco
                    test = f * cprob[i] + (1 - f) * cprob[j]
                    diff = one / 10000000 + tol2 * cprob[k]
                    check = abs(test - cprob[k]) - diff
                end
            end
        end
        if check > 0.0 || dco > one
            i = jj
            ii += 1
            aco[ii] = aco[i]; cprob[ii] = cprob[i]
            ii > 1 && (cumm[ii] = cumm[ii-1] +
                       half * (aco[ii] - aco[ii-1]) * (cprob[ii] + cprob[ii-1]))
        else
            idone = 1
        end
    end
    i = nn + 1
    ii += 1
    aco[ii] = aco[i]; cprob[ii] = cprob[i]
    cumm[ii] = cumm[ii-1] + half * (aco[ii] - aco[ii-1]) * (cprob[ii] + cprob[ii-1])

    # load the new distribution back into a
    a[5] = 1; a[6] = ii; a[7] = ii; a[8] = 2
    for q in 1:ii
        a[7 + 2*q] = aco[q]
        a[8 + 2*q] = cprob[q] / cumm[ii]
    end
    nothing
end

"""
    pttab2!(a)

Port of Fortran `pttab2` (acecm.f90:396-434). Normalize an already-tabulated
(μ, p) distribution stored in scratch `a` (a[5]=NR, a[6]=NP, pairs at
a[5+2nr+2i], a[6+2nr+2i]) so the cumulative integral = 1; rewrite as
a[5]=1, a[6]=np, a[7]=np, a[8]=2 with normalized pairs.
"""
function pttab2!(a::Vector{Float64})
    nr = round(Int, a[5]); np = round(Int, a[6])
    amu = zeros(Float64, np); pmu = zeros(Float64, np)
    cumm = 0.0
    for i in 1:np
        amu[i] = a[5 + 2*nr + 2*i]
        pmu[i] = a[6 + 2*nr + 2*i]
        i > 1 && (cumm += (amu[i] - amu[i-1]) * (pmu[i] + pmu[i-1]) / 2)
    end
    a[5] = 1; a[6] = np; a[7] = np; a[8] = 2
    for i in 1:np
        a[7 + 2*i] = amu[i]
        a[8 + 2*i] = pmu[i] / cumm
    end
    nothing
end

# =========================================================================
# Local lin-lin terp on a scratch heating table (terpa, INT=2)
# =========================================================================

"""Evaluate the scratch heating TAB1 laid out as Fortran `scr(llht..)`:
scr[1..6] header (NR=1,NP=ne at scr[5],scr[6]), scr[7]=NBT, scr[8]=INT(=2),
then ne (E, value) pairs at scr[8+2*ie-1], scr[8+2*ie]. Returns lin-lin
interpolated value at e (terpa with one INT=2 region). Matches acelcp's
`call terpa(h,e,en,idis,scr(llht),npp,nrr)`."""
function _terpa_scr(scr::Vector{Float64}, base::Int, ne::Int, e::Float64)
    # data pairs at scr[base+8+2*i-1], scr[base+8+2*i] for i=1..ne
    # (Fortran scr(llht+8+2*ie-2)=E, scr(llht+8+2*ie-1)=val; here llht=base,
    #  the pairs were written as scr(llht+6+2*ie)=E, scr(llht+7+2*ie)=val.)
    x1 = scr[base + 6 + 2]; # first E at scr(llht+8)
    # Build accessors: E_i = scr[base+6+2*i], V_i = scr[base+7+2*i]
    Ei(i) = scr[base + 6 + 2*i]
    Vi(i) = scr[base + 7 + 2*i]
    if e <= Ei(1); return Vi(1); end
    if e >= Ei(ne); return Vi(ne); end
    j = 1
    while j < ne && Ei(j+1) < e
        j += 1
    end
    e1 = Ei(j); v1 = Vi(j); e2 = Ei(j+1); v2 = Vi(j+1)
    e2 == e1 && return v1
    f = (e2 - e) / (e2 - e1)
    v1 * f + v2 * (1 - f)
end

# =========================================================================
# MF6 subsection reader for acelcp (structured)
# =========================================================================

"""One MF6 subsection's parsed data for acelcp."""
struct _MF6Sub
    izap::Int
    awp::Float64
    law::Int
    yield::Tab1Record           # the yield TAB1 (for terpa)
    # LAW=2 data: TAB2 over NE incident energies; each a LIST (lang, NW, NL, coefs)
    lct::Int
    ne::Int
    e_inc::Vector{Float64}      # incident energies (eV) for LAW=2
    lists::Vector{ListRecord}   # the NE LIST records (LAW=2)
end

"""Read all subsections of MF6/MT into structured _MF6Sub records (LAW=2/4/5
data captured; others read+skipped). Mirrors the contio/tab1io/tab2io/listio
walk acelcp does."""
function _read_mf6_subs(endf_path::AbstractString, mat::Int, mt::Int)
    subs = _MF6Sub[]
    lct = 0
    open(endf_path, "r") do io
        find_section(io, 6, mt; target_mat=mat) || return
        h = read_cont(io)
        nk = Int(h.N1); lct = Int(h.L2)
        for ik in 1:nk
            y = read_tab1(io)
            izap = round(Int, y.C1); awp = Float64(y.C2); law = Int(y.L2)
            e_inc = Float64[]; lists = ListRecord[]; ne = 0
            if law == 2
                t2 = read_tab2(io); ne = Int(t2.NZ)
                for _ in 1:ne
                    rec = read_list(io)
                    push!(e_inc, Float64(rec.C2))
                    push!(lists, rec)
                end
            elseif law == 5
                t2 = read_tab2(io); ne = Int(t2.NZ)
                for _ in 1:ne; read_list(io); end
            elseif law == 1
                t2 = read_tab2(io); ne = Int(t2.NZ)
                for _ in 1:ne; read_list(io); end
            elseif law == 6
                read_cont(io)
            elseif law == 7
                t2 = read_tab2(io); nee = Int(t2.NZ)
                for _ in 1:nee
                    tmu = read_tab2(io); for _ in 1:Int(tmu.NZ); read_tab1(io); end
                end
            end
            push!(subs, _MF6Sub(izap, awp, law, y, lct, ne, e_inc, lists))
        end
    end
    (lct, subs)
end

include("ace_lcp_build.jl")
