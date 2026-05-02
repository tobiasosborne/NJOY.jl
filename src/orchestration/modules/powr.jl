# powr module runner — EPRI-CELL / EPRI-CPM library generator.
#
# Three independent driver paths selected by `lib`:
#   lib=1  →  fast()   — GAMTAP fast-neutron library  (PHASE B in progress)
#   lib=2  →  therm()  — LIBRAR thermal-neutron lib   (Phase C — TODO)
#   lib=3  →  cpm()    — CLIB EPRI-CPM library        (Phase D — TODO)
#
# Phase B current scope (this commit): lib=1 fast() ABSORPTION-ONLY path —
# enough to handle non-fissionable, no-MF6-matrix materials like carbon
# (iwa=1, iwf=0, iwr=0, xlol(1..3)=0). Standalone test against carbon
# 68-group GAM-I oracle is BIT-IDENTICAL.
#
# Phase B remaining items (raise loud TODOs in the relevant code paths):
#   - Multi-temperature loop (rtemp > 0 walks all temps — currently only
#     reference temp is honored).
#   - MF=6 matrix accumulation: P0/P1 elastic, inelastic, n2n (loce0, loce1,
#     locin, loc2n blocks).
#   - Fission MT=18..21+38 (locsf0, locchi, locnus, locdla blocks).
#   - Self-shielding f-factors via gamff (locabs, locsf — gated on
#     nsigz > 1 + iff = 1).
#   - matd<0 (iread=1) "absorption-only direct read" path.
#
# Output format per fast() line 528-547 (carbon path):
#   L1: (6i10)            nid iwa iwf iwr        (60 cols + pad to 80)
#   L2: (20a4)            description (16 chars + pad to 80)
#   L3: (36x,1p,3e12.5)   xlol(3) xla(3) xld(3)+1
#   L4: (1p,6e12.5)       xlol(1) xla(1) xld(1)+1 xlol(2) xla(2) xld(2)+1
#   L5+: (1p,6e12.5) ×N   absorption XS, ngnd=68 values (12 lines)
# All lines right-padded to 80 chars (Fortran formatted-record default).
#
# Ref: njoy-reference/src/powr.f90:362-651 (fast),
#      :781-1020 (gamll), :1022-1399 (gamxs), :752-779 (mergt).

"""
    powr_module(tapes::TapeManager, params::PowrParams)

Dispatch on `params.lib`. Phase B handles lib=1 (carbon-style absorption-
only path); other paths (multi-temp, MF6 matrices, fission, shielding,
lib=2/3) raise a loud TODO per CLAUDE.md Rule 6.
"""
function powr_module(tapes::TapeManager, params::PowrParams)
    @info "powr: ngendf=$(params.ngendf) nout=$(params.nout) " *
          "lib=$(params.lib) iprint=$(params.iprint) iclaps=$(params.iclaps)"

    if params.lib == 1
        return _powr_lib1_fast(tapes, params)
    end

    mode = params.lib == 2 ? "thermal (LIBRAR)" : "cpm (CLIB)"
    error("powr lib=$(params.lib) [$mode] not yet ported. See " *
          "src/orchestration/modules/powr.jl header for the multi-phase plan.")
end

# ============================================================================
# lib=1 / fast() — carbon-only absorption path
# ============================================================================

# Group structure constants (Fortran fast() lines 410-412).
const _POWR_FAST_NGND = 68     # output coarse-group count
const _POWR_FAST_NGNF = 185    # input fine-group count when iclaps=0

function _powr_lib1_fast(tapes::TapeManager, params::PowrParams)
    isempty(params.fast_mats) &&
        error("powr lib=1: no materials parsed (expected at least one card 3 " *
              "with matd ≠ 0)")

    gendf_path = resolve(tapes, params.ngendf)
    isfile(gendf_path) ||
        error("powr lib=1: input GENDF (unit $(params.ngendf)) not found at $gendf_path")

    out_path = resolve(tapes, params.nout)
    open(out_path, "w") do io
        for m in params.fast_mats
            _powr_fast_one_material(io, gendf_path, m, params)
        end
    end
    register!(tapes, params.nout, out_path)
    nothing
end

function _powr_fast_one_material(io::IO, gendf_path::String,
                                  m::PowrFastMaterial, params::PowrParams)
    if m.iread != 0
        error("powr lib=1: matd<0 (iread=1) absorption-direct-read mode not " *
              "yet ported. See src/orchestration/modules/powr.jl header.")
    end

    g = _powr_read_gendf_for_fast(gendf_path, m.matd, m.rtemp)

    has_fission = any(haskey(g.mf3, mt) for mt in (18, 19, 20, 21, 38))
    iwa = 1
    iwf = has_fission ? 1 : 0
    iwr = (g.nsigz > 1 && m.iff > 0) ? 1 : 0
    nfs = g.nfs

    # MF=6 supported: MT=2 elastic; MT∈[51,91] inelastic; MT∈{6-9,16,46-49} n2n;
    # MT=18 fission (only when iwf=1).
    function _powr_mf6_kind(mt::Int)
        mt == 2                                      && return :elastic
        mt == 18                                     && return :fission
        51 <= mt <= 91                               && return :inelastic
        mt == 16 || (6 <= mt <= 9) || (46 <= mt <= 49) && return :n2n
        return :unknown
    end
    mf6_unsupported = filter(mt -> _powr_mf6_kind(mt) == :unknown, collect(keys(g.mf6)))
    if iwr != 0 || nfs > 0 || !isempty(mf6_unsupported)
        error("powr lib=1: full Phase B path not yet ported (matd=$(m.matd) " *
              "has iwr=$iwr nfs=$nfs unsupported_mf6=$mf6_unsupported). " *
              "Currently supported: absorption + MF=6 elastic / inelastic / n2n / fission.")
    end

    nid  = g.iza * 10
    ngnd = _POWR_FAST_NGND

    # cflux: coarse-group P0 + P1 flux from MF=3/MT=1 (gamxs.f90:1147-1152).
    cflux   = zeros(Float64, ngnd)
    cfluxp1 = zeros(Float64, ngnd)
    if haskey(g.mf3, 1)
        for (ig, (flux, _xs)) in g.mf3[1]
            jg = ngnd - ig + 1
            1 <= jg <= ngnd && (cflux[jg] += flux)
        end
    end
    if haskey(g.mf3p1, 1)
        for (ig, (flux1, _xs)) in g.mf3p1[1]
            jg = ngnd - ig + 1
            1 <= jg <= ngnd && (cfluxp1[jg] += flux1)
        end
    end
    @inbounds for jg in 1:ngnd
        cflux[jg]   = cflux[jg]   != 0 ? 1 / cflux[jg]   : 0.0
        cfluxp1[jg] = cfluxp1[jg] != 0 ? 1 / cfluxp1[jg] : 0.0
    end

    # Unified storage `a` mirrors Fortran's contiguous layout. Pointers:
    #   locab0 = 1                     absorption (length ngnd)
    #   locsf0 = locab0 + ngnd         sigma_f    (length ngnd)
    #   locchi = locsf0 + ngnd         chi        (length ngnd)
    #   locnus = locchi + ngnd         nu*sigf → nu (length ngnd)
    # The fission-spectrum branch in Fortran gamxs (lines 1206-1208) writes
    # to `locchi + jg2c - 1` where jg2c can be 0 or -1, silently corrupting
    # the last 2 sigf cells. We replicate that exact behavior by using a
    # single unified array and Fortran-style pointer arithmetic.
    locab0 = 1
    locsf0 = locab0 + ngnd
    locchi = locsf0 + ngnd
    locnus = locchi + ngnd
    a_total_len = locnus + ngnd
    a = zeros(Float64, a_total_len)

    # Absorption: MT ∈ [102, 150] (gamxs.f90:1157-1163).
    for mt in keys(g.mf3)
        102 <= mt <= 150 || continue
        for (ig, (flux, xs)) in g.mf3[mt]
            jg = ngnd - ig + 1
            1 <= jg <= ngnd && (a[locab0 + jg - 1] += xs * flux * cflux[jg])
        end
    end

    # MF=3 fission (gamxs.f90:1164-1175) — MT=18 first; MT∈{19,20,21,38}
    # only if MT=18 was NOT seen (i318=0).
    cspc = zeros(Float64, ngnd)   # constant fission spectrum (cspc) — set up
                                  # by MF=6/MT=18 ig=0 record.
    cnorm = 0.0
    if iwf > 0
        i318 = haskey(g.mf3, 18)
        for mt in (18, 19, 20, 21, 38)
            haskey(g.mf3, mt) || continue
            (mt > 18 && i318) && continue
            for (ig, (flux, xs)) in g.mf3[mt]
                jg = ngnd - ig + 1
                1 <= jg <= ngnd || continue
                contrib = xs * flux * cflux[jg]
                a[locsf0 + jg - 1] += contrib
                a[locab0 + jg - 1] += contrib
            end
        end

        # MF=6/MT=18 (gamxs.f90:1182-1220): nu*sigf accumulation + chi.
        if haskey(g.mf6, 18)
            i618 = false
            for r in g.mf6[18]
                ig    = r.ig
                ng2   = r.ng2
                ig2lo = r.ig2lo
                nl    = r.nl
                nz    = r.nz
                data  = r.data

                if ig == 0
                    # Save constant fission spectrum (gamxs.f90:1216-1220):
                    # cspc(ig2lo+k-1) = data[1 + nl*nz*(k-1)] for k=1..ng2.
                    for k in 1:ng2
                        idx = ig2lo + k - 1
                        1 <= idx <= ngnd || continue
                        offset = nl * nz * (k - 1)
                        cspc[idx] = data[offset + 1]
                    end
                    continue
                end

                # MT=18: i618 detection (only one MT=18 record at ig=0; further
                # ig>0 records belong to MT=18 itself, not MT∈{19,20,21,38}).
                # We don't currently merge in MT=19+ on MF=6 anyway.
                jg = ngnd - ig + 1
                1 <= jg <= ngnd || continue

                for k in 2:ng2
                    offset = nl * nz * (k - 1)
                    flux1  = data[1]                # scr(l+lz)
                    nu_sf  = data[offset + 1]       # scr(loca)
                    if ig2lo != 0
                        # matrix part (gamxs.f90:1189-1200)
                        a[locnus + jg - 1] += nu_sf * flux1 * cflux[jg]
                        ig2 = ig2lo + k - 2
                        ig2c = ig2
                        if 1 <= ig2c <= ngnd
                            jg2c = ngnd - ig2c + 1
                            a[locchi + jg2c - 1] += flux1 * nu_sf
                            cnorm += flux1 * nu_sf
                        end
                    else
                        # spectrum part (gamxs.f90:1202-1213) with the
                        # bit-faithful jg2c-1 quirk: out-of-range writes
                        # spill into a[locsf0+ngnd-2..ngnd-1] (sigf last 2).
                        flux_loca_minus_1 = data[offset]   # scr(loca-1)
                        for i in 1:ngnd
                            ig2c = i
                            (ig2c == 0 || ig2c > ngnd) && continue
                            a[locnus + jg - 1] += cspc[i] * nu_sf *
                                                  flux_loca_minus_1 * cflux[jg]
                            jg2c = ngnd - ig2c - 1   # Fortran sign: -1, NOT +1.
                            idx = locchi + jg2c - 1
                            if 1 <= idx <= a_total_len
                                a[idx] += cspc[i] * nu_sf * flux_loca_minus_1
                            end
                            cnorm += nu_sf * flux_loca_minus_1
                        end
                    end
                end
            end
        end

        # Post-loop normalization: nu = nus / sigf (gamxs.f90:1357-1363).
        for ig in 1:ngnd
            sf = a[locsf0 + ig - 1]
            if sf != 0.0
                a[locnus + ig - 1] /= sf
            end
        end

        # Normalize chi: cnorm = 1/cnorm; chi *= cnorm (gamxs.f90:1366-1370).
        if cnorm != 0.0
            cnorm = 1 / cnorm
            for ig in 1:ngnd
                a[locchi + ig - 1] *= cnorm
            end
        end
    end

    absxs = view(a, locab0:locab0 + ngnd - 1)

    # MF=6 matrix accumulation. Each kind (elastic / inelastic / n2n) follows
    # the same packed-band-diagonal storage scheme; differences:
    #   elastic   uses Legendre P0 + P1 (and applies izref / 3× / cfluxp1).
    #   inelastic uses P0 only at first sigma-zero (loca = l+lz+nl*nz*(k-1)).
    #   n2n       same as inelastic but the accumulated value is divided by 2.
    elastic_recs   = collect(filter(r -> r.ig > 0, get(g.mf6, 2, NamedTuple[])))
    inelastic_recs = NamedTuple[]
    n2n_recs       = NamedTuple[]
    for (mt, recs) in g.mf6
        kind = _powr_mf6_kind(mt)
        if kind == :inelastic
            append!(inelastic_recs, recs)
        elseif kind == :n2n
            append!(n2n_recs, recs)
        end
        # Fission (kind == :fission) handled inline above; not packed here.
    end

    e0_block, e1_block, xla3, xld3, xlol3 =
        _powr_pack_matrix(elastic_recs, ngnd, m.izref, cflux, cfluxp1; kind=:elastic)
    in_block, _,        xla1, xld1, xlol1 =
        _powr_pack_matrix(inelastic_recs, ngnd, m.izref, cflux, cfluxp1; kind=:inelastic)
    n2n_block, _,       xla2, xld2, xlol2 =
        _powr_pack_matrix(n2n_recs, ngnd, m.izref, cflux, cfluxp1; kind=:n2n)

    xla = (xla1, xla2, xla3)
    xld = (xld1, xld2, xld3)
    xlol = (xlol1, xlol2, xlol3)

    # Fission-spectrum (kscr) block, written FIRST when iwf=1 or nfs>0
    # (fast.f90:497-524). Header line is `(i6,i2,10a4)` followed by chi
    # (`(1p,6e12.5)` × 12 lines for ngnd=68).
    if iwf > 0 || nfs > 0
        # `(i6,i2,10a4)` header: nid (i6) + 0 (i2) + fsn 40 chars (10×a4).
        chi_header = @sprintf("%6d%2d", nid, 0) * _powr_a4_pack(m.fsn, 10)
        println(io, _powr_pad80(chi_header))
        chi_view = view(a, locchi:locchi + ngnd - 1)
        _powr_write_e125_block(io, chi_view)
    end

    # Main (nscr) block.
    println(io, _powr_pad80(@sprintf("%10d%10d%10d%10d", nid, iwa, iwf, iwr)))
    println(io, _powr_pad80(_powr_word16(m.word)))
    println(io, _powr_pad80(repeat(" ", 36) *
                            _powr_e125(xlol[3]) * _powr_e125(xla[3]) *
                            _powr_e125(Float64(xld[3] + 1))))
    println(io, _powr_pad80(_powr_e125(xlol[1]) * _powr_e125(xla[1]) *
                            _powr_e125(Float64(xld[1] + 1)) *
                            _powr_e125(xlol[2]) * _powr_e125(xla[2]) *
                            _powr_e125(Float64(xld[2] + 1))))
    _powr_write_e125_block(io, absxs)
    if iwf > 0
        # sigf, then nu (gamxs.f90:550-553).
        _powr_write_e125_block(io, view(a, locsf0:locsf0 + ngnd - 1))
        _powr_write_e125_block(io, view(a, locnus:locnus + ngnd - 1))
    end
    # Matrix blocks (fast.f90:572-613): elastic (P0 then P1), inelastic, n2n.
    if xlol3 != 0
        _powr_write_e125_block(io, e0_block)
        _powr_write_e125_block(io, e1_block)
    end
    if xlol1 != 0
        _powr_write_e125_block(io, in_block)
    end
    if xlol2 != 0
        _powr_write_e125_block(io, n2n_block)
    end
    nothing
end

# Pack a string into `n` a4 (4-char) words = 4n chars. Right-pad short
# strings with blanks; truncate long ones.
_powr_a4_pack(s::AbstractString, n::Int) = rpad(s, 4n)[1:4n]

# Pack one MF=6 matrix kind into the gamll/gamxs band-diagonal storage.
# Returns (block, p1_block, xla, xld, xlol). p1_block is empty Float64[]
# for non-elastic kinds.
function _powr_pack_matrix(recs::Vector{NamedTuple}, ngnd::Int, izref::Int,
                            cflux::Vector{Float64}, cfluxp1::Vector{Float64};
                            kind::Symbol)
    isempty(recs) && return (Float64[], Float64[], 0.0, 0, 0.0)

    # gamll: la = min(igc); ld = max(igc - ig2lo) (both with the
    # ngn==ngnd ⇒ icgrp=identity simplification).
    la = ngnd
    ld = 0
    for r in recs
        igc = r.ig
        ig2loc = r.ig2lo
        (igc == 0 || igc > ngnd) && continue
        la = min(la, igc)
        ld = max(ld, igc - ig2loc)
    end
    (ld == 0 && la == ngnd) && return (Float64[], Float64[], 0.0, 0, 0.0)

    xla  = Float64(ngnd - la + 1)
    xld  = ld
    lx   = round(Int, xla) + xld - ngnd - 1
    lx   = max(lx, 0)
    xlol = xla * (xld + 1) - lx * (lx - 1) / 2
    n_cells = round(Int, xlol)
    block   = zeros(Float64, n_cells)
    p1_block = kind == :elastic ? zeros(Float64, n_cells) : Float64[]
    kdp = ngnd + 1 - ld

    # Accumulation factor: n2n divides the contribution by 2 (gamxs.f90:1348).
    n2n_factor = kind == :n2n ? 0.5 : 1.0

    for r in recs
        ig = r.ig
        jg = ngnd - ig + 1
        jg < 1 && continue
        jgp = jg > kdp ? kdp : jg
        lc  = (ld + 1) * (jgp - 1)
        if jg > kdp
            for jgt in 1:(jg - kdp)
                lc += ld + 2 - jgt
            end
        end
        nl    = r.nl
        nz    = r.nz
        ng2   = r.ng2
        ig2lo = r.ig2lo
        data  = r.data
        for k in 2:ng2
            ig2  = ig2lo + k - 2
            (ig2 == 0 || ig2 > ngnd) && continue
            jg2c = max(ngnd - ig2 + 1, jg)
            loc  = lc + (jg2c - jg) + 1
            1 <= loc <= n_cells || continue
            if kind == :elastic
                # P0 uses sigma-zero index izref (gamxs.f90:1311).
                offset = nl * ((izref - 1) + nz * (k - 1))
                block[loc] += data[offset + 1] * data[1] * cflux[jg]
                if nl > 1
                    flx = cfluxp1[jg] != 0 ? cfluxp1[jg] : cflux[jg]
                    p1_block[loc] += 3 * data[offset + 2] * data[2] * flx
                end
            else
                # Inelastic + n2n use the first sigma-zero (gamxs.f90:1284, :1346).
                offset = nl * nz * (k - 1)
                block[loc] += data[offset + 1] * data[1] * cflux[jg] * n2n_factor
            end
        end
    end

    (block, p1_block, xla, xld, xlol)
end

# ============================================================================
# Minimal GENDF reader for powr's needs
# ============================================================================

# Walks a GENDF file and captures everything powr's gamll/gamxs need at the
# reference temperature: MF=3 P0 flux + xs per (mt, ig); MF=3 P1 flux + xs
# (Legendre order 2 from MT=1 only); MF=6 transfer matrix records per mt;
# plus iza, nsigz, nfs.
function _powr_read_gendf_for_fast(gendf_path::String, matd::Int, rtemp::Float64)
    mf3   = Dict{Int, Dict{Int, Tuple{Float64, Float64}}}()  # P0
    mf3p1 = Dict{Int, Dict{Int, Tuple{Float64, Float64}}}()  # P1 (when nl≥2)
    mf6   = Dict{Int, Vector{NamedTuple}}()
    iza = 0; nsigz = 1; nfs = 0

    lines = readlines(gendf_path)
    nlines = length(lines)
    eps = 1e-4

    function _line_meta(li::Int)
        p = rpad(lines[li], 80)
        (mat = _parse_int(p[67:70]),
         mf  = _parse_int(p[71:72]),
         mt  = _parse_int(p[73:75]),
         seq = _parse_int(p[76:80]),
         line = p)
    end

    # Pull iza + nsigz from the matd's MF=1/MT=451.
    idx = 1
    while idx <= nlines
        length(lines[idx]) < 75 && (idx += 1; continue)
        m = _line_meta(idx)
        if m.mat == matd && m.mf == 1 && m.mt == 451 && m.seq == 1
            iza = round(Int, parse_endf_float(m.line[1:11]))
            if idx + 1 <= nlines
                p2 = rpad(lines[idx + 1], 80)
                nsigz = _parse_int(p2[34:44])
            end
            break
        end
        idx += 1
    end

    # Helper: read one LIST record's data section into a Vector{Float64} of
    # length nw (6 values per line, 11-char fields).
    function _read_list_data(start_idx::Integer, nw::Integer)
        nw = Int(nw)
        data = Vector{Float64}(undef, nw)
        nlines_data = cld(nw, 6)
        k = 0
        for li in 0:(nlines_data - 1)
            p = rpad(lines[start_idx + li], 80)
            for col in 0:5
                k >= nw && break
                k += 1
                data[k] = parse_endf_float(p[1 + 11*col : 11 + 11*col])
            end
        end
        data
    end

    idx = 1
    while idx <= nlines
        length(lines[idx]) < 75 && (idx += 1; continue)
        m = _line_meta(idx)
        m.mat != matd && (idx += 1; continue)

        if m.mf == 3 && m.mt > 0 && m.seq == 1
            mt = m.mt
            # HEAD record: NL on column 23-33, NZ on 34-44.
            nl_section = _parse_int(m.line[23:33])
            mt_data   = get!(mf3,   mt, Dict{Int, Tuple{Float64, Float64}}())
            mt_data_p1 = nl_section >= 2 ?
                get!(mf3p1, mt, Dict{Int, Tuple{Float64, Float64}}()) : nothing
            idx += 1
            while idx <= nlines
                length(lines[idx]) < 75 && (idx += 1; continue)
                p2 = rpad(lines[idx], 80)
                m2 = (mat=_parse_int(p2[67:70]), mf=_parse_int(p2[71:72]),
                      mt=_parse_int(p2[73:75]))
                (m2.mat != matd || m2.mf != 3 || m2.mt != mt) && break
                # LIST CONT: TEMP, 0, NG2, IG2LO, NW, IG
                temp = parse_endf_float(p2[1:11])
                nw   = _parse_int(p2[45:55])
                ig   = _parse_int(p2[56:66])
                ng2  = _parse_int(p2[23:33])
                nz   = _parse_int(p2[34:44])
                if abs(temp - rtemp) <= eps && ig > 0 && nw >= 2 && idx + 1 <= nlines
                    data = _read_list_data(idx + 1, nw)
                    # MF=3 layout (NJOY GENDF): for k=1 (only k since ng2=2),
                    # values 1..nl=flux per Legendre order, then xs at
                    # positions nl+1..nl+nl*nz (P0 first sigma-zero etc).
                    nl = nl_section
                    flux_p0 = data[1]
                    xs_p0   = nl + 1 <= length(data) ? data[nl + 1] : 0.0
                    mt_data[ig] = (flux_p0, xs_p0)
                    if nl >= 2
                        flux_p1 = data[2]
                        # P1 xs is the 2nd Legendre coefficient for k=2 (group→g):
                        # position nl + 2 (l=2 of first sigma-zero of k=2 entry).
                        xs_p1 = nl + 2 <= length(data) ? data[nl + 2] : 0.0
                        mt_data_p1[ig] = (flux_p1, xs_p1)
                    end
                end
                idx += 1 + cld(nw, 6)
            end

        elseif m.mf == 6 && m.mt > 0 && m.seq == 1
            mt = m.mt
            # HEAD record: NL on cols 23-33, NZ on cols 34-44.
            nl_section = _parse_int(m.line[23:33])
            nz_section = _parse_int(m.line[34:44])
            recs = get!(mf6, mt, NamedTuple[])
            idx += 1
            while idx <= nlines
                length(lines[idx]) < 75 && (idx += 1; continue)
                p2 = rpad(lines[idx], 80)
                m2 = (mat=_parse_int(p2[67:70]), mf=_parse_int(p2[71:72]),
                      mt=_parse_int(p2[73:75]))
                (m2.mat != matd || m2.mf != 6 || m2.mt != mt) && break
                temp  = parse_endf_float(p2[1:11])
                ng2   = _parse_int(p2[23:33])
                ig2lo = _parse_int(p2[34:44])
                nw    = _parse_int(p2[45:55])
                ig    = _parse_int(p2[56:66])
                # MF=6 needs ig=0 records too (constant fission spectrum
                # cspc setup) — only filter ig < 0.
                if abs(temp - rtemp) <= eps && ig >= 0 && nw > 0 && idx + 1 <= nlines
                    data = _read_list_data(idx + 1, nw)
                    push!(recs, (ig=ig, ig2lo=ig2lo, ng2=ng2,
                                 nl=nl_section, nz=nz_section, data=data))
                end
                idx += 1 + cld(nw, 6)
            end

        elseif m.mf == 5 && m.mt == 455 && m.seq == 1
            if idx + 1 <= nlines
                p2 = rpad(lines[idx + 1], 80)
                nfs += _parse_int(p2[23:33])
            end
            idx += 1
        else
            idx += 1
        end
    end

    (mf3=mf3, mf3p1=mf3p1, mf6=mf6, iza=iza, nsigz=nsigz, nfs=nfs)
end

# ============================================================================
# Formatted-write helpers — Fortran (1p,e12.5), (20a4), 80-col padding
# ============================================================================

# (1p,e12.5): 12-char field with 1 leading digit, 5 fractional digits, E±NN.
# Fortran formats 0.0 as " 0.00000E+00". Negative/positive numbers get a
# leading sign or space.
function _powr_e125(x::Real)
    x = Float64(x)
    if x == 0.0
        return " 0.00000E+00"
    end
    s = @sprintf("%13.5E", x)   # produces "+1.67033E-01" or " 1.67033E-01"
    # %E uses 1 digit before the decimal automatically — matches 1pE12.5.
    # Width 13 includes a leading sign char; Fortran 1pE12.5 omits the
    # explicit '+' so positives get a leading space. Strip a leading '+'.
    s = replace(s, "+1." => " 1.", count = 1)
    s = replace(s, "+2." => " 2.", count = 1)   # %E always uses + or -
    # Generic fix: if the leading char is '+', replace with space.
    if startswith(s, "+")
        s = " " * s[2:end]
    end
    # Trim leading space if width is 13 (we want 12 chars).
    length(s) == 13 && (s = s[2:end])
    length(s) == 12 || error("powr: e12.5 format produced $(length(s)) chars: $s")
    s
end

# (20a4): hollerith description, 16 chars (4 a4 words) per Fortran fast().
# `word` is read as `(4a4)` from a 12-char hid (`hid=' '` then read).
# Effectively, the description is stored as up to 16 chars left-justified
# with trailing blanks. The output via (20a4) writes 20×4=80 chars.
_powr_word16(word::AbstractString) = rpad(word, 16)

# Right-pad to 80 chars (Fortran formatted-record default record length).
_powr_pad80(s::AbstractString) = rpad(s, 80)

# Write a vector as (1p,6e12.5) records — 6 values per line, padded to 80.
function _powr_write_e125_block(io::IO, vals::AbstractVector{<:Real})
    n = length(vals)
    i = 1
    while i <= n
        buf = ""
        for _ in 1:6
            i > n && break
            buf *= _powr_e125(vals[i])
            i += 1
        end
        println(io, _powr_pad80(buf))
    end
end
