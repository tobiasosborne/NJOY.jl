# powr module runner — EPRI-CELL / EPRI-CPM library generator.
#
# Three independent driver paths selected by `lib`:
#   lib=1  →  fast()   — GAMTAP fast-neutron library  (PHASES 60-65 LANDED)
#   lib=2  →  therm()  — LIBRAR thermal-neutron lib   (Phase C — TODO)
#   lib=3  →  cpm()    — CLIB EPRI-CPM library        (Phase D — TODO)
#
# Phase 60-65 scope: lib=1 fast() handles the full single-material-per-call
# path EXCEPT delayed-neutron spectra and the matd<0 iread=1 corner.
# Bit-identical against four oracles:
#   - Carbon abs-only (Phase 60).
#   - Carbon abs + MF=6 elastic + inelastic + n2n matrices (Phases 61-62).
#   - U-235 fission + MF=6/MT=18 chi/nu (Phase 63).
#   - Fe-56 multi-T multi-σ₀ self-shielding f-factors (Phases 64-65).
#
# Phase B remaining items (raise loud TODOs in the relevant code paths):
#   - Delayed-neutron spectra (`nfs > 0` with MF=5/MT=455 → locdla block).
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

    g = _powr_read_gendf_for_fast(gendf_path, m.matd, m.rtemp; nsgz_user=m.nsgz)

    # Reference temperature view: per-MT dict with the σ₀=jz_ref XS extracted.
    # gamxs uses jz=1 (= σ₀=∞) for absorption / fission accumulation; the
    # elastic matrix uses jz=izref via _powr_pack_matrix's own indexing.
    rtemp_idx = g.rtemp_idx
    @assert rtemp_idx > 0 "internal: reader returned rtemp_idx=0 with rtemp>0"
    mf3_ref   = g.mf3[rtemp_idx]
    mf3p1_ref = g.mf3p1[rtemp_idx]

    # Helper that returns a (flux, xs_at_jz1) tuple for a (mt, ig) lookup
    # at the reference temperature — preserves the pre-Phase-64 call shape.
    _xs_ref(d, mt, ig) = (d[mt][ig][1], d[mt][ig][2][1])

    has_fission = any(haskey(mf3_ref, mt) for mt in (18, 19, 20, 21, 38))
    iwa = 1
    iwf = has_fission ? 1 : 0
    iwr_initial = (g.nsigz > 1 && m.iff > 0) ? 1 : 0   # gamll iwr (binary).
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
    if nfs > 0 || !isempty(mf6_unsupported)
        error("powr lib=1: not yet ported for matd=$(m.matd) " *
              "(nfs=$nfs unsupported_mf6=$mf6_unsupported). " *
              "Currently supported: absorption + MF=6 elastic / inelastic / n2n / " *
              "fission + multi-T multi-σ₀ self-shielding (no delayed-neutron spectra).")
    end

    nid  = g.iza * 10
    ngnd = _POWR_FAST_NGND

    # cflux: coarse-group P0 + P1 flux from MF=3/MT=1 (gamxs.f90:1147-1152).
    cflux   = zeros(Float64, ngnd)
    cfluxp1 = zeros(Float64, ngnd)
    if haskey(mf3_ref, 1)
        for (ig, fxs) in mf3_ref[1]
            flux = fxs[1]
            jg = ngnd - ig + 1
            1 <= jg <= ngnd && (cflux[jg] += flux)
        end
    end
    if haskey(mf3p1_ref, 1)
        for (ig, fxs) in mf3p1_ref[1]
            flux1 = fxs[1]
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
    # Reference temperature, σ₀=jz=1 (=∞).
    for mt in keys(mf3_ref)
        102 <= mt <= 150 || continue
        for (ig, fxs) in mf3_ref[mt]
            flux, xs_vec = fxs
            xs = xs_vec[1]
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
        i318 = haskey(mf3_ref, 18)
        for mt in (18, 19, 20, 21, 38)
            haskey(mf3_ref, mt) || continue
            (mt > 18 && i318) && continue
            for (ig, fxs) in mf3_ref[mt]
                flux, xs_vec = fxs
                xs = xs_vec[1]
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

    # Self-shielding f-factors (gamff). Computed AFTER the reference-temp
    # absorption / fission accumulators are filled, since the final-pass
    # log-formula normalizes against the reference σ₀=∞ at-temp values.
    # Returns (abs_factors, sf_factors, iwr_final, jwf). When iwr_initial==0
    # this is a no-op and iwr_final = 0.
    abs_ff, sf_ff, iwr, jwf_ff = _powr_pack_ffactors(g, m, cflux, ngnd, locab0,
                                                      locsf0, a, iwr_initial)
    # nff is hardcoded to 1 in fast() at powr.f90:414. The writer block fires
    # iff nff!=0 AND iwr!=0; with nff=1 the gate reduces to iwr!=0.
    nff = 1

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

    # Main (nscr) block. iwr is the gamff-computed value (count of self-
    # shielded groups = ngmax - iglo + 1), NOT the gamll 0/1 toggle.
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
    # Self-shielding f-factor block (powr.f90:614-641) — gated on nff!=0
    # AND iwr!=0. With nff hardcoded to 1, the gate reduces to iwr!=0.
    if nff != 0 && iwr != 0
        _powr_write_ffactor_block(io, abs_ff, sf_ff, jwf_ff,
                                   g.sigz, g.tmpr, g.nsigz, length(g.tmpr), iwr)
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
# gamff port — multi-T multi-σ₀ self-shielding f-factor accumulation
# ============================================================================

# Default group range from powr.f90:411-412 — hardcoded for GAM-I 68 group.
const _POWR_FAST_NGMIN = 19
const _POWR_FAST_NGMAX = 62

# Mirror of powr.f90:1401-1542 (gamff). Returns
#   (abs_factors, sf_factors, iwr_final, jwf)
# where abs_factors / sf_factors are flat Vector{Float64} storing the log-
# transformed f-factor cells in the SAME column-major [jz, itp, jg] order as
# Fortran's `a(locabs..)` / `a(locsf..)`. iwr_final = ngmax - iglo + 1 is the
# count of self-shielded groups (iglo = first group where |1 - f| > toler).
# jwf is 1 if any fission MT contributed to sf_factors, 0 otherwise.
#
# When iwr_initial == 0 (i.e. nsigz==1 OR iff==0 OR rtemp<=0) gamff is a
# no-op and we return zero-length arrays + iwr_final=0.
#
# `cflux_ref` is the reference-temperature inverted coarse flux from gamxs.
# At itp=1 the Fortran initializes mf3mt1=1 so the first MT≥2 record uses
# the gamxs-inverted cflux as-is. The first MT=1 record (when igc<=1) zeros
# cflux and sets mf3mt1=0 which forces the next MT≥2 to invert. We replicate
# this by precomputing per-temp inverted cfluxes from each temp's MT=1 (when
# present) and falling back to `cflux_ref` for itp=1 when no MT=1 is present.
function _powr_pack_ffactors(g, m::PowrFastMaterial, cflux_ref::Vector{Float64},
                              ngnd::Int, locab0::Int, locsf0::Int,
                              a::Vector{Float64}, iwr_initial::Int)
    if iwr_initial == 0
        return (Float64[], Float64[], 0, 0)
    end

    ngmin = _POWR_FAST_NGMIN
    ngmax = _POWR_FAST_NGMAX
    ntp   = length(g.tmpr)
    nsigz = g.nsigz

    # Flat 3D accumulators in column-major [jz, itp, jg] order to match
    # Fortran's a(locabs..) layout. Linear index for (jz, itp, jg):
    #   off = (jg-1)*nsigz*ntp + (itp-1)*nsigz + (jz-1) + 1
    abs_acc = zeros(Float64, nsigz * ntp * ngnd)
    sf_acc  = zeros(Float64, nsigz * ntp * ngnd)
    @inline _idx(jz, itp, jg) = (jg-1)*nsigz*ntp + (itp-1)*nsigz + jz
    jwf = 0

    for t_idx in 1:ntp
        # Per-temp cflux: rebuild from this temp's MF=3/MT=1 if present.
        # gamff's mf3mt1 sentinel logic resolves to: at itp=1 use ref cflux
        # if no MT=1 exists; at itp>1 require MT=1 to recompute. We just
        # always rebuild if MT=1 present, else fall back to ref cflux.
        cflux_t = if haskey(g.mf3[t_idx], 1)
            ct = zeros(Float64, ngnd)
            for (ig, fxs) in g.mf3[t_idx][1]
                flux = fxs[1]
                jg = ngnd - ig + 1
                1 <= jg <= ngnd && (ct[jg] += flux)
            end
            for jg in 1:ngnd
                ct[jg] = ct[jg] != 0.0 ? 1.0 / ct[jg] : 0.0
            end
            ct
        else
            cflux_ref
        end

        i318_t = haskey(g.mf3[t_idx], 18)
        for (mt, igdict) in g.mf3[t_idx]
            if 102 <= mt <= 107
                # absorption (gamff lines 1452: mth.gt.101.and.mth.lt.108).
                for (ig, fxs) in igdict
                    flux, xs_vec = fxs
                    jg = ngnd - ig + 1
                    1 <= jg <= ngnd || continue
                    lim = min(length(xs_vec), nsigz)
                    base = _idx(0, t_idx, jg)
                    for jz in 1:lim
                        abs_acc[base + jz] += xs_vec[jz] * flux * cflux_t[jg]
                    end
                end
            elseif mt == 18 || (19 <= mt <= 21) || mt == 38
                # fission (gamff lines 1453: 17<mth<22 or mth==38).
                # i318 latches MT=18; subsequent MT∈{19,20,21,38} skip if
                # MT=18 was seen (gamff lines 1489-1490).
                (mt > 18 && i318_t) && continue
                jwf = 1
                for (ig, fxs) in igdict
                    flux, xs_vec = fxs
                    jg = ngnd - ig + 1
                    1 <= jg <= ngnd || continue
                    lim = min(length(xs_vec), nsigz)
                    base = _idx(0, t_idx, jg)
                    for jz in 1:lim
                        contrib = xs_vec[jz] * flux * cflux_t[jg]
                        sf_acc[base + jz]  += contrib
                        abs_acc[base + jz] += contrib   # gamff line 1497-1498.
                    end
                end
            end
        end
    end

    # Final pass — gamff lines 1514-1539. Convert raw accumulators to
    # log-transformed f-factors `log(σ_z) - 2*log(ratio)` where
    # ratio = a(loc) / σ_ref(jg). Track iglo = first group where the
    # f-factor differs from 1 by more than `toler` = 0.01.
    iglo  = ngmin
    toler = 0.01
    for ig in ngmin:ngmax
        ab0 = a[locab0 + ig - 1]
        ab0 = ab0 != 0.0 ? 1.0 / ab0 : 0.0
        sf0 = a[locsf0 + ig - 1]
        sf0 = sf0 != 0.0 ? 1.0 / sf0 : 0.0
        for jz in 1:nsigz
            for itp in 1:ntp
                k = _idx(jz, itp, ig)
                # absorption
                v = abs_acc[k]
                v = ab0 == 0.0 ? 1.0 : v * ab0
                v <= 0.0 && (v = 1.0)
                if abs(1.0 - v) > toler && iglo == ngmin
                    iglo = ig
                end
                abs_acc[k] = log(g.sigz[jz]) - 2.0 * log(v)
                # fission (only when jwf set)
                if jwf != 0
                    v = sf_acc[k]
                    v = sf0 == 0.0 ? 1.0 : v * sf0
                    v <= 0.0 && (v = 1.0)
                    sf_acc[k] = log(g.sigz[jz]) - 2.0 * log(v)
                end
            end
        end
    end

    iwr_final = ngmax - iglo + 1
    (abs_acc, sf_acc, iwr_final, jwf)
end

# ============================================================================
# f-factor block writer — powr.f90:614-641
# ============================================================================

# Emits, in order:
#   (4i6) nsigz ntp misc jwf                      ! misc=1 hardcoded
#   (1p,6e12.5) sigz(1..nsigz)
#   (1p,6e12.5) tmpr(1..ntp)
#   (1p,6e12.5) abs f-factors for groups [iglo, ngmax] × ntp × nsigz
#   (1p,6e12.5) fission f-factors (only if jwf != 0)
#   (1p,6e12.5) iwr*10 zero values (dummy amisc array)
#
# iwr is the gamff-computed count = ngmax - iglo + 1. The slice indices into
# the column-major [jz, itp, jg] flat array start at the first group of
# interest: jg=iglo through jg=ngmax, with itp inner-loop and jz innermost.
function _powr_write_ffactor_block(io::IO, abs_ff::Vector{Float64},
                                    sf_ff::Vector{Float64}, jwf::Int,
                                    sigz::Vector{Float64}, tmpr::Vector{Float64},
                                    nsigz::Int, ntp::Int, iwr::Int)
    misc = 1
    println(io, _powr_pad80(@sprintf("%6d%6d%6d%6d", nsigz, ntp, misc, jwf)))
    _powr_write_e125_block(io, sigz)
    _powr_write_e125_block(io, tmpr)

    ngmax = _POWR_FAST_NGMAX
    iglo  = ngmax - iwr + 1
    block_size = nsigz * ntp
    istart = (iglo - 1) * block_size + 1
    iend   = ngmax * block_size
    @assert iend == istart + iwr * block_size - 1
    _powr_write_e125_block(io, view(abs_ff, istart:iend))
    if jwf != 0
        _powr_write_e125_block(io, view(sf_ff, istart:iend))
    end
    # dummy amisc(iwr*10) zero block.
    nwmisc = iwr * 10
    _powr_write_e125_block(io, zeros(Float64, nwmisc))
end

# ============================================================================
# GENDF reader for powr's needs (multi-T, multi-σ₀)
# ============================================================================

# Walks a GENDF file and captures the metadata + per-(temperature, sigma-zero)
# data that powr lib=1 needs:
#   - iza, nsigz_file (from MAT's MF=1/MT=451 HEAD L2),
#   - sigz[1..nsigz_file]    (from MF=1/MT=451 LIST data, offset ntw),
#   - tmpr[1..ntp]           (all temps in tape order; rtemp must be FIRST,
#     matching the gamll error path "reference temperature is not on gendf
#     tape" — see powr.f90:847),
#   - rtemp_idx              (1-based; references rtemp inside tmpr; 0 means
#     no temp matched within `eps` and rtemp > 0 means the caller will error),
#   - mf3[t_idx][mt][ig]   = (flux_p0::Float64, xs_per_sigz::Vector{Float64})
#                            full per-σ₀ XS at Legendre order l=1.
#   - mf3p1[t_idx][mt][ig] = (flux_p1::Float64, xs_p1_per_sigz::Vector{Float64})
#                            for sections with nl >= 2 (MF=3/MT=1 typically),
#   - mf6[mt]              = Vector{NamedTuple} at the reference temperature
#                            only (gamff does not touch MF=6; gamxs reads
#                            MF=6 only at rtemp via `skiprz(-2)+findf(matd,
#                            1, 451)` in powr.f90:1395-1397),
#   - nfs                  delayed-neutron spectrum count.
#
# LIST record data layout for MF=3 (NJOY GENDF, ng2=2):
#   data[1..nl*nz]                     fluxes per (Legendre l, σ-zero jz).
#   data[nl*nz + 1..2*nl*nz]           cross-sections per (l, jz).
#   Indexing: data[(k-1)*nl*nz + (jz-1)*nl + l] for k ∈ {1=flux, 2=xs}.
#   The σ₀=∞ XS at l=1 lives at data[nl*nz + 1] — Fortran's gamxs uses this
#   directly (`scr(nl*jz+loca)` for jz=1 with `loca = l+lz+nl*(nz-1)` reduces
#   to `scr(7+nl*nz)`); for the multi-σ₀ f-factor accumulation gamff loops
#   jz=1..nsigz and reads `data[nl*(jz+nz-1) + 1]` (= the l=1 component at
#   σ-zero jz). We capture the full nz vector and let the call sites slice.
function _powr_read_gendf_for_fast(gendf_path::String, matd::Int, rtemp::Float64;
                                    nsgz_user::Int=0, ntpmax::Int=5, nszmax::Int=7)
    iza = 0
    nsigz_file = 1
    nfs = 0
    sigz = Float64[]
    tmpr = Float64[]
    rtemp_idx = 0

    mf3   = Vector{Dict{Int, Dict{Int, Tuple{Float64, Vector{Float64}}}}}()
    mf3p1 = Vector{Dict{Int, Dict{Int, Tuple{Float64, Vector{Float64}}}}}()
    mf6   = Dict{Int, Vector{NamedTuple}}()

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

    # Resolve (or create) the per-temp dict slot for `temp`.
    # Mirrors gamll's discovery order: rtemp must come FIRST in the tape
    # (else Fortran errors at powr.f90:847 with "reference temperature is
    # not on gendf tape"). Subsequent temps populate slots 2..ntp.
    function _temp_slot!(temp::Float64)
        if rtemp_idx == 0
            if abs(temp - rtemp) <= eps
                push!(tmpr, temp); rtemp_idx = length(tmpr)
            else
                temp > rtemp && error(
                    "powr lib=1: reference temperature $rtemp K is not on " *
                    "gendf tape (encountered $temp K first; tape must be in " *
                    "ascending temperature order). Ref: powr.f90:847.")
                # else: pre-reference temp; not stored, gamff would never see it.
                return 0
            end
        else
            slot = findfirst(t -> abs(t - temp) <= eps, tmpr)
            if slot === nothing
                length(tmpr) >= 1 + ntpmax && error(
                    "powr lib=1: too many additional temperatures " *
                    "(found $(length(tmpr) + 1), ntpmax = $ntpmax). " *
                    "Ref: powr.f90:850-852.")
                push!(tmpr, temp)
                slot = length(tmpr)
            end
        end
        # Make sure per-temp containers exist for this slot.
        while length(mf3) < length(tmpr)
            push!(mf3,   Dict{Int, Dict{Int, Tuple{Float64, Vector{Float64}}}}())
            push!(mf3p1, Dict{Int, Dict{Int, Tuple{Float64, Vector{Float64}}}}())
        end
        rtemp_idx == 0 ? 0 : findfirst(t -> abs(t - temp) <= eps, tmpr)
    end

    # Stage 1: locate matd's MF=1/MT=451 to harvest iza, nsigz_file, sigz.
    idx = 1
    while idx <= nlines
        length(lines[idx]) < 75 && (idx += 1; continue)
        m = _line_meta(idx)
        if m.mat == matd && m.mf == 1 && m.mt == 451 && m.seq == 1
            iza        = round(Int, parse_endf_float(m.line[1:11]))
            nsigz_file = max(1, _parse_int(m.line[34:44]))   # HEAD L2.
            ntw_head   = _parse_int(m.line[56:66])           # HEAD N2H.
            # First LIST record after HEAD has the ntw temp slots + nsigz σ₀'s.
            if idx + 1 <= nlines
                p_cont = rpad(lines[idx + 1], 80)
                nw_list = _parse_int(p_cont[45:55])
                if nw_list >= ntw_head + nsigz_file && idx + 2 <= nlines
                    list_data = _read_list_data(idx + 2, nw_list)
                    sigz = Float64[list_data[ntw_head + i] for i in 1:nsigz_file]
                end
            end
            break
        end
        idx += 1
    end

    # Apply user nsgz cap (gamll line 837): nsigz = min(nsigz_file, nsgz_user)
    # when nsgz_user > 0.
    nsigz = (nsgz_user > 0) ? min(nsigz_file, nsgz_user) : nsigz_file
    nsigz > nszmax && error("powr lib=1: nsigz=$nsigz exceeds nszmax=$nszmax. " *
                            "Ref: powr.f90:838-840.")

    # Stage 2: walk all (MAT, MF, MT) records, binning per-temp.
    idx = 1
    while idx <= nlines
        length(lines[idx]) < 75 && (idx += 1; continue)
        m = _line_meta(idx)
        m.mat != matd && (idx += 1; continue)

        if m.mf == 3 && m.mt > 0 && m.seq == 1
            mt = m.mt
            nl_section = _parse_int(m.line[23:33])
            nz_section = max(1, _parse_int(m.line[34:44]))
            idx += 1
            while idx <= nlines
                length(lines[idx]) < 75 && (idx += 1; continue)
                p2 = rpad(lines[idx], 80)
                m2 = (mat=_parse_int(p2[67:70]), mf=_parse_int(p2[71:72]),
                      mt=_parse_int(p2[73:75]))
                (m2.mat != matd || m2.mf != 3 || m2.mt != mt) && break
                temp = parse_endf_float(p2[1:11])
                nw   = _parse_int(p2[45:55])
                ig   = _parse_int(p2[56:66])
                if ig > 0 && nw >= nl_section * nz_section * 2 &&
                   idx + 1 <= nlines
                    t_idx = _temp_slot!(temp)
                    if t_idx > 0
                        data = _read_list_data(idx + 1, nw)
                        nl   = nl_section
                        nz   = nz_section
                        flux_p0 = data[1]
                        # XS per σ₀ at Legendre order l=1: data[nl*(jz+nz-1) + 1].
                        xs_vec = Float64[data[nl*(jz + nz - 1) + 1] for jz in 1:nz]
                        mt_data = get!(mf3[t_idx], mt,
                                       Dict{Int, Tuple{Float64, Vector{Float64}}}())
                        mt_data[ig] = (flux_p0, xs_vec)
                        if nl >= 2
                            flux_p1 = data[2]
                            xs_p1_vec = Float64[data[nl*(jz + nz - 1) + 2] for jz in 1:nz]
                            mt_data_p1 = get!(mf3p1[t_idx], mt,
                                              Dict{Int, Tuple{Float64, Vector{Float64}}}())
                            mt_data_p1[ig] = (flux_p1, xs_p1_vec)
                        end
                    end
                end
                idx += 1 + cld(nw, 6)
            end

        elseif m.mf == 6 && m.mt > 0 && m.seq == 1
            mt = m.mt
            nl_section = _parse_int(m.line[23:33])
            nz_section = max(1, _parse_int(m.line[34:44]))
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
                if abs(temp - rtemp) <= eps && ig >= 0 && nw > 0 &&
                   idx + 1 <= nlines
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

    rtemp_idx == 0 && rtemp > 0 && error(
        "powr lib=1: reference temperature $rtemp K not found on gendf tape " *
        "(tape carries temps $tmpr). Ref: powr.f90:847.")

    (mf3=mf3, mf3p1=mf3p1, mf6=mf6, iza=iza, nsigz=nsigz, nfs=nfs,
     sigz=sigz, tmpr=tmpr, rtemp_idx=rtemp_idx)
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
