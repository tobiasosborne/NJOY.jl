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

    # Walk the GENDF and accumulate per-MT per-group flux + xs at the
    # reference temperature. Returns a NamedTuple with the parsed data.
    g = _powr_read_gendf_for_fast(gendf_path, m.matd, m.rtemp)

    # gamll: detect iwa / iwf / iwr / xlol presence.
    has_fission   = any(haskey(g.mf3, mt) for mt in (18, 19, 20, 21, 38))
    has_mf6_match = !isempty(g.mf6_present)
    iwa = 1
    iwf = has_fission ? 1 : 0
    iwr = (g.nsigz > 1 && m.iff > 0) ? 1 : 0
    nfs = g.nfs    # delayed-neutron spectrum count from MF=5/MT=455

    if iwf != 0 || iwr != 0 || has_mf6_match || nfs > 0
        error("powr lib=1: full Phase B path not yet ported (matd=$(m.matd) " *
              "has iwf=$iwf iwr=$iwr mf6=$has_mf6_match nfs=$nfs). " *
              "Carbon-style absorption-only path is the only current support.")
    end

    nid = g.iza * 10
    ngnd = _POWR_FAST_NGND

    # gamxs: accumulate cflux and absorption.
    cflux = zeros(Float64, ngnd)
    if haskey(g.mf3, 1)
        for (ig, (flux, _xs)) in g.mf3[1]
            jg = ngnd - ig + 1
            1 <= jg <= ngnd && (cflux[jg] += flux)
        end
    end
    @inbounds for jg in 1:ngnd
        cflux[jg] = cflux[jg] != 0 ? 1 / cflux[jg] : 0.0
    end
    absxs = zeros(Float64, ngnd)
    for mt in keys(g.mf3)
        # Fortran gamxs.f90:1157 — `if (mth.lt.102.or.mth.gt.150) go to 280`,
        # i.e., MT ∈ [102, 150] → capture / absorption channel.
        102 <= mt <= 150 || continue
        for (ig, (flux, xs)) in g.mf3[mt]
            jg = ngnd - ig + 1
            1 <= jg <= ngnd && (absxs[jg] += xs * flux * cflux[jg])
        end
    end

    # gamll xlol: with no MF=6, xld(1..3) = 0 → xlol(1..3) = 0.
    xlol = (0.0, 0.0, 0.0)
    xla  = (0.0, 0.0, 0.0)
    xld  = (0, 0, 0)        # int form
    iwa_field = iwa
    iwf_field = iwf
    iwr_field = iwr

    # Write output records in fast()-driver order (lines 528-547).
    println(io, _powr_pad80(@sprintf("%10d%10d%10d%10d", nid, iwa_field, iwf_field, iwr_field)))
    println(io, _powr_pad80(_powr_word16(m.word)))
    println(io, _powr_pad80(repeat(" ", 36) *
                            _powr_e125(xlol[3]) * _powr_e125(xla[3]) *
                            _powr_e125(Float64(xld[3] + 1))))
    println(io, _powr_pad80(_powr_e125(xlol[1]) * _powr_e125(xla[1]) *
                            _powr_e125(Float64(xld[1] + 1)) *
                            _powr_e125(xlol[2]) * _powr_e125(xla[2]) *
                            _powr_e125(Float64(xld[2] + 1))))
    _powr_write_e125_block(io, absxs)
    nothing
end

# ============================================================================
# Minimal GENDF reader for powr's needs
# ============================================================================

# For each MT in MF=3, returns a Dict{ig => (flux, xs)} at the reference
# temperature. Also captures iza, nsigz, nfs, and any MF=6 group counts.
function _powr_read_gendf_for_fast(gendf_path::String, matd::Int, rtemp::Float64)
    mf3 = Dict{Int, Dict{Int, Tuple{Float64, Float64}}}()
    mf6_present = Set{Int}()
    iza = 0
    nsigz = 1
    nfs = 0

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

    # Walk to first record of matd
    idx = 1
    while idx <= nlines
        length(lines[idx]) < 75 && (idx += 1; continue)
        m = _line_meta(idx)
        m.mat == matd && break
        idx += 1
    end
    idx <= nlines || error("powr: MAT=$matd not found on GENDF $gendf_path")

    # Cache iza from the very first MF=1/MT=451 HEAD of this material.
    while idx <= nlines
        m = _line_meta(idx)
        m.mat != matd && break
        if m.mf == 1 && m.mt == 451 && m.seq == 1
            iza = round(Int, parse_endf_float(m.line[1:11]))
            # The 2nd CONT carries (NL, NZ, NW, NGN) in MF1/451 of GENDF.
            if idx + 1 <= nlines
                p2 = rpad(lines[idx + 1], 80)
                nsigz = _parse_int(p2[34:44])
            end
            break
        end
        idx += 1
    end

    # Walk all sections; for MF=3 LIST records at temp ≈ rtemp, capture
    # (flux, xs) per group. For MF=6 record any MTs present (we error out
    # later if any are seen — Phase B doesn't yet handle them). For
    # MF=5/MT=455 capture nfs.
    idx = 1
    while idx <= nlines
        length(lines[idx]) < 75 && (idx += 1; continue)
        m = _line_meta(idx)
        m.mat != matd && (idx += 1; continue)
        if m.mf == 3 && m.mt > 0 && m.seq == 1
            # Section HEAD; subsequent records are LISTs per group.
            mt = m.mt
            mt_data = get!(mf3, mt, Dict{Int, Tuple{Float64, Float64}}())
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
                if abs(temp - rtemp) <= eps && ig > 0 && nw >= 2 && idx + 1 <= nlines
                    p3 = rpad(lines[idx + 1], 80)
                    flux = parse_endf_float(p3[1:11])
                    xs   = parse_endf_float(p3[12:22])
                    mt_data[ig] = (flux, xs)
                end
                # Skip data lines: ceil(nw / 6).
                idx += 1 + cld(nw, 6)
            end
        elseif m.mf == 6 && m.mt > 0 && m.seq == 1
            push!(mf6_present, m.mt)
            idx += 1
        elseif m.mf == 5 && m.mt == 455 && m.seq == 1
            # NL on the 2nd CONT (or LIST) gives the # delayed groups.
            if idx + 1 <= nlines
                p2 = rpad(lines[idx + 1], 80)
                nfs += _parse_int(p2[23:33])
            end
            idx += 1
        else
            idx += 1
        end
    end

    (mf3=mf3, mf6_present=mf6_present, iza=iza, nsigz=nsigz, nfs=nfs)
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
