# mixr module runner — Mix cross sections from one or more PENDF/ENDF tapes.
#
# Reads N input tapes (PENDF or ENDF), locates a (MAT, T) on each, and writes
# a new tape with MF=1/MT=451 + MF=3 sections that are linear combinations of
# the input cross sections:
#
#     sigma_out(E, MT) = Σ wtn[i] * sigma_in[i](E, MT)
#
# The output grid for each MT is the sorted union of all input MT grids;
# outside an input's range that input contributes zero. Linear-linear
# interpolation is always used between input points (Fortran also accepts
# INT≠2 input via terp1 but emits INT=2 only — we do the same). XS values
# are rounded `sigfig(7, 0)` before write (mixr.f90:350).
#
# Output structure (Fortran mixr.f90:198-255):
#   - TPID: 66 blanks (`math=1; afend(nout, 0)`).
#   - MF=1/MT=451: HEAD + (iverf-aware CONTs) + HDATIO + description records
#     + dictionary entries (one per output MT).
#   - MF=3 sections: HEAD + TAB1 (NR=1, INT=2) + (E, σ) data + SEND.
#
# Ref: njoy-reference/src/mixr.f90:17-390 (subroutine mixr).

"""
    mixr_module(tapes::TapeManager, params::MixrParams)

Construct a new PENDF tape with mixed cross sections per `params`. The output
tape is written to `resolve(tapes, params.nout)` and registered.
"""
function mixr_module(tapes::TapeManager, params::MixrParams)
    @info "mixr: nout=$(params.nout) nin=$(params.nin) mts=$(params.mtn) " *
          "mats=$(params.matn) weights=$(params.wtn) temp=$(params.temp)"

    isempty(params.nin) && error("mixr: no input tapes specified")
    isempty(params.mtn) && error("mixr: no output MTs specified")
    length(params.matn) == length(params.nin) ||
        error("mixr: nin / matn length mismatch")

    # 1. Read each input tape and pull the requested (mat, temp) material.
    inputs = Vector{NamedTuple}(undef, length(params.nin))
    for (i, unit) in enumerate(params.nin)
        path = resolve(tapes, unit)
        isfile(path) || error("mixr: input tape unit $unit not found at $path")
        tape = read_pendf(path)
        target_mat = params.matn[i]
        material = nothing
        for m in tape.materials
            if m.mat == target_mat
                material = m; break
            end
        end
        material === nothing &&
            error("mixr: MAT=$target_mat not found on input tape unit $unit")

        meta = _mixr_parse_mf1_meta(material.mf1_lines)
        # Fortran mixr.f90:184-185 — temperature test |t - temp| < 1.
        abs(meta.temp - params.temp) < 1.0 || error(
            "mixr: MAT=$target_mat on tape unit $unit has temp=$(meta.temp), " *
            "requested temp=$(params.temp)")
        inputs[i] = (material=material, meta=meta, weight=params.wtn[i])
    end

    # 2. iverf / awi / nsub from the first input; emax = max across inputs
    #    (Fortran mixr.f90:172-180 + 171).
    iverf = inputs[1].meta.iverf
    awi   = inputs[1].meta.awi
    emax  = inputs[1].meta.emax
    for inp in inputs
        inp.meta.emax > emax && (emax = inp.meta.emax)
    end
    nsub = _mixr_compute_nsub(awi)

    # 3. For each output MT: build union grid + mix.
    out_sections = PENDFSection[]
    actual_mts = Int[]
    for mt in params.mtn
        ranges = NamedTuple[]
        for inp in inputs
            section = nothing
            for sec in inp.material.sections
                if sec.mf == 3 && sec.mt == mt
                    section = sec; break
                end
            end
            if section === nothing
                @warn "mixr: MT=$mt not present for MAT=$(inp.material.mat)"
                continue
            end
            (e, xs) = _parse_mf3_lines(section.lines)
            isempty(e) && continue
            push!(ranges, (e=e, xs=xs, w=inp.weight))
        end
        isempty(ranges) && continue

        # Sorted union of all input energies.
        all_e = sort!(unique(vcat([r.e for r in ranges]...)))

        # Mix: Σ w_i · y_i(e); y_i = 0 outside input range, lin-lin within.
        mixed = zeros(Float64, length(all_e))
        for r in ranges
            for (j, e) in enumerate(all_e)
                mixed[j] += r.w * _mixr_interp_linlin(r.e, r.xs, e)
            end
        end

        # sigfig(7, 0) per Fortran mixr.f90:350.
        @inbounds for j in eachindex(mixed)
            mixed[j] = round_sigfig(mixed[j], 7, 0)
        end

        push!(out_sections,
              PENDFSection(3, mt, _mixr_build_mf3_section(
                  all_e, mixed, params.matd, mt, params.za, params.awr)))
        push!(actual_mts, mt)
    end

    # 4. MF=1/MT=451 (uses actual_mts so dictionary matches what was written).
    mf1_lines = _mixr_build_mf1_451(params, iverf, awi, emax, nsub, actual_mts)

    # 5. Write tape directly. Mixr's tape layout has two Fortran-specific
    #    quirks vs the generic PENDF writer: TPID is 66 blanks, and there is
    #    NO SEND record between the MF=1/MT=451 dictionary and the FEND
    #    (Fortran `dictio` writes entries without auto-emitting SEND, and
    #    mixr.f90:255 calls `afend` directly — see oracle line 9).
    out_path = resolve(tapes, params.nout)
    _mixr_write_tape(out_path, params.matd, mf1_lines, out_sections)
    register!(tapes, params.nout, out_path)
    nothing
end

# ============================================================================
# Tape writer — Fortran mixr.f90:198-378 emits this exact structure
# ============================================================================

function _mixr_write_tape(path::String, matd::Int,
                           mf1_lines::Vector{String},
                           sections::Vector{PENDFSection})
    open(path, "w") do io
        # TPID: 66 blanks + (MAT=1, MF=0, MT=0, seq=0). Matches Fortran
        # `math=1; afend(nout, 0)` (mixr.f90:198-199).
        @printf(io, "%66s%4d%2d%3d%5d\n", "", 1, 0, 0, 0)

        # MF=1/MT=451: write each line, then go straight to FEND (no SEND).
        seq = 0
        for line in mf1_lines
            seq += 1
            _mixr_write_line(io, line, matd, seq)
        end
        # FEND directly (mixr.f90:255 — `afend(nout, 0)` skips the SEND).
        @printf(io, "%66s%4d%2d%3d%5d\n", "", matd, 0, 0, 0)

        # MF=3 sections: each closed with SEND (mixr.f90:370 — `asend`).
        # No per-section FEND: mixr only emits MF=1 + MF=3, and the FEND
        # after MT=451 above already separates them. The trailing FEND below
        # closes MF=3.
        for sec in sections
            seq = 0
            for line in sec.lines
                seq += 1
                _mixr_write_line(io, line, matd, seq)
            end
            @printf(io, "%66s%4d%2d%3d%5d\n", "", matd, sec.mf, 0, 99999)
        end

        # Final FEND for last MF, MEND, TEND (mixr.f90:376-378).
        @printf(io, "%66s%4d%2d%3d%5d\n", "", matd, 0, 0, 0)
        @printf(io, "%66s%4d%2d%3d%5d\n", "", 0, 0, 0, 0)
        @printf(io, "%66s%4d%2d%3d%5d\n", "", -1, 0, 0, 0)
    end
end

function _mixr_write_line(io::IO, line::AbstractString, mat::Int, seq::Int)
    data = rpad(line, 80)
    @printf(io, "%s%5d\n", data[1:75], seq)
end

# ============================================================================
# MF=1/MT=451 metadata parsing
# ============================================================================

# iverf detection mirrors Fortran mixr.f90:159-167. Layout for each version:
#   iverf=4: HEAD (with NX in N2)              + HDATIO
#   iverf=5: HEAD + 2nd CONT (all zero)        + HDATIO
#   iverf=6: HEAD + 2nd CONT (NMOD=6 in N2)    + 3rd CONT (AWI/EMAX/NSUB) + HDATIO
function _mixr_parse_mf1_meta(mf1_lines::Vector{String})
    isempty(mf1_lines) && error("mixr: MF=1/MT=451 has no records")

    iverf = 5
    awi   = 1.0
    emax  = 20.0e6
    temp  = 0.0

    if length(mf1_lines) >= 2
        p2 = rpad(mf1_lines[2], 80)
        n1 = _parse_int(p2[45:55])
        n2 = _parse_int(p2[56:66])
        if n1 != 0
            iverf = 4
        elseif n2 == 0
            iverf = 5
        else
            iverf = 6
        end
    end

    # AWI / EMAX live on the 3rd CONT only when iverf >= 6.
    if iverf >= 6 && length(mf1_lines) >= 3
        p3 = rpad(mf1_lines[3], 80)
        awi  = parse_endf_float(p3[1:11])
        emax = parse_endf_float(p3[12:22])
    end

    # TEMP lives on the HDATIO CONT — line index depends on iverf.
    hdatio_idx = iverf >= 6 ? 4 : (iverf >= 5 ? 3 : 2)
    if length(mf1_lines) >= hdatio_idx
        ph = rpad(mf1_lines[hdatio_idx], 80)
        temp = parse_endf_float(ph[1:11])
    end

    (iverf=iverf, awi=awi, emax=emax, temp=temp)
end

# Fortran mixr.f90:174-179 — incident-particle ID from AWI.
function _mixr_compute_nsub(awi::Float64)
    awi < 0.1                       && return 0      # photon
    (0.9980 < awi < 0.9990)         && return 10010  # H-1
    (1.9950 < awi < 1.9970)         && return 10020  # H-2
    (2.9895 < awi < 2.9897)         && return 10030  # H-3
    (2.9890 < awi < 2.9891)         && return 20030  # He-3
    (3.9670 < awi < 3.9680)         && return 20040  # He-4
    return 10                                        # neutron (default)
end

# ============================================================================
# MF=1/MT=451 builder — Fortran mixr.f90:198-255
# ============================================================================

function _mixr_build_mf1_451(params::MixrParams, iverf::Int, awi::Float64,
                              emax::Float64, nsub::Int, actual_mts::Vector{Int})
    nx = length(actual_mts)
    lines = String[]
    trailer = @sprintf("%4d%2d%3d", params.matd, 1, 451)
    lrp = -1   # mixr.f90:67

    # HEAD — N2 carries NX only for iverf <= 4 (Fortran mixr.f90:210).
    head_n2 = iverf <= 4 ? nx : 0
    push!(lines,
        format_endf_float(params.za) * format_endf_float(params.awr) *
        _i11(lrp) * _i11(0) * _i11(0) * _i11(head_n2) * trailer * "    0")

    # 2nd CONT for iverf > 4 — all zero except N2=6 when iverf == 6.
    if iverf > 4
        c2_n2 = iverf == 6 ? 6 : 0
        push!(lines,
            format_endf_float(0.0) * format_endf_float(0.0) *
            _i11(0) * _i11(0) * _i11(0) * _i11(c2_n2) * trailer * "    0")
    end

    # 3rd CONT for iverf >= 6 — AWI, EMAX, 0, 0, NSUB, 0.
    if iverf >= 6
        push!(lines,
            format_endf_float(awi) * format_endf_float(emax) *
            _i11(0) * _i11(0) * _i11(nsub) * _i11(0) * trailer * "    0")
    end

    # HDATIO CONT — TEMP, 0, LREL, 0, NWD, NX.
    # NWD is the count of 66-char description records (Fortran NJOY hdatio
    # reports records, not floats, even though the input buffer ndes=17 floats
    # = 1 record). Matches oracle byte-for-byte.
    lrel = iverf >= 6 ? 1 : 0
    nwd  = max(1, cld(length(params.des), 66))
    nx_field = iverf >= 5 ? nx : 0
    push!(lines,
        format_endf_float(params.temp) * format_endf_float(0.0) *
        _i11(lrel) * _i11(0) * _i11(nwd) * _i11(nx_field) * trailer * "    0")

    # Description records — pad/truncate each chunk to exactly 66 chars.
    for k in 1:nwd
        from = (k - 1) * 66 + 1
        to   = min(length(params.des), k * 66)
        chunk = from <= length(params.des) ? params.des[from:to] : ""
        push!(lines, rpad(chunk, 66) * trailer * "    0")
    end

    # Dictionary entries — Fortran `dictio` writes (22x, 4i11) per record:
    # cols 1-22 BLANK, then (MF, MT, NC, MOD). C1/C2 are not floats here.
    for mt in actual_mts
        push!(lines,
            " "^22 * _i11(3) * _i11(mt) * _i11(0) * _i11(0) *
            trailer * "    0")
    end

    lines
end

_i11(n::Integer) = lpad(string(n), 11)

# ============================================================================
# MF=3 section builder
# ============================================================================

function _mixr_build_mf3_section(es::Vector{Float64}, xss::Vector{Float64},
                                  mat::Int, mt::Int, za::Float64, awr::Float64)
    np = length(es)
    lines = String[]
    trailer = @sprintf("%4d%2d%3d", mat, 3, mt)

    # HEAD: ZA, AWR, 0, 0, 0, 0
    push!(lines,
        format_endf_float(za) * format_endf_float(awr) *
        _i11(0) * _i11(0) * _i11(0) * _i11(0) * trailer * "    0")

    # TAB1 CONT: 0, 0, 0, 0, NR=1, NP
    push!(lines,
        format_endf_float(0.0) * format_endf_float(0.0) *
        _i11(0) * _i11(0) * _i11(1) * _i11(np) * trailer * "    0")

    # Interp: NBT=NP, INT=2 (lin-lin).
    push!(lines, _i11(np) * _i11(2) * " "^44 * trailer * "    0")

    # Data: 3 (E, XS) pairs per line.
    idx = 1
    while idx <= np
        buf = ""
        for _ in 1:3
            idx > np && break
            buf *= format_endf_float(es[idx]) * format_endf_float(xss[idx])
            idx += 1
        end
        push!(lines, rpad(buf, 66) * trailer * "    0")
    end

    lines
end

# ============================================================================
# Linear-linear interpolation (returns 0 outside [e_src[1], e_src[end]])
# ============================================================================

function _mixr_interp_linlin(e_src::AbstractVector{<:Real},
                              xs_src::AbstractVector{<:Real}, e::Real)
    n = length(e_src)
    n == 0 && return 0.0
    e < e_src[1]   && return 0.0
    e > e_src[end] && return 0.0
    e == e_src[1]   && return Float64(xs_src[1])
    e == e_src[end] && return Float64(xs_src[end])
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >> 1
        e_src[mid] <= e ? (lo = mid) : (hi = mid)
    end
    x1, x2 = e_src[lo], e_src[lo + 1]
    y1, y2 = xs_src[lo], xs_src[lo + 1]
    x2 == x1 && return Float64(y1)
    y1 + (e - x1) * (y2 - y1) / (x2 - x1)
end
