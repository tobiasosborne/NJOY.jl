# resxsr module runner — Build a CCCC RESXS resonance cross-section file.
#
# Reads N PENDF input tapes; for each material, locates MT=2 (elastic),
# MT=18 (fission, optional), MT=102 (capture); builds a union energy grid
# across all reactions and temperatures; thins via 3-point linear-deviation
# test; writes a CCCC binary file with per-material control + cross-section
# blocks.
#
# Output format (resxsr.f90 lines 43-189) — Fortran unformatted records:
#   R1 file ID:    "resxsr  ", huse(1) (a8), huse(2) (a8), ivers (i32)
#   R2 file ctrl:  efirst (f32), elast (f32), nholl (i32), nmat (i32),
#                  nblok (i32 = 5000)
#   R3 hsetid:     nholl × 8 bytes — Fortran FORCES blanks (line 460-462).
#                  The user descriptive cards are NOT written here.
#   R4 file data:  (hmatn[i], i=1,nmat) (a8) ‖ (ntemp[i]) (i32) ‖ (locm[i]) (i32)
#   Per material:
#     Material control: hmat (a8), amass (f32), temp[1..ntemp] (f32),
#                       nreac (i32), nener (i32)
#     XS block(s):      packed [E, σ_R1_T1, σ_R2_T1, …, σ_RN_T1, σ_R1_T2, …]
#                       per energy, chunked into nblok-word records
#                       (nwds = floor(nblok/nn) × nn).
#
# All real data is single-precision (k4); strings are double-precision (k8 = 8
# bytes); integers are 4-byte. mult=2 (8-byte words for a8). Records use
# little-endian length markers (htol) to mirror Fortran unformatted I/O on
# x86_64 Linux.
#
# Ref: njoy-reference/src/resxsr.f90:10-502 (subroutine resxsr).

"""
    resxsr_module(tapes::TapeManager, params::ResxsrParams)

Build a CCCC RESXS file at `resolve(tapes, params.nout)`. Per material:
read the requested PENDF, gather MT=2/18/102 cross sections per temperature,
build the union grid in `[efirst, elast]`, thin to the eps-tolerance, and
write the binary record sequence above.
"""
function resxsr_module(tapes::TapeManager, params::ResxsrParams)
    @info "resxsr: nout=$(params.nout) nmat=$(params.nmat) maxt=$(params.maxt) " *
          "efirst=$(params.efirst) elast=$(params.elast) eps=$(params.eps)"

    materials_data = Vector{NamedTuple}(undef, params.nmat)
    for (im, matspec) in enumerate(params.materials)
        materials_data[im] = _resxsr_process_material(tapes, matspec, params)
    end

    out_path = resolve(tapes, params.nout)
    _resxsr_write_tape(out_path, params, materials_data)
    register!(tapes, params.nout, out_path)
    nothing
end

# ============================================================================
# Per-material processing: union grid + multi-temp interp + thinning
# ============================================================================

# Reactions resxsr collects (Fortran resxsr.f90:300 — `mth.eq.2.or.mth.eq.18.or.mth.eq.102`).
const _RESXSR_REACTIONS = (2, 18, 102)

function _resxsr_process_material(tapes::TapeManager, matspec::ResxsrMaterial,
                                    params::ResxsrParams)
    path = resolve(tapes, matspec.unit)
    isfile(path) || error("resxsr: input tape unit $(matspec.unit) not found")
    tape = read_pendf(path)

    # PENDF can hold multiple temperatures for the same MAT. Each temperature
    # appears as a separate PENDFMaterial entry (same .mat, different MF=1/451
    # TEMP value). Limit to maxt as Fortran does (line 273).
    mat_records = filter(m -> m.mat == matspec.mat, tape.materials)
    isempty(mat_records) &&
        error("resxsr: MAT=$(matspec.mat) not found on tape unit $(matspec.unit)")
    ntemp = min(length(mat_records), params.maxt)
    mat_records = mat_records[1:ntemp]

    # Per-temperature: gather (mt, energies, xs) for each present reaction.
    per_temp = Vector{NamedTuple}(undef, ntemp)
    awr = 0.0
    for (it, mat_rec) in enumerate(mat_records)
        meta = _resxsr_meta(mat_rec)
        awr  = meta.awr
        rxs  = Tuple{Int, Vector{Float64}, Vector{Float64}}[]
        for mt in _RESXSR_REACTIONS
            for sec in mat_rec.sections
                if sec.mf == 3 && sec.mt == mt
                    (e, xs) = _parse_mf3_lines(sec.lines)
                    push!(rxs, (mt, e, xs))
                    break
                end
            end
        end
        per_temp[it] = (temp=meta.temp, reactions=rxs)
    end

    # Reaction MT layout must match across temperatures (Fortran assumes this).
    mts_t1 = [r[1] for r in per_temp[1].reactions]
    for it in 2:ntemp
        mts_it = [r[1] for r in per_temp[it].reactions]
        mts_it == mts_t1 || error(
            "resxsr: MAT=$(matspec.mat) reaction MTs differ between temps " *
            "($mts_t1 vs $mts_it) — not supported")
    end
    nreac = length(mts_t1)
    nreac > 0 || error("resxsr: MAT=$(matspec.mat) has no MT=2/18/102 sections")

    # Union grid clipped to [efirst, elast] (Fortran: every gety1 walk clamps
    # against efirst/elast — lines 262-265, 339).
    all_e = Float64[]
    for tdata in per_temp
        for (_, e, _) in tdata.reactions
            append!(all_e, e)
        end
    end
    sort!(all_e); unique!(all_e)
    filter!(e -> params.efirst <= e <= params.elast, all_e)
    isempty(all_e)            && (all_e = Float64[params.efirst, params.elast])
    all_e[1] != params.efirst && pushfirst!(all_e, params.efirst)
    all_e[end] != params.elast && push!(all_e, params.elast)

    # Build pre-thin matrix: rows = energies; cols = [E, R1_T1, R2_T1, …,
    # R1_T2, …, RN_TN]. Per-temperature interpolation onto union grid.
    ne   = length(all_e)
    ncol = 1 + nreac * ntemp
    pre_thin = zeros(Float64, ne, ncol)
    @inbounds for i in 1:ne
        pre_thin[i, 1] = all_e[i]
    end
    col = 1
    for tdata in per_temp
        for (_, e_src, xs_src) in tdata.reactions
            col += 1
            @inbounds for i in 1:ne
                pre_thin[i, col] = _resxsr_interp_linlin(e_src, xs_src, all_e[i])
            end
        end
    end

    thinned = _resxsr_thin(pre_thin, params.eps)

    (matspec=matspec, awr=awr,
     temps=Float64[t.temp for t in per_temp],
     ntemp=ntemp, nreac=nreac, xs=thinned)
end

# Pull AWR (HEAD line C2) and TEMP (HDATIO line C1) from a PENDFMaterial.
function _resxsr_meta(mat_rec::PENDFMaterial)
    mf1 = mat_rec.mf1_lines
    isempty(mf1) && return (awr=0.0, temp=0.0)

    p1 = rpad(mf1[1], 80)
    awr = parse_endf_float(p1[12:22])

    iverf = 5
    if length(mf1) >= 2
        p2 = rpad(mf1[2], 80)
        n1 = _parse_int(p2[45:55])
        n2 = _parse_int(p2[56:66])
        iverf = n1 != 0 ? 4 : (n2 == 0 ? 5 : 6)
    end
    hdatio_idx = iverf >= 6 ? 4 : (iverf >= 5 ? 3 : 2)
    temp = 0.0
    if length(mf1) >= hdatio_idx
        ph = rpad(mf1[hdatio_idx], 80)
        temp = parse_endf_float(ph[1:11])
    end
    (awr=awr, temp=temp)
end

# ============================================================================
# Thinning — 3-point linear-deviation test (Fortran resxsr.f90:355-396)
# ============================================================================

# Maintain a sliding stack. After each new point, if the stack has ≥3 points
# and ANY interior point's column j (j>1) deviates from the linear interp
# between stack[1] and stack[end] by more than eps × value, emit stack[1]
# and restart the stack with [stack[end-1], stack[end]]. The very last input
# point is emitted verbatim regardless of stack contents (resxsr.f90:366
# `if (ie.eq.ne) go to 350` — the algorithm intentionally drops uncommitted
# stack contents at end-of-run; we mirror this for bit-identical output).
function _resxsr_thin(pre_thin::Matrix{Float64}, eps::Float64)
    ne, ncol = size(pre_thin)
    out = Vector{Vector{Float64}}()
    ne == 0 && return out

    if ne == 1
        push!(out, pre_thin[1, :])
        return out
    end

    stack = Vector{Vector{Float64}}()
    n = 0

    for ie in 1:ne
        sigs = pre_thin[ie, :]
        if ie == ne
            push!(out, sigs)
            return out
        end

        n += 1
        if length(stack) < n
            push!(stack, sigs)
        else
            stack[n] = sigs
        end
        n < 3 && continue

        lim = n - 1
        broken = false
        s1 = stack[1]; sn = stack[n]
        @inbounds for i in 2:lim
            si = stack[i]
            for j in 2:ncol
                test = _terp_linlin(s1[1], s1[j], sn[1], sn[j], si[1])
                if abs(test - si[j]) > eps * si[j]
                    broken = true; break
                end
            end
            broken && break
        end

        if broken
            push!(out, copy(stack[1]))
            new1 = stack[lim]
            new2 = stack[n]
            stack[1] = new1; stack[2] = new2
            n = 2
        end
    end

    out
end

# ============================================================================
# Linear interpolation helpers
# ============================================================================

function _resxsr_interp_linlin(e_src::AbstractVector{<:Real},
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

@inline function _terp_linlin(x1, y1, x2, y2, x)
    x2 == x1 && return y1
    y1 + (y2 - y1) * (x - x1) / (x2 - x1)
end

# ============================================================================
# Binary writer — Fortran unformatted records via _write_record from ccccr.jl
# ============================================================================

function _resxsr_write_tape(path::String, params::ResxsrParams,
                              materials_data::Vector{NamedTuple})
    nblok = 5000   # Fortran resxsr.f90:219 — fixed.

    # Pre-compute locm[i] (record offset at which material i's control record
    # appears, counting from after the file-header records). Fortran irec
    # increments per scratch write — locm captures the value BEFORE the new
    # material's records are emitted (line 253 `locm(im)=irec`).
    locm = Vector{Int32}(undef, length(materials_data))
    irec = 0
    for (i, md) in enumerate(materials_data)
        locm[i] = Int32(irec)
        irec += 1                                          # material control
        nn = 1 + md.nreac * md.ntemp
        nwds = (nblok ÷ nn) * nn
        nw   = length(md.xs) * nn
        nb   = nwds == 0 ? 1 : cld(nw, nwds)
        irec += nb                                          # XS block(s)
    end

    open(path, "w") do io
        # R1: file identification
        huse_padded = rpad(params.huse, 12)[1:12]
        huse1 = rpad(huse_padded[1:6], 8)
        huse2 = rpad(huse_padded[7:12], 8)
        _write_record(io, _record_buf("resxsr  ", huse1, huse2,
                                       Int32(params.ivers)))

        # R2: file control
        _write_record(io, _record_buf(
            Float32(params.efirst), Float32(params.elast),
            Int32(params.nholl), Int32(params.nmat), Int32(nblok)))

        # R3: set hollerith ID — Fortran forces blanks (resxsr.f90:460-462).
        _write_record(io, repeat(b"        ", params.nholl))

        # R4: file data — names ‖ ntemps ‖ locm
        buf = IOBuffer()
        for md in materials_data
            write(buf, codeunits(rpad(md.matspec.hmat, 8)[1:8]))
        end
        for md in materials_data
            write(buf, htol(Int32(md.ntemp)))
        end
        for v in locm
            write(buf, htol(v))
        end
        _write_record(io, take!(buf))

        # Per material: material control + xs blocks
        for md in materials_data
            mc_buf = IOBuffer()
            write(mc_buf, codeunits(rpad(md.matspec.hmat, 8)[1:8]))
            write(mc_buf, htol(Float32(md.awr)))
            for t in md.temps
                write(mc_buf, htol(Float32(t)))
            end
            write(mc_buf, htol(Int32(md.nreac)))
            write(mc_buf, htol(Int32(length(md.xs))))   # nener
            _write_record(io, take!(mc_buf))

            nn = 1 + md.nreac * md.ntemp
            nwds_per_block = (nblok ÷ nn) * nn   # word capacity per record

            block_buf = IOBuffer()
            words_in_block = 0
            for row in md.xs
                length(row) == nn || error(
                    "resxsr: thinned row width $(length(row)) ≠ expected $nn")
                for v in row
                    write(block_buf, htol(Float32(v)))
                    words_in_block += 1
                end
                if words_in_block >= nwds_per_block
                    _write_record(io, take!(block_buf))
                    words_in_block = 0
                end
            end
            if words_in_block > 0
                _write_record(io, take!(block_buf))
            end
        end
    end
end
