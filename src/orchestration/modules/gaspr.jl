# gaspr module runner — Gas production cross sections.
#
# Reads the heatr-output PENDF, computes MT=203/204/205/206/207 (proton,
# deuteron, triton, He-3, alpha production) by summing partial reactions
# weighted by particle multiplicities (`accumulate_gas` in
# src/processing/gaspr.jl), then writes a new PENDF with these sections
# spliced in and the MF1/MT451 directory updated.
#
# Ref: njoy-reference/src/gaspr.f90 (subroutine gaspr).
#
# Current coverage:
#  - Constant per-MT yields (MT=11,22-45,103-117,154-200): YES
#  - MF6/MT5 emission spectra (energy-dependent yields, izap=1001..2004): NO
#  - MT=51-91 LR=22..36 competitive flags: NO (Fortran lr field; defer)
#  - Multi-temperature loop: NO (single-temperature only; matches
#    99% of NJOY usage)
#  - Existing MT=203..207 deletion: YES (filtered out before re-emit)

"""
    gaspr_module(tapes::TapeManager, params::GasprParams)

Add gas production cross sections (MT=203..207) to the PENDF tape.

Reads `params.npendf_in`, accumulates gas production from existing MF3
sections via the multiplicity table in `accumulate_gas`, builds new MT=203..207
sections on a subset of MT=1's energy grid (starting at the gas threshold,
backed up one point per channel), and writes the augmented tape to
`params.npendf_out`. The MF1/MT451 directory is updated to list the new
sections (MOD=1).

Defers MF6/MT5 and MT=51-91 LR-flag yields — see header comment.
"""
function gaspr_module(tapes::TapeManager, params::GasprParams)
    @info "gaspr: nendf=$(params.nendf) npendf_in=$(params.npendf_in) " *
          "npendf_out=$(params.npendf_out)"

    if params.npendf_out <= 0
        @warn "gaspr: no output unit — nothing to do"
        return nothing
    end

    out_path = resolve(tapes, params.npendf_out)

    if params.npendf_in <= 0
        touch(out_path)
        register!(tapes, params.npendf_out, out_path)
        return nothing
    end

    in_path = resolve(tapes, params.npendf_in)
    isfile(in_path) ||
        error("gaspr: input PENDF (unit $(params.npendf_in)) not found at $in_path")

    tape = read_pendf(in_path)

    new_materials = PENDFMaterial[]
    for material in tape.materials
        push!(new_materials, _gaspr_one_material(material))
    end

    write_pendf_tape(out_path, PENDFTape(tape.tpid, new_materials))
    register!(tapes, params.npendf_out, out_path)
    nothing
end

# =========================================================================
# One material: extract MF3, accumulate gas, splice MT=203..207
# =========================================================================

function _gaspr_one_material(material::PENDFMaterial)::PENDFMaterial
    mat = material.mat

    # Pull every MF3 section's (energies, xs)
    mf3 = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
    for sec in material.sections
        sec.mf == 3 || continue
        e, xs = _parse_mf3_lines(sec.lines)
        isempty(e) && continue
        mf3[sec.mt] = (e, xs)
    end

    # Need MT=1 to provide the gas-grid skeleton (Fortran gaspr.f90:435 —
    # "use the energy grid of mt1 starting at thrg")
    haskey(mf3, 1) || begin
        @info "gaspr: no MT=1 in MAT=$mat — copying material unchanged"
        return material
    end

    # Identify gas-producing MTs and the threshold (min first-energy among them).
    # Fortran sets thrg = min(enext) across all MTs that pass the gas filter
    # (gaspr.f90:285-424).
    gas_mts_present = Int[]
    thrg = Inf
    for (mt, (e, _)) in mf3
        mt in _GASPR_SKIP_MTS && continue
        any(!iszero, gas_yield(mt)) || continue
        push!(gas_mts_present, mt)
        e[1] < thrg && (thrg = e[1])
    end

    if isempty(gas_mts_present)
        @info "gaspr: no gas-producing reactions in MAT=$mat — copying unchanged"
        return material
    end

    # Build gas grid = MT=1 grid restricted to e >= thrg
    e1, _ = mf3[1]
    gas_grid = Float64[e for e in e1 if e >= thrg]
    ngas = length(gas_grid)
    ngas == 0 && return material

    # Accumulate gas production: sgas[ip,d,t,h3,a × ngas]
    # Each gas-producing MT is interpolated onto gas_grid and weighted.
    # PENDF after broadr is linearized, so LinLin matches Fortran terpa default.
    sgas = zeros(5, ngas)
    for mt in gas_mts_present
        p, d, t, h3, a = gas_yield(mt)
        e_src, xs_src = mf3[mt]
        for (i, eg) in enumerate(gas_grid)
            y = _interp_linlin_at(e_src, xs_src, eg)
            y == 0.0 && continue
            p  != 0 && (sgas[1, i] += p  * y)
            d  != 0 && (sgas[2, i] += d  * y)
            t  != 0 && (sgas[3, i] += t  * y)
            h3 != 0 && (sgas[4, i] += h3 * y)
            a  != 0 && (sgas[5, i] += a  * y)
        end
    end

    # Per-channel write: find first nonzero, back up by 1 (gaspr.f90:1055-1057),
    # sigfig(7,0) the values (gaspr.f90:1094).
    gas_mts = (203, 204, 205, 206, 207)
    new_sections = Dict{Int, Vector{String}}()
    np_per_mt = Dict{Int, Int}()
    for (jg, mt) in enumerate(gas_mts)
        ii = 0
        @inbounds for i in 1:ngas
            if sgas[jg, i] != 0
                ii = i
                break
            end
        end
        ii == 0 && continue        # all-zero channel: skip section
        i_start = max(ii - 1, 1)
        np = ngas - i_start + 1
        es  = gas_grid[i_start:end]
        xss = [round_sigfig(sgas[jg, i], 7, 0) for i in i_start:ngas]

        # Carry ZA/AWR from material's MF1/MT451 head line for these new sections
        za, awr = _gaspr_za_awr(material.mf1_lines)
        new_sections[mt] = _gaspr_build_mf3_section(es, xss, mat, mt, za, awr)
        np_per_mt[mt] = np
    end

    # Drop any pre-existing MT=203..207 sections (Fortran "delete old gas sections")
    kept = filter(s -> !(s.mf == 3 && s.mt in 203:207), material.sections)

    # Insert new gas sections in MT order, then resort by (MF, MT)
    out_sections = PENDFSection[]
    append!(out_sections, kept)
    for mt in sort!(collect(keys(new_sections)))
        push!(out_sections, PENDFSection(3, mt, new_sections[mt]))
    end
    sort!(out_sections, by = s -> (s.mf, s.mt))

    # MF1/MT451 directory: bump NXC, splice (MF=3, MT, NC=(NP+2)÷3, MOD=1) entries
    # in MT-sorted order. Fortran gaspr.f90:961-1005, 1019.
    new_mf1 = _gaspr_update_mf1_directory(material.mf1_lines, mat, np_per_mt)

    PENDFMaterial(mat, new_mf1, out_sections)
end

# =========================================================================
# Helpers
# =========================================================================

# Linear interpolation (LinLin) of (xs_src) on grid (e_src) at point e.
# Returns 0 outside the source range. Matches Fortran gety1/terpa default
# for PENDF after broadr (everything is INT=2 by then).
function _interp_linlin_at(e_src::AbstractVector{<:Real},
                            xs_src::AbstractVector{<:Real}, e::Real)
    n = length(e_src)
    n == 0 && return 0.0
    e < e_src[1] && return 0.0
    e > e_src[end] && return 0.0
    e == e_src[1] && return Float64(xs_src[1])
    e == e_src[end] && return Float64(xs_src[end])
    # Binary search for interval [i, i+1] with e_src[i] <= e < e_src[i+1]
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

# Pull ZA, AWR from MF1/MT451 head line (line 1 of mf1_lines).
function _gaspr_za_awr(mf1_lines::Vector{String})
    isempty(mf1_lines) && return (0.0, 0.0)
    p = rpad(mf1_lines[1], 80)
    (parse_endf_float(p[1:11]), parse_endf_float(p[12:22]))
end

# Build an MF3 section (HEAD + TAB1 header + interp + data + SEND) for a new
# gas reaction. ZA/AWR carried from MF1/451 head; QM=QI=0; LR=0; INT=2.
# Matches the layout written by Fortran gaspr.f90:1064-1109.
function _gaspr_build_mf3_section(es::Vector{Float64}, xss::Vector{Float64},
                                    mat::Int, mt::Int, za::Float64, awr::Float64)
    np = length(es)
    lines = String[]
    trailer = @sprintf("%4d%2d%3d", mat, 3, mt)

    # HEAD: ZA, AWR, 0, 0, 0, 0
    push!(lines, format_endf_float(za) * format_endf_float(awr) *
                 lpad("0", 11) * lpad("0", 11) * lpad("0", 11) * lpad("0", 11) *
                 trailer * "    0")

    # TAB1 header: QM=0, QI=0, 0, LR=0, NR=1, NP
    push!(lines, format_endf_float(0.0) * format_endf_float(0.0) *
                 lpad("0", 11) * lpad("0", 11) *
                 lpad("1", 11) * lpad(string(np), 11) *
                 trailer * "    0")

    # Interp: NBT=NP, INT=2 (linear-linear)
    push!(lines, lpad(string(np), 11) * lpad("2", 11) * " "^44 *
                 trailer * "    0")

    # Data: 3 (E,XS) pairs per line (SEND added by writer)
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

# Splice new gas-production directory entries into MF1/MT451 and bump NXC.
# `np_per_mt`: Dict{MT => NP} (number of (E,XS) points written for that MT).
# Fortran gaspr.f90:889-1028 — set NC = (NP+2)÷3, MOD=1, insert in MT order.
function _gaspr_update_mf1_directory(mf1_lines::Vector{String}, mat::Int,
                                       np_per_mt::Dict{Int, Int})
    isempty(np_per_mt) && return mf1_lines

    # Locate the line carrying NXC (the 3rd CONT line of MF1/MT451 in
    # ENDF-6 = TEMP/ERR/LDRV/0/NWD/NXC). For PENDF this is line 4 of mf1_lines
    # (line 1=HEAD ZA/AWR, line 2=blank CONT, line 3=AWI/EMAX CONT [iverf>=6],
    # line 4=TEMP/ERR CONT). Detect by finding the line whose seq position
    # matches and which precedes the NWD descriptive lines.
    out = copy(mf1_lines)
    nxc_line_idx = _gaspr_find_nxc_line(out)
    nxc_line_idx == 0 && error("gaspr: cannot locate MF1/MT451 NXC line in MAT=$mat")

    # Parse current NXC, NWD from that CONT line (cols 45-55=N1=NWD, 56-66=N2=NXC).
    cont = rpad(out[nxc_line_idx], 80)
    nwd = _parse_int(cont[45:55])
    nxc = _parse_int(cont[56:66])

    # Existing directory entries follow the NWD description lines.
    # Each is a CONT line with cols 23-33=MF, 34-44=MT, 45-55=NC, 56-66=MOD.
    dir_start = nxc_line_idx + nwd + 1
    dir_end   = dir_start + nxc - 1
    dir_end <= length(out) ||
        error("gaspr: MF1/MT451 directory truncated for MAT=$mat (need lines " *
              "$dir_start..$dir_end, have $(length(out)))")

    # Strip any pre-existing MF=3, MT=203..207 entries (Fortran "delete old gas")
    new_entries = String[]
    for li in dir_start:dir_end
        l = rpad(out[li], 80)
        emf = _parse_int(l[23:33])
        emt = _parse_int(l[34:44])
        (emf == 3 && emt in 203:207) && continue
        push!(new_entries, l)
    end

    # Build new entries for the added gas MTs. NC = (NP+2)÷3, MOD = 1.
    for mt in sort!(collect(keys(np_per_mt)))
        np = np_per_mt[mt]
        nc = (np + 2) ÷ 3
        push!(new_entries,
            " "^22 *
            lpad("3", 11) * lpad(string(mt), 11) *
            lpad(string(nc), 11) * lpad("1", 11) *
            @sprintf("%4d%2d%3d", mat, 1, 451) * "    0")
    end

    # Re-sort directory by (MF, MT) so MT=203..207 land in MF=3 cluster
    # (between MT=111 and MT=301-style entries).
    sort!(new_entries, by = l -> begin
        p = rpad(l, 80)
        (_parse_int(p[23:33]), _parse_int(p[34:44]))
    end)

    # Update NXC field on the CONT line. Format-preserving rewrite of cols 56-66.
    new_nxc = length(new_entries)
    out[nxc_line_idx] = string(cont[1:55], lpad(string(new_nxc), 11), cont[67:80])

    # Reassemble: pre-directory + description + new directory + post-directory
    pre  = out[1:nxc_line_idx + nwd]
    post = out[dir_end + 1:end]
    vcat(pre, new_entries, post)
end

# Locate the MF1/MT451 CONT line that holds (TEMP, ERR, LDRV, 0, NWD, NXC).
# Two PENDF dialects exist: Fortran-NJOY uses HEAD + 2 extra CONTs (blank,
# AWI/EMAX) before TEMP/ERR (so NXC is on line 4); Julia heatr's writer
# emits HEAD + TEMP/ERR directly (NXC on line 2). Detect by structural fit:
# the NXC line is the unique CONT where `idx + NWD + NXC == n` — everything
# after it is exactly NWD description lines + NXC directory entries.
function _gaspr_find_nxc_line(mf1_lines::Vector{String})
    n = length(mf1_lines)
    for idx in 2:min(n, 6)
        p = rpad(mf1_lines[idx], 80)
        nwd = _parse_int(p[45:55])
        nxc = _parse_int(p[56:66])
        nwd >= 0 && nxc >= 1 || continue
        idx + nwd + nxc == n && return idx
    end
    0
end

