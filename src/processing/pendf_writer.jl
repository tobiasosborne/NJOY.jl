# PENDF writer -- streaming output of pointwise cross sections
#
# PROPOSAL B DESIGN: The writer is designed for streaming I/O:
#   - Each section is written independently, never buffering the whole file
#   - The `write_pendf` function works with `PointwiseMaterial` directly
#   - MF2 is copied verbatim from the original ENDF source if available
#   - Provides both PointwiseMaterial and NamedTuple interfaces
#
# Correspondence to NJOY2016 reconr.f90:
#   recout (4984-5441) -> write_pendf, write_pendf_file

using Printf

# ==========================================================================
# Main entry points
# ==========================================================================

"""
    write_pendf(io::IO, material::PointwiseMaterial;
                endf_source=nothing, temperature=0.0, err=0.001,
                description=String[])

Write a PENDF tape for `material` to the output stream `io`.

The output is in standard ENDF-6 format with:
- MF1/MT451: descriptive data with dictionary
- MF2/MT151: resonance parameters (copied from source or minimal placeholder)
- MF3/MT{N}: pointwise cross sections with lin-lin interpolation
- Proper SEND/FEND/MEND/TEND delimiters
"""
function write_pendf(io::IO, material::PointwiseMaterial;
                     endf_source::Union{IO, Nothing} = nothing,
                     temperature::Float64 = 0.0,
                     err::Float64 = 0.001,
                     description::Vector{String} = String[],
                     za::Float64 = 0.0,
                     awr::Float64 = 0.0)
    mat = Int(material.mat)
    n_pts = length(material.energies)
    ns = Ref(1)  # line sequence counter

    # TPID record
    _write_tpid_line(io, "PENDF tape produced by NJOY.jl", mat)

    # MF1/MT451
    _write_mf1_section(io, mat, material, temperature, err, description, ns)

    # MF2/MT151
    if endf_source !== nothing
        _copy_mf2_from_source(io, endf_source, mat)
    else
        _write_minimal_mf2(io, mat, ns)
    end

    # MF3 for each reaction -- continuous sequence numbering across sections
    ns[] = 1
    for (col, mt) in enumerate(material.mt_list)
        _write_mf3_from_matrix(io, mat, mt, material.energies,
                               material.cross_sections, col, ns;
                               za=za, awr=awr)
    end
    _write_fend(io, mat; ns=ns)

    # MEND + TEND
    _write_record_sep(io, 0, 0, 0)
    _write_record_sep(io, -1, 0, 0)
end

"""
    write_pendf(io::IO, result::NamedTuple; mat=0, label="",
                err=0.001, tempr=0.0)

Write PENDF from the legacy reconr() NamedTuple result format.
"""
function write_pendf(io::IO, result::NamedTuple;
                     mat::Integer = 0,
                     label::AbstractString = "reconstructed data",
                     err::Float64 = 0.001,
                     tempr::Float64 = 0.0)
    mf2 = result.mf2
    actual_mat = mat > 0 ? Int32(mat) : Int32(max(1, round(Int, mf2.ZA / 10)))
    ns = Ref(1)

    # Determine reactions
    reactions = _collect_reactions(result)

    # TPID
    _write_tpid_line(io, label, Int(actual_mat))

    # MF1/MT451
    _write_legacy_mf1(io, mf2, actual_mat, err, tempr, length(result.energies), reactions, ns)

    # MF2/MT151
    _write_legacy_mf2(io, mf2, actual_mat, ns)

    # MF3 sections
    _write_legacy_mf3(io, result, actual_mat, reactions, ns)

    # MEND and TEND (NS=0 matching Fortran)
    _write_fend_zero(io, 0)   # MEND: MAT=0, MF=0, MT=0, NS=0
    blanks = repeat(" ", 66)
    @printf(io, "%s%4d%2d%3d%5d\n", blanks, -1, 0, 0, 0)  # TEND
end

"""
    write_pendf_file(filename, data; kwargs...)

Write PENDF to a file. `data` can be `PointwiseMaterial` or `NamedTuple`.
"""
function write_pendf_file(filename::AbstractString, data; kwargs...)
    open(filename, "w") do io
        write_pendf(io, data; kwargs...)
    end
end

# ==========================================================================
# Low-level ENDF formatting
# ==========================================================================

function _write_tpid_line(io::IO, text::AbstractString, mat::Int)
    t = rpad(text, 66)[1:66]
    # TPID always uses MAT=1, MF=0, MT=0, NS=0 per NJOY convention
    @printf(io, "%s%4d%2d%3d%5d\n", t, 1, 0, 0, 0)
end

function _write_record_sep(io::IO, mat::Int, mf::Int, mt::Int;
                           ns::Union{Ref{Int}, Nothing}=nothing)
    blanks = repeat(" ", 66)
    seq = ns === nothing ? 99999 : ns[]
    @printf(io, "%s%4d%2d%3d%5d\n", blanks, mat, mf, mt, seq)
    if ns !== nothing
        ns[] += 1
    end
end

function _write_fend(io::IO, mat::Int; ns::Union{Ref{Int}, Nothing}=nothing)
    # FEND only: each MT section writes its own SEND; no extra SEND here
    _write_record_sep(io, mat, 0, 0; ns=ns)
end

"""Write SEND record with NS=99999 (matching Fortran convention)."""
function _write_send(io::IO, mat::Int, mf::Int)
    blanks = repeat(" ", 66)
    @printf(io, "%s%4d%2d%3d%5d\n", blanks, mat, mf, 0, 99999)
end

"""Write FEND record with NS=0 (matching Fortran convention)."""
function _write_fend_zero(io::IO, mat::Int)
    blanks = repeat(" ", 66)
    @printf(io, "%s%4d%2d%3d%5d\n", blanks, mat, 0, 0, 0)
end

function _write_cont_line(io::IO, c1, c2, l1, l2, n1, n2,
                          mat::Int, mf::Int, mt::Int, ns::Ref{Int})
    s1 = format_endf_float(Float64(c1))
    s2 = format_endf_float(Float64(c2))
    @printf(io, "%s%s%11d%11d%11d%11d%4d%2d%3d%5d\n",
            s1, s2, l1, l2, n1, n2, mat, mf, mt, ns[])
    ns[] += 1
end

function _write_data_values(io::IO, vals, mat::Int, mf::Int, mt::Int, ns::Ref{Int};
                            pair_data::Bool=false)
    n = length(vals)
    i = 1
    while i <= n
        for j in 1:6
            if i <= n
                # For (E, XS) pair data: odd positions are energies (may use 9-sigfig
                # extended format), even positions are XS (always 7-sigfig scientific).
                # Matching Fortran a11 behavior where only energies from round_sigfig(x,9)
                # have genuine 9-sigfig precision.
                ext = !pair_data || isodd(i)
                write(io, format_endf_float(Float64(vals[i]); extended=ext))
                i += 1
            else
                write(io, "           ")
            end
        end
        @printf(io, "%4d%2d%3d%5d\n", mat, mf, mt, ns[])
        ns[] += 1
    end
end

# ==========================================================================
# MF3 writer for PointwiseMaterial
# ==========================================================================

function _write_mf3_from_matrix(io::IO, mat::Int, mt::Int,
                                 energies::Vector{Float64},
                                 xs::Matrix{Float64}, col::Int,
                                 ns::Ref{Int};
                                 za::Float64=0.0, awr::Float64=0.0)
    n = length(energies)
    n >= 2 || return

    # HEAD: ZA, AWR, 0, LR=99, 0, 0 (matching reference tape format)
    _write_cont_line(io, za, awr, 0, 99, 0, 0, mat, 3, mt, ns)

    # TAB1 header: QM=0, QI=0, 0, 0, NR=1, NP=n
    _write_cont_line(io, 0.0, 0.0, 0, 0, 1, n, mat, 3, mt, ns)

    # Interpolation table: NBT=n, INT=2
    vals = Float64[Float64(n), 2.0]
    _write_data_values(io, vals, mat, 3, mt, ns)

    # Data pairs
    data = Float64[]
    sizehint!(data, 2 * n)
    for i in 1:n
        push!(data, energies[i])
        push!(data, xs[i, col])
    end
    _write_data_values(io, data, mat, 3, mt, ns; pair_data=true)

    # SEND with continuous sequence numbering
    _write_record_sep(io, mat, 3, 0; ns=ns)
end

# ==========================================================================
# MF1 section writer for PointwiseMaterial
# ==========================================================================

function _write_mf1_section(io::IO, mat::Int, material::PointwiseMaterial,
                            temperature::Float64, err::Float64,
                            description::Vector{String}, ns::Ref{Int})
    ns[] = 1
    nxc = length(material.mt_list) + 2  # MF1 + MF2 + MF3 entries
    nwd = length(description)

    # HEAD
    _write_cont_line(io, 0.0, 0.0, 0, 0, 0, 0, mat, 1, 451, ns)
    # Control
    _write_cont_line(io, temperature, err, 0, 0, nwd, nxc, mat, 1, 451, ns)

    # Description cards
    for desc in description
        padded = rpad(desc, 66)[1:66]
        @printf(io, "%s%4d%2d%3d%5d\n", padded, mat, 1, 451, ns[])
        ns[] += 1
    end

    # Dictionary entries
    n_pts = length(material.energies)
    nc_per_section = 3 + cld(n_pts, 3)
    _write_cont_line(io, 0.0, 0.0, 1, 451, 3, 0, mat, 1, 451, ns)
    _write_cont_line(io, 0.0, 0.0, 2, 151, 4, 0, mat, 1, 451, ns)
    for mt in material.mt_list
        _write_cont_line(io, 0.0, 0.0, 3, mt, nc_per_section, 0, mat, 1, 451, ns)
    end

    # SEND + FEND
    _write_record_sep(io, mat, 1, 0)
    _write_record_sep(io, mat, 0, 0)
end

# ==========================================================================
# MF2 handling
# ==========================================================================

function _write_minimal_mf2(io::IO, mat::Int, ns::Ref{Int})
    ns[] = 1
    _write_cont_line(io, 0.0, 0.0, 0, 0, 1, 0, mat, 2, 151, ns)
    _write_cont_line(io, 0.0, 0.0, 0, 0, 0, 0, mat, 2, 151, ns)
    _write_record_sep(io, mat, 2, 0)
    _write_record_sep(io, mat, 0, 0)
end

function _copy_mf2_from_source(io_out::IO, io_src::IO, mat::Int)
    seekstart(io_src)
    copying = false
    while !eof(io_src)
        line = readline(io_src)
        p = rpad(line, 80)
        mf = _pendf_parse_int(p[71:72])
        mat_line = _pendf_parse_int(p[67:70])

        if mat_line == mat && mf == 2
            copying = true
        end

        if copying
            println(io_out, rstrip(p))
            if mf == 0 && copying
                break
            end
        end
    end

    if !copying
        ns = Ref(1)
        _write_minimal_mf2(io_out, mat, ns)
    end
end

# ==========================================================================
# Legacy interface (NamedTuple result)
# ==========================================================================

function _collect_reactions(result)
    # Collect reactions in ENDF file order (matching Fortran emerge).
    # MT=1 (total) always first. Then non-redundant MTs from MF3 in file order.
    # MT=4 is output as a computed redundant sum of MT=51-91.
    reactions = Tuple{Int32, String}[]
    push!(reactions, (Int32(1), "total"))

    has_inelastic = any(s -> Int(s.mt) >= 51 && Int(s.mt) <= 91, result.mf3_sections)

    for sec in result.mf3_sections
        mt = Int(sec.mt)
        # Skip truly redundant MTs (never output by Fortran)
        (mt == 1 || mt == 3 || mt == 101 || mt == 120 ||
         mt == 151 || mt == 27 || (mt >= 251 && mt <= 300 && mt != 261)) && continue
        # MT=4: output as redundant sum if inelastic levels present
        if mt == 4
            has_inelastic && push!(reactions, (Int32(4), "inelastic"))
            continue
        end
        push!(reactions, (Int32(mt), "MT$mt"))
    end
    return reactions
end

function _write_legacy_mf1(io::IO, mf2::MF2Data, mat::Int32, err, tempr,
                            n_pts, reactions, ns::Ref{Int})
    ns[] = 1
    nxc = 2 + length(reactions)

    _write_cont_line(io, mf2.ZA, mf2.AWR, 2, 0, 0, 0,
                     Int(mat), 1, 451, ns)
    _write_cont_line(io, tempr, err, 0, 0, 0, nxc,
                     Int(mat), 1, 451, ns)

    nc_per = 3 + cld(n_pts, 3)
    _write_cont_line(io, 0.0, 0.0, 1, 451, 3, 0, Int(mat), 1, 451, ns)
    _write_cont_line(io, 0.0, 0.0, 2, 151, 4, 0, Int(mat), 1, 451, ns)
    for (mt, _) in reactions
        _write_cont_line(io, 0.0, 0.0, 3, Int(mt), nc_per, 0,
                         Int(mat), 1, 451, ns)
    end

    _write_send(io, Int(mat), 1)
    _write_fend_zero(io, Int(mat))
end

function _write_legacy_mf2(io::IO, mf2::MF2Data, mat::Int32, ns::Ref{Int})
    ns[] = 1
    nis = length(mf2.isotopes)
    _write_cont_line(io, mf2.ZA, mf2.AWR, 0, 0, nis, 0, Int(mat), 2, 151, ns)

    for iso in mf2.isotopes
        _write_cont_line(io, iso.ZAI, iso.ABN, 0, 0, 1, 0, Int(mat), 2, 151, ns)
        el = length(iso.ranges) > 0 ? iso.ranges[1].EL : 1.0e-5
        eh = length(iso.ranges) > 0 ? iso.ranges[end].EH : 2.0e7
        _write_cont_line(io, el, eh, 0, 0, 0, 0, Int(mat), 2, 151, ns)
        spi = 1.0; ap = 0.0
        if length(iso.ranges) > 0
            p = iso.ranges[1].parameters
            if hasproperty(p, :SPI); spi = p.SPI; end
            if hasproperty(p, :AP); ap = p.AP; end
        end
        _write_cont_line(io, spi, ap, 0, 0, 0, 0, Int(mat), 2, 151, ns)
    end

    _write_send(io, Int(mat), 2)
    _write_fend_zero(io, Int(mat))
end

function _write_legacy_mf3(io::IO, result, mat::Int32, reactions, ns::Ref{Int})
    energies = result.energies
    awr = result.mf2.AWR
    za = result.mf2.ZA

    for (mt, _) in reactions
        ns[] = 1  # Restart sequence per section (matching Fortran)

        # Get xs data and determine threshold
        sec_data = _get_legacy_section(result, Int(mt))
        sec_data === nothing && continue

        sec_e, sec_xs, qm, qi = sec_data

        # HEAD: ZA, AWR, 0, LR, 0, 0
        # LR=99 only for MT=1 (total), matching Fortran emerge
        lr = mt == 1 ? 99 : 0
        _write_cont_line(io, za, awr, 0, lr, 0, 0,
                         Int(mat), 3, Int(mt), ns)

        # TAB1: QM, QI, 0, 0, NR=1, NP
        np = length(sec_e)
        _write_cont_line(io, qm, qi, 0, 0, 1, np,
                         Int(mat), 3, Int(mt), ns)

        # Interpolation table: NBT, INT as integers (Fortran format)
        @printf(io, "%11d%11d%44s%4d%2d%3d%5d\n",
                np, 2, "", Int(mat), 3, Int(mt), ns[])
        ns[] += 1

        # Data pairs
        data = Float64[]
        sizehint!(data, 2 * np)
        for i in 1:np
            push!(data, sec_e[i])
            push!(data, sec_xs[i])
        end
        _write_data_values(io, data, Int(mat), 3, Int(mt), ns; pair_data=true)

        # SEND (NS=99999 per ENDF convention)
        _write_send(io, Int(mat), 3)
    end

    # FEND (MAT, MF=0, MT=0, NS=0)
    _write_fend_zero(io, Int(mat))
end

"""Get energies, XS, QM, QI for a given MT, handling thresholds."""
function _get_legacy_section(result, mt::Int)
    energies = result.energies
    awr = result.mf2.AWR

    # Look up QM/QI from the MF3 section data
    qm, qi = 0.0, 0.0
    for sec in result.mf3_sections
        if Int(sec.mt) == mt
            qm, qi = sec.QM, sec.QI
            break
        end
    end

    if mt == 1
        return energies, result.total, 0.0, 0.0
    elseif mt == 4
        # Redundant sum: MT=4 = sum(MT=51-91)
        # Use the grid from the first inelastic level
        first_inel = findfirst(s -> Int(s.mt) >= 51 && Int(s.mt) <= 91, result.mf3_sections)
        first_inel === nothing && return nothing
        sec_first = result.mf3_sections[first_inel]
        qi_4 = sec_first.QI  # Use QI from first inelastic level
        thrx_4 = awr > 0.0 ? -qi_4 * (awr + 1) / awr : -qi_4
        thrxx_4 = round_sigfig(thrx_4, 7, +1)
        # Build grid: all energies at or above threshold
        sec_e = Float64[]
        sec_xs = Float64[]
        for e in energies
            thrxx_4 > 0.0 && thrxx_4 - e > 1.0e-9 * thrxx_4 && continue
            # Sum all MT=51-91 at this energy
            total_inel = 0.0
            for sec in result.mf3_sections
                smt = Int(sec.mt)
                (smt < 51 || smt > 91) && continue
                bg = interpolate(sec.tab, e)
                # Threshold-adjusted interpolation for each level
                s_qi = sec.QI
                if s_qi < 0.0
                    s_thrx = awr > 0.0 ? -s_qi * (awr + 1) / awr : -s_qi
                    s_thrxx = round_sigfig(s_thrx, 7, +1)
                    if s_thrxx > 0.0 && e >= s_thrxx && length(sec.tab.x) >= 2 &&
                       e < sec.tab.x[2] && sec.tab.x[1] < s_thrxx
                        bg = _threshold_interp(sec.tab, e, s_thrxx)
                    end
                    if s_thrxx > 0.0 && abs(s_thrxx - e) < 1.0e-9 * s_thrxx
                        bg = 0.0
                    end
                end
                total_inel += round_sigfig(bg, 7)
            end
            push!(sec_e, e)
            push!(sec_xs, round_sigfig(total_inel, 7))
        end
        isempty(sec_e) && return nothing
        # Pseudo-threshold: skip leading zero values (matching Fortran
        # recout ith tracking — start at the point before first nonzero)
        first_nz = findfirst(>(0.0), sec_xs)
        if first_nz !== nothing && first_nz > 1
            sec_e = sec_e[first_nz-1:end]
            sec_xs = sec_xs[first_nz-1:end]
        end
        isempty(sec_e) && return nothing
        return sec_e, sec_xs, 0.0, qi_4
    elseif mt == 2
        return energies, result.elastic, qm, qi
    elseif mt == 18 || mt == 19
        return energies, result.fission, qm, qi
    elseif mt == 102
        return energies, result.capture, qm, qi
    end

    # Non-primary MT: find the MF3 section
    for sec in result.mf3_sections
        Int(sec.mt) != mt && continue

        qm = sec.QM
        qi = sec.QI

        # Compute threshold for reactions with Q < 0
        thrxx = 0.0
        if qi < 0.0
            thrx = awr > 0.0 ? -qi * (awr + 1) / awr : -qi
            thrxx = round_sigfig(thrx, 7, +1)
        end

        # Build the section output matching Fortran emerge:
        # 1. Skip points below threshold (line 4792)
        # 2. Pseudo-threshold skip: skip zero-XS panels (lines 1973-1976)
        # 3. Set xs=0 at threshold energy (line 4795)
        sec_e = Float64[]
        sec_xs = Float64[]
        found_nonzero = false
        for (idx, e) in enumerate(energies)
            if thrxx > 0.0 && thrxx - e > 1.0e-9 * thrxx
                continue  # below threshold
            end
            # Near threshold, modify first breakpoint to (thrxx, 0) and interpolate
            # using the MF3 interpolation law (Fortran emerge:1936 + gety1)
            bg = if thrxx > 0.0 && e >= thrxx && e < sec.tab.x[2] && sec.tab.x[1] < thrxx
                _threshold_interp(sec.tab, e, thrxx)
            else
                interpolate(sec.tab, e)
            end
            bg = round_sigfig(bg, 7)

            # Pseudo-threshold skip: if current xs ≈ 0, check if we've
            # reached non-zero data yet. Skip leading zero-XS panels.
            if !found_nonzero
                if abs(bg) < 1.0e-30
                    # Look ahead: is the NEXT point also zero?
                    next_idx = idx + 1
                    while next_idx <= length(energies)
                        next_e = energies[next_idx]
                        if thrxx > 0.0 && thrxx - next_e > 1.0e-9 * thrxx
                            next_idx += 1; continue
                        end
                        break
                    end
                    if next_idx <= length(energies)
                        next_bg = round_sigfig(interpolate(sec.tab, energies[next_idx]), 7)
                        if abs(next_bg) < 1.0e-30
                            continue  # both zero → skip
                        end
                    end
                end
                found_nonzero = true
            end

            # Set xs=0 at threshold (line 4795)
            if thrxx > 0.0 && abs(thrxx - e) < 1.0e-9 * thrxx
                bg = 0.0
            end

            push!(sec_e, e)
            push!(sec_xs, bg)
        end

        isempty(sec_e) && return nothing
        return sec_e, sec_xs, qm, qi
    end
    return nothing
end

function _pendf_parse_int(s::AbstractString)
    t = strip(s)
    isempty(t) && return 0
    return parse(Int, t)
end

"""
    _threshold_interp(tab, e, thrxx)

Threshold-adjusted interpolation: temporarily modify the first breakpoint
to (thrxx, 0.0) and interpolate using the MF3's own interpolation law.
Matches Fortran emerge line 1936 (modify ex(1)) + gety1.
"""
function _threshold_interp(tab::TabulatedFunction, e::Float64, thrxx::Float64)
    saved_x, saved_y = tab.x[1], tab.y[1]
    tab.x[1] = thrxx
    tab.y[1] = 0.0
    result = interpolate(tab, e)
    tab.x[1] = saved_x
    tab.y[1] = saved_y
    return result
end
