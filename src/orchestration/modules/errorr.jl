# errorr module runner -- Covariance processing from ENDF MF33/MF31
#
# Matches Fortran errorr.f90: read covariance data from ENDF, compute
# group-averaged cross sections, expand covariances onto union grid,
# write output tape in GENDF-like format.

"""
    errorr_module(tapes::TapeManager, params::ErrorrParams)

Run ERRORR: read ENDF covariance data, compute multigroup covariances,
write output tape. Matches Fortran errorr.f90 interface.

Input tapes: nendf (ENDF evaluation), npend (PENDF, optional), ngout (GENDF, optional)
Output tapes: nout (covariance output)
"""
function errorr_module(tapes::TapeManager, params::ErrorrParams)
    @info "errorr: MAT=$(params.mat) mfcov=$(params.mfcov) → tape $(params.nout)"

    endf_path = resolve(tapes, params.nendf)
    nout_path = resolve(tapes, params.nout)

    # Determine group structure
    egn = _errorr_group_structure(params, endf_path)
    ngn = length(egn) - 1
    @info "errorr: $ngn groups"

    # Read ZA and AWR from ENDF
    za, awr = _read_za_awr(endf_path, params.mat)

    # Compute group-averaged cross sections
    group_xs = Dict{Int, Vector{Float64}}()
    if params.npend > 0 && params.iread == 0
        pendf_path = resolve(tapes, params.npend)
        group_xs = _errorr_group_average(pendf_path, params.mat, egn, params.iwt)
    elseif params.ngout != 0 && params.iread != 0
        gendf_path = resolve(tapes, abs(params.ngout))
        group_xs = _errorr_read_gendf_xs(gendf_path, params.mat, params.mfcov)
    end

    # Read and process covariance data
    cov_results = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    reaction_mts = Int[]

    if params.nin > 0
        # Read covariance from input tape (second errorr call pattern)
        nin_path = resolve(tapes, params.nin)
        cov_results, reaction_mts, group_xs, egn =
            _errorr_read_input_cov(nin_path, params.mat, params.mfcov)
        ngn = length(egn) - 1
    else
        # Read covariance from ENDF tape
        cov_results, reaction_mts = _errorr_compute_cov(endf_path, params.mat,
                                                         params.mfcov, egn)
    end

    # Write output tape
    open(nout_path, "w") do io
        _write_errorr_tape(io, params.mat, za, awr, egn, group_xs,
                          cov_results, reaction_mts, params.mfcov)
    end
    register!(tapes, params.nout, nout_path)

    lines = countlines(nout_path)
    @info "errorr: wrote $nout_path ($lines lines, $ngn groups, $(length(reaction_mts)) reactions)"
    nothing
end

# =========================================================================
# Group structure determination
# =========================================================================

function _errorr_group_structure(params::ErrorrParams, endf_path::String)
    if params.ign < 0 && !isempty(params.user_egn)
        # User-defined structure — union with covariance grid
        return _errorr_union_grid(endf_path, params.mat, params.mfcov, params.user_egn)
    elseif params.ign == 3
        return collect(Float64, LANL_30)
    else
        # Default: build from covariance file
        return _errorr_union_grid(endf_path, params.mat,
                                  params.mfcov > 0 ? params.mfcov : 33, Float64[])
    end
end

"""Build union energy grid from MF33/MF31 covariance data + user boundaries."""
function _errorr_union_grid(endf_path::String, mat::Int, mfcov::Int,
                            user_egn::Vector{Float64})
    pts = Set{Float64}()
    for e in user_egn; push!(pts, e); end

    open(endf_path) do io
        seekstart(io)
        # Scan for all MFcov sections and extract energy grids
        while !eof(io)
            line = readline(io); length(line) < 75 && continue
            p = rpad(line, 80)
            mf_val = _parse_int(p[71:72]); mt_val = _parse_int(p[73:75])
            mat_val = _parse_int(p[67:70])
            if mat_val == mat && mf_val == mfcov && mt_val > 0
                # Found a covariance section — read its energy grids
                _extract_cov_energies!(pts, io, mat, mfcov, mt_val)
            end
        end
    end

    grid = sort!(collect(pts))
    # Clamp to user boundaries if provided
    if !isempty(user_egn)
        elow, ehigh = extrema(user_egn)
        filter!(e -> e >= elow && e <= ehigh, grid)
        elow in grid || pushfirst!(grid, elow)
        ehigh in grid || push!(grid, ehigh)
    end
    grid
end

"""Extract energy values from one MFcov section's NI-type sub-subsections."""
function _extract_cov_energies!(pts::Set{Float64}, io::IO, mat::Int,
                                 mfcov::Int, mt::Integer)
    seekstart(io)
    find_section(io, mfcov, mt; target_mat=mat) || return
    head = read_cont(io); nl = Int(head.N2)
    for _ in 1:nl
        sh = read_cont(io); nc = Int(sh.N1); ni = Int(sh.N2)
        for _ in 1:nc; read_list(io); end  # skip NC-type
        for _ in 1:ni
            lst = read_list(io)
            lb = Int(lst.L2); ne = Int(lst.N2)
            if lb in (0,1,2,3,4)
                nk = div(Int(lst.N1), 2)
                for k in 1:nk; push!(pts, lst.data[2k-1]); push!(pts, lst.data[2k-1]); end
                # Add upper boundary of last interval
                nk >= 1 && push!(pts, lst.data[2nk-1])
            elseif lb in (5,6)
                for k in 1:ne; push!(pts, lst.data[k]); end
            end
        end
    end
end

# =========================================================================
# Group-averaged cross sections from PENDF
# =========================================================================

function _errorr_group_average(pendf_path::String, mat::Int,
                               egn::Vector{Float64}, iwt::Int)
    result = Dict{Int, Vector{Float64}}()
    tape = read_pendf(pendf_path)
    mf3 = extract_mf3_all(tape, mat)
    ngn = length(egn) - 1
    wfn = iwt == 3 ? inv_e_weight : constant_weight

    for (mt, (energies, xs)) in mf3
        mt in (1, 451) && continue  # skip total and directory
        flux = group_integrate(energies, [wfn(e) for e in energies], egn)
        sig_flux = group_integrate(energies, [xs[i] * wfn(energies[i]) for i in eachindex(energies)], egn)
        avg = [flux[g] > 0 ? sig_flux[g] / flux[g] : 0.0 for g in 1:ngn]
        result[mt] = avg
    end
    result
end

# =========================================================================
# Read group-averaged XS from GENDF tape
# =========================================================================

function _errorr_read_gendf_xs(gendf_path::String, mat::Int, mfcov::Int)
    result = Dict{Int, Vector{Float64}}()
    isfile(gendf_path) || return result

    lines = readlines(gendf_path)
    idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = _parse_int(p[71:72]); mt = _parse_int(p[73:75])
        mat_val = _parse_int(p[67:70])
        if mat_val == mat && mf == 3 && mt > 0
            # Read HEAD to get number of groups
            ng = _parse_int(p[56:66])
            idx += 1
            # Read the group-averaged values
            vals = Float64[]
            while idx <= length(lines) && length(vals) < ng
                length(lines[idx]) < 66 && (idx += 1; continue)
                p2 = rpad(lines[idx], 80)
                mf2 = _parse_int(p2[71:72]); mt2 = _parse_int(p2[73:75])
                (mf2 != 3 || mt2 != mt) && break
                for col in 0:5
                    s = p2[1+11*col:11+11*col]
                    v = tryparse(Float64, strip(replace(s, r"([0-9])([+-])(\d)" => s"\1e\2\3")))
                    v !== nothing && push!(vals, v)
                    length(vals) >= ng && break
                end
                idx += 1
            end
            result[mt] = vals
        else
            idx += 1
        end
    end
    result
end

# =========================================================================
# Covariance computation
# =========================================================================

function _errorr_compute_cov(endf_path::String, mat::Int, mfcov::Int,
                              egn::Vector{Float64})
    cov_results = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    reaction_mts = Int[]

    open(endf_path) do io
        seekstart(io)
        # Find all MTs in MFcov
        avail_mts = Int[]
        while !eof(io)
            line = readline(io); length(line) < 75 && continue
            p = rpad(line, 80)
            mf = _parse_int(p[71:72]); mt = _parse_int(p[73:75])
            mat_val = _parse_int(p[67:70])
            if mat_val == mat && mf == mfcov && mt > 0 && !(mt in avail_mts)
                push!(avail_mts, mt)
            end
        end

        for mt in avail_mts
            try
                cd = read_mf33(io, mat, mt)
                if !isempty(cd.ni_subsections)
                    cm = multigroup_covariance(cd.ni_subsections, egn;
                                               is_relative=true)
                    mt1 = cd.ni_subsections[1].mt1
                    mt2 = cd.ni_subsections[1].mt2
                    cov_results[(mt, mt1)] = cm.matrix
                    mt in reaction_mts || push!(reaction_mts, mt)
                    mt1 != mt && !(mt1 in reaction_mts) && push!(reaction_mts, mt1)
                end
            catch e
                @warn "errorr: skipping MF$mfcov/MT$mt: $e"
            end
        end
    end

    sort!(reaction_mts)
    cov_results, reaction_mts
end

function _errorr_read_input_cov(nin_path::String, mat::Int, mfcov::Int)
    # Read covariance data from a previous errorr output tape
    # For the second errorr call, this reads from tape23 (first errorr output)
    cov_results = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    reaction_mts = Int[]
    group_xs = Dict{Int, Vector{Float64}}()
    egn = Float64[]

    lines = readlines(nin_path)
    idx = 1

    # Parse MF1/MT451 to get group boundaries
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = _parse_int(p[71:72]); mt = _parse_int(p[73:75])
        if mf == 1 && mt == 451
            idx += 1  # Move to next line (N1, N2, NW fields)
            length(lines[idx]) < 66 && (idx += 1; continue)
            p2 = rpad(lines[idx], 80)
            ngn = _parse_int(p2[23:33])  # N1 = number of groups
            nw = _parse_int(p2[45:55])   # N2 = number of words in directory
            idx += 1
            # Read group boundaries
            while idx <= length(lines) && length(egn) < ngn + 1
                length(lines[idx]) < 66 && (idx += 1; continue)
                p3 = rpad(lines[idx], 80)
                mt3 = _parse_int(p3[73:75])
                mt3 != 451 && break
                for col in 0:5
                    s = p3[1+11*col:11+11*col]
                    v = tryparse(Float64, strip(replace(s, r"([0-9])([+-])(\d)" => s"\1e\2\3")))
                    v !== nothing && v != 0.0 && push!(egn, v)
                    length(egn) >= ngn + 1 && break
                end
                idx += 1
            end
            break
        end
        idx += 1
    end

    # Parse MF3 sections for group-averaged XS
    idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = _parse_int(p[71:72]); mt = _parse_int(p[73:75])
        mat_val = _parse_int(p[67:70])
        if mat_val == mat && mf == 3 && mt > 0
            ng = _parse_int(p[56:66])
            idx += 1; vals = Float64[]
            while idx <= length(lines) && length(vals) < ng
                p2 = rpad(lines[idx], 80)
                mf2 = _parse_int(p2[71:72]); mt2 = _parse_int(p2[73:75])
                (mf2 != 3 || mt2 != mt) && break
                for col in 0:5
                    s = p2[1+11*col:11+11*col]
                    v = tryparse(Float64, strip(replace(s, r"([0-9])([+-])(\d)" => s"\1e\2\3")))
                    v !== nothing && push!(vals, v)
                    length(vals) >= ng && break
                end
                idx += 1
            end
            group_xs[mt] = vals
            mt in reaction_mts || push!(reaction_mts, mt)
        else
            idx += 1
        end
    end

    # Parse MFcov sections for covariance matrices
    idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = _parse_int(p[71:72]); mt = _parse_int(p[73:75])
        mat_val = _parse_int(p[67:70])
        if mat_val == mat && mf == mfcov && mt > 0
            ngn = length(egn) - 1
            # Read covariance matrix from the sub-sections
            idx += 1  # skip HEAD
            # Read sub-section headers and data
            while idx <= length(lines)
                p2 = rpad(lines[idx], 80)
                mf2 = _parse_int(p2[71:72]); mt2 = _parse_int(p2[73:75])
                (mf2 != mfcov || mt2 != mt) && break
                # Sub-section header: mt1 from L2
                mt1 = _parse_int(p2[23:33])
                ni = _parse_int(p2[56:66])
                idx += 1
                # Read NI-type blocks (each row of the matrix)
                matrix = zeros(Float64, ngn, ngn)
                row = 0
                for _ in 1:ni
                    idx > length(lines) && break
                    p3 = rpad(lines[idx], 80)
                    ne = _parse_int(p3[12:22]); lb = _parse_int(p3[23:33])
                    nvals = _parse_int(p3[45:55])
                    row_idx = _parse_int(p3[56:66])
                    idx += 1
                    # Read values
                    vals = Float64[]
                    while idx <= length(lines) && length(vals) < nvals
                        p4 = rpad(lines[idx], 80)
                        mf4 = _parse_int(p4[71:72]); mt4 = _parse_int(p4[73:75])
                        (mf4 != mfcov || mt4 != mt) && break
                        for col in 0:5
                            s = p4[1+11*col:11+11*col]
                            v = tryparse(Float64, strip(replace(s, r"([0-9])([+-])(\d)" => s"\1e\2\3")))
                            v !== nothing && push!(vals, v)
                            length(vals) >= nvals && break
                        end
                        idx += 1
                    end
                    if row_idx >= 1 && row_idx <= ngn && length(vals) >= ngn
                        matrix[row_idx, :] .= vals[1:ngn]
                    end
                end
                mt1_val = mt1 > 0 ? mt1 : mt
                cov_results[(mt, mt1_val)] = matrix
                mt in reaction_mts || push!(reaction_mts, mt)
            end
        else
            idx += 1
        end
    end

    sort!(reaction_mts)
    cov_results, reaction_mts, group_xs, egn
end

# =========================================================================
# Helpers
# =========================================================================

function _read_za_awr(endf_path::String, mat::Int)
    za = 0.0; awr = 0.0
    open(endf_path) do io
        while !eof(io)
            line = readline(io); length(line) < 75 && continue
            p = rpad(line, 80)
            mat_val = _parse_int(p[67:70]); mf = _parse_int(p[71:72])
            if mat_val == mat && mf > 0
                za = parse_endf_float(p[1:11])
                awr = parse_endf_float(p[12:22])
                return za, awr
            end
        end
    end
    za, awr
end

# =========================================================================
# Output tape writer
# =========================================================================

"""Write ERRORR output tape in GENDF-like format matching Fortran covout."""
function _write_errorr_tape(io::IO, mat::Int, za::Float64, awr::Float64,
                            egn::Vector{Float64}, group_xs::Dict{Int,Vector{Float64}},
                            cov_results::Dict{Tuple{Int,Int},Matrix{Float64}},
                            reaction_mts::Vector{Int}, mfcov::Int)
    ngn = length(egn) - 1
    nw = length(egn)

    # TPID record (blank)
    @printf(io, "%66s%4d%2d%3d%5d\n", "", 0, 0, 0, 0)

    # MF1/MT451 — directory with group boundaries
    seq = 1
    # HEAD: ZA, AWR, LRP=5, 0, NFC=-nw, 0
    _write_cont_line(io, za, awr, 5, 0, -nw, 0, mat, 1, 451, seq); seq += 1
    # CONT: 0, 0, NG=ngn, 0, NW=nw, 0
    _write_cont_line(io, 0.0, 0.0, ngn, 0, nw, 0, mat, 1, 451, seq); seq += 1
    # Write group boundaries (6 per line)
    idx = 1
    while idx <= nw
        buf = ""
        for col in 1:6
            idx > nw && break
            buf *= format_endf_float(egn[idx]); idx += 1
        end
        _write_data_line(io, buf, mat, 1, 451, seq); seq += 1
    end
    # SEND for MF1
    _write_send_line(io, mat, 1); seq = 1
    # FEND
    _write_fend_line(io, mat)

    # MF3 — group-averaged cross sections for each reaction
    for mt in reaction_mts
        haskey(group_xs, mt) || continue
        xs = group_xs[mt]
        seq = 1
        # HEAD: ZA, 0, 0, 0, NG, 0
        _write_cont_line(io, za, 0.0, 0, 0, ngn, 0, mat, 3, mt, seq); seq += 1
        # Write values (6 per line)
        idx_v = 1
        while idx_v <= length(xs)
            buf = ""
            for col in 1:6
                idx_v > length(xs) && break
                # Use free-format float (like Fortran's F10.7 for XS)
                buf *= _fmt_errorr_float(xs[idx_v]); idx_v += 1
            end
            _write_data_line(io, buf, mat, 3, mt, seq); seq += 1
        end
        _write_send_line(io, mat, 3)
    end
    _write_fend_line(io, mat)

    # MFcov — covariance matrices
    for mt in reaction_mts
        # Count sub-sections for this MT
        sub_keys = [(mt, mt2) for (mt1, mt2) in keys(cov_results) if mt1 == mt]
        isempty(sub_keys) && continue

        seq = 1
        nl = length(sub_keys)
        # HEAD: ZA, AWR, 0, 0, 0, NL
        _write_cont_line(io, za, awr, 0, 0, 0, nl, mat, mfcov, mt, seq); seq += 1

        for (mt1, mt2) in sub_keys
            matrix = cov_results[(mt1, mt2)]
            # Sub-section header: 0, 0, 0, MT2, 0, NG
            _write_cont_line(io, 0.0, 0.0, 0, mt2, 0, ngn, mat, mfcov, mt, seq); seq += 1

            # Write each row as an NI-type block (LB=1, one row per block)
            for row in 1:ngn
                # Block header: 0, 0, NE=ngn, LB=1, NP=ngn, row_idx
                _write_cont_line(io, 0.0, 0.0, ngn, 1, ngn, row, mat, mfcov, mt, seq); seq += 1
                # Write row values
                idx_v = 1
                while idx_v <= ngn
                    buf = ""
                    for col in 1:6
                        idx_v > ngn && break
                        buf *= _fmt_errorr_float(matrix[row, idx_v]); idx_v += 1
                    end
                    _write_data_line(io, buf, mat, mfcov, mt, seq); seq += 1
                end
            end
        end
        _write_send_line(io, mat, mfcov)
    end
    _write_fend_line(io, mat)

    # MEND + TEND
    _write_fend_line(io, 0)  # MEND
end

# =========================================================================
# ENDF line formatting helpers
# =========================================================================

function _write_cont_line(io::IO, c1, c2, l1, l2, n1, n2, mat, mf, mt, seq)
    s1 = format_endf_float(Float64(c1))
    s2 = format_endf_float(Float64(c2))
    @printf(io, "%s%s%11d%11d%11d%11d%4d%2d%3d%5d\n",
            s1, s2, l1, l2, n1, n2, mat, mf, mt, seq)
end

function _write_data_line(io::IO, data::String, mat, mf, mt, seq)
    @printf(io, "%-66s%4d%2d%3d%5d\n", data, mat, mf, mt, seq)
end

function _write_send_line(io::IO, mat, mf)
    @printf(io, "%66s%4d%2d%3d%5d\n", "", mat, mf, 0, 99999)
end

function _write_fend_line(io::IO, mat)
    @printf(io, "%66s%4d%2d%3d%5d\n", "", mat, 0, 0, 0)
end

"""Format a float for errorr output — matches Fortran's free-format style."""
function _fmt_errorr_float(x::Float64)
    if x == 0.0
        return "0.000000+0 "
    end
    # Use standard ENDF float format but with slight adjustment for errorr
    format_endf_float(x)
end
