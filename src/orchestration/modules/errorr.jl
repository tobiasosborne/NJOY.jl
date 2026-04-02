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

    # Read ZA and AWR from ENDF
    za, awr = _read_za_awr(endf_path, params.mat)

    # Find available MTs in the covariance file
    mfcov = params.mfcov > 0 ? params.mfcov : 33
    avail_mts = _find_mfcov_mts(endf_path, params.mat, mfcov)
    @info "errorr: MF$mfcov MTs available: $avail_mts"

    if params.nin > 0
        # Second errorr call: read covariance from previous output tape
        nin_path = resolve(tapes, params.nin)
        _errorr_passthrough(nin_path, nout_path, params)
        register!(tapes, params.nout, nout_path)
        @info "errorr: passthrough from tape $(params.nin) → tape $(params.nout)"
        return nothing
    end

    # Build union energy grid from covariance energies + user boundaries
    egn = _build_errorr_grid_from_endf(endf_path, params.mat, mfcov, params)
    ngn = length(egn) - 1
    @info "errorr: $ngn groups"

    # Determine which reaction MTs are involved
    reaction_mts = sort(avail_mts)
    @info "errorr: reactions: $reaction_mts"

    # Compute group-averaged cross sections from PENDF
    group_xs = Dict{Int, Vector{Float64}}()
    if params.npend > 0 && params.iread == 0
        pendf_path = resolve(tapes, params.npend)
        group_xs = _errorr_group_average(pendf_path, params.mat, egn, params.iwt)
    end

    # Compute multigroup covariance matrices using process_covariance
    cov_results = process_covariance(endf_path, egn; mts=avail_mts)
    cov_matrices = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    for cm in cov_results
        cov_matrices[(cm.mt1, cm.mt2)] = cm.matrix
    end

    # Write output tape
    open(nout_path, "w") do io
        _write_errorr_tape(io, params.mat, za, awr, egn, group_xs,
                          cov_matrices, reaction_mts, mfcov)
    end
    register!(tapes, params.nout, nout_path)

    lines = countlines(nout_path)
    @info "errorr: wrote $nout_path ($lines lines, $ngn groups, $(length(reaction_mts)) reactions)"
    nothing
end

# =========================================================================
# MF33/MF31 discovery
# =========================================================================

function _find_mfcov_mts(endf_path::String, mat::Int, mfcov::Int)
    mts = Int[]
    open(endf_path) do io
        while !eof(io)
            line = readline(io); length(line) < 75 && continue
            p = rpad(line, 80)
            mat_val = _parse_int(p[67:70])
            mf_val = _parse_int(p[71:72])
            mt_val = Int(_parse_int(p[73:75]))
            if mat_val == mat && mf_val == mfcov && mt_val > 0 && !(mt_val in mts)
                push!(mts, mt_val)
            end
        end
    end
    sort!(mts)
end

# =========================================================================
# Group structure
# =========================================================================

"""Build the errorr union energy grid directly from raw ENDF MFcov energy values."""
function _build_errorr_grid_from_endf(endf_path::String, mat::Int, mfcov::Int,
                                       params::ErrorrParams)
    pts = Set{Float64}()
    for e in params.user_egn; push!(pts, e); end

    # Scan ENDF MFcov sections for energy values in LIST records
    open(endf_path) do io
        while !eof(io)
            line = readline(io); length(line) < 75 && continue
            p = rpad(line, 80)
            mat_val = _parse_int(p[67:70])
            mf_val = _parse_int(p[71:72])
            mt_val = _parse_int(p[73:75])

            if mat_val == mat && mf_val == mfcov && mt_val > 0
                # We're in an MFcov section — scan for energy-like values
                # The LIST records in NI sub-subsections start with energies
                # Just extract all floats from this section that look like energies
                # (positive, in the nuclear physics range 1e-5 to 2e7)
                val = parse_endf_float(p[1:11])
                if val > 0 && val <= 2.1e7 && val != Float64(mat_val)
                    # Additional heuristic: energies in LB=5 blocks are in the first NE items
                    push!(pts, val)
                end
            end
        end
    end

    # Better approach: read MF33 and extract energies from the first NI sub-section of each MT
    for mt in _find_mfcov_mts(endf_path, mat, mfcov)
        try
            open(endf_path) do io
                seekstart(io)
                find_section(io, mfcov, mt; target_mat=mat) || return
                head = read_cont(io); nl = Int(head.N2)
                for il in 1:nl
                    sh = read_cont(io)
                    nc = Int(sh.N1); ni = Int(sh.N2)
                    # Skip NC sub-subsections
                    for _ in 1:nc
                        try; read_list(io); catch; break; end
                    end
                    # Read NI sub-subsections (energies are here)
                    for _ in 1:ni
                        try
                            lst = read_list(io)
                            lb = Int(lst.L2)
                            if lb in (0,1,2,3,4)
                                nk = div(Int(lst.N1), 2)
                                for k in 1:nk
                                    push!(pts, lst.data[2k-1])
                                end
                            elseif lb in (5,6)
                                ne = Int(lst.N2)
                                for k in 1:ne
                                    push!(pts, lst.data[k])
                                end
                            end
                        catch
                            break
                        end
                    end
                    break  # Only read first sub-section for energies
                end
            end
        catch e
            @warn "errorr: grid extraction failed for MT=$mt: $e"
        end
    end

    grid = sort!(collect(pts))

    # Clamp to user boundaries if provided
    if !isempty(params.user_egn) && length(params.user_egn) >= 2
        elow, ehigh = extrema(params.user_egn)
        filter!(e -> e >= elow && e <= ehigh, grid)
        elow in grid || pushfirst!(grid, elow)
        ehigh in grid || push!(grid, ehigh)
        sort!(unique!(grid))
    end

    grid
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
        mt in (1, 451) && continue

        # Build product arrays at data points
        n = length(energies)
        w_vals = [wfn(energies[i]) for i in 1:n]
        sig_w = [xs[i] * w_vals[i] for i in 1:n]

        # Exact panel integration
        flux = group_integrate(energies, w_vals, egn)
        sig_int = group_integrate(energies, sig_w, egn)

        avg = [flux[g] > 0 ? sig_int[g] / flux[g] : 0.0 for g in 1:ngn]
        result[mt] = avg
    end
    result
end

# =========================================================================
# Passthrough for second errorr call
# =========================================================================

"""Copy covariance data from one errorr output to another (for chained errorr calls)."""
function _errorr_passthrough(nin_path::String, nout_path::String, params::ErrorrParams)
    # For the second errorr call, we pass through the first errorr's output
    # with potential modifications (different mfcov, different group structure)
    # For now, just copy the input tape
    cp(nin_path, nout_path; force=true)
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
                            cov_matrices::Dict{Tuple{Int,Int},Matrix{Float64}},
                            reaction_mts::Vector{Int}, mfcov::Int)
    ngn = length(egn) - 1
    nw = length(egn)

    # TPID record (blank)
    @printf(io, "%66s%4d%2d%3d%5d\n", "", 0, 0, 0, 0)

    # MF1/MT451 — directory with group boundaries
    seq = 1
    _write_cont_line(io, za, awr, 5, 0, -nw, 0, mat, 1, 451, seq); seq += 1
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
    _write_send_line(io, mat, 1)
    _write_fend_line(io, mat)

    # MF3 — group-averaged cross sections for each reaction
    for mt in reaction_mts
        haskey(group_xs, mt) || continue
        xs = group_xs[mt]
        length(xs) == ngn || continue
        seq = 1
        _write_cont_line(io, za, 0.0, 0, 0, ngn, 0, mat, 3, mt, seq); seq += 1
        idx_v = 1
        while idx_v <= length(xs)
            buf = ""
            for col in 1:6
                idx_v > length(xs) && break
                buf *= _fmt_errorr_xs(xs[idx_v]); idx_v += 1
            end
            _write_data_line(io, buf, mat, 3, mt, seq); seq += 1
        end
        _write_send_line(io, mat, 3)
    end
    _write_fend_line(io, mat)

    # MFcov — covariance matrices
    for mt in reaction_mts
        # Find all sub-sections for this MT
        sub_keys = [(mt1, mt2) for (mt1, mt2) in keys(cov_matrices) if mt1 == mt]
        isempty(sub_keys) && continue

        seq = 1
        nl = length(sub_keys)
        _write_cont_line(io, za, awr, 0, 0, 0, nl, mat, mfcov, mt, seq); seq += 1

        for (mt1, mt2) in sub_keys
            matrix = cov_matrices[(mt1, mt2)]
            # Sub-section header: 0, 0, 0, MT2, 0, NG
            _write_cont_line(io, 0.0, 0.0, 0, mt2, 0, ngn, mat, mfcov, mt, seq); seq += 1

            # Each row as an NI-type block (LB=1 symmetric, one row per block)
            for row in 1:ngn
                _write_cont_line(io, 0.0, 0.0, ngn, 1, ngn, row, mat, mfcov, mt, seq); seq += 1
                idx_v = 1
                while idx_v <= ngn
                    buf = ""
                    for col in 1:6
                        idx_v > ngn && break
                        buf *= _fmt_errorr_cov(matrix[row, idx_v]); idx_v += 1
                    end
                    _write_data_line(io, buf, mat, mfcov, mt, seq); seq += 1
                end
            end
        end
        _write_send_line(io, mat, mfcov)
    end
    _write_fend_line(io, mat)

    # MEND
    _write_fend_line(io, 0)
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

"""Format a float for errorr XS output — Fortran free-format style (F10.7)."""
function _fmt_errorr_xs(x::Float64)
    # The errorr XS output uses free-format floats, not ENDF a11 format
    s = @sprintf("%.7g", x)
    rpad(s, 11)
end

"""Format a float for errorr covariance output — scientific notation."""
function _fmt_errorr_cov(x::Float64)
    if x == 0.0
        return " 0.0       "
    end
    format_endf_float(x)
end
