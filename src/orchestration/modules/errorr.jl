# errorr module runner -- Covariance processing from ENDF MF33/MF31
#
# Matches Fortran errorr.f90: read covariance data from ENDF, compute
# group-averaged cross sections, expand covariances onto union grid,
# write output tape in GENDF-like format.

"""
    errorr_dummy_mf33_stub!(tapes::TapeManager, nin_unit::Int, nout_unit::Int)

STUB for Fortran errorr's "999 option" — inserts dummy MF33 covariance
sections into an ENDF tape. Full implementation would synthesize
placeholder MF33 records for each requested MT (see Fortran errorr.f90
line ~900). For now, copy the input tape unchanged so downstream
reconr/broadr find a readable file.
"""
function errorr_dummy_mf33_stub!(tapes::TapeManager, nin_unit::Int, nout_unit::Int)
    if nin_unit <= 0 || nout_unit <= 0
        @warn "errorr 999-mode: bad units nin=$nin_unit nout=$nout_unit — skipping"
        return nothing
    end
    nin_path  = resolve(tapes, nin_unit)
    nout_path = resolve(tapes, nout_unit)
    @info "errorr: 999-mode STUB (cp tape$(nin_unit) → tape$(nout_unit); " *
          "dummy MF33 insertion not implemented)"
    if isfile(nin_path)
        nin_path != nout_path && cp(nin_path, nout_path; force=true)
        register!(tapes, nout_unit, nout_path)
    else
        @warn "errorr 999-mode: input tape$(nin_unit) at $nin_path not found"
    end
    nothing
end

"""
    errorr_module(tapes::TapeManager, params::ErrorrParams)

Run ERRORR: read ENDF covariance data, compute multigroup covariances,
write output tape. Matches Fortran errorr.f90 interface.

Input tapes: nendf (ENDF evaluation), npend (PENDF, optional), ngout (GENDF, optional)
Output tapes: nout (covariance output)
"""
function errorr_module(tapes::TapeManager, params::ErrorrParams)
    @info "errorr: MAT=$(params.mat) mfcov=$(params.mfcov) → tape $(params.nout)"

    # Graceful skip for malformed/unsupported deck shapes that left nout=0.
    # The dedicated 999-mode stub is dispatched separately in pipeline.jl.
    if params.nout <= 0
        @warn "errorr: nout=0 — skipping (unsupported deck shape)"
        return nothing
    end

    endf_path = resolve(tapes, params.nendf)
    nout_path = resolve(tapes, params.nout)

    # Read ZA and AWR from ENDF
    za, awr = _read_za_awr(endf_path, params.mat)

    # Find available MTs in the covariance file
    mfcov = params.mfcov > 0 ? params.mfcov : 33
    avail_mts = _find_mfcov_mts(endf_path, params.mat, mfcov)
    @info "errorr: MF$mfcov MTs available: $avail_mts"

    if params.nin > 0
        # Second errorr call: copy previous covariance tape + append new MT
        nin_path = resolve(tapes, params.nin)
        _errorr_second_call(tapes, params, endf_path, nin_path, nout_path, za, awr, mfcov)
        register!(tapes, params.nout, nout_path)
        @info "errorr: second call, tape $(params.nin) + new cov → tape $(params.nout)"
        return nothing
    end

    # Output group structure — matches Fortran egngpn (errorr.f90:9716).
    # ign=-1: union of user_egn + MFcov breakpoints (replaces egn).
    # ign=1:  user_egn exactly (2306 / 9 / 7 lines typical).
    # ign>=2: library structure (LANL-30 for ign=3, etc.).
    egn = if params.ign == -1
        _build_errorr_grid_from_endf(endf_path, params.mat, mfcov, params)
    else
        _errorr_output_grid(params)
    end
    ngn = length(egn) - 1
    @info "errorr: ign=$(params.ign) → $ngn groups"

    # Determine which reaction MTs are involved
    reaction_mts = sort(avail_mts)
    @info "errorr: reactions: $reaction_mts"

    # Compute group-averaged cross sections from PENDF
    group_xs = Dict{Int, Vector{Float64}}()
    if params.npend > 0 && params.iread == 0
        pendf_path = resolve(tapes, params.npend)
        group_xs = _errorr_group_average(pendf_path, params.mat, egn, params.iwt)
    end

    # Compute multigroup covariance matrices
    # Read MF33 sub-sections and expand NI blocks onto the group grid.
    # Only use NI blocks from standalone sub-sections (NC=0).
    # Sub-sections with NC>0 contain derived covariance that requires
    # NC processing (linear combinations of other MTs) — skip for now.
    cov_matrices = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    for mt in avail_mts
        try
            open(endf_path) do io
                seekstart(io)
                find_section(io, mfcov, mt; target_mat=params.mat) || return
                head = read_cont(io); nl = Int(head.N2)
                for _ in 1:nl
                    sh = read_cont(io)
                    mt2 = Int(sh.L2); mt2 == 0 && (mt2 = mt)
                    nc = Int(sh.N1); ni = Int(sh.N2)
                    # Read NC sub-subsections
                    for _ in 1:nc
                        try; read_list(io); catch; break; end
                    end
                    # Only expand NI blocks from standalone sub-sections (nc==0)
                    for _ in 1:ni
                        try
                            lst = read_list(io)
                            nc > 0 && continue  # Skip NI in derived sub-sections
                            lb = Int(lst.L2)
                            lb in (0,1,2,5,6) || continue
                            ne = Int(lst.N2); np = Int(lst.N1)
                            if lb in (0,1,2,3,4)
                                nk = div(np, 2)
                                nk == 0 && continue
                                block = CovarianceBlock(mt, mt2, lb, Int(lst.L1),
                                    lst.data[1:2:2nk], lst.data[2:2:2nk])
                            elseif lb == 5
                                ek = lst.data[1:ne]
                                fvals = lst.data[ne+1:np]
                                block = CovarianceBlock(mt, mt2, lb, Int(lst.L1), ek, fvals)
                            elseif lb == 6
                                nek = ne; nel = Int(lst.L1)
                                ek = vcat(lst.data[1:nek], lst.data[nek+1:nek+nel])
                                fvals = lst.data[nek+nel+1:np]
                                block = CovarianceBlock(mt, mt2, lb, Int(lst.L1), ek, fvals)
                            else
                                continue
                            end
                            C = expand_covariance_block(block, egn)
                            key = (mt, mt2)
                            if haskey(cov_matrices, key)
                                cov_matrices[key] .+= C
                            else
                                cov_matrices[key] = C
                            end
                        catch e
                            @debug "errorr: skipping NI block: $e"
                        end
                    end
                end
            end
        catch e
            @warn "errorr: covariance failed for MT=$mt: $e"
        end
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

"""Resolve the errorr output group structure from params.ign.

Mirrors Fortran gengpn (groupr.f90:1557): user-defined grid for
ign<0 or ign==1 (read from params.user_egn), library structure
otherwise (LANL-30 for ign=3, WIMS-69 for ign=9, RRD-50 for ign=5,
and whatever `get_group_structure` supports for the rest)."""
function _errorr_output_grid(params::ErrorrParams)
    ign = params.ign
    if ign <= 1
        isempty(params.user_egn) && error(
            "errorr: ign=$ign requires user_egn, got empty — check input deck parser")
        return collect(Float64, params.user_egn)
    end
    try
        return collect(Float64, get_group_structure(ign))
    catch
        @warn "errorr: unsupported ign=$ign, falling back to LANL-30"
        return collect(Float64, LANL_30)
    end
end

"""Build the errorr union energy grid directly from raw ENDF MFcov energy values.

DEPRECATED for output grid — use `_errorr_output_grid` instead. Kept
because the union of MF33 breakpoints is still needed internally for
covcal-style covariance expansion (see T03_phase7 T04 tape25 residual)."""
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

    for (mt, (energies, xs)) in mf3
        mt in (1, 451) && continue
        if iwt == 3
            # 1/E weighting: use analytical integration matching Fortran epanel
            # For piecewise-linear σ(E) on each panel [E1,E2]:
            #   ∫σ/E dE = (σ1·E2 - σ2·E1)/(E2-E1) · ln(E2/E1) + (σ2 - σ1)
            #   ∫1/E dE = ln(E2/E1)
            avg = _group_average_inv_e(energies, xs, egn)
        else
            # Flat weighting: standard panel integration
            n = length(energies)
            w_vals = [constant_weight(energies[i]) for i in 1:n]
            sig_w = [xs[i] * w_vals[i] for i in 1:n]
            flux = group_integrate(energies, w_vals, egn)
            sig_int = group_integrate(energies, sig_w, egn)
            avg = [flux[g] > 0 ? sig_int[g] / flux[g] : 0.0 for g in 1:ngn]
        end
        result[mt] = avg
    end
    result
end

"""Group averaging with 1/E weight matching Fortran errorr epanel.
Uses trapezoidal rule with 1% energy stepping (enext = 1.01*E),
replicating Fortran's egtwtf for iwt=3 which returns enext = s101*e
where s101 = 1.01. Between PENDF data points, σ is interpolated linearly
and 1/E is evaluated exactly at each sub-panel boundary."""
function _group_average_inv_e(energies::Vector{Float64}, xs::Vector{Float64},
                               egn::Vector{Float64})
    ngn = length(egn) - 1
    ne = length(energies)
    avg = zeros(Float64, ngn)
    ip = 1
    s101 = 1.01  # Fortran errorr stepping factor

    for g in 1:ngn
        elo, ehi = egn[g], egn[g+1]
        flux_g = 0.0; sigflux_g = 0.0

        # Advance to first panel overlapping this group
        while ip < ne - 1 && energies[ip+1] <= elo
            ip += 1
        end

        jp = ip
        while jp < ne && energies[jp] < ehi
            e1, e2 = energies[jp], energies[jp+1]
            s1, s2 = xs[jp], xs[jp+1]

            # Clip to group boundaries
            ea = max(e1, elo); eb = min(e2, ehi)
            ea >= eb && (jp += 1; continue)

            # Trapezoidal integration with 1% sub-panels matching Fortran epanel
            e_lo = ea
            while e_lo < eb
                e_hi = min(s101 * e_lo, eb)  # next weight function point or group boundary

                # Interpolate σ at sub-panel boundaries from PENDF data
                frac_lo = (e_lo - e1) / (e2 - e1)
                frac_hi = (e_hi - e1) / (e2 - e1)
                sig_lo = s1 + frac_lo * (s2 - s1)
                sig_hi = s1 + frac_hi * (s2 - s1)

                # Trapezoidal: ∫f dE ≈ (f(lo) + f(hi))/2 · ΔE
                bq = (e_hi - e_lo) / 2.0
                flux_lo = 1.0 / e_lo; flux_hi = 1.0 / e_hi
                flux_g += (flux_lo + flux_hi) * bq
                sigflux_g += (sig_lo * flux_lo + sig_hi * flux_hi) * bq

                e_lo = e_hi
            end

            jp += 1
        end

        avg[g] = flux_g > 0 ? sigflux_g / flux_g : 0.0
    end
    avg
end

# =========================================================================
# Passthrough for second errorr call
# =========================================================================

"""Second errorr call: copy previous covariance + append new MT covariance."""
function _errorr_second_call(tapes::TapeManager, params::ErrorrParams,
                              endf_path::String, nin_path::String, nout_path::String,
                              za::Float64, awr::Float64, mfcov::Int)
    mat = params.mat

    # Build second group structure from user boundaries
    egn2 = params.user_egn
    isempty(egn2) && error("errorr second call requires user energy grid")
    ngn2 = length(egn2) - 1

    # Read nubar from GENDF tape (ngout)
    nubar_values = Float64[]
    if params.ngout != 0
        gendf_path = resolve(tapes, abs(params.ngout))
        nubar_values = _read_gendf_nubar(gendf_path, mat, egn2)
    end

    # Read MF31 covariance from ENDF and collapse onto second group grid.
    # Fortran errorr's covcal works on a UNION grid (all block breakpoints +
    # group boundaries) then collapses — direct per-group expansion skips
    # sub-group structure from finer blocks and diverges from Fortran.
    cov_mt = 452
    cov_matrix = zeros(Float64, ngn2, ngn2)
    blocks = CovarianceBlock[]
    try
        open(endf_path) do io
            seekstart(io)
            find_section(io, mfcov, cov_mt; target_mat=mat) || return
            head = read_cont(io); nl = Int(head.N2)
            for _ in 1:nl
                sh = read_cont(io)
                nc = Int(sh.N1); ni = Int(sh.N2)
                for _ in 1:nc
                    try; read_list(io); catch; break; end
                end
                nc > 0 && continue
                for _ in 1:ni
                    try
                        lst = read_list(io)
                        lb = Int(lst.L2); lt = Int(lst.L1)
                        np = Int(lst.N1); ne = Int(lst.N2)
                        block = if lb in (0,1,2,3,4)
                            nk = div(np, 2)
                            CovarianceBlock(cov_mt, cov_mt, lb, lt,
                                lst.data[1:2:2nk], lst.data[2:2:2nk])
                        elseif lb == 5
                            CovarianceBlock(cov_mt, cov_mt, lb, lt,
                                lst.data[1:ne], lst.data[ne+1:np])
                        elseif lb == 6
                            CovarianceBlock(cov_mt, cov_mt, lb, lt,
                                lst.data[1:ne+lt], lst.data[ne+lt+1:np])
                        else
                            nothing
                        end
                        block !== nothing && push!(blocks, block)
                    catch; end
                end
                break
            end
        end
        if !isempty(blocks)
            cm = multigroup_covariance(blocks, egn2)
            cov_matrix = cm.matrix
        end
    catch e
        @warn "errorr: MF$mfcov/MT$cov_mt read failed: $e"
    end

    # Read nin content
    nin_lines = readlines(nin_path)
    # Find last content line (before MEND/TEND)
    last_content = length(nin_lines)
    while last_content > 0
        p = rpad(nin_lines[last_content], 80)
        mat_val = _parse_int(p[67:70])
        mat_val <= 0 || break
        last_content -= 1
    end

    open(nout_path, "w") do io
        # Phase A (continuous seq): TPID + copied nin content + MEND separator
        # + 2nd material's MF1/MT451. Matches referenceTape25 lines 0-84.
        seq = 0
        for i in 1:last_content
            p = rpad(nin_lines[i], 80)
            @printf(io, "%s%5d\n", p[1:75], seq)
            seq += 1
        end
        @printf(io, "%66s%4d%2d%3d%5d\n", "", 0, 0, 0, seq); seq += 1  # MEND separator

        _write_cont_line(io, za, awr, 5, 0, -11, 0, mat, 1, 451, seq); seq += 1
        _write_cont_line(io, 0.0, 0.0, ngn2, 0, length(egn2), 0, mat, 1, 451, seq); seq += 1
        idx = 1
        while idx <= length(egn2)
            buf = ""
            for col in 1:6
                idx > length(egn2) && break
                buf *= format_endf_float(egn2[idx]); idx += 1
            end
            _write_data_line(io, buf, mat, 1, 451, seq); seq += 1
        end

        # Phase B (standard ENDF seq): SEND=99999, FEND=0, section seq resets to 1.
        _write_send_line(io, mat, 1)
        _write_fend_line(io, mat)

        seq = 1
        _write_cont_line(io, za, 0.0, 0, 0, ngn2, 0, mat, 3, cov_mt, seq); seq += 1
        idx_v = 1
        while idx_v <= length(nubar_values)
            buf = ""
            for col in 1:6
                idx_v > length(nubar_values) && break
                buf *= _fmt_errorr_xs(nubar_values[idx_v]); idx_v += 1
            end
            _write_data_line(io, buf, mat, 3, cov_mt, seq); seq += 1
        end
        _write_send_line(io, mat, 3)
        _write_fend_line(io, mat)

        seq = 1
        _write_cont_line(io, za, awr, 0, 0, 0, 1, mat, 33, cov_mt, seq); seq += 1
        _write_cont_line(io, 0.0, 0.0, 0, cov_mt, 0, ngn2, mat, 33, cov_mt, seq); seq += 1
        for row in 1:ngn2
            _write_cont_line(io, 0.0, 0.0, ngn2, 1, ngn2, row, mat, 33, cov_mt, seq); seq += 1
            idx_v = 1
            while idx_v <= ngn2
                buf = ""
                for col in 1:6
                    idx_v > ngn2 && break
                    buf *= _fmt_errorr_cov(cov_matrix[row, idx_v]); idx_v += 1
                end
                _write_data_line(io, buf, mat, 33, cov_mt, seq); seq += 1
            end
        end
        _write_send_line(io, mat, 33)
        _write_fend_line(io, mat)

        _write_fend_line(io, 0)
        @printf(io, "%66s%4d%2d%3d%5d\n", "", -1, 0, 0, 0)
    end
end

"""Read group-averaged nubar from GENDF tape for specified energy groups."""
function _read_gendf_nubar(gendf_path::String, mat::Int, egn::Vector{Float64})
    # Read nubar tabulated data from the GENDF MF3/MT452 records
    ngn = length(egn) - 1
    values = Float64[]

    lines = readlines(gendf_path)
    idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = _parse_int(p[71:72]); mt = _parse_int(p[73:75])
        mat_val = _parse_int(p[67:70])
        if mat_val == mat && mf == 3 && mt == 452
            # Found MT452 section — read group records
            idx += 1
            while idx <= length(lines)
                p2 = rpad(lines[idx], 80)
                mf2 = _parse_int(p2[71:72]); mt2 = _parse_int(p2[73:75])
                (mf2 != 3 || mt2 != 452) && break
                # Group record header: check if it's data (not CONT header)
                nw = _parse_int(p2[45:55])
                ig = _parse_int(p2[56:66])
                if ig > 0 && nw > 0
                    idx += 1
                    # Next line has the data values (flux, nubar, sigf_avg)
                    length(lines[idx]) < 22 && (idx += 1; continue)
                    p3 = rpad(lines[idx], 80)
                    # Second value is nubar
                    nubar = parse_endf_float(p3[12:22])
                    push!(values, nubar)
                end
                idx += 1
            end
            break
        end
        idx += 1
    end

    # Now interpolate the GENDF nubar (on LANL-30 grid) onto the user group grid
    # by finding which GENDF group each user group falls into
    if length(values) > 0
        # Read the GENDF group boundaries from MF1
        gendf_egn = _read_gendf_group_bounds(gendf_path, mat)
        if length(gendf_egn) > 1
            return _regroup_nubar(values, gendf_egn, egn)
        end
    end
    values
end

"""Read group boundaries from a GENDF tape's MF1/MT451."""
function _read_gendf_group_bounds(gendf_path::String, mat::Int)
    egn = Float64[]
    lines = readlines(gendf_path)
    idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = _parse_int(p[71:72]); mt = _parse_int(p[73:75])
        if mf == 1 && mt == 451
            idx += 1  # Skip first HEAD
            p2 = rpad(lines[idx], 80)
            ngn = _parse_int(p2[23:33])
            nw = _parse_int(p2[45:55])
            idx += 1
            # Read boundaries (skip first 2 values: elow, ehigh)
            all_vals = Float64[]
            while idx <= length(lines) && length(all_vals) < nw
                p3 = rpad(lines[idx], 80)
                _parse_int(p3[73:75]) != 451 && break
                for col in 0:5
                    s = p3[1+11*col:11+11*col]
                    v = tryparse(Float64, strip(replace(s, r"([0-9])([+-])(\d)" => s"\1e\2\3")))
                    v !== nothing && push!(all_vals, v)
                end
                idx += 1
            end
            # Skip elow (0.0) and ehigh (1e10), then take ngn+1 boundaries
            if length(all_vals) >= 3
                egn = all_vals[3:min(end, 2+ngn+1)]
            end
            break
        end
        idx += 1
    end
    egn
end

"""Regroup nubar from GENDF groups to user groups by flux weighting."""
function _regroup_nubar(gendf_nubar::Vector{Float64}, gendf_egn::Vector{Float64},
                         user_egn::Vector{Float64})
    ngn_user = length(user_egn) - 1
    result = Vector{Float64}(undef, ngn_user)

    for g in 1:ngn_user
        elo, ehi = user_egn[g], user_egn[g+1]
        # Find weighted average of GENDF nubar values overlapping this user group
        wt_sum = 0.0; nu_wt_sum = 0.0
        for k in 1:length(gendf_nubar)
            k >= length(gendf_egn) && break
            glo, ghi = gendf_egn[k], gendf_egn[k+1]
            # Overlap
            olo = max(elo, glo); ohi = min(ehi, ghi)
            olo >= ohi && continue
            # Weight by lethargy width of overlap
            wt = log(ohi / olo)
            wt_sum += wt
            nu_wt_sum += wt * gendf_nubar[k]
        end
        result[g] = wt_sum > 0 ? nu_wt_sum / wt_sum : 0.0
    end
    result
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
    # N1=-11 is the errorr convention (errorr.f90:5927), not a computed count.
    seq = 1
    _write_cont_line(io, za, awr, 5, 0, -11, 0, mat, 1, 451, seq); seq += 1
    _write_cont_line(io, 0.0, 0.0, ngn, 0, nw, 0, mat, 1, 451, seq); seq += 1
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
    # Sequence resets to 1 at each SEND (new MT within same MF).
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
        sub_keys = [(mt, mt2) for (mt1, mt2) in keys(cov_matrices) if mt1 == mt]
        for mt2 in reaction_mts
            mt2 <= mt && continue
            (mt, mt2) in sub_keys || push!(sub_keys, (mt, mt2))
        end
        sort!(sub_keys, by = x -> x[2])

        isempty(sub_keys) && continue
        nl = length(sub_keys)
        seq = 1
        _write_cont_line(io, za, awr, 0, 0, 0, nl, mat, mfcov, mt, seq); seq += 1

        for (mt1, mt2) in sub_keys
            matrix = get(cov_matrices, (mt1, mt2), nothing)
            _write_cont_line(io, 0.0, 0.0, 0, mt2, 0, ngn, mat, mfcov, mt, seq); seq += 1

            if matrix !== nothing && any(!iszero, matrix)
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
            else
                _write_cont_line(io, 0.0, 0.0, 1, 9, 1, ngn, mat, mfcov, mt, seq); seq += 1
                _write_data_line(io, format_endf_float(0.0), mat, mfcov, mt, seq); seq += 1
            end
        end
        _write_send_line(io, mat, mfcov)
    end
    _write_fend_line(io, mat)

    # MEND + TEND
    _write_fend_line(io, 0)
    @printf(io, "%66s%4d%2d%3d%5d\n", "", -1, 0, 0, 0)
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

"""Format a float for errorr XS output — uses ENDF a11 format (same as Fortran)."""
function _fmt_errorr_xs(x::Float64)
    format_endf_float(x)
end

"""Format a float for errorr covariance output — scientific notation."""
function _fmt_errorr_cov(x::Float64)
    if x == 0.0
        return " 0.0       "
    end
    format_endf_float(x)
end
