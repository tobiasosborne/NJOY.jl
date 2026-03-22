# RECONR pipeline -- top-level reconstruction and legacy interface
#
# PROPOSAL B DESIGN: The pipeline is a chain of composable functions:
#   build_evaluator -> build_grid -> adaptive_reconstruct -> merge_background!
#
# Key differences from Fortran / direct translation:
#   1. `build_evaluator` returns a closure -- no global state, fully testable
#   2. `build_grid` is a pure function: MF2 + MF3 data in, grid out
#   3. The evaluator is a separate concern from the adaptive algorithm
#   4. `PointwiseMaterial` uses a Matrix for cross sections (columnar layout)
#   5. merge_background! works in-place on the matrix to avoid copies
#   6. Multiple dispatch on formalism types for the cross section evaluator
#
# Correspondence to NJOY2016 reconr.f90:
#   sigma (2571-2667) -> build_evaluator (returns a closure)
#   lunion (1771-2238) -> build_grid (union of MF2 nodes + MF3 breakpoints)
#   emerge (4646-4982) -> merge_background! (in-place matrix merge)
#   reconr (main)      -> reconstruct (top-level pipeline)

# ==========================================================================
# Step 4: Top-level pipeline
# ==========================================================================

"""
    reconstruct(endf_file; mat=0, err=0.001, errmax=nothing, errint=nothing,
                temperature=0.0) -> PointwiseMaterial

Full RECONR pipeline: read ENDF, reconstruct, merge, return pointwise data.

The pipeline is:
1. Read MF2 (resonance parameters) and MF3 (background cross sections)
2. Build cross section evaluator closure from MF2
3. Build initial energy grid from MF2 nodes + MF3 breakpoints
4. Run adaptive reconstruction (generic algorithm from adaptive_grid.jl)
5. Merge MF3 backgrounds into reconstructed data
6. Return `PointwiseMaterial`

# Arguments
- `endf_file`: path to ENDF-format file
- `mat`: MAT number (0 = auto-detect from file)
- `err`: fractional reconstruction tolerance
- `errmax`: relaxed tolerance (default 10*err)
- `errint`: integral tolerance (default err/20000)
- `temperature`: reconstruction temperature [K] (0.0 = zero temperature)

# Returns
`PointwiseMaterial` with columns [total, elastic, fission, capture].
"""
function reconstruct(endf_file::AbstractString;
                     mat::Integer = 0,
                     err::Real = 0.001,
                     errmax::Union{Nothing, Real} = nothing,
                     errint::Union{Nothing, Real} = nothing,
                     temperature::Real = 0.0)
    open(endf_file, "r") do io
        # Detect MAT number
        actual_mat = _detect_mat(io, mat)

        # Read MF2 (filter by MAT for multi-material tapes)
        seekstart(io)
        found = find_section(io, 2, 151; target_mat=actual_mat)
        found || error("reconstruct: MF2/MT151 not found for MAT=$actual_mat in $endf_file")
        mf2 = read_mf2(io)

        # Read MF3 sections
        mf3_sections = read_mf3_sections(io, actual_mat)

        # Build evaluator
        table = temperature > 0.0 ? build_faddeeva_table() : nothing
        xs_eval = build_evaluator(mf2; temperature=temperature, table=table)

        # Check if we need the RML fallback evaluator (LRF=7)
        rml_data = read_rml_data(io)
        if rml_data !== nothing
            rml_eval = build_rml_evaluator(rml_data)
            base_eval = xs_eval
            xs_eval = function(E::Float64)
                r1 = base_eval(E)
                r2 = rml_eval(E)
                return (r1[1] + r2[1], r1[2] + r2[2],
                        r1[3] + r2[3], r1[4] + r2[4])
            end
        end

        # Build grid
        initial_grid = build_grid(mf2, mf3_sections)

        # Filter grid to resonance range
        eresl, eresh, eresr = _resonance_bounds(mf2)

        # Fallback: scan raw MF2 CONT records for range boundaries when the
        # reader could not parse the resonance formalism (e.g. LRF=7).
        if isinf(eresl) || eresh == 0.0
            eresl, eresh, eresr = _resonance_bounds_from_file(io)
        end

        # Build config with step_guard_limit set to eresr (resolved upper bound)
        config = AdaptiveConfig(Float64(err);
                                errmax = Float64(something(errmax, 10 * err)),
                                errint = Float64(something(errint, err / 20000)),
                                step_guard_limit = eresr > 0.0 ? eresr : Inf)

        # Handle materials with no resonance parameters
        if isinf(eresl) || eresh == 0.0
            energies = sort(unique(initial_grid))
            if length(energies) < 2
                energies = [1.0e-5, 2.0e7]
            end
            # Linearize for 1/v: add midpoints, decade points, thermal point
            linearize_one_over_v!(energies, Float64(err))
            n_pts = length(energies)
            values = zeros(n_pts, 4)
            merge_background!(energies, values, mf3_sections, mf2)
            mt_list = [1, 2, 18, 102]
            return PointwiseMaterial(Int32(actual_mat), energies, values, mt_list)
        end

        res_grid = filter(e -> e >= eresl && e <= eresh, initial_grid)
        if isempty(res_grid) || res_grid[1] > eresl
            pushfirst!(res_grid, eresl)
        end
        if res_grid[end] < eresh
            push!(res_grid, eresh)
        end
        sort!(res_grid)
        unique!(res_grid)

        # Adaptive reconstruction (uses the generic algorithm)
        energies, values = adaptive_reconstruct(xs_eval, res_grid, config)

        # Extend grid to include MF3 energies outside the resonance range.
        # Fortran RECONR covers the full MF3 energy span, not just the
        # resonance region. Points outside [eresl, eresh] get zero resonance
        # cross sections but still receive MF3 backgrounds.
        mf3_extra = Float64[]
        for sec in mf3_sections
            mt = Int(sec.mt)
            (mt == 1 || mt == 3 || mt == 101) && continue
            for e in sec.tab.x
                if e > 0.0 && (e < eresl || e > eresh)
                    push!(mf3_extra, e)
                end
            end
        end
        if !isempty(mf3_extra)
            sort!(mf3_extra)
            unique!(mf3_extra)
            filter!(e -> !(e >= eresl && e <= eresh), mf3_extra)
            # Linearize the extension grid for 1/v representation
            linearize_one_over_v!(mf3_extra, Float64(err))
            filter!(e -> !(e >= eresl && e <= eresh), mf3_extra)
            n_extra = length(mf3_extra)
            extra_values = zeros(n_extra, 4)
            all_energies = vcat(mf3_extra, energies)
            all_values = vcat(extra_values, values)
            perm = sortperm(all_energies)
            energies = all_energies[perm]
            values = all_values[perm, :]
            mask = trues(length(energies))
            for i in 2:length(energies)
                if energies[i] == energies[i-1]
                    mask[i] = false
                end
            end
            if !all(mask)
                energies = energies[mask]
                values = values[mask, :]
            end
        end

        # Merge backgrounds
        merge_background!(energies, values, mf3_sections, mf2)

        mt_list = [1, 2, 18, 102]  # total, elastic, fission, capture
        return PointwiseMaterial(Int32(actual_mat), energies, values, mt_list)
    end
end

# ==========================================================================
# Convenience: legacy interface compatibility
# ==========================================================================

# Provide the old `reconr` name as an alias to `reconstruct` for compatibility
# with the existing test suite. Returns a NamedTuple instead of PointwiseMaterial.
function reconr(endf_file::AbstractString;
                mat::Integer = 0,
                err::Float64 = 0.001,
                errmax::Union{Nothing, Float64} = nothing,
                errint::Union{Nothing, Float64} = nothing)
    errmax_val = something(errmax, 10.0 * err)
    errint_val = something(errint, err / 20000.0)

    open(endf_file, "r") do io
        actual_mat = _detect_mat(io, mat)

        seekstart(io)
        found = find_section(io, 2, 151; target_mat=actual_mat)
        found || error("reconr: MF2/MT151 not found for MAT=$actual_mat in $endf_file")
        mf2 = read_mf2(io)

        mf3_sections = read_mf3_sections(io, actual_mat)

        # Check for RML (LRF=7) data and build fallback evaluator
        rml_data = read_rml_data(io)

        initial_grid = build_grid(mf2, mf3_sections)
        eresl, eresh, eresr = _resonance_bounds(mf2)

        # Fallback: if no parsed ranges were found, scan the raw MF2 CONT
        # records for resonance range boundaries.  This handles materials
        # whose formalism (e.g. LRF=7 R-Matrix Limited) is not yet fully
        # parsed by the reader but still have LRU=1 resonance data.
        if isinf(eresl) || eresh == 0.0
            eresl, eresh, eresr = _resonance_bounds_from_file(io)
        end

        # Handle materials with no resonances (LRU=0 only, e.g. H-2)
        if (isinf(eresl) || eresh == 0.0) && rml_data === nothing
            # For LRU=0, Fortran sets eresr/eresh = ehigh = 20e6 (constant,
            # reconr.f90:144,308-309). NOT from MF2 EH.
            ehigh_const = 20.0e6  # Fortran parameter ehigh=20.e6_kr

            # Build union grid matching Fortran lunion (reconr.f90:1771-2238).
            # elim = min(0.99e6, eresr) = min(0.99e6, 20e6) = 0.99e6
            all_energies = lunion_grid(mf3_sections, err;
                                       awr=mf2.AWR)

            n_pts = length(all_energies)
            res_xs = [CrossSections() for _ in 1:n_pts]
            merged_xs = merge_background_legacy(all_energies, res_xs, mf3_sections;
                                                awr=mf2.AWR)

            total_arr = [xs.total for xs in merged_xs]
            elastic_arr = [xs.elastic for xs in merged_xs]
            fission_arr = [xs.fission for xs in merged_xs]
            capture_arr = [xs.capture for xs in merged_xs]

            return (energies=all_energies, total=total_arr, elastic=elastic_arr,
                    fission=fission_arr, capture=capture_arr,
                    mf2=mf2, mf3_sections=mf3_sections)
        end

        # Build full union grid from MF3 sections (matching Fortran lunion,
        # which runs BEFORE resxs — reconr.f90:358)
        mf2_nodes = Float64[]
        _add_mf2_nodes!(mf2_nodes, mf2)
        bg_grid = lunion_grid(mf3_sections, err;
                              nodes=mf2_nodes, awr=mf2.AWR,
                              eresl=eresl, eresr=eresr, eresh=eresh)

        # Filter union grid to resonance range for adaptive reconstruction.
        # Fortran resxs skips below eresl, stops at eresh (reconr.f90:2335,2370).
        # The lunion grid already has shaded boundary points (e.g. sigfig(EL,7,+1)),
        # so we do NOT add exact eresl/eresh — those were removed by boundary filtering.
        res_grid = filter(e -> e >= eresl && e < eresh, bg_grid)
        isempty(res_grid) && error("reconr: no grid points in resonance range [$eresl, $eresh)")

        # Build cross section function: combine MF2 + RML evaluators
        xs_fn = if rml_data !== nothing
            rml_eval = build_rml_evaluator(rml_data)
            function(E)
                sig = sigma_mf2(E, mf2)
                rml_xs = rml_eval(Float64(E))
                CrossSections(sig.total + rml_xs[1], sig.elastic + rml_xs[2],
                              sig.fission + rml_xs[3], sig.capture + rml_xs[4])
            end
        else
            E -> sigma_mf2(E, mf2)
        end

        # For adaptive convergence, test only partials (elastic, fission, capture)
        # matching Fortran resxs which tests j=1..nsig-1 = partials, NOT total.
        # Total makes convergence stricter, producing extra grid points.
        xs_partials = E -> let xs = xs_fn(E)
            (xs.elastic, xs.fission, xs.capture)
        end

        # Adaptive reconstruction in resonance range (matching Fortran resxs)
        config = AdaptiveConfig(err; errmax=errmax_val, errint=errint_val,
                                step_guard_limit = eresr > 0.0 ? eresr : Inf)
        res_energies, _ = adaptive_reconstruct(xs_partials, res_grid, config)

        # Merge: lunion grid (outside resonance) + adaptive grid (inside resonance)
        # Fortran emerge merges bg grid (outside [eresl, eresh)) with resxs grid
        bg_outside = filter(e -> e < eresl || e >= eresh, bg_grid)
        all_energies = sort!(vcat(bg_outside, res_energies))
        _dedup_tol!(all_energies)

        # Evaluate resonance XS at all grid points
        n_pts = length(all_energies)
        res_xs = Vector{CrossSections{Float64}}(undef, n_pts)
        for i in 1:n_pts
            e = all_energies[i]
            res_xs[i] = (e >= eresl && e < eresh) ? xs_fn(e) : CrossSections()
        end

        merged_xs = merge_background_legacy(all_energies, res_xs, mf3_sections;
                                            awr=mf2.AWR)

        total_arr = [xs.total for xs in merged_xs]
        elastic_arr = [xs.elastic for xs in merged_xs]
        fission_arr = [xs.fission for xs in merged_xs]
        capture_arr = [xs.capture for xs in merged_xs]

        return (energies=all_energies, total=total_arr, elastic=elastic_arr,
                fission=fission_arr, capture=capture_arr,
                mf2=mf2, mf3_sections=mf3_sections)
    end
end

# ==========================================================================
# Raw MF2 resonance-range boundary scanner
# ==========================================================================

"""
    _resonance_bounds_from_file(io) -> (eresl, eresh, eresr)

Scan the MF2/MT151 section directly from the ENDF file to extract resonance
range boundaries (EL, EH, LRU) without requiring the formalism to be fully
parsed.  This handles materials like Sr-88 (LRF=7, R-Matrix Limited) whose
resonance parameters are skipped by `read_mf2` but still need their energy
bounds for the RECONR pipeline.

Returns the same `(eresl, eresh, eresr)` triple as `_resonance_bounds`.
"""
function _resonance_bounds_from_file(io::IO)
    eresl = Inf
    eresh = 0.0
    eresr = 0.0

    saved_pos = position(io)
    seekstart(io)

    # Find MF2/MT151
    if !find_section(io, 2, 151)
        seek(io, saved_pos)
        return (eresl, eresh, eresr)
    end

    # HEAD record: ZA, AWR, 0, 0, NIS, 0
    head = read_head(io)
    NIS = Int(head.N1)

    for _ in 1:NIS
        # Isotope CONT: ZAI, ABN, 0, LFW, NER, 0
        iso_cont = read_cont(io)
        NER = Int(iso_cont.N1)

        for _ in 1:NER
            # Range CONT: EL, EH, LRU, LRF, NRO, NAPS
            range_cont = read_cont(io)
            EL = range_cont.C1
            EH = range_cont.C2
            LRU = Int(range_cont.L1)

            if LRU > 0
                eresl = min(eresl, EL)
                eresh = max(eresh, EH)
                if LRU == 1
                    eresr = max(eresr, EH)
                end
            end

            # Skip the rest of this range's data (we only need the CONT line)
            while !eof(io)
                pos = position(io)
                line = readline(io)
                p = rpad(line, 80)
                mf = _parse_int(p[71:72])
                mt = _parse_int(p[73:75])
                if mf == 0 || mt == 0 || mf != 2 || mt != 151
                    seek(io, pos)
                    break
                end
            end
        end
    end

    eresr = clamp(eresr, eresl, eresh)
    seek(io, saved_pos)
    return (eresl, eresh, eresr)
end
