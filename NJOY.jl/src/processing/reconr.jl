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

        # Read MF2
        seekstart(io)
        found = find_section(io, 2, 151)
        found || error("reconstruct: MF2/MT151 not found in $endf_file")
        mf2 = read_mf2(io)

        # Read MF3 sections
        mf3_sections = read_mf3_sections(io, actual_mat)

        # Build evaluator
        table = temperature > 0.0 ? build_faddeeva_table() : nothing
        xs_eval = build_evaluator(mf2; temperature=temperature, table=table)

        # Build grid
        initial_grid = build_grid(mf2, mf3_sections)

        # Filter grid to resonance range
        eresl, eresh, eresr = _resonance_bounds(mf2)

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
        found = find_section(io, 2, 151)
        found || error("reconr: MF2/MT151 not found in $endf_file")
        mf2 = read_mf2(io)

        mf3_sections = read_mf3_sections(io, actual_mat)

        initial_grid = build_grid(mf2, mf3_sections)
        eresl, eresh, eresr = _resonance_bounds(mf2)

        # Handle materials with no resonances (LRU=0 only, e.g. H-2)
        if eresl == Inf || eresh == 0.0
            # No resonance range -- use MF3 energies directly with zero resonance XS
            # Build grid from MF3 breakpoints only
            all_energies = Float64[]
            for sec in mf3_sections
                append!(all_energies, sec.tab.x)
            end
            if isempty(all_energies)
                push!(all_energies, 1.0e-5)
                push!(all_energies, 2.0e7)
            end
            sort!(all_energies)
            unique!(all_energies)
            filter!(e -> e > 0.0, all_energies)

            n_pts = length(all_energies)
            res_xs = [CrossSections() for _ in 1:n_pts]
            merged_xs = merge_background_legacy(all_energies, res_xs, mf3_sections)

            total_arr = [xs.total for xs in merged_xs]
            elastic_arr = [xs.elastic for xs in merged_xs]
            fission_arr = [xs.fission for xs in merged_xs]
            capture_arr = [xs.capture for xs in merged_xs]

            return (energies=all_energies, total=total_arr, elastic=elastic_arr,
                    fission=fission_arr, capture=capture_arr,
                    mf2=mf2, mf3_sections=mf3_sections)
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

        xs_fn = E -> sigma_mf2(E, mf2)

        config = AdaptiveConfig(err; errmax=errmax_val, errint=errint_val,
                                step_guard_limit = eresr > 0.0 ? eresr : Inf)
        energies, values = adaptive_reconstruct(xs_fn, res_grid, config)

        n_pts = length(energies)
        res_xs = Vector{CrossSections}(undef, n_pts)
        for i in 1:n_pts
            res_xs[i] = CrossSections(values[i, 1], values[i, 2],
                                       values[i, 3], values[i, 4])
        end

        merged_xs = merge_background_legacy(energies, res_xs, mf3_sections)

        total_arr = [xs.total for xs in merged_xs]
        elastic_arr = [xs.elastic for xs in merged_xs]
        fission_arr = [xs.fission for xs in merged_xs]
        capture_arr = [xs.capture for xs in merged_xs]

        return (energies=energies, total=total_arr, elastic=elastic_arr,
                fission=fission_arr, capture=capture_arr,
                mf2=mf2, mf3_sections=mf3_sections)
    end
end
