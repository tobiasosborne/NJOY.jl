# RECONR pipeline -- composable pointwise cross section reconstruction
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
# Data types
# ==========================================================================

"""
    MF3Section

A single MF3/MF13/MF23 cross section section from ENDF.
Stores the reaction type and the tabulated cross section data.
"""
struct MF3Section
    mt::Int32
    QM::Float64
    QI::Float64
    tab::TabulatedFunction
end

"""
    ENDFMaterial

All ENDF data needed for RECONR processing of one material.
"""
struct ENDFMaterial
    mat::Int32                              # MAT number
    awr::Float64                            # atomic weight ratio
    mf2::MF2Data                            # File 2 resonance parameters
    mf3_sections::Vector{MF3Section}        # File 3 cross section sections
    redundant_mts::Set{Int}                 # MTs that are sums (1, 3, 4, 18, 101, 27)
end

"""
    PointwiseMaterial

Result of RECONR processing: pointwise cross sections on a linearized grid.
"""
struct PointwiseMaterial
    mat::Int32
    energies::Vector{Float64}
    cross_sections::Matrix{Float64}  # (n_energies, n_reactions)
    mt_list::Vector{Int}             # which MT numbers are stored in each column
end

# Column indices for the primary 4 resonance channels
const _COL_TOTAL   = 1
const _COL_ELASTIC = 2
const _COL_FISSION = 3
const _COL_CAPTURE = 4

# ==========================================================================
# Step 1: Build evaluator (closure over MF2 data)
# ==========================================================================

"""
    build_evaluator(mf2::MF2Data; temperature=0.0, table=nothing)

Build a cross section evaluator closure from MF2 resonance data.
Returns a callable `f(E::Float64) -> NTuple{4, Float64}` that computes
`(total, elastic, fission, capture)` at energy `E`.

The returned function:
- Sums over all isotope sections weighted by abundance
- Dispatches to the correct formalism via `cross_section`
- Clamps negative components to zero
- Is fully type-stable (returns `NTuple{4, Float64}`)

This replaces NJOY2016's `sigma` subroutine (reconr.f90:2571-2667).
"""
function build_evaluator(mf2::MF2Data;
                         temperature::Real = 0.0,
                         table::Union{Nothing, FaddeevaTable} = nothing)
    # Pre-collect the isotope/range/abundance triples to avoid
    # iterating nested structures in the hot path
    sections = _collect_sections(mf2)

    function evaluator(E::Float64)
        sig_t = 0.0
        sig_e = 0.0
        sig_f = 0.0
        sig_c = 0.0

        @inbounds for (rng, abn) in sections
            if E >= rng.EL && E < rng.EH && Int(rng.LRU) > 0
                sigp = cross_section(E, rng; temperature=temperature, table=table)
                sig_t += max(0.0, sigp.total) * abn
                sig_e += max(0.0, sigp.elastic) * abn
                sig_f += max(0.0, sigp.fission) * abn
                sig_c += max(0.0, sigp.capture) * abn
            end
        end
        return (sig_t, sig_e, sig_f, sig_c)
    end

    return evaluator
end

# Flatten isotope/range hierarchy into a flat vector of (range, abundance) pairs
function _collect_sections(mf2::MF2Data)
    result = Tuple{ResonanceRange, Float64}[]
    for iso in mf2.isotopes
        abn = iso.ABN
        for rng in iso.ranges
            push!(result, (rng, abn))
        end
    end
    return result
end

# ==========================================================================
# Step 2: Build initial energy grid
# ==========================================================================

"""
    build_grid(mf2::MF2Data, mf3_sections::Vector{MF3Section}) -> Vector{Float64}

Build the initial energy grid for adaptive reconstruction by taking the
union of:
1. MF2 resonance energy nodes (range boundaries, resonance peaks, half-widths)
2. MF3 energy breakpoints from all non-redundant cross section tables
3. Standard anchor points (1e-5 eV, 0.0253 eV thermal)

This replaces the node-building logic from NJOY2016's `rdfil2` and `lunion`.
"""
function build_grid(mf2::MF2Data, mf3_sections::Vector{MF3Section})
    nodes = Float64[]
    sizehint!(nodes, 2048)

    # Add MF2 resonance nodes
    _add_mf2_nodes!(nodes, mf2)

    # Add MF3 breakpoints (filtered: skip redundant MT=1, MT=3)
    for sec in mf3_sections
        mt = Int(sec.mt)
        (mt == 1 || mt == 3 || mt == 101) && continue
        for e in sec.tab.x
            push!(nodes, e)
        end
    end

    # Add standard anchor points
    push!(nodes, 1.0e-5)
    push!(nodes, 0.0253)   # thermal point

    # Sort and deduplicate
    sort!(nodes)
    unique!(nodes)

    # Remove any non-positive energies
    filter!(e -> e > 0.0, nodes)

    return nodes
end

"""
    _add_mf2_nodes!(nodes, mf2)

Add resonance energy nodes from MF2 data:
- Range boundaries [EL, EH] with sigfig shading
- Resonance peak energies Er
- Half-width offsets Er +/- Gamma_total/2
"""
function _add_mf2_nodes!(nodes::Vector{Float64}, mf2::MF2Data)
    elow = 1.0e-5
    small = 1.0e-6

    for iso in mf2.isotopes
        for rng in iso.ranges
            Int(rng.LRU) == 0 && continue
            el, eh = rng.EL, rng.EH

            # Range boundary nodes with shading
            if abs(el - elow) > small
                push!(nodes, round_sigfig(el, 7, -1))
                push!(nodes, round_sigfig(el, 7, +1))
            else
                push!(nodes, el)
            end
            push!(nodes, round_sigfig(eh, 7, -1))
            push!(nodes, round_sigfig(eh, 7, +1))

            # Thermal point if range spans it
            if el < 0.0253
                push!(nodes, 0.0253)
            end

            # Resonance peak nodes
            _add_peak_nodes!(nodes, rng.parameters, el, eh)
        end
    end
end

# Dispatch on formalism type for resonance peak nodes
function _add_peak_nodes!(nodes, params::SLBWParameters, el, eh)
    _add_bw_peaks!(nodes, params.Er, params.Gn, params.Gg,
                   params.Gf, params.NLS, el, eh)
end

function _add_peak_nodes!(nodes, params::MLBWParameters, el, eh)
    _add_bw_peaks!(nodes, params.Er, params.Gn, params.Gg,
                   params.Gf, params.NLS, el, eh)
end

function _add_peak_nodes!(nodes, params::ReichMooreParameters, el, eh)
    for il in 1:Int(params.NLS)
        for ir in eachindex(params.Er[il])
            er = params.Er[il][ir]
            (er <= el || er > eh) && continue
            gn = params.Gn[il][ir]
            gg = params.Gg[il][ir]
            gfa = params.Gfa[il][ir]
            gfb = params.Gfb[il][ir]
            hw = (abs(gn) + abs(gg) + abs(gfa) + abs(gfb)) / 2.0
            _push_peak_triplet!(nodes, er, hw, el, eh)
        end
    end
end

# Fallback for unsupported formalisms
function _add_peak_nodes!(nodes, ::AbstractResonanceFormalism, el, eh)
    # No peak nodes for unsupported formalisms
end

# Shared logic for Breit-Wigner peak nodes
function _add_bw_peaks!(nodes, Er, Gn, Gg, Gf, NLS, el, eh)
    for il in 1:Int(NLS)
        for ir in eachindex(Er[il])
            er = Er[il][ir]
            (er <= el || er > eh) && continue
            gn = Gn[il][ir]
            gg = Gg[il][ir]
            gf = Gf[il][ir]
            hw = (abs(gn) + abs(gg) + abs(gf)) / 2.0
            _push_peak_triplet!(nodes, er, hw, el, eh)
        end
    end
end

# Push (Er - hw, Er, Er + hw) with appropriate rounding
function _push_peak_triplet!(nodes, er, hw, el, eh)
    ndig = 5
    if er > 0.0 && hw > 0.0
        ndig = clamp(2 + round(Int, log10(er / max(hw / 10, 1e-30))), 5, 9)
    end
    push!(nodes, round_sigfig(er, ndig, 0))
    er_lo = er - hw
    er_hi = er + hw
    if er_lo > el
        push!(nodes, round_sigfig(er_lo, ndig, 0))
    end
    if er_hi < eh
        push!(nodes, round_sigfig(er_hi, ndig, 0))
    end
end

# ==========================================================================
# Step 3: Merge MF3 background into pointwise result
# ==========================================================================

"""
    merge_background!(energies, values, mf3_sections, mf2)

Merge MF3 background cross sections into the resonance cross section
matrix `values` in-place. This replaces NJOY2016's `emerge`.

For each energy and matching MF3 reaction:
- MT 2 -> add to elastic column
- MT 18, 19 -> add to fission column
- MT 102 -> add to capture column
Then recompute total = elastic + fission + capture.
"""
function merge_background!(energies::Vector{Float64},
                            values::Matrix{Float64},
                            mf3_sections::Vector{MF3Section},
                            mf2::MF2Data)
    n = length(energies)
    small = 1.0e-8

    # Find resonance range boundaries for overlap detection
    eresl, eresh, eresr = _resonance_bounds(mf2)

    for i in 1:n
        e = energies[i]

        for sec in mf3_sections
            mt = Int(sec.mt)
            # Skip redundant reactions
            (mt == 1 || mt == 3 || mt == 101) && continue

            bg = interpolate(sec.tab, e)
            bg == 0.0 && continue

            # In the unresolved-resolved overlap region, backgrounds for
            # primary channels are assigned to the unresolved component
            # (matching NJOY's emerge logic)
            if e >= eresr && e < eresh && (mt == 2 || mt == 18 || mt == 19 || mt == 102)
                continue
            end

            if mt == 2
                values[i, _COL_ELASTIC] += bg
            elseif mt == 18 || mt == 19
                values[i, _COL_FISSION] += bg
            elseif mt == 102
                values[i, _COL_CAPTURE] += bg
            end
        end

        # Clamp elastic to a small positive value
        if values[i, _COL_ELASTIC] <= small
            values[i, _COL_ELASTIC] = small
        end

        # Round to 7 significant figures (matching NJOY's emerge:sigfig call)
        for j in _COL_ELASTIC:_COL_CAPTURE
            values[i, j] = round_sigfig(values[i, j], 7, 0)
        end

        # Recompute total
        values[i, _COL_TOTAL] = values[i, _COL_ELASTIC] +
                                 values[i, _COL_FISSION] +
                                 values[i, _COL_CAPTURE]
    end
end

# Extract resonance range boundaries from MF2 data
function _resonance_bounds(mf2::MF2Data)
    eresl = Inf
    eresh = 0.0
    eresr = 0.0
    for iso in mf2.isotopes
        for rng in iso.ranges
            Int(rng.LRU) == 0 && continue
            eresl = min(eresl, rng.EL)
            eresh = max(eresh, rng.EH)
            if Int(rng.LRU) == 1
                eresr = max(eresr, rng.EH)
            end
        end
    end
    eresr = clamp(eresr, eresl, eresh)
    return (eresl, eresh, eresr)
end

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
# MF3 reader
# ==========================================================================

"""
    read_mf3_sections(io::IO, mat::Integer) -> Vector{MF3Section}

Read all MF3 cross section sections for material `mat`.
"""
function read_mf3_sections(io::IO, mat::Integer)
    sections = MF3Section[]
    seekstart(io)

    while !eof(io)
        pos = position(io)
        line = readline(io)
        p = rpad(line, 80)
        mf = _parse_int(p[71:72])
        mt = _parse_int(p[73:75])
        mat_line = _parse_int(p[67:70])

        mat_line != mat && continue

        if mf == 3 && mt > 0
            seek(io, pos)
            try
                head = read_cont(io)
                tab1 = read_tab1(io)
                tf = TabulatedFunction(tab1)
                push!(sections, MF3Section(Int32(mt), head.C1, tab1.C2, tf))
                # Skip to SEND
                _skip_to_send(io)
            catch e
                @warn "read_mf3_sections: skipping MF3/MT=$mt due to parse error" exception=(e, catch_backtrace())
                _skip_to_send(io)
            end
        end
    end
    return sections
end

function _skip_to_send(io::IO)
    while !eof(io)
        line = readline(io)
        p = rpad(line, 80)
        mt = _parse_int(p[73:75])
        mt == 0 && break
    end
end

function _detect_mat(io::IO, requested::Integer)
    requested > 0 && return Int32(requested)
    seekstart(io)
    while !eof(io)
        line = readline(io)
        p = rpad(line, 80)
        mf_val = _parse_int(p[71:72])
        mat_val = _parse_int(p[67:70])
        if mat_val > 0 && mf_val > 0
            return Int32(mat_val)
        end
    end
    error("reconstruct: could not detect MAT number from file")
end

# _parse_int is already defined in endf/io.jl -- no need to redefine here

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

"""
    sigma_mf2(E, mf2) -> CrossSections

Evaluate resonance cross sections at energy E by summing over all isotopes.
Legacy interface returning CrossSections struct.
"""
function sigma_mf2(E::Real, mf2::MF2Data)
    sig = CrossSections()
    E_f = Float64(E)
    for iso in mf2.isotopes
        abn = iso.ABN
        for rng in iso.ranges
            if E_f >= rng.EL && E_f < rng.EH && Int(rng.LRU) > 0
                try
                    sigp = cross_section(E_f, rng)
                    total = max(0.0, sigp.total)
                    elastic = max(0.0, sigp.elastic)
                    fission = max(0.0, sigp.fission)
                    capture = max(0.0, sigp.capture)
                    sig = sig + abn * CrossSections(total, elastic, fission, capture)
                catch e
                    @warn "sigma_mf2: cross_section failed at E=$E_f for range [$(rng.EL), $(rng.EH)]" exception=(e, catch_backtrace())
                end
            end
        end
    end
    return sig
end

"""
    merge_background_legacy(energies, res_xs, mf3_sections) -> Vector{CrossSections}

Legacy merge returning Vector{CrossSections}.
"""
function merge_background_legacy(energies::Vector{Float64},
                                  res_xs::Vector{CrossSections},
                                  mf3_sections::Vector{MF3Section})
    n = length(energies)
    result = Vector{CrossSections}(undef, n)
    for i in 1:n
        e = energies[i]
        elastic = res_xs[i].elastic
        fission = res_xs[i].fission
        capture = res_xs[i].capture
        for sec in mf3_sections
            mt = Int(sec.mt)
            (mt == 1 || mt == 3 || mt == 101) && continue
            bg = interpolate(sec.tab, e)
            bg == 0.0 && continue
            if mt == 2
                elastic += bg
            elseif mt == 18 || mt == 19
                fission += bg
            elseif mt == 102
                capture += bg
            end
        end
        if elastic <= 1.0e-8
            elastic = 1.0e-8
        end
        total = elastic + fission + capture
        result[i] = CrossSections(total, elastic, fission, capture)
    end
    return result
end
