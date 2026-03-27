# RECONR grid builder -- construct initial energy grid from MF2/MF3 data
#
# Step 2 of the RECONR pipeline: build_grid creates the initial energy grid
# for adaptive reconstruction by taking the union of MF2 resonance nodes
# and MF3 energy breakpoints.
#
# The lunion_grid function matches Fortran RECONR's `lunion` subroutine
# (reconr.f90:1771-2238): for each MF3 section, panels are processed with a
# unified DFS that simultaneously forces decade points, checks ratios
# (linear sections), and checks interpolation errors (nonlinear sections).

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
function build_grid(mf2::MF2Data, mf3_sections::Vector{MF3Section};
                    tol::Float64 = 0.001)
    nodes = Float64[]
    sizehint!(nodes, 2048)

    # Add MF2 resonance nodes
    _add_mf2_nodes!(nodes, mf2)

    # Add MF3 breakpoints (filtered: skip redundant MT=1, MT=3)
    # Also linearize non-linear interpolation regions (log-log, lin-log)
    # matching Fortran lunion's 1/v linearisation behaviour.
    for sec in mf3_sections
        mt = Int(sec.mt)
        (mt == 1 || mt == 3 || mt == 101) && continue
        for e in sec.tab.x
            push!(nodes, e)
        end
        # Linearize sections with non-linear interpolation (law != LinLin)
        _linearize_mf3!(nodes, sec, tol)
    end

    # Add standard anchor points
    push!(nodes, 1.0e-5)
    push!(nodes, 0.0253)   # thermal point

    # Add 1-2-5-10 decade points within the energy range (matching
    # Fortran lunion's forced decade grid)
    _add_decade_points!(nodes)

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
    _add_bw_peaks!(nodes, params.Er, params.GT, params.NLS, el, eh)
end

function _add_peak_nodes!(nodes, params::MLBWParameters, el, eh)
    _add_bw_peaks!(nodes, params.Er, params.GT, params.NLS, el, eh)
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

# SAMMY/RML (LRF=7) peak nodes — matching Fortran rdsammy (samm.f90:1151-1177)
function _add_peak_nodes!(nodes, params::SAMMYParameters, el, eh)
    for sg in params.spin_groups
        nchan = Int(sg.nchan)
        for ires in 1:Int(sg.nres)
            er = sg.eres[ires]
            (er <= el || er >= eh) && continue
            # hw = gamgam/2 + Σ|gamma_j|/2 (samm.f90:1152-1155)
            # gamgam raw (no abs), channel widths with abs — matches Fortran exactly
            hw = sg.gamgam[ires] / 2.0
            for ich in 1:nchan
                hw += abs(sg.gamma[ich][ires]) / 2.0
            end
            _push_peak_triplet!(nodes, er, hw, el, eh)
        end
    end
end

# Fallback for unsupported formalisms
function _add_peak_nodes!(nodes, ::AbstractResonanceFormalism, el, eh)
    # No peak nodes for unsupported formalisms
end

# Shared logic for Breit-Wigner peak nodes
# Uses GT/2 (total width from ENDF) for half-width, matching Fortran rdf2bw
function _add_bw_peaks!(nodes, Er, GT, NLS, el, eh)
    for il in 1:Int(NLS)
        for ir in eachindex(Er[il])
            er = Er[il][ir]
            (er <= el || er > eh) && continue
            hw = GT[il][ir] / 2.0
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
# Decade points -- 1-2-5-10 grid matching Fortran lunion
# ==========================================================================

"""
    _add_decade_points!(nodes)

Add standard 1-2-5-10 decade energy points from 1e-5 to 2e7 eV.
These ensure that the grid covers all orders of magnitude, matching
Fortran lunion's forced decade grid (reconr.f90:2112-2174).
"""
function _add_decade_points!(nodes::Vector{Float64})
    # Decade multipliers within each order of magnitude
    mults = (1.0, 2.0, 5.0)
    # Cover from 1e-5 to 1e7 eV in decade steps
    for exp_val in -5:7
        base = exp10(exp_val)
        for m in mults
            e = m * base
            if e >= 1.0e-5 && e <= 2.0e7
                push!(nodes, e)
            end
        end
    end
    push!(nodes, 2.0e7)
end

# ==========================================================================
# lunion_grid -- unified grid builder matching Fortran reconr.f90:1771-2238
# ==========================================================================

"""
    lunion_grid(mf3_sections, err; nodes=Float64[], elim=1e6) -> Vector{Float64}

Build the union energy grid for RECONR by processing each MF3 section
through a unified bisection loop that simultaneously:
1. Forces decade points (1, 2, 5 × 10^n eV and thermal 0.0253 eV)
2. Checks energy ratio for linear interpolation (ratio ≤ stpmax)
3. Checks interpolation error for nonlinear sections (|true-linear| ≤ errn·|true|)

This matches Fortran RECONR's `lunion` subroutine (reconr.f90:1771-2238).
Each MF3 section is processed sequentially; its panels are checked with
a stack-based DFS and accepted points are added to the cumulative grid.
"""
function lunion_grid(mf3_sections::Vector{MF3Section}, err::Float64;
                     nodes::Vector{Float64} = Float64[],
                     elim::Float64 = 0.99e6,  # Fortran: reconr.f90:1797
                     eresl::Float64 = 0.0,
                     eresr::Float64 = 0.0,
                     eresh::Float64 = 0.0,
                     awr::Float64 = 0.0)  # AWR for threshold computation
    stpmax = 1.0 + sqrt(5.3 * err)

    # Cumulative grid starts from initial nodes
    grid = copy(nodes)
    sort!(grid); unique!(grid)

    # Detect redundant charged-particle groups (Fortran anlyzd lines 567-590)
    # Only check MF=3 sections for partial reactions
    _has_mt_range(lo, hi) = any(s -> Int(s.mf) == 3 && lo <= Int(s.mt) <= hi, mf3_sections)
    skip_mt103 = _has_mt_range(600, 649)
    skip_mt104 = _has_mt_range(650, 699)
    skip_mt105 = _has_mt_range(700, 749)
    skip_mt106 = _has_mt_range(750, 799)
    skip_mt107 = _has_mt_range(800, 849)

    for sec in mf3_sections
        mt = Int(sec.mt)
        # Skip redundant reactions (matching Fortran reconr.f90:1881-1901)
        # Only apply to MF=3 sections — MF=12/13 bypass these checks
        # (Fortran line 1881: if (mfh.ne.3) go to 180)
        if Int(sec.mf) == 3
            (mt == 1 || mt == 3 || mt == 101) && continue
            mt == 4 && continue   # inelastic total (redundant)
            (mt >= 251 && mt <= 300 && mt != 261) && continue
            mt == 120 && continue
            mt == 151 && continue
            # Skip MT=103-107 when partials exist (Fortran lunion lines 1884-1888)
            (mt == 103 && skip_mt103) && continue
            (mt == 104 && skip_mt104) && continue
            (mt == 105 && skip_mt105) && continue
            (mt == 106 && skip_mt106) && continue
            (mt == 107 && skip_mt107) && continue
        end

        tab = sec.tab
        npts = length(tab.x)
        npts < 2 && continue

        # Compute physical threshold for reactions with Q < 0
        # (matching Fortran lunion lines 1911-1921)
        # Fortran reads AWR from each MF3 section's HEAD record (awrx=c2h/awin).
        # Use per-section AWR when available; fall back to global awr.
        sec_awr = sec.awr > 0.0 ? sec.awr : awr
        thrx = 0.0
        qx = sec.QI
        if qx < 0.0 && mt != 2 && mt != 18 && mt != 19 && mt != 102
            thrx = sec_awr > 0.0 ? -qx * (sec_awr + 1) / sec_awr : -qx
        end

        # Threshold replacement: Fortran REPLACES the first MF3 breakpoint
        # with the kinematic threshold in-place (reconr.f90:1925-1943).
        work_x = tab.x  # Default: alias original x data
        start_k = 1

        if thrx > 0.0
            thrxx = round_sigfig(thrx, 7, +1)
            if tab.x[1] < thrxx
                work_x = copy(tab.x)
                work_x[1] = thrxx
                # Fix subsequent breakpoints now <= new first energy
                # (Fortran lines 1937-1943)
                l = 1
                while l < npts && work_x[l + 1] <= work_x[l]
                    (l > 11) && break
                    work_x[l + 1] = round_sigfig(work_x[l], 7, +1)
                    l += 1
                end
            end
        end

        # Pseudo-threshold advancement: skip leading zero-XS panels
        # for non-primary reactions (Fortran label 205, lines 1968-1976)
        if mt != 2 && mt != 18 && mt != 19 && mt != 102
            ir = start_k
            while ir < npts - 1
                if abs(tab.y[ir]) < _SSMALL && abs(tab.y[ir + 1]) < _SSMALL
                    ir += 1
                else
                    break
                end
            end
            start_k = ir
        end

        # Process this section's panels through the bisection loop.
        # Split panels at existing grid nodes (matching Fortran lunion
        # labels 210-280 which merges existing grid into panel processing).
        new_points = Float64[]
        sizehint!(new_points, 1024)

        for k in 1:(npts - 1)
            e_lo = tab.x[k]
            e_hi = tab.x[k + 1]
            e_hi <= e_lo && continue

            law = _law_for_interval(tab, k)
            is_nonlinear = (law != LinLin && law != Histogram)

            # Find existing grid points within this panel
            i_start = searchsortedfirst(grid, e_lo)
            i_end = searchsortedlast(grid, e_hi)

            # Build sub-panels by splitting at grid nodes
            panel_lo = e_lo
            panel_y_lo = tab.y[k]
            for gi in i_start:i_end
                gp = grid[gi]
                (gp <= panel_lo || gp >= e_hi) && continue
                # Sub-panel [panel_lo, gp]
                gp_y = is_nonlinear ? interpolate(tab, gp) :
                       tab.y[k] + (tab.y[k+1] - tab.y[k]) * (gp - e_lo) / (e_hi - e_lo)
                _bisect_panel!(new_points, tab, panel_lo, panel_y_lo,
                               gp, gp_y, is_nonlinear, law, stpmax, err, elim)
                panel_lo = gp
                panel_y_lo = gp_y
            end
            # Final sub-panel [panel_lo, e_hi]
            _bisect_panel!(new_points, tab, panel_lo, panel_y_lo,
                           e_hi, tab.y[k+1], is_nonlinear, law, stpmax, err, elim)
        end

        # Add section breakpoints with singularity shading.
        # Fortran lunion (lines 1930-1937): when consecutive MF3 breakpoints
        # are very close (duplicates encoding discontinuities), shade them to
        # sigfig(x, 7, -1) and sigfig(x, 7, +1) to prevent singularities.
        #
        # Histogram shading (reconr.f90:2050-2064): for sections with histogram
        # interpolation (like MF=12 photon multiplicities), each interior
        # breakpoint is replaced by a shaded pair {sigfig(E,7,-1), sigfig(E,7,+1)}.
        # This is a two-pass mechanism: first pass creates the lower shade,
        # second pass (triggered when et==enl) creates the upper shade. Any
        # existing grid point at the exact breakpoint energy is consumed.
        is_histogram = any(l == Histogram for l in tab.interp.law)
        shaded_energies = Float64[]

        for k in start_k:npts
            e = work_x[k]
            if k < npts && abs(work_x[k+1] - e) < 1.0e-9 * abs(e)
                # Duplicate pair detected. Fortran lunion (line 2092) checks
                # if the y-values differ: if same, skip shading (go to 260).
                # Only shade when y-values actually differ (real discontinuity).
                if abs(tab.y[k+1] - tab.y[k]) < 1.0e-9 * abs(tab.y[k] + 1.0e-30)
                    # Same y-values → not a real discontinuity. Just keep
                    # one copy (Fortran advances ir without shading).
                    push!(grid, e)
                else
                    # First of a duplicate pair: Fortran uses different shading
                    # for initial (label 207, line 1979: sigfig(er,7,0)) vs
                    # mid-data (label 270, line 2029: sigfig(er,7,-1)) duplicates
                    bias = (k == start_k) ? 0 : -1
                    push!(grid, round_sigfig(e, 7, bias))
                end
            elseif k > start_k && abs(e - work_x[k-1]) < 1.0e-9 * abs(e)
                # Second of a duplicate pair: check if it was a real discontinuity
                if abs(tab.y[k] - tab.y[k-1]) < 1.0e-9 * abs(tab.y[k-1] + 1.0e-30)
                    continue  # Same y-values → skip (already inserted one copy)
                end
                shaded_up = round_sigfig(e, 7, +1)
                push!(grid, shaded_up)
                # Fortran modifies breakpoints in-place (line 1980/2030:
                # scr(ibase+ir*2+1)=ernext). This creates cascaded duplicates:
                # if the shaded-up value matches the NEXT breakpoint, the next
                # iteration will detect a new duplicate. Only update work_x
                # when a cascade is actually present (next breakpoint matches).
                if k < npts && abs(work_x[k+1] - shaded_up) < 1.0e-9 * shaded_up
                    work_x[k] = shaded_up
                end
            elseif is_histogram && k > start_k && k < npts
                # Histogram interior breakpoint → shaded pair
                push!(grid, round_sigfig(e, 7, -1))
                push!(grid, round_sigfig(e, 7, +1))
                push!(shaded_energies, e)
            else
                # Cross-section coincidence check (Fortran label 220, lines
                # 1996-2003): when a breakpoint coincides with an existing
                # grid point AND has nonzero XS, shade to a pair
                # {sigfig(e,7,0), sigfig(e,7,+1)}. Fortran applies this to
                # ALL breakpoints in the section, not just the first one.
                if k == start_k && e >= round_sigfig(1.0e-5, 7, +1) && abs(tab.y[k]) > 0.0
                    # Fortran lunion label 220 (lines 1992-2005):
                    # Find eg = first old grid point >= er*(1-small)
                    gi = searchsortedfirst(grid, e * (1 - 1.0e-9))
                    if gi <= length(grid)
                        eg = abs(grid[gi])
                        # Fortran line 1996: compare er with sigfig(|eg|,7,-1)
                        eg_dn = round_sigfig(eg, 7, -1)
                        if (e - eg_dn) <= 1.0e-8 * e
                            # Coincident with existing shaded grid point — shade
                            push!(grid, round_sigfig(e, 7, 0))
                            e = round_sigfig(e, 7, +1)
                        end
                    end
                end
                push!(grid, e)
            end
        end

        # Remove exact values that were replaced by histogram shading.
        # The Fortran consumes old-grid points at these energies (label 280).
        if !isempty(shaded_energies)
            filter!(grid) do e
                !any(abs(e - s) <= 1.0e-8 * abs(s) for s in shaded_energies)
            end
        end
        append!(grid, new_points)

        sort!(grid)
        # Tolerance-based deduplication matching Fortran lunion (small=1e-9).
        # The sigfig bias creates near-duplicates that exact unique! misses.
        _dedup_tol!(grid)
        filter!(>(0.0), grid)
    end

    # Post-processing: remove interior points that coincide with resonance
    # range boundaries (matching Fortran lunion lines 2196-2226).
    # First and last points and points above emax=19 MeV are always kept.
    small = 1.0e-9
    emax = 19.0e6
    boundaries = Float64[eresl, eresr, eresh]
    filter!(>(0.0), boundaries)
    if !isempty(boundaries) && length(grid) > 2
        mask = trues(length(grid))
        for i in 2:(length(grid) - 1)
            eg = grid[i]
            eg > emax && continue  # keep points above emax
            for eb in boundaries
                if abs(eg - eb) <= small * eb
                    mask[i] = false
                    break
                end
            end
        end
        grid = grid[mask]
    end

    # Post-processing: step-ratio enforcement on the cumulative grid.
    # Fortran lunion processes each section cumulatively — panel bisection
    # sees grid points from ALL previous sections. Julia's per-section
    # _bisect_panel! only sees its own breakpoints. After combining,
    # consecutive grid points may have step ratios > stpmax that need
    # bisection. Matches Fortran lunion lines 2131-2140 (step ratio check).
    if elim > 0.0
        changed = true
        while changed
            changed = false
            i = 1
            while i < length(grid)
                xl, xh = grid[i], grid[i+1]
                if xl > 0.0 && xl <= elim && xh / xl > stpmax
                    xm = round_sigfig((xl + xh) / 2, 7)
                    if xm > xl && xm < xh
                        insert!(grid, i+1, xm)
                        changed = true
                    end
                end
                i += 1
            end
        end
    end

    grid
end

"""Remove near-duplicate energies from a sorted vector (tolerance=1e-9 relative)."""
function _dedup_tol!(v::Vector{Float64})
    isempty(v) && return v
    j = 1
    for i in 2:length(v)
        if abs(v[i] - v[j]) > 1.0e-9 * max(abs(v[j]), abs(v[i]))
            j += 1
            v[j] = v[i]
        end
    end
    resize!(v, j)
    v
end

# Fortran constants (reconr.f90:1796-1806)
const _STPMIN = 1.001
const _UP = 1.001
const _DN = 0.999
const _THERM = 0.0253
const _TRANGE = 0.4999
const _SSMALL = 1.0e-30

# Precomputed decade points: 1, 0.5, 0.2 × 10^ipwr for ipwr in -4:5
const _DECADE_POINTS = let
    pts = Float64[]
    for ipwr in -4:5
        base = 10.0^ipwr
        push!(pts, base)
        push!(pts, round_sigfig(base / 2, 1))
        push!(pts, round_sigfig(base * 2 / 10, 1))
    end
    push!(pts, _THERM)
    sort!(unique(pts))
end

"""
    _bisect_panel!(out, tab, e_lo, y_lo, e_hi, y_hi, is_nonlinear, law,
                   stpmax, err, elim)

Stack-based DFS matching Fortran lunion label 310 (reconr.f90:2106-2174).

For each sub-panel [lo, hi] on the stack:
1. Skip if ratio < stpmin (panel too narrow)
2. If lo < elim: check for forced decade points (1, 0.5, 0.2 × 10^n, thermal)
3. For linear sections (lo < elim): check ratio ≤ stpmax
4. For nonlinear sections: check |true - linear| ≤ errn·|true| + 1e-30
5. Accept or bisect with arithmetic midpoint rounded to 7 significant figures

Points are pushed into `out`; the caller adds them to the grid.
"""
function _bisect_panel!(out::Vector{Float64}, tab::TabulatedFunction,
                        e_lo::Float64, y_lo::Float64,
                        e_hi::Float64, y_hi::Float64,
                        is_nonlinear::Bool, law,
                        stpmax::Float64, err::Float64, elim::Float64)
    # Stack entries: (x_lower, y_lower, x_upper, y_upper)
    # Fortran stack has x(i)=lower, x(i-1)=upper (grows upward)
    # We use a Julia vector with explicit tuples.
    stack = NTuple{4, Float64}[(e_lo, y_lo, e_hi, y_hi)]
    max_stack = 50  # Fortran ndim=50

    while !isempty(stack)
        xl, yl, xh, yh = pop!(stack)
        length(stack) >= max_stack && continue

        # Label 310: check if panel is too narrow
        xh / xl < _STPMIN && continue

        # Force decade points if lower < elim (reconr.f90:2113-2127)
        forced = false
        if xl <= elim
            for dp in _DECADE_POINTS
                if xh > _UP * dp && xl < _DN * dp
                    xm = dp
                    ym = _eval_at(tab, is_nonlinear, xl, yl, xh, yh, xm, law)
                    push!(out, xm)
                    # Push both sub-panels (upper first so lower is processed first)
                    push!(stack, (xm, ym, xh, yh))
                    push!(stack, (xl, yl, xm, ym))
                    forced = true
                    break
                end
            end
        end
        forced && continue

        # Label 320: linear vs nonlinear check
        if !is_nonlinear
            # LINEAR (reconr.f90:2131-2140): ratio test
            xl > elim && continue  # accept above elim
            xh / xl <= stpmax && continue  # accept if ratio small enough
            xm = round_sigfig((xl + xh) / 2, 7)
            ym = _eval_at(tab, false, xl, yl, xh, yh, xm, law)
        else
            # NONLINEAR (reconr.f90:2142-2152): interpolation error test
            xm = round_sigfig((xl + xh) / 2, 7)
            ym = _eval_at(tab, true, xl, yl, xh, yh, xm, law)
            yl_lin = yl + (yh - yl) * (xm - xl) / (xh - xl)
            errn = xh < _TRANGE ? err / 5 : err
            test = errn * abs(ym) + _SSMALL
            abs(ym - yl_lin) <= test && continue  # accept
        end

        # Label 330: push midpoint and both sub-panels
        xm <= xl && continue
        xm >= xh && continue
        push!(out, xm)
        push!(stack, (xm, ym, xh, yh))
        push!(stack, (xl, yl, xm, ym))
    end
end

"""Evaluate cross section at xm: true interpolation for nonlinear, linear for linear."""
function _eval_at(tab::TabulatedFunction, is_nonlinear::Bool,
                  xl::Float64, yl::Float64, xh::Float64, yh::Float64,
                  xm::Float64, law)
    if is_nonlinear
        # Use the original tabulated function to get the true value
        interpolate(tab, xm)
    else
        # Linear interpolation between panel endpoints
        yl + (yh - yl) * (xm - xl) / (xh - xl)
    end
end

# ==========================================================================
# 1/v linearization -- bisect step-ratio gaps and force decade points
# ==========================================================================

"""
    linearize_one_over_v!(energies, err; elim=1e6)

Refine `energies` so that the step ratio between consecutive points is small
enough to represent 1/v behaviour to fractional tolerance `err`.

Matches Fortran lunion's LINEAR branch (reconr.f90:2131-2140) with forced
decade points (lines 2113-2127). Uses arithmetic midpoints rounded to 7
significant figures.
"""
function linearize_one_over_v!(energies::Vector{Float64}, err::Float64;
                               elim::Float64 = 0.99e6)
    _inject_decade_and_thermal!(energies, elim)
    length(energies) < 2 && return energies

    stpmax = 1.0 + sqrt(5.3 * err)

    i = 1
    while i < length(energies)
        lo, hi = energies[i], energies[i + 1]
        if lo >= elim
            i += 1; continue
        end
        hi_eff = min(hi, elim)

        # Check for decade points first (matching Fortran label 310)
        forced = false
        for dp in _DECADE_POINTS
            if hi_eff > _UP * dp && lo < _DN * dp
                if dp > lo && dp < hi
                    insert!(energies, i + 1, dp)
                    forced = true
                    break
                end
            end
        end
        forced && continue  # re-check from same i

        if hi_eff / lo > stpmax
            mid = round_sigfig((lo + hi_eff) / 2, 7)
            if mid > lo && mid < hi
                insert!(energies, i + 1, mid)
            else
                i += 1
            end
        else
            i += 1
        end
    end
    energies
end

"""Insert decade points (1,2,5 x 10^n) and thermal (0.0253 eV) into sorted grid."""
function _inject_decade_and_thermal!(energies::Vector{Float64}, elim::Float64)
    isempty(energies) && return
    emin, emax = first(energies), min(last(energies), elim)
    if 0.0253 >= emin && 0.0253 <= emax
        push!(energies, 0.0253)
    end
    for exp_val in -5:7
        base = exp10(exp_val)
        for m in (1.0, 2.0, 5.0)
            e = m * base
            e >= emin && e <= emax && push!(energies, e)
        end
    end
    sort!(energies)
    unique!(energies)
    filter!(>(0.0), energies)
end

# ==========================================================================
# MF3 linearization -- kept for backward compatibility
# ==========================================================================

"""
    _linearize_mf3!(nodes, sec, tol)

For MF3 sections with non-linear interpolation (log-log, lin-log, etc.),
add midpoints between existing breakpoints until the cross section is
within `tol` of a linear interpolation.
"""
function _linearize_mf3!(nodes::Vector{Float64}, sec::MF3Section, tol::Float64)
    tab = sec.tab
    npts = length(tab.x)
    npts < 2 && return

    has_nonlinear = false
    for law in tab.interp.law
        if law != LinLin && law != Histogram
            has_nonlinear = true
            break
        end
    end
    has_nonlinear || return

    max_depth = 20
    stack = Tuple{Float64, Float64, Float64, Float64, Int}[]

    for k in 1:(npts - 1)
        x1 = tab.x[k]
        x2 = tab.x[k + 1]
        x2 <= x1 && continue

        law = _law_for_interval(tab, k)
        (law == LinLin || law == Histogram) && continue

        y1 = tab.y[k]
        y2 = tab.y[k + 1]

        push!(stack, (x1, y1, x2, y2, 0))

        while !isempty(stack)
            xa, ya, xb, yb, depth = pop!(stack)
            depth >= max_depth && continue
            xb - xa <= xa * 1.0e-10 && continue

            xm = round_sigfig((xa + xb) / 2, 7)
            ym_true = interpolate(tab, xm)
            ym_linear = ya + (yb - ya) * (xm - xa) / (xb - xa)

            errn = xm < 0.4999 ? tol / 5 : tol
            test = errn * abs(ym_true) + 1.0e-30
            if abs(ym_true - ym_linear) > test
                push!(nodes, xm)
                push!(stack, (xa, ya, xm, ym_true, depth + 1))
                push!(stack, (xm, ym_true, xb, yb, depth + 1))
            end
        end
    end
end

"""
    _law_for_interval(tab, idx)

Return the interpolation law for interval `idx` in the tabulated function.
"""
function _law_for_interval(tab::TabulatedFunction, idx::Int)
    for r in eachindex(tab.interp.nbt)
        idx < tab.interp.nbt[r] && return tab.interp.law[r]
    end
    return tab.interp.law[end]
end
