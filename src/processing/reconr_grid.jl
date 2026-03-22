# RECONR grid builder -- construct initial energy grid from MF2/MF3 data
#
# Step 2 of the RECONR pipeline: build_grid creates the initial energy grid
# for adaptive reconstruction by taking the union of MF2 resonance nodes
# and MF3 energy breakpoints.

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
    # matching Fortran lunion's 1/v linearization behaviour.
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
# 1/v linearization -- bisect step-ratio gaps and force decade points
# ==========================================================================

"""
    linearize_one_over_v!(energies, err; elim=1e6)

Refine `energies` so that the step ratio between consecutive points is small
enough to represent 1/v behaviour to fractional tolerance `err`.

This matches the logic of Fortran RECONR's `lunion` subroutine
(reconr.f90:1771-2238):

1. Force decade points (1, 2, 5 x 10^n) and the thermal point 0.0253 eV.
2. Walk through sorted energies; where `e[i+1]/e[i] > stpmax`, insert the
   geometric midpoint rounded to 7 significant figures.
3. Apply 5x tighter tolerance below the thermal range (0.4999 eV), matching
   Fortran's `trange` logic (line 2150).
4. Only linearize below `elim` (default 1 MeV), matching Fortran's `elim`
   (line 1811).

The vector is sorted and deduplicated on exit.
"""
function linearize_one_over_v!(energies::Vector{Float64}, err::Float64;
                               elim::Float64 = 1.0e6)
    _inject_decade_and_thermal!(energies, elim)
    length(energies) < 2 && return energies

    # Fortran constants (reconr.f90:1821, 2150)
    trange     = 0.4999
    stpmax     = 1.0 + sqrt(5.3 * err)
    stpmax_lo  = 1.0 + sqrt(5.3 * err / 5.0)  # 5x tighter below trange

    # Single-pass scan: inserting a midpoint re-checks the left sub-interval
    # before advancing, so one pass suffices.
    i = 1
    while i < length(energies)
        lo, hi = energies[i], energies[i + 1]
        if lo >= elim
            i += 1; continue
        end
        hi_eff = min(hi, elim)
        threshold = lo < trange ? stpmax_lo : stpmax
        if hi_eff / lo > threshold
            mid = round_sigfig(sqrt(lo * hi_eff), 7)
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
# MF3 linearization -- add midpoints for non-linear interpolation
# ==========================================================================

"""
    _linearize_mf3!(nodes, sec, tol)

For MF3 sections with non-linear interpolation (log-log, lin-log, etc.),
add midpoints between existing breakpoints until the cross section is
within `tol` of a linear interpolation. This ensures the 1/v shape is
represented in the initial grid, matching Fortran lunion's behaviour.
"""
function _linearize_mf3!(nodes::Vector{Float64}, sec::MF3Section, tol::Float64)
    tab = sec.tab
    npts = length(tab.x)
    npts < 2 && return

    # Check if any interpolation region uses a non-linear law
    has_nonlinear = false
    for law in tab.interp.law
        if law != LinLin && law != Histogram
            has_nonlinear = true
            break
        end
    end
    has_nonlinear || return

    # Process each interval, adding midpoints where needed
    max_depth = 20  # prevent infinite recursion
    stack = Tuple{Float64, Float64, Float64, Float64, Int}[]

    for k in 1:(npts - 1)
        x1 = tab.x[k]
        x2 = tab.x[k + 1]
        x2 <= x1 && continue

        # Determine the interpolation law for this interval
        law = _law_for_interval(tab, k)
        (law == LinLin || law == Histogram) && continue

        y1 = tab.y[k]
        y2 = tab.y[k + 1]

        # Push the interval onto the stack for iterative subdivision
        push!(stack, (x1, y1, x2, y2, 0))

        while !isempty(stack)
            xa, ya, xb, yb, depth = pop!(stack)
            depth >= max_depth && continue
            xb - xa <= xa * 1.0e-10 && continue

            xm = 0.5 * (xa + xb)
            ym_true = interpolate(tab, xm)
            ym_linear = 0.5 * (ya + yb)

            # Check if linear interpolation is within tolerance
            ref = max(abs(ym_true), abs(ym_linear), 1.0e-30)
            if abs(ym_true - ym_linear) > tol * ref
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
