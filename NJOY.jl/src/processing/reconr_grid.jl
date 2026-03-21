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
