# readers.jl -- Additional ENDF section readers for module orchestration
#
# These readers extract data needed by heatr and thermr that isn't
# available through the existing reconr/MF2 readers.

# =========================================================================
# MF12: Photon multiplicity data (for heatr gamma recoil)
# =========================================================================

"""
    read_mf12_gammas(filename::AbstractString, mat::Integer; mt::Integer=102)
        -> Vector{Tuple{Float64, Float64}}

Read MF12 photon yield data for a given material and MT. Returns a vector
of (E_gamma, yield) tuples, one per discrete gamma line.

The first TAB1 subsection (total yield) is skipped; per-gamma subsections
follow, each starting with E_gamma in the TAB1 C1 field.

Replaces the manual 60-line parsing block in t01_pipeline.jl.
"""
function read_mf12_gammas(filename::AbstractString, mat::Integer; mt::Integer=102)
    gammas = Tuple{Float64, Float64}[]

    open(filename, "r") do io
        found = find_section(io, 12, mt; target_mat=Int(mat))
        found || return gammas

        # HEAD record: ZA, AWR, LO, 0, NK, 0
        head = read_cont(io)
        nk = Int(head.N1)
        nk <= 0 && return gammas

        # First subsection: total yield TAB1 — skip it
        _discard_tab1(io)

        # NK per-gamma subsections follow the total
        for k in 1:nk
            tab1 = read_tab1(io)
            eg = Float64(tab1.C1)              # gamma energy [eV]
            yield_val = isempty(tab1.y) ? 0.0 : Float64(tab1.y[1])  # yield at first energy
            push!(gammas, (eg, yield_val))
        end
    end

    gammas
end

"""Skip one TAB1 record (interp params + data pairs) from an IO stream."""
function _discard_tab1(io::IO)
    read_tab1(io)  # just read and discard
    nothing
end

# =========================================================================
# MF5: Energy distribution data (for heatr evaporation)
# =========================================================================

"""
    read_mf5_evaporation(filename::AbstractString, mat::Integer; mt::Integer=91)
        -> Union{Nothing, NamedTuple{(:u, :theta, :lf), Tuple{Float64, Float64, Int}}}

Read MF5 evaporation spectrum parameters for a given MT.
Returns (u=threshold, theta=nuclear_temperature, lf=distribution_law)
or `nothing` if MF5/MT not found.

Only handles LF=9 (simple evaporation spectrum) and LF=5 (general evaporation).
For LF=9, theta is the constant nuclear temperature.
"""
function read_mf5_evaporation(filename::AbstractString, mat::Integer; mt::Integer=91)
    open(filename, "r") do io
        found = find_section(io, 5, mt; target_mat=Int(mat))
        found || return nothing

        # HEAD: ZA, AWR, 0, 0, NK, 0
        head = read_cont(io)
        nk = Int(head.N1)
        nk <= 0 && return nothing

        # First subsection: probability TAB1 with LF in L2
        prob_tab1 = read_tab1(io)
        u = Float64(prob_tab1.C1)     # threshold energy [eV]
        lf = Int(prob_tab1.L2)        # distribution law

        if lf == 9  # simple fission / evaporation spectrum
            # Next TAB1: theta(E)
            theta_tab1 = read_tab1(io)
            theta = isempty(theta_tab1.y) ? 0.0 : Float64(theta_tab1.y[1])
            return (u=u, theta=theta, lf=lf)
        elseif lf == 5  # general evaporation
            # Next TAB1: theta(E), same format
            theta_tab1 = read_tab1(io)
            theta = isempty(theta_tab1.y) ? 0.0 : Float64(theta_tab1.y[1])
            return (u=u, theta=theta, lf=lf)
        end

        # Other LF values: return what we have
        return (u=u, theta=0.0, lf=lf)
    end
end
