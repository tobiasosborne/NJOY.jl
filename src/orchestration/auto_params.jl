# auto_params.jl -- Auto-computation of module parameters from ENDF/PENDF data
#
# Replaces hardcoded values in test pipeline scripts. Every parameter is
# derived from the input deck + ENDF data, matching how Fortran NJOY computes them.

# =========================================================================
# Broadr: thnmax computation
# =========================================================================

"""
    compute_thnmax(mf3_sections, awr::Float64) -> Float64

Compute the maximum energy for broadening/thinning from MF3 section thresholds.
Matches Fortran broadr.f90 logic: thnmax = lowest threshold of any reaction
with negative QI (inelastic).

Returns `Inf` if no inelastic reactions found (i.e., broaden everything).
"""
function compute_thnmax(mf3_sections, awr::Float64)
    thnmax = Inf
    for sec in mf3_sections
        sec.QI >= 0 && continue  # skip non-threshold reactions
        # Physical threshold: E_thr = -QI * (AWR + 1) / AWR
        s_awr = sec.awr > 0 ? Float64(sec.awr) : awr
        thrx = -Float64(sec.QI) * (s_awr + 1) / s_awr
        thrx > 0 && (thnmax = min(thnmax, thrx))
    end
    thnmax
end

"""
    resolve_thnmax(user_thnmax::Float64, mf3_sections, awr::Float64) -> Float64

Apply the Fortran broadr convention for the user-supplied thnmax value:
- thnmax == 0: auto-compute from reaction thresholds
- thnmax < 0: use abs(thnmax) directly
- thnmax > 0: use user value directly
"""
function resolve_thnmax(user_thnmax::Float64, mf3_sections, awr::Float64)
    if user_thnmax == 0.0
        return compute_thnmax(mf3_sections, awr)
    elseif user_thnmax < 0.0
        return abs(user_thnmax)
    else
        return user_thnmax
    end
end

# =========================================================================
# Broadr: partials selection
# =========================================================================

"""
    select_broadr_partials(mf3_data::Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}})
        -> (xs_matrix, mt_ids)

Select which MF3 reactions to use for broadn_grid convergence testing.
Matches Fortran broadr nreac: uses non-redundant primary partials only
(elastic, fission, capture), NOT the total.

Returns (Matrix{Float64} with one column per partial, Vector{Int} of MT numbers).
"""
function select_broadr_partials(mf3_data::Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}})
    # Primary partials in Fortran order
    primary_mts = Int[]
    haskey(mf3_data, 2)   && push!(primary_mts, 2)    # elastic
    haskey(mf3_data, 18)  && push!(primary_mts, 18)   # fission
    haskey(mf3_data, 102) && push!(primary_mts, 102)  # capture

    isempty(primary_mts) && error("No broadr partials found (need MT=2, 18, or 102)")

    # All partials must share the same energy grid (reconr grid)
    ref_e = mf3_data[primary_mts[1]][1]
    n = length(ref_e)

    xs_matrix = Matrix{Float64}(undef, n, length(primary_mts))
    for (col, mt) in enumerate(primary_mts)
        e, xs = mf3_data[mt]
        if length(e) != n
            error("MT=$mt energy grid length $(length(e)) != reference $(n)")
        end
        xs_matrix[:, col] = xs
    end

    (xs_matrix, primary_mts)
end

# =========================================================================
# Thermr: Bragg lattice parameters
# =========================================================================

"""
Bragg lattice parameters for coherent elastic scattering.
Matches Fortran thermr.f90 `sigcoh` select-case over material numbers.

Fields:
- `a`, `c`: lattice constants [cm]
- `sigma_coh`: coherent scattering cross section [barn]
- `A_mass`: atomic mass [amu]
- `lat`: lattice type (1=graphite, 2=beryllium, 3=BeO, 4=FCC, 5=BCC, 6=HCP)
"""
const BRAGG_LATTICE_PARAMS = Dict{Int, NamedTuple{(:a, :c, :sigma_coh, :A_mass, :lat),
                                                   NTuple{5, Float64}}}(
    # Graphite (MAT 1065)
    1065 => (a=2.4573e-8, c=6.7e-8, sigma_coh=5.50, A_mass=12.011, lat=1.0),
    # Beryllium (MAT 1064)
    1064 => (a=2.2856e-8, c=3.5832e-8, sigma_coh=7.63, A_mass=9.0122, lat=2.0),
    # Beryllium oxide - Be part (MAT 1066)
    1066 => (a=2.6950e-8, c=4.3900e-8, sigma_coh=7.63, A_mass=9.0122, lat=3.0),
    # Beryllium oxide - O part (MAT 1067)
    1067 => (a=2.6950e-8, c=4.3900e-8, sigma_coh=4.232, A_mass=15.999, lat=3.0),
    # Aluminum (MAT 1097, FCC a=4.0495 Angstrom)
    1097 => (a=4.0495e-8, c=4.0495e-8, sigma_coh=1.495, A_mass=26.982, lat=4.0),
    # Iron (MAT 1028, BCC a=2.8665 Angstrom)
    1028 => (a=2.8665e-8, c=2.8665e-8, sigma_coh=11.22, A_mass=55.845, lat=5.0),
    # Silicon dioxide - Si part (MAT 1070, alpha quartz)
    1070 => (a=4.9134e-8, c=5.4052e-8, sigma_coh=2.163, A_mass=28.086, lat=6.0),
    # Silicon dioxide - O part (MAT 1071)
    1071 => (a=4.9134e-8, c=5.4052e-8, sigma_coh=4.232, A_mass=15.999, lat=6.0),
)

"""
    lookup_bragg_params(mat_thermal::Int) -> NamedTuple

Look up Bragg lattice parameters for a thermal scattering material.
Matches Fortran sigcoh select-case.
"""
function lookup_bragg_params(mat_thermal::Int)
    haskey(BRAGG_LATTICE_PARAMS, mat_thermal) ||
        error("No Bragg lattice parameters for MAT=$mat_thermal. " *
              "Add to BRAGG_LATTICE_PARAMS in auto_params.jl.")
    BRAGG_LATTICE_PARAMS[mat_thermal]
end

# =========================================================================
# General: Z extraction from ZA
# =========================================================================

"""
    extract_Z(za::Real) -> Int

Extract atomic number Z from ENDF ZA identifier (ZA = Z*1000 + A).
"""
extract_Z(za::Real) = div(round(Int, za), 1000)
