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
    resolve_thnmax(user_thnmax, mf3_sections, awr; eresh=Inf) -> Float64

Apply the Fortran broadr convention for the user-supplied thnmax value:
- thnmax == 0: auto-compute (min of eresh and lowest inelastic threshold)
- thnmax < 0: use abs(thnmax) directly
- thnmax > 0: use min(user_value, eresh)

`eresh` is the resolved resonance range upper energy from MF2/MT151.
Fortran broadr.f90 lines 427-431.
"""
function resolve_thnmax(user_thnmax::Float64, mf3_sections, awr::Float64;
                        eresh::Float64=Inf)
    if user_thnmax == 0.0
        thnmax = eresh
        thresh = compute_thnmax(mf3_sections, awr)
        return min(thnmax, thresh)
    elseif user_thnmax < 0.0
        return abs(user_thnmax)
    else
        return min(user_thnmax, eresh)
    end
end

# =========================================================================
# Broadr: partials selection
# =========================================================================

"""
    select_broadr_partials(mf3_data, mt1_e; emin, lrp, eresh, iverf, nppmt)
        -> (xs_matrix, mt_ids)

Select which MF3 reactions broadr broadens (the `mtr[]` set) and project each
onto the MT=1 master energy grid `mt1_e`.

Mirrors the membership loop in broadr.f90:500-524. For each MF3 reaction `mth`
with first grid energy `enext = mf3_data[mth][1][1]`:

  HARD EXCLUDE (never broadened) — broadr.f90:500-512:
    - iverf≥6: 201≤mth≤599 or mth>850;  iverf<6: mth>150
    - mth ∈ {3,4,19};  46≤mth≤49
  INCLUDE — broadr.f90:516-522:
    - mth==18 (always)                                  (f90:512, "go to 170")
    - mth is an LRF7 channel partial in `nppmt`          (f90:513-519)
    - enext ≤ emin                                       (f90:521)
    - lrp==1 and enext < eresh                           (f90:522)
  Else SKIP (high-threshold reaction for lrp==0).

`emin` is 1.0 for lrp==0, `eresh` for lrp==1 (broadr.f90:356,432). Exothermic
reactions (QI>0, threshold ≈ 1e-5 eV) satisfy `enext ≤ emin` and are therefore
INCLUDED — the bug fixed by bead e5n was selecting only {2,18,102}, which
dropped the exothermic partials that dominate the thermal total for nuclides
like B-10 ({2,102,103,107,113}) and Co-58m1 ({2,51,102,103,107}).

All selected partials are broadened TOGETHER on one shared union grid
(broadr.f90:620, bfile3). MT1's grid is the union of all reaction grids, so a
lin-lin interpolation of a (linearized) partial onto `mt1_e` is exact and is the
faithful "union grid" — for {2,18,102} already on the MT1 grid it is the
identity. Non-mtr[] reactions are copied verbatim on their own grids by the
writer and are not returned here.

Returns (Matrix{Float64} with one column per partial on `mt1_e`, sorted Vector{Int}
of MT numbers in Fortran MF3 dictionary (ascending) order).
"""
function select_broadr_partials(mf3_data::Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}},
                                mt1_e::AbstractVector{Float64};
                                emin::Float64, lrp::Int, eresh::Float64,
                                iverf::Int, nppmt::Set{Int}=Set{Int}())
    selected = Int[]
    for mth in sort(collect(keys(mf3_data)))
        mth == 1 && continue  # total is never a broadened partial

        # HARD EXCLUDE (broadr.f90:500-511)
        if iverf >= 6
            (mth > 200 && mth < 600) && continue
            mth > 850 && continue
        else
            mth > 150 && continue
        end
        (mth == 3 || mth == 4 || mth == 19) && continue
        (mth >= 46 && mth <= 49) && continue

        enext = mf3_data[mth][1][1]  # first grid energy of this reaction

        # INCLUDE tests (broadr.f90:512-522)
        if mth == 18 || (mth in nppmt) || (enext <= emin) ||
           (lrp == 1 && enext < eresh)
            push!(selected, mth)
        end
        # else: high-threshold reaction → SKIP (copied verbatim by the writer)
    end

    isempty(selected) && error("No broadr partials selected (mtr[] predicate empty)")

    n = length(mt1_e)
    xs_matrix = Matrix{Float64}(undef, n, length(selected))
    for (col, mt) in enumerate(selected)
        e, xs = mf3_data[mt]
        xs_matrix[:, col] = _interp_onto_grid(mt1_e, e, xs)
    end

    (xs_matrix, selected)
end

"""
    _interp_onto_grid(target_e, src_e, src_xs) -> Vector{Float64}

Lin-lin interpolation of a linearized cross section `(src_e, src_xs)` onto the
master grid `target_e` (the MF3/MT1 union grid). Because every reaction's grid
is a subset of the union grid and the data are linearized, this is exact on
shared points and an exact in-segment evaluation elsewhere. Outside the source
range the value is zeroed below the first source energy (threshold reactions:
no cross section below threshold) and held flat above the last (mirrors gety1's
out-of-range behaviour for the broadened union grid).

When `src_e === target_e` (the common {2,18,102}-on-MT1-grid case) the result
is the identity copy of `src_xs`.
"""
function _interp_onto_grid(target_e::AbstractVector{Float64},
                           src_e::AbstractVector{Float64},
                           src_xs::AbstractVector{Float64})
    if length(src_e) == length(target_e) && src_e == target_e
        return copy(Vector{Float64}(src_xs))
    end
    out = Vector{Float64}(undef, length(target_e))
    for (i, e) in enumerate(target_e)
        if e <= src_e[1]
            out[i] = e == src_e[1] ? src_xs[1] : 0.0
        elseif e >= src_e[end]
            out[i] = src_xs[end]
        else
            idx = searchsortedfirst(src_e, e)
            if src_e[idx] == e
                out[i] = src_xs[idx]
            else
                f = (e - src_e[idx-1]) / (src_e[idx] - src_e[idx-1])
                out[i] = src_xs[idx-1] + f * (src_xs[idx] - src_xs[idx-1])
            end
        end
    end
    return out
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
