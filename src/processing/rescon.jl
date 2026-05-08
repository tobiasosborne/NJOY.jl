# =========================================================================
# rescon — MF=32 → MF=33 resonance-parameter covariance sandwich
#
# Mirrors Fortran covout (errorr.f90:7465) → rescon (errorr.f90:8513-8819).
# When MF=32 is present, Fortran propagates resonance-parameter
# uncertainty into MF=33 group-cov by:
#   1. resprx (errorr.f90:3011) reads MF=32 covariance into module-level
#      cff/cfg/cgg/cee/cef/ceg/ctt arrays (resolved) via rpxlc0/rpxlc12/
#      rpxsamm, and uff/.../utt (URR) via rpxunr.
#   2. The fill itself happens inside the format-specific reader, e.g.
#      rpxlc12 (errorr.f90:4541-4593) computes finite-difference
#      sensitivities ∂σ_X(ig)/∂RP_p via paired rpendf/rpxgrp calls and
#      accumulates `c**[igind] += cov[i,j]·sens(channel,i,ig)·sens(channel,j,ig2)`.
#   3. covout calls rescon (errorr.f90:8513) which simply ADDS the
#      precomputed c**[igind] entries into cova[ig, ig2] for the right
#      (mt, mt2) pair (the dispatch table at errorr.f90:8531-8539).
#
# Data flow: MF=32 cov → finite-difference sensitivities → triangular
# packed cgg/cff/.../ctt (and full-matrix cfg/cef/ceg) → cova additive.
#
# Status (Phase 72, 2026-05-08): MF=32 reader landed (read_mf32 in
# mf32_reader.jl). Sensitivity builder + cov fill is the next-session
# work. `apply_rescon!` below parses the data, validates the dispatch,
# and exits with a structured summary — it does NOT fabricate values.
# Output covariances for the seven (mt, mt2) pairs remain unchanged
# until the sensitivity port lands. The Phase-71 RED canary in
# `test_errorr_covcal_lb5.jl` stays @test_broken.
# =========================================================================

# Seven (mt, mt2) reaction pairs that receive an RP-cov contribution.
# Mirrors errorr.f90:8531-8539 (rescon itp dispatch). Order matches the
# Fortran branch order so future sensitivity-builder code can dispatch
# by index 1..7 directly.
const RESCON_PAIRS = (
    (18,  18 ),  # itp=1: cff (triangular)
    (18,  102),  # itp=2: cfg (full)
    (102, 102),  # itp=3: cgg (triangular) — U-238 MT=102 row 1..14 canary
    (2,   2  ),  # itp=4: cee (triangular)
    (2,   18 ),  # itp=5: cef (full)
    (2,   102),  # itp=6: ceg (full)
    (1,   1  ),  # itp=7: ctt (triangular)
)

"""
    rescon_pair_index(mt::Int, mt2::Int) -> Int

Map a `(mt, mt2)` pair to its rescon `itp` index 1..7, or 0 if the pair
does not receive an RP-cov contribution. Mirrors the seven-way `if`
chain at errorr.f90:8531-8539. Order is `mt ≤ mt2` (canonical).
"""
function rescon_pair_index(mt::Int, mt2::Int)
    canonical = mt <= mt2 ? (mt, mt2) : (mt2, mt)
    for (i, p) in enumerate(RESCON_PAIRS)
        canonical == p && return i
    end
    return 0
end

"""
    apply_rescon!(cov_matrices, endf_path, mat, egn, group_xs)

Add MF=32 (resonance-parameter) covariance contributions to the seven
relevant `(mt, mt2)` entries of `cov_matrices`. Mirrors Fortran covout
(errorr.f90:7465) → rescon (errorr.f90:8513-8819).

Currently this function:
  1. Reads MF=32 via `read_mf32` (the foundation landed in Phase 72).
  2. Validates the data (range count, MPAR/NRB consistency).
  3. Logs the pending-sensitivity status.
  4. Returns without modifying `cov_matrices`.

The sensitivity builder (Fortran rpendf + rpxgrp + the perturbation
loop at errorr.f90:4249-4523) is the next-session port. Until it
lands, output MF=33 covariances for the seven (mt, mt2) pairs remain
unchanged — exactly the current behaviour, but now with the data
parsed and validated upstream.

Returns the parsed `MF32Data` (or `nothing` if MF=32 is absent).
"""
function apply_rescon!(
    cov_matrices::Dict{Tuple{Int,Int},Matrix{Float64}},
    endf_path::AbstractString,
    mat::Integer,
    egn::AbstractVector{<:Real},
    group_xs::Dict{Int,Vector{Float64}},
)
    # Re-check presence — caller already gates on _mf32_present, but
    # making this self-contained means apply_rescon! is unit-testable
    # without the orchestrator.
    open(endf_path, "r") do io
        find_section(io, 32, 151; target_mat=Int(mat)) || return nothing
    end

    data = read_mf32(endf_path, mat)

    nres = mf32_resonance_count(data)
    npar = mf32_param_count(data)
    n_ranges = sum(length(iso.resolved_ranges) for iso in data.isotopes)

    @info "apply_rescon!: parsed MF=32 — MAT=$mat, $(length(data.isotopes)) \
          isotope(s), $n_ranges resolved range(s), $nres resonance(s), \
          $npar uncertain parameter(s). Sensitivity-builder port pending \
          (rpendf + rpxgrp + perturbation loop, errorr.f90:4249-4523); \
          MF=33 cov for (1,1)/(2,2)/(2,18)/(2,102)/(18,18)/(18,102)/(102,102) \
          unchanged this run."

    # Sensitivity builder + sandwich fill go here. Spec:
    #   for each isotope, range:
    #     for each subsection:
    #       sens = build_sensitivity_jacobian(range, group_xs, egn)
    #         # 4-D: sens[channel, p, ig] for channel ∈ 1:4 (tot,el,fis,cap),
    #         # p ∈ 1:npar, ig ∈ 1:ngn. Per errorr.f90:4500-4514.
    #       for each pair (mt, mt2) ∈ RESCON_PAIRS:
    #         apply_sandwich!(cov_matrices, mt, mt2, sens, subsection.cov)
    # Channel mapping (per errorr.f90:4550-4593):
    #   sens[1,*,*] = total   → ctt          (1,1)
    #   sens[2,*,*] = elastic → cee/cef/ceg  (2,2)/(2,18)/(2,102)
    #   sens[3,*,*] = fission → cff/cef/cfg  (18,18)/(2,18)/(18,102)
    #   sens[4,*,*] = capture → cgg/ceg/cfg  (102,102)/(2,102)/(18,102)
    return data
end

"""
    rescon_supports_pair(mt::Int, mt2::Int) -> Bool

True iff `(mt, mt2)` (in either order) is one of the seven pairs that
rescon contributes to. See `RESCON_PAIRS`.
"""
rescon_supports_pair(mt::Int, mt2::Int) = rescon_pair_index(mt, mt2) != 0
