# RED→GREEN test for errorr MF33 NC-block (derived covariance) expansion.
#
# T15 tape26 MF33 MT=2 and MT=4 are NC-derived: their covariance is a
# linear combination (LTY=0) of other MTs' covariances. Previously the
# Julia errorr port read NC LIST records and discarded them
# (errorr.jl:131-134, pre-NJOY.jl-km1), so MT=2 / MT=4 sub-sections
# emitted only zero-stub rows.
#
# Reference T15 tape26 (line counts via `_mf_mt_line_counts`):
#   MT=2: 1400 lines  (NK=35: self + 34 cross-pairs vs all reaction MTs)
#   MT=4: 1106 lines  (NK=34: self + 33 cross-pairs)
# Pre-fix Julia: ~106 / ~103 (just empty stubs).
#
# Fortran reference: errorr.f90 gridd (lines 1091-1483) builds the
# `akxy` coefficient array; covout (lines 7431-7438) accumulates
# `cov_out[ig,igp] += akxy[iy,ix,k] * akxy[iyp,ixp,kp] * cov_in[jg,jgp]`.
# For the U-238 case (no input cross-covariances among referenced MTs)
# this collapses to:
#   Cov(mt, mt)        = sum_i c_i^2 * Cov(ref_i, ref_i)   over E in [E1,E2]
#   Cov(mt, ref_j)     = c_j     * Cov(ref_j, ref_j)       (for ref_j > mt)
# (Cross-pairs where BOTH endpoints are NC-derived — e.g. Cov(2,4) — are
# v2; v1 emits them as zero stubs, costing ~40 lines vs reference.)
#
# Run:  julia --project=. test/validation/test_errorr_nc_expansion.jl

using Test
using NJOY

const REPO      = normpath(joinpath(@__DIR__, "..", ".."))
const T15_CACHE = joinpath(REPO, "test", "validation", "oracle_cache", "test15")
const T15_INPUT = joinpath(REPO, "njoy-reference", "tests", "15", "input")
const T15_REF   = joinpath(REPO, "njoy-reference", "tests", "15",
                           "referenceTape26")

"""Per-MT line count for a given MAT/MF on an ENDF-format tape."""
function _mf_mt_line_counts(path::AbstractString, mat::Int, mf::Int)
    counts = Dict{Int, Int}()
    for line in eachline(path)
        length(line) < 75 && continue
        p = rpad(line, 80)
        m = tryparse(Int, strip(p[67:70])); m === nothing && continue
        f = tryparse(Int, strip(p[71:72])); f === nothing && continue
        t = tryparse(Int, strip(p[73:75])); t === nothing && continue
        m == mat && f == mf && t > 0 || continue
        counts[t] = get(counts, t, 0) + 1
    end
    return counts
end

"""
Per-(mt, mt2) sub-section line count for MAT/MF on an ENDF-format tape.
Identifies sub-section boundaries by the CONT header `0 0 MAT1 MT1 NC NI`
emitted between the per-MT HEAD record and each NI sub-subsection. The
header line is the one with `C1 == 0.0 && C2 == 0.0 && L1 == 0 && N1 == 0`
within an MF/MT block (after the leading HEAD record).
"""
function _mf_mt_mt2_line_counts(path::AbstractString, mat::Int, mf::Int)
    counts = Dict{Tuple{Int,Int}, Int}()
    cur_mt = 0
    cur_mt2 = 0
    cur_seq_after_head = 0
    for line in eachline(path)
        length(line) < 75 && continue
        p = rpad(line, 80)
        m = tryparse(Int, strip(p[67:70])); m === nothing && continue
        f = tryparse(Int, strip(p[71:72])); f === nothing && continue
        t = tryparse(Int, strip(p[73:75])); t === nothing && continue
        if m != mat || f != mf || t <= 0
            cur_mt = 0; cur_mt2 = 0; cur_seq_after_head = 0
            continue
        end
        if t != cur_mt
            cur_mt = t; cur_mt2 = 0; cur_seq_after_head = 0
        end
        cur_seq_after_head += 1
        # First record in MF/MT block is the HEAD; skip it.
        cur_seq_after_head == 1 && continue
        c1 = strip(p[1:11]); c2 = strip(p[12:22])
        l1 = tryparse(Int, strip(p[23:33]))
        n1 = tryparse(Int, strip(p[45:55]))
        is_subsection_header = (c1 == "0.000000+0" && c2 == "0.000000+0" &&
                                l1 !== nothing && l1 == 0 &&
                                n1 !== nothing && n1 == 0)
        if is_subsection_header
            l2 = tryparse(Int, strip(p[34:44]))
            l2 === nothing && continue
            cur_mt2 = l2
        end
        cur_mt2 == 0 && continue
        key = (cur_mt, cur_mt2)
        counts[key] = get(counts, key, 0) + 1
    end
    return counts
end

function _run_t15_errorr_mfcov33(; work_dir::AbstractString = mktempdir())
    txt = read(T15_INPUT, String)
    calls = NJOY.tokenise_njoy_input(txt)
    gcall = calls[findfirst(c -> c.name === :groupr, calls)]
    gparams = NJOY.parse_groupr(gcall)
    tapes = NJOY.TapeManager(Dict{Int,String}(), work_dir)
    NJOY.register!(tapes, 21, joinpath(T15_CACHE, "run_broadr", "tape20"))
    NJOY.register!(tapes, 23, joinpath(T15_CACHE, "after_broadr.pendf"))
    NJOY.groupr_module(tapes, gparams)
    errorrs = findall(c -> c.name === :errorr, calls)
    ecall = calls[errorrs[2]]
    eparams = NJOY.parse_errorr(ecall)
    NJOY.errorr_module(tapes, eparams)
    return work_dir, NJOY.resolve(tapes, abs(eparams.nout))
end

@testset "errorr MF33 NC-block expansion (T15 tape26 MT=2/MT=4)" begin
    work_dir, tape26 = _run_t15_errorr_mfcov33()
    @test isfile(tape26)

    jul = _mf_mt_line_counts(tape26, 9237, 33)
    ref = _mf_mt_line_counts(T15_REF, 9237, 33)

    # MT=2 self+cross-cov derived from NC. Reference NK=35 → 1400 lines.
    # Tolerance band accounts for the known per-MT row-sparsity drift
    # (NJOY.jl-f8k). The v2 fix (double-NC Cov(2,4)) brings line count
    # closer to reference but still leaves the f8k content drift.
    @info "T15 MF33 MT=2: Julia=$(get(jul, 2, 0)), ref=$(get(ref, 2, 0))"
    @test get(jul, 2, 0) >= 1300
    @test get(jul, 2, 0) <= get(ref, 2, 0) + 150

    # MT=4 self+cross-cov derived from NC. Reference NK=34 → 1106 lines.
    @info "T15 MF33 MT=4: Julia=$(get(jul, 4, 0)), ref=$(get(ref, 4, 0))"
    @test get(jul, 4, 0) >= 1050
    @test get(jul, 4, 0) <= get(ref, 4, 0) + 50

    # Total MF33 should grow by at least ~2300 lines from pre-fix 1556.
    # Remaining gap to reference 5655 is f8k content drift on the
    # non-NC-derived MTs (independent issue tracked separately).
    total_jul = sum(values(jul)); total_ref = sum(values(ref))
    @info "T15 MF33 total — Julia: $total_jul, ref: $total_ref"
    @test total_jul >= 3700

    # Total tape26 should grow well past pre-fix 1859. Same f8k gap
    # caps the achievable maximum below reference 5958.
    total_lines = countlines(tape26)
    ref_total = countlines(T15_REF)
    @info "T15 tape26 total — Julia: $total_lines, ref: $ref_total"
    @test total_lines >= 4000

    # ---- v2: double-NC-derived cross-pair Cov(2, 4) ----
    # Both MT=2 and MT=4 are NC-derived in T15 U-238 JENDL. Reference
    # tape26 emits a real ~69-line LB-block sub-section for (mt=2, mt2=4).
    # v1 emitted only a 3-line zero stub (writer's fallback when
    # cov_matrices has no entry for the pair). v2 implements the full
    # double sum from Fortran covout (errorr.f90:7431-7438) so the pair
    # gets real values:  Cov(2, 4) = -Σ_{iz ∈ {51..77, 91}} Cov_in(iz, iz)
    # (refs(2) ∩ refs(4) = {51..77, 91}; coefficients are -1 in MT=2's
    # formula and +1 in MT=4's, so the product is -1 per term).
    jul_pairs = _mf_mt_mt2_line_counts(tape26, 9237, 33)
    ref_pairs = _mf_mt_mt2_line_counts(T15_REF, 9237, 33)
    jul_24 = get(jul_pairs, (2, 4), 0)
    ref_24 = get(ref_pairs, (2, 4), 0)
    @info "T15 MF33 Cov(2, 4): Julia=$jul_24, ref=$ref_24"
    # Ref is 69 lines (one full LB-block sub-section). Pre-v2 Julia
    # was 3 lines (zero-stub). Real expansion should be at least 30
    # lines (depends on group-grid sparsity); upper bound allows for
    # f8k row-sparsity drift on the constituent partials.
    @test jul_24 >= 30
    @test jul_24 <= ref_24 + 50
end
