# RED→GREEN test for errorr MF33 sparse per-row emission.
#
# T15 tape26 MF33 per-MT line counts: Julia dense-NGN-row dumps produce
# ~16-180 extra lines per MT pair with real data (hardcoded all-col per
# row), vs Fortran's covout which emits only the nonzero column range
# [ig2lo, ng2] per row and skips all-zero rows (except the last).
#
# Ref: errorr.f90:7530-7605 (covout row emission).
#
# Pre-fix: MT=1 MF33 has 287 lines (Julia) vs 271 (ref). Delta +16.
# Post-fix: closer to ref (dense rows become variable-width).
#
# Orthogonal issue NJOY.jl-km1: MT=2/MT=4 MF33 are UNDER-emitted (Julia
# 103/100 vs ref 1400/1106) because Julia's cov_matrices skips NC-derived
# cross-MT blocks. Separate phase; this test guards only the sparse-row
# emission outcome on MTs where Julia already has populated matrices.
#
# Run:  julia --project=. test/validation/test_errorr_mf33_sparse.jl

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

@testset "errorr MF33 sparse row emission (T15 tape26)" begin
    work_dir, tape26 = _run_t15_errorr_mfcov33()
    @test isfile(tape26)

    jul = _mf_mt_line_counts(tape26, 9237, 33)
    ref = _mf_mt_line_counts(T15_REF, 9237, 33)
    total_jul = sum(values(jul)); total_ref = sum(values(ref))
    @info "T15 tape26 MF33 line counts — Julia total: $total_jul, \
           reference total: $total_ref"

    # MTs where Julia's cov_matrices is populated (self-cov present in the
    # ENDF as standalone NI blocks). Assert Julia ≤ ref + slack for
    # these — sparse emission should match or undercut the dense reference.
    #
    # Excluded: MT=2, MT=4 (under-emission from NJOY.jl-km1, NC-derived
    # cross-pair blocks not yet expanded) and MT=0 (tape delimiters).
    populated_self_cov_mts = [1, 16, 17, 18, 37, 51, 52, 53, 54, 55, 56,
                              57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
                              68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
                              91, 102]

    # Hard assertion: Julia must not emit more lines than reference for
    # these MTs (allow +15 slack to tolerate MT=77-class covcal content
    # drift; NJOY.jl-f8k tracks the underlying matrix-extent mismatch).
    for mt in populated_self_cov_mts
        jl = get(jul, mt, 0); rf = get(ref, mt, 0)
        jl == 0 && continue  # MT not emitted by Julia — skip
        @test jl <= rf + 15
    end

    # Total tape26 should drop well below post-Phase-46 baseline (8205).
    # Sparse emission alone drives this under the reference (5958) because
    # MT=2/4 cross-MT data is still missing (NJOY.jl-km1). Not asserting
    # ≥ ref — just that the raw over-expansion is gone.
    total_lines = countlines(tape26)
    ref_total = countlines(T15_REF)
    @info "T15 tape26 total — Julia: $total_lines, ref: $ref_total"
    @test total_lines < 4000  # post-Phase-46 was 8205

    # Document the remaining under-emission gap for MT=2/MT=4 (known).
    mt2_gap = get(ref, 2, 0) - get(jul, 2, 0)
    mt4_gap = get(ref, 4, 0) - get(jul, 4, 0)
    @info "Known gap (NJOY.jl-km1): MT=2 under by $mt2_gap lines, \
           MT=4 under by $mt4_gap lines"
end
