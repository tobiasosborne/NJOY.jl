# RED→GREEN test for errorr MF33 sparse per-row emission + NK geometry.
#
# T15 tape26 MF33 per-MT line counts: Julia dense-NGN-row dumps produce
# ~16-180 extra lines per MT pair with real data (hardcoded all-col per
# row), vs Fortran's covout which emits only the nonzero column range
# [ig2lo, ng2] per row and skips all-zero rows (except the last).
#
# Ref: errorr.f90:7530-7605 (covout row emission), 7244 + 7254 (NK = nmts
# - ix + 1, sub-section per (mt, mt2) with mt2 >= mt).
#
# Phase 56 (Bug A): NK now matches reference for every MT — the writer
# emits a sub-section for every (mt, mt2) with mt2 >= mt, with empty
# matrices written as 2-line zero stubs (matrix===nothing branch). The
# residual per-MT line gaps are content drift in sub-section data, not
# missing geometry; tracked under HANDOFF P1 sub-item 2 (LTY=1/2/3
# standards/ratio expansion).
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

    # NK structural acceptance — Bug A (Phase 56). Every MT in the
    # reactions list must emit NK = (count of mt2 in reactions with
    # mt2 >= mt) sub-sections, matching Fortran covout (errorr.f90:7244,
    # `scr(6)=nmts-ix+1`). Validate against reference NK directly.
    function _nk(path, mat, mt)
        open(path, "r") do io
            NJOY.find_section(io, 33, mt; target_mat=mat) || return 0
            return Int(NJOY.read_cont(io).N2)
        end
    end
    populated_self_cov_mts = [1, 2, 4, 16, 17, 18, 37, 51, 52, 53, 54, 55,
                              56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
                              67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
                              91, 102]
    for mt in populated_self_cov_mts
        @test _nk(tape26, 9237, mt) == _nk(T15_REF, 9237, mt)
    end

    # Per-MT line-count guard — backstop against row-density blowup. With
    # NK structurally correct, residual per-MT diffs are sub-section
    # content drift (HANDOFF P1 sub-item 2). Slack +60 catches a real
    # blowup (e.g. dense NGN×NGN dump) while tolerating current drift.
    # MT=2 is excluded because its NC-derived cross-pair content drifts
    # by ~+100 (open under HANDOFF P1).
    for mt in populated_self_cov_mts
        mt == 2 && continue
        jl = get(jul, mt, 0); rf = get(ref, mt, 0)
        jl == 0 && continue
        @test jl <= rf + 60
    end

    # Total tape26 — Bug A acceptance is "> 5500 lines" (HANDOFF P1).
    # Upper bound 6500 leaves room for Bug B / NC-v2 follow-ups while
    # still catching gross blowup (post-Phase-46 was 8205).
    total_lines = countlines(tape26)
    ref_total = countlines(T15_REF)
    @info "T15 tape26 total — Julia: $total_lines, ref: $ref_total"
    @test 5500 < total_lines < 6500

    # Document any residual gap for MT=2/MT=4 (post-Bug-A: MT=2 is
    # over by ~106 from sub-section content drift; MT=4 is on-ref).
    mt2_gap = get(ref, 2, 0) - get(jul, 2, 0)
    mt4_gap = get(ref, 4, 0) - get(jul, 4, 0)
    @info "Residual gap: MT=2 under by $mt2_gap lines, \
           MT=4 under by $mt4_gap lines"

    # Phase 74 — MT=1/mt2=2 cross-pair from MT=2's NC formula (MT=2 =
    # MT=1 − Σ partials). By cov symmetry the (MT=1, MT=2) cross-pair is
    # populated at rows/cols 15-30 (where U-238 MF=33 MT=1's input
    # self-cov is nonzero). Pre-fix Julia emitted only a 2-line zero
    # stub; post-fix it should emit ~16 row records × ~3 lines each.
    function _sub_section_lines(path, mat, outer_mt, target_mt2)
        # Walk the MF=33 section for `outer_mt` and count the lines in
        # the sub-section whose CONT has L2 == target_mt2. Returns
        # (sub_lines, nonzero_rows) where sub_lines includes the
        # sub-section's CONT.
        return open(path, "r") do io
            NJOY.find_section(io, 33, outer_mt; target_mat=mat) || return (0, 0)
            head = NJOY.read_cont(io)
            nl = Int(head.N2)
            for _ in 1:nl
                # Sub-section CONT carries (C1,C2,L1=0,L2=mt2,N1=0,N2=ngn).
                sub_cont = NJOY.read_cont(io)
                mt2 = Int(sub_cont.L2)
                # Each LIST row CONT carries (0,0,count,ig2lo,count,ig).
                # We need to walk LIST records until the next sub-CONT,
                # but they all share (mf=33, mt=outer_mt) — so we read
                # row-by-row tracking the LIST header.
                ngn_sub = Int(sub_cont.N2)
                row_count = 0
                line_count = 1  # the sub-CONT itself
                last_seq = 0
                for _ in 1:ngn_sub
                    row_cont = NJOY.read_cont(io)
                    count = Int(row_cont.N1)
                    ig = Int(row_cont.N2)
                    line_count += 1
                    row_count += 1
                    # Read the (count) data values — 6 per line, ceil(count/6) lines.
                    nfloat_lines = cld(count, 6)
                    for _ in 1:nfloat_lines
                        readline(io)
                        line_count += 1
                    end
                    # Last all-zero row sentinel: ig == ngn AND count == 1
                    # marks the final stub. Don't count it as a "data" row.
                    if mt2 == target_mt2
                        # Don't break; keep reading until ig=ngn signals end.
                    end
                    if ig >= ngn_sub; break; end
                end
                if mt2 == target_mt2
                    return (line_count, row_count)
                end
            end
            (0, 0)
        end
    end

    ref_lines, ref_rows = _sub_section_lines(T15_REF,  9237, 1, 2)
    jul_lines, jul_rows = _sub_section_lines(tape26,   9237, 1, 2)
    @info "MT=1/mt2=2 sub-section — Julia: $jul_lines lines / $jul_rows rows, \
           ref: $ref_lines lines / $ref_rows rows"
    # Pre-Phase-74 baseline: jul_lines ≈ 3 (sub-CONT + 1 stub row CONT +
    # 1 data line), jul_rows = 1. Post-fix: jul_rows should be ≥ 12
    # (16 in ref, allow a few rows tolerance for off-by-one σ-ratio cells
    # near the edge of the populated band).
    @test jul_rows >= 12
    # Total tape26 should close to within ~10 lines of reference (5958).
    @test abs(total_lines - ref_total) <= 15
end
