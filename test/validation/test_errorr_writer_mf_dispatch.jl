# RED→GREEN test for `_write_errorr_tape` body-MF / MT dispatch.
#
# Mirrors Fortran `covout` (errorr.f90:7018-7787) — the writer must remap
# the body MF column based on `mfcov`, NOT emit it literally:
#
#   mfcov=31 (ν̄)     → body MF = 33  (errorr.f90:7214 + 7472)
#   mfcov=33 (XS)     → body MF = 33  (no remap)
#   mfcov=34 (mubar)  → body MF = 34, but MT collapsed to 251
#                       (errorr.f90:7245-7250 + 7480-7485 + 7585)
#                       Per-pair info encoded as scr(3)=251, scr(4)=ld,
#                       scr(5)=ld1 in the sub-section CONT.
#   mfcov=35 (spec)   → body MF = 35  (no remap)
#   mfcov=40 (act)    → body MF = 40  (no remap)
#
# Plus: errorr's MF=3 echo must include MT=1 when present in the input
# (Fortran `colaps` does NOT filter MT=1; only `rdgout` synthesises a
# missing MT=1 from partials, which is unrelated to writer-side echo).
#
# The reader (`covard`, src/processing/covr_io.jl:231-289, mirroring
# covr.f90:820-823) sets:
#   mfflg=-14            → mf3x=40
#   mfflg=-12            → mf3x=35
#   mfflg=-11, mt=251    → mf3x=34
#   mfflg=-11, otherwise → mf3x=33
# and looks up `tape.sections[(mf3x, mt)]`. The current writer emits
# body sections keyed under literal mfcov, so for mfcov∈{31,34} covard
# fails with `MF$mf3x/MT=$mt not present`.
#
# Run:  julia --project=. test/validation/test_errorr_writer_mf_dispatch.jl

using Test
using NJOY
using LinearAlgebra: I

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const _MAT = 9237
const _ZA  = 92238.0
const _AWR = 235.984
const _EGN = [1.0e-5, 0.0253, 1.0, 100.0, 1.0e3, 1.0e6, 2.0e7]   # 6 groups

"""Round-trip: build a synthetic errorr tape via `_write_errorr_tape` for
`mfcov`, parse with `read_errorr_tape`, return the parsed `ErrorrTape`."""
function _writer_roundtrip(mfcov::Int, reaction_mts::Vector{Int};
                           pairs::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[],
                           extra_xs_mts::Vector{Int} = Int[])
    ngn = length(_EGN) - 1
    group_xs = Dict{Int, Vector{Float64}}()
    for mt in reaction_mts
        group_xs[mt] = fill(2.5, ngn)
    end
    for mt in extra_xs_mts
        group_xs[mt] = fill(1.7, ngn)
    end

    cov_matrices = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    requested_pairs = isempty(pairs) ?
        [(m, m) for m in reaction_mts] : pairs
    for (m1, m2) in requested_pairs
        cov_matrices[(m1, m2)] = 0.001 * Matrix{Float64}(I, ngn, ngn)
    end

    path = tempname() * ".errorr"
    open(path, "w") do io
        NJOY._write_errorr_tape(io, _MAT, _ZA, _AWR, _EGN,
                                group_xs, cov_matrices, reaction_mts, mfcov)
    end
    return NJOY.read_errorr_tape(path; mat=_MAT)
end

@testset "errorr writer — body MF/MT dispatch (Fortran covout)" begin

    @testset "mfcov=31 (ν̄) → body MF=33 (errorr.f90:7214,7472)" begin
        tape = _writer_roundtrip(31, [452, 455, 456])
        @test tape.mfflg == -11
        @test haskey(tape.sections, (33, 452))   # remapped to MF=33
        @test haskey(tape.sections, (33, 455))
        @test haskey(tape.sections, (33, 456))
        @test !haskey(tape.sections, (31, 452))  # NOT under literal mfcov
        @test !haskey(tape.sections, (31, 455))
        @test !haskey(tape.sections, (31, 456))

        # covard must succeed without "MF33/MT=452 not present"
        result = NJOY.covard(tape, _MAT, 452, _MAT, 452)
        @test result !== nothing
        @test result.izero >= 0
    end

    @testset "mfcov=33 (cross-section) → body MF=33 unchanged" begin
        tape = _writer_roundtrip(33, [1, 2, 4, 102])
        @test tape.mfflg == -11
        for mt in (1, 2, 4, 102)
            @test haskey(tape.sections, (33, mt))
        end
        result = NJOY.covard(tape, _MAT, 2, _MAT, 2)
        @test result !== nothing
    end

    @testset "mfcov=35 (spectrum) → body MF=35 unchanged" begin
        tape = _writer_roundtrip(35, [18])
        @test tape.mfflg == -12       # 35 → -12 sentinel
        @test haskey(tape.sections, (35, 18))
        @test !haskey(tape.sections, (33, 18))
    end

    @testset "mfcov=35 MF=3 echo limited to reaction_mts (T34 regression)" begin
        # T34 regression: when GENDF readback puts every MF=3 MT into
        # `group_xs` (not just the cov reactions), emitting MF=3 for ALL
        # of them inflates the downstream covr boxer-format output by
        # 100×+. Only reaction_mts (and MT=251 for mubar) should appear
        # in MF=3 — every other group_xs entry stays out.
        tape = _writer_roundtrip(35, [18];
                                 extra_xs_mts=[1, 2, 4, 16, 102, 251])
        @test haskey(tape.mf3_xs, 18)
        for mt_extra in (1, 2, 4, 16, 102, 251)
            @test !haskey(tape.mf3_xs, mt_extra)
        end
    end

    @testset "mfcov=33 MF=3 echo limited to reaction_mts" begin
        # Same guard for mfcov=33: extra group_xs entries (e.g. nubar
        # MTs from the GENDF) must not bleed into MF=3 unless they were
        # in the MF=33 cov reactions list.
        tape = _writer_roundtrip(33, [1, 2, 102];
                                 extra_xs_mts=[251, 452, 455, 456])
        for mt_in in (1, 2, 102)
            @test haskey(tape.mf3_xs, mt_in)
        end
        for mt_out in (251, 452, 455, 456)
            @test !haskey(tape.mf3_xs, mt_out)
        end
    end

    @testset "mfcov=40 (activation) → body MF=40 unchanged" begin
        tape = _writer_roundtrip(40, [102])
        @test tape.mfflg == -14       # 40 → -14 sentinel
        @test haskey(tape.sections, (40, 102))
        @test !haskey(tape.sections, (33, 102))
    end

    @testset "mfcov=34 (mubar) → body MF=34/MT=251 (errorr.f90:7246+7481)" begin
        # Mubar input has per-reaction MTs (e.g. T65: just MT=2; T15 mubar
        # tape: MT=2/4/16/...). Fortran covout collapses every reaction
        # into MF=34/MT=251 sections — sub-section CONT carries L1=251,
        # L2=ld, N1=ld1, with per-pair differentiation via Legendre
        # order. covard with (mt=251, mt1=251) finds the first match.
        # Writer must mirror this and ALSO emit MF=3/MT=251 (musigc,
        # errorr.f90:5897-6036) so covard's sandwich-rule lookup works.
        tape = _writer_roundtrip(34, [2]; extra_xs_mts=[251])

        @test tape.mfflg == -11       # mfcov=34 keeps -11 sentinel
        @test haskey(tape.sections, (34, 251))     # mubar lives at (34, 251)
        @test !haskey(tape.sections, (34, 2))      # per-reaction MT NOT used
        @test !haskey(tape.sections, (33, 2))      # NOT under MF=33

        # MF=3/MT=251 must be present so covard's xs lookup succeeds:
        # `covard` calls `tape.mf3_xs[251]` for relative-cov normalisation
        # (covr.f90:786-812 / covr_io.jl:255-265).
        @test haskey(tape.mf3_xs, 251)

        # `covard(mat, 251, mat, 251)` must succeed: `mf3x = 34` because
        # mt=251, sub-section L1=251 matches mt1=251 caller (covr_io.jl
        # :332, mirroring covr.f90:849).
        result = NJOY.covard(tape, _MAT, 251, _MAT, 251)
        @test result !== nothing
    end
end

@testset "errorr writer — MF=3 echo must NOT filter MT=1 (T16 regression)" begin
    # T16 crash: `covard: MF3/MT=1 missing for MAT=9237`. PENDF MF=3 path
    # in `_errorr_group_average` (errorr.jl:364) skipped MT=1 with
    # `mt in (1, 451) && continue`. Fortran `colaps` (errorr.f90:9097+)
    # does NOT filter MT=1. The skip is wrong for any caller that needs
    # cov involving MT=1 (e.g. T15 MF=33 cov pair (1, 1) and T16's full
    # U-238 MF=33 cov set).
    tape = _writer_roundtrip(33, [1, 2, 102])
    @test haskey(tape.mf3_xs, 1)      # MT=1 echoed
    @test haskey(tape.mf3_xs, 2)
    @test haskey(tape.mf3_xs, 102)

    # And covard must work for the (1, 1) self-cov pair.
    result = NJOY.covard(tape, _MAT, 1, _MAT, 1)
    @test result !== nothing
end

@testset "_errorr_group_average no longer filters MT=1 (T16 PENDF path)" begin
    # The synthetic round-trip above passes group_xs directly — it
    # cannot exercise the actual T16 bug location. T16 takes the PENDF
    # path (`_errorr_group_average`, errorr.jl:356-383) because its
    # input deck has no groupr (npend>0, ngout=0). The pre-fix filter
    # `mt in (1, 451) && continue` silently dropped MT=1 from the
    # returned dict, so MF=3/MT=1 was missing from tape26 and covard
    # crashed with `MF3/MT=1 missing`. We exercise the function directly
    # against the T15/T16 cached broadr PENDF (same material, same
    # broadr params).
    cache = joinpath(REPO, "test", "validation", "oracle_cache",
                     "test15", "after_broadr.pendf")
    if isfile(cache)
        # T15/T16 use ign=3 (LANL 30-group). Use the standard structure
        # boundaries (any valid grid will do — we only assert MT=1
        # presence in the returned dict, not the values).
        egn = [1.0e-5, 1.0, 1.0e3, 1.0e6, 2.0e7]
        result = NJOY._errorr_group_average(cache, 9237, egn, 1)
        @test haskey(result, 1)        # was filtered pre-fix → T16 crash
        @test !haskey(result, 451)     # MT=451 (directory) still skipped
    else
        @info "skipping T16 PENDF path test — cache $cache absent"
    end
end
