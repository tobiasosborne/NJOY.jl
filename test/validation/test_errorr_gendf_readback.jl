# RED→GREEN test for errorr GENDF MF3 readback (Fortran colaps).
#
# T15/T16/T17 errorr decks run with `npend=0, ngout=91` — group-averaged
# cross sections must come from the groupr GENDF, not a PENDF. Fortran
# handles this in `colaps` (errorr.f90:9097-9532): walk MF=3 on ngout,
# emit one LIST per group with (flux, [nubar,] sigma).
#
# Pre-fix: Julia gates group_xs population on `params.npend > 0`. With
# npend=0, the dict stays empty → the MF3 writer skips every MT → tape26
# is short by ~216 lines (HANDOFF Phase 44).
#
# Post-fix: `_errorr_read_gendf_xs` reads per-group sigma (position 2 of
# each GENDF LIST body) when npend=0 and ngout>0.
#
# Scope: same-ign case — groupr ign matches errorr ign. T15 uses ign=3
# in both. Cross-ign collapse (full Fortran colaps behaviour) deferred.
#
# Fast loop: reuses the T15 post-broadr cache + a single groupr run to
# produce a real Julia GENDF, then drives errorr_module directly.
#
# Run:  julia --project=. test/validation/test_errorr_gendf_readback.jl

using Test
using NJOY

const REPO      = normpath(joinpath(@__DIR__, "..", ".."))
const T15_CACHE = joinpath(REPO, "test", "validation", "oracle_cache", "test15")
const T15_INPUT = joinpath(REPO, "njoy-reference", "tests", "15", "input")

# Auto-expanded MTs (Fortran nextr yield on T15 MF=3). The errorr MF=3
# output block should contain all of these + explicit nubar MTs.
const EXPECTED_MF3_MTS = Int[1, 2, 4, 16, 17, 18, 51, 52, 53, 54, 55, 56,
    57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
    74, 75, 76, 77, 91, 102]

"""Return sorted unique MTs for (mat, mf) on an ENDF/GENDF tape."""
function _tape_mts(path::AbstractString, mat::Int; mf::Int = 3)
    mts = Set{Int}()
    for line in eachline(path)
        length(line) < 75 && continue
        p = rpad(line, 80)
        m = tryparse(Int, strip(p[67:70])); m === nothing && continue
        f = tryparse(Int, strip(p[71:72])); f === nothing && continue
        t = tryparse(Int, strip(p[73:75])); t === nothing && continue
        m == mat && f == mf && t > 0 && push!(mts, t)
    end
    return sort!(collect(mts))
end

"""Run T15's groupr first (produces Julia GENDF tape91), then the first
T15 errorr call (mfcov=31, nubar). Returns (work_dir, tape25_path)."""
function _run_t15_errorr_mfcov31(; work_dir::AbstractString = mktempdir())
    txt = read(T15_INPUT, String)
    calls = NJOY.tokenise_njoy_input(txt)

    # 1. groupr → tape91
    gcall = calls[findfirst(c -> c.name === :groupr, calls)]
    gparams = NJOY.parse_groupr(gcall)

    tapes = NJOY.TapeManager(Dict{Int,String}(), work_dir)
    NJOY.register!(tapes, 21, joinpath(T15_CACHE, "run_broadr", "tape20"))
    NJOY.register!(tapes, 23, joinpath(T15_CACHE, "after_broadr.pendf"))
    NJOY.groupr_module(tapes, gparams)

    # 2. first errorr (mfcov=31, nubar cov → tape25)
    ecall = calls[findfirst(c -> c.name === :errorr, calls)]
    eparams = NJOY.parse_errorr(ecall)
    NJOY.errorr_module(tapes, eparams)

    return work_dir, NJOY.resolve(tapes, abs(eparams.nout))
end

"""Run T15's second errorr call (mfcov=33, XS cov → tape26).  Much bigger
output — this is the tape where the 35-MT MF3 readback is most visible."""
function _run_t15_errorr_mfcov33(; work_dir::AbstractString = mktempdir())
    txt = read(T15_INPUT, String)
    calls = NJOY.tokenise_njoy_input(txt)

    gcall = calls[findfirst(c -> c.name === :groupr, calls)]
    gparams = NJOY.parse_groupr(gcall)

    tapes = NJOY.TapeManager(Dict{Int,String}(), work_dir)
    NJOY.register!(tapes, 21, joinpath(T15_CACHE, "run_broadr", "tape20"))
    NJOY.register!(tapes, 23, joinpath(T15_CACHE, "after_broadr.pendf"))
    NJOY.groupr_module(tapes, gparams)

    # Pick the errorr call with mfcov=33 (second errorr in the deck).
    errorrs = findall(c -> c.name === :errorr, calls)
    @assert length(errorrs) >= 2 "T15 deck should have 3 errorr calls"
    ecall = calls[errorrs[2]]
    eparams = NJOY.parse_errorr(ecall)
    NJOY.errorr_module(tapes, eparams)
    return work_dir, NJOY.resolve(tapes, abs(eparams.nout))
end

@testset "errorr GENDF MF3 readback (T15)" begin
    @testset "tape26 (mfcov=33, XS cov) — auto-expanded MTs present" begin
        work_dir, tape26 = _run_t15_errorr_mfcov33()
        @test isfile(tape26)

        mts = _tape_mts(tape26, 9237; mf = 3)
        @info "T15 tape26 MF=3 MTs ($(length(mts))): $mts"

        # Every auto-expanded MT should now have an MF=3 section on tape26.
        for mt in EXPECTED_MF3_MTS
            @test mt in mts
        end
    end

    @testset "tape25 (mfcov=31, nubar cov) — MT=452 MF3 present" begin
        work_dir, tape25 = _run_t15_errorr_mfcov31()
        @test isfile(tape25)
        mts = _tape_mts(tape25, 9237; mf = 3)
        @info "T15 tape25 MF=3 MTs: $mts"
        # Nubar covariance output always has MT=452 MF3 (nubar xs), via existing
        # _read_gendf_nubar path — should NOT regress from the new logic.
        @test 452 in mts
    end
end
