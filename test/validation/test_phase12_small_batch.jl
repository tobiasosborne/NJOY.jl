# Phase 12 RED→GREEN tests — small-fix batch:
#   T06 plotr dispatch stub
#   T20 errorr 999-mode cp stub
#   T43 broadr T=0 pass-through
#
# Each test asserts run_ok (no crash). Content still diffs — these are
# structural stubs/guards, not bit-identical fixes.
#
# T15/T17 were originally grouped into this batch as BoundsError fixes
# (per T10 §Recommendations item 2), but the failure mode shifted
# between the sweep report and this session — broadr now runs to
# completion (471s) then trips the harness hard-timeout. Tracked as a
# separate performance bead rather than a crash fix.
#
# Run:  julia --project=. test/validation/test_phase12_small_batch.jl

using Test
using NJOY

include(joinpath(@__DIR__, "reference_test.jl"))

const TIMEOUT = 600.0

function _quiet(n::Int)
    run_reference_test(n; verbose=false, timeout_hard_sec=TIMEOUT,
                       tolerances=[1e-9])
end

@testset "Phase 12 small-fix batch" begin
    @testset "T06 — plotr stub" begin
        r = _quiet(6)
        @test r.run_ok
        @test !occursin("No such file", r.run_error)
    end

    @testset "T20 — errorr 999-mode stub" begin
        r = _quiet(20)
        @test r.run_ok
        @test !occursin("tape unit 0 is invalid", r.run_error)
        @test !occursin("tape21", r.run_error) || !occursin("No such file", r.run_error)
    end

    @testset "T43 — broadr T=0 pass-through" begin
        r = _quiet(43)
        @test r.run_ok
        @test !occursin("InexactError", r.run_error)
        @test !occursin("NaN", r.run_error)
    end
end
