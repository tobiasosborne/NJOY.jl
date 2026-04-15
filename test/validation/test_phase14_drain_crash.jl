# Phase 14 RED→GREEN tests — drain the structural CRASH bucket.
#
#   T24 heatr plot-tape stub (4th card-1 int)
#   T60 broadr MT=1 fallback for dosimetry/IRDFF evaluations
#   T09 thermr free-gas fallback when SAB tape is empty (leapr stub)
#
# Run:  julia --project=. test/validation/test_phase14_drain_crash.jl

using Test
using NJOY

include(joinpath(@__DIR__, "reference_test.jl"))

const TIMEOUT = 900.0

function _quiet(n::Int)
    run_reference_test(n; verbose=false, timeout_hard_sec=TIMEOUT,
                       tolerances=[1e-9])
end

@testset "Phase 14 drain-CRASH" begin
    @testset "T24 — heatr plot-tape stub" begin
        # Pre-P14: CRASH "tape26: No such file" at viewr (heatr plot
        # tape never written because HeatrParams had no nplot field).
        r = _quiet(24)
        @test r.run_ok
        @test !occursin("tape26", r.run_error) ||
              !occursin("No such file", r.run_error)
    end

    @testset "T60 — broadr MT=1 fallback (dosimetry/IRDFF)" begin
        # Pre-P14: CRASH "broadr: MT=1 (total) not found on input PENDF"
        # Fe-nat IRDFF-II has dosimetry-only MTs, no elastic/total.
        r = _quiet(60)
        @test r.run_ok
        @test !occursin("MT=1 (total) not found", r.run_error)
    end

    @testset "T09 — thermr SAB fallback to free-gas" begin
        # Pre-P14: CRASH "MF7/MT4 not found for MAT=101" because leapr
        # stub outputs an empty tape and thermr's SAB path needs real
        # MF7. With iinc=2 + empty SAB tape, thermr now falls back to
        # iinc=1 (free-gas) with a warning.
        r = _quiet(9)
        @test r.run_ok
        @test !occursin("MF7/MT4 not found", r.run_error)
    end
end
