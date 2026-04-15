# Phase 11 RED→GREEN test for covr dispatch.
#
# Pre-fix: T05 CRASHes with
#   SystemError: opening file "…/tape34": No such file or directory
# because viewr runs after an undispatched `covr` and finds no tape34.
#
# Post-fix: covr_module touches tape34 so the chain reaches viewr; run_ok
# flips to true. Output content is still 0 vs 78796 lines (stub writes an
# empty file) — that's DIFFS-class, not CRASH.
#
# T16 is another covr consumer but takes ~20 min (JENDL U-238 errorr is
# slow) so it's excluded from the automated suite; verify manually via
#   julia --project=. test/validation/reference_test.jl 16
#
# Run:  julia --project=. test/validation/test_phase11_covr_dispatch.jl

using Test
using NJOY

include(joinpath(@__DIR__, "reference_test.jl"))

const TIMEOUT = 600.0

function _quiet(n::Int)
    run_reference_test(n; verbose=false, timeout_hard_sec=TIMEOUT,
                       tolerances=[1e-9])
end

@testset "Phase 11 RED→GREEN — covr dispatch" begin
    @testset "T05 — covr writes tape34 (stub)" begin
        r = _quiet(5)
        @test r.run_ok
        @test !occursin("tape34", r.run_error)
        @test !occursin("No such file", r.run_error)
    end
end
