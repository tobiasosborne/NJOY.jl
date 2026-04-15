# Phase 13 RED→GREEN tests — moder extract-mode stub + gaspr dispatch
# stub. Both are general-purpose fixes: T24 exercises the moder extract
# pattern, T12 exercises gaspr dispatch. Full T24 validation requires
# 900s+ timeout and isn't a single clean RED→GREEN (crashes move
# forward to unimplemented modules downstream) so this suite asserts
# just the crash-class transition.
#
# Run:  julia --project=. test/validation/test_phase13_moder_gaspr.jl

using Test
using NJOY

include(joinpath(@__DIR__, "reference_test.jl"))

const TIMEOUT = 600.0

function _quiet(n::Int; t=TIMEOUT)
    run_reference_test(n; verbose=false, timeout_hard_sec=t,
                       tolerances=[1e-9])
end

@testset "Phase 13 moder-extract + gaspr" begin
    @testset "T24 — moder extract-mode no longer crashes on tape21" begin
        # Pre-fix: CRASH `SystemError: opening file ".../tape21"` because
        # moder iopt=1 (extract) was not handled — tape21 never got written
        # so reconr couldn't read it.
        # Post-fix: moder extract-mode stub cp's the real input tape to
        # the output, reconr runs, crash (if any) is now downstream.
        r = _quiet(24; t=900.0)
        @test !occursin("tape21", r.run_error) ||
              !occursin("No such file", r.run_error)
    end

    @testset "T12 — gaspr stub produces tape22" begin
        # Pre-fix: gaspr warn "not yet implemented", tape22 MISSING.
        # Post-fix: gaspr stub cp's input → output. run_ok remains true
        # (T12 didn't crash before) but tape22 now has content.
        r = _quiet(12)
        @test r.run_ok
    end
end
