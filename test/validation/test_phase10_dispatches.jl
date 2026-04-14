# Phase 10 RED→GREEN tests for items 1-5 from worklog/T09_acer_dispatch.md
#
# Each test asserts `run_ok` (no crash). Pre-fix: all RED (specific crash
# signatures from reports/REFERENCE_SWEEP.md). Post-fix: all GREEN.
#
# Run:  julia --project=. test/validation/test_phase10_dispatches.jl

using Test
using NJOY

# Pull in the reference-test runner (NJOY-internal helper file)
include(joinpath(@__DIR__, "reference_test.jl"))

const TIMEOUT = 600.0   # most affected tests are < 60s; URR-heavy ones up to ~10min

"""Run a reference test and return its NamedTuple result."""
function _quiet(n::Int)
    run_reference_test(n; verbose=false, timeout_hard_sec=TIMEOUT,
                       tolerances=[1e-9])
end

@testset "Phase 10 RED→GREEN" begin

    @testset "Item 3 — heatr h scope bug (T08, T26, T49, T79)" begin
        # Pre-fix crash: UndefVarError: `h` not defined in local scope
        # All four use heatr; T08 is fastest.
        for n in (8, 26, 49, 79)
            r = _quiet(n)
            @test r.run_ok      # no crash
            @test !occursin("UndefVarError", r.run_error)
        end
    end

    @testset "Item 4 — Bragg lookup gating + entries (T25, T32, T67-70, T74)" begin
        # T32, T68 have icoh=0 → Bragg should be skipped (gating fix)
        # T25, T67, T69, T70, T74 have icoh>0 → need real Bragg entries
        for n in (25, 32, 67, 68, 69, 70, 74)
            r = _quiet(n)
            @test r.run_ok
            @test !occursin("No Bragg lattice parameters", r.run_error)
        end
    end

    @testset "Item 5 — INT=0 in MF TAB1 reader (T15, T17)" begin
        # JENDL-3.3 U-238 emits INT=0; Fortran NJOY silently treats as linear
        for n in (15, 17)
            r = _quiet(n)
            @test r.run_ok
            @test !occursin("InterpolationLaw: 0", r.run_error)
        end
    end

    @testset "Item 1 — purr dispatch (T35, T36, T63, T28)" begin
        # Pre-fix crash: SystemError opening tape43 (purr never wrote it)
        # Stub: copy npendf_in → npendf_out so downstream finds something
        for n in (35, 36, 63, 28)
            r = _quiet(n)
            @test r.run_ok
        end
    end

    @testset "Item 2 — leapr dispatch (T22, T23, T80)" begin
        # Pre-fix: T22/T23/T80 ran (no leapr call) but produced no tape20.
        # Post-fix: leapr_module writes a (possibly empty) tape20 →
        #          test framework finds a file (DIFFS or STRUCTURAL_FAIL),
        #          run_ok stays true. Asserts only no crash.
        for n in (22, 23, 80)
            r = _quiet(n)
            @test r.run_ok
        end
    end

end
