# Fast isolated oracle for build_thermal_grid on T67 (D-in-7LiD, LTHR=3).
#
# Reproduces the NJOY_jl-km2 bug at the build_thermal_grid level WITHOUT running
# the ~timeout full T67 thermr pipeline. The reference merged MT221 grid lives on
# njoy-reference/tests/67/referenceTape60 (7060 points). The Fortran `coh`
# subroutine (thermr.f90:748-922) builds it by an adaptive LIFO-stack merge of
# the broadened elastic grid (MT2) with the Bragg edges of MF7/MT2 (3374 edges),
# terminating not on edge-count but on the emax cutoff (thermr.f90:870-875) after
# `sigcoh` returns the emax sentinel (thermr.f90:1192-1195) once the real edges
# are exhausted.
#
# Inputs (T67 card-2 `15 128 …`, card-4 `0.001 10.`):
#   tsl  : tsl-DinLiD-ENDF8.1-beta.endf, MAT=15, T=400  -> lthr=3, n_edges=3374
#   pendf: referenceTape60, MAT=128  -> MT2 = broadened elastic grid (input),
#                                       MT221 = reference merged grid (7060 pts)
#   emax = 10.0, tol = 0.001
#
# Run:  julia --project=. test/validation/test_thermr_build_grid_t67.jl

using Test
using NJOY

const REFDIR = joinpath(@__DIR__, "..", "..", "njoy-reference", "tests")
const RES    = joinpath(REFDIR, "resources")

@testset "build_thermal_grid T67 oracle (NJOY_jl-km2)" begin
    tsl_path = joinpath(RES, "tsl-DinLiD-ENDF8.1-beta.endf")
    lthr, bragg, incoh = NJOY.read_mf7_mt2(tsl_path, 15, 400.0)
    @test lthr == 3
    @test bragg !== nothing
    @test bragg.n_edges == 3374

    # searchsorted prefix-sum precondition (Defect-2 perf fix relies on this).
    @test issorted(bragg.tau_sq)

    tape = NJOY.read_pendf(joinpath(REFDIR, "67", "referenceTape60"))
    mf3 = NJOY.extract_mf3_all(tape, 128)
    @test haskey(mf3, 2)    # broadened elastic input grid
    @test haskey(mf3, 221)  # reference merged thermal grid

    mt2_energies = mf3[2][1]
    ref_grid     = mf3[221][1]
    @test length(ref_grid) == 7060

    t0 = time()
    grid = NJOY.build_thermal_grid(bragg, mt2_energies, 10.0; tol=0.001)
    NJOY._append_emax_sentinels!(grid, 10.0)
    dt = time() - t0
    @info "build_thermal_grid T67: $(length(grid)) pts (ref $(length(ref_grid))) in $(round(dt; digits=2))s"

    # PRIMARY: structural oracle. Length is exact; values match the written
    # reference to the `round_sigfig` upward bias (×1.0000000000001, CLAUDE.md
    # trap #2). The reference MT221 grid is parsed from the a11-formatted PENDF,
    # so the bias has been stripped by format_endf_float on write; build_thermal_grid
    # returns the internal biased Float64. isapprox at the bias tolerance is the
    # faithful structural comparison (the production write path re-applies the same
    # sigfig formatting, so the tape is bit-identical).
    @test length(grid) == 7060
    @test all(isapprox.(grid, ref_grid; rtol=2e-13))

    # Regression guards: the input points the old edge-count loop dropped, and the
    # prefix that always agreed (the divergence used to start at index 2314).
    @test any(x -> isapprox(x, 1.5; rtol=2e-13), grid)
    @test any(x -> isapprox(x, 5.0; rtol=2e-13), grid)
    n_common = min(2312, length(grid), length(ref_grid))
    @test all(isapprox.(grid[1:n_common], ref_grid[1:n_common]; rtol=2e-13))

    # Tail structure: the [last_real_edge, emax] interval is now processed —
    # the binary-subdivision sequence 5.0/5.3125/.../9.6875 and the emax sentinels.
    @test all(any(x -> isapprox(x, v; rtol=2e-13), grid)
              for v in (5.3125, 5.625, 5.9375, 6.25, 9.687493))
end
