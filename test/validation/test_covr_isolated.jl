# test_covr_isolated.jl -- Bit-identity check for covr against the
# Fortran reference errorr tapes.
#
# Why isolated? Julia's errorr (an upstream module) emits a sparse subset
# of the full Fortran reference covariance tape (see HANDOFF P1 —
# "covcal content drift"). Running the full T05/T16 pipeline through
# Julia's errorr → covr therefore diffs against ref because the *errorr*
# output diverges. To test covr in isolation we feed the Fortran
# reference errorr tape (`referenceTapeNN`) directly to `covr_module` and
# compare against `referenceTape{plot_unit}`.
#
# When errorr lands its full multigroup covariance (Bug A and Bug B
# under HANDOFF "Covcal content drift"), the full-pipeline T05/T16
# reference tests will pass too — without further changes to covr.

using NJOY
using Test

const REPO = joinpath(@__DIR__, "..", "..")

function _run_covr_iso(work_dir::AbstractString,
                       in_unit::Int, in_tape_name::AbstractString,
                       params::NJOY.CovrParams)
    mkpath(work_dir)
    cp(joinpath(REPO, "njoy-reference", "tests", in_tape_name),
       joinpath(work_dir, "tape$in_unit"); force=true)
    tm = NJOY.TapeManager(Dict{Int,String}(), work_dir)
    NJOY.register!(tm, in_unit, joinpath(work_dir, "tape$in_unit"))
    NJOY.covr_module(tm, params)
    tm
end

function _check_bit_identical(label::AbstractString,
                              julia_path::AbstractString,
                              ref_path::AbstractString)
    @assert isfile(julia_path) "$label: Julia output missing at $julia_path"
    @assert isfile(ref_path)  "$label: reference missing at $ref_path"
    j = read(julia_path)
    r = read(ref_path)
    if j == r
        @info "$label: BIT-IDENTICAL ($(length(j)) bytes)"
        return true
    else
        @warn "$label: differs — Julia=$(length(j))B vs Ref=$(length(r))B"
        return false
    end
end

@testset "covr isolation — Fortran reference errorr tape → bit-identical plot tape" begin
    # ---- T05: C-12 (MAT=1306) MF33, ncase=1 imt=0 expansion ----
    @testset "T05 referenceTape34 (C-12 MF33 plot)" begin
        wd = mktempdir(; cleanup=true)
        params = NJOY.CovrParams(
            33, 0, 34, 1,
            [0.001, 0.1, 0.2, 0.3, 0.6, 1.0],
            0.0, 1, 1, 0, 1, 1, 3, "", "",
            [NJOY.CovrCase(1306, 0, 0, 0)])
        tm = _run_covr_iso(wd, 33, "05/referenceTape33", params)
        @test _check_bit_identical(
            "T05 tape34", joinpath(wd, "tape34"),
            joinpath(REPO, "njoy-reference", "tests", "05", "referenceTape34"))
    end

    # ---- T16: U-238 (MAT=9237) MF33 with explicit pair list ----
    @testset "T16 referenceTape36 (U-238 MF33 plot, 14 explicit pairs)" begin
        wd = mktempdir(; cleanup=true)
        # 14 (mat,mt,mat1,mt1) pairs from njoy-reference/tests/16/input.
        cases = NJOY.CovrCase[
            NJOY.CovrCase(9237,   1, 9237,   1),
            NJOY.CovrCase(9237,   1, 9237,   2),
            NJOY.CovrCase(9237,   2, 9237,   2),
            NJOY.CovrCase(9237,   2, 9237,   4),
            NJOY.CovrCase(9237,   2, 9237,  16),
            NJOY.CovrCase(9237,   2, 9237,  17),
            NJOY.CovrCase(9237,   2, 9237,  18),
            NJOY.CovrCase(9237,   2, 9237, 102),
            NJOY.CovrCase(9237,   4, 9237,   4),
            NJOY.CovrCase(9237,  16, 9237,  16),
            NJOY.CovrCase(9237,  17, 9237,  17),
            NJOY.CovrCase(9237,  18, 9237,  18),
            NJOY.CovrCase(9237,  18, 9237, 102),
            NJOY.CovrCase(9237, 102, 9237, 102),
        ]
        params = NJOY.CovrParams(
            26, 0, 36, 1,
            [0.001, 0.1, 0.2, 0.3, 0.6, 1.0],
            0.0, 1, 14, 0, 1, 1, 3, "", "", cases)
        tm = _run_covr_iso(wd, 26, "16/referenceTape26", params)
        @test _check_bit_identical(
            "T16 tape36", joinpath(wd, "tape36"),
            joinpath(REPO, "njoy-reference", "tests", "16", "referenceTape36"))
    end

    # ---- T16 second covr card: U-238 MF34 (mubar) ----
    @testset "T16 referenceTape37 (U-238 MF34 mubar plot)" begin
        wd = mktempdir(; cleanup=true)
        params = NJOY.CovrParams(
            27, 0, 37, 1,
            [0.001, 0.1, 0.2, 0.3, 0.6, 1.0],
            0.0, 1, 1, 0, 1, 1, 3, "", "",
            [NJOY.CovrCase(9237, 0, 0, 0)])
        tm = _run_covr_iso(wd, 27, "16/referenceTape27", params)
        @test _check_bit_identical(
            "T16 tape37", joinpath(wd, "tape37"),
            joinpath(REPO, "njoy-reference", "tests", "16", "referenceTape37"))
    end
end
