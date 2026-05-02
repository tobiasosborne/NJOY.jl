# powr lib=1 (fast/GAMTAP) full-matrix test — Phase 62 (2026-05-02)
#
# Adds MF=6 inelastic (MT=51, MT=91) on top of the elastic matrix from
# Phase 61. Exercises the gamll lap/ldp + lain/ldin detection paths and
# the gamxs accumulation for elastic + inelastic together.
#
# Output: 4 header + 12 absorption + 34 P0-elastic + 34 P1-elastic +
# 29 inelastic = 113 lines.
#
# Acceptance: bit-identical to Fortran NJOY's tape50.

using Test
using NJOY
using NJOY: TapeManager, register!, resolve
using NJOY: ModuleCall, parse_powr, powr_module

const POWR_M_FIXTURE = joinpath(@__DIR__, "oracle_cache", "powr_lib1_full_matrices")
const POWR_M_GENDF   = joinpath(POWR_M_FIXTURE, "tape23")
const POWR_M_REF     = joinpath(POWR_M_FIXTURE, "tape50")

if !isfile(POWR_M_GENDF) || !isfile(POWR_M_REF)
    @info "powr lib=1 full-matrix standalone: oracle fixtures missing — skipping."
    exit(0)
end

@testset "powr lib=1 (carbon, MF6 elastic + inelastic)" begin
    work_dir = mktempdir()
    tapes = TapeManager(Dict{Int,String}(), work_dir)
    register!(tapes, 23, POWR_M_GENDF)

    raw_cards = [
        ["23", "50"],
        ["1", "1", "1"],
        ["1306", "296.", "1", "0", "1"],
        ["'carbon'"],
        ["'thermal'"],
        ["0"],
    ]
    mc = ModuleCall(:powr, raw_cards)
    params = parse_powr(mc)

    powr_module(tapes, params)
    out_path = resolve(tapes, 50)
    @test isfile(out_path)

    ref_bytes = read(POWR_M_REF)
    jul_bytes = read(out_path)
    @test length(ref_bytes) == length(jul_bytes)
    bit_identical = ref_bytes == jul_bytes
    @test bit_identical

    if !bit_identical
        ref_lines = readlines(POWR_M_REF)
        jul_lines = readlines(out_path)
        println("\n=== powr lib=1 full-matrix DIFF ===")
        println("ref lines: $(length(ref_lines))   jul lines: $(length(jul_lines))")
        ndiff = 0
        for i in 1:min(length(ref_lines), length(jul_lines))
            if ref_lines[i] != jul_lines[i]
                ndiff += 1
                if ndiff <= 8
                    println("L$i:")
                    println("  R: \"$(ref_lines[i])\"")
                    println("  J: \"$(jul_lines[i])\"")
                end
            end
        end
        println("matching lines: $(min(length(ref_lines), length(jul_lines)) - ndiff) / $(length(ref_lines))")
    end
end
