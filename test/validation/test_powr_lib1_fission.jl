# powr lib=1 (fast/GAMTAP) fission path test — Phase 63 (2026-05-02)
#
# Drives Julia's powr_module against a Fortran-generated GENDF for U-235
# (MAT=9228, T=296K, 68-group GAM-I) WITH MF=6/MT=18 fission spectrum +
# nu*sigma_f. Exercises the gamxs MF=3 fission accumulation, the MF=6/MT=18
# spectrum-only branch (cspc setup at ig=0; spectrum-distribution at ig>0
# with the bit-faithful jg2c-1 quirk), the post-loop nu = nus/sigf
# normalization, and the chi-norm pass.
#
# Output structure: 1 chi header + 12 chi values + 4 main header + 12 abs +
# 12 sigf + 12 nu = 53 lines.
#
# Acceptance: bit-identical to Fortran NJOY's tape50.

using Test
using NJOY
using NJOY: TapeManager, register!, resolve
using NJOY: ModuleCall, parse_powr, powr_module

const POWR_F_FIXTURE = joinpath(@__DIR__, "oracle_cache", "powr_lib1_fission")
const POWR_F_GENDF   = joinpath(POWR_F_FIXTURE, "tape23")
const POWR_F_REF     = joinpath(POWR_F_FIXTURE, "tape50")

if !isfile(POWR_F_GENDF) || !isfile(POWR_F_REF)
    @info "powr lib=1 fission standalone: oracle fixtures missing — skipping."
    exit(0)
end

@testset "powr lib=1 (U-235, MF6 fission spectrum + nu)" begin
    work_dir = mktempdir()
    tapes = TapeManager(Dict{Int,String}(), work_dir)
    register!(tapes, 23, POWR_F_GENDF)

    raw_cards = [
        ["23", "50"],
        ["1", "1", "1"],
        ["9228", "296.", "1", "0", "1"],
        ["'u-235'"],
        ["'thermal-fission'"],
        ["0"],
    ]
    mc = ModuleCall(:powr, raw_cards)
    params = parse_powr(mc)

    powr_module(tapes, params)
    out_path = resolve(tapes, 50)
    @test isfile(out_path)

    ref_bytes = read(POWR_F_REF)
    jul_bytes = read(out_path)
    @test length(ref_bytes) == length(jul_bytes)
    bit_identical = ref_bytes == jul_bytes
    @test bit_identical

    if !bit_identical
        ref_lines = readlines(POWR_F_REF)
        jul_lines = readlines(out_path)
        println("\n=== powr lib=1 fission DIFF ===")
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
