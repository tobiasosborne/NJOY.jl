# powr lib=1 (fast/GAMTAP) standalone test — Phase B initial port (2026-05-02)
#
# Drives Julia's powr_module against a Fortran-generated GENDF for carbon
# (MAT=1306, T=296K, 68-group GAM-I structure ign=6, single sigma-zero).
# The carbon path exercises the absorption-only branch (iwa=1, iwf=0,
# iwr=0, no MF=6 matrices), which is the simplest of fast()'s code paths.
# Bit-identical to Fortran NJOY's tape50.
#
# Oracle fixtures (regen recipe at the bottom):
#   oracle_cache/powr_lib1_standalone/tape23 — input GENDF
#   oracle_cache/powr_lib1_standalone/tape50 — Fortran NJOY GAMTAP reference
#
# Phase B follow-ups (separate sessions):
#   - Multi-temperature accumulation (ntp > 1).
#   - MF=6 matrix paths (P0/P1 elastic, inelastic, n2n).
#   - Fission paths (locsf0/locchi/locnus + MT=18..21+38).
#   - gamff self-shielding factors (iff=1 with nsigz > 1).
#   - matd<0 iread=1 absorption-direct read.

using Test
using NJOY
using NJOY: TapeManager, register!, resolve
using NJOY: ModuleCall, parse_powr, powr_module

const POWR_FIXTURE_DIR = joinpath(@__DIR__, "oracle_cache", "powr_lib1_standalone")
const POWR_GENDF       = joinpath(POWR_FIXTURE_DIR, "tape23")
const POWR_REF_TAPE    = joinpath(POWR_FIXTURE_DIR, "tape50")

if !isfile(POWR_GENDF) || !isfile(POWR_REF_TAPE)
    @info "powr lib=1 standalone: oracle fixtures missing — skipping. " *
          "Regenerate via the comment block at the end of this file."
    exit(0)
end

@testset "powr lib=1 standalone (carbon 68-group, absorption-only)" begin
    work_dir = mktempdir()
    tapes = TapeManager(Dict{Int,String}(), work_dir)
    register!(tapes, 23, POWR_GENDF)

    # Mirror oracle deck (the powr block of test/validation/oracle_cache/
    # powr_lib1_standalone/input):
    #   powr
    #    23 50/
    #    1 1 1/
    #    1306 296. 1 0 1/
    #    'carbon'/
    #    'thermal'/
    #    0/
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

    @test params.lib == 1
    @test params.ngendf == 23
    @test params.nout == 50
    @test length(params.fast_mats) == 1
    @test params.fast_mats[1].matd == 1306
    @test params.fast_mats[1].rtemp == 296.0
    @test params.fast_mats[1].iff == 1
    @test params.fast_mats[1].word == "carbon"
    @test params.fast_mats[1].fsn == "thermal"

    powr_module(tapes, params)
    out_path = resolve(tapes, 50)
    @test isfile(out_path)

    if isfile(out_path)
        ref_bytes = read(POWR_REF_TAPE)
        jul_bytes = read(out_path)
        @test length(ref_bytes) == length(jul_bytes)
        bit_identical = ref_bytes == jul_bytes
        @test bit_identical

        if !bit_identical
            ref_lines = readlines(POWR_REF_TAPE)
            jul_lines = readlines(out_path)
            println("\n=== powr lib=1 DIFF ===")
            println("ref lines: $(length(ref_lines))   jul lines: $(length(jul_lines))")
            for i in 1:min(length(ref_lines), length(jul_lines))
                if ref_lines[i] != jul_lines[i]
                    println("L$i:")
                    println("  R: \"$(ref_lines[i])\"")
                    println("  J: \"$(jul_lines[i])\"")
                end
            end
        end
    end
end

#=  ORACLE REGENERATION RECIPE
#
#   cd test/validation/oracle_cache/powr_lib1_standalone
#   cp ../../../../njoy-reference/tests/resources/t511 tape20
#   cat > input <<'EOF'
#   reconr
#    20 21/
#    'pendf for c-nat'/
#    1306 0/
#    .005/
#    0/
#   broadr
#    20 21 22/
#    1306 1/
#    .005/
#    296./
#    0/
#   groupr
#    20 22 0 23/
#    1306 6 0 3 1 1 1 1/
#    'c-nat 68-group GAM-I'/
#    296.
#    1.e10/
#    3/
#    0/
#    0/
#   powr
#    23 50/
#    1 1 1/
#    1306 296. 1 0 1/
#    'carbon'/
#    'thermal'/
#    0/
#   stop
#   EOF
#   ../../../../njoy-reference/build/njoy < input > output
=#
