"""
powr lib=1 multi-temperature + self-shielding (gamff) standalone test
— Phase 64+65 (2026-05-03).

Drives Julia's powr_module against a Fortran-generated multi-temperature
multi-σ₀ GENDF for Fe-56 (MAT=2631, JEFF-3.3, T = {296, 600, 1200} K,
σ₀ = {1e10, 1e3, 1e2, 1e1} barns, 68-group GAM-I structure ign=6).

This is the canonical exercise for the gamff f-factor block emitted by
powr lib=1 at powr.f90:614-641 (gated on `nff!=0 && iwr!=0`). Without
self-shielding the multi-temp paths would be unobservable; the Phase 60
carbon oracle deliberately lives in the iwr=0 corner that bypasses this
block. Fe-56 has resolved resonances and is groupr's σ₀-weighted-flux
treatment exposes per-σ₀ XS differences in the resonance range.

Oracle fixtures (regen recipe at the bottom):
  oracle_cache/powr_lib1_selfshielding/tape24 — input GENDF (multi-T, multi-σ₀)
  oracle_cache/powr_lib1_selfshielding/tape50 — Fortran NJOY GAMTAP reference
                                                (177 lines: 4 header + 12
                                                abs XS + 1 (4i6) + sigz +
                                                tmpr + 86 f-factors +
                                                72 dummy amisc zeros)
"""

using Test
using NJOY
using NJOY: TapeManager, register!, resolve
using NJOY: ModuleCall, parse_powr, powr_module

const POWR_FIXTURE_DIR = joinpath(@__DIR__, "oracle_cache", "powr_lib1_selfshielding")
const POWR_GENDF       = joinpath(POWR_FIXTURE_DIR, "tape24")
const POWR_REF_TAPE    = joinpath(POWR_FIXTURE_DIR, "tape50")

if !isfile(POWR_GENDF) || !isfile(POWR_REF_TAPE)
    @info "powr lib=1 self-shielding: oracle fixtures missing — skipping. " *
          "Regenerate via the comment block at the end of this file."
    exit(0)
end

@testset "powr lib=1 self-shielding (Fe-56 68-group, ntemp=3, nsigz=4)" begin
    work_dir = mktempdir()
    tapes = TapeManager(Dict{Int,String}(), work_dir)
    register!(tapes, 24, POWR_GENDF)

    raw_cards = [
        ["24", "50"],
        ["1", "1", "1"],
        ["2631", "296.", "1", "4", "1"],
        ["'iron-56'"],
        ["'resonant'"],
        ["0"],
    ]
    mc = ModuleCall(:powr, raw_cards)
    params = parse_powr(mc)

    @test params.lib == 1
    @test params.ngendf == 24
    @test params.nout == 50
    @test length(params.fast_mats) == 1
    @test params.fast_mats[1].matd == 2631
    @test params.fast_mats[1].rtemp == 296.0
    @test params.fast_mats[1].iff == 1
    @test params.fast_mats[1].nsgz == 4
    @test params.fast_mats[1].word == "iron-56"

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
            println("\n=== powr lib=1 self-shielding DIFF ===")
            println("ref lines: $(length(ref_lines))   jul lines: $(length(jul_lines))")
            n_diffs = 0
            for i in 1:min(length(ref_lines), length(jul_lines))
                if ref_lines[i] != jul_lines[i]
                    n_diffs += 1
                    if n_diffs <= 30
                        println("L$i:")
                        println("  R: \"$(ref_lines[i])\"")
                        println("  J: \"$(jul_lines[i])\"")
                    end
                end
            end
            n_diffs > 30 && println("... ($(n_diffs - 30) more diffs)")
            println("Total differing lines: $n_diffs")
        end
    end
end

#=  ORACLE REGENERATION RECIPE
#
#   cd test/validation/oracle_cache/powr_lib1_selfshielding
#   cp ../../../../njoy-reference/tests/resources/n-026_Fe_056-JEFF3.3.endf tape20
#   cat > input <<'EOF'
#   reconr
#    20 21/
#    'pendf for fe-56 (jeff3.3)'/
#    2631 0/
#    .005/
#    0/
#   broadr
#    20 21 22/
#    2631 3/
#    .005/
#    296. 600. 1200./
#    0/
#   unresr
#    20 22 23/
#    2631 3 4 1/
#    296. 600. 1200./
#    1.e10 1.e3 1.e2 1.e1/
#    0/
#   groupr
#    20 23 0 24/
#    2631 6 0 3 1 3 4 1/
#    'fe-56 68-group GAM-I (multi-T multi-sigma_z)'/
#    296. 600. 1200./
#    1.e10 1.e3 1.e2 1.e1/
#    3/
#    0/
#    3/
#    0/
#    3/
#    0/
#    0/
#   powr
#    24 50/
#    1 1 1/
#    2631 296. 1 4 1/
#    'iron-56'/
#    'resonant'/
#    0/
#   stop
#   EOF
#   ../../../../njoy-reference/build/njoy < input > output
#
#   Note: Fe-56 has no unresolved-resonance parameters. unresr will report
#   "copy as is to nout" — that is expected. The σ₀-weighted flux that
#   produces the per-σ₀ XS differences in the resolved range is computed by
#   groupr's narrow-resonance treatment.
=#
