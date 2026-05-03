"""
powr lib=1 delayed-neutron spectra (MF=5/MT=455) standalone test
— Phase 66 (2026-05-03).

Drives Julia's powr_module against a Fortran-generated GENDF for U-235
(MAT=9228, ENDF/B-VIII, T=296K, single sigma-zero, 68-group GAM-I) with
explicit MF=5/MT=455 inclusion in the groupr deck (groupr's auto-finder
excludes MF=5 — `5 455/` must be requested explicitly).

The Fortran reference tape50 has 131 lines:
  L1-13:   prompt chi block (i2=0 header + 12 data lines)
  L14-91:  6 delayed-chi blocks (each: i2=k header + 12 data lines, k=1..6)
  L92-131: nscr block (header iwf=1, abs XS, sigf, nu, MF=6 chi/nu carry)

This exercises the locdla pointer (between locchi and locnus, sized
ngnd*nfs), the MT=455 dispatch in gamxs at powr.f90:1135 (both mfh=3 and
mfh=5 branches), and the delayed-chi normalization at powr.f90:1371-1393
(rnorm = 1 / nu_total at jgdnu, with the sumd[il]*rnorm replacement at
ig=ngnd-1 and the locb-replacement at ig=ngnd).
"""

using Test
using NJOY
using NJOY: TapeManager, register!, resolve
using NJOY: ModuleCall, parse_powr, powr_module

const POWR_FIXTURE_DIR = joinpath(@__DIR__, "oracle_cache", "powr_lib1_delayed")
const POWR_GENDF       = joinpath(POWR_FIXTURE_DIR, "tape23")
const POWR_REF_TAPE    = joinpath(POWR_FIXTURE_DIR, "tape50")

if !isfile(POWR_GENDF) || !isfile(POWR_REF_TAPE)
    @info "powr lib=1 delayed: oracle fixtures missing — skipping. " *
          "Regenerate via the comment block at the end of this file."
    exit(0)
end

@testset "powr lib=1 delayed neutrons (U-235 68-group, nfs=6)" begin
    work_dir = mktempdir()
    tapes = TapeManager(Dict{Int,String}(), work_dir)
    register!(tapes, 23, POWR_GENDF)

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

    @test params.lib == 1
    @test params.fast_mats[1].matd == 9228
    @test params.fast_mats[1].rtemp == 296.0

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
            println("\n=== powr lib=1 delayed-neutron DIFF ===")
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
#   cd test/validation/oracle_cache/powr_lib1_delayed
#   cp ../../../../njoy-reference/tests/resources/n-092_U_235-ENDF8.0.endf tape20
#   cat > input <<'EOF'
#   reconr
#    20 21/
#    'pendf for u-235 (delayed neutrons)'/
#    9228 0/
#    .005/
#    0/
#   broadr
#    20 21 22/
#    9228 1/
#    .005/
#    296./
#    0/
#   groupr
#    20 22 0 23/
#    9228 6 0 3 1 1 1 1/
#    'u-235 68-group GAM-I + fission + delayed-chi'/
#    296.
#    1.e10/
#    3/
#    5 455 'delayed neutron spectra'/
#    6 18/
#    0/
#    0/
#   powr
#    23 50/
#    1 1 1/
#    9228 296. 1 0 1/
#    'u-235'/
#    'thermal-fission'/
#    0/
#   stop
#   EOF
#   ../../../../njoy-reference/build/njoy < input > output
=#
