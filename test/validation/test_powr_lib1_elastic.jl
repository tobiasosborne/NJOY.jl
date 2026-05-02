# powr lib=1 (fast/GAMTAP) MF=6 elastic matrix test — Phase 61 (2026-05-02)
#
# Drives Julia's powr_module against a Fortran-generated GENDF for carbon
# (MAT=1306, T=296K) that includes BOTH MF=3 cross sections AND MF=6/MT=2
# elastic transfer matrix (added via `6 2/` to the groupr deck). This
# exercises the gamll lap/ldp/xlol(3) detection AND the gamxs P0+P1 elastic
# matrix accumulation paths inside fast().
#
# Output is 84 × 80-char lines: 4 header + 12 absorption + 34 P0 elastic
# matrix + 34 P1 elastic matrix.
#
# Acceptance: bit-identical to Fortran NJOY's tape50.

using Test
using NJOY
using NJOY: TapeManager, register!, resolve
using NJOY: ModuleCall, parse_powr, powr_module

const POWR_E_FIXTURE_DIR = joinpath(@__DIR__, "oracle_cache", "powr_lib1_elastic")
const POWR_E_GENDF       = joinpath(POWR_E_FIXTURE_DIR, "tape23")
const POWR_E_REF_TAPE    = joinpath(POWR_E_FIXTURE_DIR, "tape50")

if !isfile(POWR_E_GENDF) || !isfile(POWR_E_REF_TAPE)
    @info "powr lib=1 elastic standalone: oracle fixtures missing — skipping. " *
          "Regenerate via the comment block at the end of this file."
    exit(0)
end

@testset "powr lib=1 elastic (carbon 68-group, +MF6 elastic matrix)" begin
    work_dir = mktempdir()
    tapes = TapeManager(Dict{Int,String}(), work_dir)
    register!(tapes, 23, POWR_E_GENDF)

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

    ref_bytes = read(POWR_E_REF_TAPE)
    jul_bytes = read(out_path)
    @test length(ref_bytes) == length(jul_bytes)
    bit_identical = ref_bytes == jul_bytes
    @test bit_identical

    if !bit_identical
        ref_lines = readlines(POWR_E_REF_TAPE)
        jul_lines = readlines(out_path)
        println("\n=== powr lib=1 elastic DIFF ===")
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

#=  ORACLE REGENERATION RECIPE
#
#   cd test/validation/oracle_cache/powr_lib1_elastic
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
#    'c-nat 68-group GAM-I + MF6 elastic'/
#    296.
#    1.e10/
#    3/
#    6 2/
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
