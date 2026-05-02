# mixr standalone reference test (initial port — 2026-05-02)
#
# No reference test in njoy-reference/tests/ exercises mixr, so the oracle is
# generated locally by running Fortran NJOY's mixr against the carbon (T01,
# MAT=1306) post-reconr PENDF with weight 1.0. A single-input weight-1.0 mix
# round-trips: output MF=3 = input MF=3 (modulo MF=1/MT=451 rebuild + sigfig).
#
# Oracle fixtures (regenerated via the recipe at the bottom):
#   oracle_cache/mixr_standalone/tape21  — input PENDF (carbon after_reconr)
#   oracle_cache/mixr_standalone/tape50  — Fortran NJOY mixr reference
#
# Acceptance: bit-identical to tape50 (single-input, weight=1.0, sigfig=7
# pass-through). If anything diffs it's a real bug in the port, not FP noise.

using Test
using NJOY
using NJOY: TapeManager, register!, resolve
using NJOY: ModuleCall, parse_mixr, mixr_module

const MIXR_FIXTURE_DIR  = joinpath(@__DIR__, "oracle_cache", "mixr_standalone")
const MIXR_INPUT_TAPE   = joinpath(MIXR_FIXTURE_DIR, "tape21")
const MIXR_REF_TAPE     = joinpath(MIXR_FIXTURE_DIR, "tape50")
const MIXR_SELFMIX_TAPE = joinpath(MIXR_FIXTURE_DIR, "tape51")

if !isfile(MIXR_INPUT_TAPE) || !isfile(MIXR_REF_TAPE)
    @info "mixr standalone: oracle fixtures missing — skipping. Regenerate " *
          "via the comment block at the end of this file."
    exit(0)
end

@testset "mixr standalone (carbon, weight=1.0 round-trip)" begin
    work_dir = mktempdir()
    tapes = TapeManager(Dict{Int,String}(), work_dir)
    register!(tapes, 21, MIXR_INPUT_TAPE)

    # Mirror oracle deck:
    #   mixr
    #    50 21/
    #    1 2 102 0/
    #    1306 1.0 0/
    #    0./
    #    1306 6.000000+3 1.189690+1/
    #    'mixr standalone test for c-nat'/
    raw_cards = [
        ["50", "21"],
        ["1", "2", "102", "0"],
        ["1306", "1.0", "0"],
        ["0."],
        ["1306", "6.000000+3", "1.189690+1"],
        ["'mixr standalone test for c-nat'"],
    ]
    mc = ModuleCall(:mixr, raw_cards)

    params = parse_mixr(mc)
    @test params.nout == 50
    @test params.nin == [21]
    @test params.mtn == [1, 2, 102]
    @test params.matn == [1306]
    @test params.wtn == [1.0]
    @test params.matd == 1306
    @test params.des == "mixr standalone test for c-nat"

    mixr_module(tapes, params)
    out_path = resolve(tapes, 50)
    @test isfile(out_path)

    if isfile(out_path)
        ref_bytes = read(MIXR_REF_TAPE)
        jul_bytes = read(out_path)
        bit_identical = ref_bytes == jul_bytes
        @test bit_identical

        if !bit_identical
            ref_lines = readlines(MIXR_REF_TAPE)
            jul_lines = readlines(out_path)
            println("\n=== mixr DIFF ===")
            println("ref lines: $(length(ref_lines))   jul lines: $(length(jul_lines))")
            ndiff = 0
            for i in 1:min(length(ref_lines), length(jul_lines))
                if ref_lines[i] != jul_lines[i]
                    ndiff += 1
                    if ndiff <= 10
                        println("L$i:")
                        println("  R: \"$(ref_lines[i])\"")
                        println("  J: \"$(jul_lines[i])\"")
                    end
                end
            end
            n_match = min(length(ref_lines), length(jul_lines)) - ndiff
            println("matching lines: $n_match / $(length(ref_lines))")
        end
    end
end

# Second test: two-input self-mix with weights 0.5 + 0.5 = 1.0. Exercises the
# multi-input union-grid + weighted-sum path. 0.5*x is exact in IEEE 754 and
# 0.5*x + 0.5*x = x exactly, so the result should still be bit-identical.
if isfile(MIXR_SELFMIX_TAPE)
@testset "mixr standalone (self-mix, 0.5 + 0.5 = 1.0)" begin
    work_dir = mktempdir()
    tapes = TapeManager(Dict{Int,String}(), work_dir)
    register!(tapes, 21, MIXR_INPUT_TAPE)
    register!(tapes, 22, MIXR_INPUT_TAPE)   # same tape, different unit

    raw_cards = [
        ["51", "21", "22"],
        ["1", "2", "102", "0"],
        ["1306", "0.5", "1306", "0.5", "0"],
        ["0."],
        ["1306", "6.000000+3", "1.189690+1"],
        ["'mixr self-mix 0.5+0.5 c-nat'"],
    ]
    mc = ModuleCall(:mixr, raw_cards)
    params = parse_mixr(mc)
    @test params.nin == [21, 22]
    @test params.matn == [1306, 1306]
    @test params.wtn == [0.5, 0.5]

    mixr_module(tapes, params)
    out_path = resolve(tapes, 51)
    @test isfile(out_path)
    bit_identical = read(MIXR_SELFMIX_TAPE) == read(out_path)
    @test bit_identical
end
end

#=  ORACLE REGENERATION RECIPE
#
#   cd test/validation/oracle_cache/mixr_standalone
#   cp ../test01/after_reconr.pendf tape21
#   cat > input <<'EOF'
#   mixr
#    50 21/
#    1 2 102 0/
#    1306 1.0 0/
#    0./
#    1306 6.000000+3 1.189690+1/
#    'mixr standalone test for c-nat'/
#   stop
#   EOF
#   ../../../../njoy-reference/build/njoy < input > output
#
#   tape50 is the reference. tape21 is the input. Commit both.
=#
