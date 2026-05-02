# resxsr standalone reference test (initial port — 2026-05-02)
#
# No reference test in njoy-reference/tests/ exercises resxsr, so the oracle
# is generated locally by running Fortran NJOY's resxsr against the broadened
# carbon (T01, MAT=1306, T=296K) PENDF. The output is a CCCC-IV binary RESXS
# file containing per-temp resonance cross sections for elastic + capture
# (no fission for C-nat).
#
# Oracle fixtures (regenerable via the recipe at the bottom):
#   oracle_cache/resxsr_standalone/tape21  — input PENDF (after_broadr.pendf)
#   oracle_cache/resxsr_standalone/tape50  — Fortran NJOY resxsr reference
#
# Acceptance: bit-identical to tape50.

using Test
using NJOY
using NJOY: TapeManager, register!, resolve
using NJOY: ModuleCall, parse_resxsr, resxsr_module

const RESXSR_FIXTURE_DIR = joinpath(@__DIR__, "oracle_cache", "resxsr_standalone")
const RESXSR_INPUT_TAPE  = joinpath(RESXSR_FIXTURE_DIR, "tape21")
const RESXSR_REF_TAPE    = joinpath(RESXSR_FIXTURE_DIR, "tape50")

if !isfile(RESXSR_INPUT_TAPE) || !isfile(RESXSR_REF_TAPE)
    @info "resxsr standalone: oracle fixtures missing — skipping. Regenerate " *
          "via the comment block at the end of this file."
    exit(0)
end

@testset "resxsr standalone (carbon @296K, eps=0.001)" begin
    work_dir = mktempdir()
    tapes = TapeManager(Dict{Int,String}(), work_dir)
    register!(tapes, 21, RESXSR_INPUT_TAPE)

    # Mirror the oracle deck:
    #   resxsr
    #    50/
    #    1 1 1 1.e-5 2.e7 0.001/
    #    'NJOY.jl' 1/
    #    'resxsr standalone test for c-nat'/
    #    'cnat' 1306 21/
    raw_cards = [
        ["50"],
        ["1", "1", "1", "1.e-5", "2.e7", "0.001"],
        ["'NJOY.jl'", "1"],
        ["'resxsr standalone test for c-nat'"],
        ["'cnat'", "1306", "21"],
    ]
    mc = ModuleCall(:resxsr, raw_cards)
    params = parse_resxsr(mc)

    @test params.nout == 50
    @test params.nmat == 1
    @test params.maxt == 1
    @test params.nholl == 1
    @test params.efirst == 1.e-5
    @test params.elast  == 2.e7
    @test params.eps    == 0.001
    @test params.huse   == "NJOY.jl"
    @test params.ivers  == 1
    @test length(params.materials) == 1
    @test params.materials[1].hmat == "cnat"
    @test params.materials[1].mat  == 1306
    @test params.materials[1].unit == 21

    resxsr_module(tapes, params)
    out_path = resolve(tapes, 50)
    @test isfile(out_path)

    if isfile(out_path)
        ref_bytes = read(RESXSR_REF_TAPE)
        jul_bytes = read(out_path)
        @test length(ref_bytes) == length(jul_bytes)
        bit_identical = ref_bytes == jul_bytes
        @test bit_identical

        if !bit_identical
            n = min(length(ref_bytes), length(jul_bytes))
            first_diff = findfirst(i -> ref_bytes[i] != jul_bytes[i], 1:n)
            println("\n=== resxsr DIFF ===")
            println("ref bytes: $(length(ref_bytes))   jul bytes: $(length(jul_bytes))")
            if first_diff !== nothing
                lo = max(first_diff - 16, 1); hi = min(first_diff + 47, n)
                println("first byte diff at offset $first_diff (1-based)")
                println("  R[$lo..$hi]: ", bytes2hex(ref_bytes[lo:hi]))
                println("  J[$lo..$hi]: ", bytes2hex(jul_bytes[lo:hi]))
            end
        end
    end
end

#=  ORACLE REGENERATION RECIPE
#
#   cd test/validation/oracle_cache/resxsr_standalone
#   cp ../test01/after_broadr.pendf tape21
#   cat > input <<'EOF'
#   resxsr
#    50/
#    1 1 1 1.e-5 2.e7 0.001/
#    'NJOY.jl' 1/
#    'resxsr standalone test for c-nat'/
#    'cnat' 1306 21/
#   stop
#   EOF
#   ../../../../njoy-reference/build/njoy < input > output
=#
