# T11 wimsr standalone reference test (Phase 58, 2026-05-02)
#
# Drives Julia's wimsr_module with a Fortran-generated full-precision
# ASCII GENDF (njoy-reference/tests/11 chain output, dumped via moder
# `-26 29/`) and asserts the output matches Fortran's referenceTape27
# byte-for-byte (or NUMERIC_PASS @1e-7 if FP order differs).
#
# Oracle fixtures (regenerable via the helper at the bottom):
#   gendf_input.tape   — Fortran tape29 (ASCII GENDF, 10379 lines)
#   wims_reference.tape — Fortran tape27 (= njoy-reference/tests/11/referenceTape27)
#
# Expected to FAIL until Phase 58 implementation lands (red bar before code,
# CLAUDE.md Law 1).

using Test
using NJOY
using NJOY: TapeManager, register!, resolve
using NJOY: parse_njoy_input, ModuleCall

const FIXTURE_DIR = joinpath(@__DIR__, "oracle_cache", "test11", "wimsr_standalone")
const GENDF_PATH  = joinpath(FIXTURE_DIR, "gendf_input.tape")
const WIMS_REF    = joinpath(FIXTURE_DIR, "wims_reference.tape")

if !isfile(GENDF_PATH) || !isfile(WIMS_REF)
    @info "T11 wimsr standalone: oracle fixtures missing — skipping. " *
          "Regenerate via the comment block at the end of this file " *
          "(needs Fortran NJOY binary + njoy-reference/tests/resources/t404)."
    exit(0)
end

@testset "T11 wimsr standalone (Phase 58)" begin
    work_dir = mktempdir()
    tapes = TapeManager(Dict{Int,String}(), work_dir)

    # Pre-register the GENDF input on unit 26 (matches T11 deck `-26 27/`).
    register!(tapes, 26, GENDF_PATH)

    # Synthesize a ModuleCall mirroring T11 deck lines 68-73:
    #   wimsr
    #    -26 27/
    #    1/
    #    1050 1 1050./
    #    3 7 1e10 3 10.890 221 0/
    #    1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1./
    raw_cards = [
        ["-26", "27"],
        ["1"],
        ["1050", "1", "1050."],
        ["3", "7", "1e10", "3", "10.890", "221", "0"],
        fill("1.", 13),
    ]
    mc = ModuleCall(:wimsr, raw_cards)

    # Drive Julia wimsr — this currently does NOT exist; that is the red bar.
    has_module = isdefined(NJOY, :wimsr_module) && isdefined(NJOY, :parse_wimsr)
    @test has_module
    if has_module
        params = NJOY.parse_wimsr(mc)
        NJOY.wimsr_module(tapes, params)

        out_path = resolve(tapes, 27)
        @test isfile(out_path)

        if isfile(out_path)
            # Bit-identical first; if not, classify the diff.
            ref_bytes = read(WIMS_REF)
            jul_bytes = read(out_path)
            bit_identical = ref_bytes == jul_bytes
            @test bit_identical

            if !bit_identical
                # Quick diff diagnostics — line-by-line until first divergence.
                ref_lines = readlines(WIMS_REF)
                jul_lines = readlines(out_path)
                println("\n=== T11 wimsr DIFF ===")
                println("ref lines: $(length(ref_lines)), jul lines: $(length(jul_lines))")
                ndiff = 0
                for i in 1:min(length(ref_lines), length(jul_lines))
                    if ref_lines[i] != jul_lines[i]
                        ndiff += 1
                        if ndiff <= 5
                            println("  L$i ref: |$(ref_lines[i])|")
                            println("  L$i jul: |$(jul_lines[i])|")
                        end
                    end
                end
                println("total differing lines: $ndiff (showing first 5)")
            end
        end
    end
end

# ============================================================================
# Oracle regeneration (run by hand if fixtures need refreshing)
# ============================================================================
#
#   mkdir -p /tmp/t11_oracle_v2
#   cp njoy-reference/tests/resources/t404 /tmp/t11_oracle_v2/tape20
#   # (paste njoy-reference/tests/11/input with `moder \n -26 29/` injected
#   #  immediately before the wimsr call to dump groupr's binary tape26 as
#   #  ASCII tape29 — see worklog/T11_phase58_wimsr_wiring.md)
#   cd /tmp/t11_oracle_v2 && njoy-reference/build/njoy < input
#   cmp tape27 njoy-reference/tests/11/referenceTape27   # must be 0
#   cp tape29 test/validation/oracle_cache/test11/wimsr_standalone/gendf_input.tape
#   cp tape27 test/validation/oracle_cache/test11/wimsr_standalone/wims_reference.tape
