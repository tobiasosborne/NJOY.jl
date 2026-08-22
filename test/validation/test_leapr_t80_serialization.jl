using Test
using NJOY

# NJOY_jl-sc5: exact MF7/MT4 B(4) serialization for T80.
# Ref: njoy-reference/src/leapr.f90:3291-3323 (endout), especially
#      scr(10)=sigfig(therm*beta(nbeta),7,0) before listio.
# ENDF-6 Formats Manual, Section 7.4: B(4) is EMAX for the principal scatterer.

const LEAPR_T80_DIR = joinpath(
    @__DIR__, "..", "..", "njoy-reference", "tests", "80")
const LEAPR_T80_WORK = "/tmp/njoy_sc5_t80_serialization"

@testset "T80 tape24 line 19 MF7/MT4 B(4) is byte-identical" begin
    rm(LEAPR_T80_WORK; recursive=true, force=true)
    run_njoy(joinpath(LEAPR_T80_DIR, "input");
             work_dir=LEAPR_T80_WORK, verbose=false)

    trial_path = joinpath(LEAPR_T80_WORK, "tape24")
    @test isfile(trial_path)

    reference = readlines(joinpath(LEAPR_T80_DIR, "referenceTape24"))
    trial = readlines(trial_path)
    @test trial[19] == reference[19]
end
