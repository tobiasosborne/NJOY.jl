using Test
using NJOY

# NJOY_jl-cgy: BROADR must preserve RECONR's complete ENDF-6 MF1/MT451
# metadata and canonical blank Hollerith description.
# Ref: njoy-reference/src/broadr.f90:638-680 (copy header and dictionary),
#      njoy-reference/src/reconr.f90:419-453 (blank card for ncards=0).

const BROADR_HEADER_T44_DIR = joinpath(
    @__DIR__, "..", "..", "njoy-reference", "tests", "44")
const BROADR_HEADER_T44_ORACLE = joinpath(BROADR_HEADER_T44_DIR, "referenceTape35")
const BROADR_HEADER_T44_WORK = "/tmp/njoy_cgy_t44_header"

@testset "T44 tape35 MF1 header lines 1-13 are byte-identical" begin
    rm(BROADR_HEADER_T44_WORK; recursive=true, force=true)
    run_njoy(joinpath(BROADR_HEADER_T44_DIR, "input");
             work_dir=BROADR_HEADER_T44_WORK, verbose=false)

    trial_path = joinpath(BROADR_HEADER_T44_WORK, "tape35")
    @test isfile(trial_path)

    reference = readlines(BROADR_HEADER_T44_ORACLE)
    trial = readlines(trial_path)
    @test trial[1:13] == reference[1:13]
end
