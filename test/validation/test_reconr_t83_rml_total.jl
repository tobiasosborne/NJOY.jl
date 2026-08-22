using Test
using NJOY

# NJOY_jl-fod: a zero MF3 smooth background must not suppress an RML reaction
# contribution before redundant MT1 accumulation.  Fortran initializes sn=0,
# then unconditionally adds res(1+itype) before sigfig and total accumulation.
# Ref: njoy-reference/src/reconr.f90:4791-4809,4832-4893 (emerge).

const RECONR_T83_DIR = joinpath(
    @__DIR__, "..", "..", "njoy-reference", "tests", "83")
const RECONR_T83_ORACLE = joinpath(RECONR_T83_DIR, "referenceTape50")
const RECONR_T83_WORK = "/tmp/njoy_fod_t83_rml_total"

@testset "T83 MT152 and first resolved MT1 record are byte-identical" begin
    rm(RECONR_T83_WORK; recursive=true, force=true)
    run_njoy(joinpath(RECONR_T83_DIR, "input");
             work_dir=RECONR_T83_WORK, verbose=false)

    trial_path = joinpath(RECONR_T83_WORK, "tape50")
    @test isfile(trial_path)

    reference = readlines(RECONR_T83_ORACLE)
    trial = readlines(trial_path)
    @test trial[101:130] == reference[101:130]
    @test trial[134] == reference[134]
end
