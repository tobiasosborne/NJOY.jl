using Test
using NJOY

# NJOY_jl-4lw: exact LEAPR endout serialization boundaries.
# Ref: njoy-reference/src/leapr.f90:3354-3417 (sigfig before TAB1 output),
#      njoy-reference/src/endf.f90:882-981 (canonical a11 extended form).

const LEAPR_T33_DIR = joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "33")
const LEAPR_T33_WORK = "/tmp/njoy_4lw_t33_serialization"
const LEAPR_T33_LINES = Dict(
    24 => [5886, 9668, 15683, 53147],
    34 => [2, 76, 78, 2850, 4388, 6469, 13584, 53147],
)

@testset "T33 LEAPR residual tape24/tape34 records are byte-identical" begin
    rm(LEAPR_T33_WORK; recursive=true, force=true)
    run_njoy(joinpath(LEAPR_T33_DIR, "input");
             work_dir=LEAPR_T33_WORK, verbose=false)

    for tape in (24, 34)
        trial_path = joinpath(LEAPR_T33_WORK, "tape$tape")
        @test isfile(trial_path)

        reference = readlines(joinpath(LEAPR_T33_DIR, "referenceTape$tape"))
        trial = readlines(trial_path)
        lines = LEAPR_T33_LINES[tape]
        @test trial[lines] == reference[lines]
    end
end
