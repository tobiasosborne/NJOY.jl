using Test
using NJOY

# Direct-RECONR PENDF serialization regression for NJOY_jl-gkw.
#
# T84 is the smallest complete oracle for this path.  Its referenceTape100
# lines 1-1115 exercise the TPID, the four-record ENDF-6 MF1/MT451 header, the
# canonical blank Hollerith card inserted for `ncards=0`, the full dictionary,
# continuous per-file sequence numbers, MF2/MF3, and reconstructed MF12.
#
# Ref: njoy-reference/src/reconr.f90:419-453 (ruina inserts the blank card),
#      :5028-5174 (recout writes MF1 and its dictionary),
#      :5178-5441 (recout writes MF2 and the remaining files).

const RECONR_SERIAL_T84_DIR = joinpath(
    @__DIR__, "..", "..", "njoy-reference", "tests", "84")
const RECONR_SERIAL_T84_ORACLE = joinpath(RECONR_SERIAL_T84_DIR, "referenceTape100")
const RECONR_SERIAL_T84_WORK = "/tmp/njoy_gkw_t84_serialization"

function first_record_difference(reference::Vector{String}, trial::Vector{String})
    common = min(length(reference), length(trial))
    mismatch = findfirst(i -> reference[i] != trial[i], 1:common)
    mismatch !== nothing && return mismatch
    length(reference) == length(trial) ? nothing : common + 1
end

@testset "T84 referenceTape100 lines 1-1115 are byte-identical" begin
    rm(RECONR_SERIAL_T84_WORK; recursive=true, force=true)
    run_njoy(joinpath(RECONR_SERIAL_T84_DIR, "input");
             work_dir=RECONR_SERIAL_T84_WORK, verbose=false)

    trial_path = joinpath(RECONR_SERIAL_T84_WORK, "tape100")
    @test isfile(trial_path)

    reference = readlines(RECONR_SERIAL_T84_ORACLE)
    trial = readlines(trial_path)
    @test length(trial) == length(reference) == 1115
    @test first_record_difference(reference, trial) === nothing
end
