using Test
using NJOY

# T84 H-2 is the bundled oracle for MF12 TAB1 metadata preservation.
# Ref: referenceTape100 line 853 and reconr.f90:4924-4927: emerge rewrites
# NR/NP and interpolation metadata, but preserves the TAB1 L1/L2 fields.
const T84_DIR = joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "84")
const T84_ORACLE = joinpath(T84_DIR, "referenceTape100")
const T84_WORK_DIR = "/tmp/njoy_t84_photon_metadata"

function t84_mf12_columns(path::AbstractString)
    records = String[]
    for line in eachline(path)
        padded = rpad(line, 80)
        mat = NJOY._parse_int(padded[67:70])
        mf = NJOY._parse_int(padded[71:72])
        mt = NJOY._parse_int(padded[73:75])
        mat == 128 && mf == 12 && mt == 102 && push!(records, padded[1:66])
    end
    records
end

@testset "T84 MF12/MT102 preserves TAB1 metadata at oracle line 853" begin
    rm(T84_WORK_DIR; recursive=true, force=true)
    run_njoy(joinpath(T84_DIR, "input"); work_dir=T84_WORK_DIR, verbose=false)
    trial = joinpath(T84_WORK_DIR, "tape100")
    @test isfile(trial)
    @test t84_mf12_columns(trial) == t84_mf12_columns(T84_ORACLE)

    tab1 = rpad(readlines(trial)[findfirst(line ->
        length(line) >= 75 && NJOY._parse_int(rpad(line, 80)[71:72]) == 12 &&
        NJOY._parse_int(rpad(line, 80)[73:75]) == 102,
        readlines(trial)) + 1], 80)
    @test (NJOY._parse_int(tab1[23:33]), NJOY._parse_int(tab1[34:44])) == (2, 2)
end
