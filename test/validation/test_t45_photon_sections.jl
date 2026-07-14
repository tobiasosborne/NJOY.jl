using Test
using NJOY

# T45 regression for NJOY_jl-2k4.
#
# Oracle slices (one-based referenceTape40 lines, excluding SEND):
#   MF12/MT102  6639-6979
#   MF13/MT4    6982-7104
#   MF13/MT103  7106-7184
#
# Ref: njoy-reference/src/reconr.f90:1771-2197 (lunion), 4726-4960
# (emerge). RECONR retains only each section's total TAB1, evaluates it on
# the union grid, and deliberately omits MF14.

const T45_DIR = joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "45")
const T45_ORACLE = joinpath(T45_DIR, "referenceTape40")
const T45_WORK_DIR = "/tmp/njoy_t45_photon_regression"

function t45_section_columns(path::AbstractString, mf_target::Int, mt_target::Int)
    records = String[]
    for line in eachline(path)
        padded = rpad(line, 80)
        mat = NJOY._parse_int(padded[67:70])
        mf = NJOY._parse_int(padded[71:72])
        mt = NJOY._parse_int(padded[73:75])
        mat > 0 && mf == mf_target && mt == mt_target && push!(records, padded[1:66])
    end
    records
end

function t45_has_mf(path::AbstractString, mf_target::Int)
    any(eachline(path)) do line
        padded = rpad(line, 80)
        NJOY._parse_int(padded[67:70]) > 0 &&
            NJOY._parse_int(padded[71:72]) == mf_target
    end
end

function t45_directory_entries(path::AbstractString)
    entries = NTuple{4,Int}[]
    for line in eachline(path)
        padded = rpad(line, 80)
        NJOY._parse_int(padded[67:70]) == 525 || continue
        NJOY._parse_int(padded[71:72]) == 1 || continue
        NJOY._parse_int(padded[73:75]) == 451 || continue
        mf = NJOY._parse_int(padded[23:33])
        mt = NJOY._parse_int(padded[34:44])
        mf in (12, 13, 14) || continue
        push!(entries, (mf, mt, NJOY._parse_int(padded[45:55]),
                        NJOY._parse_int(padded[56:66])))
    end
    entries
end

function t45_photon_structure(path::AbstractString)
    events = Tuple{Symbol,Int,Int}[]
    last_section = (0, 0)
    current_mf = 0
    for line in eachline(path)
        padded = rpad(line, 80)
        mat = NJOY._parse_int(padded[67:70])
        mf = NJOY._parse_int(padded[71:72])
        mt = NJOY._parse_int(padded[73:75])
        if mat == 525 && mf in (12, 13) && mt > 0
            current_mf = mf
            if (mf, mt) != last_section
                push!(events, (:section, mf, mt))
                last_section = (mf, mt)
            end
        elseif mat == 525 && mf in (12, 13) && mt == 0
            push!(events, (:send, mf, 0))
        elseif mat == 525 && mf == 0 && mt == 0 && current_mf in (12, 13)
            push!(events, (:fend, current_mf, 0))
            current_mf = 0
        end
    end
    events
end

@testset "T45 photon totals match oracle fixed columns" begin
    rm(T45_WORK_DIR; recursive=true, force=true)
    run_njoy(joinpath(T45_DIR, "input"); work_dir=T45_WORK_DIR, verbose=false)
    trial = joinpath(T45_WORK_DIR, "tape40")
    @test isfile(trial)

    @testset "MF12/MT102 oracle lines 6639-6979" begin
        @test t45_section_columns(trial, 12, 102) ==
              t45_section_columns(T45_ORACLE, 12, 102)
    end
    @testset "MF13/MT4 oracle lines 6982-7104" begin
        @test t45_section_columns(trial, 13, 4) ==
              t45_section_columns(T45_ORACLE, 13, 4)
    end
    @testset "MF13/MT103 oracle lines 7106-7184" begin
        @test t45_section_columns(trial, 13, 103) ==
              t45_section_columns(T45_ORACLE, 13, 103)
    end

    # MF14 is present in the source ENDF, but anlyzd/lunion exclude it from
    # the reconstructed PENDF. Its absence is canonical, not missing work.
    @test !t45_has_mf(trial, 14)
    @test t45_directory_entries(trial) ==
          [(12, 102, 341, 0), (13, 4, 123, 0), (13, 103, 79, 0)]
    @test t45_photon_structure(trial) == [
        (:section, 12, 102), (:send, 12, 0), (:fend, 12, 0),
        (:section, 13, 4), (:send, 13, 0),
        (:section, 13, 103), (:send, 13, 0), (:fend, 13, 0),
    ]
end
