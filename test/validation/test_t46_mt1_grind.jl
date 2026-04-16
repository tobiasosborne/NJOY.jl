using Test
using NJOY

# TDD grind for T46 MT=1 (Fe-56 JEFF-3.3 total XS)
#
# Pre-grind state: 72/73 PERFECT.  MT=1 total differs from Fortran oracle
# at 2 high-energy lines (±1-3 in 7th sigfig).  Per HANDOFF this is
# "IEEE 754 summation order" — Fortran accumulates each sigfig'd section
# in tape order, Julia sums first then rounds.
#
# Goal: MT=1 bit-identical.

const T46_ORACLE = joinpath(@__DIR__, "oracle_cache", "test46", "after_reconr.pendf")
const T46_INPUT  = joinpath(@__DIR__, "oracle_cache", "test46", "run_reconr", "tape20")
const T46_JULIA  = "/tmp/t46_julia.pendf"

function t46_parse_mf3(fn, mt_target)
    lines = readlines(fn); idx = 1
    data = String[]
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = NJOY._parse_int(p[71:72]); mt = NJOY._parse_int(p[73:75])
        mat_line = NJOY._parse_int(p[67:70])
        if mf == 3 && mt == mt_target && mat_line > 0
            idx += 3
            while idx <= length(lines)
                pd = rpad(lines[idx], 80)
                mfc = NJOY._parse_int(pd[71:72]); mtc = NJOY._parse_int(pd[73:75])
                (mfc != 3 || mtc != mt) && break
                push!(data, pd[1:66]); idx += 1
            end
            return data
        else; idx += 1; end
    end
    return data
end

function ensure_t46()
    src = joinpath(@__DIR__, "..", "..", "src")
    if !isfile(T46_JULIA) ||
       mtime(T46_JULIA) < mtime(T46_INPUT) ||
       mtime(T46_JULIA) < mtime(joinpath(src, "processing", "reconr.jl")) ||
       mtime(T46_JULIA) < mtime(joinpath(src, "processing", "reconr_evaluator.jl")) ||
       mtime(T46_JULIA) < mtime(joinpath(src, "processing", "pendf_writer.jl"))
        r = reconr(T46_INPUT; mat=2631, err=0.001)
        write_pendf_file(T46_JULIA, r; mat=2631, err=0.001)
    end
end

@testset "T46 MT=1 (Fe-56 JEFF total) — bit-identical grind" begin
    ensure_t46()
    j = t46_parse_mf3(T46_JULIA, 1)
    f = t46_parse_mf3(T46_ORACLE, 1)
    @test length(j) == length(f)
    ndiff = 0
    first_diff = 0
    for i in 1:min(length(j), length(f))
        if j[i] != f[i]
            ndiff += 1
            first_diff == 0 && (first_diff = i)
        end
    end
    println("MT=1: ", ndiff, " / ", length(j), " lines differ; first at line ", first_diff)
    @test ndiff == 0
end
