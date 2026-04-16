using Test
using NJOY

# TDD grind for T20 MT=600 (Cl-35 n,p proton channel, SAMMY/RML)
#
# Pre-grind state (2026-04-16):
#   MT=600: 3577 lines, 1973 differ. Pattern — diffs concentrated in
#   [100 keV, 1.2 MeV] with systematic -1.8e-5 relative bias peaking
#   near 1 MeV and vanishing above 1.5 MeV. NOT FP noise.
#
# Goal: bit-identical MT=600 with Fortran oracle.

const ORACLE = joinpath(@__DIR__, "oracle_cache", "test20", "after_reconr.pendf")
const INPUT  = joinpath(@__DIR__, "oracle_cache", "test20", "run_reconr", "tape20")
const JULIA_OUT = "/tmp/t20_julia.pendf"

function parse_mf3_lines(fn::AbstractString, mt_target::Int)
    lines = readlines(fn); idx = 1
    data = String[]
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = NJOY._parse_int(p[71:72]); mt = NJOY._parse_int(p[73:75])
        mat = NJOY._parse_int(p[67:70])
        if mf == 3 && mt == mt_target && mat > 0
            idx += 3
            while idx <= length(lines)
                pd = rpad(lines[idx], 80)
                mfc = NJOY._parse_int(pd[71:72]); mtc = NJOY._parse_int(pd[73:75])
                (mfc != 3 || mtc != mt) && break
                push!(data, pd[1:66])
                idx += 1
            end
            return data
        else
            idx += 1
        end
    end
    return data
end

function ensure_julia_output()
    if !isfile(JULIA_OUT) ||
       mtime(JULIA_OUT) < mtime(INPUT) ||
       mtime(JULIA_OUT) < mtime(joinpath(@__DIR__, "..", "..", "src", "resonances", "sammy.jl"))
        r = reconr(INPUT; mat=1725, err=0.01)
        write_pendf_file(JULIA_OUT, r; mat=1725, err=0.01)
    end
end

@testset "T20 MT=600 (Cl-35 n,p) — bit-identical grind" begin
    ensure_julia_output()
    j = parse_mf3_lines(JULIA_OUT, 600)
    f = parse_mf3_lines(ORACLE,    600)
    @test length(j) == length(f)
    ndiff = 0
    first_diff = 0
    for i in 1:min(length(j), length(f))
        if j[i] != f[i]
            ndiff += 1
            first_diff == 0 && (first_diff = i)
        end
    end
    println("MT=600: ", ndiff, " / ", length(j), " lines differ; first at line ", first_diff)
    @test ndiff == 0
end
