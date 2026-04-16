using Test
using NJOY

# TDD grind for T34 MT=102 (Pu-240 n,γ capture, Reich-Moore + URR LSSF=1)
#
# Pre-grind state: 52/53 PERFECT.  1 MT differs (MT=102) with 3 lines
# showing ±1 in the 7th sigfig at E=630.04, 2089.07, 4526.46 eV.  Phase
# 14 gdb investigation labelled this "CONFIRMED IRREDUCIBLE" FP
# accumulation noise in the l=1 non-fissile Reich-Moore small-parameter
# path.  Per project insight #4 every prior irreducible claim has
# turned out to be a real bug.
#
# Goal: MT=102 bit-identical with Fortran oracle.

const T34_ORACLE = joinpath(@__DIR__, "oracle_cache", "test34", "after_reconr.pendf")
const T34_INPUT  = joinpath(@__DIR__, "oracle_cache", "test34", "run_reconr", "tape20")
const T34_JULIA  = "/tmp/t34_julia.pendf"

function t34_parse_mf3(fn, mt_target)
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

function ensure_t34()
    sammy_path = joinpath(@__DIR__, "..", "..", "src", "resonances", "sammy.jl")
    rm_path    = joinpath(@__DIR__, "..", "..", "src", "resonances", "reich_moore.jl")
    if !isfile(T34_JULIA) ||
       mtime(T34_JULIA) < mtime(T34_INPUT) ||
       mtime(T34_JULIA) < mtime(sammy_path) ||
       mtime(T34_JULIA) < mtime(rm_path)
        r = reconr(T34_INPUT; mat=9440, err=0.001)
        write_pendf_file(T34_JULIA, r; mat=9440, err=0.001)
    end
end

@testset "T34 MT=102 (Pu-240 n,γ) — bit-identical grind" begin
    ensure_t34()
    j = t34_parse_mf3(T34_JULIA, 102)
    f = t34_parse_mf3(T34_ORACLE, 102)
    @test length(j) == length(f)
    ndiff = 0
    first_diff = 0
    for i in 1:min(length(j), length(f))
        if j[i] != f[i]
            ndiff += 1
            first_diff == 0 && (first_diff = i)
        end
    end
    println("MT=102: ", ndiff, " / ", length(j), " lines differ; first at line ", first_diff)
    @test ndiff == 0
end
