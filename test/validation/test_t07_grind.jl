using Test
using NJOY

# TDD grind for T07 (U-235, SLBW + URR mode=12)
#
# Pre-grind state: 24/27 PERFECT.  3 MTs (18, 19, 102) differ at the
# URR-boundary grid points sigfig(82, 7, ±1) and sigfig(25000, 7, -1).
# HANDOFF labels this "FP accumulation in _gnrl Gauss-Laguerre".
# Per project insight #4, every prior such label has been a real bug.

const T07_ORACLE = joinpath(@__DIR__, "oracle_cache", "test07", "after_reconr.pendf")
const T07_INPUT  = joinpath(@__DIR__, "oracle_cache", "test07", "run_reconr", "tape20")
const T07_JULIA  = "/tmp/t07_julia.pendf"

function t07_parse_mf3(fn, mt_target)
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

function ensure_t07()
    src = joinpath(@__DIR__, "..", "..", "src")
    if !isfile(T07_JULIA) ||
       mtime(T07_JULIA) < mtime(T07_INPUT) ||
       mtime(T07_JULIA) < mtime(joinpath(src, "resonances", "unresolved.jl")) ||
       mtime(T07_JULIA) < mtime(joinpath(src, "processing", "reconr.jl"))
        r = reconr(T07_INPUT; mat=1395, err=0.005)
        write_pendf_file(T07_JULIA, r; mat=1395, err=0.005)
    end
end

@testset "T07 U-235 — bit-identical grind" begin
    ensure_t07()
    for mt in (18, 19, 102)
        j = t07_parse_mf3(T07_JULIA, mt)
        f = t07_parse_mf3(T07_ORACLE, mt)
        @test length(j) == length(f)
        ndiff = 0
        for i in 1:min(length(j), length(f))
            j[i] != f[i] && (ndiff += 1)
        end
        println("MT=$mt: $ndiff / $(length(j)) lines differ")
        @test ndiff == 0
    end
end
