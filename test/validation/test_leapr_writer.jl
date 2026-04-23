using Test
using NJOY

# Fortran ref: njoy-reference/src/leapr.f90:2972-3623 (endout).

const T22_INPUT = read(joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "22", "input"), String)

@testset "write_leapr_tape structural (T22)" begin
    calls = NJOY.tokenise_njoy_input(T22_INPUT)
    leapr_mc = first(c for c in calls if c.name == :leapr)
    p = NJOY.parse_leapr(leapr_mc)

    # Synthesize plausible S(α,β) for smoke test — values don't have to be
    # physically correct at this stage; the test checks tape STRUCTURE is
    # readable by the existing NJOY.jl MF7 reader.
    ssm = zeros(p.nbeta, p.nalpha, p.ntempr)
    ssp = zeros(p.nbeta, p.nalpha, p.ntempr)
    for it in 1:p.ntempr, j in 1:p.nalpha, i in 1:p.nbeta
        ssm[i, j, it] = exp(-p.alpha[j] - 0.1*p.beta[i])
        ssp[i, j, it] = exp(-p.alpha[j] - 0.1*p.beta[i])
    end

    tempf = copy(p.temperatures)
    dwpix = ones(p.ntempr)

    tmppath = tempname() * ".endf"
    NJOY.write_leapr_tape(tmppath, p, ssm, ssp;
                          tempf=tempf, dwpix=dwpix, isym=1)
    @test isfile(tmppath)

    # File has reasonable line count (T22 reference: 4636 lines)
    lines = readlines(tmppath)
    # Per-β block = 22 lines (TAB1 header+interp+ceil(59/3)=20). 209 rows → 4598.
    # Plus MF1/MT451 ~24 lines + MF7/MT4 header/LIST/TAB2 ~6 + tail ~5 + SEND/FEND/MEND/TEND = ~4640
    @test 4500 < length(lines) < 5000

    # TPID present on line 1
    @test length(lines[1]) >= 66
    @test lstrip(lines[1][67:75]) in ("1 0  0", "1 0 0")  # MAT=1 MF=0 MT=0

    # MF1/MT451 HEAD1 on line 2: ZA=1001
    @test occursin("1.001000+3", lines[2])
    @test occursin(" 2 1451", lines[2])

    # MF7/MT4 section present (MF at cols 71-72, MT at 73-75)
    mf7_lines = [l for l in lines
                 if length(l) >= 75 && strip(l[71:72]) == "7" && strip(l[73:75]) == "4"]
    @test length(mf7_lines) > 100

    # MEND + TEND at the end (last non-empty 2 lines)
    @test occursin("   0 0  0", lines[end-1]) || occursin("0 0  0    0", lines[end-1])
    @test occursin("  -1 0  0", lines[end]) || occursin("-1 0  0    0", lines[end])
end
