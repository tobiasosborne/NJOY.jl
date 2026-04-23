using Test
using NJOY

# Fortran ref: njoy-reference/src/leapr.f90:230-400 (deck read, card-by-card)
# T22 fixture: njoy-reference/tests/22/input (para-H2 at 20 K)

const T22_INPUT = read(joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "22", "input"), String)

@testset "parse_leapr T22 (para-H2 at 20K)" begin
    calls = NJOY.tokenise_njoy_input(T22_INPUT)
    leapr_mc = first(c for c in calls if c.name == :leapr)
    p = NJOY.parse_leapr(leapr_mc)

    # Card 1: global output unit
    @test p.nout == 20

    # Card 2: title
    @test occursin("para", lowercase(p.title))
    @test occursin("20k",  lowercase(p.title))

    # Card 3: run control (iprint=2 provided; nphon defaults to 100)
    @test p.ntempr == 1
    @test p.iprint == 2
    @test p.nphon  == 100

    # Card 4: ENDF metadata
    @test p.mat   == 2
    @test p.za    == 1001
    @test p.isabt == 0
    @test p.ilog  == 0
    @test p.smin  ≈ 2e-38 rtol=1e-12

    # Card 5: principal scatterer (nsk omitted -> 0)
    @test p.awr   ≈ 0.99917 rtol=1e-6
    @test p.spr   ≈ 20.478  rtol=1e-6
    @test p.npr   == 2
    @test p.iel   == 0
    @test p.ncold == 2
    @test p.nsk   == 0

    # Card 6: no secondary scatterer
    @test p.nss == 0

    # Card 7: alpha/beta control
    @test p.nalpha == 59
    @test p.nbeta  == 105
    @test p.lat    == 0

    # Cards 8-9: grids
    @test length(p.alpha) == 59
    @test length(p.beta)  == 105
    @test p.alpha[1]  ≈ 0.001 rtol=1e-9
    @test p.alpha[22] ≈ 10.0  rtol=1e-9
    @test p.alpha[59] ≈ 600.0 rtol=1e-9
    @test p.beta[1]   == 0.0
    @test p.beta[2]   ≈ 0.01  rtol=1e-9
    @test p.beta[105] ≈ 300.0 rtol=1e-9

    # Per-temperature data (1 entry for ntempr=1)
    @test length(p.temperatures) == 1
    @test p.temperatures[1] ≈ 20.0 rtol=1e-9

    # Card 11: delta1, ni
    @test length(p.delta1) == 1
    @test length(p.ni)     == 1
    @test p.delta1[1] ≈ 2.5e-4 rtol=1e-9
    @test p.ni[1] == 48

    # Card 12: phonon DOS (48 values)
    @test length(p.p1_dos)    == 1
    @test length(p.p1_dos[1]) == 48
    @test p.p1_dos[1][1]  == 0.0
    @test p.p1_dos[1][2]  ≈ 0.01563 rtol=1e-6
    @test p.p1_dos[1][48] == 0.0

    # Card 13: translational params
    @test p.twt[1]   ≈ 0.025 rtol=1e-9
    @test p.c[1]     ≈ 40.0  rtol=1e-9
    @test p.tbeta[1] ≈ 0.475 rtol=1e-9

    # Card 14: no discrete oscillators
    @test p.nd[1] == 0
    @test isempty(p.bdel[1])
    @test isempty(p.adel[1])

    # Cards 17-18: pair-correlation s(kappa) (read because ncold>0)
    @test p.nka[1] == 301
    @test p.dka[1] ≈ 0.05 rtol=1e-9
    @test length(p.ska[1])  == 301
    @test p.ska[1][1]   ≈ 9.29815e-2 rtol=1e-6
    @test p.ska[1][301] ≈ 9.99400e-1 rtol=1e-6

    # Card 19: cfrac (NOT read because nsk=0)
    @test isempty(p.cfrac) || p.cfrac[1] == 0.0

    # Comments
    @test length(p.comments) >= 1
    @test any(occursin("para-h", lowercase(c)) for c in p.comments)
end

const T23_INPUT = read(joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "23", "input"), String)

@testset "parse_leapr T23 (BeO, 8 temps, secondary scatterer)" begin
    calls = NJOY.tokenise_njoy_input(T23_INPUT)
    leapr_mc = first(c for c in calls if c.name == :leapr)
    p = NJOY.parse_leapr(leapr_mc)

    @test p.nout == 20
    @test occursin("beo", lowercase(p.title))

    # Card 3: 8 temperatures, iprint=1, nphon default 100
    @test p.ntempr == 8
    @test p.iprint == 1
    @test p.nphon  == 100

    # Card 4
    @test p.mat == 27
    @test p.za  == 127
    @test p.smin ≈ 2e-38 rtol=1e-12

    # Card 5: principal = Be in BeO (ncold/nsk omitted → 0)
    @test p.awr   ≈ 8.93478 rtol=1e-6
    @test p.spr   ≈ 6.15    rtol=1e-6
    @test p.npr   == 1
    @test p.iel   == 3      # BeO coherent elastic model
    @test p.ncold == 0
    @test p.nsk   == 0

    # Card 6: secondary = O in BeO
    @test p.nss == 1
    @test p.b7  ≈ 0.0
    @test p.aws ≈ 15.858 rtol=1e-6
    @test p.sps ≈ 3.7481 rtol=1e-6
    @test p.mss == 1

    # Card 7: nalpha=90, nbeta=117, lat=1
    @test p.nalpha == 90
    @test p.nbeta  == 117
    @test p.lat    == 1
    @test length(p.alpha) == 90
    @test length(p.beta)  == 117
    @test p.alpha[1] ≈ 0.01008 rtol=1e-6
    @test p.alpha[90] ≈ 80.0 rtol=1e-9
    @test p.beta[1] == 0.0
    @test p.beta[117] ≈ 80.0 rtol=1e-9

    # Temperatures: 296, -400, -500, -600, -700, -800, -1000, -1200 (abs in struct)
    @test length(p.temperatures) == 8
    @test p.temperatures[1] == 296.0
    @test p.temperatures[2] == 400.0
    @test p.temperatures[8] == 1200.0

    # First temperature carries the real data; subsequent temps reuse
    @test p.ni[1] == 84
    @test length(p.p1_dos[1]) == 84
    @test p.ni[2] == 84  # reused from temp 1
    @test p.p1_dos[2] == p.p1_dos[1]  # reused

    @test p.twt[1] == 0.0
    @test p.c[1]   == 0.0
    @test p.tbeta[1] ≈ 1.0 rtol=1e-9

    @test p.nd[1] == 0
end

const T80_INPUT = read(joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "80", "input"), String)

@testset "parse_leapr T80 (H_HF, large grids, bare title)" begin
    calls = NJOY.tokenise_njoy_input(T80_INPUT)
    leapr_mc = first(c for c in calls if c.name == :leapr)
    p = NJOY.parse_leapr(leapr_mc)

    @test p.nout == 24
    # "H_HF/" title — bare word, no quote delimiter
    @test occursin("h_hf", lowercase(p.title))

    # Card 3: ntempr=1, iprint=2
    @test p.ntempr == 1
    @test p.iprint == 2

    # Card 4
    @test p.mat == 8
    @test p.za  == 108

    # Card 5: H in HF, ncold/nsk omitted → 0 (no cold-H convolution or s(kappa))
    @test p.awr   ≈ 0.9991673 rtol=1e-6
    @test p.spr   ≈ 20.43608  rtol=1e-6
    @test p.npr   == 1
    @test p.iel   == 0
    @test p.ncold == 0
    @test p.nsk   == 0

    @test p.nss == 0

    # Card 7: very fine grids
    @test p.nalpha == 200
    @test p.nbeta  == 1325
    @test p.lat    == 1
    @test length(p.alpha) == 200
    @test length(p.beta)  == 1325

    # Single temperature
    @test length(p.temperatures) == 1
    @test p.temperatures[1] ≈ 343.0 rtol=1e-9

    # Card 11: fine DOS
    @test p.delta1[1] ≈ 1.240688e-3 rtol=1e-6
    @test p.ni[1]     == 551
    @test length(p.p1_dos[1]) == 551

    # No cold-H convolution, no s(kappa) cards
    @test p.nka[1] == 0
    @test isempty(p.ska[1])
end
