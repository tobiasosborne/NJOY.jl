# RED→GREEN unit test for thermr incoherent-elastic (iel / LTHR=2) — Stage 1.
#
# Wires the deterministic pure-function layer of NJOY_jl-69j:
#   - incoh_elastic_xs(E, sigma_b, dwp, natom): the integrated incoherent-elastic
#     cross section (Fortran thermr.f90:1322,1379-1384, subroutine iel).
#   - read_mf7_mt2 LTHR=2/LTHR=3 incoherent branch returning an
#     IncoherentElasticData object (SB + W'(T) table).
#   - incoh_dwp_at(iel_data, T): W'(T) interpolation (thermr.f90:1303-1318).
#   - build_incoh_elastic_records(esi, sigma_b, dwp, nbin, natom): the iel
#     record builder (xie + LANG=3 equiprobable-cosine MF6 records, :1377-1406).
#
# Fortran ground truth (subroutine iel, thermr.f90:1244-1425):
#   c1 = sb/(2*natom)                       (line 1322)
#   c2 = 2*e*dwa                            (line 1379)
#   x1 = exp(-2*c2)                         (line 1382)
#   rc2 = 1/c2                              (line 1381)
#   xsec = c1*rc2*(1-x1)                    (line 1384)
#     => σ(E) = sb/(2*natom) * 1/(2*E*W') * (1 - exp(-4*E*W'))
#
# Exact anchors parsed from the ENDF resources (natom=1 for T74/T69, =2 for T67):
#   T74 H-in-ZrH  (MAT=7,  LTHR=2): SB=81.98006, W'(296)=8.486993
#       σ(1e-5, 296) = 81.96615 barn
#   T69 Zr-in-ZrH (MAT=58, LTHR=2): SB=6.337872, W'(400)=2.677764
#       σ(1e-5, 400) = 6.337533 barn
#   T67 D-in-LiD  (MAT=15, LTHR=3): SB(incoh)=2.054202, natom=2, W'(400)=10.47985
#       σ(1e-5, 400) = 1.026886 barn
#
# Run:  julia --project=. test/validation/test_thermr_incoherent_elastic.jl

using Test
using NJOY

const RES = joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "resources")

@testset "thermr incoherent elastic (iel / LTHR=2) — Stage 1" begin

    @testset "incoh_elastic_xs integrated formula (thermr.f90:1322,1379-1384)" begin
        # T74 H-in-ZrH at T=296, E=1e-5
        @test isapprox(NJOY.incoh_elastic_xs(1e-5, 81.98006, 8.486993, 1),
                       81.96615; rtol=1e-5)
        # T69 Zr-in-ZrH at T=400, E=1e-5
        @test isapprox(NJOY.incoh_elastic_xs(1e-5, 6.337872, 2.677764, 1),
                       6.337533; rtol=1e-5)
        # T67 D-in-LiD at T=400, E=1e-5, natom=2
        @test isapprox(NJOY.incoh_elastic_xs(1e-5, 2.054202, 10.47985, 2),
                       1.026886; rtol=1e-5)
        # E<=0 guard
        @test NJOY.incoh_elastic_xs(0.0, 81.98006, 8.486993, 1) == 0.0
        @test NJOY.incoh_elastic_xs(-1.0, 81.98006, 8.486993, 1) == 0.0
        # High-E saturation: σ -> sb/(2*natom)*1/(2*E*W') for large E*W'
        big = NJOY.incoh_elastic_xs(100.0, 81.98006, 8.486993, 1)
        @test isapprox(big, 81.98006 / (2 * 1) / (2 * 100.0 * 8.486993); rtol=1e-6)
    end

    @testset "read_mf7_mt2 LTHR=2 incoherent (T74 H-in-ZrH, MAT=7, T=296)" begin
        f = joinpath(RES, "tsl-HinZrH-ENDF8.0.endf")
        lthr, bragg, incoh = NJOY.read_mf7_mt2(f, 7, 296.0)
        @test lthr == 2
        @test bragg === nothing
        @test incoh !== nothing
        @test isapprox(incoh.sigma_b, 81.98006; rtol=1e-6)
        # W'(296) is a table point -> exact
        @test isapprox(NJOY.incoh_dwp_at(incoh, 296.0), 8.486993; rtol=1e-6)
        # σ via the full chain at the anchor
        σ = NJOY.incoh_elastic_xs(1e-5, incoh.sigma_b,
                                  NJOY.incoh_dwp_at(incoh, 296.0), 1)
        @test isapprox(σ, 81.96615; rtol=1e-5)
    end

    @testset "incoh_dwp_at interpolation (thermr.f90:1303-1318)" begin
        f = joinpath(RES, "tsl-HinZrH-ENDF8.0.endf")
        _, _, incoh = NJOY.read_mf7_mt2(f, 7, 296.0)
        # exact table points
        @test isapprox(NJOY.incoh_dwp_at(incoh, 296.0), 8.486993; rtol=1e-7)
        @test isapprox(NJOY.incoh_dwp_at(incoh, 400.0), 9.093191; rtol=1e-7)
        @test isapprox(NJOY.incoh_dwp_at(incoh, 1200.0), 17.14050; rtol=1e-6)
        # lin-lin between 400 and 500 (9.093191 -> 9.828159) at T=450
        expected = 9.093191 + (450.0 - 400.0) * (9.828159 - 9.093191) /
                   (500.0 - 400.0)
        @test isapprox(NJOY.incoh_dwp_at(incoh, 450.0), expected; rtol=1e-7)
    end

    @testset "read_mf7_mt2 LTHR=2 (T69 Zr-in-ZrH, MAT=58, T=400)" begin
        f = joinpath(RES, "tsl-ZrinZrH-ENDF8.0.endf")
        lthr, bragg, incoh = NJOY.read_mf7_mt2(f, 58, 400.0)
        @test lthr == 2
        @test bragg === nothing
        @test incoh !== nothing
        @test isapprox(incoh.sigma_b, 6.337872; rtol=1e-6)
        @test isapprox(NJOY.incoh_dwp_at(incoh, 400.0), 2.677764; rtol=1e-6)
        σ = NJOY.incoh_elastic_xs(1e-5, incoh.sigma_b,
                                  NJOY.incoh_dwp_at(incoh, 400.0), 1)
        @test isapprox(σ, 6.337533; rtol=1e-5)
    end

    @testset "read_mf7_mt2 LTHR=3 mixed (T67 D-in-LiD, MAT=15, T=400)" begin
        f = joinpath(RES, "tsl-DinLiD-ENDF8.1-beta.endf")
        lthr, bragg, incoh = NJOY.read_mf7_mt2(f, 15, 400.0)
        @test lthr == 3
        @test bragg !== nothing            # coherent part still parsed
        @test incoh !== nothing            # incoherent W'(T) part parsed
        @test isapprox(incoh.sigma_b, 2.054202; rtol=1e-6)
        @test isapprox(NJOY.incoh_dwp_at(incoh, 400.0), 10.47985; rtol=1e-6)
        σ = NJOY.incoh_elastic_xs(1e-5, incoh.sigma_b,
                                  NJOY.incoh_dwp_at(incoh, 400.0), 2)
        @test isapprox(σ, 1.026886; rtol=1e-5)
    end

    @testset "build_incoh_elastic_records (iel record builder, :1377-1406)" begin
        nbin = 8
        natom = 1
        sigma_b = 81.98006
        dwp = 8.486993
        esi = Float64[1e-5, 1e-3, 1e-2, 1e-1, 1.0]
        xie, recs = NJOY.build_incoh_elastic_records(esi, sigma_b, dwp, nbin, natom)

        # xie matches the σ anchor at esi[1]=1e-5 (raw, no sigfig — line 1400)
        @test length(xie) == length(esi)
        @test isapprox(xie[1], 81.96615; rtol=1e-5)
        @test all(xie[i] == NJOY.incoh_elastic_xs(esi[i], sigma_b, dwp, natom)
                  for i in eachindex(esi))

        # one MF6 record per incident energy
        @test length(recs) == length(esi)
        for (i, r) in enumerate(recs)
            @test r isa NJOY.MF6ListRecord
            @test isapprox(r.E_incident, esi[i]; rtol=0, atol=0)
            @test length(r.entries) == 1               # one secondary energy (elastic)
            entry = r.entries[1]
            @test length(entry) == nbin + 2            # [E', weight, μ₁...μ_nbin]
            @test entry[1] == esi[i]                   # E' = E (elastic)
            @test entry[2] == 1.0                      # weight word (scr(8)=1)
            cosines = entry[3:end]
            @test length(cosines) == nbin
            @test all(-1.0 - 1e-9 .<= cosines .<= 1.0 + 1e-9)
        end
    end
end
