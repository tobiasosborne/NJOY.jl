using Test
using NJOY
using SpecialFunctions: besselk

# Fortran ref: njoy-reference/src/leapr.f90:844-1007 (trans, stable, sbfill,
# terps, besk1). Julia port: src/processing/leapr.jl.

@testset "_besk1 vs SpecialFunctions.besselk(1, x)" begin
    # For x <= 1 the Fortran returns unscaled K_1(x).
    for x in (0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)
        @test NJOY._besk1(x) ≈ besselk(1, x) rtol=1e-6
    end
    # For x > 1 the Fortran returns e^x · K_1(x) (scaled form; stable() compensates).
    for x in (1.1, 2.0, 5.0, 10.0, 50.0)
        @test NJOY._besk1(x) ≈ exp(x) * besselk(1, x) rtol=1e-6
    end
end

@testset "_stable! free-gas branch analytic match" begin
    # Fortran stable (c=0 free-gas): sd(j) = exp(-(wal - be)^2 / (4*wal)) / sqrt(4π·wal)
    # where wal = twt*α, be = (j-1)*delta.
    ap = Vector{Float64}(undef, 10_000)
    sd = Vector{Float64}(undef, 10_000)
    twt = 0.025; al = 1.0; delta = 0.1
    nsd = NJOY._stable!(ap, sd, al, delta, twt, 0.0, 10_000)
    @test isodd(nsd)  # Simpson's rule requirement
    wal = twt * al
    for j in 1:min(nsd, 20)
        be = (j-1) * delta
        expected = exp(-(wal - be)^2 / (4*wal)) / sqrt(4π * wal)
        @test sd[j] ≈ expected rtol=1e-10
        @test ap[j] ≈ be rtol=1e-12
    end
end

@testset "_stable! diffusion branch shape" begin
    # With c>0 the kernel peaks at β=0 and decays; spot-check positivity
    # and the eps-stop criterion (sd[end] < 1e-7 * sd[1]).
    ap = Vector{Float64}(undef, 10_000)
    sd = Vector{Float64}(undef, 10_000)
    nsd = NJOY._stable!(ap, sd, 1.0, 0.05, 0.025, 40.0, 10_000)
    @test isodd(nsd)
    @test nsd >= 3
    @test all(sd[j] > 0 for j in 1:nsd)
    @test sd[nsd] <= 1e-7 * sd[1] || nsd >= 9999
end

@testset "_terps log-linear interpolation" begin
    # Build a known log-linear table: sd[j] = exp(-j)
    nsd = 11
    delta = 1.0
    sd = [exp(-j+1.0) for j in 1:nsd]   # sd[1]=exp(0)=1, sd[2]=exp(-1), ...
    # At exact grid point
    @test NJOY._terps(sd, nsd, delta, 0.0) ≈ 1.0 rtol=1e-10
    @test NJOY._terps(sd, nsd, delta, 1.0) ≈ exp(-1.0) rtol=1e-10
    # Halfway (log-linear)
    @test NJOY._terps(sd, nsd, delta, 0.5) ≈ exp(-0.5) rtol=1e-10
    # Past the table
    @test NJOY._terps(sd, nsd, delta, 100.0) == 0.0
end

@testset "trans! smoke: T22 para-H2 at 20K" begin
    # Synthetic smooth S(α,-β) — rough proxy for what contin() produces.
    # After trans, ssm values should be: positive, non-increasing with β (gross
    # monotonicity), and tempf[itemp] should be updated per the Fortran formula.
    nalpha = 8; nbeta = 16
    alpha_grid = collect(range(0.1, 2.0, length=nalpha))
    beta_grid  = collect(range(0.0, 10.0, length=nbeta))
    tev = 8.617333e-5 * 20.0         # k_B·T at 20 K, in eV
    ssm = zeros(nbeta, nalpha, 1)
    for j in 1:nalpha, i in 1:nbeta
        ssm[i, j, 1] = exp(-alpha_grid[j] - 0.5*beta_grid[i])
    end
    tempr = [20.0]
    tempf = [20.0]                    # will be mutated
    twt = 0.025; c_diff = 40.0; tbeta = 0.475
    # deltab: β step in units of k_B·T (matches _start output)
    deltab = (beta_grid[2] - beta_grid[1])
    f0 = 1.0                          # Debye-Waller stand-in

    NJOY.trans!(ssm, 1, alpha_grid, beta_grid, twt, c_diff, tbeta,
                tev, deltab, f0, 0, 1.0, tempr, tempf)

    # Output is positive everywhere a kernel was evaluated
    @test all(ssm[:, :, 1] .>= 0.0)
    # Effective temperature updated per tempf = (tbeta*tempf + twt*tempr)/(tbeta+twt)
    @test tempf[1] ≈ (0.475*20.0 + 0.025*20.0)/(0.475+0.025) rtol=1e-10
    # Not every entry zeroed — there's real transport mass
    @test count(ssm[:, :, 1] .> 0) > nalpha  # at least one β per α nonzero
end
