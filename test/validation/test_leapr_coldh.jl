using Test
using NJOY

# Fortran ref: njoy-reference/src/leapr.f90:1798-2466 (bfill, exts, sint,
# coldh, bt, sumh, cn, sjbes, terpk). Julia port: src/processing/leapr.jl.

@testset "_bt statistical weight (partition normalization)" begin
    # At x → 0: p_j = (2j+1) / [2·Σ_k (2k'+1)] where k' has same parity as j
    # and ranges over 10 values (k ∈ 0..18 stepped by 2 from matching start).
    for j in (0, 1, 2, 3, 4)
        pj = NJOY._bt(j, 0.0)
        # Σ over 10 parity-matched k ranges: k' = parity_of_j + 0, 2, ..., 18
        start = isodd(j) ? 1 : 0
        s = sum((2*(start + 2*i) + 1) for i in 0:9)
        @test pj ≈ (2j + 1) / (2 * s) rtol=1e-12
    end
    # All weights positive and sum bounded
    @test all(NJOY._bt(j, 0.1) > 0 for j in 0:10)
end

@testset "_cn Clebsch-Gordon edge cases" begin
    # (jj+ll+nn) odd → 0
    @test NJOY._cn(1, 1, 1) == 0.0
    @test NJOY._cn(0, 1, 0) == 0.0
    @test NJOY._cn(2, 1, 0) == 0.0
    # (0,0,0) symmetric → sqrt(1) = 1
    @test NJOY._cn(0, 0, 0) ≈ 1.0 rtol=1e-12
    # Stable for moderate j
    for (jj, ll, nn) in [(1,1,0), (1,1,2), (2,2,0), (2,2,2), (2,2,4), (3,3,0)]
        @test isfinite(NJOY._cn(jj, ll, nn))
    end
end

@testset "_sjbes vs analytic spherical Bessel" begin
    # j_0(x) = sin(x)/x, j_1(x) = sin(x)/x² - cos(x)/x
    for x in (0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0)
        j0 = NJOY._sjbes(0, x)
        j1 = NJOY._sjbes(1, x)
        j0_exact = x < 1e-3 ? 1 - x^2/6 : sin(x)/x
        j1_exact = x < 1e-3 ? x/3 - x^3/30 : sin(x)/x^2 - cos(x)/x
        @test j0 ≈ j0_exact rtol=1e-5
        @test j1 ≈ j1_exact rtol=1e-4
    end
    # Higher-order: j_2(x) = (3/x² - 1)·sin(x)/x - 3·cos(x)/x²
    for x in (0.5, 1.0, 2.0)
        j2 = NJOY._sjbes(2, x)
        j2_exact = (3/x^2 - 1) * sin(x)/x - 3*cos(x)/x^2
        @test j2 ≈ j2_exact rtol=1e-4
    end
end

@testset "_sumh boundary cases" begin
    # j=0, jp=0: (sjbes(0,y)·cn(0,0,0))² = (sin(y)/y)²
    for y in (0.1, 0.5, 1.0, 2.0)
        @test NJOY._sumh(0, 0, y) ≈ (sin(y)/y)^2 rtol=1e-5
    end
    # All positive (squared sums)
    for j in 0:3, jp in 0:3
        for y in (0.1, 1.0, 3.0)
            @test NJOY._sumh(j, jp, y) >= 0
        end
    end
end

@testset "_terpk linear interpolation (Fortran-faithful off-by-one)" begin
    # Fortran's terpk (leapr.f90:2444) does `i=int(be/delta); i=i+1; terpk=ska(i)+...`
    # so at be = k·dka it returns ska[k+1], not ska[k]. Our port preserves this.
    nka = 5; dka = 0.1
    ska = [1.0, 2.0, 4.0, 3.0, 2.0]
    @test NJOY._terpk(ska, nka, dka, 0.1) ≈ 2.0 rtol=1e-12   # → ska[2]
    @test NJOY._terpk(ska, nka, dka, 0.2) ≈ 4.0 rtol=1e-12   # → ska[3]
    @test NJOY._terpk(ska, nka, dka, 0.15) ≈ 3.0 rtol=1e-12  # halfway: (2+4)/2
    @test NJOY._terpk(ska, nka, dka, 0.25) ≈ 3.5 rtol=1e-12  # halfway: (4+3)/2
    @test NJOY._terpk(ska, nka, dka, 10.0) == 1.0            # past table → 1
end

@testset "_bfill! monotone extended β grid" begin
    nbeta = 5; maxbb = 2*nbeta + 1
    betan = [0.0, 1.0, 2.0, 4.0, 8.0]
    bex   = zeros(maxbb)
    rdbex = zeros(maxbb)
    nbx = NJOY._bfill!(bex, rdbex, betan, nbeta, maxbb)
    @test nbx >= 9    # (nbeta-1 negatives) + 0 + (nbeta-1 positives) = 9
    @test all(bex[i+1] > bex[i] for i in 1:nbx-1)
    @test 0.0 in bex[1:nbx]
    # rdbex reciprocal spacing
    for i in 1:nbx-1
        @test rdbex[i] ≈ 1/(bex[i+1] - bex[i]) rtol=1e-12
    end
end

@testset "_exts! detailed balance" begin
    nbeta = 4; maxbb = 2*nbeta + 1
    betan = [0.0, 1.0, 2.0, 3.0]
    sexpb = [1.0, 0.5, 0.2, 0.1]        # S(α,-β) values (decreasing in β)
    exb   = [exp(-betan[i]/2) for i in 1:nbeta]
    sex   = zeros(maxbb)
    NJOY._exts!(sexpb, sex, exb, betan, nbeta, maxbb)
    # For β[1]=0, no duplicate → sex[nbeta-1+1] = sexpb[1]
    # For i>1, sex maps positive side with detailed balance: exp(-β)·sexpb
    # Specifically sex[nbeta+1] = sexpb[2]·exp(-β₂), sex[nbeta+2] = sexpb[3]·exp(-β₃), ...
    @test sex[nbeta+1] ≈ sexpb[2] * exb[2]^2 rtol=1e-12   # = sexpb[2]·exp(-β₂)
    @test sex[nbeta+2] ≈ sexpb[3] * exb[3]^2 rtol=1e-12
    @test sex[nbeta+3] ≈ sexpb[4] * exb[4]^2 rtol=1e-12
end

@testset "_sint interpolation + SCT fallback" begin
    nbeta = 5; maxbb = 2*nbeta + 1
    betan = [0.0, 1.0, 2.0, 3.0, 4.0]
    exb   = [exp(-betan[i]/2) for i in 1:nbeta]
    sexpb = [1.0, 0.5, 0.3, 0.2, 0.1]
    bex   = zeros(maxbb); rdbex = zeros(maxbb); sex = zeros(maxbb)
    nbx = NJOY._bfill!(bex, rdbex, betan, nbeta, maxbb)
    NJOY._exts!(sexpb, sex, exb, betan, nbeta, maxbb)
    # At β=0 (grid point) → sexpb[1] (since first β is zero, no duplicate)
    s0 = NJOY._sint(0.0, bex, rdbex, sex, nbx, 1.0, 1.0, 1.0, betan, nbeta, maxbb)
    @test s0 > 0
    # SCT fallback beyond table: uses Gaussian form, positive
    s_sct = NJOY._sint(10.0, bex, rdbex, sex, nbx, 1.0, 0.5, 1.0, betan, nbeta, maxbb)
    @test s_sct >= 0
end

@testset "coldh! smoke test (para-H2, ncold=2)" begin
    # Minimal synthetic input — verify ssm/ssp populated with finite,
    # non-negative values; no NaN/Inf; energy-grid basics.
    nalpha = 4; nbeta = 8
    alpha_grid = collect(range(0.1, 2.0, length=nalpha))
    beta_grid  = collect(range(0.0, 10.0, length=nbeta))
    nka = 20; dka = 0.05
    ska = [1.0 - 0.3*exp(-(i*dka - 0.5)^2) for i in 1:nka]   # bumpy positive

    ssm = zeros(nbeta, nalpha, 1)
    ssp = zeros(nbeta, nalpha, 1)
    # Seed ssm with a smooth S(α,-β) approximation
    for j in 1:nalpha, i in 1:nbeta
        ssm[i, j, 1] = exp(-alpha_grid[j] - 0.5*beta_grid[i])
    end

    temp = 20.0
    tev  = NJOY.PhysicsConstants.bk * temp
    twt = 0.025; tbeta = 0.475
    tempr = [20.0]; tempf = [20.0]

    NJOY.coldh!(ssm, ssp, 1, alpha_grid, beta_grid, ska, nka, dka,
                2, 0, temp, tev, twt, tbeta, 0, 1.0, tempr, tempf)

    @test all(isfinite.(ssm[:, :, 1]))
    @test all(isfinite.(ssp[:, :, 1]))
    @test all(ssm[:, :, 1] .>= 0)
    @test all(ssp[:, :, 1] .>= 0)
    # Should see real transport mass in both directions
    @test count(ssm[:, :, 1] .> 0) > nalpha
    @test count(ssp[:, :, 1] .> 0) > nalpha
end
