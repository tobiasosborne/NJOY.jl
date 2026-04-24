# Structural tests for the purr MT152/MT153 writer.
#
# These tests exercise _write_mt152_purr and _write_mt153_purr in
# isolation: minimal URRForUnresr fixture + synthetic PurrBlock data,
# check record layout, line counts, DICT NC formulas against the
# Fortran purr.f90 reference formulas.

using Test
using Printf
using NJOY

@testset "purr DICT NC counts (Fortran formulas)" begin
    # Fortran purr.f90:365-368 canonical formulas.
    @test NJOY._purr_mt152_nc(1, 13) == 2 + div(1 + 13 * 6 - 1, 6)
    @test NJOY._purr_mt152_nc(1, 13) == 15          # matches T21 DICT
    @test NJOY._purr_mt152_nc(1, 36) == 38          # matches T38 DICT
    @test NJOY._purr_mt153_nc(20, 13) == 264        # matches T21 DICT
    @test NJOY._purr_mt153_nc(20, 36) == 727        # matches T38 DICT
end

# Minimal URR stand-in sharing URRForUnresr's interface expected by the writer.
@testset "MT152 record layout" begin
    urr = (za=2.6058e4, awr=57.4356, lssf=1, intunr=5, energies=[3.5e5, 4e5, 5e5])
    sigz = [1e10]
    block = NJOY.PurrBlock(
        293.6,
        reshape(Float64[4.567886, 4.556627, 0.0, 0.01126, 4.567886,
                         4.382575, 4.372347, 0.0, 0.01023, 4.382575,
                         4.117175, 4.108556, 0.0, 0.00862, 4.117175],
                 3, 5),
        zeros(20, 5, 3),
        zeros(20, 3),
    )

    buf = IOBuffer()
    NJOY._write_mt152_purr(buf, 2637, urr, sigz, block)
    lines = split(String(take!(buf)), '\n'; keepempty=false)

    # 2 header lines + ceil((1 + 3*(1+5))/6)=4 data lines + 1 SEND = 7
    @test length(lines) == 7

    # CONT header: cols 71-73 should read "2", col 74 should be "152" prefix
    @test occursin("2637 2152", lines[1])
    # LIST header should have NPL=19 and NE=3
    @test occursin(lpad("19", 11), lines[2])
    @test occursin(lpad("3", 11) * "2637 2152", lines[2])
    # First data line starts with "1.00000+10" (sigz[1])
    @test startswith(lines[3], " 1.00000+10")
    # Last line is SEND (MT=0)
    @test occursin("2637 2  0", lines[end])
end

@testset "MT153 record layout" begin
    urr = (za=3.6083e4, awr=82.202, lssf=0, intunr=5, energies=[2.72e2, 3e2, 3.5e2])
    nbin = 20
    block = NJOY.PurrBlock(293.6, zeros(3, 5), zeros(nbin, 5, 3), zeros(nbin, 3))

    buf = IOBuffer()
    NJOY._write_mt153_purr(buf, 3640, urr, nbin, -1, -1, block)
    lines = split(String(take!(buf)), '\n'; keepempty=false)

    # Per energy: 1 + 6*nbin = 121. Three energies => 363. Lines = 2 (hdr) +
    # ceil(363/6) = 61 + 1 (SEND) = 64
    @test length(lines) == 2 + div(121 * 3 + 5, 6) + 1
    @test occursin("3640 2153", lines[1])
    # LIST header NPL = (1+6*20)*3 = 363
    @test occursin(lpad("363", 11), lines[2])
    @test occursin("3640 2  0", lines[end])
end
