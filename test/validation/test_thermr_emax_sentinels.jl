# RED→GREEN unit test for the emax-cutoff sentinel sequence on the shared
# thermr MT221/MT222 energy grid (Fortran tpend / sigcoh / save-elastic).
#
# Bug: T70 tape60 MF3/MT221 and MF3/MT222 emit NP=1484, but the reference
# (njoy-reference/tests/70/referenceTape60) has NP=1485 — one fewer energy
# point. The missing point is sigfig(emax,7,-1) (=2.276999 for emax=2.277),
# the shaded-DOWN held duplicate at the thermal cutoff.
#
# The reference grid near emax is the 4-point sequence:
#   sigfig(emax,7,-1) [held thermal value]  — coherent grid's last point
#   emax              [held value]
#   sigfig(emax,7,+1) [σ=0 above cutoff]
#   etop=2e7          [σ=0]
#
# Fortran ground truth:
#   thermr.f90:1194 (sigcoh)        enext = sigfig(elim,7,-1)  -> last coh point
#   thermr.f90:393  (save-elastic)  enext = emax               -> held at emax
#   thermr.f90:394  (save-elastic)  enext = up*emax            -> sigfig(emax,7,+1)
#   thermr.f90:3211-3213 (tpend)    ex(1) = etop  (=2e7) at ib = ne+1
#
# Run:  julia --project=. test/validation/test_thermr_emax_sentinels.jl

using Test
using NJOY

@testset "thermr emax sentinel sequence (T70 MT221/MT222 NP)" begin
    emax = 2.277
    em1 = NJOY.round_sigfig(emax, 7, -1)   # 2.276999 — held thermal value
    ep1 = NJOY.round_sigfig(emax, 7, 1)    # 2.277001 — step to zero
    etop = 2e7

    # Simulate the shared thermal grid as produced upstream WITHOUT the
    # cutoff sentinels: lower-energy points only, ending below emax. The
    # helper is responsible for injecting the full 4-point cutoff sequence
    # (sigcoh's sigfig(emax,7,-1) is NOT already present here).
    grid = Float64[1e-5, 1.0, 2.0]
    NJOY._append_emax_sentinels!(grid, emax)

    # ---- THE FIX: shaded-down held point must survive (FAILS before fix) ----
    @test em1 in grid

    # ---- other boundary points present ----
    @test emax in grid
    @test ep1 in grid
    @test etop in grid

    # ---- grid is sorted ----
    @test issorted(grid)

    # ---- the four boundary points appear in order -1, 0, +1, etop ----
    i_em1 = findfirst(==(em1), grid)
    i_em  = findfirst(==(emax), grid)
    i_ep1 = findfirst(==(ep1), grid)
    i_top = findfirst(==(etop), grid)
    @test i_em1 < i_em < i_ep1 < i_top
end
