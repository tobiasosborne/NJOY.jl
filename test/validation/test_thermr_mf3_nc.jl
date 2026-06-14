# RED→GREEN unit test for the MF1/MT451 directory NC (line-count) formula
# used by the PENDF writer for thermr-augmented MF3 sections.
#
# Bug NJOY_jl-h61: T70 PENDF MF1/MT451 directory wrote NC=497 for MF3/MT221
# and MF3/MT222 (NP=1485), but the reference (njoy-reference/tests/70/
# referenceTape60 lines 69-70) says NC=498 and the emitted sections ARE 498
# lines. The writer substituted a stale, smaller `thermr_coh_ne` for the
# actual emitted np. The faithful fix mirrors Fortran thermr tpend.
#
# Fortran ground truth:
#   thermr.f90:3087  nc = 3 + (ne+2)/3   (integer division)
#   thermr.f90:3198  scr(6) = ne+1       (written NP includes the etop sentinel)
#   thermr.f90:3211  ex(1) = etop        (the appended sentinel point)
# So in terms of emitted np = ne+1:  nc = 3 + div((np-1)+2, 3) = 3 + div(np+1, 3).
#
# This DIVERGES from the reconr/ceil form  3 + cld(np,3) = 3 + div(np+2,3)
# exactly when np ≡ 1 (mod 3): there Fortran undercounts the true record
# count by 1. We reproduce that quirk (GROUND-TRUTH PRINCIPLE), not "fix" it.
# reconr tpend (reconr.f90:5115-5116) has no such sentinel and uses cld(np,3).
#
# Run:  julia --project=. test/validation/test_thermr_mf3_nc.jl

using Test
using NJOY

@testset "thermr/reconr MF3 directory NC formulas (NJOY_jl-h61)" begin
    # ---- Reference values from njoy-reference reference tapes ----
    # T70 MF3/MT221 & MT222: NP=1485 → NC=498  (referenceTape60 lines 69-70)
    @test NJOY._thermr_mf3_dir_nc(1485) == 498
    # T01 MF3/MT221 free-gas: NP=146 → NC=52
    @test NJOY._thermr_mf3_dir_nc(146) == 52
    # T68 MF3/MT221 coherent (nbin=20): NP=378 → NC=129
    @test NJOY._thermr_mf3_dir_nc(378) == 129

    # ---- np ≡ 1 (mod 3): Fortran undercount, NOT 499 ----
    # 1486 mod 3 == 1.  thermr (sentinel) form gives 3 + div(1487,3) = 3+495 = 498.
    @test NJOY._thermr_mf3_dir_nc(1486) == 498
    # Contrast: reconr/ceil form gives 3 + cld(1486,3) = 3 + div(1488,3) = 3+496 = 499.
    @test NJOY._reconr_mf3_dir_nc(1486) == 499

    # ---- The two formulas agree for np ≢ 1 (mod 3) ----
    for np in (1485, 146, 378, 1487, 1488)  # ≡ 0 or 2 mod 3
        @test NJOY._thermr_mf3_dir_nc(np) == NJOY._reconr_mf3_dir_nc(np)
    end
end
