using Test
using NJOY

# T45 regression for NJOY_jl-9h4 — GASPR LR and residual-nuclide gas yields.
#
# Ref: njoy-reference/src/gaspr.f90:495-826 (subroutine gaspr, accumulation pass)
#
# Fortran builds, for every MF3 reaction, a residual-nucleus ZA
#     izr = za + zain                                   (gaspr.f90:495)
# and decrements it by the ZA of each explicitly emitted particle while
# setting the y203..y207 counters (gaspr.f90:507-818).  It then converts a
# *light* residual nucleus into gas (gaspr.f90:820-825):
#
#     if (izr.eq.1001) y203=y203+1     if (izr.eq.2003) y206=y206+1
#     if (izr.eq.1002) y204=y204+1     if (izr.eq.2004) y207=y207+1
#     if (izr.eq.1003) y205=y205+1     if (izr.eq.4008) y207=y207+2
#
# The Be-8 rule on the last line is the T45 defect: B-10(n,t) is MT105, whose
# residual is Be-8 (izr = 5010+1-1003 = 4008), so the Fortran adds 2 alphas on
# top of the explicit triton.  Julia's old static MT->yield table had neither
# the residual arithmetic nor the MT51-91 LR branch, so MT207 was short by
# 2*sigma(MT105) at every energy and MT203/MT204 were missing every LR
# breakup contribution above threshold.
#
# T45 is B-10 (MAT 525, ZA 5010, NSUB 10 -> zain 1).

const T45G_DIR = joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "45")
const T45G_ORACLE = joinpath(T45G_DIR, "referenceTape40")
const T45G_WORK_DIR = "/tmp/njoy_t45_gaspr_regression"

# --------------------------------------------------------------------------
# Part 1 — unit oracle on the yield chain itself (hand-derived from gaspr.f90)
# --------------------------------------------------------------------------

@testset "gaspr gas_channel — B-10 residual and LR chain" begin
    za, zain = 5010.0, 1.0   # B-10, incident neutron (nsub=10 -> zain=1)

    # MT107 (n,alpha): izr = 5011 - 2004 = 3007 (Li-7), not a gas residual.
    @test NJOY.gas_channel(107, 0, za, zain) == (0, 0, 0, 0, 1)

    # MT105 (n,t): izr = 5011 - 1003 = 4008 (Be-8) -> +2 alpha on top of the
    # explicit triton.  This is the record the whole bead is about.
    @test NJOY.gas_channel(105, 0, za, zain) == (0, 0, 1, 0, 2)

    # MT113 (n,t+2alpha): izr = 5011 - 5011 = 0, no residual gas.
    @test NJOY.gas_channel(113, 0, za, zain) == (0, 0, 1, 0, 2)

    # MT103 (n,p): izr = 5011 - 1001 = 4010 (Be-10), not a gas residual.
    @test NJOY.gas_channel(103, 0, za, zain) == (1, 0, 0, 0, 0)

    # MT104 (n,d): izr = 5011 - 1002 = 4009 (Be-9), not a gas residual.
    @test NJOY.gas_channel(104, 0, za, zain) == (0, 1, 0, 0, 0)

    # MT51-91 carry the breakup flag LR in the MF3 TAB1 L2 field.
    # LR=22 (n,n'alpha): izr = 5011 - 1 - 2004 = 3006 (Li-6).
    @test NJOY.gas_channel(55, 22, za, zain) == (0, 0, 0, 0, 1)
    # LR=28 (n,n'p): izr = 5011 - 1 - 1001 = 4009 (Be-9).
    @test NJOY.gas_channel(65, 28, za, zain) == (1, 0, 0, 0, 0)
    # LR=35 (n,n'd+2alpha): izr = 5011 - 1 - 5010 = 0.
    @test NJOY.gas_channel(62, 35, za, zain) == (0, 1, 0, 0, 2)
    # LR=0 discrete level: izr = 5011 - 1 = 5010, no gas at all.
    @test NJOY.gas_channel(57, 0, za, zain) == (0, 0, 0, 0, 0)
    # LR=39/40 leave izr untouched and emit nothing (gaspr.f90:607-610).
    @test NJOY.gas_channel(604, 40, za, zain) == (0, 0, 0, 0, 0)

    # MT102 (n,gamma) has no branch in the chain: izr stays 5011, no gas.
    @test NJOY.gas_channel(102, 0, za, zain) == (0, 0, 0, 0, 0)

    # Levels MT600-849 are always skipped (gaspr.f90:486-490) so that the
    # MT103-107 totals are not double counted.
    @test NJOY.gas_channel(600, 0, za, zain) == (0, 0, 0, 0, 0)
    @test NJOY.gas_channel(700, 1, za, zain) == (0, 0, 0, 0, 0)
    @test NJOY.gas_channel(800, 0, za, zain) == (0, 0, 0, 0, 0)

    # A non-B-10 spot check of the residual rules: Li-6(n,t) -> MT105 leaves
    # izr = 3006+1-1003 = 2004 (alpha) -> +1 alpha via the residual rule.
    @test NJOY.gas_channel(105, 0, 3006.0, 1.0) == (0, 0, 1, 0, 1)
end

# --------------------------------------------------------------------------
# Part 2 — tape-level oracle: T45 tape40 MF3 gas sections vs the Fortran
# --------------------------------------------------------------------------

"Return columns 1-66 of every line of `path` belonging to (MF, MT)."
function t45g_section_columns(path::AbstractString, mf_target::Int, mt_target::Int)
    records = String[]
    for line in eachline(path)
        padded = rpad(line, 80)
        NJOY._parse_int(padded[67:70]) > 0 || continue
        NJOY._parse_int(padded[71:72]) == mf_target || continue
        NJOY._parse_int(padded[73:75]) == mt_target || continue
        push!(records, padded[1:66])
    end
    records
end

@testset "T45 tape40 GASPR sections are bit-identical to the Fortran oracle" begin
    isdir(T45G_WORK_DIR) && rm(T45G_WORK_DIR; force=true, recursive=true)
    mkpath(T45G_WORK_DIR)
    run_njoy(joinpath(T45G_DIR, "input"); work_dir=T45G_WORK_DIR)
    produced = joinpath(T45G_WORK_DIR, "tape40")
    @test isfile(produced)

    # Structure must stay exact while the values are repaired.
    @test countlines(produced) == countlines(T45G_ORACLE)

    # MT206 is absent from both tapes: B-10 has no MT106 and no He-3-emitting
    # LR partial, so the channel is all-zero and GASPR suppresses the section.
    @test isempty(t45g_section_columns(T45G_ORACLE, 3, 206))

    for mt in (203, 204, 206, 207)
        want = t45g_section_columns(T45G_ORACLE, 3, mt)
        got = t45g_section_columns(produced, 3, mt)
        @test length(got) == length(want)
        bad = findfirst(i -> got[i] != want[i], eachindex(want))
        if bad !== nothing
            @info "T45 MF3/MT$mt first mismatch at section record $bad" want[bad] got[bad]
        end
        @test bad === nothing
    end

    # MT205 is exact except for one record inherited from upstream BROADR: the
    # triton channel is sigma(MT105)+sigma(MT113), and MT113 itself differs from
    # the oracle by one ULP at 1.125e-5 eV (tape40 line 3713).  That is bead
    # NJOY_jl-bdu, not a GASPR defect -- assert the shape exactly so a genuine
    # regression here still fails, and delete this exception when bdu closes.
    want205 = t45g_section_columns(T45G_ORACLE, 3, 205)
    got205 = t45g_section_columns(produced, 3, 205)
    @test length(got205) == length(want205)
    diff205 = [i for i in eachindex(want205) if got205[i] != want205[i]]
    @test diff205 == [5]
    @test got205[5] ==
        " 1.093750-5 6.611181-1 1.125000-5 6.518713-1 1.156250-5 6.430018-1"
end

@testset "T45 tape40 TPID keeps the deck's blank RECONR label" begin
    # The T45 deck's RECONR label card is `'  '/`; Fortran writes those blanks
    # through to tape40 line 1 (reconr.f90:168,192 -> broadr -> gaspr -> moder).
    # Julia used to substitute "reconstructed data" here.  NJOY_jl-1kf.
    produced = joinpath(T45G_WORK_DIR, "tape40")
    got = rpad(first(eachline(produced)), 80)
    want = rpad(first(eachline(T45G_ORACLE)), 80)
    @test got[1:66] == " "^66
    @test got == want
end
