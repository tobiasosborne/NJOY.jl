using Test
using NJOY

# T13 regression for NJOY_jl-xpf — GASPR MT203-207 survive final_assembly!.
#
# T13's deck is moder -> reconr -> broadr -> heatr -> gaspr -> moder(-25 28).
# heatr populates ctx.extra_mf3, which makes `needs_assembly` true, so
# `final_assembly!` rebuilds the last moder output from the RunContext via
# write_full_pendf instead of leaving moder's copy of gaspr's tape in place.
# GASPR was never collected into the context (`_collect_gaspr!` existed at
# pipeline.jl:471 with zero call sites), so the rebuild silently dropped the
# gas sections that gaspr had correctly written.
#
# Measured before the fix: gaspr's own tape25 carried MF3 MT203 and MT207, but
# tape28 came back with only the 23 sections heatr had contributed.
#
# The symmetric hazard: collecting gas sections must NOT by itself make
# `needs_assembly` true, or a plain moder->reconr->broadr->gaspr->moder chain
# (T45) would start reconstructing a tape that was already exact.

const T13C_DIR = joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "13")
const T13C_WORK_DIR = "/tmp/njoy_t13_gaspr_collection"

"MF3 MT numbers present for `mat` in an ENDF/PENDF tape."
function t13c_mf3_mts(path::AbstractString)
    mts = Set{Int}()
    for line in eachline(path)
        p = rpad(line, 80)
        NJOY._parse_int(p[67:70]) > 0 || continue
        NJOY._parse_int(p[71:72]) == 3 || continue
        mt = NJOY._parse_int(p[73:75])
        mt > 0 && push!(mts, mt)
    end
    mts
end

@testset "T13 final tape keeps GASPR's MT203/MT207 through final_assembly!" begin
    isdir(T13C_WORK_DIR) && rm(T13C_WORK_DIR; force=true, recursive=true)
    mkpath(T13C_WORK_DIR)
    run_njoy(joinpath(T13C_DIR, "input"); work_dir=T13C_WORK_DIR)

    gaspr_out = joinpath(T13C_WORK_DIR, "tape25")   # gaspr's own output
    final_out = joinpath(T13C_WORK_DIR, "tape28")   # after moder + assembly
    @test isfile(gaspr_out)
    @test isfile(final_out)

    # Ni-61 produces protons and alphas but no deuterons/tritons/He-3, so
    # gaspr emits exactly MT203 and MT207 here.
    gas = t13c_mf3_mts(gaspr_out)
    @test 203 in gas
    @test 207 in gas

    # The whole point of the bead: they must still be there at the end.
    final = t13c_mf3_mts(final_out)
    @test 203 in final
    @test 207 in final

    # And the reference agrees they belong there.
    ref = t13c_mf3_mts(joinpath(T13C_DIR, "referenceTape28"))
    @test 203 in ref
    @test 207 in ref

    # Assembly must not lose anything else gaspr passed through either.
    @test issubset(gas, final)
end
