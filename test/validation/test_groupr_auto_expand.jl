# RED→GREEN test for groupr `3 /` auto-expand sentinel.
#
# T15/T16/T17 groupr decks use `3 /` (mfd=3, mtd absent) as a Fortran
# sentinel meaning "process every MF=3 MT on the PENDF passing the filter
# mt<=200 OR 203<=mt<=207 OR mt>300". See groupr.f90:612-634 + nextr at
# 1087-1123.
#
# Pre-fix: Julia parse_groupr reads missing mtd as 0; groupr_module drops
# the entry. T15 tape91 gets 6 MTs {251, 252, 452, 455, 456, MF=5/MT=18}
# instead of the reference's 40.
#
# Post-fix: parser emits a sentinel, groupr_module expands to the full
# auto-set from the PENDF tape.
#
# Fast loop: uses the cached post-broadr PENDF at
# test/validation/oracle_cache/test15/after_broadr.pendf — skips the
# ~500s U-238 broadening.
#
# Run:  julia --project=. test/validation/test_groupr_auto_expand.jl

using Test
using NJOY

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const T15_CACHE = joinpath(REPO, "test", "validation", "oracle_cache", "test15")
const T15_INPUT = joinpath(REPO, "njoy-reference", "tests", "15", "input")

# Auto-expanded MTs present in the Fortran reference tape91 (9237).
# Derived from referenceTape91 MF=3 scan:
#   awk 'length($0)>=75 && substr($0,67,4)+0==9237 &&
#        substr($0,71,2)+0==3 && substr($0,73,3)+0>0'
# minus explicit {251, 252, 452, 455, 456}.
const AUTO_EXPANDED_MTS = sort!(collect(
    Int[1, 2, 4, 16, 17, 18, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
        61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75,
        76, 77, 91, 102]
))

const EXPLICIT_MTS = sort!(Int[251, 252, 452, 455, 456])

"""Scan an ENDF/GENDF tape, return sorted unique (mf, mt) pairs for `mat`."""
function _tape_mts(path::AbstractString, mat::Int; mf::Int = 3)
    mts = Set{Int}()
    for line in eachline(path)
        length(line) < 75 && continue
        p = rpad(line, 80)
        m  = tryparse(Int, strip(p[67:70])); m === nothing && continue
        f  = tryparse(Int, strip(p[71:72])); f === nothing && continue
        t  = tryparse(Int, strip(p[73:75])); t === nothing && continue
        m == mat && f == mf && t > 0 && push!(mts, t)
    end
    return sort!(collect(mts))
end

"""Run Julia groupr for T15 against cached broadr PENDF. Returns output path."""
function _run_t15_groupr(; work_dir::AbstractString = mktempdir())
    # Parse T15 input deck, pluck the first groupr call.
    txt = read(T15_INPUT, String)
    calls = NJOY.tokenise_njoy_input(txt)
    gcall = calls[findfirst(c -> c.name === :groupr, calls)]
    params = NJOY.parse_groupr(gcall)

    # Stage tapes: unit 21 → real ENDF, unit 23 → cached broadr PENDF.
    tapes = NJOY.TapeManager(Dict{Int,String}(), work_dir)
    NJOY.register!(tapes, 21,
        joinpath(T15_CACHE, "run_broadr", "tape20"))        # ENDF source
    NJOY.register!(tapes, 23,
        joinpath(T15_CACHE, "after_broadr.pendf"))           # broadr PENDF

    NJOY.groupr_module(tapes, params)
    return NJOY.resolve(tapes, abs(params.nout))
end

@testset "groupr `3 /` auto-expand (T15)" begin
    out = _run_t15_groupr()
    @test isfile(out)

    mts = _tape_mts(out, 9237; mf = 3)
    @info "T15 Julia groupr MF=3 MTs ($(length(mts))): $mts"

    # Every auto-expanded MT must be present (the sentinel bug — fails pre-fix).
    for mt in AUTO_EXPANDED_MTS
        @test mt in mts
    end

    # Nubar MTs already worked pre-fix; guard against regression.
    for mt in (452, 455, 456)
        @test mt in mts
    end

    # Known out-of-scope gaps, filed as separate beads:
    #   NJOY.jl-cdy: MT=251/252 need MF=4/MF=6 derivation (not MF=3 lookup).
    #   NJOY.jl-5oi: MT=37 emitted with all-zero groups (threshold 17.82 MeV
    #                above LANL-30 top 17.0 MeV). Fortran skips such MTs.
    @test_broken 251 in mts
    @test_broken 252 in mts
    @test_broken !(37 in mts)
end
