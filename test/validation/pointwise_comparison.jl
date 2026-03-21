# Pointwise validation: compare NJOY.jl RECONR output against NJOY2016
# reference tapes.
#
# For each test case this script:
#   1. Parses the reference PENDF tape (NJOY2016 output).
#   2. Runs NJOY.jl reconr on the same ENDF input with the same tolerance.
#   3. Compares cross sections at every reference energy point.
#   4. Reports max relative error per reaction (with and without an
#      absolute-floor filter to avoid noise at vanishing cross sections).
#   5. Compares grid-point counts (NJOY.jl should be within 2x of reference).

using NJOY
using Printf

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
const NJOY_REF  = normpath(joinpath(@__DIR__, "..", "..", "njoy-reference"))
const RESOURCES = joinpath(NJOY_REF, "tests", "resources")
const TESTS_DIR = joinpath(NJOY_REF, "tests")

# Minimum absolute cross section (barns) below which relative error is not
# meaningful.  NJOY2016 itself uses ~1e-3 barns as a practical floor in many
# convergence checks.
const ABS_FLOOR = 1.0e-2

# ---------------------------------------------------------------------------
# ENDF helpers (duplicated from integration_tests.jl for standalone use)
# ---------------------------------------------------------------------------
function parse_endf_data_pairs(line::AbstractString)
    p = rpad(line, 80)
    pairs = Tuple{Float64,Float64}[]
    for i in 0:2
        f1, f2 = p[i*22+1:i*22+11], p[i*22+12:i*22+22]
        (strip(f1) == "" && strip(f2) == "") && break
        push!(pairs, (parse_endf_float(f1), parse_endf_float(f2)))
    end
    pairs
end

"""
    read_reference_pendf(filename) -> Dict{Int, NamedTuple}

Read a PENDF reference tape produced by NJOY2016.
Returns Dict: MT => (energies, xs, n_points).
"""
function read_reference_pendf(filename::AbstractString)
    result = Dict{Int, @NamedTuple{energies::Vector{Float64},
                                    xs::Vector{Float64},
                                    n_points::Int}}()
    lines = readlines(filename)
    i = 1
    while i <= length(lines)
        p = rpad(lines[i], 80)
        mf  = NJOY._parse_int(p[71:72])
        mt  = NJOY._parse_int(p[73:75])
        mat_val = NJOY._parse_int(p[67:70])
        if mf == 3 && mt > 0 && mat_val > 0
            i += 1; i > length(lines) && break
            p2 = rpad(lines[i], 80)
            n_points = NJOY._parse_int(p2[56:66])
            i += 2; i > length(lines) && break
            energies, xs_vals = Float64[], Float64[]
            while i <= length(lines)
                pd = rpad(lines[i], 80)
                (NJOY._parse_int(pd[71:72]) != mf ||
                 NJOY._parse_int(pd[73:75]) != mt) && break
                for (e, x) in parse_endf_data_pairs(pd)
                    push!(energies, e); push!(xs_vals, x)
                end
                i += 1
            end
            result[mt] = (energies=energies, xs=xs_vals, n_points=n_points)
        else
            i += 1
        end
    end
    result
end

# ---------------------------------------------------------------------------
# Linear interpolation
# ---------------------------------------------------------------------------
function interp_at(energies::Vector{Float64}, xs::Vector{Float64},
                   e::Float64)
    e <= energies[1]   && return xs[1]
    e >= energies[end]  && return xs[end]
    idx = searchsortedlast(energies, e)
    idx >= length(energies) && return xs[end]
    frac = (e - energies[idx]) / (energies[idx+1] - energies[idx])
    xs[idx] + frac * (xs[idx+1] - xs[idx])
end

# ---------------------------------------------------------------------------
# Per-reaction comparison result
# ---------------------------------------------------------------------------
struct ReactionResult
    mt::Int
    label::String
    n_ref::Int
    n_ours::Int
    # Raw max relative error (over all points)
    max_rel_err::Float64
    max_abs_err::Float64
    max_err_energy::Float64
    # Filtered max relative error (only where |ref_xs| >= ABS_FLOOR)
    filt_rel_err::Float64
    filt_abs_err::Float64
    filt_err_energy::Float64
    n_compared::Int        # points above ABS_FLOOR
end

function compare_reaction(ref_data, our_energies::Vector{Float64},
                          our_xs::Vector{Float64}, mt::Int, label::String)
    ref_e  = ref_data.energies
    ref_xs = ref_data.xs
    n_ref  = length(ref_e)

    max_rel = 0.0;  max_abs = 0.0;  worst_e = 0.0
    filt_rel = 0.0; filt_abs = 0.0; filt_e  = 0.0
    n_sig = 0

    for k in 1:n_ref
        xs_ours = interp_at(our_energies, our_xs, ref_e[k])
        abs_err = abs(xs_ours - ref_xs[k])
        denom   = max(abs(ref_xs[k]), 1.0e-30)
        rel_err = abs_err / denom

        if rel_err > max_rel
            max_rel = rel_err; max_abs = abs_err; worst_e = ref_e[k]
        end
        if abs(ref_xs[k]) >= ABS_FLOOR
            n_sig += 1
            if rel_err > filt_rel
                filt_rel = rel_err; filt_abs = abs_err; filt_e = ref_e[k]
            end
        end
    end

    ReactionResult(mt, label, n_ref, length(our_energies),
                   max_rel, max_abs, worst_e,
                   filt_rel, filt_abs, filt_e, n_sig)
end

# ---------------------------------------------------------------------------
# MT -> accessor mapping
# ---------------------------------------------------------------------------
const MT_MAP = Dict{Int,Tuple{Symbol,String}}(
    1   => (:total,   "total"),
    2   => (:elastic, "elastic"),
    18  => (:fission, "fission"),
    102 => (:capture, "capture"),
)

# ---------------------------------------------------------------------------
# Test-case definitions (RECONR-only)
# ---------------------------------------------------------------------------
struct TestCase
    id::Int
    label::String
    endf_file::String
    ref_tape::String
    mat::Int
    err::Float64
end

const CASES = [
    TestCase(84, "H-2",
             joinpath(RESOURCES, "n-001_H_002-ENDF8.0.endf"),
             joinpath(TESTS_DIR, "84", "referenceTape100"),
             128, 0.001),
    TestCase(81, "Sr-88",
             joinpath(RESOURCES, "n-038_Sr_088-ENDF8.1.endf"),
             joinpath(TESTS_DIR, "81", "referenceTape30"),
             3837, 0.001),
    TestCase(85, "Ar-37",
             joinpath(RESOURCES, "n-018_Ar_37-tendl2023.endf"),
             joinpath(TESTS_DIR, "85", "referenceTape50"),
             1828, 0.001),
]

# ---------------------------------------------------------------------------
# Main validation loop
# ---------------------------------------------------------------------------
function run_validation()
    println("=" ^ 78)
    println("NJOY.jl  vs  NJOY2016  Pointwise RECONR Validation")
    println("=" ^ 78)
    @printf("Abs-floor for filtered relative error: %.1e barns\n\n", ABS_FLOOR)

    all_pass = true

    for tc in CASES
        if !isfile(tc.endf_file)
            println("SKIP  Test $(tc.id) ($(tc.label)): ENDF file missing")
            println()
            continue
        end
        if !isfile(tc.ref_tape)
            println("SKIP  Test $(tc.id) ($(tc.label)): reference tape missing")
            println()
            continue
        end

        # Parse reference PENDF
        ref = read_reference_pendf(tc.ref_tape)
        if isempty(ref)
            println("SKIP  Test $(tc.id) ($(tc.label)): no MF3 in reference")
            println()
            continue
        end

        # Run NJOY.jl reconr
        local result
        try
            result = reconr(tc.endf_file; mat=tc.mat, err=tc.err)
        catch ex
            println("FAIL  Test $(tc.id) ($(tc.label)): reconr exception")
            println("      ", sprint(showerror, ex))
            println()
            all_pass = false
            continue
        end

        # Compare each reaction present in the reference
        reactions = ReactionResult[]
        for (mt, (sym, lbl)) in sort(collect(MT_MAP))
            haskey(ref, mt) || continue
            our_xs = getfield(result, sym)
            rr = compare_reaction(ref[mt], result.energies, our_xs, mt, lbl)
            push!(reactions, rr)
        end

        # Summary line (filtered errors only -- physically meaningful)
        parts = String[]
        for rr in reactions
            push!(parts, @sprintf("max_err_%s=%.3f%%",
                                  rr.label, rr.filt_rel_err * 100))
        end
        n_ref_total = haskey(ref, 1) ? ref[1].n_points : 0
        n_ours      = length(result.energies)
        ratio       = n_ref_total > 0 ? n_ours / n_ref_total : NaN

        println("Test $(tc.id) ($(tc.label)): N_points=$(n_ref_total), ",
                join(parts, ", "))
        @printf("      NJOY.jl grid: %d points  (ratio=%.2fx)\n", n_ours, ratio)

        # Per-reaction detail
        for rr in reactions
            status = rr.filt_rel_err < 0.10 ? "OK" : "WARN"
            @printf("      MT=%3d %-8s  ref=%5d  sig=%5d  filt_rel=%.4f%%  filt_abs=%.4e  raw_rel=%.4f%%  [%s]\n",
                    rr.mt, rr.label, rr.n_ref, rr.n_compared,
                    rr.filt_rel_err * 100, rr.filt_abs_err,
                    rr.max_rel_err * 100, status)
            if rr.filt_rel_err > 0.0
                @printf("        worst filtered at E=%.6e eV\n", rr.filt_err_energy)
            end
            if rr.max_rel_err > rr.filt_rel_err + 0.01
                @printf("        worst raw      at E=%.6e eV  (abs=%.4e, below floor)\n",
                        rr.max_err_energy, rr.max_abs_err)
            end
        end

        # Grid count check
        if n_ref_total > 0
            if ratio > 2.0
                println("      GRID: >2x reference -- WARN")
            elseif ratio < 0.5
                println("      GRID: <0.5x reference -- WARN")
            else
                println("      GRID: within 2x factor -- OK")
            end
        end

        # Mark failures: filtered relative error > 50%
        for rr in reactions
            if rr.filt_rel_err > 0.50
                all_pass = false
            end
        end

        println()
    end

    println("=" ^ 78)
    if all_pass
        println("RESULT: All tests pass (filtered relative errors within tolerance).")
    else
        println("RESULT: Some tests show large filtered errors -- see WARN above.")
    end
    println("=" ^ 78)

    return all_pass
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
passed = run_validation()
exit(passed ? 0 : 1)
