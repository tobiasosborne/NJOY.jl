# run_all.jl -- Full 85-test NJOY validation pipeline
#
# For each test: parse input deck, execute reconr/broadr, compare against
# reference PENDF tapes, report PASS/MARGINAL/FAIL/SKIP/ERROR.
#
# Usage:
#   julia --project test/validation/run_all.jl [--verbose] [--tests 81,84,85]

using NJOY
using Printf

include("test_executor.jl")
include("reference_comparator.jl")

const ALL_TESTS = vcat(collect(1:76), collect(78:85))

struct ValidationResult
    test_number::Int; category::Symbol
    execution::ExecutionResult
    comparison::Union{ComparisonReport, Nothing}
    overall_status::Symbol
end

# Full validation pipeline for one test.
function validate_test(test_number::Int; verbose::Bool=false)
    tc = parse_njoy_test(test_number)
    cat = deck_category(tc.deck)
    if isempty(tc.deck.calls)
        return ValidationResult(test_number, cat,
            ExecutionResult(test_number, :skip, "no input deck",
                           StepResult[], nothing, 0.0),
            nothing, :skip)
    end
    exec = execute_test(tc)
    comp = nothing
    if exec.pendf !== nothing && !isempty(tc.reference_tapes)
        comp = compare_best_ref(exec.pendf, tc.reference_tapes, test_number)
    end
    overall = exec.status == :skip ? :skip :
              exec.status == :error ? :error :
              comp === nothing ? :skip : comp.overall
    if verbose
        st = uppercase(string(overall))
        ms = join(tc.modules, ",")
        if comp !== nothing && !isempty(comp.reactions)
            w = maximum(r -> r.filt_rel_err, comp.reactions)
            @printf("  %-8s Test %02d [%-12s] mods=%-40s worst=%.3f%%\n",
                    st, test_number, string(cat), ms, w * 100)
        else
            @printf("  %-8s Test %02d [%-12s] mods=%s\n", st, test_number, string(cat), ms)
        end
    end
    ValidationResult(test_number, cat, exec, comp, overall)
end

# Run the full validation pipeline for all specified tests.
function validate_all(; tests::AbstractVector{Int}=ALL_TESTS, verbose::Bool=true)
    results = ValidationResult[]
    if verbose
        println("=" ^ 78)
        println("NJOY.jl Full Validation Pipeline -- $(length(tests)) tests")
        println("=" ^ 78)
        println()
    end
    for tn in tests
        push!(results, validate_test(tn; verbose=verbose))
    end
    verbose && (println(); print_summary(results))
    results
end

# Print a comprehensive summary of all validation results.
function print_summary(results::Vector{ValidationResult})
    println("=" ^ 78)
    println("VALIDATION SUMMARY")
    println("=" ^ 78)
    np = count(r -> r.overall_status == :pass, results)
    nm = count(r -> r.overall_status == :marginal, results)
    nf = count(r -> r.overall_status == :fail, results)
    ns = count(r -> r.overall_status == :skip, results)
    ne = count(r -> r.overall_status == :error, results)
    nd = count(r -> r.overall_status == :no_data, results)
    @printf("\n  PASS:     %3d  (filtered rel error < 1%%)\n", np)
    @printf("  MARGINAL: %3d  (1%% - 10%%)\n", nm)
    @printf("  FAIL:     %3d  (> 10%%)\n", nf)
    @printf("  ERROR:    %3d  (execution error)\n", ne)
    @printf("  SKIP:     %3d  (not runnable / no reference)\n", ns)
    @printf("  NO_DATA:  %3d  (no MF3 in reference)\n", nd)
    @printf("  TOTAL:    %3d\n\n", length(results))

    # Category breakdown
    cats = Dict{Symbol,Int}()
    for r in results; cats[r.category] = get(cats, r.category, 0) + 1; end
    println("By category:")
    for (c, n) in sort(collect(cats))
        @printf("  %-20s %3d\n", string(c), n)
    end
    println()

    # Module frequency
    all_mods = Symbol[]
    for r in results
        tc = parse_njoy_test(r.test_number)
        append!(all_mods, tc.modules)
    end
    mc = Dict{Symbol,Int}()
    for m in all_mods; mc[m] = get(mc, m, 0) + 1; end
    println("Module frequency:")
    for (mod, n) in sort(collect(mc); by=x->-x[2])
        av = mod in AVAILABLE_MODULES ? "available" :
             (mod in SKIP_MODULES ? "skip(viz)" : "MISSING")
        @printf("  %-10s %3d  [%s]\n", string(mod), n, av)
    end
    println()

    # FAIL/MARGINAL details
    bad = filter(r -> r.overall_status in (:fail, :marginal), results)
    if !isempty(bad)
        println("FAIL/MARGINAL details:")
        println("-" ^ 60)
        for r in bad
            r.comparison !== nothing ? println(format_comparison(r.comparison)) :
                @printf("  Test %02d: %s\n", r.test_number, uppercase(string(r.overall_status)))
        end
    end

    # Error details
    errs = filter(r -> r.overall_status == :error, results)
    if !isempty(errs)
        println("ERROR details:")
        println("-" ^ 60)
        for r in errs
            @printf("  Test %02d: %s\n", r.test_number, r.execution.message)
            for s in r.execution.steps
                s.status == :error && @printf("    %s: %s\n", s.module_name, s.message)
            end
        end
        println()
    end

    # Skip summary
    skipped = filter(r -> r.overall_status == :skip, results)
    if !isempty(skipped)
        println("Skipped tests:")
        for r in skipped
            tc = parse_njoy_test(r.test_number)
            @printf("  Test %02d: %-40s [%s]\n", r.test_number,
                    r.execution.message, join(tc.modules, ","))
        end
        println()
    end

    println("=" ^ 78)
    if nf == 0 && ne == 0
        println("VERDICT: All runnable tests pass or are marginal.")
    else
        @printf("VERDICT: %d FAIL, %d ERROR.\n", nf, ne)
    end
    println("=" ^ 78)
end

# =========================================================================
# CLI entry point
# =========================================================================

function main()
    verbose = true; test_nums = ALL_TESTS
    i = 1
    while i <= length(ARGS)
        ARGS[i] == "--quiet" && (verbose = false)
        ARGS[i] == "--verbose" && (verbose = true)
        if ARGS[i] == "--tests" && i < length(ARGS)
            i += 1; test_nums = [parse(Int, strip(s)) for s in split(ARGS[i], ',')]
        end
        i += 1
    end
    results = validate_all(; tests=test_nums, verbose=verbose)
    nf = count(r -> r.overall_status == :fail, results)
    exit(nf > 0 ? 1 : 0)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
