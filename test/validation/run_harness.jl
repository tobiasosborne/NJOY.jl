# run_harness.jl — Batch diagnostic harness runner with dashboard
#
# Usage:
#   julia --project=. test/validation/run_harness.jl
#   julia --project=. test/validation/run_harness.jl --tests 1,2,7,12,84
#
# Runs the per-module diagnostic harness across multiple tests and produces
# a summary dashboard showing module culprit ranking and formalism error patterns.

using Printf

include("diagnose_harness.jl")
include("test_classifier.jl")

const ALL_TESTS = [i for i in 1:85 if i != 77]

# =========================================================================
# Batch result
# =========================================================================

struct BatchResult
    test_number::Int
    classification::TestClassification
    ab_results::Vector{ModuleABResult}
    culprit::Union{Symbol, Nothing}
    culprit_error::Float64
    elapsed::Float64
end

# =========================================================================
# Runner
# =========================================================================

function run_harness(; tests::Vector{Int}=ALL_TESTS, verbose::Bool=false)
    results = BatchResult[]
    t_total = time()

    println("═" ^ 78)
    println("  NJOY.jl Diagnostic Harness — Batch Run ($(length(tests)) tests)")
    println("═" ^ 78)

    for (i, tn) in enumerate(tests)
        t0 = time()
        @printf("\n[%d/%d] Test %02d ... ", i, length(tests), tn)
        flush(stdout)

        tc_class = classify_test(tn)

        # Skip empty tests
        if tc_class.chain_pattern == :empty
            @printf("SKIP (empty)\n")
            push!(results, BatchResult(tn, tc_class, ModuleABResult[], nothing, 0.0, 0.0))
            continue
        end

        try
            tc = parse_njoy_test(tn)
            oracle_refs = Dict{String, String}()
            try
                njoy_bin = ensure_njoy_binary()
                oracle_refs = generate_module_references(tn; njoy_binary=njoy_bin)
            catch; end

            ab_results = if !isempty(tc.deck.calls) && tc.mat > 0
                # Simplified A/B comparison (no verbose output)
                _run_ab_quiet(tc, oracle_refs)
            else
                ModuleABResult[]
            end

            # Find culprit
            failing = filter(r -> r.verdict == :fail, ab_results)
            culprit = nothing
            culprit_err = 0.0
            if !isempty(failing)
                worst = argmax(r -> r.overall_error, failing)
                culprit = worst.module_name
                culprit_err = worst.overall_error
            end

            elapsed = time() - t0
            push!(results, BatchResult(tn, tc_class, ab_results, culprit, culprit_err, elapsed))

            status = culprit !== nothing ? "FAIL ($(culprit): $(round(culprit_err*100, digits=1))%)" :
                     isempty(ab_results) ? "SKIP" : "PASS"
            @printf("%-50s (%.1fs)\n", status, elapsed)

        catch ex
            elapsed = time() - t0
            msg = first(split(sprint(showerror, ex), '\n'))
            @printf("ERROR: %s (%.1fs)\n", msg[1:min(45, length(msg))], elapsed)
            push!(results, BatchResult(tn, tc_class, ModuleABResult[], nothing, 0.0, elapsed))
        end
    end

    print_dashboard(results)
    @printf("\nTotal time: %.1f seconds\n", time() - t_total)
    results
end

"""Run A/B comparison without verbose output."""
function _run_ab_quiet(tc::NJOYTestCase, oracle_refs::Dict{String, String})
    endf_file = find_endf_file(tc)
    endf_file === nothing && return ModuleABResult[]

    results = ModuleABResult[]
    pendf = nothing
    step = 0
    mod_counts = Dict{Symbol, Int}()

    for mc in tc.deck.calls
        mc.name == :moder && continue
        mc.name in SKIP_MODULES && continue
        step += 1
        mod_counts[mc.name] = get(mod_counts, mc.name, 0) + 1
        step_label = mod_counts[mc.name] > 1 ?
            "$(mc.name)_$(mod_counts[mc.name])" : string(mc.name)

        # Run Julia module (quiet)
        diag, new_pendf, _ = try
            if mc.name == :reconr
                d, r, s = inspect_reconr(endf_file, mc, tc.mat, tc.err)
                d, r, s
            elseif pendf === nothing
                continue
            else
                pm = ensure_pendf(pendf, tc.mat)
                if mc.name == :broadr
                    inspect_broadr(pm, mc; awr=tc.awr)
                elseif mc.name == :heatr
                    inspect_heatr(pm, mc; awr=tc.awr)
                elseif mc.name == :thermr
                    inspect_thermr(pm, mc)
                else
                    continue
                end
            end
        catch
            continue
        end

        new_pendf !== nothing && (pendf = new_pendf)

        # Compare vs oracle
        mt_errors = Dict{Int, Float64}()
        overall = 0.0
        verdict = :no_ref

        oracle_path = get(oracle_refs, step_label, "")
        if !isempty(oracle_path) && isfile(oracle_path)
            try
                ref_data = read_ref_pendf(oracle_path)
                result_nt = pendf isa PointwiseMaterial ? pendf_to_result(pendf) : pendf
                if result_nt !== nothing && !isempty(ref_data)
                    for (mt, rd) in ref_data
                        fi = get(MT_FIELD_MAP, mt, nothing)
                        fi === nothing && continue
                        sym, label = fi
                        hasproperty(result_nt, sym) || continue
                        rc = compare_reaction(rd, result_nt.energies,
                                              getproperty(result_nt, sym), mt, label)
                        mt_errors[mt] = rc.filt_rel_err
                        rc.filt_rel_err > overall && (overall = rc.filt_rel_err)
                    end
                    verdict = overall <= THRESH_PASS ? :pass :
                              overall <= THRESH_MARGINAL ? :marginal : :fail
                end
            catch; end
        end

        push!(results, ModuleABResult(mc.name, step, diag, mt_errors, overall, verdict))
    end
    results
end

# =========================================================================
# Dashboard
# =========================================================================

function print_dashboard(results::Vector{BatchResult})
    println("\n", "═" ^ 78)
    println("  DASHBOARD")
    println("═" ^ 78)

    # Overall stats
    n_total = length(results)
    n_pass = count(r -> r.culprit === nothing && !isempty(r.ab_results), results)
    n_fail = count(r -> r.culprit !== nothing, results)
    n_skip = count(r -> isempty(r.ab_results), results)
    @printf("\n  Tests: %d total, %d pass, %d fail, %d skip\n", n_total, n_pass, n_fail, n_skip)

    # Module culprit ranking
    println("\n  MODULE CULPRIT RANKING")
    println("  ", "─" ^ 65)
    @printf("  %-10s  %8s  %10s  %s\n", "Module", "Culprit#", "Avg_err%", "Worst_test")
    println("  ", "─" ^ 65)

    module_stats = Dict{Symbol, Vector{Tuple{Int, Float64}}}()
    for r in results
        r.culprit === nothing && continue
        push!(get!(module_stats, r.culprit, Tuple{Int, Float64}[]),
              (r.test_number, r.culprit_error))
    end
    for (mod, entries) in sort(collect(module_stats); by=p -> -length(p.second))
        n = length(entries)
        avg = sum(e -> e[2], entries) / n
        worst_tn, worst_err = argmax(e -> e[2], entries)
        @printf("  %-10s  %8d  %10.1f  Test %02d (%.1f%%)\n",
                mod, n, avg*100, worst_tn, worst_err*100)
    end

    # Formalism error summary
    println("\n  FORMALISM ERROR SUMMARY")
    println("  ", "─" ^ 65)
    @printf("  %-12s  %6s  %10s  %s\n", "Formalism", "Tests", "Avg_err%", "Worst_test")
    println("  ", "─" ^ 65)

    form_stats = Dict{Symbol, Vector{Tuple{Int, Float64}}}()
    for r in results
        isempty(r.ab_results) && continue
        worst_err = maximum((ab.overall_error for ab in r.ab_results); init=0.0)
        push!(get!(form_stats, r.classification.formalism, Tuple{Int, Float64}[]),
              (r.test_number, worst_err))
    end
    for (form, entries) in sort(collect(form_stats); by=first)
        n = length(entries)
        avg = sum(e -> e[2], entries) / n
        worst_tn, worst_err = argmax(e -> e[2], entries)
        @printf("  %-12s  %6d  %10.1f  Test %02d (%.1f%%)\n",
                form, n, avg*100, worst_tn, worst_err*100)
    end

    # Detail for failing tests
    fails = filter(r -> r.culprit !== nothing, results)
    if !isempty(fails)
        println("\n  FAILING TESTS DETAIL")
        println("  ", "─" ^ 65)
        for r in sort(fails; by=r -> -r.culprit_error)
            @printf("  Test %02d  %-10s  %.1f%%  (%s, %s)\n",
                    r.test_number, r.culprit, r.culprit_error*100,
                    r.classification.formalism, r.classification.material_class)
        end
    end
end

# =========================================================================
# CLI
# =========================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    tests = ALL_TESTS
    for (i, arg) in enumerate(ARGS)
        if arg == "--tests" && i < length(ARGS)
            tests = [parse(Int, s) for s in split(ARGS[i+1], ",")]
        end
    end
    run_harness(; tests=tests)
end
