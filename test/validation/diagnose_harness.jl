# diagnose_harness.jl — Per-module diagnostic harness with Fortran oracle
#
# Usage:
#   julia --project=. test/validation/diagnose_harness.jl [test_number]
#   julia --project=. test/validation/diagnose_harness.jl [test_number] --isolate
#
# Five-phase diagnostic:
#   1. Deck Analysis       — parse input, classify test
#   2. Fortran Oracle      — run Fortran NJOY with truncated decks (cached)
#   3. A/B Comparison      — Julia vs Fortran at each module boundary
#   4. Deep Inspection     — zoom into the culprit module
#   5. Enhanced Comparison — all MTs, error classification, regional analysis

using NJOY
using Printf

include("module_inspector.jl")
include("reference_comparator.jl")
include("error_classifier.jl")
include("fortran_oracle.jl")

const SEP_H = "─" ^ 78
const THICK_H = "═" ^ 78

banner(t) = println("\n", THICK_H, "\n  ", t, "\n", THICK_H)
section(t) = println("\n", SEP_H, "\n  ", t, "\n", SEP_H)

# =========================================================================
# Phase 1: Deck Analysis
# =========================================================================

function phase_deck_analysis(test_number::Int)
    section("Phase 1: Deck Analysis")

    tc = parse_njoy_test(test_number)
    println("  Modules: ", join(tc.modules, " → "))
    println("  MAT: $(tc.mat)  err: $(tc.err)  AWR: $(tc.awr)")
    println("  Status: $(tc.status)")

    # ENDF files
    println("\n  ENDF files:")
    for (unit, path) in sort(collect(tc.endf_files))
        sz = isfile(path) ? filesize(path) : 0
        println("    tape$unit → $(basename(path))  ($(round(sz/1024, digits=1)) KB)")
    end

    # Reference tapes
    println("\n  Reference tapes:")
    for (name, path) in sort(collect(tc.reference_tapes))
        sz = isfile(path) ? filesize(path) : 0
        groupr = is_groupr_tape(path) ? " [GROUPR — skip]" : ""
        println("    $name  ($(round(sz/1024, digits=1)) KB)$groupr")
    end

    # Module params
    println("\n  Parsed params:")
    for mc in tc.deck.calls
        mc.name == :moder && continue
        mc.name in SKIP_MODULES && continue
        try
            if mc.name == :reconr
                p = parse_reconr(mc)
                @printf("    reconr: MAT=%d, err=%.4f\n", p.mat, p.err)
            elseif mc.name == :broadr
                p = parse_broadr(mc)
                @printf("    broadr: ntemp=%d, tol=%.4f, temps=%s\n",
                        p.ntemp, p.tol, p.temperatures)
            elseif mc.name == :heatr
                p = parse_heatr(mc)
                @printf("    heatr: mts=%s\n", p.mts)
            elseif mc.name == :thermr
                p = parse_thermr(mc)
                @printf("    thermr: model=%s, temps=%s, emax=%.2f\n",
                        p.nin_thermal > 0 ? "S(a,b)" : "free_gas",
                        p.temperatures, p.emax)
            else
                @printf("    %s: %d cards\n", mc.name, length(mc.raw_cards))
            end
        catch; @printf("    %s: (parse error)\n", mc.name); end
    end

    tc
end

# =========================================================================
# Phase 2: Fortran Oracle
# =========================================================================

function phase_fortran_oracle(test_number::Int)
    section("Phase 2: Fortran Oracle References")

    njoy_bin = ""
    try
        njoy_bin = ensure_njoy_binary()
    catch ex
        println("  NJOY2016 binary not available: $ex")
        println("  Falling back to final reference tapes only.")
        return Dict{String, String}()
    end

    println("  Binary: $njoy_bin")
    refs = generate_module_references(test_number; njoy_binary=njoy_bin)

    println("\n  Generated per-module references:")
    for (mod, path) in sort(collect(refs); by=first)
        sz = isfile(path) ? filesize(path) : 0
        println("    after $mod → $(round(sz/1024, digits=1)) KB")
    end
    refs
end

# =========================================================================
# Phase 3: Instrumented Pipeline + A/B Comparison
# =========================================================================

struct ModuleABResult
    module_name::Symbol
    step_index::Int
    julia_diag::ModuleDiagnostics
    # Per-MT error vs Fortran oracle at this exact stage
    mt_errors::Dict{Int, Float64}    # MT -> max filtered rel error
    overall_error::Float64
    verdict::Symbol                  # :pass, :marginal, :fail, :no_ref
end

function phase_ab_comparison(tc::NJOYTestCase, oracle_refs::Dict{String, String})
    section("Phase 3: A/B Comparison (Julia vs Fortran at each stage)")

    endf_file = find_endf_file(tc)
    if endf_file === nothing
        println("  No ENDF file found")
        return ModuleABResult[]
    end

    results = ModuleABResult[]
    pendf = nothing
    step = 0
    # Track module occurrence count for step labels (e.g. thermr_2)
    mod_counts = Dict{Symbol, Int}()

    for mc in tc.deck.calls
        mc.name == :moder && continue
        mc.name in SKIP_MODULES && continue

        step += 1
        mod_counts[mc.name] = get(mod_counts, mc.name, 0) + 1
        # Build step label matching fortran_oracle's _step_label
        step_label = mod_counts[mc.name] > 1 ?
            "$(mc.name)_$(mod_counts[mc.name])" : string(mc.name)

        # Run Julia module
        diag, new_pendf, step_result = if mc.name == :reconr
            d, r, s = inspect_reconr(endf_file, mc, tc.mat, tc.err)
            d, r, s
        elseif pendf === nothing
            println("  [$step] $(mc.name) ... SKIP (no upstream)")
            continue
        else
            pm = ensure_pendf(pendf, tc.mat)
            if mc.name == :broadr
                d, r, s = inspect_broadr(pm, mc; awr=tc.awr)
                d, r, s
            elseif mc.name == :heatr
                d, r, s = inspect_heatr(pm, mc; awr=tc.awr)
                d, r, s
            elseif mc.name == :thermr
                d, r, s = inspect_thermr(pm, mc)
                d, r, s
            else
                # Pass-through modules
                println("  [$step] $(mc.name) ... pass-through")
                continue
            end
        end

        new_pendf !== nothing && (pendf = new_pendf)

        # Compare vs Fortran oracle at this stage (lookup by step label)
        mt_errors = Dict{Int, Float64}()
        overall = 0.0
        verdict = :no_ref

        oracle_path = get(oracle_refs, step_label, "")
        if !isempty(oracle_path) && isfile(oracle_path)
            ref_data = read_ref_pendf(oracle_path)
            if !isempty(ref_data)
                result_nt = if pendf isa PointwiseMaterial
                    pendf_to_result(pendf)
                else
                    pendf
                end

                if result_nt !== nothing
                    for (mt, rd) in ref_data
                        field_info = get(MT_FIELD_MAP, mt, nothing)
                        field_info === nothing && continue
                        sym, label = field_info
                        hasproperty(result_nt, sym) || continue
                        our_xs = getproperty(result_nt, sym)
                        rc = compare_reaction(rd, result_nt.energies, our_xs, mt, label)
                        mt_errors[mt] = rc.filt_rel_err
                        rc.filt_rel_err > overall && (overall = rc.filt_rel_err)
                    end
                    verdict = overall <= THRESH_PASS ? :pass :
                              overall <= THRESH_MARGINAL ? :marginal : :fail
                end
            end
        end

        push!(results, ModuleABResult(mc.name, step, diag, mt_errors, overall, verdict))
    end

    # Print A/B table
    println()
    @printf("  %-4s  %-8s  %7s  %7s  %10s  %10s  %10s  %10s  %s\n",
            "Step", "Module", "Jul_pts", "F_pts",
            "MT1_err%", "MT2_err%", "MT18_err%", "MT102_err%", "Verdict")
    @printf("  %s\n", "─" ^ 90)

    for ab in results
        jul_pts = ab.julia_diag.n_grid_out

        # Try to get Fortran point count from oracle
        f_pts = "-"
        # Find matching oracle ref by module name prefix
        for (label, path) in oracle_refs
            if startswith(label, string(ab.module_name))
                ref_data = read_ref_pendf(path)
                if haskey(ref_data, 1)
                    f_pts = string(length(ref_data[1].energies))
                end
                break
            end
        end

        mt1  = haskey(ab.mt_errors, 1)   ? @sprintf("%9.2f%%", ab.mt_errors[1]*100)   : "    N/A  "
        mt2  = haskey(ab.mt_errors, 2)   ? @sprintf("%9.2f%%", ab.mt_errors[2]*100)   : "    N/A  "
        mt18 = haskey(ab.mt_errors, 18)  ? @sprintf("%9.2f%%", ab.mt_errors[18]*100)  : "    N/A  "
        mt102= haskey(ab.mt_errors, 102) ? @sprintf("%9.2f%%", ab.mt_errors[102]*100) : "    N/A  "

        v = uppercase(string(ab.verdict))
        @printf("  [%2d]  %-8s  %7d  %7s  %s  %s  %s  %s  %s\n",
                ab.step_index, ab.module_name, jul_pts, f_pts,
                mt1, mt2, mt18, mt102, v)
    end

    # Identify culprit
    failing = filter(r -> r.verdict == :fail, results)
    if !isempty(failing)
        worst = argmax(r -> r.overall_error, failing)
        println("\n  CULPRIT: $(worst.module_name) ($(round(worst.overall_error*100, digits=1))% max error)")
    elseif !isempty(results)
        println("\n  No clear culprit — all modules within tolerance or no oracle data")
    end

    results
end

# =========================================================================
# Phase 4: Deep Module Inspection
# =========================================================================

function phase_deep_inspection(results::Vector{ModuleABResult})
    section("Phase 4: Deep Module Inspection")

    failing = filter(r -> r.verdict == :fail, results)
    if isempty(failing)
        println("  No failing modules — skipping deep inspection")
        return
    end

    worst = argmax(r -> r.overall_error, failing)
    println("  Inspecting: $(worst.module_name)")
    println()

    # Print full diagnostics
    print(format_diagnostics(worst.julia_diag; indent="  "))

    # Fortran cross-reference
    xref = Dict(
        :reconr => ("reconr.f90", "emerge:4646 (merge bg), lunion:1771 (grid), resxs:2240 (adaptive), sigma:2571 (evaluate)"),
        :broadr => ("broadr.f90", "bsigma:1510 (kernel), broadn:1256 (main), thinb:1677 (thinning)"),
        :heatr  => ("heatr.f90",  "heatr:128 (main)"),
        :thermr => ("thermr.f90", "thermr (main)"),
    )
    if haskey(xref, worst.module_name)
        file, subs = xref[worst.module_name]
        println("\n  Fortran reference: $file")
        println("    Key subroutines: $subs")
    end

    # Per-MT error detail
    if !isempty(worst.mt_errors)
        println("\n  Per-MT errors:")
        for (mt, err) in sort(collect(worst.mt_errors))
            label = get(EXTENDED_MT_NAMES, mt, "MT$mt")
            status = err <= THRESH_PASS ? "PASS" : (err <= THRESH_MARGINAL ? "MARG" : "FAIL")
            @printf("    MT=%-3d %-8s  %8.2f%%  %s\n", mt, label, err*100, status)
        end
    end
end

# =========================================================================
# Phase 5: Enhanced Comparison
# =========================================================================

function phase_enhanced_comparison(tc::NJOYTestCase, pendf_result,
                                   oracle_refs::Dict{String, String})
    section("Phase 5: Enhanced Comparison (all MTs)")

    if pendf_result === nothing
        println("  No pipeline output")
        return
    end

    result_nt = pendf_result isa PointwiseMaterial ? pendf_to_result(pendf_result) : pendf_result

    # Compare against all available reference sources
    all_refs = Dict{String, String}()

    # Oracle per-module references
    for (mod, path) in oracle_refs
        all_refs["oracle_$(mod)"] = path
    end

    # Original test reference tapes
    for (name, path) in tc.reference_tapes
        is_groupr_tape(path) && continue
        all_refs[name] = path
    end

    for (label, path) in sort(collect(all_refs))
        ref_data = read_ref_pendf(path)
        isempty(ref_data) && continue

        println("\n  --- vs $label ($(length(ref_data)) MTs) ---")

        # Show all available MTs in reference
        ref_mts = sort(collect(keys(ref_data)))
        println("    Reference MTs: $ref_mts")

        # Compare primary MTs with error classification
        for mt in [1, 2, 18, 102]
            haskey(ref_data, mt) || continue
            field_info = get(MT_FIELD_MAP, mt, nothing)
            field_info === nothing && continue
            sym, mt_label = field_info
            hasproperty(result_nt, sym) || continue

            our_xs = getproperty(result_nt, sym)
            prof = build_error_profile(result_nt.energies, our_xs,
                                       ref_data[mt].energies, ref_data[mt].xs,
                                       mt, mt_label)
            print(format_error_profile(prof; indent="    "))
        end
    end
end

# =========================================================================
# Main orchestrator
# =========================================================================

function diagnose_deep(test_number::Int; isolate::Bool=false)
    t_start = time()
    banner("NJOY.jl Diagnostic Harness: Test $test_number")

    # Phase 1
    tc = phase_deck_analysis(test_number)

    # Phase 2
    oracle_refs = phase_fortran_oracle(test_number)

    # Phase 3
    ab_results = phase_ab_comparison(tc, oracle_refs)

    # Phase 4
    phase_deep_inspection(ab_results)

    # Phase 5: get final pipeline result for enhanced comparison
    endf_file = find_endf_file(tc)
    pendf_result = nothing
    if endf_file !== nothing
        try
            er = execute_test(tc)
            pendf_result = er.pendf
        catch; end
    end
    phase_enhanced_comparison(tc, pendf_result, oracle_refs)

    # Summary
    section("Summary")
    @printf("  Total time: %.1f seconds\n", time() - t_start)
    @printf("  Oracle references: %d\n", length(oracle_refs))
    @printf("  A/B comparisons: %d\n", length(ab_results))

    n_pass = count(r -> r.verdict == :pass, ab_results)
    n_fail = count(r -> r.verdict == :fail, ab_results)
    n_marg = count(r -> r.verdict == :marginal, ab_results)
    n_noref = count(r -> r.verdict == :no_ref, ab_results)
    @printf("  Verdicts: %d pass, %d marginal, %d fail, %d no_ref\n",
            n_pass, n_marg, n_fail, n_noref)

    println()
end

# =========================================================================
# CLI
# =========================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    tn = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
    isolate = "--isolate" in ARGS
    diagnose_deep(tn; isolate=isolate)
end
