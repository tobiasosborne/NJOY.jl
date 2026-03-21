# diagnose_test.jl — Verbose diagnostic runner for individual NJOY validation tests
#
# Usage: julia --project=. test/validation/diagnose_test.jl [test_number]
#
# Prints detailed info at every step: parsed params, grid sizes, sample XS,
# comparison metrics, timing. Designed for debugging pipeline failures.

using NJOY
using Printf

include("test_executor.jl")
include("reference_comparator.jl")

const SEP = "─" ^ 78
const THICK = "═" ^ 78

# =========================================================================
# Pretty-printing helpers
# =========================================================================

function banner(title::String)
    println("\n", THICK)
    println("  ", title)
    println(THICK)
end

function section(title::String)
    println("\n", SEP)
    println("  ", title)
    println(SEP)
end

function show_pendf_summary(label::String, result; indent="  ")
    e = result.energies
    println("$(indent)$(label): $(length(e)) energy points")
    println("$(indent)  E range: $(e[1]) — $(e[end]) eV")

    # Sample at key energies: thermal, 1eV, 1keV, 100keV, 1MeV, 10MeV
    probes = [0.0253, 1.0, 1e3, 1e5, 1e6, 1e7]
    mt_names = [:total, :elastic, :fission, :capture]

    # Build interpolation arrays for each MT
    cols = Dict{Symbol, Vector{Float64}}()
    for sym in mt_names
        if hasproperty(result, sym)
            cols[sym] = getproperty(result, sym)
        elseif result isa PointwiseMaterial
            mt_num = Dict(:total=>1, :elastic=>2, :fission=>18, :capture=>102)[sym]
            idx = findfirst(==(mt_num), result.mt_list)
            cols[sym] = idx !== nothing ? result.cross_sections[:, idx] :
                                          zeros(length(e))
        end
    end

    # Header
    @printf("%s  %12s  %12s  %12s  %12s  %12s\n", indent,
            "E (eV)", "total", "elastic", "fission", "capture")
    for ep in probes
        ep < e[1] || ep > e[end] && continue
        vals = Float64[]
        for sym in mt_names
            haskey(cols, sym) || (push!(vals, NaN); continue)
            push!(vals, _interp_at(e, cols[sym], ep))
        end
        @printf("%s  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n",
                indent, ep, vals...)
    end
end

function show_comparison_detail(report::ComparisonReport; indent="  ")
    println("$(indent)Reference tape: $(report.ref_tape)")
    println("$(indent)Points: ref=$(report.n_ref_total_points)  ours=$(report.n_ours_points)  ratio=$(round(report.grid_ratio, digits=2))x")
    println("$(indent)Overall: $(uppercase(string(report.overall)))")
    println()
    if isempty(report.reactions)
        println("$(indent)  (no MF3 reactions found in reference)")
        return
    end
    @printf("%s  %-5s %-8s  %6s  %6s  %10s  %10s  %14s  %s\n", indent,
            "MT", "name", "n_ref", "n_sig", "filt_err%", "mean_err%", "worst_E(eV)", "verdict")
    for rc in report.reactions
        v = uppercase(string(rc.classification))
        @printf("%s  MT=%-3d %-8s  %6d  %6d  %10.4f  %10.4f  %14.6e  %s\n",
                indent, rc.mt, rc.label, rc.n_ref, rc.n_significant,
                rc.filt_rel_err * 100, rc.mean_rel_err * 100,
                rc.filt_worst_energy, v)
    end
end

# =========================================================================
# Main diagnostic flow
# =========================================================================

function diagnose(test_number::Int)
    t_start = time()

    banner("NJOY.jl Diagnostic: Test $test_number")

    # --- Parse test case ---
    section("1. Test Case Parsing")
    tc = parse_njoy_test(test_number)
    println("  Modules: ", join(tc.modules, " → "))
    println("  Category: ", deck_category(tc.deck))
    println("  MAT: $(tc.mat)  err: $(tc.err)  AWR: $(tc.awr)")
    println("  Status: $(tc.status)")

    println("\n  ENDF files:")
    for (unit, path) in sort(collect(tc.endf_files))
        sz = isfile(path) ? filesize(path) : 0
        println("    tape$unit → $(basename(path))  ($(round(sz/1024, digits=1)) KB)")
    end

    println("\n  Reference tapes:")
    for (name, path) in sort(collect(tc.reference_tapes))
        sz = isfile(path) ? filesize(path) : 0
        println("    $name  ($(round(sz/1024, digits=1)) KB)")
    end

    println("\n  Parsed module params:")
    for mc in tc.deck.calls
        if mc.name == :reconr
            p = parse_reconr(mc)
            println("    reconr: MAT=$(p.mat), err=$(p.err), nendf=$(p.nendf), npend=$(p.npend)")
        elseif mc.name == :broadr
            p = parse_broadr(mc)
            println("    broadr: MAT=$(p.mat), ntemp=$(p.ntemp), tol=$(p.tol), thnmax=$(p.thnmax)")
            println("            temps=$(p.temperatures)")
        elseif mc.name == :heatr
            p = parse_heatr(mc)
            println("    heatr: MAT=$(p.mat), nqa=$(p.nqa), mts=$(p.mts)")
        elseif mc.name == :thermr
            p = parse_thermr(mc)
            println("    thermr: MAT=$(p.mat), temps=$(p.temperatures), emax=$(p.emax)")
        elseif mc.name == :unresr
            p = parse_unresr(mc)
            println("    unresr: MAT=$(p.mat), ntemp=$(p.ntemp), nsigz=$(p.nsigz)")
            println("            temps=$(p.temperatures), sigz=$(p.sigz)")
        elseif mc.name == :purr
            p = parse_purr(mc)
            println("    purr: MAT=$(p.mat), ntemp=$(p.ntemp), nsigz=$(p.nsigz)")
        elseif mc.name == :acer
            p = parse_acer(mc)
            println("    acer: MAT=$(p.mat), iopt=$(p.iopt), temp=$(p.temp)")
        elseif mc.name == :gaspr
            p = parse_gaspr(mc)
            println("    gaspr: nendf=$(p.nendf), nin=$(p.npendf_in), nout=$(p.npendf_out)")
        else
            println("    $(mc.name): $(length(mc.raw_cards)) cards")
        end
    end

    # --- Execute step by step ---
    section("2. Pipeline Execution")

    endf_file = find_endf_file(tc)
    if endf_file === nothing
        println("  ✘ No ENDF file found — cannot proceed")
        return
    end
    println("  ENDF file: $(endf_file)")

    pendf = nothing
    step_num = 0

    for mc in tc.deck.calls
        step_num += 1
        t0 = time()
        mod = mc.name

        if mod == :moder
            @printf("  [%d] %-8s ... pass-through (%.1fms)\n", step_num, mod, (time()-t0)*1000)
            continue
        end

        if mod in SKIP_MODULES
            @printf("  [%d] %-8s ... skipped (visualization)\n", step_num, mod)
            continue
        end

        @printf("  [%d] %-8s ... ", step_num, mod)
        flush(stdout)

        if mod == :reconr
            r, s = execute_reconr(endf_file, mc; mat_override=tc.mat, err_override=tc.err)
            dt = (time() - t0) * 1000
            if s.status == :ok
                pendf = r
                @printf("OK  %s  (%.0fms)\n", s.message, dt)
                show_pendf_summary("RECONR output", pendf)
            else
                @printf("FAIL  %s  (%.0fms)\n", s.message, dt)
            end
        elseif pendf === nothing
            @printf("SKIP  (no upstream pendf)\n")
        else
            pm = ensure_pendf(pendf, tc.mat)
            r, s = _execute_downstream(mc, pm, tc.mat; awr=tc.awr)
            dt = (time() - t0) * 1000
            status_str = uppercase(string(s.status))
            @printf("%-4s  %s  (%.0fms)\n", status_str, s.message, dt)
            if r !== nothing && r !== pm
                pendf = r
                show_pendf_summary("After $(mod)", pendf isa PointwiseMaterial ? pendf_to_result(pendf) : pendf)
            end
        end
    end

    # --- Comparison ---
    section("3. Reference Comparison")

    result = if pendf isa PointwiseMaterial
        pendf_to_result(pendf)
    else
        pendf
    end

    if result === nothing
        println("  ✘ No pipeline output to compare")
        return
    end

    if isempty(tc.reference_tapes)
        println("  ✘ No reference tapes available")
        return
    end

    for (name, path) in sort(collect(tc.reference_tapes))
        println("\n  --- vs $name ---")
        ref_data = read_ref_pendf(path)
        if isempty(ref_data)
            println("    (no MF3 data in reference tape)")
            continue
        end

        println("    Reference MTs found: ", sort(collect(keys(ref_data))))
        for (mt, rd) in sort(collect(ref_data))
            @printf("      MT=%3d: %d points, E=[%.3e, %.3e]\n",
                    mt, length(rd.energies), rd.energies[1], rd.energies[end])
        end

        report = compare_all(result, path, test_number)
        println()
        show_comparison_detail(report; indent="    ")

        # Show worst-point details for any failing reactions
        for rc in report.reactions
            rc.filt_rel_err <= THRESH_PASS && continue
            haskey(ref_data, rc.mt) || continue
            rd = ref_data[rc.mt]
            field = get(MT_FIELD_MAP, rc.mt, nothing)
            field === nothing && continue
            our_xs = getproperty(result, field[1])

            println("\n    Worst points for MT=$(rc.mt) ($(rc.label)):")
            # Find top-5 worst points
            errs = Tuple{Float64,Float64,Float64,Float64}[]
            for k in 1:length(rd.energies)
                re = rd.energies[k]
                (re <= 0 || !isfinite(re)) && continue
                xo = _interp_at(result.energies, our_xs, re)
                abs(rd.xs[k]) < ABS_FLOOR && continue
                rel = abs(xo - rd.xs[k]) / abs(rd.xs[k])
                push!(errs, (rel, re, xo, rd.xs[k]))
            end
            sort!(errs; rev=true)
            @printf("    %14s  %14s  %14s  %10s\n", "E (eV)", "ours", "ref", "rel_err%")
            for (rel, e, xo, xr) in errs[1:min(10, length(errs))]
                @printf("    %14.6e  %14.6e  %14.6e  %10.3f%%\n", e, xo, xr, rel*100)
            end
        end
    end

    # --- Summary ---
    dt_total = time() - t_start
    section("4. Summary")
    @printf("  Total time: %.1f seconds\n", dt_total)
    println("  Pipeline steps: $step_num")
    if result !== nothing
        println("  Final grid: $(length(result.energies)) points")
    end
    println()
end

# =========================================================================
# CLI
# =========================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    tn = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
    diagnose(tn)
end
