# reference_test.jl -- Faithful Julia port of njoy-reference/tests/execute.py.
#
# For each Fortran reference test NN:
#   1. run_njoy("njoy-reference/tests/NN/input"; work_dir=<tmp>)
#      (build_tape_manager stages resource tapes exactly like CMakeLists.txt does)
#   2. For every `referenceTapeNN` in the test dir, find the matching produced
#      `tapeNN` in work_dir and compare line-by-line with execute.py semantics:
#      date wildcarded to XX/XX/XX, floats compared with rtol/atol, non-numeric
#      residue must match.
#
# Verbose output: per-module start/done with elapsed time, plus a heartbeat
# that prints "[STILL RUNNING after Ns] last: <msg>" every 10s when no
# module transition has happened. Soft timeout at 600s (prints a loud
# warning — does NOT kill the run; Ctrl-C to stop).

using NJOY
using Printf

const NJOY_REF_ROOT = normpath(joinpath(@__DIR__, "..", "..", "njoy-reference"))
const TESTS_ROOT    = joinpath(NJOY_REF_ROOT, "tests")

# =========================================================================
# execute.py::lineEquivalence port
# =========================================================================

# Same regex as execute.py's floatPattern (including optional period+digits and
# optional designator-less [+-]\d+ exponent for ENDF a11 format).
const FLOAT_RE = r"([-+]?\d+(\.\d+)?(([eE])?[-+]\d+)?)"
const DATE_RE  = r"\d{2}/\d{2}/\d{2}"

"""
    _parse_endf_or_normal(s::AbstractString) -> Float64 or Nothing

Parse a float that may be in ENDF compact format (`1.234567+8`) or regular
Julia format. Matches execute.py::makeFloat semantics, including Python's
`float()` overflow behavior (returns ±Inf for exponents >= 309, 0.0 for
<= -324) — Julia's `tryparse` returns nothing on overflow, which would
spuriously fail lines where the greedy regex absorbs trailing ENDF trailer
digits into the exponent (e.g. "6.706487+11306" where 1306 is the MAT).
"""
function _parse_endf_or_normal(s::AbstractString)
    v = tryparse(Float64, s)
    v !== nothing && return v
    # Try ENDF compact: insert E before the last unprefixed +/-
    for i in length(s)-1:-1:2
        c = s[i]
        if (c == '+' || c == '-') && !(s[i-1] in ('e', 'E')) && isdigit(s[i-1])
            candidate = string(s[1:i-1], 'E', s[i:end])
            v = tryparse(Float64, candidate)
            v !== nothing && return v
            # Overflow path — match Python float() semantics.
            mant   = tryparse(Float64, s[1:i-1])
            expval = tryparse(Int,     s[i:end])
            if mant !== nothing && expval !== nothing
                if expval >= 309
                    return mant == 0.0 ? 0.0 : sign(mant) * Inf
                elseif expval <= -324
                    return 0.0
                end
            end
        end
    end
    nothing
end

"""
    line_equivalence(ref, trial; rel_tol=1e-9, abs_tol=1e-10) -> Bool

Port of execute.py::lineEquivalence. True if `ref` and `trial` are equivalent
under the NJOY reference-comparison semantics:

  - If byte-equal, true.
  - Otherwise: extract floats from both via FLOAT_RE. If count differs → false.
    All floats must pass `isapprox(a, b; rtol, atol)`. The non-numeric residue
    (input with floats substituted to empty) must also match byte-for-byte.
  - If ref has no floats → false (text differs and nothing numeric to reconcile).
"""
function line_equivalence(ref::AbstractString, trial::AbstractString;
                           rel_tol::Float64=1e-9, abs_tol::Float64=1e-10)
    ref == trial && return true

    ref_matches   = collect(eachmatch(FLOAT_RE, ref))
    trial_matches = collect(eachmatch(FLOAT_RE, trial))

    # No floats → text-only comparison already failed.
    isempty(ref_matches) && return false

    length(ref_matches) != length(trial_matches) && return false

    for (rm, tm) in zip(ref_matches, trial_matches)
        rv = _parse_endf_or_normal(rm.match)
        tv = _parse_endf_or_normal(tm.match)
        (rv === nothing || tv === nothing) && return false
        isapprox(rv, tv; rtol=rel_tol, atol=abs_tol) || return false
    end

    # Non-numeric residue must match.
    ref_residue   = replace(ref,   FLOAT_RE => "")
    trial_residue = replace(trial, FLOAT_RE => "")
    ref_residue == trial_residue
end

"""
    tape_compare(ref_path, trial_path; rel_tol=1e-9, abs_tol=1e-10, max_show=5)

Compare two ENDF-format tape files line-by-line using execute.py semantics.
Returns a NamedTuple with per-line stats and up to `max_show` example diffs.
"""
function tape_compare(ref_path::AbstractString, trial_path::AbstractString;
                       rel_tol::Float64=1e-9, abs_tol::Float64=1e-10,
                       max_show::Int=5)
    ref_lines   = readlines(ref_path)
    trial_lines = readlines(trial_path)

    # Normalize date strings (execute.py does this prior to comparison)
    ref_lines   = [replace(l, DATE_RE => "XX/XX/XX") for l in ref_lines]
    trial_lines = [replace(l, DATE_RE => "XX/XX/XX") for l in trial_lines]

    n_ref   = length(ref_lines)
    n_trial = length(trial_lines)
    n_cmp   = min(n_ref, n_trial)
    line_count_match = (n_ref == n_trial)

    n_pass = 0
    diffs  = Tuple{Int,String,String}[]
    for i in 1:n_cmp
        if line_equivalence(ref_lines[i], trial_lines[i]; rel_tol=rel_tol, abs_tol=abs_tol)
            n_pass += 1
        elseif length(diffs) < max_show
            push!(diffs, (i, ref_lines[i], trial_lines[i]))
        end
    end

    (; n_ref, n_trial, n_cmp, n_pass, line_count_match,
       n_diff = n_cmp - n_pass,
       diffs,
       all_pass = line_count_match && n_pass == n_cmp)
end

# =========================================================================
# Heartbeat: prints STILL RUNNING when run_njoy has no module transition
# =========================================================================

"""
    with_heartbeat(f, progress; test_id="", heartbeat_sec=10.0, timeout_warn_sec=600.0) -> f()

Run `f()` while a background Timer watches `progress::RunProgress`. Every
`heartbeat_sec` seconds, if the current module has not changed (i.e. it's been
running for more than `heartbeat_sec`), print a heartbeat line.

At `timeout_warn_sec` the heartbeat gets loud ("⚠ SOFT TIMEOUT EXCEEDED").
Does NOT kill — Ctrl-C to stop. Cache safety requires in-process execution.
"""
function with_heartbeat(f, progress::NJOY.RunProgress;
                         test_id::AbstractString="",
                         heartbeat_sec::Float64=10.0,
                         timeout_warn_sec::Float64=600.0)
    run_start = time()
    last_heartbeat = Ref(run_start)
    warned_timeout = Ref(false)

    hb = Timer(heartbeat_sec; interval=heartbeat_sec) do _
        now = time()
        elapsed_total  = now - run_start
        elapsed_module = now - progress.module_start
        # Only heartbeat if a module is running AND has been >= heartbeat_sec
        if progress.status == :running && elapsed_module >= heartbeat_sec
            prefix = isempty(test_id) ? "" : "[$test_id] "
            if elapsed_total >= timeout_warn_sec && !warned_timeout[]
                @printf("%s⚠ SOFT TIMEOUT at %.0fs — still in %s (%.1fs)  last: %s\n",
                        prefix, elapsed_total, progress.current_module,
                        elapsed_module, progress.last_message)
                warned_timeout[] = true
            else
                @printf("%s… still in %s (%.1fs, total %.1fs)  last: %s\n",
                        prefix, progress.current_module,
                        elapsed_module, elapsed_total, progress.last_message)
            end
            flush(stdout)
            last_heartbeat[] = now
        end
    end

    try
        return f()
    finally
        close(hb)
    end
end

# =========================================================================
# Runner
# =========================================================================

"""
    run_reference_test(n::Int; tolerances, work_dir, verbose, heartbeat_sec, timeout_warn_sec)

Run test number `n` via `run_njoy(njoy-reference/tests/NN/input)` and compare
every produced `tape{U}` against every `referenceTape{U}` in the test dir.

Returns a NamedTuple:
    (; test, input, work_dir, run_ok, run_error, elapsed,
       tape_results::Vector{NamedTuple(unit, ref_path, trial_path, exists, compare, tolerance_pass)},
       summary::String, all_pass::Bool)

`tolerances` controls which rel_tol levels to report per tape. `all_pass`
reflects the STRICTEST tolerance (first entry).
"""
function run_reference_test(n::Int;
                             tolerances::Vector{Float64} = [1e-9, 1e-7, 1e-5],
                             work_dir::Union{Nothing,String} = nothing,
                             verbose::Bool = true,
                             heartbeat_sec::Float64 = 10.0,
                             timeout_warn_sec::Float64 = 600.0,
                             max_show::Int = 3)
    test_id   = @sprintf("T%02d", n)
    test_dir  = joinpath(TESTS_ROOT, lpad(n, 2, '0'))
    input_pth = joinpath(test_dir, "input")

    if !isfile(input_pth)
        return (; test=n, input=input_pth, work_dir=nothing, run_ok=false,
                  run_error="no input file", elapsed=0.0,
                  tape_results=NamedTuple[], summary="NO_INPUT", all_pass=false)
    end

    wdir = work_dir === nothing ? mktempdir(prefix="njoy_$(test_id)_") : work_dir
    mkpath(wdir)

    verbose && (@printf("╔ %s start  deck=%s  work_dir=%s\n", test_id, input_pth, wdir); flush(stdout))

    progress = NJOY.RunProgress()
    progress.run_start = time()

    run_ok    = false
    run_error = ""
    t0        = time()
    try
        with_heartbeat(progress; test_id=test_id,
                       heartbeat_sec=heartbeat_sec,
                       timeout_warn_sec=timeout_warn_sec) do
            NJOY.run_njoy(input_pth; work_dir=wdir, verbose=verbose, progress=progress)
        end
        run_ok = true
    catch ex
        run_error = first(split(sprint(showerror, ex), '\n'))
        verbose && (@printf("║ %s CRASH: %s\n", test_id, run_error); flush(stdout))
    end
    elapsed = time() - t0

    # Discover reference tapes (referenceTape{U} → unit U)
    ref_files = Tuple{Int,String}[]
    if isdir(test_dir)
        for f in sort(readdir(test_dir))
            m = match(r"^referenceTape(\d+)$", f)
            m === nothing && continue
            push!(ref_files, (parse(Int, m.captures[1]), joinpath(test_dir, f)))
        end
    end

    tape_results = NamedTuple[]
    for (unit, ref_path) in ref_files
        trial_path = joinpath(wdir, "tape$unit")
        if !isfile(trial_path)
            push!(tape_results, (;
                unit, ref_path, trial_path, exists=false,
                compare=nothing,
                tolerance_pass=Dict{Float64,Bool}(rtol => false for rtol in tolerances),
                status=:MISSING))
            verbose && (@printf("║ %s tape%d  MISSING (not produced)\n", test_id, unit); flush(stdout))
            continue
        end

        # Compare at each tolerance (strictest first)
        cmp_strict = tape_compare(ref_path, trial_path;
                                   rel_tol=first(tolerances), abs_tol=1e-10,
                                   max_show=max_show)
        tol_pass = Dict{Float64,Bool}()
        for rtol in tolerances
            if rtol == first(tolerances)
                tol_pass[rtol] = cmp_strict.all_pass
            else
                c = tape_compare(ref_path, trial_path; rel_tol=rtol, abs_tol=1e-10, max_show=0)
                tol_pass[rtol] = c.all_pass
            end
        end

        status = cmp_strict.all_pass ? :BIT_IDENTICAL :
                 !cmp_strict.line_count_match ? :STRUCTURAL_FAIL :
                 any(last, sort(collect(tol_pass), by=first)) ? :NUMERIC_PASS :
                 :DIFFS

        push!(tape_results, (;
            unit, ref_path, trial_path, exists=true,
            compare=cmp_strict, tolerance_pass=tol_pass, status))

        if verbose
            passed_tols = [rtol for rtol in tolerances if tol_pass[rtol]]
            tol_str = isempty(passed_tols) ? "no-tolerance-passes" :
                      "passes at " * join([@sprintf("%.0e", t) for t in passed_tols], ", ")
            @printf("║ %s tape%d  %s   %d/%d lines (%s)",
                    test_id, unit, String(status),
                    cmp_strict.n_pass, cmp_strict.n_cmp, tol_str)
            if !cmp_strict.line_count_match
                @printf("   [LINE COUNT: %d vs %d]", cmp_strict.n_ref, cmp_strict.n_trial)
            end
            println()
            for (ln, r, t) in cmp_strict.diffs
                @printf("║   first-diff line %d:\n║     R: %s\n║     J: %s\n", ln,
                        length(r) > 72 ? r[1:72] * "…" : r,
                        length(t) > 72 ? t[1:72] * "…" : t)
            end
            flush(stdout)
        end
    end

    # Overall summary
    all_pass_strict = run_ok && !isempty(tape_results) &&
                      all(tr -> get(tr.tolerance_pass, first(tolerances), false), tape_results)
    summary = if !run_ok
        "CRASH: $run_error"
    elseif isempty(ref_files)
        "NO_REFERENCE"
    elseif all_pass_strict
        @sprintf("ALL PASS @ rtol=%.0e", first(tolerances))
    else
        n_ok = count(tr -> get(tr.tolerance_pass, first(tolerances), false), tape_results)
        @sprintf("%d/%d tapes pass @ rtol=%.0e", n_ok, length(tape_results), first(tolerances))
    end

    verbose && (@printf("╚ %s done  %.1fs  %s\n\n", test_id, elapsed, summary); flush(stdout))

    (; test=n, input=input_pth, work_dir=wdir, run_ok, run_error, elapsed,
       tape_results, summary, all_pass=all_pass_strict)
end

# Script entry: `julia --project=. test/validation/reference_test.jl NN [NN ...]`
if abspath(PROGRAM_FILE) == @__FILE__
    nums = isempty(ARGS) ? [3, 4] : parse.(Int, ARGS)
    results = [run_reference_test(n) for n in nums]
    println("\n========================================")
    println("Summary")
    println("========================================")
    for r in results
        println(@sprintf("  T%02d  %-40s  (%.1fs)", r.test, r.summary, r.elapsed))
    end
end
