# njoy_test_runner.jl -- Legacy runner, delegates to the new pipeline files.
#
# This file is kept for backward compatibility. New code should use:
#   include("run_all.jl") for the full pipeline
#   include("test_executor.jl") for execution only
#   include("reference_comparator.jl") for comparison only
#   include("input_parser.jl") for parsing only

using NJOY
using Printf

include("test_executor.jl")
include("reference_comparator.jl")

# Re-export key types and functions for backward compatibility
const NJOY_REF_COMPAT  = NJOY_REF
const RESOURCES_COMPAT = RESOURCES
const TESTS_DIR_COMPAT = TESTS_DIR

# Compatibility: run_njoy_test returns a TestResult-like structure
struct TestResult
    number::Int; status::Symbol; message::String
    comparisons::Vector{ReactionComparison}; elapsed::Float64
end

function run_njoy_test(tc::NJOYTestCase; rtol::Float64=0.05)
    exec = execute_test(tc)
    comps = ReactionComparison[]
    ref_tape = ""
    if exec.pendf !== nothing && !isempty(tc.reference_tapes)
        report = compare_best_ref(exec.pendf, tc.reference_tapes, tc.number)
        comps = report.reactions
        ref_tape = report.ref_tape
    end
    if isempty(comps)
        return TestResult(tc.number, exec.status == :skip ? :skip : :error,
                          exec.message, comps, exec.elapsed)
    end
    all_pass = all(c -> c.filt_rel_err <= rtol, comps)
    worst = maximum(c -> c.filt_rel_err, comps)
    msg = @sprintf("max_err=%.2f%% ref=%s", worst * 100, ref_tape)
    TestResult(tc.number, all_pass ? :pass : :fail, msg, comps, exec.elapsed)
end

# Format a TestResult
function format_result(r::TestResult)
    buf = IOBuffer()
    @printf(buf, "%-6s Test %02d  (%.1fs)  %s\n",
            uppercase(string(r.status)), r.number, r.elapsed, r.message)
    for c in r.comparisons
        fl = c.filt_rel_err <= 0.01 ? "OK" : (c.filt_rel_err <= 0.10 ? "WARN" : "FAIL")
        @printf(buf, "       MT=%-3d %-8s  ref=%5d  ours=%6d  filt=%.3f%%  [%s]\n",
                c.mt, c.label, c.n_ref, c.n_ours, c.filt_rel_err*100, fl)
    end
    String(take!(buf))
end
