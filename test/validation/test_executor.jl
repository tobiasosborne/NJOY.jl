# test_executor.jl -- Execute parsed NJOY test cases through NJOY.jl
#
# Maps NJOY module names to Julia functions, chains outputs, handles
# graceful degradation for modules we cannot fully run.

using NJOY
using Printf

include("input_parser.jl")

const NJOY_REF  = normpath(joinpath(@__DIR__, "..", "..", "njoy-reference"))
const RESOURCES = joinpath(NJOY_REF, "tests", "resources")
const TESTS_DIR = joinpath(NJOY_REF, "tests")

# =========================================================================
# Types
# =========================================================================

struct NJOYTestCase
    number::Int
    modules::Vector{Symbol}
    runnable::Vector{Symbol}
    missing_mods::Vector{Symbol}
    endf_files::Dict{Int, String}
    reference_tapes::Dict{String, String}
    deck::NJOYInputDeck
    mat::Int; err::Float64; status::Symbol
end

struct StepResult
    module_name::Symbol; status::Symbol; message::String; elapsed::Float64
end

struct ExecutionResult
    test_number::Int; status::Symbol; message::String
    steps::Vector{StepResult}
    reconr_output::Any; broadr_output::Any; elapsed::Float64
end

# =========================================================================
# Test case builder
# =========================================================================

function parse_njoy_test(test_number::Int)
    dir = joinpath(TESTS_DIR, lpad(string(test_number), 2, '0'))
    input_file = joinpath(dir, "input")
    cmake_file = joinpath(dir, "CMakeLists.txt")

    deck = isfile(input_file) ? parse_njoy_input(input_file) : NJOYInputDeck(ModuleCall[])
    mods = modules_in_deck(deck)
    run = runnable_modules(deck); miss = missing_modules(deck)

    endf_files = Dict{Int,String}(); ref_tapes = Dict{String,String}()
    if isfile(cmake_file)
        ct = read(cmake_file, String)
        for (unit, fname) in parse_cmake_tapes(ct)
            fp = joinpath(RESOURCES, fname)
            !isfile(fp) && (fp = joinpath(dir, fname))
            isfile(fp) && (endf_files[unit] = fp)
        end
        for rn in parse_cmake_references(ct)
            rp = joinpath(dir, rn); isfile(rp) && (ref_tapes[rn] = rp)
        end
    end
    if isdir(dir)
        for f in readdir(dir)
            startswith(f, "referenceTape") && !haskey(ref_tapes, f) &&
                (ref_tapes[f] = joinpath(dir, f))
        end
    end

    mat = 0; err = 0.001
    for c in deck.calls
        if c.name == :reconr && length(c.raw_cards) >= 3
            rp = parse_reconr(c); mat = rp.mat; err = rp.err; break
        end
    end

    status = isempty(mods) ? :empty :
             (isempty(miss) || all(m -> m in SKIP_MODULES, miss)) ? :ready : :partial

    NJOYTestCase(test_number, mods, run, miss, endf_files, ref_tapes,
                 deck, mat, err, status)
end

parse_njoy_test(dir::AbstractString) = parse_njoy_test(parse(Int, basename(dir)))

# =========================================================================
# Module execution
# =========================================================================

function find_endf_file(tc::NJOYTestCase)
    for u in [20, 21, 22, 23, 30]
        haskey(tc.endf_files, u) && return tc.endf_files[u]
    end
    for (_, p) in sort(collect(tc.endf_files)); isfile(p) && return p; end
    nothing
end

function execute_reconr(endf_file::String, mc::ModuleCall;
                        mat_override::Int=0, err_override::Float64=0.0)
    t0 = time()
    p = parse_reconr(mc)
    mat = mat_override > 0 ? mat_override : p.mat
    err = err_override > 0 ? err_override : p.err
    mat <= 0 && return nothing, StepResult(:reconr, :error, "MAT not determined", time()-t0)
    try
        r = reconr(endf_file; mat=mat, err=err)
        return r, StepResult(:reconr, :ok, "MAT=$mat pts=$(length(r.energies))", time()-t0)
    catch ex
        msg = first(split(sprint(showerror, ex), '\n'))
        return nothing, StepResult(:reconr, :error, "reconr: $msg", time()-t0)
    end
end

function execute_broadr(reconr_result, mc::ModuleCall)
    t0 = time()
    p = parse_broadr(mc)
    if isempty(p.temperatures) || all(t -> t == 0.0, p.temperatures)
        return reconr_result, StepResult(:broadr, :ok, "T=0 pass-through", time()-t0)
    end
    temp = p.temperatures[1]
    try
        xs = hcat(reconr_result.total, reconr_result.elastic,
                  reconr_result.fission, reconr_result.capture)
        pendf = PointwiseMaterial(Int32(max(p.mat, 0)), copy(reconr_result.energies),
                                 xs, [1, 2, 18, 102])
        thnmax = p.thnmax < 0 ? abs(p.thnmax) : (p.thnmax > 0 ? p.thnmax : Inf)
        br = doppler_broaden(pendf, temp; T_old=0.0, awr=1.0,
                             tol=p.tol, thnmax=Float64(thnmax))
        return br, StepResult(:broadr, :ok, "T=$(temp)K pts=$(length(br.energies))", time()-t0)
    catch ex
        msg = first(split(sprint(showerror, ex), '\n'))
        return nothing, StepResult(:broadr, :error, "broadr: $msg", time()-t0)
    end
end

skip_step(name::Symbol, reason::String="not implemented") =
    StepResult(name, :skipped, reason, 0.0)

# =========================================================================
# Full test execution
# =========================================================================

function execute_test(tc::NJOYTestCase)
    t0 = time()
    steps = StepResult[]
    endf_file = find_endf_file(tc)
    endf_file === nothing && return ExecutionResult(tc.number, :skip,
        "No ENDF file", steps, nothing, nothing, time()-t0)
    isempty(tc.deck.calls) && return ExecutionResult(tc.number, :skip,
        "Empty deck", steps, nothing, nothing, time()-t0)
    tc.mat == 0 && any(c -> c.name == :reconr, tc.deck.calls) &&
        return ExecutionResult(tc.number, :skip, "No MAT", steps, nothing, nothing, time()-t0)

    reconr_out = nothing; broadr_out = nothing
    had_error = false; has_reconr = false

    for mc in tc.deck.calls
        if mc.name == :moder
            push!(steps, StepResult(:moder, :ok, "pass-through", 0.0))
        elseif mc.name == :reconr
            has_reconr = true
            r, s = execute_reconr(endf_file, mc; mat_override=tc.mat, err_override=tc.err)
            push!(steps, s); s.status == :ok ? (reconr_out = r) : (had_error = true)
        elseif mc.name == :broadr
            if reconr_out === nothing
                push!(steps, skip_step(:broadr, "no reconr output"))
            else
                r, s = execute_broadr(reconr_out, mc)
                push!(steps, s); s.status == :ok ? (broadr_out = r) : (had_error = true)
            end
        elseif mc.name in SKIP_MODULES
            push!(steps, skip_step(mc.name, "visualization"))
        else
            push!(steps, skip_step(mc.name, "skipped in pipeline"))
        end
    end

    status = had_error ? :error : (!has_reconr ? :skip :
             (reconr_out !== nothing ? :pass : :error))
    n_ok = count(s -> s.status == :ok, steps)
    n_skip = count(s -> s.status == :skipped, steps)
    msg = status == :skip ? "No reconr [$(join(unique(tc.modules), ","))]" :
          status == :pass ? "ran $n_ok, skipped $n_skip" : "errors occurred"
    ExecutionResult(tc.number, status, msg, steps, reconr_out, broadr_out, time()-t0)
end

# Format a single execution result for display.
function format_execution(er::ExecutionResult)
    buf = IOBuffer()
    @printf(buf, "%-6s Test %02d  (%.1fs)  %s\n",
            uppercase(string(er.status)), er.test_number, er.elapsed, er.message)
    for s in er.steps
        @printf(buf, "       %-8s %-6s  %s\n",
                s.module_name, uppercase(string(s.status)), s.message)
    end
    String(take!(buf))
end
