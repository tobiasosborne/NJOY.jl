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
    mat::Int; err::Float64; awr::Float64; status::Symbol
end

struct StepResult
    module_name::Symbol; status::Symbol; message::String; elapsed::Float64
end

struct ExecutionResult
    test_number::Int; status::Symbol; message::String
    steps::Vector{StepResult}
    pendf::Any; elapsed::Float64
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

    mat = 0; err = 0.001; awr = 1.0
    for c in deck.calls
        if c.name == :reconr && length(c.raw_cards) >= 3
            rp = parse_reconr(c); mat = rp.mat; err = rp.err; break
        end
    end

    # Extract AWR from the ENDF file's MF1/MT451 HEAD record
    if mat > 0
        for u in [20, 21, 22, 23, 30]
            haskey(endf_files, u) || continue
            try
                open(endf_files[u]) do io
                    while !eof(io)
                        line = readline(io)
                        p = rpad(line, 80)
                        mf = NJOY._parse_int(p[71:72])
                        mt = NJOY._parse_int(p[73:75])
                        if mf == 1 && mt == 451
                            awr = parse_endf_float(p[12:22])
                            break
                        end
                    end
                end
                awr > 0 && break
            catch; end
        end
    end

    status = isempty(mods) ? :empty :
             (isempty(miss) || all(m -> m in SKIP_MODULES, miss)) ? :ready : :partial

    NJOYTestCase(test_number, mods, run, miss, endf_files, ref_tapes,
                 deck, mat, err, awr, status)
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

"""Convert reconr NamedTuple to PointwiseMaterial for downstream modules."""
function reconr_to_pendf(r, mat::Int)
    xs = hcat(r.total, r.elastic, r.fission, r.capture)
    PointwiseMaterial(Int32(mat), copy(r.energies), xs, [1, 2, 18, 102])
end

"""Convert PointwiseMaterial to a NamedTuple the comparator can consume."""
function pendf_to_result(pendf::PointwiseMaterial)
    mt_to_col = Dict(zip(pendf.mt_list, axes(pendf.cross_sections, 2)))
    col(mt) = haskey(mt_to_col, mt) ? pendf.cross_sections[:, mt_to_col[mt]] :
                                       zeros(length(pendf.energies))
    (energies=pendf.energies, total=col(1), elastic=col(2),
     fission=col(18), capture=col(102))
end

"""Ensure pipeline state is a PointwiseMaterial."""
ensure_pendf(p::PointwiseMaterial, ::Int) = p
ensure_pendf(r, mat::Int) = reconr_to_pendf(r, mat)

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

function execute_broadr(pendf::PointwiseMaterial, mc::ModuleCall; awr::Float64=1.0)
    t0 = time()
    p = parse_broadr(mc)
    if isempty(p.temperatures) || all(t -> t == 0.0, p.temperatures)
        return pendf, StepResult(:broadr, :ok, "T=0 pass-through", time()-t0)
    end
    temp = p.temperatures[1]
    try
        thnmax = p.thnmax < 0 ? abs(p.thnmax) : (p.thnmax > 0 ? p.thnmax : Inf)
        br = doppler_broaden(pendf, temp; T_old=0.0, awr=awr,
                             tol=p.tol, thnmax=Float64(thnmax))
        return br, StepResult(:broadr, :ok,
            @sprintf("T=%.0fK AWR=%.1f pts=%d", temp, awr, length(br.energies)), time()-t0)
    catch ex
        msg = first(split(sprint(showerror, ex), '\n'))
        return nothing, StepResult(:broadr, :error, "broadr: $msg", time()-t0)
    end
end

function execute_heatr(pendf::PointwiseMaterial, mc::ModuleCall; awr::Float64=1.0)
    t0 = time()
    p = parse_heatr(mc)
    try
        kr = compute_kerma(pendf; awr=awr, Z=0)
        msg = @sprintf("KERMA %d pts, peak=%.2e", length(kr.energies), maximum(abs, kr.total_kerma))
        return pendf, StepResult(:heatr, :ok, msg, time()-t0)
    catch ex
        msg = first(split(sprint(showerror, ex), '\n'))
        return pendf, StepResult(:heatr, :warn, "heatr: $msg (pendf unchanged)", time()-t0)
    end
end

function execute_thermr(pendf::PointwiseMaterial, mc::ModuleCall)
    t0 = time()
    p = parse_thermr(mc)
    temp = isempty(p.temperatures) ? 296.0 : p.temperatures[1]
    try
        result = compute_thermal(pendf, temp, 1.0; emax=p.emax)
        return result, StepResult(:thermr, :ok, "T=$(temp)K", time()-t0)
    catch ex
        msg = first(split(sprint(showerror, ex), '\n'))
        return pendf, StepResult(:thermr, :warn, "thermr: $msg (pendf unchanged)", time()-t0)
    end
end

function execute_gaspr(pendf::PointwiseMaterial, mc::ModuleCall)
    t0 = time()
    try
        gr = compute_gas_production(pendf)
        n = length(gr.energies)
        return pendf, StepResult(:gaspr, :ok, "$n pts", time()-t0)
    catch ex
        msg = first(split(sprint(showerror, ex), '\n'))
        return pendf, StepResult(:gaspr, :warn, "gaspr: $msg (pendf unchanged)", time()-t0)
    end
end

function execute_acer(pendf::PointwiseMaterial, mc::ModuleCall)
    t0 = time()
    p = parse_acer(mc)
    p.iopt != 1 && return pendf, StepResult(:acer, :skipped, "iopt=$(p.iopt)", time()-t0)
    try
        ace = build_ace_from_pendf(pendf; suffix=p.suffix, temp_kelvin=p.temp)
        return pendf, StepResult(:acer, :ok, "NES=$(length(pendf.energies))", time()-t0)
    catch ex
        msg = first(split(sprint(showerror, ex), '\n'))
        return pendf, StepResult(:acer, :warn, "acer: $msg", time()-t0)
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
        "No ENDF file", steps, nothing, time()-t0)
    isempty(tc.deck.calls) && return ExecutionResult(tc.number, :skip,
        "Empty deck", steps, nothing, time()-t0)
    tc.mat == 0 && any(c -> c.name == :reconr, tc.deck.calls) &&
        return ExecutionResult(tc.number, :skip, "No MAT", steps, nothing, time()-t0)

    pendf = nothing        # latest pipeline state (NamedTuple or PointwiseMaterial)
    had_error = false; has_reconr = false

    for mc in tc.deck.calls
        if mc.name == :moder
            push!(steps, StepResult(:moder, :ok, "pass-through", 0.0))

        elseif mc.name == :reconr
            has_reconr = true
            r, s = execute_reconr(endf_file, mc; mat_override=tc.mat, err_override=tc.err)
            push!(steps, s)
            s.status == :ok ? (pendf = r) : (had_error = true)

        elseif mc.name in (:broadr, :heatr, :thermr, :gaspr, :acer,
                           :unresr, :purr, :groupr, :errorr, :covr)
            if pendf === nothing
                push!(steps, skip_step(mc.name, "no upstream pendf"))
            else
                r, s = _execute_downstream(mc, pendf, tc.mat; awr=tc.awr)
                push!(steps, s)
                r !== nothing && (pendf = r)
                s.status == :error && (had_error = true)
            end

        elseif mc.name in SKIP_MODULES
            push!(steps, skip_step(mc.name, "visualization"))
        else
            push!(steps, skip_step(mc.name, "unknown module"))
        end
    end

    # For comparison, convert to NamedTuple if needed
    result = pendf isa PointwiseMaterial ? pendf_to_result(pendf) : pendf

    status = had_error ? :error : (!has_reconr ? :skip :
             (pendf !== nothing ? :pass : :error))
    n_ok = count(s -> s.status == :ok, steps)
    n_skip = count(s -> s.status == :skipped, steps)
    msg = status == :skip ? "No reconr [$(join(unique(tc.modules), ","))]" :
          status == :pass ? "ran $n_ok, skipped $n_skip" : "errors occurred"
    ExecutionResult(tc.number, status, msg, steps, result, time()-t0)
end

"""Dispatch a downstream module call on the current pipeline state."""
function _execute_downstream(mc::ModuleCall, pendf, mat::Int; awr::Float64=1.0)
    pm = ensure_pendf(pendf, mat)
    if mc.name == :broadr;  return execute_broadr(pm, mc; awr=awr)
    elseif mc.name == :heatr;   return execute_heatr(pm, mc; awr=awr)
    elseif mc.name == :thermr;  return execute_thermr(pm, mc)
    elseif mc.name == :gaspr;   return execute_gaspr(pm, mc)
    elseif mc.name == :acer;    return execute_acer(pm, mc)
    elseif mc.name == :unresr;  return pm, skip_step(:unresr, "self-shielding (pass-through)")
    elseif mc.name == :purr;    return pm, skip_step(:purr, "prob tables (pass-through)")
    elseif mc.name == :groupr;  return pm, skip_step(:groupr, "multigroup (pass-through)")
    elseif mc.name == :errorr;  return pm, skip_step(:errorr, "covariance (pass-through)")
    elseif mc.name == :covr;    return pm, skip_step(:covr, "cov output (pass-through)")
    else;                       return pm, skip_step(mc.name)
    end
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
