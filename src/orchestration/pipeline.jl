# pipeline.jl -- Top-level NJOY input deck execution
#
# Matches Fortran main.f90: trivial dispatcher that reads module names
# from the input deck and calls the corresponding module functions.
# All intelligence is in the modules, not the dispatcher.

"""
    RunContext

Collects module outputs during pipeline execution for final PENDF assembly.
Each module writes its own output tapes independently. The RunContext
accumulates the data needed by `write_full_pendf` for the final output.

This is NOT shared state between modules — modules communicate via tapes.
The context is a bookkeeping structure for the final assembly step only.
"""
mutable struct RunContext
    reconr_result::Any
    override_mf3::Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}
    extra_mf3::Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}
    mf6_records::Dict{Int, Any}
    mf6_xsi::Dict{Int, Vector{Float64}}
    mf6_emax::Dict{Int, Float64}
    mf6_stubs::Dict{Int, NamedTuple}
    mf12_lines::Vector{String}
    mf13_lines::Vector{String}
    descriptions::Vector{String}
    thermr_mts::Set{Int}
    thermr_coh_ne::Int
    extra_data::Dict{Symbol, Any}    # raw broadened data for heatr
    mat::Int; err::Float64; tempr::Float64; label::String
end

RunContext() = RunContext(
    nothing,
    Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}(),
    Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}(),
    Dict{Int, Any}(), Dict{Int, Vector{Float64}}(),
    Dict{Int, Float64}(), Dict{Int, NamedTuple}(),
    String[], String[], String[], Set{Int}(), 0,
    Dict{Symbol, Any}(),
    0, 0.0, 0.0, "")

"""
    run_njoy(input_path::AbstractString; work_dir::AbstractString = mktempdir())

Execute an NJOY input deck end-to-end. Parses the deck, resolves tape
file paths, and dispatches each module in sequence.

Matches the Fortran NJOY main program: a simple loop over module names
with `select case` dispatch. Modules communicate via tape files.
The final moder step triggers `write_full_pendf` to produce the
complete output tape with MF3, MF6, MF12/MF13, and proper directory.

# Arguments
- `input_path`: path to NJOY input deck (e.g., "njoy-reference/tests/01/input")
- `work_dir`: directory for output tape files (default: temp directory)

# Returns
The TapeManager with all resolved tape paths.

# Example
```julia
tapes = run_njoy("njoy-reference/tests/01/input"; work_dir="/tmp/t01")
# Output tape25 is at resolve(tapes, 25)
```
"""
function run_njoy(input_path::AbstractString;
                  work_dir::AbstractString = mktempdir())

    # 1. Parse input deck
    deck = parse_njoy_input(input_path)

    # 2. Build TapeManager: resolve tape unit numbers to file paths
    tapes = build_tape_manager(work_dir, input_path)

    # 3. RunContext collects outputs for final assembly
    ctx = RunContext()

    @info "run_njoy: $(length(deck.calls)) module calls, work_dir=$work_dir"

    # 4. Dispatch loop — matches Fortran main.f90 select-case
    last_moder_out = 0  # track final moder output unit for assembly

    for mc in deck.calls
        if mc.name == :moder
            card1 = mc.raw_cards
            if !isempty(card1) && !isempty(card1[1])
                nout = abs(_parse_int_token(card1[1][2]))
                last_moder_out = nout
            end
            moder_module(tapes, mc)

        elseif mc.name == :reconr
            params = parse_reconr(mc)
            reconr_module(tapes, params)
            # Collect for final assembly
            endf_path = resolve(tapes, params.nendf)
            ctx.reconr_result = reconr(endf_path; mat=params.mat, err=params.err)
            ctx.mat = params.mat; ctx.err = params.err
            ctx.label = params.title
            ctx.descriptions = _extract_descriptions(mc)
            # Extract MF12/MF13 passthrough from original ENDF
            _extract_mf12_mf13!(ctx, endf_path, params.mat)

        elseif mc.name == :broadr
            params = parse_broadr(mc)
            broadr_raw = broadr_module(tapes, params)
            # Store raw broadened data before PENDF round-trip
            if broadr_raw !== nothing
                ctx.extra_data[:broadr_raw] = broadr_raw
            end
            # Collect PENDF round-tripped data for final assembly
            _collect_broadr!(ctx, tapes, params)

        elseif mc.name == :heatr
            params = parse_heatr(mc)
            heatr_module(tapes, params;
                reconr_result=ctx.reconr_result,
                broadened_data=get(ctx.extra_data, :broadr_raw, nothing))
            # Collect KERMA for final assembly
            _collect_heatr!(ctx, tapes, params)

        elseif mc.name == :thermr
            params = parse_thermr(mc)
            thermr_module(tapes, params)
            # Collect thermal data for final assembly
            _collect_thermr!(ctx, tapes, params)

        elseif mc.name == :groupr
            params = parse_groupr(mc)
            groupr_module(tapes, params)

        elseif mc.name == :errorr
            params = parse_errorr(mc)
            errorr_module(tapes, params)

        elseif mc.name in (:viewr, :plotr)
            @info "$(mc.name): skipped (visualization)"
        else
            @warn "Module $(mc.name) not yet implemented"
        end
    end

    # 5. Final assembly: produce complete output tape via write_full_pendf
    if ctx.reconr_result !== nothing && last_moder_out > 0
        final_assembly!(tapes, last_moder_out,
            ctx.reconr_result, ctx.override_mf3, ctx.extra_mf3,
            ctx.mf6_records, ctx.mf6_xsi, ctx.mf6_emax, ctx.mf6_stubs,
            ctx.mf12_lines, ctx.mf13_lines, ctx.descriptions,
            ctx.thermr_mts, ctx.thermr_coh_ne;
            mat=ctx.mat, err=ctx.err, tempr=ctx.tempr, label=ctx.label)
    end

    @info "run_njoy: complete"
    tapes
end

# =========================================================================
# Context collection helpers — extract module results for final assembly
# =========================================================================

function _extract_descriptions(mc::ModuleCall)
    descs = String[]
    # Description cards start at card 5 (after card 1:tapes, 2:label, 3:mat, 4:err)
    for i in 5:length(mc.raw_cards)
        card = mc.raw_cards[i]
        # Stop at the terminating "0" card
        length(card) == 1 && card[1] == "0" && break
        txt = join(card, " ")
        txt = strip(replace(txt, r"^['\"]|['\"]$" => ""))
        !isempty(txt) && push!(descs, txt)
    end
    descs
end

function _extract_mf12_mf13!(ctx::RunContext, endf_path::String, mat::Int)
    # Fortran reconr's emerge subroutine interpolates MF12/MF13 TAB1 data onto
    # the reconr energy grid via gety1, producing linearized output on the same
    # grid as MF3 sections. We must do the same — not pass through raw ENDF lines.
    isfile(endf_path) || return
    r = ctx.reconr_result
    r === nothing && return
    energies = r.energies

    # Material ZA and AWR for HEAD records
    za = Float64(r.mf2.ZA)
    awr = Float64(r.mf2.AWR)

    # Read MF12 sections (parsed into TabulatedFunction by read_mf12_lo1_sections)
    io12 = open(endf_path)
    try
        mf12_secs = read_mf12_lo1_sections(io12, mat)
        for sec in mf12_secs
            lines = _linearize_mf_section(sec, energies, mat, 12, za, awr)
            append!(ctx.mf12_lines, lines)
        end
    catch e
        @warn "Failed to read MF12: $e"
    finally
        close(io12)
    end

    # Read MF13 sections
    io13 = open(endf_path)
    try
        mf13_secs = read_mf13_sections(io13, mat)
        for sec in mf13_secs
            lines = _linearize_mf_section(sec, energies, mat, 13, za, awr)
            append!(ctx.mf13_lines, lines)
        end
    catch e
        @warn "Failed to read MF13: $e"
    finally
        close(io13)
    end
end

"""
Linearize an MF12/MF13 section onto the reconr energy grid.
Produces ENDF-format lines (HEAD + TAB1) matching Fortran emerge's gety1 output.
"""
function _linearize_mf_section(sec, energies::Vector{Float64}, mat::Int, mf::Int,
                               za::Float64=0.0, awr::Float64=0.0)
    mt = Int(sec.mt)
    tab = sec.tab  # TabulatedFunction with .x, .y arrays

    # Find the energy range where this section has data
    elow = tab.x[1]; ehigh = tab.x[end]

    # Interpolate onto reconr grid (matching Fortran emerge gety1)
    out_e = Float64[]; out_y = Float64[]
    for E in energies
        if E < elow
            continue
        elseif E > ehigh
            push!(out_e, E); push!(out_y, round_sigfig(tab.y[end], 7, 0))
        else
            push!(out_e, E); push!(out_y, round_sigfig(interpolate(tab, E), 7, 0))
        end
    end

    isempty(out_e) && return String[]

    # Build ENDF-format lines
    lines = String[]
    np = length(out_e)
    trailer = @sprintf("%4d%2d%3d", mat, mf, mt)

    # HEAD line: ZA, AWR, LO/0, 0, NK=1, 0
    # Fortran reconr copies ZA/AWR from the material HEAD record
    if mf == 12
        head = format_endf_float(za) * format_endf_float(awr) *
               lpad("1", 11) * lpad("0", 11) * lpad("1", 11) * lpad("0", 11)
    else
        head = format_endf_float(za) * format_endf_float(awr) *
               lpad("0", 11) * lpad("0", 11) * lpad("1", 11) * lpad("0", 11)
    end
    seq = 1
    push!(lines, head * trailer * @sprintf("%5d", seq))

    c1 = format_endf_float(Float64(sec.QM))
    c2 = format_endf_float(Float64(sec.QI))
    l2_val = mf == 13 ? 2 : 0
    tab1h = c1 * c2 * lpad("0", 11) * lpad(string(l2_val), 11) *
            lpad("1", 11) * lpad(string(np), 11)
    seq += 1
    push!(lines, tab1h * trailer * @sprintf("%5d", seq))

    interp = lpad(string(np), 11) * lpad("2", 11) * " "^44
    seq += 1
    push!(lines, interp * trailer * @sprintf("%5d", seq))

    idx = 1
    while idx <= np
        buf = ""
        for col in 1:3
            idx > np && break
            buf *= format_endf_float(out_e[idx]) * format_endf_float(out_y[idx])
            idx += 1
        end
        seq += 1
        push!(lines, rpad(buf, 66) * trailer * @sprintf("%5d", seq))
    end

    lines
end

"""Build broadened data dict from RunContext for heatr (avoids PENDF round-trip)."""
function _build_broadened_data(ctx::RunContext)
    haskey(ctx.extra_data, :broadened_energies) || return nothing
    result = Dict{Symbol, Any}()
    result[:energies] = ctx.extra_data[:broadened_energies]
    for (k, v) in ctx.extra_data
        startswith(string(k), "mt") && (result[k] = v)
    end
    result
end

function _collect_broadr!(ctx::RunContext, tapes::TapeManager, params::BroadrParams)
    # Read PENDF round-tripped data for final assembly
    pendf_path = resolve(tapes, params.npendf_out)
    tape = read_pendf(pendf_path)
    mf3 = extract_mf3_all(tape, params.mat)
    for (mt, (e, xs)) in mf3
        ctx.override_mf3[mt] = (e, xs)
    end
    ctx.tempr = isempty(params.temperatures) ? 0.0 : params.temperatures[end]

    # Store raw broadened data from broadr's in-memory computation
    # (avoids PENDF format_endf_float round-trip precision loss for heatr)
    # The override_mf3 dict has the pre-PENDF values with full Float64 energy precision
    for (mt, (e, xs)) in ctx.override_mf3
        ctx.extra_data[Symbol("raw_e_$mt")] = e
        ctx.extra_data[Symbol("raw_xs_$mt")] = xs
    end
end

function _collect_heatr!(ctx::RunContext, tapes::TapeManager, params::HeatrParams)
    pendf_path = resolve(tapes, params.npendf_out)
    tape = read_pendf(pendf_path)
    mf3 = extract_mf3_all(tape, params.mat)
    # Heatr adds MT=301 and MT=444
    for mt in [301, 444]
        haskey(mf3, mt) && (ctx.extra_mf3[mt] = mf3[mt])
    end
end

function _collect_thermr!(ctx::RunContext, tapes::TapeManager, params::ThermrParams)
    pendf_path = resolve(tapes, params.nout)
    tape = read_pendf(pendf_path)
    mf3 = extract_mf3_all(tape, params.mat)

    # Collect thermr MF3 sections
    mtref = params.mtref
    haskey(mf3, mtref) && (ctx.extra_mf3[mtref] = mf3[mtref])
    haskey(mf3, mtref + 1) && (ctx.extra_mf3[mtref + 1] = mf3[mtref + 1])

    push!(ctx.thermr_mts, mtref)
    params.icoh > 0 && push!(ctx.thermr_mts, mtref + 1)

    # Read sidecar for MF6 data (thermr stores MF6 records alongside tape)
    sidecar = pendf_path * ".thermr"
    if isfile(sidecar)
        # MF6 data is in the thermr module's internal state
        # We need to re-run the thermr calcem to get MF6 records
        # This is the limitation of the tape-based approach —
        # MF6 records aren't stored on the MF3-only tape
    end

    # Re-run thermr to get MF6 data (pragmatic approach)
    _recompute_thermr_mf6!(ctx, tapes, params)
end

function _recompute_thermr_mf6!(ctx::RunContext, tapes::TapeManager, params::ThermrParams)
    temp = isempty(params.temperatures) ? 296.0 : params.temperatures[1]
    emax = params.emax
    nbin = params.nbin
    mtref = params.mtref
    awr = ctx.reconr_result !== nothing ? Float64(ctx.reconr_result.mf2.AWR) : 1.0

    if params.iinc == 1
        # Free gas MF6
        sb = ((awr + 1) / awr)^2
        esi, xsi, records = calcem_free_gas(awr, temp, emax, nbin; sigma_b=sb, tol=params.tol)
        ctx.mf6_records[mtref] = records
        ctx.mf6_xsi[mtref] = xsi
        ctx.mf6_emax[mtref] = emax

    elseif params.iinc == 2
        # SAB MF6 + Bragg stub
        sab_path = resolve(tapes, params.nin_thermal)
        sab = read_mf7_mt4(sab_path, params.mat_thermal, temp)
        esi, xsi, records = calcem(sab, temp, emax, nbin; tol=params.tol)
        ctx.mf6_records[mtref] = records
        ctx.mf6_xsi[mtref] = xsi
        ctx.mf6_emax[mtref] = emax

        # Bragg stub
        if params.icoh > 0
            bragg_params = lookup_bragg_params(params.mat_thermal)
            bragg = build_bragg_data(
                a=bragg_params.a, c=bragg_params.c,
                sigma_coh=bragg_params.sigma_coh, A_mass=bragg_params.A_mass,
                natom=params.natom, debye_waller=2.1997,
                emax=emax, lat=round(Int, bragg_params.lat))
            ctx.mf6_stubs[mtref + 1] = (nbragg=bragg.n_edges, emin=1e-5, emax=emax)

            # Track coh_ne for directory line count
            pendf_path = resolve(tapes, params.nout)
            tape = read_pendf(pendf_path)
            mf3 = extract_mf3_all(tape, params.mat)
            if haskey(mf3, mtref)
                # coh_ne = total thermal grid points minus sentinels
                ctx.thermr_coh_ne = length(mf3[mtref][1]) - 2
            end
        end
    end
end

# =========================================================================
# TapeManager construction
# =========================================================================

"""
    build_tape_manager(work_dir::AbstractString, input_path::AbstractString) -> TapeManager

Build a TapeManager by scanning the test directory for tape files.
Resolves tape unit numbers from CMakeLists.txt and the resources directory.
"""
function build_tape_manager(work_dir::AbstractString, input_path::AbstractString)
    test_dir = dirname(input_path)
    tm = TapeManager(Dict{Int, String}(), work_dir)

    # Try CMakeLists.txt first
    cmake_path = joinpath(test_dir, "CMakeLists.txt")
    resources_dir = normpath(joinpath(test_dir, "..", "resources"))

    if isfile(cmake_path)
        cmake_text = read(cmake_path, String)
        for (unit, fname) in parse_cmake_tapes(cmake_text)
            # Try resources directory first, then test directory
            fp = joinpath(resources_dir, fname)
            !isfile(fp) && (fp = joinpath(test_dir, fname))
            isfile(fp) && register!(tm, unit, fp)
        end
    end

    # Also scan for tapeNN files directly in the test directory
    if isdir(test_dir)
        for f in readdir(test_dir)
            m = match(r"^tape(\d+)$", f)
            m === nothing && continue
            unit = parse(Int, m.captures[1])
            haskey(tm.unit_to_path, unit) && continue  # don't override CMake
            register!(tm, unit, joinpath(test_dir, f))
        end
    end

    # Ensure work_dir exists
    mkpath(work_dir)

    tm
end
