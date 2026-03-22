# fortran_oracle.jl — Run the real NJOY2016 Fortran binary with truncated decks
#
# Produces per-module intermediate reference PENDFs by running the Fortran
# binary with input decks truncated at each module boundary. Each truncated
# deck ends with `moder -XX YY` to convert the result to ASCII for comparison.
#
# Caches results in oracle_cache/ to avoid re-running Fortran on every invocation.

using Printf

include("input_parser.jl")

const NJOY_BINARY = joinpath(@__DIR__, "..", "..", "njoy-reference", "build", "njoy")
const ORACLE_CACHE = joinpath(@__DIR__, "oracle_cache")
const NJOY_REF_ORACLE = joinpath(@__DIR__, "..", "..", "njoy-reference")
const RESOURCES_ORACLE = joinpath(NJOY_REF_ORACLE, "tests", "resources")

# Modules that produce PENDF output worth comparing
const PENDF_MODULES = Set([:reconr, :broadr, :heatr, :thermr, :unresr, :purr, :gaspr])

# =========================================================================
# Tape unit tracking
# =========================================================================

"""
Track which tape unit holds the 'current PENDF' after each module step,
by inspecting the module's I/O card (first raw_card line with unit numbers).
Returns (input_unit, output_unit) where negative = binary.
"""
function _module_io_units(mc::ModuleCall)
    isempty(mc.raw_cards) && return (0, 0)
    tokens = mc.raw_cards[1]
    length(tokens) < 2 && return (0, 0)
    units = Int[]
    for t in tokens
        v = tryparse(Int, replace(t, "/" => ""))
        v !== nothing && push!(units, v)
    end
    if mc.name == :reconr
        # reconr: nendf npend → output is npend
        length(units) >= 2 && return (units[1], units[2])
    elseif mc.name == :broadr
        # broadr: nendf nin nout → output is nout
        length(units) >= 3 && return (units[2], units[3])
    elseif mc.name == :heatr
        # heatr: nendf nin nout → output is nout (or nin if nout=0)
        length(units) >= 3 && return (units[2], units[3] == 0 ? units[2] : units[3])
    elseif mc.name == :thermr
        # thermr: nin_thermal nin_pendf nout
        length(units) >= 3 && return (units[2], units[3])
    elseif mc.name in (:unresr, :purr)
        # unresr/purr: nendf nin nout
        length(units) >= 3 && return (units[2], units[3])
    elseif mc.name == :gaspr
        # gaspr: nendf nin nout
        length(units) >= 3 && return (units[2], units[3])
    elseif mc.name == :groupr
        # groupr: nendf npend ngout ngg
        length(units) >= 3 && return (units[2], units[3])
    elseif mc.name == :moder
        # moder: nin nout
        length(units) >= 2 && return (units[1], units[2])
    elseif mc.name == :acer
        # acer: nendf npend nace ndir
        length(units) >= 2 && return (units[2], 0)
    end
    (0, 0)
end

# =========================================================================
# Truncated deck generation
# =========================================================================

"""
    build_truncated_deck(input_file::String, deck::NJOYInputDeck, stop_after_idx::Int) -> (String, Int)

Build a truncated input deck by reading the raw input file and truncating at
the module boundary. Preserves exact NJOY formatting (slashes, quotes, etc.).
Appends a `moder` call to convert the last PENDF output to ASCII tape 90.
Returns (deck_text, output_tape_unit).
"""
function build_truncated_deck(input_file::String, deck::NJOYInputDeck,
                              stop_after_idx::Int)
    raw_lines = readlines(input_file)

    # Find the line in the raw file where module stop_after_idx+1 starts (or "stop")
    # by scanning for module keywords
    mod_idx = 0
    cut_line = length(raw_lines) + 1  # default: keep everything

    for (i, line) in enumerate(raw_lines)
        word = lowercase(strip(replace(line, r"\s*/.*" => "")))
        word = strip(replace(word, r"['\"]" => ""))
        word = replace(word, r"^\[x\d+\]\s*" => "")
        if Symbol(word) in NJOY_MODULES || word == "stop"
            mod_idx += 1
            if mod_idx > stop_after_idx + 1  # +1 because we want to include stop_after_idx
                # But we need to be more careful: stop_after_idx counts only
                # non-moder modules. Let me count ALL modules.
                break
            end
        end
    end

    # Actually, simpler: count module calls in order, matching deck.calls
    # The deck.calls includes ALL modules (moder, reconr, etc.)
    call_idx = 0
    cut_line = length(raw_lines) + 1

    for (i, line) in enumerate(raw_lines)
        word = lowercase(strip(replace(line, r"\s*/.*" => "")))
        word = strip(replace(word, r"['\"]" => ""))
        word = replace(word, r"^\[x\d+\]\s*" => "")
        if Symbol(word) in NJOY_MODULES
            call_idx += 1
            if call_idx > stop_after_idx
                cut_line = i
                break
            end
        elseif word == "stop"
            cut_line = i
            break
        end
    end

    # Keep lines up to (but not including) cut_line
    kept = raw_lines[1:min(cut_line - 1, length(raw_lines))]

    # Track the last PENDF output unit
    last_pendf_unit = 0
    for i in 1:min(stop_after_idx, length(deck.calls))
        mc = deck.calls[i]
        _, out = _module_io_units(mc)
        if mc.name in PENDF_MODULES && out != 0
            last_pendf_unit = out
        end
    end

    # Append moder to convert to ASCII
    if last_pendf_unit != 0
        push!(kept, "moder")
        push!(kept, " $(last_pendf_unit) 90/")
    end
    push!(kept, "stop")

    (join(kept, "\n") * "\n", 90)
end

# =========================================================================
# Running the Fortran binary
# =========================================================================

"""
    run_fortran_njoy(deck_text, tape_files; workdir, njoy_binary) -> String

Run the Fortran NJOY binary with the given input deck.
`tape_files` maps tape unit numbers to source file paths (will be symlinked).
Returns the path to the ASCII output PENDF (tape90), or empty string on failure.
"""
function run_fortran_njoy(deck_text::String, tape_files::Dict{Int,String};
                          workdir::String, njoy_binary::String=NJOY_BINARY)
    mkpath(workdir)

    # Write input deck
    input_path = joinpath(workdir, "input")
    write(input_path, deck_text)

    # Symlink tape files
    for (unit, src) in tape_files
        dst = joinpath(workdir, "tape$unit")
        isfile(dst) && rm(dst)
        islink(dst) && rm(dst)
        cp(src, dst; force=true)
    end

    # Run NJOY
    output_path = joinpath(workdir, "output")
    error_path = joinpath(workdir, "error")

    env = copy(ENV)
    env["LD_LIBRARY_PATH"] = get(ENV, "LD_LIBRARY_PATH", "") *
        ":" * dirname(njoy_binary)

    cmd = pipeline(Cmd(`$njoy_binary`; dir=workdir, env=env);
                   stdin=input_path, stdout=output_path, stderr=error_path)

    try
        run(cmd)
    catch ex
        @warn "Fortran NJOY failed" workdir exception=ex
        if isfile(error_path)
            err_text = read(error_path, String)
            !isempty(err_text) && @warn "NJOY stderr" err_text[1:min(500, length(err_text))]
        end
        return ""
    end

    # Return the output tape path
    tape90 = joinpath(workdir, "tape90")
    isfile(tape90) ? tape90 : ""
end

# =========================================================================
# Per-test oracle generation
# =========================================================================

"""
    generate_module_references(test_number; njoy_binary, cache_dir)

Generate per-module reference PENDFs for one test by running Fortran NJOY
with truncated input decks. Caches results.

Returns Dict{String, String} mapping step label (e.g. "reconr", "thermr_2") to ASCII PENDF path.
"""
function generate_module_references(test_number::Int;
                                     njoy_binary::String=NJOY_BINARY,
                                     cache_dir::String=ORACLE_CACHE)
    test_dir = joinpath(NJOY_REF_ORACLE, "tests", lpad(string(test_number), 2, '0'))
    isdir(test_dir) || error("Test directory not found: $test_dir")

    input_file = joinpath(test_dir, "input")
    isfile(input_file) || error("No input file: $input_file")
    cmake_file = joinpath(test_dir, "CMakeLists.txt")

    # Parse deck
    deck = parse_njoy_input(input_file)

    # Parse tape mappings from CMakeLists.txt
    tape_files = Dict{Int,String}()
    if isfile(cmake_file)
        ct = read(cmake_file, String)
        for (unit, fname) in parse_cmake_tapes(ct)
            fp = joinpath(RESOURCES_ORACLE, fname)
            !isfile(fp) && (fp = joinpath(test_dir, fname))
            isfile(fp) && (tape_files[unit] = fp)
        end
    end

    # Cache directory for this test
    test_cache = joinpath(cache_dir, @sprintf("test%02d", test_number))
    mkpath(test_cache)

    results = Dict{String, String}()

    # Find truncation points: after each PENDF-producing module
    for (idx, mc) in enumerate(deck.calls)
        mc.name in PENDF_MODULES || continue

        label = _step_label(deck.calls, idx)
        cache_file = joinpath(test_cache, "after_$(label).pendf")

        if isfile(cache_file) && filesize(cache_file) > 100
            results[label] = cache_file
            @info "  Cached: $label" path=cache_file
            continue
        end

        # Build truncated deck and run
        deck_text, _ = build_truncated_deck(input_file, deck, idx)
        workdir = joinpath(test_cache, "run_$(label)")

        @info "  Running Fortran NJOY (after $label)..."
        tape90 = run_fortran_njoy(deck_text, tape_files;
                                  workdir=workdir, njoy_binary=njoy_binary)

        if !isempty(tape90) && isfile(tape90) && filesize(tape90) > 100
            cp(tape90, cache_file; force=true)
            results[label] = cache_file
            @info "  Generated: $label" size=filesize(cache_file)
        else
            @warn "  Failed to generate reference for $label"
        end
    end

    results
end

"""Build a unique step label for the module at index `idx`."""
function _step_label(calls::Vector{ModuleCall}, idx::Int)
    name = string(calls[idx].name)
    # Count occurrences of this module name before idx
    n = count(i -> calls[i].name == calls[idx].name, 1:idx)
    n > 1 ? "$(name)_$(n)" : name
end

# =========================================================================
# Binary availability check
# =========================================================================

"""Check if the NJOY2016 Fortran binary is available."""
function njoy_binary_available(;njoy_binary::String=NJOY_BINARY)
    isfile(njoy_binary) && return true
    # Try building
    build_dir = joinpath(NJOY_REF_ORACLE, "build")
    isfile(joinpath(build_dir, "njoy")) && return true
    false
end

"""Return the path to the NJOY binary, building if needed."""
function ensure_njoy_binary(;njoy_binary::String=NJOY_BINARY)
    isfile(njoy_binary) && return njoy_binary

    # Check alternate location
    alt = joinpath(NJOY_REF_ORACLE, "build", "njoy")
    isfile(alt) && return alt

    error("NJOY2016 binary not found. Build with:\n" *
          "  cd njoy-reference && mkdir -p build && cd build\n" *
          "  cmake .. -DCMAKE_BUILD_TYPE=Release && make -j\$(nproc)")
end
