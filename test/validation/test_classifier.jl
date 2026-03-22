# test_classifier.jl — Classify all 85 NJOY tests for batch diagnostic prioritization
#
# Groups tests by module chain pattern, resonance formalism, material weight
# class, and expected difficulty.

using NJOY
using Printf

include("input_parser.jl")

const NJOY_REF_CL = normpath(joinpath(@__DIR__, "..", "..", "njoy-reference"))
const RESOURCES_CL = joinpath(NJOY_REF_CL, "tests", "resources")
const TESTS_DIR_CL = joinpath(NJOY_REF_CL, "tests")

# =========================================================================
# Types
# =========================================================================

struct TestClassification
    test_number::Int
    chain_modules::Vector{Symbol}   # ordered, excluding moder
    chain_pattern::Symbol           # canonical pattern
    formalism::Symbol               # :slbw, :mlbw, :rm, :sammy, :none, :mixed
    material_class::Symbol          # :light, :medium, :heavy, :actinide
    mat::Int
    za::Float64
    has_thermal::Bool
    has_urr::Bool
    has_fission::Bool
    n_ref_mf3_mts::Int              # MTs in best reference tape
    difficulty::Symbol              # :simple, :moderate, :complex
end

# =========================================================================
# Classification
# =========================================================================

function classify_test(test_number::Int)
    dir = joinpath(TESTS_DIR_CL, lpad(string(test_number), 2, '0'))
    input_file = joinpath(dir, "input")
    cmake_file = joinpath(dir, "CMakeLists.txt")

    # Defaults for missing/empty tests
    default = TestClassification(test_number, Symbol[], :empty, :none, :light,
                                  0, 0.0, false, false, false, 0, :simple)
    isfile(input_file) || return default

    deck = parse_njoy_input(input_file)
    mods = [mc.name for mc in deck.calls if mc.name != :moder && mc.name != :stop]
    isempty(mods) && return default

    # Chain pattern: unique ordered module names
    unique_mods = unique(filter(m -> m != :moder && !(m in SKIP_MODULES), mods))
    pattern = isempty(unique_mods) ? :empty : Symbol(join(string.(unique_mods), "_"))

    # Extract MAT from reconr
    mat = 0; err = 0.001
    for mc in deck.calls
        if mc.name == :reconr && length(mc.raw_cards) >= 3
            p = parse_reconr(mc); mat = p.mat; err = p.err; break
        end
    end

    # Find ENDF file and extract ZA
    za = 0.0
    endf_file = ""
    if isfile(cmake_file)
        ct = read(cmake_file, String)
        for (unit, fname) in parse_cmake_tapes(ct)
            fp = joinpath(RESOURCES_CL, fname)
            !isfile(fp) && (fp = joinpath(dir, fname))
            if isfile(fp) && unit in [20, 21, 30]
                endf_file = fp
                break
            end
        end
    end

    if !isempty(endf_file) && mat > 0
        try
            open(endf_file) do io
                while !eof(io)
                    line = readline(io)
                    p = rpad(line, 80)
                    mf = NJOY._parse_int(p[71:72])
                    mt = NJOY._parse_int(p[73:75])
                    if mf == 1 && mt == 451
                        za = parse_endf_float(p[1:11])
                        break
                    end
                end
            end
        catch; end
    end

    # Material class from Z
    z = floor(Int, za / 1000)
    mat_class = z <= 0 ? :unknown : z <= 20 ? :light : z <= 56 ? :medium :
                z <= 88 ? :heavy : :actinide

    # Formalism from MF2
    formalism = :none
    if !isempty(endf_file) && mat > 0
        formalism = detect_formalism(endf_file, mat)
    end

    # Flags
    has_thermal = :thermr in mods || :leapr in mods
    has_urr = :unresr in mods || :purr in mods
    has_fission = za > 0 && z >= 89

    # Count MF3 MTs in reference tapes
    n_ref_mts = 0
    if isdir(dir)
        for f in readdir(dir)
            startswith(f, "referenceTape") || continue
            ref_data = _count_mf3_mts(joinpath(dir, f))
            ref_data > n_ref_mts && (n_ref_mts = ref_data)
        end
    end

    # Difficulty
    n_real_mods = length(filter(m -> m in [:reconr, :broadr, :heatr, :thermr,
                                            :unresr, :purr, :groupr, :errorr,
                                            :acer, :gaspr, :leapr, :gaminr], mods))
    difficulty = n_real_mods <= 2 ? :simple :
                 n_real_mods <= 5 ? :moderate : :complex

    TestClassification(test_number, unique_mods, pattern, formalism, mat_class,
                       mat, za, has_thermal, has_urr, has_fission, n_ref_mts, difficulty)
end

"""Classify all tests. Returns sorted vector."""
function classify_all(; tests::AbstractVector{Int}=1:85)
    [classify_test(t) for t in tests if isdir(joinpath(TESTS_DIR_CL, lpad(string(t), 2, '0')))]
end

"""Detect the primary resonance formalism from MF2/MT151."""
function detect_formalism(endf_file::String, mat::Int)
    try
        open(endf_file) do io
            in_mf2 = false
            while !eof(io)
                line = readline(io)
                length(line) < 75 && continue
                p = rpad(line, 80)
                mf = NJOY._parse_int(p[71:72])
                mt = NJOY._parse_int(p[73:75])
                if mf == 2 && mt == 151
                    in_mf2 = true
                elseif in_mf2
                    # Read LRF from first range CONT record
                    # LRF is in L2 position (column 34-44 area)
                    lrf = NJOY._parse_int(p[34:44])
                    return Dict(0 => :none, 1 => :slbw, 2 => :mlbw,
                                3 => :reich_moore, 7 => :sammy)[get(
                                Dict(0=>0, 1=>1, 2=>2, 3=>3, 7=>7), lrf, 0)]
                end
                mf > 2 && break
            end
        end
    catch; end
    :none
end

"""Count distinct MF3 MT numbers in a reference tape."""
function _count_mf3_mts(filepath::String)
    mts = Set{Int}()
    try
        for line in eachline(filepath)
            length(line) < 75 && continue
            p = rpad(line, 80)
            mf = NJOY._parse_int(p[71:72])
            mt = NJOY._parse_int(p[73:75])
            mf == 3 && mt > 0 && push!(mts, mt)
        end
    catch; end
    length(mts)
end

"""Group classifications by a field."""
function group_by(classifications::Vector{TestClassification}, field::Symbol)
    groups = Dict{Any, Vector{TestClassification}}()
    for tc in classifications
        key = getfield(tc, field)
        push!(get!(groups, key, TestClassification[]), tc)
    end
    groups
end

"""Print a summary table of all classifications."""
function print_classification_table(classifications::Vector{TestClassification})
    @printf("%-5s  %-6s  %-12s  %-8s  %-10s  %5s  %5s  %s\n",
            "Test", "MAT", "Formalism", "Material", "Difficulty",
            "Therm", "URR", "Chain")
    println("─" ^ 80)
    for tc in classifications
        @printf("%-5d  %-6d  %-12s  %-8s  %-10s  %5s  %5s  %s\n",
                tc.test_number, tc.mat, tc.formalism, tc.material_class,
                tc.difficulty,
                tc.has_thermal ? "yes" : "-",
                tc.has_urr ? "yes" : "-",
                join(string.(tc.chain_modules[1:min(4, length(tc.chain_modules))]), "→"))
    end
end
