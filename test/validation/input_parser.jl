# input_parser.jl -- Parse NJOY2016 input decks into structured ModuleCall objects
#
# Handles all 23 NJOY module types. Extracts enough information to drive
# NJOY.jl API calls: module sequence, MAT numbers, tolerances, temperatures.

# =========================================================================
# Types
# =========================================================================

struct ModuleCall
    name::Symbol
    raw_cards::Vector{Vector{String}}
end

struct NJOYInputDeck
    calls::Vector{ModuleCall}
end

const NJOY_MODULES = Set([
    :moder, :reconr, :broadr, :heatr, :thermr, :unresr, :purr,
    :groupr, :errorr, :acer, :viewr, :plotr, :leapr, :gaminr,
    :dtfr, :ccccr, :matxsr, :wimsr, :gaspr, :mixr, :powr, :resxsr, :covr,
])

const AVAILABLE_MODULES = Set([
    :moder, :reconr, :broadr, :heatr, :thermr, :unresr, :purr,
    :groupr, :errorr, :acer, :gaspr, :leapr, :gaminr, :covr,
])

const SKIP_MODULES = Set([:viewr, :plotr])

# =========================================================================
# Tokeniser
# =========================================================================

function tokenise_njoy_input(text::AbstractString)::Vector{ModuleCall}
    calls = ModuleCall[]
    lines = split(text, '\n')
    i = 1
    while i <= length(lines)
        line = strip(lines[i])
        if isempty(line) || startswith(line, "--") || startswith(line, "*")
            i += 1; continue
        end
        lowercase(strip(line)) == "stop" && break
        word = lowercase(replace(line, r"\s*/.*" => ""))
        word = strip(replace(word, r"['\"]" => ""))
        word = replace(word, r"^\[x\d+\]\s*" => "")
        mod_sym = Symbol(word)
        if mod_sym in NJOY_MODULES
            i += 1
            cards = Vector{String}[]
            while i <= length(lines)
                cline = strip(lines[i])
                if startswith(cline, "--") || startswith(cline, "*")
                    i += 1; continue
                end
                occursin(r"^\[x\d+\]\s*$", cline) && (i += 1; continue)
                cline = replace(cline, r"^\[x\d+\]\s*" => "")
                trimmed = lstrip(cline)
                if !(startswith(trimmed, "'") || startswith(trimmed, "\""))
                    cword = lowercase(replace(cline, r"\s*/.*" => ""))
                    cword = strip(replace(cword, r"['\"]" => ""))
                    if Symbol(cword) in NJOY_MODULES || cword == "stop"
                        break
                    end
                end
                raw = cline
                si = findfirst('/', raw)
                si !== nothing && (raw = raw[1:si-1])
                tokens = _tokenise_card(raw)
                !isempty(tokens) && push!(cards, tokens)
                i += 1
            end
            push!(calls, ModuleCall(mod_sym, cards))
        else
            i += 1
        end
    end
    return calls
end

function _tokenise_card(raw::AbstractString)::Vector{String}
    tokens = String[]
    s = strip(raw)
    isempty(s) && return tokens
    i = 1
    while i <= length(s)
        c = s[i]
        if c == '\'' || c == '"'
            j = findnext(c, s, i + 1)
            j === nothing && (j = length(s))
            push!(tokens, s[i:j]); i = j + 1
        elseif c in (' ', '\t', ',')
            i += 1
        else
            j = i
            while j <= length(s) && !(s[j] in (' ', '\t', ',')); j += 1; end
            push!(tokens, s[i:j-1]); i = j
        end
    end
    return tokens
end

# =========================================================================
# Numeric helpers
# =========================================================================

function _parse_num(s::AbstractString)::Float64
    t = strip(s)
    isempty(t) && return 0.0
    t = replace(t, r"([0-9])([+-])(\d)" => s"\1e\2\3")
    parse(Float64, t)
end

function _parse_int_token(s::AbstractString)::Int
    t = strip(s)
    isempty(t) && return 0
    occursin('.', t) && return round(Int, _parse_num(t))
    parse(Int, t)
end

_fnum(card, n; default=0.0) = n <= length(card) ? _parse_num(card[n]) : default
_fint(card, n; default=0)   = n <= length(card) ? _parse_int_token(card[n]) : default

# =========================================================================
# Module-specific parameter extractors
# =========================================================================

struct ReconrParams
    nendf::Int; npend::Int; mat::Int; err::Float64; title::String
end

function parse_reconr(mc::ModuleCall)::ReconrParams
    cards = mc.raw_cards
    isempty(cards) && return ReconrParams(0, 0, 0, 0.001, "")
    nendf = _fint(cards[1], 1); npend = _fint(cards[1], 2)
    title = length(cards) >= 2 ? join(cards[2], " ") : ""
    mat = length(cards) >= 3 ? _fint(cards[3], 1) : 0
    err = length(cards) >= 4 ? _fnum(cards[4], 1; default=0.001) : 0.001
    ReconrParams(abs(nendf), abs(npend), mat, err, title)
end

struct BroadrParams
    nendf::Int; npendf_in::Int; npendf_out::Int
    mat::Int; ntemp::Int; tol::Float64; thnmax::Float64
    temperatures::Vector{Float64}
end

function parse_broadr(mc::ModuleCall)::BroadrParams
    cards = mc.raw_cards
    isempty(cards) && return BroadrParams(0,0,0,0,1,0.001,0.0,Float64[])
    nendf = abs(_fint(cards[1], 1)); nin = abs(_fint(cards[1], 2))
    nout = length(cards[1]) >= 3 ? abs(_fint(cards[1], 3)) : 0
    mat = length(cards) >= 2 ? _fint(cards[2], 1) : 0
    ntemp = length(cards) >= 2 ? _fint(cards[2], 2; default=1) : 1
    tol = length(cards) >= 3 ? _fnum(cards[3], 1; default=0.001) : 0.001
    thnmax = length(cards) >= 3 && length(cards[3]) >= 2 ?
             _fnum(cards[3], 2; default=0.0) : 0.0
    temps = Float64[]
    for ci in 4:min(4 + max(ntemp,1), length(cards))
        for t in cards[ci]
            startswith(t, "'") && continue; t == "0" && break
            push!(temps, _parse_num(t))
        end
        length(temps) >= max(ntemp,1) && break
    end
    BroadrParams(nendf, nin, nout, mat, max(ntemp, length(temps)), tol, thnmax, temps)
end

struct HeatrParams
    nendf::Int; npendf_in::Int; npendf_out::Int
    mat::Int; nqa::Int; mts::Vector{Int}
end

function parse_heatr(mc::ModuleCall)::HeatrParams
    cards = mc.raw_cards
    isempty(cards) && return HeatrParams(0,0,0,0,0,Int[])
    nendf = abs(_fint(cards[1], 1)); nin = abs(_fint(cards[1], 2))
    nout = length(cards[1]) >= 3 ? abs(_fint(cards[1], 3)) : 0
    mat = length(cards) >= 2 ? _fint(cards[2], 1) : 0
    nqa = length(cards) >= 2 ? _fint(cards[2], 2; default=0) : 0
    mts = Int[]
    if length(cards) >= 3 && nqa > 0
        for t in cards[3]; startswith(t, "'") && continue; push!(mts, _parse_int_token(t)); end
    end
    HeatrParams(nendf, nin, nout, mat, nqa, mts)
end

# =========================================================================
# CMakeLists.txt parser
# =========================================================================

function parse_cmake_tapes(cmake_text::AbstractString)::Dict{Int, String}
    mapping = Dict{Int, String}()
    for m in eachmatch(
        r"configure_file\(\s*\"([^\"]+)\"\s*\n?\s*\"[^\"]*/(tape(\d+))\"\s*COPYONLY\s*\)",
        cmake_text)
        fname = replace(m.captures[1], r".*/" => "")
        fname = replace(fname, r"\$\{RESOURCES\}" => "")
        fname = replace(fname, r"\$\{CMAKE_CURRENT_SOURCE_DIR\}" => "")
        mapping[parse(Int, m.captures[3])] = strip(fname, '/')
    end
    mapping
end

function parse_cmake_references(cmake_text::AbstractString)::Vector{String}
    unique([m.match for m in eachmatch(r"referenceTape\d+", cmake_text)])
end

# =========================================================================
# Top-level interface
# =========================================================================

function parse_njoy_input(input_file::AbstractString)::NJOYInputDeck
    NJOYInputDeck(tokenise_njoy_input(read(input_file, String)))
end

modules_in_deck(deck::NJOYInputDeck) = [c.name for c in deck.calls]

function runnable_modules(deck::NJOYInputDeck)
    [c.name for c in deck.calls if c.name in AVAILABLE_MODULES]
end

function missing_modules(deck::NJOYInputDeck)
    [c.name for c in deck.calls
     if !(c.name in AVAILABLE_MODULES) && !(c.name in SKIP_MODULES)]
end

function deck_category(deck::NJOYInputDeck)::Symbol
    mods = modules_in_deck(deck)
    isempty(mods) && return :empty
    data_mods = filter(m -> m != :moder, mods)
    isempty(data_mods) && return :moder_only
    fm = data_mods[1]
    fm == :leapr && return :leapr
    fm == :reconr && return :reconr_chain
    fm == :acer && return :acer_only
    fm == :errorr && return :errorr_chain
    fm == :plotr && return :plotr_only
    fm == :gaminr && return :gaminr_chain
    return :other
end
