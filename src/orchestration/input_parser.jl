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
                # Strip '/' terminator, but only outside quotes
                si = _find_slash_outside_quotes(raw)
                si > 0 && (raw = raw[1:si-1])
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

"""Find the first '/' outside single/double quotes. Returns 0 if none found."""
function _find_slash_outside_quotes(s::AbstractString)::Int
    in_single = false; in_double = false
    for (i, c) in enumerate(s)
        if c == '\'' && !in_double; in_single = !in_single
        elseif c == '"' && !in_single; in_double = !in_double
        elseif c == '/' && !in_single && !in_double; return i
        end
    end
    return 0
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

struct ReconrMatSpec
    mat::Int; err::Float64; description::String
end

struct ReconrParams
    nendf::Int; npend::Int; mat::Int; err::Float64; title::String
    materials::Vector{ReconrMatSpec}
end
ReconrParams(nendf, npend, mat, err, title) =
    ReconrParams(nendf, npend, mat, err, title, [ReconrMatSpec(mat, err, "")])

function parse_reconr(mc::ModuleCall)::ReconrParams
    cards = mc.raw_cards
    isempty(cards) && return ReconrParams(0, 0, 0, 0.001, "")
    nendf = abs(_fint(cards[1], 1)); npend = abs(_fint(cards[1], 2))
    title = length(cards) >= 2 ? strip(replace(join(cards[2], " "), r"^['\"]|['\"]$" => "")) : ""
    materials = ReconrMatSpec[]
    ci = 3
    while ci <= length(cards)
        mat_val = _fint(cards[ci], 1)
        mat_val == 0 && break
        ncards = length(cards[ci]) >= 2 ? _fint(cards[ci], 2; default=0) : 0
        ngrid = length(cards[ci]) >= 3 ? _fint(cards[ci], 3; default=0) : 0
        err_val = ci + 1 <= length(cards) ? _fnum(cards[ci+1], 1; default=0.001) : 0.001
        desc = ""
        if ncards > 0 && ci + 2 <= length(cards)
            desc = strip(replace(join(cards[ci+2], " "), r"^['\"]|['\"]$" => ""))
        end
        push!(materials, ReconrMatSpec(abs(mat_val), err_val, desc))
        ci += 2 + max(ncards, 0) + max(ngrid, 0)
    end
    first_mat = isempty(materials) ? 0 : materials[1].mat
    first_err = isempty(materials) ? 0.001 : materials[1].err
    ReconrParams(nendf, npend, first_mat, first_err, title, materials)
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

struct ThermrParams
    nin_thermal::Int; nin_pendf::Int; nout::Int
    mat_thermal::Int; mat::Int; nbin::Int; ntemp::Int; iinc::Int; icoh::Int
    iform::Int; natom::Int; mtref::Int; iprint::Int
    temperatures::Vector{Float64}; tol::Float64; emax::Float64
end

function parse_thermr(mc::ModuleCall)::ThermrParams
    cards = mc.raw_cards
    isempty(cards) && return ThermrParams(0,0,0,0,0,8,1,0,0,0,1,221,0,Float64[],0.05,10.0)
    nin_th = _fint(cards[1], 1); nin_p = abs(_fint(cards[1], 2))
    nout = length(cards[1]) >= 3 ? abs(_fint(cards[1], 3)) : 0
    # Card 2: matde, matdp, nbin, ntemp, iinc, icoh, iform, natom, mtref, iprint
    mat_th = length(cards) >= 2 ? _fint(cards[2], 1) : 0
    mat    = length(cards) >= 2 ? _fint(cards[2], 2) : 0
    nbin   = length(cards) >= 2 ? _fint(cards[2], 3; default=8) : 8
    ntemp  = length(cards) >= 2 ? _fint(cards[2], 4; default=1) : 1
    iinc   = length(cards) >= 2 ? _fint(cards[2], 5; default=0) : 0
    icoh   = length(cards) >= 2 ? _fint(cards[2], 6; default=0) : 0
    iform  = length(cards) >= 2 ? _fint(cards[2], 7; default=0) : 0
    natom  = length(cards) >= 2 ? _fint(cards[2], 8; default=1) : 1
    mtref  = length(cards) >= 2 ? _fint(cards[2], 9; default=221) : 221
    iprint = length(cards) >= 2 ? _fint(cards[2], 10; default=0) : 0
    temps = Float64[]
    if length(cards) >= 3
        for t in cards[3]; startswith(t, "'") && continue; push!(temps, _parse_num(t)); end
    end
    tol = length(cards) >= 4 ? _fnum(cards[4], 1; default=0.05) : 0.05
    emax = length(cards) >= 4 ? _fnum(cards[4], 2; default=10.0) : 10.0
    ThermrParams(nin_th, nin_p, nout, mat_th, mat, nbin, max(ntemp, length(temps)),
                 iinc, icoh, iform, natom, mtref, iprint, temps, tol, emax)
end

struct AcerParams
    nendf::Int; npendf::Int; nace::Int; ndir::Int
    mat::Int; iopt::Int; temp::Float64; suffix::String
end

function parse_acer(mc::ModuleCall)::AcerParams
    cards = mc.raw_cards
    isempty(cards) && return AcerParams(0,0,0,0,0,1,300.0,"80c")
    nendf = abs(_fint(cards[1], 1)); npendf = abs(_fint(cards[1], 2))
    nace = length(cards[1]) >= 3 ? abs(_fint(cards[1], 3)) : 0
    ndir = length(cards[1]) >= 4 ? abs(_fint(cards[1], 4)) : 0
    iopt = length(cards) >= 2 ? _fint(cards[2], 1; default=1) : 1
    # Card 3 is a title string (hz); Card 4 has mat, temp for iopt=1
    mat = 0; temp = 300.0; suffix = "80c"
    if iopt == 1 && length(cards) >= 4
        # Card 4: matd, tempd, [local], [iprint]
        card4 = cards[4]
        mat = _fint(card4, 1)
        temp = _fnum(card4, 2; default=300.0)
    end
    AcerParams(nendf, npendf, nace, ndir, mat, iopt, temp, suffix)
end

struct GasprParams
    nendf::Int; npendf_in::Int; npendf_out::Int
end

function parse_gaspr(mc::ModuleCall)::GasprParams
    cards = mc.raw_cards
    isempty(cards) && return GasprParams(0,0,0)
    GasprParams(abs(_fint(cards[1], 1)), abs(_fint(cards[1], 2)),
                length(cards[1]) >= 3 ? abs(_fint(cards[1], 3)) : 0)
end

struct UnresrParams
    nendf::Int; npendf_in::Int; npendf_out::Int
    mat::Int; ntemp::Int; nsigz::Int
    temperatures::Vector{Float64}; sigz::Vector{Float64}
end

function parse_unresr(mc::ModuleCall)::UnresrParams
    cards = mc.raw_cards
    isempty(cards) && return UnresrParams(0,0,0,0,1,1,Float64[],Float64[])
    nendf = abs(_fint(cards[1], 1)); nin = abs(_fint(cards[1], 2))
    nout = length(cards[1]) >= 3 ? abs(_fint(cards[1], 3)) : 0
    mat = length(cards) >= 2 ? _fint(cards[2], 1) : 0
    ntemp = length(cards) >= 2 ? _fint(cards[2], 2; default=1) : 1
    nsigz = length(cards) >= 2 ? _fint(cards[2], 3; default=1) : 1
    temps = Float64[]
    if length(cards) >= 3
        for t in cards[3]; startswith(t, "'") && continue
            v = _parse_num(t); v == 0 && break; push!(temps, v); end
    end
    sigz = Float64[]
    if length(cards) >= 4
        for t in cards[4]; startswith(t, "'") && continue
            v = _parse_num(t); v == 0 && break; push!(sigz, v); end
    end
    UnresrParams(nendf, nin, nout, mat, max(ntemp, length(temps)),
                 max(nsigz, length(sigz)), temps, sigz)
end

struct PurrParams
    nendf::Int; npendf_in::Int; npendf_out::Int
    mat::Int; ntemp::Int; nsigz::Int; nbin::Int; nladr::Int
    temperatures::Vector{Float64}; sigz::Vector{Float64}
end

struct ErrorrParams
    nendf::Int; npend::Int; ngout::Int; nout::Int; nin::Int; nstan::Int
    mat::Int; ign::Int; iwt::Int; iprint::Int; irelco::Int
    mprint::Int; tempin::Float64
    iread::Int; mfcov::Int; irespr::Int; legord::Int
    user_egn::Vector{Float64}  # user energy group boundaries (when ign<0)
end

function parse_errorr(mc::ModuleCall)::ErrorrParams
    cards = mc.raw_cards
    isempty(cards) && return ErrorrParams(0,0,0,0,0,0,0,1,6,1,1,0,0.0,0,33,1,1,Float64[])
    # Card 1: nendf, npend, ngout, nout, nin, nstan
    nendf = abs(_fint(cards[1], 1)); npend = abs(_fint(cards[1], 2))
    ngout = _fint(cards[1], 3; default=0); nout = abs(_fint(cards[1], 4; default=0))
    nin = abs(_fint(cards[1], 5; default=0)); nstan = abs(_fint(cards[1], 6; default=0))
    # Card 2: matd, ign, iwt, iprint, irelco
    mat = length(cards) >= 2 ? _fint(cards[2], 1) : 0
    ign = length(cards) >= 2 ? _fint(cards[2], 2; default=1) : 1
    iwt = length(cards) >= 2 ? _fint(cards[2], 3; default=6) : 6
    iprint = length(cards) >= 2 ? _fint(cards[2], 4; default=1) : 1
    irelco = length(cards) >= 2 ? _fint(cards[2], 5; default=1) : 1
    # Card 3: mprint, tempin
    mprint = length(cards) >= 3 ? _fint(cards[3], 1; default=0) : 0
    tempin = length(cards) >= 3 ? _fnum(cards[3], 2; default=0.0) : 0.0
    # Card 4: iread, mfcov, irespr, legord
    iread = length(cards) >= 4 ? _fint(cards[4], 1; default=0) : 0
    mfcov = length(cards) >= 4 ? _fint(cards[4], 2; default=33) : 33
    irespr = length(cards) >= 4 ? _fint(cards[4], 3; default=1) : 1
    legord = length(cards) >= 4 ? _fint(cards[4], 4; default=1) : 1
    # User group structure (ign < 0 or ign == 1): read ngn then boundaries
    user_egn = Float64[]
    if ign < 0 || ign == 1
        ci = 5
        ngn = ci <= length(cards) ? _fint(cards[ci], 1; default=0) : 0
        ci += 1
        while ci <= length(cards) && length(user_egn) < ngn + 1
            for t in cards[ci]
                startswith(t, "'") && continue
                push!(user_egn, _parse_num(t))
                length(user_egn) >= ngn + 1 && break
            end
            ci += 1
        end
    end
    ErrorrParams(nendf, npend, ngout, nout, nin, nstan,
                 mat, ign, iwt, iprint, irelco, mprint, tempin,
                 iread, mfcov, irespr, legord, user_egn)
end

struct GrouprParams
    nendf::Int; npend::Int; ngout_in::Int; nout::Int
    mat::Int; ign::Int; igg::Int; iwt::Int; lord::Int
    ntemp::Int; nsigz::Int; iprint::Int
    title::String
    temperatures::Vector{Float64}; sigz::Vector{Float64}
    mt_list::Vector{Tuple{Int,Int,String}}  # (mfd, mtd, name) triplets
end

function parse_groupr(mc::ModuleCall)::GrouprParams
    cards = mc.raw_cards
    isempty(cards) && return GrouprParams(0,0,0,0,0,3,0,6,0,1,1,0,"",Float64[],Float64[],Tuple{Int,Int,String}[])
    # Card 1: nendf, npend, ngout, nout
    nendf = abs(_fint(cards[1], 1)); npend = abs(_fint(cards[1], 2))
    ngout_in = _fint(cards[1], 3; default=0); nout = abs(_fint(cards[1], 4; default=0))
    # Card 2: matd, ign, igg, iwt, lord, ntemp, nsigz, iprint
    mat   = length(cards) >= 2 ? _fint(cards[2], 1) : 0
    ign   = length(cards) >= 2 ? _fint(cards[2], 2; default=3) : 3
    igg   = length(cards) >= 2 ? _fint(cards[2], 3; default=0) : 0
    iwt   = length(cards) >= 2 ? _fint(cards[2], 4; default=6) : 6
    lord  = length(cards) >= 2 ? _fint(cards[2], 5; default=0) : 0
    ntemp = length(cards) >= 2 ? _fint(cards[2], 6; default=1) : 1
    nsigz = length(cards) >= 2 ? _fint(cards[2], 7; default=1) : 1
    iprint= length(cards) >= 2 ? _fint(cards[2], 8; default=0) : 0
    # Card 3: title string
    title = ""
    ci = 3
    if ci <= length(cards)
        title = strip(replace(join(cards[ci], " "), r"^['\"]|['\"]$" => ""))
        ci += 1
    end
    # Card 4+: temperatures (ntemp values)
    temps = Float64[]
    while ci <= length(cards) && length(temps) < ntemp
        for t in cards[ci]
            startswith(t, "'") && continue
            push!(temps, _parse_num(t))
            length(temps) >= ntemp && break
        end
        ci += 1
    end
    # Card 5+: sigz values (nsigz values)
    sigz = Float64[]
    while ci <= length(cards) && length(sigz) < nsigz
        for t in cards[ci]
            startswith(t, "'") && continue
            push!(sigz, _parse_num(t))
            length(sigz) >= nsigz && break
        end
        ci += 1
    end
    # MT list: read (mfd, mtd, name) until mfd=0
    mt_list = Tuple{Int,Int,String}[]
    while ci <= length(cards)
        mfd = _fint(cards[ci], 1; default=0)
        mfd == 0 && break
        mtd = _fint(cards[ci], 2; default=0)
        name = length(cards[ci]) >= 3 ?
            strip(replace(join(cards[ci][3:end], " "), r"^['\"]|['\"]$" => "")) : ""
        push!(mt_list, (mfd, mtd, name))
        ci += 1
    end
    GrouprParams(nendf, npend, ngout_in, nout, mat, ign, igg, iwt, lord,
                 ntemp, nsigz, iprint, title, temps, sigz, mt_list)
end

function parse_purr(mc::ModuleCall)::PurrParams
    cards = mc.raw_cards
    isempty(cards) && return PurrParams(0,0,0,0,1,1,20,64,Float64[],Float64[])
    nendf = abs(_fint(cards[1], 1)); nin = abs(_fint(cards[1], 2))
    nout = length(cards[1]) >= 3 ? abs(_fint(cards[1], 3)) : 0
    mat = length(cards) >= 2 ? _fint(cards[2], 1) : 0
    ntemp = length(cards) >= 2 ? _fint(cards[2], 2; default=1) : 1
    nsigz = length(cards) >= 2 ? _fint(cards[2], 3; default=1) : 1
    nbin = length(cards) >= 2 ? _fint(cards[2], 4; default=20) : 20
    nladr = length(cards) >= 2 ? _fint(cards[2], 5; default=64) : 64
    temps = Float64[]
    if length(cards) >= 3
        for t in cards[3]; startswith(t, "'") && continue
            v = _parse_num(t); v == 0 && break; push!(temps, v); end
    end
    sigz = Float64[]
    if length(cards) >= 4
        for t in cards[4]; startswith(t, "'") && continue
            v = _parse_num(t); v == 0 && break; push!(sigz, v); end
    end
    PurrParams(nendf, nin, nout, mat, max(ntemp, length(temps)),
               max(nsigz, length(sigz)), nbin, nladr, temps, sigz)
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
