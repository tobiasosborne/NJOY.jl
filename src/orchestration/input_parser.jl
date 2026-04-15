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
    :dtfr, :matxsr, :viewr,
])

const SKIP_MODULES = Set([:plotr])

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
    # Card 7 (optional): weight function parameters
    # iwt=1: tabulated weight function → skip 1 card (plus tab1 data, simplified here)
    # iwt=4: ehi, sigpot, ebreak, tb[, tc, eb, ef] → skip 1 card
    # iwt=5: EPRI-cell weight function → skip 1 card
    if iwt in (1, 4, 5) && ci <= length(cards)
        ci += 1
    end
    # MT list: read (mfd, mtd, name) until mfd=0 — first temperature block only
    # (Fortran groupr uses the first block's MT list for all temperatures)
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

# =========================================================================
# GAMINR parser — matches gaminr.f90 ruing subroutine
# Card 1: nendf npend ngam1 ngam2
# Card 2: matb igg iwt lord iprint
# Card 3: title
# Card 4: (igg=1 only) ngg egg(1)..egg(ngg+1)
# Card 5: (iwt=1 only) weight function data
# Then repeat: card 6 (mfd mtd name) until mfd=0/-1, card 7 (matd, 0=stop)
# =========================================================================

struct GaminrParams
    nendf::Int; npend::Int; ngam1::Int; ngam2::Int
    matb::Int; igg::Int; iwt::Int; lord::Int; iprint::Int
    title::String
    materials::Vector{Tuple{Int,Vector{Tuple{Int,Int,String}}}}
    # (mat, [(mfd,mtd,name)...])
end

function parse_gaminr(mc::ModuleCall)::GaminrParams
    cards = mc.raw_cards
    isempty(cards) && return GaminrParams(0,0,0,0,0,3,3,0,1,"",[])
    # Card 1: unit numbers
    nendf = abs(_fint(cards[1], 1)); npend = abs(_fint(cards[1], 2))
    ngam1 = _fint(cards[1], 3; default=0); ngam2 = abs(_fint(cards[1], 4; default=0))
    # Card 2: control parameters
    matb   = length(cards) >= 2 ? _fint(cards[2], 1) : 0
    igg    = length(cards) >= 2 ? _fint(cards[2], 2; default=3) : 3
    iwt    = length(cards) >= 2 ? _fint(cards[2], 3; default=3) : 3
    lord   = length(cards) >= 2 ? _fint(cards[2], 4; default=0) : 0
    iprint = length(cards) >= 2 ? _fint(cards[2], 5; default=1) : 1
    # Card 3: title
    title = ""
    ci = 3
    if ci <= length(cards)
        title = strip(replace(join(cards[ci], " "), r"^['\"]|['\"]$" => ""))
        ci += 1
    end
    # Card 4: custom group structure (igg=1 only) — skip
    if igg == 1 && ci <= length(cards)
        ci += 1  # ngg line
        ci <= length(cards) && (ci += 1)  # egg values
    end
    # Card 5: custom weight function (iwt=1 only) — skip
    if iwt == 1 && ci <= length(cards)
        ci += 1
    end
    # Material loop: card 6 (mfd mtd name) + card 7 (matd)
    # First material uses matb from card 2
    materials = Tuple{Int,Vector{Tuple{Int,Int,String}}}[]
    cur_mat = matb
    while cur_mat != 0
        reactions = Tuple{Int,Int,String}[]
        while ci <= length(cards)
            mfd = _fint(cards[ci], 1; default=0)
            if mfd == 0
                ci += 1; break
            elseif mfd == -1
                push!(reactions, (-1, 0, ""))
                ci += 1; break
            else
                mtd = _fint(cards[ci], 2; default=0)
                name = length(cards[ci]) >= 3 ?
                    strip(replace(join(cards[ci][3:end], " "), r"^['\"]|['\"]$" => "")) : ""
                push!(reactions, (mfd, mtd, name))
                ci += 1
            end
        end
        push!(materials, (cur_mat, reactions))
        # Card 7: next matd (0=stop)
        if ci <= length(cards)
            cur_mat = _fint(cards[ci], 1; default=0)
            ci += 1
        else
            cur_mat = 0
        end
    end
    GaminrParams(nendf, npend, ngam1, ngam2, matb, igg, iwt, lord, iprint,
                 title, materials)
end

# =========================================================================
# DTFR parser — matches dtfr.f90 ruin subroutine
# Card 1: nin nout npend nplot
# Card 2: iprint ifilm iedit
# iedit=0: Card 3 (nlmax ng iptotl ipingp itabl ned ntherm)
#          Card 3a (ntherm≠0: mti mtc nlc)
#          Card 4 (edit names, iptotl-3 entries)
#          Card 5 (edit specs: ned triplets jpos mt mult)
# iedit=1: Card 6 (nlmax ng)
# Card 7: nptabl ngp
# Card 8: (repeat) hisnam mat jsigz dtemp — blank terminates
# =========================================================================

struct DtfrEditSpec
    jpos::Int; mt::Int; mult::Int
end

struct DtfrMaterial
    name::String; mat::Int; jsigz::Int; dtemp::Float64
end

struct DtfrParams
    nin::Int; nout::Int; npend::Int; nplot::Int
    iprint::Int; ifilm::Int; iedit::Int
    nlmax::Int; ng::Int; iptotl::Int; ipingp::Int; itabl::Int
    ned::Int; ntherm::Int; mti::Int; mtc::Int; nlc::Int
    edit_names::Vector{String}; edit_specs::Vector{DtfrEditSpec}
    nptabl::Int; ngp::Int
    materials::Vector{DtfrMaterial}
end

function parse_dtfr(mc::ModuleCall)::DtfrParams
    cards = mc.raw_cards
    isempty(cards) && return DtfrParams(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,[],[],0,0,[])
    # Card 1: units
    nin   = abs(_fint(cards[1], 1)); nout  = abs(_fint(cards[1], 2; default=0))
    npend = abs(_fint(cards[1], 3; default=0)); nplot = abs(_fint(cards[1], 4; default=0))
    # Card 2: options
    iprint = length(cards) >= 2 ? _fint(cards[2], 1; default=0) : 0
    ifilm  = length(cards) >= 2 ? _fint(cards[2], 2; default=0) : 0
    iedit  = length(cards) >= 2 ? _fint(cards[2], 3; default=0) : 0
    ci = 3
    nlmax = 0; ng = 0; iptotl = 0; ipingp = 0; itabl = 0; ned = 0; ntherm = 0
    mti = 0; mtc = 0; nlc = 0
    edit_names = String[]; edit_specs = DtfrEditSpec[]
    if iedit == 0
        # Card 3: table parameters
        if ci <= length(cards)
            nlmax  = _fint(cards[ci], 1; default=0)
            ng     = _fint(cards[ci], 2; default=0)
            iptotl = _fint(cards[ci], 3; default=0)
            ipingp = _fint(cards[ci], 4; default=0)
            itabl  = _fint(cards[ci], 5; default=0)
            ned    = _fint(cards[ci], 6; default=0)
            ntherm = _fint(cards[ci], 7; default=0)
            ci += 1
        end
        # Card 3a: thermal MTs (ntherm≠0)
        if ntherm != 0 && ci <= length(cards)
            mti = _fint(cards[ci], 1; default=0)
            mtc = _fint(cards[ci], 2; default=0)
            nlc = _fint(cards[ci], 3; default=0)
            ci += 1
        end
        # Card 4: edit names (nedpos = iptotl - 3 entries)
        nedpos = max(iptotl - 3, 0)
        if nedpos > 0 && ci <= length(cards)
            for j in 1:min(nedpos, length(cards[ci]))
                push!(edit_names, strip(replace(cards[ci][j], r"^['\"]|['\"]$" => "")))
            end
            ci += 1
        end
        # Card 5: edit specs (ned triplets)
        if ned > 0 && ci <= length(cards)
            j = 1
            for k in 1:ned
                jpos = _fint(cards[ci], j; default=0)
                mt_  = _fint(cards[ci], j+1; default=0)
                mult = _fint(cards[ci], j+2; default=1)
                push!(edit_specs, DtfrEditSpec(jpos, mt_, mult))
                j += 3
            end
            ci += 1
        end
    else
        # iedit=1: Card 6
        if ci <= length(cards)
            nlmax = _fint(cards[ci], 1; default=5)
            ng    = _fint(cards[ci], 2; default=30)
            ci += 1
        end
    end
    # Card 7: photon table parameters
    nptabl = 0; ngp = 0
    if ci <= length(cards)
        nptabl = _fint(cards[ci], 1; default=0)
        ngp    = _fint(cards[ci], 2; default=0)
        ci += 1
    end
    # Card 8: materials (repeat until blank/empty)
    materials = DtfrMaterial[]
    while ci <= length(cards)
        isempty(cards[ci]) && break
        name_tok = length(cards[ci]) >= 1 ? strip(replace(cards[ci][1], r"^['\"]|['\"]$" => "")) : ""
        matd = _fint(cards[ci], 2; default=0)
        matd == 0 && isempty(name_tok) && break
        jsigz = _fint(cards[ci], 3; default=1)
        dtemp = _fnum(cards[ci], 4; default=300.0)
        push!(materials, DtfrMaterial(name_tok, matd, jsigz, dtemp))
        ci += 1
    end
    DtfrParams(nin, nout, npend, nplot, iprint, ifilm, iedit,
               nlmax, ng, iptotl, ipingp, itabl, ned, ntherm, mti, mtc, nlc,
               edit_names, edit_specs, nptabl, ngp, materials)
end

# =========================================================================
# MATXSR parser — matches matxsr.f90 ruinm subroutine
# Card 1: ngen1 ngen2 nmatx (from main matxsr, not ruinm)
# Card 2: ivers huse  (from main matxsr)
# Card 3: npart ntype nholl nmat  (from ruinm)
# Card 4: (nholl cards) hollerith set ID
# Card 5: particle names (npart)
# Card 6: group counts (npart)
# Card 7: data type names (ntype)
# Card 8: jinp (ntype)
# Card 9: joutp (ntype)
# Card 10: (nmat cards) hmat matno matgg
# =========================================================================

struct MatxsrMaterial
    name::String; matno::Int; matgg::Int
end

struct MatxsrParams
    ngen1::Int; ngen2::Int; nmatx::Int
    ivers::Int; huse::String
    npart::Int; ntype::Int; nholl::Int; nmat::Int
    hsetid::Vector{String}
    particles::Vector{String}; ngrp::Vector{Int}
    dtypes::Vector{String}; jinp::Vector{Int}; joutp::Vector{Int}
    materials::Vector{MatxsrMaterial}
end

function parse_matxsr(mc::ModuleCall)::MatxsrParams
    cards = mc.raw_cards
    isempty(cards) && return MatxsrParams(0,0,0,0,"",0,0,0,0,[],[],[],[],[],[],[])
    # Card 1: unit numbers
    ngen1 = _fint(cards[1], 1; default=0)
    ngen2 = _fint(cards[1], 2; default=0)
    nmatx = abs(_fint(cards[1], 3; default=0))
    # Card 2: ivers, huse
    ivers = length(cards) >= 2 ? _fint(cards[2], 1; default=0) : 0
    huse  = length(cards) >= 2 && length(cards[2]) >= 2 ?
        strip(replace(join(cards[2][2:end], " "), r"^['\"]|['\"]$" => "")) : ""
    ci = 3
    # Card 3: npart, ntype, nholl, nmat
    npart = ci <= length(cards) ? _fint(cards[ci], 1; default=1) : 1
    ntype = ci <= length(cards) ? _fint(cards[ci], 2; default=1) : 1
    nholl = ci <= length(cards) ? _fint(cards[ci], 3; default=1) : 1
    nmat  = ci <= length(cards) ? _fint(cards[ci], 4; default=1) : 1
    ci += 1
    # Card 4: hollerith set ID (nholl cards)
    hsetid = String[]
    for k in 1:nholl
        ci > length(cards) && break
        push!(hsetid, strip(replace(join(cards[ci], " "), r"^['\"]|['\"]$" => "")))
        ci += 1
    end
    # Card 5: particle names
    particles = String[]
    if ci <= length(cards)
        for j in 1:min(npart, length(cards[ci]))
            push!(particles, strip(replace(cards[ci][j], r"^['\"]|['\"]$" => "")))
        end
        ci += 1
    end
    # Card 6: group counts
    ngrp = Int[]
    if ci <= length(cards)
        for j in 1:min(npart, length(cards[ci]))
            push!(ngrp, _fint(cards[ci], j))
        end
        ci += 1
    end
    # Card 7: data type names
    dtypes = String[]
    if ci <= length(cards)
        for j in 1:min(ntype, length(cards[ci]))
            push!(dtypes, strip(replace(cards[ci][j], r"^['\"]|['\"]$" => "")))
        end
        ci += 1
    end
    # Card 8: jinp
    jinp = Int[]
    if ci <= length(cards)
        for j in 1:min(ntype, length(cards[ci]))
            push!(jinp, _fint(cards[ci], j))
        end
        ci += 1
    end
    # Card 9: joutp
    joutp = Int[]
    if ci <= length(cards)
        for j in 1:min(ntype, length(cards[ci]))
            push!(joutp, _fint(cards[ci], j))
        end
        ci += 1
    end
    # Card 10: materials (nmat cards)
    materials = MatxsrMaterial[]
    for k in 1:nmat
        ci > length(cards) && break
        name = length(cards[ci]) >= 1 ?
            strip(replace(cards[ci][1], r"^['\"]|['\"]$" => "")) : ""
        matno = _fint(cards[ci], 2; default=0)
        matgg = _fint(cards[ci], 3; default=0)
        push!(materials, MatxsrMaterial(name, matno, matgg))
        ci += 1
    end
    MatxsrParams(ngen1, ngen2, nmatx, ivers, huse,
                 npart, ntype, nholl, nmat, hsetid,
                 particles, ngrp, dtypes, jinp, joutp, materials)
end

# =========================================================================
# VIEWR parser — matches viewr.f90
# Card 1: infile nps  (from nsysi — the NJOY input)
# All other data comes from infile (the plot tape)
# =========================================================================

struct ViewrParams
    infile::Int; nps::Int
end

function parse_viewr(mc::ModuleCall)::ViewrParams
    cards = mc.raw_cards
    isempty(cards) && return ViewrParams(0, 0)
    infile = abs(_fint(cards[1], 1; default=0))
    nps    = abs(_fint(cards[1], 2; default=0))
    ViewrParams(infile, nps)
end

struct LeaprParams
    nout::Int
end

function parse_leapr(mc::ModuleCall)::LeaprParams
    cards = mc.raw_cards
    isempty(cards) && return LeaprParams(0)
    LeaprParams(abs(_fint(cards[1], 1)))
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

struct CovrParams
    nin::Int      # input covariance tape from errorr (or zero)
    npend::Int    # input pendf tape (or zero; unused by stub)
    nout::Int     # output plot-tape unit consumed by viewr
end

function parse_covr(mc::ModuleCall)::CovrParams
    cards = mc.raw_cards
    isempty(cards) && return CovrParams(0, 0, 0)
    nin   = abs(_fint(cards[1], 1))
    npend = abs(_fint(cards[1], 2; default=0))
    nout  = abs(_fint(cards[1], 3; default=0))
    CovrParams(nin, npend, nout)
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
