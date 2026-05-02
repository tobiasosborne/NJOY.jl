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
    # Verbatim deck lines belonging to this module call (in order, including
    # comment-only lines and blank-card markers). Fortran-faithful parsers
    # like `parse_covr` need this because the tokeniser collapses empty
    # default-cards (`/` with no tokens) — covr treats every `/` as a
    # distinct card terminator (including `//` for two consecutive defaults).
    raw_lines::Vector{String}
end
ModuleCall(name::Symbol, raw_cards::Vector{Vector{String}}) =
    ModuleCall(name, raw_cards, String[])

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
    :dtfr, :matxsr, :viewr, :wimsr, :mixr, :resxsr, :powr,
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
            raw_kept = String[]
            while i <= length(lines)
                cline = strip(lines[i])
                # Treat `*…*` as a string card (NJOY title/comment-string convention)
                # rather than a Julia-side pure comment. Lines with only one `*`
                # (e.g. "* explanatory note") are still skipped as comments.
                if startswith(cline, "--") || (startswith(cline, "*") && !occursin(r"\*.*\*", cline))
                    i += 1; continue
                end
                occursin(r"^\[x\d+\]\s*$", cline) && (i += 1; continue)
                cline = replace(cline, r"^\[x\d+\]\s*" => "")
                trimmed = lstrip(cline)
                if !(startswith(trimmed, "'") || startswith(trimmed, "\"") || startswith(trimmed, "*"))
                    cword = lowercase(replace(cline, r"\s*/.*" => ""))
                    cword = strip(replace(cword, r"['\"]" => ""))
                    if Symbol(cword) in NJOY_MODULES || cword == "stop"
                        break
                    end
                end
                raw = cline
                push!(raw_kept, String(cline))
                # Strip '/' terminator, but only outside quotes
                si = _find_slash_outside_quotes(raw)
                si > 0 && (raw = raw[1:si-1])
                tokens = _tokenise_card(raw)
                !isempty(tokens) && push!(cards, tokens)
                i += 1
            end
            push!(calls, ModuleCall(mod_sym, cards, raw_kept))
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
        elseif c == '*'
            # `*…*` alternative string delimiter (NJOY convention, e.g. T23 title).
            j = findnext('*', s, i + 1)
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
    nplot::Int           # card 1, 4th int — plot tape unit (0 = no plot)
    mat::Int; nqa::Int; mts::Vector{Int}
end

function parse_heatr(mc::ModuleCall)::HeatrParams
    cards = mc.raw_cards
    isempty(cards) && return HeatrParams(0,0,0,0,0,0,Int[])
    nendf = abs(_fint(cards[1], 1)); nin = abs(_fint(cards[1], 2))
    nout = length(cards[1]) >= 3 ? abs(_fint(cards[1], 3)) : 0
    nplot = length(cards[1]) >= 4 ? abs(_fint(cards[1], 4)) : 0
    mat = length(cards) >= 2 ? _fint(cards[2], 1) : 0
    nqa = length(cards) >= 2 ? _fint(cards[2], 2; default=0) : 0
    mts = Int[]
    if length(cards) >= 3 && nqa > 0
        for t in cards[3]; startswith(t, "'") && continue; push!(mts, _parse_int_token(t)); end
    end
    HeatrParams(nendf, nin, nout, nplot, mat, nqa, mts)
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
    # Card 1 units. Fortran acer.f90 reads 5 slots: nendf, npend, ngend, nace, ndir.
    # For iopt=1: nendf=ENDF input, npend=PENDF input, nace=ACE out, ndir=xsdir out.
    # For iopt=7: npend=input ACE, ngend=output viewr plot, nace=output ACE, ndir=summary.
    nendf::Int; npendf::Int; ngend::Int; nace::Int; ndir::Int
    mat::Int; iopt::Int; iprint::Int; itype::Int; suffix::String
    temp::Float64; title::String
    nplot::Int  # iopt=7: plot option (-1 = no plot)
end

function parse_acer(mc::ModuleCall)::AcerParams
    cards = mc.raw_cards
    isempty(cards) && return AcerParams(0,0,0,0,0, 0,1,0,1,"80c", 300.0, "", 0)

    # Card 1: nendf npend ngend nace ndir
    c1 = cards[1]
    nendf  = length(c1) >= 1 ? abs(_fint(c1, 1)) : 0
    npendf = length(c1) >= 2 ? abs(_fint(c1, 2)) : 0
    ngend  = length(c1) >= 3 ? abs(_fint(c1, 3)) : 0
    nace   = length(c1) >= 4 ? abs(_fint(c1, 4)) : 0
    ndir   = length(c1) >= 5 ? abs(_fint(c1, 5)) : 0

    # Card 2: iopt iprint itype suff ...
    iopt   = length(cards) >= 2 ? _fint(cards[2], 1; default=1) : 1
    iprint = length(cards) >= 2 && length(cards[2]) >= 2 ? _fint(cards[2], 2; default=1) : 1
    itype  = length(cards) >= 2 && length(cards[2]) >= 3 ? _fint(cards[2], 3; default=1) : 1
    # Suffix appears in card 2 position 4 as a decimal like `.10` or `.80` — read as string
    # to preserve leading dot and trailing zeros, then build a suffix tag like `10c`, `80c`.
    # iopt=1 → suffix `<dd>c`; iopt=7 → output suffix preserved from the ACE header.
    suffix = "80c"
    if length(cards) >= 2 && length(cards[2]) >= 4
        tok = String(strip(cards[2][4]))
        if startswith(tok, ".") && length(tok) >= 3
            suffix = tok[2:end]  # ".10" → "10", ".80" → "80"
            # For iopt=1 (fast neutron/particle), suffix is .XXc/.XXh etc.; charged-particle
            # incident uses .XXa (alpha), .XXh (proton), .XXd (deuteron), .XXt (triton),
            # .XXs (He-3). Default to `c` if not specified by card 2 position 5.
            suffix = length(suffix) == 2 ? suffix * "c" : suffix
        end
    end

    # Card 3: title string (hz, up to 70 chars). Tokenizer keeps surrounding
    # quotes; strip them here. Same pattern used by parse_reconr for its
    # description cards.
    title = length(cards) >= 3 && !isempty(cards[3]) ?
            String(strip(replace(join(cards[3], " "), r"^['\"]|['\"]$" => ""))) : ""

    # Card 4 (iopt=1): matd, tempd, [local], [iprint]
    mat  = 0
    temp = 0.0  # default per Fortran acer.f90 — no broadening when tempd=0
    if iopt == 1 && length(cards) >= 4 && !isempty(cards[4])
        mat  = _fint(cards[4], 1)
        temp = _fnum(cards[4], 2; default=0.0)
    end

    # iopt=7 card 2 position 4 is nplot (e.g. `-1` → no plot). For iopt=7 suffix is not
    # on the card — reused slot is nplot.
    nplot = 0
    if iopt == 7 && length(cards) >= 2 && length(cards[2]) >= 4
        nplot = _fint(cards[2], 4; default=0)
    end

    AcerParams(nendf, npendf, ngend, nace, ndir, mat, iopt, iprint, itype, suffix,
               temp, title, nplot)
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

# mixr: linear-combination of cross sections from N input PENDF/ENDF tapes.
# Card 1: nout, nin1..nin10 (trailing zeros end list)
# Card 2: mtn(1..20)        (output MT list, trailing zeros end list)
# Card 3: (matn, wtn) pairs (one per input tape; matn=0 ends list)
# Card 4: temp              (target temperature)
# Card 5: matd, za, awr     (output material id / ZA / AWR)
# Card 6: 'description'     (66-char text written to MF=1/MT=451)
# Ref: njoy-reference/src/mixr.f90:31-58 (input description),
#      :96-122 (parser),    :150-195 (mat/temp lookup).
struct MixrParams
    nout::Int
    nin::Vector{Int}
    mtn::Vector{Int}
    matn::Vector{Int}
    wtn::Vector{Float64}
    temp::Float64
    matd::Int
    za::Float64
    awr::Float64
    des::String
end

function parse_mixr(mc::ModuleCall)::MixrParams
    cards = mc.raw_cards
    length(cards) >= 6 || error(
        "mixr: need 6 cards (units, MTs, mats+weights, temp, output id, " *
        "description); got $(length(cards))")

    # Card 1: nout, nin1..
    c1 = cards[1]
    nout = abs(_fint(c1, 1))
    nin = Int[]
    for i in 2:length(c1)
        v = abs(_fint(c1, i))
        v == 0 && break
        push!(nin, v)
    end

    # Card 2: MT list
    c2 = cards[2]
    mtn = Int[]
    for i in 1:length(c2)
        v = _fint(c2, i)
        v == 0 && break
        push!(mtn, v)
    end

    # Card 3: (matn, wtn) pairs
    c3 = cards[3]
    matn = Int[]; wtn = Float64[]
    i = 1
    while i + 1 <= length(c3)
        m = _fint(c3, i)
        m == 0 && break
        push!(matn, m); push!(wtn, _fnum(c3, i + 1))
        i += 2
    end
    length(matn) == length(nin) || error(
        "mixr: number of (mat, weight) pairs ($(length(matn))) must equal " *
        "number of input tapes ($(length(nin)))")

    # Card 4: temperature
    temp = _fnum(cards[4], 1)

    # Card 5: matd, za, awr
    c5 = cards[5]
    matd = _fint(c5, 1)
    za   = _fnum(c5, 2)
    awr  = _fnum(c5, 3)

    # Card 6: description (66 chars max). Strip surrounding quotes/asterisks.
    c6 = cards[6]
    des = ""
    if !isempty(c6)
        des = String(strip(c6[1], ['\'', '"', '*', ' ']))
    end
    length(des) > 66 && (des = des[1:66])

    MixrParams(nout, nin, mtn, matn, wtn, temp, matd, za, awr, des)
end

# resxsr: build a CCCC RESXS resonance cross-section binary file from PENDF.
# Card 1: nout
# Card 2: nmat, maxt, nholl, efirst, elast, eps
# Card 3: 'huse' (12-char) ivers
# Card 4: 'holl' (repeat nholl times) — descriptive text
# Card 5: 'hmat' mat unit (repeat nmat times)
# Ref: njoy-reference/src/resxsr.f90:16-40 (input description),
#      :237-248 (parser).
struct ResxsrMaterial
    hmat::String   # 8-char hollerith name
    mat::Int       # ENDF MAT number
    unit::Int      # PENDF input tape unit
end

struct ResxsrParams
    nout::Int
    nmat::Int
    maxt::Int
    nholl::Int
    efirst::Float64
    elast::Float64
    eps::Float64
    huse::String   # 12-char user identifier
    ivers::Int
    holl::Vector{String}
    materials::Vector{ResxsrMaterial}
end

function parse_resxsr(mc::ModuleCall)::ResxsrParams
    cards = mc.raw_cards
    length(cards) >= 3 || error("resxsr: need at least 3 control cards")

    nout = abs(_fint(cards[1], 1))

    c2 = cards[2]
    nmat   = _fint(c2, 1)
    maxt   = _fint(c2, 2)
    nholl  = _fint(c2, 3)
    efirst = _fnum(c2, 4)
    elast  = _fnum(c2, 5)
    eps    = _fnum(c2, 6)

    c3 = cards[3]
    huse  = isempty(c3) ? "" : String(strip(c3[1], ['\'', '"', '*', ' ']))
    ivers = length(c3) >= 2 ? _fint(c3, 2) : 0

    length(cards) >= 3 + nholl + nmat || error(
        "resxsr: insufficient cards — need 3 + nholl ($nholl) + nmat ($nmat) " *
        "= $(3 + nholl + nmat); got $(length(cards))")

    holl = String[]
    idx = 4
    for _ in 1:nholl
        ck = cards[idx]
        text = isempty(ck) ? "" : String(strip(ck[1], ['\'', '"', '*', ' ']))
        push!(holl, text)
        idx += 1
    end

    materials = ResxsrMaterial[]
    for _ in 1:nmat
        ck = cards[idx]
        hmat = String(strip(ck[1], ['\'', '"', '*', ' ']))
        mat  = _fint(ck, 2)
        unit = abs(_fint(ck, 3))
        push!(materials, ResxsrMaterial(hmat, mat, unit))
        idx += 1
    end

    ResxsrParams(nout, nmat, maxt, nholl, efirst, elast, eps,
                  huse, ivers, holl, materials)
end

# powr: produce input for EPRI-CELL (GAMTAP, LIBRAR) and EPRI-CPM (CLIB).
# Card 1: ngendf nout
# Card 2: lib iprint iclaps        (lib: 1=fast/GAMTAP, 2=thermal/LIBRAR, 3=cpm/CLIB)
# Cards 3+ (lib=1): per material: (matd, rtemp, iff, nsgz, izref) + word(16) +
#                   fsn(40)   [or, if matd<0, abs(ngnd) values via card 6].
#                   Terminate with `matd=0/`.
# Cards 3+ (lib=2/3): structured differently (see powr.f90:111-232). Phase A
# keeps those raw_cards for later phases.
# Ref: njoy-reference/src/powr.f90:63-296.
struct PowrFastMaterial
    matd::Int             # ENDF MAT (positive)
    iread::Int            # 0=normal, 1=read absorption directly (matd<0 in input)
    rtemp::Float64        # reference temperature (K)
    iff::Int              # f-factor option (0/1)
    nsgz::Int             # max sigma-zeros (0=all)
    izref::Int            # reference sigma-zero index for elastic
    word::String          # 16-char description
    fsn::String           # 40-char fission spectrum title (unused for non-fissile)
    abs_in::Vector{Float64}  # length-ngnd absorption values when iread=1
end

struct PowrParams
    ngendf::Int
    nout::Int
    lib::Int               # 1=fast, 2=thermal, 3=cpm
    iprint::Int
    iclaps::Int
    fast_mats::Vector{PowrFastMaterial}    # populated only when lib=1
    raw_cards::Vector{Vector{String}}      # full deck for phases C/D
end

function parse_powr(mc::ModuleCall)::PowrParams
    cards = mc.raw_cards
    length(cards) >= 2 || error(
        "powr: need at least 2 control cards (ngendf nout / lib iprint iclaps)")

    ngendf = abs(_fint(cards[1], 1))
    nout   = abs(_fint(cards[1], 2))

    c2     = cards[2]
    lib    = _fint(c2, 1)
    iprint = length(c2) >= 2 ? _fint(c2, 2) : 0
    iclaps = length(c2) >= 3 ? _fint(c2, 3) : 0

    lib in (1, 2, 3) ||
        error("powr: lib must be 1 (fast), 2 (thermal), or 3 (cpm); got $lib")

    fast_mats = PowrFastMaterial[]
    if lib == 1
        ngnd = 68    # Fortran fast() line 410.
        idx  = 3
        while idx <= length(cards)
            ck = cards[idx]
            matd_signed = _fint(ck, 1)
            matd_signed == 0 && break
            iread = matd_signed < 0 ? 1 : 0
            matd  = abs(matd_signed)
            rtemp = length(ck) >= 2 ? _fnum(ck, 2)        : 300.0
            iff   = length(ck) >= 3 ? _fint(ck, 3)        : 1
            nsgz  = length(ck) >= 4 ? _fint(ck, 4)        : 0
            izref = length(ck) >= 5 ? _fint(ck, 5)        : 1
            idx  += 1

            word = ""
            fsn  = ""
            abs_in = Float64[]
            if iread == 0
                idx <= length(cards) || error(
                    "powr lib=1: missing description card after matd=$matd")
                word = isempty(cards[idx]) ? "" :
                       String(strip(cards[idx][1], ['\'', '"', '*', ' ']))
                idx += 1
                idx <= length(cards) || error(
                    "powr lib=1: missing fission spectrum title card after matd=$matd")
                fsn = isempty(cards[idx]) ? "" :
                      String(strip(cards[idx][1], ['\'', '"', '*', ' ']))
                idx += 1
            else
                # iread=1: read ngnd absorption values directly (one card)
                idx <= length(cards) || error(
                    "powr lib=1: missing absorption card after matd=$matd (iread=1)")
                abs_in = [_fnum(cards[idx], k) for k in 1:min(length(cards[idx]), ngnd)]
                length(abs_in) < ngnd && append!(abs_in, zeros(ngnd - length(abs_in)))
                idx += 1
            end

            push!(fast_mats, PowrFastMaterial(matd, iread, rtemp, iff, nsgz,
                                                izref, word, fsn, abs_in))
        end
    end

    PowrParams(ngendf, nout, lib, iprint, iclaps, fast_mats, cards)
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
    # (Fortran groupr uses the first block's MT list for all temperatures).
    #
    # Ref: njoy-reference/src/groupr.f90:610-634. Fortran inits mtdp=-1000
    # before the list-directed `read(nsysi,*) mfd,mtdp,strng`, so a card with
    # only mfd (e.g. `3 /`) leaves mtdp at -1000, which triggers auto-expand
    # via `nextr` (label 382, iauto=1). We preserve the sentinel verbatim;
    # groupr_module expands it against the PENDF MF=3 MT set.
    mt_list = Tuple{Int,Int,String}[]
    while ci <= length(cards)
        mfd = _fint(cards[ci], 1; default=0)
        mfd == 0 && break
        mtd = length(cards[ci]) >= 2 ? _parse_int_token(cards[ci][2]) : -1000
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

"""
    LeaprParams

Full leapr input-deck contract. Mirrors the Fortran `leapr` reads at
`njoy-reference/src/leapr.f90:230-400`.

Per-temperature vectors have length `ntempr`. When the Fortran deck uses a
negative temperature (`temp < 0`) to signal "reuse previous", we copy the
prior values into the matching slot — the struct is always fully populated.
"""
struct LeaprParams
    # Card 1: output ENDF tape unit
    nout::Int
    # Card 2: title (80-char text; quotes stripped)
    title::String
    # Card 3: run control
    ntempr::Int           # number of temperatures (default 1)
    iprint::Int           # print flag (default 1)
    nphon::Int            # phonon-expansion order (default 100)
    # Card 4: ENDF metadata
    mat::Int
    za::Float64
    isabt::Int            # default 0 (isabt=1 unsupported; Fortran warns)
    ilog::Int             # default 0 (1 → log10(S))
    smin::Float64         # default 1.0e-75
    # Card 5: principal scatterer
    awr::Float64
    spr::Float64
    npr::Int
    iel::Int              # coherent-elastic option; default 0
    ncold::Int            # cold moderator option; default 0
    nsk::Int              # s(kappa) option; default 0
    # Card 6: secondary scatterer (defaults for nss=0)
    nss::Int
    b7::Float64
    aws::Float64
    sps::Float64
    mss::Int
    # Card 7: alpha/beta control
    nalpha::Int
    nbeta::Int
    lat::Int              # default 0
    # Cards 8-9: grids
    alpha::Vector{Float64}
    beta::Vector{Float64}
    # Per-temperature data (len = ntempr; mirrors Fortran copy-on-reread)
    temperatures::Vector{Float64}
    delta1::Vector{Float64}
    ni::Vector{Int}
    p1_dos::Vector{Vector{Float64}}
    twt::Vector{Float64}
    c::Vector{Float64}
    tbeta::Vector{Float64}
    nd::Vector{Int}
    bdel::Vector{Vector{Float64}}   # oscillator energies (empty if nd=0)
    adel::Vector{Vector{Float64}}   # oscillator weights (empty if nd=0)
    nka::Vector{Int}
    dka::Vector{Float64}
    ska::Vector{Vector{Float64}}    # pair-correlation s(kappa) (empty if nsk=0 && ncold=0)
    cfrac::Vector{Float64}          # coherent fraction (empty unless nsk>0)
    # Card 20: MF1/MT451 descriptive lines
    comments::Vector{String}
end

# ---------- leapr deck cursor helpers ----------
# The Fortran `read(nsysi,*)` free-format reader treats `/` as a soft terminator:
# unread variables keep their pre-initialised defaults. The Julia tokeniser
# splits on `/`, so each `cards[i]` is exactly one Fortran read's worth of
# explicit tokens (remaining variables default). For array reads like
# `(alpha(i), i=1, nalpha)` values flow across successive cards until `nalpha`
# are consumed — matching Fortran's "keep reading until count satisfied".

"Consume one card (for a scalar Fortran read). Caller supplies defaults."
function _leapr_scalars!(cards::Vector{Vector{String}}, ci::Base.RefValue{Int})
    ci[] > length(cards) && return String[]
    card = cards[ci[]]
    ci[] += 1
    return card
end

"Consume N tokens across consecutive cards (for array Fortran reads)."
function _leapr_array!(cards::Vector{Vector{String}}, ci::Base.RefValue{Int}, n::Int)
    vals = Float64[]
    while length(vals) < n
        ci[] > length(cards) && error("leapr deck ran out of tokens (needed $n, have $(length(vals)))")
        for t in cards[ci[]]
            (startswith(t, "'") || startswith(t, "\"") || startswith(t, "*")) && continue
            push!(vals, _parse_num(t))
            length(vals) >= n && break
        end
        ci[] += 1
    end
    return vals
end

"Strip surrounding quotes (`'…'`, `\"…\"`, `*…*`). Preserves internal
whitespace (the NJOY convention keeps leading/trailing spaces inside
quotes for proper MF1/MT451 hollerith formatting)."
function _leapr_strip_title(s::AbstractString)
    t = String(strip(s))                      # outer-only strip
    (startswith(t, "'") && endswith(t, "'")) && return String(t[2:end-1])
    (startswith(t, "\"") && endswith(t, "\"")) && return String(t[2:end-1])
    (startswith(t, "*") && endswith(t, "*"))  && return String(t[2:end-1])
    return t
end

function parse_leapr(mc::ModuleCall)::LeaprParams
    cards = mc.raw_cards
    if isempty(cards)
        return LeaprParams(
            0, "", 1, 1, 100, 0, 0.0, 0, 0, 1e-75,
            0.0, 0.0, 0, 0, 0, 0,
            0, 0.0, 0.0, 0.0, 0,
            0, 0, 0,
            Float64[], Float64[],
            Float64[], Float64[], Int[], Vector{Float64}[],
            Float64[], Float64[], Float64[], Int[],
            Vector{Float64}[], Vector{Float64}[],
            Int[], Float64[], Vector{Float64}[], Float64[],
            String[],
        )
    end

    ci = Ref(1)

    # Card 1: nout
    c = _leapr_scalars!(cards, ci)
    nout = abs(_fint(c, 1))

    # Card 2: title (one quoted token, typically)
    c = _leapr_scalars!(cards, ci)
    title = isempty(c) ? "" : _leapr_strip_title(c[1])

    # Card 3: ntempr, iprint, nphon (defaults 1, 1, 100)
    c = _leapr_scalars!(cards, ci)
    ntempr = _fint(c, 1; default=1)
    iprint = _fint(c, 2; default=1)
    nphon  = _fint(c, 3; default=100)

    # Card 4: mat, za, isabt, ilog, smin
    c = _leapr_scalars!(cards, ci)
    mat   = _fint(c, 1)
    za    = _fnum(c, 2)
    isabt = _fint(c, 3; default=0)
    ilog  = _fint(c, 4; default=0)
    smin  = _fnum(c, 5; default=1e-75)

    # Card 5: awr, spr, npr, iel, ncold, nsk
    c = _leapr_scalars!(cards, ci)
    awr   = _fnum(c, 1)
    spr   = _fnum(c, 2)
    npr   = _fint(c, 3)
    iel   = _fint(c, 4; default=0)
    ncold = _fint(c, 5; default=0)
    nsk   = _fint(c, 6; default=0)

    # Card 6: nss, b7, aws, sps, mss
    c = _leapr_scalars!(cards, ci)
    nss = _fint(c, 1; default=0)
    b7  = _fnum(c, 2; default=0.0)
    aws = _fnum(c, 3; default=0.0)
    sps = _fnum(c, 4; default=0.0)
    mss = _fint(c, 5; default=0)

    # Card 7: nalpha, nbeta, lat
    c = _leapr_scalars!(cards, ci)
    nalpha = _fint(c, 1)
    nbeta  = _fint(c, 2)
    lat    = _fint(c, 3; default=0)

    # Cards 8-9: alpha and beta grids
    alpha = _leapr_array!(cards, ci, nalpha)
    beta  = _leapr_array!(cards, ci, nbeta)

    # Per-temperature loop (Fortran `do itemp = 1, ntempr`)
    temps  = Float64[]; delta1s = Float64[]; nis = Int[]
    p1s    = Vector{Vector{Float64}}()
    twts   = Float64[]; cs = Float64[]; tbetas = Float64[]
    nds    = Int[]; bdels = Vector{Vector{Float64}}(); adels = Vector{Vector{Float64}}()
    nkas   = Int[]; dkas = Float64[]; skas = Vector{Vector{Float64}}()
    cfracs = Float64[]

    for itemp in 1:ntempr
        c = _leapr_scalars!(cards, ci)
        temp = _fnum(c, 1)
        push!(temps, abs(temp))

        # Reread cards 11-19 only when itemp==1 OR temp>=0
        if itemp == 1 || temp >= 0
            # Card 11: delta1, ni
            c = _leapr_scalars!(cards, ci)
            delta1 = _fnum(c, 1)
            ni = _fint(c, 2)
            push!(delta1s, delta1); push!(nis, ni)
            # Card 12: p1 DOS (ni values)
            p1 = _leapr_array!(cards, ci, ni)
            push!(p1s, p1)
            # Card 13: twt, c, tbeta
            c = _leapr_scalars!(cards, ci)
            push!(twts,   _fnum(c, 1; default=0.0))
            push!(cs,     _fnum(c, 2; default=0.0))
            push!(tbetas, _fnum(c, 3; default=0.0))
            # Card 14: nd
            c = _leapr_scalars!(cards, ci)
            nd = _fint(c, 1; default=0)
            push!(nds, nd)
            if nd > 0
                # Cards 15-16: oscillator data
                push!(bdels, _leapr_array!(cards, ci, nd))
                push!(adels, _leapr_array!(cards, ci, nd))
            else
                push!(bdels, Float64[]); push!(adels, Float64[])
            end
            # Cards 17-18: s(kappa) if needed
            if nsk > 0 || ncold > 0
                c = _leapr_scalars!(cards, ci)
                nka = _fint(c, 1)
                dka = _fnum(c, 2)
                push!(nkas, nka); push!(dkas, dka)
                push!(skas, _leapr_array!(cards, ci, nka))
            else
                push!(nkas, 0); push!(dkas, 0.0); push!(skas, Float64[])
            end
            # Card 19: cfrac (only if nsk>0)
            if nsk > 0
                c = _leapr_scalars!(cards, ci)
                push!(cfracs, _fnum(c, 1; default=0.0))
            end
        else
            # Reuse prior temp's data
            push!(delta1s, delta1s[end]); push!(nis, nis[end]); push!(p1s, p1s[end])
            push!(twts, twts[end]); push!(cs, cs[end]); push!(tbetas, tbetas[end])
            push!(nds, nds[end]); push!(bdels, bdels[end]); push!(adels, adels[end])
            push!(nkas, nkas[end]); push!(dkas, dkas[end]); push!(skas, skas[end])
            nsk > 0 && push!(cfracs, cfracs[end])
        end
    end

    # Card 20: free-form comment block (typically quoted strings) until EOF/empty
    comments = String[]
    while ci[] <= length(cards)
        card = cards[ci[]]
        ci[] += 1
        isempty(card) && continue
        push!(comments, _leapr_strip_title(join(card, " ")))
    end

    return LeaprParams(
        nout, title, ntempr, iprint, nphon,
        mat, za, isabt, ilog, smin,
        awr, spr, npr, iel, ncold, nsk,
        nss, b7, aws, sps, mss,
        nalpha, nbeta, lat,
        alpha, beta,
        temps, delta1s, nis, p1s,
        twts, cs, tbetas, nds,
        bdels, adels,
        nkas, dkas, skas, cfracs,
        comments,
    )
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

"""
One mat/mt/mat1/mt1 case from a covr card 4. Matches the
`(imat,imt,imat1,imt1)` quadruple in covr.f90 line 259. Defaults match
Fortran lines 256-258: imt=imat1=imt1=0 if omitted on the card.
"""
struct CovrCase
    mat::Int
    mt::Int       # 0 = all MTs for this MAT; negative = strip-list semantics
    mat1::Int     # 0 = use mat
    mt1::Int      # 0 = use mt
end

"""
Full input parameters for a single Fortran covr invocation. Mirrors the
free-format card sequence in covr.f90 lines 168-260:

  card 1:  nin, nout, nplot
  cards 2/2'/2a/3a   (plot mode, nout ≤ 0)
    card 2:   icolor (0=mono, 1=color, 2=color+custom-tlev)
    card 2':  nlev, tlev[1..nlev]    (only when icolor=2)
    card 2a:  epmin                  (multiplied by rdn=0.999998 like Fortran)
    card 3a:  irelco, ncase, noleg, nstart, ndiv
  cards 2b/3b/3c     (library mode, nout > 0)
    card 2b:  matype, ncase
    card 3b:  hlibid (≤ 6 chars, single quoted token)
    card 3c:  hdescr (≤ 21 chars, single quoted token)
  card 4:    mat, mt, mat1, mt1   (repeated ncase times)

Fortran defaults retained verbatim:
  icolor=0, nlev=6, tlev=(.001,.1,.2,.3,.6,1.0), epmin=0,
  irelco=1, ncase=1, noleg=0, nstart=1, ndiv=1, matype=3.
"""
struct CovrParams
    nin::Int
    nout::Int
    nplot::Int
    # Plot-mode (used when nout ≤ 0)
    icolor::Int
    tlev::Vector{Float64}    # length nlev, last entry must equal 1.0
    epmin::Float64
    irelco::Int
    ncase::Int
    noleg::Int
    nstart::Int
    ndiv::Int
    # Library-mode (used when nout > 0)
    matype::Int              # 3=cov, 4=corr
    hlibid::String           # ≤ 6 chars
    hdescr::String           # ≤ 21 chars
    # Cases (length = ncase; valid in either mode)
    cases::Vector{CovrCase}
end

# Default tlev matches covr.f90 line 161 (drop trailing zero placeholders).
const _COVR_DEFAULT_TLEV =
    Float64[0.001, 0.1, 0.2, 0.3, 0.6, 1.0]
const _COVR_RDN = 0.999998     # covr.f90 line 164

# Split a raw covr deck line into N cards, one per `/` terminator.
# Tokens on the line precede their `/`; consecutive `/` (e.g. `//`) yield
# trailing empty cards. Tokens following the last `/` are discarded
# (Fortran free-format ignores trailing comments after the terminator).
function _covr_split_card_line(line::AbstractString)::Vector{Vector{String}}
    out = Vector{String}[]
    s = line
    while true
        si = _find_slash_outside_quotes(s)
        si == 0 && break
        push!(out, _tokenise_card(s[1:si-1]))
        s = s[si+1:end]
    end
    # No `/` at all: treat the whole line as one card (continuation form).
    isempty(out) && !isempty(strip(line)) && push!(out, _tokenise_card(line))
    out
end

# Flatten covr raw_lines into a flat sequence of cards. Each `/` produces
# one card (possibly empty). Comment-only lines (already filtered upstream)
# do not appear here. Quoted strings span one card.
function _covr_collect_cards(raw_lines::Vector{String})::Vector{Vector{String}}
    cards = Vector{String}[]
    for line in raw_lines
        for c in _covr_split_card_line(line)
            push!(cards, c)
        end
    end
    cards
end

# Strip surrounding ' or " or * delimiters from a single token.
function _covr_unquote(tok::AbstractString)::String
    t = String(tok)
    isempty(t) && return ""
    if (t[1] == '\'' && t[end] == '\'' && length(t) >= 2) ||
       (t[1] == '"'  && t[end] == '"'  && length(t) >= 2) ||
       (t[1] == '*'  && t[end] == '*'  && length(t) >= 2)
        return t[2:end-1]
    end
    t
end

function parse_covr(mc::ModuleCall)::CovrParams
    cards = isempty(mc.raw_lines) ? mc.raw_cards : _covr_collect_cards(mc.raw_lines)
    if isempty(cards)
        return CovrParams(0, 0, 0,
                          0, copy(_COVR_DEFAULT_TLEV), 0.0,
                          1, 1, 0, 1, 1,
                          3, "", "", CovrCase[])
    end

    # Card 1: nin, nout, nplot.  Fortran reads `nout=0, nplot=0` defaults,
    # then `read(nsysi,*) nin,nout,nplot`.
    c1 = cards[1]
    nin   = abs(_fint(c1, 1))
    nout  = _fint(c1, 2; default=0)   # signed in Fortran; sign chooses ASCII/binary
    nplot = _fint(c1, 3; default=0)
    idx = 2

    icolor = 0
    tlev   = copy(_COVR_DEFAULT_TLEV)
    epmin  = 0.0
    irelco = 1
    ncase  = 1
    noleg  = 0
    nstart = 1
    ndiv   = 1
    matype = 3
    hlibid = ""
    hdescr = ""

    if nout <= 0
        # ----- plot mode -----
        # card 2: icolor
        if idx <= length(cards)
            icolor = _fint(cards[idx], 1; default=0)
            idx += 1
        end
        # card 2' (only when icolor=2): nlev, tlev[1..nlev]
        if icolor == 2 && idx <= length(cards)
            c = cards[idx]
            nlev = _fint(c, 1; default=6)
            nlev > 9 && error("covr: nlev=$nlev exceeds nlevmx=9 (covr.f90:154)")
            tlev = Float64[_fnum(c, 1+i; default=0.0) for i in 1:nlev]
            for i in 2:nlev
                tlev[i] > tlev[i-1] ||
                    error("covr: tlev must be strictly increasing (covr.f90:206-211)")
            end
            tlev[end] != 1.0 && (tlev[end] = 1.0)   # covr.f90:213-218
            idx += 1
        end
        # card 2a: epmin
        if idx <= length(cards)
            epmin = _fnum(cards[idx], 1; default=0.0) * _COVR_RDN
            idx += 1
        end
        # card 3a: irelco, ncase, noleg, nstart, ndiv
        if idx <= length(cards)
            c = cards[idx]
            irelco = _fint(c, 1; default=1)
            ncase  = _fint(c, 2; default=1)
            noleg  = _fint(c, 3; default=0)
            nstart = _fint(c, 4; default=1)
            ndiv   = _fint(c, 5; default=1)
            ndiv == 0 && (ndiv = 1)
            idx += 1
        end
    else
        # ----- library mode -----
        # card 2b: matype, ncase
        if idx <= length(cards)
            c = cards[idx]
            matype = _fint(c, 1; default=3)
            matype != 4 && (matype = 3)
            ncase  = _fint(c, 2; default=1)
            ncase <= 0 && (ncase = 1)
            idx += 1
        end
        # card 3b: hlibid (one quoted token)
        if idx <= length(cards) && !isempty(cards[idx])
            hlibid = _covr_unquote(cards[idx][1])
            length(hlibid) > 6 && (hlibid = hlibid[1:6])
            idx += 1
        end
        # card 3c: hdescr (one quoted token)
        if idx <= length(cards) && !isempty(cards[idx])
            hdescr = _covr_unquote(cards[idx][1])
            length(hdescr) > 21 && (hdescr = hdescr[1:21])
            idx += 1
        end
    end

    # ncase × card 4
    cases = CovrCase[]
    for k in 1:ncase
        if idx + k - 1 > length(cards)
            error("covr: only $(idx-1+k-1)/$ncase case cards present (covr.f90:255-260)")
        end
        c = cards[idx + k - 1]
        m   = _fint(c, 1; default=0)
        mt  = _fint(c, 2; default=0)
        m1  = _fint(c, 3; default=0)
        mt1 = _fint(c, 4; default=0)
        push!(cases, CovrCase(m, mt, m1, mt1))
    end

    CovrParams(abs(nin), nout, nplot,
               icolor, tlev, epmin,
               irelco, ncase, noleg, nstart, ndiv,
               matype, hlibid, hdescr, cases)
end

struct PlotrParams
    nplt::Int     # output plot-command tape (card 1, first int)
end

function parse_plotr(mc::ModuleCall)::PlotrParams
    cards = mc.raw_cards
    isempty(cards) && return PlotrParams(0)
    PlotrParams(abs(_fint(cards[1], 1)))
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
# wimsr — WIMS-D/E reactor library generator
# Ref: njoy-reference/src/wimsr.f90:48-307 (subroutine wimsr, card layout).
# Card layout (with conditionals):
#   1: ngendf nout                        (input GENDF unit, output WIMS unit)
#   2: iprint iverw igroup                (verbosity; WIMS version 4|5; group struct)
#   2a: ngnd nfg nrg igref                (ONLY if igroup==9; user group structure)
#   3: mat nfid rdfid iburn               (ENDF MAT, WIMS ID, resonance ID, burnup flag)
#   4: ntemp nsigz sgref ires sigp mti mtc ip1opt inorf isof ifprod jp1
#                                         (12 fields with defaults, free-form)
#   5: ntis efiss                         (ONLY if iburn>0)
#   6a-c: burnup chain pairs              (ONLY if iburn>0)
#   7: glam(1:nrg)                        (Goldstein lambdas)
#   8: p1flx(1:jp1)                       (ONLY if jp1>0)
# =========================================================================

struct WimsrBurnupCard
    ntis::Int
    efiss::Float64
    pairs::Vector{Tuple{Int,Float64}}  # (ident, yield/decay)
end

struct WimsrParams
    ngendf::Int                 # card 1: GENDF input unit (always abs() — sign is binary flag)
    nout::Int                   # card 1: WIMS output unit (positive = ASCII)
    iprint::Int                 # card 2 default 0
    iverw::Int                  # card 2 default 4 (WIMS-D)
    igroup::Int                 # card 2 default 0 (69-group default)
    ngnd::Int                   # card 2a (only if igroup==9; else 69)
    nfg::Int                    # card 2a (only if igroup==9; else 14)
    nrg::Int                    # card 2a (only if igroup==9; else 13)
    igref::Int                  # card 2a (only if igroup==9; else 0)
    mat::Int                    # card 3
    nfid::Int                   # card 3 (set from nint(rdfid))
    rdfid::Float64              # card 3
    iburn::Int                  # card 3 default 0
    ntemp::Int                  # card 4 default 0 (means "use all on tape")
    nsigz::Int                  # card 4 default 0
    sgref::Float64              # card 4 default 1e10 (infinite dilution)
    ires::Int                   # card 4 default 0 (no resonance tables)
    sigp::Float64               # card 4 default 0
    mti::Int                    # card 4 thermal inelastic MT (e.g. 221)
    mtc::Int                    # card 4 thermal coherent MT (0 if absent)
    ip1opt::Int                 # card 4 default 1 (1=skip P1, 0=include)
    inorf::Int                  # card 4 default 0
    isof::Int                   # card 4 default 0 (skip fission spectrum)
    ifprod::Int                 # card 4 default 0
    jp1::Int                    # card 4 default 0 (no current spectrum)
    burnup::Union{Nothing,WimsrBurnupCard}    # cards 5/6 if iburn>0
    glam::Vector{Float64}       # card 7 (length nrg)
    p1flx::Vector{Float64}      # card 8 (length jp1; empty if jp1==0)
end

function parse_wimsr(mc::ModuleCall)::WimsrParams
    cards = mc.raw_cards
    isempty(cards) && error("wimsr: empty card list")

    # Card 1
    ngendf = abs(_fint(cards[1], 1))
    nout   = abs(_fint(cards[1], 2))

    # Card 2 (defaults from wimsr.f90:170-172)
    iprint = length(cards) >= 2 ? _fint(cards[2], 1; default=0) : 0
    iverw  = length(cards) >= 2 ? _fint(cards[2], 2; default=4) : 4
    igroup = length(cards) >= 2 ? _fint(cards[2], 3; default=0) : 0

    # Card 2a — only if igroup==9 (wimsr.f90:186-189). Default 69-group structure.
    next_ci = 3
    ngnd, nfg, nrg, igref = 69, 14, 13, 0
    if igroup == 9
        c2a = cards[next_ci]
        ngnd  = _fint(c2a, 1; default=69)
        nfg   = _fint(c2a, 2; default=14)
        nrg   = _fint(c2a, 3; default=13)
        igref = _fint(c2a, 4; default=0)
        next_ci += 1
    end

    # Card 3 (wimsr.f90:202-205)
    c3 = cards[next_ci]
    mat   = _fint(c3, 1; default=0)
    nfid  = _fint(c3, 2; default=0)
    rdfid = _fnum(c3, 3; default=0.0)
    iburn = _fint(c3, 4; default=0)
    # Fortran wimsr.f90:204 unconditionally does `nfid=nint(rdfid)` — the
    # deck's second field is read into `nfid` first but then *overwritten*
    # by the third field. T11 deck `1050 1 1050.` ends up with nfid=1050,
    # not 1. Fall back to nfid (float'd) only if rdfid is absent (==0).
    if rdfid != 0.0
        nfid = round(Int, rdfid)
    elseif nfid != 0
        rdfid = Float64(nfid)
    end
    next_ci += 1

    # Card 4 (wimsr.f90:217-218) — 12 fields with defaults
    c4 = cards[next_ci]
    ntemp  = _fint(c4, 1;  default=0)
    nsigz  = _fint(c4, 2;  default=0)
    sgref  = _fnum(c4, 3;  default=1e10)
    ires   = _fint(c4, 4;  default=0)
    sigp   = _fnum(c4, 5;  default=0.0)
    mti    = _fint(c4, 6;  default=0)
    mtc    = _fint(c4, 7;  default=0)
    ip1opt = _fint(c4, 8;  default=1)
    inorf  = _fint(c4, 9;  default=0)
    isof   = _fint(c4, 10; default=0)
    ifprod = _fint(c4, 11; default=0)
    jp1    = _fint(c4, 12; default=0)
    next_ci += 1

    # Cards 5/6: burnup data — only if iburn>0 (wimsr.f90:242-255)
    burnup = nothing
    if iburn > 0
        c5 = cards[next_ci]
        ntis = _fint(c5, 1; default=0)
        efiss = _fnum(c5, 2; default=0.0)
        next_ci += 1
        pairs = Tuple{Int,Float64}[]
        # Fortran reads ntis (ident, yield) pairs across one or more lines
        for _ in 1:ntis
            next_ci > length(cards) && break
            row = cards[next_ci]
            for j in 1:2:length(row)
                j+1 > length(row) && break
                push!(pairs, (_fint(row, j), _fnum(row, j+1)))
            end
            next_ci += 1
        end
        burnup = WimsrBurnupCard(ntis, efiss, pairs)
    end

    # Card 7: glam (nrg lambdas)
    glam = Float64[]
    if next_ci <= length(cards)
        for tok in cards[next_ci]
            startswith(tok, "'") && continue
            push!(glam, _parse_num(tok))
            length(glam) >= nrg && break
        end
        next_ci += 1
    end

    # Card 8: p1flx (only if jp1>0)
    p1flx = Float64[]
    if jp1 > 0 && next_ci <= length(cards)
        for tok in cards[next_ci]
            startswith(tok, "'") && continue
            push!(p1flx, _parse_num(tok))
            length(p1flx) >= jp1 && break
        end
    end

    WimsrParams(ngendf, nout, iprint, iverw, igroup,
                ngnd, nfg, nrg, igref, mat, nfid, rdfid, iburn,
                ntemp, nsigz, sgref, ires, sigp, mti, mtc,
                ip1opt, inorf, isof, ifprod, jp1,
                burnup, glam, p1flx)
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
