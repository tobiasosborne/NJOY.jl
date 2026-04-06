#=
    NJOY plot tape parser — reads Fortran free-format plot cards.
    Handles trailing / terminators, missing values with -99 sentinels,
    Fortran-style exponents (1.0+3 -> 1.0e+3), quoted strings,
    comments after /, and blank lines.

    Used by viewr_render! to parse the plot tape produced by plotr/covr/dtfr.
=#

"""
    read_plot_reals(line) -> Vector{Float64}

Parse Fortran free-format real values from a slash-terminated line.
Returns an empty vector for blank lines or bare `/`.

Handles:
- Trailing `/` terminators (everything after `/` is ignored)
- Fortran-style exponents: `1.0+3` -> `1.0e+3`, `2.5-4` -> `2.5e-4`
- Missing values become -99.0 sentinel
- Standard floating point: `1.23`, `-4.5e6`, `.001`
"""
function read_plot_reals(line::AbstractString)
    s = strip(line)

    # empty or bare slash
    if isempty(s) || s == "/"
        return Float64[]
    end

    # truncate at first unquoted /
    in_quote = false
    slash_pos = 0
    for i in eachindex(s)
        c = s[i]
        if c == '\''
            in_quote = !in_quote
        elseif c == '/' && !in_quote
            slash_pos = i
            break
        end
    end
    if slash_pos > 0
        s = s[1:prevind(s, slash_pos)]
    end
    s = strip(s)
    isempty(s) && return Float64[]

    # remove any quoted strings (they are not numeric)
    s = replace(s, r"'[^']*'" => "")
    s = strip(s)
    isempty(s) && return Float64[]

    vals = Float64[]
    for tok in split(s)
        isempty(tok) && continue
        # Fix Fortran-style exponents: 1.0+3 -> 1.0e+3, 2.5-4 -> 2.5e-4
        # Pattern: a digit or dot followed by + or - then digits, without 'e'
        fixed = replace(String(tok), r"([0-9\.])([+-])(\d)" => s"\1e\2\3")
        try
            push!(vals, parse(Float64, fixed))
        catch
            # skip unparseable tokens
        end
    end
    return vals
end

"""
    read_plot_string(line) -> (str, len)

Extract a string from a plot tape card.
Handles quoted strings (`'text'`), unquoted text, and trailing `/`.
Returns the string and its effective length (trailing blanks stripped).
"""
function read_plot_string(line::AbstractString)
    s = strip(line)

    # remove trailing / and anything after it
    slash_pos = 0
    in_quote = false
    for i in eachindex(s)
        c = s[i]
        if c == '\''
            in_quote = !in_quote
        elseif c == '/' && !in_quote
            slash_pos = i
            break
        end
    end
    if slash_pos > 0
        s = s[1:prevind(s, slash_pos)]
    end
    s = strip(s)

    # extract quoted string
    m = match(r"^'(.*)'$", s)
    if m !== nothing
        txt = m.captures[1]
    else
        m2 = match(r"'([^']*)'", s)
        if m2 !== nothing
            txt = m2.captures[1]
        else
            txt = s
        end
    end

    # compute length = last non-blank character position
    len = 0
    for i in 1:length(txt)
        if txt[i] != ' '
            len = i
        end
    end

    return (txt, len)
end

"""
    read_plot_ints(line) -> Vector{Int}

Parse integer values from a plot tape card.
Uses read_plot_reals and rounds to Int.
"""
function read_plot_ints(line::AbstractString)
    vals = read_plot_reals(line)
    return [round(Int, v) for v in vals]
end

# ---- Fortran-compatible card readers for the viewr main loop ----

"""
    _read_card_reals(lines, li, n, defaults) -> (values, li_new)

Read one card from the plot tape at line index `li`, parse up to `n` real
values with given defaults. Returns the values vector and the next line index.
"""
function _read_card_reals(lines::Vector{String}, li::Int, n::Int,
                          defaults::Vector{Float64})
    vals = copy(defaults)
    if li <= length(lines)
        parsed = read_plot_reals(lines[li])
        for i in 1:min(length(parsed), n)
            vals[i] = parsed[i]
        end
        li += 1
    end
    return (vals, li)
end

"""
    _read_card_string(lines, li) -> (str, len, li_new)

Read a string card from the plot tape. Returns the string, its length,
and the next line index.
"""
function _read_card_string(lines::Vector{String}, li::Int)
    if li <= length(lines)
        str, len = read_plot_string(lines[li])
        return (str, len, li + 1)
    end
    return ("", 0, li)
end
