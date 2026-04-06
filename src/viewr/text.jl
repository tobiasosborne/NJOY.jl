#=
    Text rendering engine with DISSPLA markup for the NJOY graph engine.
    Faithful port of graph.f90 subroutines: charw, stripv, iget, rget,
    ssum, txtlen, dchr, text3, text2.

    All format strings match Fortran exactly for bit-identical PS output.
=#

using Printf

# ---- Character width lookup ----

"""
    charw(ch::Char, itable::Int) -> Float64

Return the width of character `ch` in font `itable` (1=Times-Roman,
2=Helvetica, 3=Symbol).  Width is a fraction of character height.
Matches graph.f90 `charw` (lines 2850-2863).
"""
function charw(ch::Char, itable::Int)
    i = 1 + Int(ch)
    return CW[i, itable]
end

# ---- Strip leading/trailing blanks ----

"""
    stripv(str::String) -> (stripped::String, length::Int)

Strip leading and trailing non-printable / blank characters, keeping
only ASCII 33-126 (printable non-space).  Returns the stripped string
and its length.
Matches graph.f90 `stripv` (lines 2865-2896).
"""
function stripv(str::String)
    n = length(str)
    j = 0  # first printable position
    k = 0  # last printable position
    for i in 1:n
        ic = Int(str[i])
        if ic > 32 && ic < 127
            if j == 0
                j = i
            end
            k = i
        end
    end
    if j == 0
        return ("", 0)
    end
    len = k - j + 1
    return (str[j:k], len)
end

# ---- Parse integer from string ----

"""
    iget(str::String, n::Int) -> Int

Decode an integer from the first `n` characters of `str`.
Handles optional leading minus sign.
Matches graph.f90 `iget` (lines 2033-2059).
"""
function iget(str::String, n::Int)
    digits = "0123456789"
    neg = 0
    ival = 0
    for i in 1:n
        c = str[i]
        j = findfirst(==(c), digits)
        if j !== nothing
            ival = 10 * ival + j - 1
        elseif c == '-'
            neg = 1
        end
    end
    if neg == 1
        ival = -ival
    end
    return ival
end

# ---- Parse real from string ----

"""
    rget(str::String, n::Int) -> Float64

Decode a real number from the first `n` characters of `str`.
Handles optional sign and decimal point; no exponent support.
Matches graph.f90 `rget` (lines 2061-2095).
"""
function rget(str::String, n::Int)
    digits = "0123456789"
    neg = 0
    ival = 0
    ipwr = 0
    for i in 1:n
        c = str[i]
        j = findfirst(==(c), digits)
        if j !== nothing
            ival = 10 * ival + j - 1
            if ipwr != 0
                ipwr += 1
            end
        elseif c == '.'
            ipwr = 1
        elseif c == '-'
            neg = 1
        end
    end
    val = Float64(ival)
    if ipwr != 0
        val = val * 10.0^(1 - ipwr)
    end
    if neg == 1
        val = -val
    end
    return val
end

# ---- Sum character widths ----

"""
    ssum(str::String, l::Int, ifont::Int, ialf::Int, scale::Float64) -> Float64

Sum the widths of the first `l` characters in `str` using font `ifont`,
scaled by `scale`.  Used by `txtlen` to measure text width.
Matches graph.f90 `ssum` (lines 2013-2031).
"""
function ssum(str::String, l::Int, ifont::Int, ialf::Int, scale::Float64)
    itabl = ifont
    s = 0.0
    for i in 1:l
        s += scale * charw(str[i], itabl)
    end
    return s
end

# ---- Measure text width including markup ----

"""
    txtlen(gs::GraphState, txt::String, ntxt::Int, ifont::Int, ht::Float64) -> Float64

Measure the total rendered width of text string `txt` (length `ntxt`),
including DISSPLA markup processing.  Returns width in the same units
as `ht`.  Does not produce any PostScript output.
Matches graph.f90 `txtlen` (lines 1838-2011).
"""
function txtlen(gs::GraphState, txt::String, ntxt::Int, ifont::Int, ht::Float64)
    siz = 0.0
    cbase = 0.0
    celev = 0.0
    icase = 0
    ialf = 0
    scale = 1.0
    imode = 0
    iskip = 0
    jfont = ifont
    i = 1
    j = 0
    l = 0
    temp = Vector{Char}(undef, 80)
    comm = ' '

    while i <= ntxt
        c = txt[i]

        if imode == 1
            # ---- Command argument accumulation ----
            if (c >= '0' && c <= '9') || c == '.' || c == '-'
                j += 1
                temp[j] = c
            elseif c == 'X'
                j += 1
                temp[j] = c
            else
                if j == 0
                    temp[1] = 'D'
                    j = 1
                end
                if comm == 'E'
                    if temp[1] == 'D'
                        temp[1] = '.'
                        temp[2] = '5'
                        j = 2
                    end
                    if temp[1] == 'X'
                        celev = 0.0
                    else
                        y = rget(String(temp[1:j]), j)
                        if y > 0
                            celev = ht * y
                        else
                            celev = -ht * y
                        end
                    end
                elseif comm == 'L'
                    if temp[1] == 'D'
                        temp[1] = '.'
                        temp[2] = '5'
                        j = 2
                    end
                    if temp[1] == 'X'
                        celev = 0.0
                    else
                        y = rget(String(temp[1:j]), j)
                        if y > 0
                            celev = -ht * y
                        else
                            celev = ht * y
                        end
                    end
                elseif comm == 'H'
                    if temp[1] == 'D'
                        temp[1] = '.'
                        temp[2] = '7'
                        j = 2
                    end
                    if temp[1] == 'X'
                        scale = 1.0
                    else
                        scale = rget(String(temp[1:j]), j)
                    end
                elseif comm == 'O'
                    y = rget(String(temp[1:j]), j)
                    cbase += y
                elseif comm == 'F'
                    jfont = iget(String(temp[1:j]), j)
                elseif comm == 'M'
                    ialf = iget(String(temp[1:j]), j)
                end
                comm = c
                j = 0
                if c == '<' || c == '>' || c == '[' || c == ']'
                    imode = 0
                end
            end
        else
            # ---- Normal text mode ----
            ins = 0
            iskip = 0

            if i < ntxt
                if c == '#'
                    if txt[i+1] == '#'
                        i += 1
                    else
                        ins = 1
                    end
                end
                if c == '<'
                    if txt[i+1] == '<'
                        i += 1
                    else
                        icase = 1
                        jfont = ifont
                        iskip = 1
                    end
                end
                if c == '>'
                    if txt[i+1] == '>'
                        i += 1
                    else
                        icase = 0
                        jfont = ifont
                        iskip = 1
                    end
                end
                if c == '['
                    if txt[i+1] == '['
                        i += 1
                    else
                        icase = 1
                        jfont = 3
                        iskip = 1
                    end
                end
                if c == ']'
                    if txt[i+1] == ']'
                        i += 1
                    else
                        icase = 0
                        jfont = 3
                        iskip = 1
                    end
                end
            end

            if ins == 1
                # Flush accumulated characters before entering command mode
                imode = 1
                siz += ssum(String(temp[1:l]), l, jfont, ialf, scale)
                comm = c
                j = 0
            elseif iskip == 1
                iskip = 0
            else
                if ialf == 0
                    j += 1
                    if icase == 1
                        ii = Int(c)
                        if ii >= 96
                            ii -= 32
                        end
                        c = Char(ii)
                    end
                    temp[j] = c
                    if c != ' '
                        l = j
                    end
                else
                    j += 1
                    # Encode as octal escape: \ooo (4 chars)
                    num = @sprintf("%c%03o", '\\', Int(c) + 128)
                    for k in 1:4
                        temp[j+k-1] = num[k]
                    end
                    j += 3
                    l = j
                end
            end
        end

        i += 1
    end

    # Flush any remaining accumulated characters
    if j > 0
        siz += ssum(String(temp[1:l]), l, jfont, ialf, scale)
    end

    return siz * ht
end

# ---- 3D perspective transforms (used by text3!) ----

"""
    _trans3(gs, x, y, z) -> (u, v)

3D perspective transformation.
Matches graph.f90 `trans3` (lines 571-583).
"""
@inline function _trans3(gs::GraphState, x::Float64, y::Float64, z::Float64)
    xo = gs.ct * (gs.cp * x + gs.sp * y) + gs.st * z
    s = gs.rs / (gs.ro - xo)
    u = s * (-gs.sp * x + gs.cp * y) + gs.du
    v = s * (-gs.st * (gs.cp * x + gs.sp * y) + gs.ct * z) + gs.dv
    return (u, v)
end

"""
    _transw(gs, x, y) -> (u, v)

Window coordinate transformation.
Matches graph.f90 `transw` (lines 363-384).
"""
@inline function _transw(gs::GraphState, x::Float64, y::Float64)
    if gs.www == 0.0
        return (x, y)
    else
        ct = cos(2 * pi * gs.wwr / 360)
        st = sin(2 * pi * gs.wwr / 360)
        u = gs.xwll + x * ct - y * st
        v = gs.ywll + x * st + y * ct
        return (u, v)
    end
end

# ---- Render single character ----

"""
    dchr!(gs, io, x, y, c, ifont, ialf, ht, xx, yx, xy, yy)

Render a single character `c` at position `(x,y)` in PostScript,
with font `ifont`, size `ht`, and direction vectors `(xx,yx)`, `(xy,yy)`.
Handles landscape rotation, font selection, paren escaping, and
alternate character encoding (`ialf != 0`).
Matches graph.f90 `dchr` (lines 2795-2848).
"""
function dchr!(gs::GraphState, io::IO, x::Float64, y::Float64,
               c::Char, ifont::Int, ialf::Int, ht::Float64,
               xx::Float64, yx::Float64, xy::Float64, yy::Float64)
    if gs.land == 1
        u  = gs.uwidth - 72 * y + gs.ushift
        v  = 72 * x + gs.vshift
        ux = -72 * yx
        vx =  72 * xx
        uy = -72 * yy
        vy =  72 * xy
    else
        u  = 72 * x + gs.ushift
        v  = 72 * y + gs.vshift
        ux = 72 * xx
        vx = 72 * yx
        uy = 72 * xy
        vy = 72 * yy
    end

    # Clamp coordinates
    test = 2000.0
    if u > test; u = 2000.0; end
    test = -1000.0
    if u < test; u = -1000.0; end
    test = 2000.0
    if v > test; v = 2000.0; end
    test = -1000.0
    if v < test; v = -1000.0; end

    # moveto
    @printf(io, "%9.2f%9.2f moveto\n", u, v)

    # Font selection: strip the font name, then emit findfont/makefont
    fname = FONT_NAMES[ifont]
    sname, lname = stripv(fname)
    @printf(io, "/%s findfont [%7.2f%7.2f%7.2f%7.2f 0. 0.] makefont setfont\n",
            sname, ht * ux, ht * vx, ht * uy, ht * vy)

    # Character output
    if ialf == 0
        if c == '(' || c == ')'
            @printf(io, "(\\%c) show\n", c)
        else
            @printf(io, "(%c) show\n", c)
        end
    else
        @printf(io, "%c%03o) show\n", '\\', Int(c) + 128)
    end
end

# ---- Full 3D markup text engine ----

"""
    text3!(gs, io, txt, ntxt, ifont, ht, xo, yo, zo, xx, yx, zx, xy, yy, zy)

Draw a text string with DISSPLA markup in 3D perspective.
`ifont` and `ht` are the default font and height.
`(xo,yo,zo)` is the origin (lower-left of first character).
`(xx,yx,zx)` is the unit vector in the character advance direction.
`(xy,yy,zy)` is the unit vector in the character up direction.

Markup commands:
- `<` uppercase + restore ifont, `>` lowercase + restore ifont
- `[` uppercase + Symbol font, `]` lowercase + Symbol font
- Doubled `<<`, `>>`, `[[`, `]]` produce literal characters
- `#` enters command mode; letter = command, digits/./- = argument
- Commands: E(elevate), L(lower), H(height), O(offset), F(font), M(mode)

Matches graph.f90 `text3` (lines 1649-1836).
"""
function text3!(gs::GraphState, io::IO, txt::String, ntxt::Int,
                ifont::Int, ht::Float64,
                xo::Float64, yo::Float64, zo::Float64,
                xx::Float64, yx::Float64, zx::Float64,
                xy::Float64, yy::Float64, zy::Float64)
    scale = 1.0
    x = xo
    y = yo
    z = zo
    delta = 0.0
    cbase = 0.0
    celev = 0.0
    icase = 0
    ialf = 0
    imode = 0
    jfont = ifont
    i = 1
    j = 0
    temp = Vector{Char}(undef, 80)
    comm = ' '

    while i <= ntxt
        c = txt[i]

        if imode == 1
            # ---- Command argument accumulation ----
            if (c >= '0' && c <= '9') || c == '.' || c == '-'
                j += 1
                temp[j] = c
            elseif c == 'X'
                j += 1
                temp[j] = c
            else
                if j == 0
                    temp[1] = 'D'
                    j = 1
                end
                if comm == 'E'
                    if temp[1] == 'D'
                        temp[1] = '.'
                        temp[2] = '5'
                        j = 2
                    end
                    if temp[1] == 'X'
                        celev = 0.0
                    else
                        dy = rget(String(temp[1:j]), j)
                        if dy > 0
                            celev = ht * dy
                        else
                            celev = -ht * dy
                        end
                    end
                    delta = celev
                elseif comm == 'L'
                    if temp[1] == 'D'
                        temp[1] = '.'
                        temp[2] = '5'
                        j = 2
                    end
                    if temp[1] == 'X'
                        celev = 0.0
                    else
                        dy = rget(String(temp[1:j]), j)
                        if dy > 0
                            celev = -ht * dy
                        else
                            celev = ht * dy
                        end
                    end
                    delta = celev
                elseif comm == 'H'
                    if temp[1] == 'D'
                        temp[1] = '.'
                        temp[2] = '7'
                        j = 2
                    end
                    if temp[1] == 'X'
                        scale = 1.0
                    else
                        scale = rget(String(temp[1:j]), j)
                    end
                elseif comm == 'O'
                    dy = rget(String(temp[1:j]), j)
                    cbase += dy
                elseif comm == 'F'
                    jfont = iget(String(temp[1:j]), j)
                elseif comm == 'M'
                    ialf = iget(String(temp[1:j]), j)
                end
                comm = c
                j = 0
                if c == '<' || c == '>' || c == '[' || c == ']'
                    imode = 0
                    i -= 1
                end
            end
        else
            # ---- Normal text mode ----
            ins = 0
            iskip = 0

            if c == '#'
                if i < ntxt && txt[i+1] == '#'
                    i += 1
                else
                    ins = 1
                end
            end
            if c == '<'
                if i < ntxt && txt[i+1] == '<'
                    i += 1
                else
                    jfont = ifont
                    icase = 1
                    iskip = 1
                end
            end
            if c == '>'
                if i < ntxt && txt[i+1] == '>'
                    i += 1
                else
                    jfont = ifont
                    icase = 0
                    iskip = 1
                end
            end
            if c == '['
                if i < ntxt && txt[i+1] == '['
                    i += 1
                else
                    jfont = 3
                    icase = 1
                    iskip = 1
                end
            end
            if c == ']'
                if i < ntxt && txt[i+1] == ']'
                    i += 1
                else
                    jfont = 3
                    icase = 0
                    iskip = 1
                end
            end

            if ins == 1
                imode = 1
                comm = c
                j = 0
            elseif iskip == 1
                iskip = 0
            else
                if icase == 1
                    ii = Int(c)
                    if ii >= 96
                        ii -= 32
                    end
                    c = Char(ii)
                end
                siz = scale * ht
                wi = ssum(String([c]), 1, jfont, ialf, siz)
                xn = x + delta * xy
                yn = y + delta * yy
                zn = z + delta * zy
                wu, wv = _trans3(gs, xn, yn, zn)
                u, v = _transw(gs, wu, wv)
                wu, wv = _trans3(gs, xn + xx, yn + yx, zn + zx)
                ux, vx = _transw(gs, wu, wv)
                wu, wv = _trans3(gs, xn + xy, yn + yy, zn + zy)
                uy, vy = _transw(gs, wu, wv)
                ux -= u
                vx -= v
                uy -= u
                vy -= v
                dchr!(gs, io, u, v, c, jfont, ialf, siz, ux, vx, uy, vy)
                x += wi * xx
                y += wi * yx
                z += wi * zx
            end
        end

        i += 1
    end
end

# ---- 2D text wrapper ----

"""
    text2!(gs, io, txt, ntxt, ifont, ht, xo, yo, xx, yx, xy, yy)

Draw a 2D text string using the 3D text engine with z-components set to zero.
Matches graph.f90 `text2` (lines 1634-1647).
"""
function text2!(gs::GraphState, io::IO, txt::String, ntxt::Int,
                ifont::Int, ht::Float64,
                xo::Float64, yo::Float64,
                xx::Float64, yx::Float64,
                xy::Float64, yy::Float64)
    text3!(gs, io, txt, ntxt, ifont, ht,
           xo, yo, 0.0, xx, yx, 0.0, xy, yy, 0.0)
end
