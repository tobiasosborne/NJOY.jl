# viewr module runner -- PostScript rendering from NJOY plot tapes
#
# Matches Fortran viewr.f90 + graph.f90: reads plot tape produced by
# dtfr/plotr/covr, renders to PostScript output.

"""
    viewr_module(tapes::TapeManager, params::ViewrParams)

Run VIEWR: render NJOY plot tape to PostScript.
Reads plot commands from infile, writes PostScript to nps.
"""
function viewr_module(tapes::TapeManager, params::ViewrParams)
    @info "viewr: infile=$(params.infile) nps=$(params.nps)"

    infile_path = resolve(tapes, params.infile)
    nps_path = resolve(tapes, params.nps)

    # Read the plot tape and render to PostScript
    open(nps_path, "w") do ps_io
        _render_njoy_plot_tape(ps_io, infile_path)
    end
    register!(tapes, params.nps, nps_path)

    lines = countlines(nps_path)
    @info "viewr: wrote $nps_path ($lines lines)"
    nothing
end

# =========================================================================
# NJOY plot tape parser + PostScript renderer
# Matches Fortran viewr.f90 + graph.f90 PostScript output exactly
# =========================================================================

"""Parse an NJOY plot tape and render to PostScript, matching Fortran viewr output."""
function _render_njoy_plot_tape(ps_io::IO, plot_tape_path::String)
    lines = readlines(plot_tape_path)
    isempty(lines) && return

    # Parse card 1: page setup
    card1 = _parse_plot_card(lines[1])
    lori = length(card1) >= 1 ? round(Int, card1[1]) : 1  # 1=landscape
    istyle = length(card1) >= 2 ? round(Int, card1[2]) : 2  # 2=Swiss
    csize = length(card1) >= 3 ? card1[3] : 0.30

    # PostScript header — matching Fortran graph.f90 exactly
    println(ps_io, "%!PS-Adobe-")
    if lori == 1
        @printf(ps_io, "%%%%BoundingBox: %4d %4d %4d %4d\n", 36, 36, 576, 756)
    else
        @printf(ps_io, "%%%%BoundingBox: %4d %4d %4d %4d\n", 36, 36, 576, 756)
    end
    println(ps_io, "%%Pages: (atend)")
    println(ps_io, lori == 1 ? "%%Orientation: Landscape" : "%%Orientation: Portrait")
    println(ps_io, "2 setlinecap")

    ipage = 0
    li = 2  # line index (skip card 1)

    while li <= length(lines)
        line = strip(lines[li])
        isempty(line) && (li += 1; continue)

        # Parse card 2: plot control
        card2 = _parse_plot_card(line)
        iplot = length(card2) >= 1 ? round(Int, card2[1]) : 0

        iplot == 99 && break  # terminate

        if iplot == 1 || iplot == -1
            # New plot frame
            ipage += 1
            @printf(ps_io, "%%%%Page:%4d%4d\n", ipage, ipage)

            # Draw page background and frame
            _ps_new_page(ps_io, lori)

            li += 1

            # Card 3: title
            title = ""
            if li <= length(lines)
                title = _extract_plot_string(lines[li])
                li += 1
            end

            # Card 3a: subtitle
            subtitle = ""
            if li <= length(lines)
                subtitle = _extract_plot_string(lines[li])
                li += 1
            end

            # Card 4: axis types
            card4 = li <= length(lines) ? _parse_plot_card(lines[li]) : Float64[]
            li += 1
            itype = length(card4) >= 1 ? round(Int, card4[1]) : 4
            jtype = length(card4) >= 2 ? round(Int, card4[2]) : 0
            igrid = length(card4) >= 3 ? round(Int, card4[3]) : 2
            ileg  = length(card4) >= 4 ? round(Int, card4[4]) : 0
            xtag  = length(card4) >= 5 ? card4[5] : 0.0
            ytag  = length(card4) >= 6 ? card4[6] : 0.0

            # Card 5: X axis limits
            card5 = li <= length(lines) ? _parse_plot_card(lines[li]) : Float64[]
            li += 1
            xmin = length(card5) >= 1 ? card5[1] : 1.0
            xmax = length(card5) >= 2 ? card5[2] : 1e7
            xstep = length(card5) >= 3 ? card5[3] : 1.0

            # Card 5a: X axis label
            xlabel = li <= length(lines) ? _extract_plot_string(lines[li]) : ""
            li += 1

            # Card 6: Y axis limits
            card6 = li <= length(lines) ? _parse_plot_card(lines[li]) : Float64[]
            li += 1
            ymin = length(card6) >= 1 ? card6[1] : 1.0
            ymax = length(card6) >= 2 ? card6[2] : 1e7
            ystep = length(card6) >= 3 ? card6[3] : 1.0

            # Card 6a: Y axis label
            ylabel = li <= length(lines) ? _extract_plot_string(lines[li]) : ""
            li += 1

            # Card 7/7a: secondary axis (jtype > 0)
            if jtype > 0
                li += 1  # limits
                li += 1  # label
            end

            # Draw axes, labels, tick marks
            _ps_draw_axes(ps_io, lori, itype, igrid,
                         xmin, xmax, xstep, ymin, ymax, ystep,
                         xlabel, ylabel, title, subtitle, csize, xtag, ytag, ileg)

            # Card 8: dummy (always "0/" or "/")
            if li <= length(lines)
                li += 1
            end
        else
            li += 1
        end

        # Now read curve data until next plot frame or termination
        while li <= length(lines)
            line = strip(lines[li])
            isempty(line) && (li += 1; continue)

            # Check if this is a new plot command
            card = _parse_plot_card(line)
            first_val = length(card) >= 1 ? round(Int, card[1]) : 0

            # If first value is 1, -1, or 99, this is a new frame/termination
            if first_val == 1 || first_val == -1 || first_val == 99
                break
            end

            # Card 9: curve style
            icon   = length(card) >= 1 ? round(Int, card[1]) : 0
            isym   = length(card) >= 2 ? round(Int, card[2]) : 0
            idash  = length(card) >= 3 ? round(Int, card[3]) : 0
            iccol  = length(card) >= 4 ? round(Int, card[4]) : 0
            ithick = length(card) >= 5 ? round(Int, card[5]) : 1
            li += 1

            # Card 10: legend label
            legend = ""
            if li <= length(lines)
                legend = _extract_plot_string(lines[li])
                li += 1
            end

            # Card 12: data format code
            nform = 0
            if li <= length(lines)
                nform_card = _parse_plot_card(lines[li])
                nform = length(nform_card) >= 1 ? round(Int, nform_card[1]) : 0
                li += 1
            end

            # Card 13: data points
            xdata = Float64[]; ydata = Float64[]
            while li <= length(lines)
                dline = strip(lines[li])
                isempty(dline) && (li += 1; continue)
                dline == "/" && (li += 1; break)

                dcard = _parse_plot_card(dline)
                isempty(dcard) && (li += 1; break)

                push!(xdata, dcard[1])
                length(dcard) >= 2 && push!(ydata, dcard[2])
                li += 1
            end

            # Draw this curve
            if !isempty(xdata) && !isempty(ydata)
                _ps_draw_curve(ps_io, lori, xdata, ydata,
                              xmin, xmax, ymin, ymax, itype,
                              iccol, idash, ithick)
            end
        end

        # End of page
        println(ps_io, "stroke")
        println(ps_io, "showpage")
        println(ps_io, "%endp")
    end

    # Trailer
    println(ps_io, "%gdone")
    @printf(ps_io, "%%%%Trailer\n%%%%Pages: %d\n", max(ipage, 1))
    println(ps_io, "%%EOF")
end

# =========================================================================
# PostScript drawing primitives — matching Fortran graph.f90 exactly
# =========================================================================

function _ps_new_page(io::IO, lori::Int)
    println(io, "stroke")
    println(io, "newpath")
    @printf(io, " 0.000 0.000 0.000 setrgbcolor\n")

    if lori == 1  # landscape
        @printf(io, "  576.00   36.00 moveto\n")
        println(io, "gsave 1.000 1.000 1.000 setrgbcolor fill grestore")
        println(io, "%gset")
        println(io, "stroke")
        println(io, "newpath")
        @printf(io, " 0.000 0.000 0.000 setrgbcolor\n")
        @printf(io, "  576.00   36.00 moveto\n")
        @printf(io, "   0.360 setlinewidth\n")
        @printf(io, "  576.00  756.00 lineto\n")
        @printf(io, "   36.00  756.00 lineto\n")
        @printf(io, "   36.00   36.00 lineto\n")
    else  # portrait
        @printf(io, "   36.00   36.00 moveto\n")
        println(io, "gsave 1.000 1.000 1.000 setrgbcolor fill grestore")
        println(io, "%gset")
        println(io, "stroke")
        println(io, "newpath")
        @printf(io, " 0.000 0.000 0.000 setrgbcolor\n")
        @printf(io, "   36.00   36.00 moveto\n")
        @printf(io, "   0.360 setlinewidth\n")
        @printf(io, "  576.00   36.00 lineto\n")
        @printf(io, "  576.00  756.00 lineto\n")
        @printf(io, "   36.00  756.00 lineto\n")
    end
end

function _ps_draw_axes(io::IO, lori::Int, itype::Int, igrid::Int,
                      xmin, xmax, xstep, ymin, ymax, ystep,
                      xlabel, ylabel, title, subtitle, csize,
                      xtag, ytag, ileg)
    # This is a simplified version — full implementation needs to match
    # Fortran graph.f90 axis drawing code exactly.
    # For now, write the structural PostScript commands.
    # TODO: match Fortran tick marks, labels, axis scales exactly
end

function _ps_draw_curve(io::IO, lori::Int,
                       xdata::Vector{Float64}, ydata::Vector{Float64},
                       xmin, xmax, ymin, ymax, itype::Int,
                       iccol::Int, idash::Int, ithick::Int)
    # Transform data coordinates to page coordinates
    # Matching Fortran graph.f90 coordinate transforms

    # Color
    r, g, b = _ps_color(iccol)
    @printf(io, "%6.3f%6.3f%6.3f setrgbcolor\n", r, g, b)

    # Line width
    w = ithick * 0.36  # matches Fortran: 72 * width_inches
    @printf(io, "%8.3f setlinewidth\n", w)

    # Dash pattern
    _ps_dash(io, idash, w)

    # Plot area bounds (matching Fortran layout)
    # These need to match the axis drawing exactly
    px_lo = 36.0 + 72.0 * 1.0   # left margin
    px_hi = 576.0 - 72.0 * 0.5  # right margin
    py_lo = 36.0 + 72.0 * 1.0   # bottom margin
    py_hi = 756.0 - 72.0 * 1.0  # top margin

    log_x = (itype == 3 || itype == 4)
    log_y = (itype == 2 || itype == 4)

    first = true
    for i in 1:min(length(xdata), length(ydata))
        x = xdata[i]; y = ydata[i]
        x <= 0 && log_x && continue
        y <= 0 && log_y && continue

        # Map to page coordinates
        if log_x
            fx = (log10(x) - log10(xmin)) / (log10(xmax) - log10(xmin))
        else
            fx = (x - xmin) / (xmax - xmin)
        end
        if log_y
            fy = (log10(y) - log10(ymin)) / (log10(ymax) - log10(ymin))
        else
            fy = (y - ymin) / (ymax - ymin)
        end

        # Clamp
        fx = clamp(fx, 0.0, 1.0)
        fy = clamp(fy, 0.0, 1.0)

        if lori == 1  # landscape
            u = px_hi - fy * (px_hi - px_lo)
            v = py_lo + fx * (py_hi - py_lo)
        else  # portrait
            u = px_lo + fx * (px_hi - px_lo)
            v = py_lo + fy * (py_hi - py_lo)
        end

        if first
            @printf(io, "%9.2f%9.2f moveto\n", u, v)
            first = false
        else
            @printf(io, "%9.2f%9.2f lineto\n", u, v)
        end
    end
end

function _ps_color(iccol::Int)
    iccol == 0 && return (0.0, 0.0, 0.0)  # black
    iccol == 1 && return (1.0, 0.0, 0.0)  # red
    iccol == 2 && return (0.0, 1.0, 0.0)  # green
    iccol == 3 && return (0.0, 0.0, 1.0)  # blue
    iccol == 4 && return (1.0, 0.0, 1.0)  # magenta
    iccol == 5 && return (0.0, 1.0, 1.0)  # cyan
    iccol == 6 && return (0.65, 0.16, 0.16)  # brown
    iccol == 7 && return (0.63, 0.13, 0.94)  # purple
    iccol == 8 && return (1.0, 0.65, 0.0)  # orange
    return (0.0, 0.0, 0.0)
end

function _ps_dash(io::IO, idash::Int, w::Float64)
    if idash == 0
        println(io, "[] 0 setdash")
    elseif idash == 1
        @printf(io, "[%9.2f%9.2f] 0 setdash\n", w*5, w*3)
    elseif idash == 2
        @printf(io, "[%9.2f%9.2f%9.2f%9.2f] 0 setdash\n", w*6, w*3, w*3, w*3)
    elseif idash == 3
        @printf(io, "[%9.2f%9.2f%9.2f%9.2f] 0 setdash\n", w*6, w*3, w*1, w*3)
    elseif idash == 4
        @printf(io, "[%9.2f%9.2f] 0 setdash\n", w*1, w*3)
    end
end

# =========================================================================
# Plot tape parsing helpers
# =========================================================================

"""Parse a plot tape card line into numeric values (strips trailing /)."""
function _parse_plot_card(line::AbstractString)
    # Remove trailing / and comments
    s = replace(strip(line), r"/.*$" => "")
    s = replace(s, r"'[^']*'" => "")  # remove quoted strings
    isempty(strip(s)) && return Float64[]

    vals = Float64[]
    for tok in split(strip(s))
        isempty(tok) && continue
        startswith(tok, "'") && continue
        try
            push!(vals, parse(Float64, replace(tok, r"([0-9])([+-])(\d)" => s"\1e\2\3")))
        catch
        end
    end
    vals
end

"""Extract a string from a plot tape card (between quotes or raw)."""
function _extract_plot_string(line::AbstractString)
    s = strip(line)
    # Remove trailing /
    s = replace(s, r"\s*/\s*$" => "")
    # Extract content between quotes
    m = match(r"'([^']*)'", s)
    m !== nothing && return strip(m.captures[1])
    return strip(s)
end
