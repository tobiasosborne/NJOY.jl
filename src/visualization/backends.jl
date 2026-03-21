# Visualization: Rendering backends (Proposal B)
# Consumes PlotSpec (from plotting.jl), renders to ASCII / PostScript / recipe Dict.
# ASCII = zero deps. PostScript = NJOY viewr compat. Recipe = CairoMakie/Plots.jl.

# --- Shared helpers ---
function _fmt_num(v::Float64)
    av = abs(v)
    av == 0.0 ? "0" : av >= 1e4 || av < 0.01 ? @sprintf("%.1e", v) :
    av >= 1.0 ? @sprintf("%.1f", v) : @sprintf("%.3f", v)
end

function _bresenham!(canvas::Matrix{Char}, r1::Int, c1::Int, r2::Int, c2::Int, ch::Char)
    dr, dc = abs(r2 - r1), abs(c2 - c1)
    sr, sc, err, r, c = r1<r2 ? 1 : -1, c1<c2 ? 1 : -1, dc - dr, r1, c1
    for _ in 1:(dr+dc+1)
        1<=r<=size(canvas,1) && 1<=c<=size(canvas,2) && canvas[r,c]==' ' && (canvas[r,c]=ch)
        r==r2 && c==c2 && break
        e2 = 2*err
        e2 > -dr && (err -= dr; c += sc)
        e2 < dc && (err += dc; r += sr)
    end
end

function _make_mapper(lo::Float64, hi::Float64, use_log::Bool)
    if use_log && lo > 0 && hi > lo
        ll, sp = log10(lo), log10(hi) - log10(lo)
        v -> (log10(clamp(v, lo, hi)) - ll) / sp
    else
        sp = hi - lo; sp == 0.0 && (sp = 1.0)
        v -> (clamp(v, lo, hi) - lo) / sp
    end
end

# --- ASCII backend (zero dependencies) ---

"""    render_ascii(p::PlotSpec; width=80, height=24) -> String
Render a PlotSpec as ASCII art for terminal display. Zero external dependencies."""
function render_ascii(p::PlotSpec; width::Int=80, height::Int=24)
    width >= 20 || throw(ArgumentError("width must be >= 20"))
    height >= 8 || throw(ArgumentError("height must be >= 8"))
    xlo, xhi, ylo, yhi = resolve_axes(p)
    xlo >= xhi && (xhi = xlo + 1.0); ylo >= yhi && (yhi = ylo + 1.0)
    lm = 12; pw = max(4, width - lm - 1); ph = max(4, height - 4)
    canvas = fill(' ', ph, pw)
    mks = ['.', '+', 'x', 'o', '*', '#', '@', '~']
    lx, ly = p.scale in (LOGLIN, LOGLOG), p.scale in (LINLOG, LOGLOG)
    xm, ym = _make_mapper(xlo, xhi, lx), _make_mapper(ylo, yhi, ly)
    for (ci, cv) in enumerate(p.curves)
        mk = mks[mod1(ci, length(mks))]
        for i in eachindex(cv.x)
            c = clamp(round(Int, xm(cv.x[i])*(pw-1))+1, 1, pw)
            r = clamp(ph - round(Int, ym(cv.y[i])*(ph-1)), 1, ph)
            canvas[r, c] = mk
        end
        if length(cv.x) >= 2 && cv.linestyle != LINE_INVISIBLE
            for i in 1:(length(cv.x)-1)
                c1 = clamp(round(Int, xm(cv.x[i])*(pw-1))+1, 1, pw)
                r1 = clamp(ph - round(Int, ym(cv.y[i])*(ph-1)), 1, ph)
                c2 = clamp(round(Int, xm(cv.x[i+1])*(pw-1))+1, 1, pw)
                r2 = clamp(ph - round(Int, ym(cv.y[i+1])*(ph-1)), 1, ph)
                _bresenham!(canvas, r1, c1, r2, c2, mk)
            end
        end
    end
    lines = String[]
    !isempty(p.title) && push!(lines, " "^max(0,div(width-length(p.title),2)) * p.title)
    push!(lines, " "^lm * "+" * "-"^pw * "+")
    ylr = Dict{Int,String}()
    for (row, f) in [(1, 1.0), (div(ph,2), 0.5), (ph, 0.0)]
        v = ly && ylo > 0 ? 10.0^(log10(ylo)+f*(log10(yhi)-log10(ylo))) : ylo+f*(yhi-ylo)
        ylr[row] = _fmt_num(v)
    end
    for r in 1:ph
        push!(lines, lpad(get(ylr, r, ""), lm-1) * " |" * String(canvas[r,:]) * "|")
    end
    push!(lines, " "^lm * "+" * "-"^pw * "+")
    tw = lm + 1 + pw + 1; xl = collect(' '^tw)
    lo_s, hi_s = _fmt_num(xlo), _fmt_num(xhi)
    mv = lx && xlo > 0 ? 10.0^((log10(xlo)+log10(xhi))/2) : (xlo+xhi)/2
    ms = _fmt_num(mv)
    for (s, p0) in [(lo_s, lm+2), (ms, lm+1+div(pw,2)-div(length(ms),2)+1),
                     (hi_s, lm+1+pw-length(hi_s)+1)]
        for (i, ch) in enumerate(s); q = p0+i-1; 1<=q<=tw && (xl[q]=ch); end
    end
    push!(lines, String(xl))
    !isempty(p.xaxis.label) && push!(lines, " "^max(0,div(width-length(p.xaxis.label),2)) * p.xaxis.label)
    if any(c -> !isempty(c.label), p.curves)
        push!(lines, "")
        for (ci, c) in enumerate(p.curves)
            !isempty(c.label) && push!(lines, "  $(mks[mod1(ci,length(mks))])  $(c.label)")
        end
    end
    join(lines, "\n")
end

# --- PostScript backend (NJOY viewr compatibility) ---

_ps_escape(s::AbstractString) = replace(replace(replace(s, "\\"=>"\\\\"), "("=>"\\("), ")"=>"\\)")

const _PS_COLORS = Dict{CurveColor,NTuple{3,Float64}}(
    COLOR_BLACK=>(0.,0.,0.), COLOR_RED=>(1.,0.,0.), COLOR_GREEN=>(0.,.6,0.),
    COLOR_BLUE=>(0.,0.,1.), COLOR_MAGENTA=>(1.,0.,1.), COLOR_CYAN=>(0.,.7,.7),
    COLOR_BROWN=>(.6,.3,0.), COLOR_PURPLE=>(.5,0.,.5), COLOR_ORANGE=>(1.,.5,0.))

const _PS_DASHES = Dict{LineStyle,String}(
    LINE_SOLID=>"[] 0", LINE_DASHED=>"[6 4] 0", LINE_CHAIN_DASH=>"[6 3 2 3] 0",
    LINE_CHAIN_DOT=>"[1 3] 0", LINE_DOT=>"[2 3] 0", LINE_INVISIBLE=>"[] 0")

"""Write a self-contained EPS file from a PlotSpec. NJOY viewr compatible."""
render_postscript(p::PlotSpec, fn::AbstractString) = open(io -> render_postscript(p, io), fn, "w")

"""Write PostScript rendering of a PlotSpec to an IO stream (EPS format)."""
function render_postscript(p::PlotSpec, io::IO)
    xlo, xhi, ylo, yhi = resolve_axes(p)
    xlo >= xhi && (xhi = xlo + 1.0); ylo >= yhi && (yhi = ylo + 1.0)
    ox, oy, pw, ph = 72.0, 72.0, 468.0, 360.0
    println(io, "%!PS-Adobe-3.0 EPSF-3.0\n%%BoundingBox: 54 54 558 468")
    println(io, "%%Title: ", p.title, "\n%%Creator: NJOY.jl visualization")
    println(io, "%%EndComments\n/Helvetica findfont 12 scalefont setfont")
    if !isempty(p.title)
        @printf(io, "%.1f %.1f moveto (%s) dup stringwidth pop 2 div neg 0 rmoveto show\n",
                ox+pw/2, oy+ph+24, _ps_escape(p.title))
    end
    if !isempty(p.subtitle)
        println(io, "/Helvetica findfont 10 scalefont setfont")
        @printf(io, "%.1f %.1f moveto (%s) dup stringwidth pop 2 div neg 0 rmoveto show\n",
                ox+pw/2, oy+ph+10, _ps_escape(p.subtitle))
        println(io, "/Helvetica findfont 12 scalefont setfont")
    end
    println(io, "0.5 setlinewidth")
    @printf(io, "%.1f %.1f moveto %.1f %.1f lineto %.1f %.1f lineto %.1f %.1f lineto closepath stroke\n",
            ox, oy, ox+pw, oy, ox+pw, oy+ph, ox, oy+ph)
    @printf(io, "%.1f %.1f moveto (%s) dup stringwidth pop 2 div neg 0 rmoveto show\n",
            ox+pw/2, oy-30, _ps_escape(p.xaxis.label))
    @printf(io, "gsave\n%.1f %.1f translate 90 rotate 0 0 moveto (%s) dup stringwidth pop 2 div neg 0 rmoveto show\ngrestore\n",
            ox-40, oy+ph/2, _ps_escape(p.yaxis.label))
    lx, ly = p.scale in (LOGLIN, LOGLOG), p.scale in (LINLOG, LOGLOG)
    xm, ym = _make_mapper(xlo, xhi, lx), _make_mapper(ylo, yhi, ly)
    println(io, "gsave")
    @printf(io, "%.1f %.1f moveto %.1f %.1f lineto %.1f %.1f lineto %.1f %.1f lineto closepath clip newpath\n",
            ox, oy, ox+pw, oy, ox+pw, oy+ph, ox, oy+ph)
    for curve in p.curves
        isempty(curve.x) && continue; curve.linestyle == LINE_INVISIBLE && continue
        r, g, b = _PS_COLORS[curve.color]
        @printf(io, "gsave\n%.2f %.2f %.2f setrgbcolor\n%.1f setlinewidth\n%s setdash\nnewpath\n",
                r, g, b, max(0.5, curve.thickness*0.5), _PS_DASHES[curve.linestyle])
        for i in eachindex(curve.x)
            @printf(io, "%.2f %.2f %s\n", ox+pw*xm(curve.x[i]), oy+ph*ym(curve.y[i]),
                    i==1 ? "moveto" : "lineto")
        end
        println(io, "stroke\ngrestore")
    end
    println(io, "grestore\nshowpage\n%%EOF")
end

# --- Recipe backend (CairoMakie / Plots.jl interop) ---

const _RECIPE_LS = Dict{LineStyle,Symbol}(LINE_SOLID=>:solid, LINE_DASHED=>:dash,
    LINE_CHAIN_DASH=>:dashdot, LINE_CHAIN_DOT=>:dashdot, LINE_DOT=>:dot, LINE_INVISIBLE=>:solid)
const _RECIPE_COLORS = Dict{CurveColor,NTuple{3,Float64}}(
    COLOR_BLACK=>(0.,0.,0.), COLOR_RED=>(.8,0.,0.), COLOR_GREEN=>(0.,.6,0.),
    COLOR_BLUE=>(0.,0.,.8), COLOR_MAGENTA=>(.8,0.,.8), COLOR_CYAN=>(0.,.6,.6),
    COLOR_BROWN=>(.6,.3,0.), COLOR_PURPLE=>(.5,0.,.5), COLOR_ORANGE=>(1.,.5,0.))

"""Convert PlotSpec to a Dict for Plots.jl / CairoMakie.
Keys: title, subtitle, xlabel, ylabel, xscale, yscale, xlim, ylim, series, etc."""
function to_plot_recipe(p::PlotSpec)
    lx, ly = p.scale in (LOGLIN, LOGLOG), p.scale in (LINLOG, LOGLOG)
    series = [Dict{String,Any}("x"=>copy(c.x), "y"=>copy(c.y), "label"=>c.label,
        "linestyle"=>_RECIPE_LS[c.linestyle], "color"=>_RECIPE_COLORS[c.color],
        "linewidth"=>max(1, c.thickness)) for c in p.curves]
    ax, ay = p.xaxis, p.yaxis
    Dict{String,Any}("title"=>p.title, "subtitle"=>p.subtitle,
        "xlabel"=>ax.label, "ylabel"=>ay.label,
        "xscale"=>lx ? :log10 : :identity, "yscale"=>ly ? :log10 : :identity,
        "xlim"=>needs_autoscale(ax) ? :auto : (ax.lo, ax.hi),
        "ylim"=>needs_autoscale(ay) ? :auto : (ay.lo, ay.hi),
        "grid"=>p.grid, "orientation"=>p.orientation, "series"=>series)
end

"""Convert HeatmapSpec to a recipe Dict for heatmap rendering."""
function to_plot_recipe(hm::HeatmapSpec)
    Dict{String,Any}("type"=>:heatmap, "title"=>hm.title, "xlabel"=>hm.xlabel,
        "ylabel"=>hm.ylabel, "xbins"=>copy(hm.xbins), "ybins"=>copy(hm.ybins),
        "values"=>copy(hm.values), "clabel"=>hm.clabel)
end
