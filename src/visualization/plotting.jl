# Visualization: Plot specification types and data extraction (Proposal B)
#
# Plot specifications are PURE DATA (structs). No rendering happens here.
# This mirrors the NJOY plotr/viewr separation: plotr emits plot commands,
# viewr consumes them. Rendering is delegated to backends.jl.

# --- Enums matching NJOY plotr/viewr conventions ---

"""Axis scaling mode (NJOY itype: 1=linlin, 2=linlog, 3=loglin, 4=loglog)."""
@enum AxisScale LINLIN=1 LINLOG=2 LOGLIN=3 LOGLOG=4

"""Grid/tick mark control (NJOY igrid)."""
@enum GridStyle GRID_NONE=0 GRID_LINES=1 GRID_TICKS_OUT=2 GRID_TICKS_IN=3

"""Line dash pattern (NJOY idash)."""
@enum LineStyle LINE_SOLID=0 LINE_DASHED=1 LINE_CHAIN_DASH=2 LINE_CHAIN_DOT=3 LINE_DOT=4 LINE_INVISIBLE=5

"""Plot marker shapes (NJOY isym subset)."""
@enum MarkerShape MARKER_NONE=-1 MARKER_SQUARE=0 MARKER_OCTAGON=1 MARKER_TRIANGLE=2 MARKER_CROSS=3 MARKER_EX=4 MARKER_DIAMOND=5 MARKER_CIRCLE=16

"""Curve colors (NJOY iccol)."""
@enum CurveColor COLOR_BLACK=0 COLOR_RED=1 COLOR_GREEN=2 COLOR_BLUE=3 COLOR_MAGENTA=4 COLOR_CYAN=5 COLOR_BROWN=6 COLOR_PURPLE=7 COLOR_ORANGE=8

# --- Plot specification types (pure data, no rendering logic) ---

"""A single curve with x,y data and styling. All fields are plain data."""
struct CurveData
    x::Vector{Float64}
    y::Vector{Float64}
    label::String
    color::CurveColor
    linestyle::LineStyle
    marker::MarkerShape
    marker_interval::Int  # 0=no markers, n=marker every nth point (NJOY icon)
    thickness::Int        # line thickness (NJOY ithick, 0=invisible)
    function CurveData(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                       label::AbstractString="", color::CurveColor=COLOR_BLACK,
                       linestyle::LineStyle=LINE_SOLID, marker::MarkerShape=MARKER_NONE,
                       marker_interval::Integer=0, thickness::Integer=1)
        length(x) == length(y) || throw(ArgumentError(
            "x and y must have equal length (got $(length(x)) vs $(length(y)))"))
        new(collect(Float64, x), collect(Float64, y), String(label),
            color, linestyle, marker, Int(marker_interval), Int(thickness))
    end
end

"""Axis limits and labeling. NaN limits trigger auto-scaling (NJOY default)."""
struct AxisSpec
    lo::Float64; hi::Float64; step::Float64; label::String
    function AxisSpec(; lo::Real=NaN, hi::Real=NaN, step::Real=0.0, label::AbstractString="")
        new(Float64(lo), Float64(hi), Float64(step), String(label))
    end
end

needs_autoscale(a::AxisSpec) = isnan(a.lo) || isnan(a.hi)

"""Complete specification of a 2-D plot. Purely descriptive, no rendering state.
Mirrors NJOY plotr output: titles, axis config, curve data."""
struct PlotSpec
    title::String; subtitle::String
    xaxis::AxisSpec; yaxis::AxisSpec
    curves::Vector{CurveData}
    scale::AxisScale; grid::GridStyle; orientation::Symbol
    function PlotSpec(; title::AbstractString="", subtitle::AbstractString="",
            xaxis::AxisSpec=AxisSpec(label="Energy (eV)"),
            yaxis::AxisSpec=AxisSpec(label="Cross section (barns)"),
            curves::AbstractVector{CurveData}=CurveData[],
            scale::AxisScale=LOGLOG, grid::GridStyle=GRID_TICKS_OUT,
            orientation::Symbol=:landscape)
        orientation in (:landscape, :portrait) || throw(ArgumentError(
            "orientation must be :landscape or :portrait, got :$orientation"))
        new(String(title), String(subtitle), xaxis, yaxis,
            collect(CurveData, curves), scale, grid, orientation)
    end
end

"""Return a new PlotSpec with curve `c` appended (functional style)."""
function add_curve(p::PlotSpec, c::CurveData)
    PlotSpec(title=p.title, subtitle=p.subtitle, xaxis=p.xaxis, yaxis=p.yaxis,
             curves=vcat(p.curves, [c]), scale=p.scale, grid=p.grid,
             orientation=p.orientation)
end

"""2-D heatmap specification for covariance/correlation matrices."""
struct HeatmapSpec
    title::String; subtitle::String; xlabel::String; ylabel::String
    xbins::Vector{Float64}; ybins::Vector{Float64}
    values::Matrix{Float64}; clabel::String
    function HeatmapSpec(; title::AbstractString="", subtitle::AbstractString="",
            xlabel::AbstractString="Energy group", ylabel::AbstractString="Energy group",
            xbins::AbstractVector{<:Real}=Float64[], ybins::AbstractVector{<:Real}=Float64[],
            values::AbstractMatrix{<:Real}=Matrix{Float64}(undef,0,0), clabel::AbstractString="")
        new(String(title), String(subtitle), String(xlabel), String(ylabel),
            collect(Float64, xbins), collect(Float64, ybins),
            collect(Float64, values), String(clabel))
    end
end

# --- Auto-scaling (translation of NJOY viewr ascalv) ---

"""Auto-scale linear axis. Returns (lo, hi, step, major). From NJOY ascalv."""
function auto_scale_linear(data_min::Real, data_max::Real; min_divisions::Int=4)
    zmin, zmax = Float64(data_min), Float64(data_max)
    if zmax <= zmin || min_divisions <= 0
        return (0.0, 2.0 * max(abs(zmax), 1.0), 1.0, 1)
    end
    if zmin != 0.0 && zmax != 0.0
        ratio = zmax / zmin
        abs(ratio) / 1000.0 >= 1.0 && (zmin = 0.0)
        abs(ratio) * 1000.0 <= 1.0 && (zmax = 0.0)
    end
    zmax -= abs(zmax) / 1e6; zmin += abs(zmin) / 1e6
    p = (zmax - zmin) / min_divisions
    tenk = 1.0
    if p < 1.0
        p = 1.0 / p
        while p >= 10.0; p /= 10.0; tenk *= 10.0; end
        p = 10.0 / p; tenk = 1.0 / (10.0 * tenk)
    else
        while p >= 10.0; p /= 10.0; tenk *= 10.0; end
    end
    p_nice = p <= 1.99 ? 1.0 : (p < 4.99 ? 2.0 : 5.0)
    dz = p_nice * tenk
    n1 = floor(Int, data_min / dz)
    n1 * dz > data_min && (n1 -= 1)
    n2 = ceil(Int, data_max / dz)
    n2 * dz < data_max && (n2 += 1)
    return (n1 * dz, n2 * dz, dz, n2 - n1)
end

"""Auto-scale log axis. Returns (lo, hi) where lo <= data_min, hi >= data_max.
Snaps to whole decades. Uses Base.Math.exp10 for exact powers of 10."""
function auto_scale_log(data_min::Real, data_max::Real; max_decades::Int=13)
    dmin, dmax = Float64(data_min), Float64(data_max)
    dmin > 0 || throw(ArgumentError("log scale requires positive data_min, got $dmin"))
    dmax > dmin || throw(ArgumentError("data_max must exceed data_min"))
    _pow10(n::Int) = Base.Math.exp10(n)  # exact for small integer exponents
    top_exp = log10(dmax)
    hi_exp = top_exp >= 0 ? floor(Int, top_exp + 0.99) : floor(Int, top_exp)
    hi = _pow10(hi_exp)
    bot_exp = floor(Int, log10(dmin))
    lo = _pow10(bot_exp)
    lo > dmin && (lo = _pow10(bot_exp - 1))  # guard against float rounding
    hi / lo > _pow10(max_decades) && (lo = _pow10(floor(Int, top_exp - max_decades)))
    return (lo, hi)
end

"""Compute effective axis limits, auto-scaling from curve data where needed."""
function resolve_axes(p::PlotSpec)
    isempty(p.curves) && return (0.0, 1.0, 0.0, 1.0)
    xlo, xhi, ylo, yhi = Inf, -Inf, Inf, -Inf
    for c in p.curves
        for v in c.x; isfinite(v) || continue; v < xlo && (xlo = v); v > xhi && (xhi = v); end
        for v in c.y; isfinite(v) || continue; v < ylo && (ylo = v); v > yhi && (yhi = v); end
    end
    log_x, log_y = p.scale in (LOGLIN, LOGLOG), p.scale in (LINLOG, LOGLOG)
    _auto_one(ax, lo, hi, curves, use_log, field) = begin
        if needs_autoscale(ax)
            if use_log
                pos = minimum((v for c in curves for v in getfield(c, field) if v > 0); init=Inf)
                pos == Inf && (pos = 1e-5)
                return auto_scale_log(pos, max(hi, pos * 10))
            else
                r = auto_scale_linear(lo, hi); return (r[1], r[2])
            end
        else
            return (ax.lo, ax.hi)
        end
    end
    xlo_f, xhi_f = _auto_one(p.xaxis, xlo, xhi, p.curves, log_x, :x)
    ylo_f, yhi_f = _auto_one(p.yaxis, ylo, yhi, p.curves, log_y, :y)
    return (xlo_f, xhi_f, ylo_f, yhi_f)
end

# --- MT name lookup ---
const _MT_NAMES = Dict{Int,String}(
    1=>"Total", 2=>"Elastic", 3=>"Nonelastic", 4=>"Inelastic",
    16=>"(n,2n)", 17=>"(n,3n)", 18=>"Fission", 27=>"Absorption",
    51=>"Inelastic L1", 101=>"Disappearance", 102=>"Capture (n,g)",
    103=>"(n,p)", 104=>"(n,d)", 105=>"(n,t)", 106=>"(n,He3)", 107=>"(n,a)")
_mt_label(mt::Int) = get(_MT_NAMES, mt, "MT$mt")

# --- Data extraction ---

"""Extract a single MT cross section from a PointwiseMaterial as CurveData."""
function extract_curve(pm, mt::Integer; label::AbstractString="",
                       color::CurveColor=COLOR_BLACK, kwargs...)
    col_idx = findfirst(==(Int(mt)), pm.mt_list)
    col_idx === nothing && throw(ArgumentError(
        "MT=$mt not found in material (available: $(pm.mt_list))"))
    xs = pm.cross_sections[:, col_idx]; energies = pm.energies
    mask = xs .> 0
    lbl = isempty(label) ? _mt_label(Int(mt)) : label
    CurveData(energies[mask], xs[mask]; label=lbl, color=color, kwargs...)
end

"""Create a PlotSpec from a PointwiseMaterial. Auto-assigns colors per curve."""
function plot_material(pm; mts::AbstractVector{<:Integer}=Int[],
                       title::AbstractString="", kwargs...)
    mt_plot = isempty(mts) ? pm.mt_list : collect(Int, mts)
    colors = [COLOR_BLACK, COLOR_RED, COLOR_BLUE, COLOR_GREEN,
              COLOR_MAGENTA, COLOR_CYAN, COLOR_BROWN, COLOR_PURPLE, COLOR_ORANGE]
    curves = CurveData[]
    for (i, mt) in enumerate(mt_plot)
        mt in pm.mt_list || continue
        push!(curves, extract_curve(pm, mt; color=colors[mod1(i, length(colors))]))
    end
    auto_title = isempty(title) ? "MAT $(pm.mat)" : title
    PlotSpec(title=auto_title, xaxis=AxisSpec(label="Energy (eV)"),
             yaxis=AxisSpec(label="Cross section (barns)"),
             curves=curves, scale=LOGLOG; kwargs...)
end

"""Create a step-function PlotSpec from multigroup cross sections."""
function plot_multigroup(mg; mt::Int=1, title::AbstractString="")
    col = findfirst(==(mt), mg.mt_list)
    col === nothing && throw(ArgumentError("MT $mt not found in MultiGroupXS"))
    ng = length(mg.group_bounds) - 1
    xvals, yvals = Float64[], Float64[]
    for g in 1:ng
        push!(xvals, mg.group_bounds[g]); push!(yvals, mg.xs[g, col])
        push!(xvals, mg.group_bounds[g+1]); push!(yvals, mg.xs[g, col])
    end
    auto_title = isempty(title) ? "Multigroup $(_mt_label(mt))" : title
    PlotSpec(title=auto_title, xaxis=AxisSpec(label="Energy (eV)"),
             yaxis=AxisSpec(label="Group-averaged XS (barns)"),
             curves=[CurveData(xvals, yvals; label=_mt_label(mt))], scale=LOGLOG)
end

"""Create a HeatmapSpec from a covariance matrix."""
function plot_covariance(cov; title::AbstractString="")
    auto_title = isempty(title) ?
        "$(cov.is_relative ? "Relative " : "")Covariance MT$(cov.mt1)/MT$(cov.mt2)" : title
    HeatmapSpec(title=auto_title, xbins=cov.energy_groups, ybins=cov.energy_groups,
                values=cov.matrix,
                clabel=cov.is_relative ? "Relative covariance" : "Absolute covariance")
end

"""Plot probability table at energy index `ie`."""
function plot_probability_table(ptable, ie::Int; title::AbstractString="")
    1 <= ie <= length(ptable.energies) || throw(BoundsError(ptable.energies, ie))
    E = ptable.energies[ie]
    auto_title = isempty(title) ? @sprintf("Probability table at %.3g eV", E) : title
    nb = ptable.nbins
    cprob = cumsum(ptable.prob[1:nb, ie])
    cprob[end] > 0 && (cprob ./= cprob[end])
    colors = [COLOR_BLACK, COLOR_RED, COLOR_BLUE, COLOR_GREEN]
    styles = [LINE_SOLID, LINE_DASHED, LINE_DOT, LINE_CHAIN_DASH]
    curves = CurveData[]
    for (i, (name, mat)) in enumerate([("Total", ptable.total), ("Elastic", ptable.elastic),
                                        ("Fission", ptable.fission), ("Capture", ptable.capture)])
        vals = mat[1:nb, ie]
        any(v -> v > 0, vals) && push!(curves, CurveData(cprob, vals;
            label=name, color=colors[i], linestyle=styles[i]))
    end
    PlotSpec(title=auto_title, xaxis=AxisSpec(label="Cumulative probability"),
             yaxis=AxisSpec(label="Cross section (barns)"), curves=curves, scale=LINLOG)
end
