# ENDF interpolation and integration -- translation of NJOY2016 endf.f90:
#   terp1, terpa (-> interpolate), intega (-> integrate), gral (-> panel_integral)

"""
    terp1(x1, y1, x2, y2, x, law; coulomb_threshold=0.0) -> Float64

Interpolate one point between `(x1,y1)` and `(x2,y2)` using ENDF interpolation `law`.
The `coulomb_threshold` keyword is the kinematic threshold for Coulomb penetrability
interpolation (law 6). Defaults to 0.0.
"""
function terp1(x1::Real, y1::Real, x2::Real, y2::Real,
               x::Real, law::InterpolationLaw;
               coulomb_threshold::Real=0.0)
    x2 == x1 && return Float64(y1)
    (law == Histogram || y2 == y1 || x == x1) && return Float64(y1)
    if law == LinLin
        y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    elseif law == LinLog
        y1 + log(x / x1) * (y2 - y1) / log(x2 / x1)
    elseif law == LogLin
        y1 * exp((x - x1) * log(y2 / y1) / (x2 - x1))
    elseif law == LogLog
        y1 == 0.0 ? 0.0 : y1 * exp(log(x / x1) * log(y2 / y1) / log(x2 / x1))
    elseif law == CoulombPen
        y1 == 0.0 && return 0.0
        thr = Float64(coulomb_threshold)
        t = sqrt(x1 - thr)
        b = log((x2 * y2) / (x1 * y1)) / (1.0 / t - 1.0 / sqrt(x2 - thr))
        a = exp(b / t) * x1 * y1
        (a / x) * exp(-b / sqrt(x - thr))
    else
        error("Unknown interpolation law: $law")
    end
end

terp1(x1, y1, x2, y2, x, i::Integer; coulomb_threshold::Real=0.0) =
    terp1(x1, y1, x2, y2, x, InterpolationLaw(i); coulomb_threshold=coulomb_threshold)

# --- Interval lookup ---

function _find_interval(tab::TabulatedFunction, x::Real)
    n = length(tab)
    x <= tab.x[1] && return (1, tab.interp.law[1])
    x >= tab.x[n] && return (n - 1, tab.interp.law[end])
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >> 1
        tab.x[mid] <= x ? (lo = mid) : (hi = mid)
    end
    law = _law_for_point(tab.interp, lo)
    return (lo, law)
end

function _law_for_point(interp::InterpolationTable, idx::Int)
    for r in eachindex(interp.nbt)
        idx < interp.nbt[r] && return interp.law[r]
    end
    interp.law[end]
end

"""
    interpolate(tab::TabulatedFunction, x; coulomb_threshold=0.0) -> Float64

Interpolate the tabulated function at `x`. Returns 0.0 outside range.
The `coulomb_threshold` keyword is passed through to `terp1` for law 6.
"""
function interpolate(tab::TabulatedFunction, x::Real;
                     coulomb_threshold::Real=0.0)
    n = length(tab)
    n >= 1 || return 0.0
    (x < tab.x[1] || x > tab.x[n]) && return 0.0
    x == tab.x[1] && return Float64(tab.y[1])
    x == tab.x[n] && return Float64(tab.y[n])
    i, law = _find_interval(tab, x)
    terp1(tab.x[i], tab.y[i], tab.x[i+1], tab.y[i+1], x, law;
          coulomb_threshold=coulomb_threshold)
end

# --- Analytical panel integral (NJOY's gral) ---
# Helpers for fallback chains used by log-log law

function _gral_linlin(yl, yh, xl, xh, x1, x2)
    b = (yh - yl) / (xh - xl); a = yl - b * xl
    (x2 - x1) * (a + b * (x2 + x1) / 2)
end

function _gral_loglin(yl, yh, xl, xh, x1, x2, brk)
    b = log(yh / yl) / (xh - xl); a = log(yl) - b * xl
    z = (x2 - x1) * b
    abs(z) <= brk ? exp(a + b * x1) * (x2 - x1) * (1 + z * (0.5 + z / 6)) :
                    exp(a + b * x1) * (exp(z) - 1) / b
end

function _gral_linlog(yl, yh, xl, xh, x1, x2, brk)
    b = (yh - yl) / log(xh / xl); z = (x2 - x1) / x1
    base = (x2 - x1) * (yl + b * log(x1 / xl))
    abs(z) <= brk ? base + b * x1 * z^2 * (1 + z * (-1/3 + z * (1/6 - z / 10))) / 2 :
                    base + b * x1 * (1 + (x2 / x1) * (log(x2 / x1) - 1))
end

"""
    panel_integral(xl, yl, xh, yh, x1, x2, law) -> Float64

Analytical integral from `x1` to `x2` of function on `(xl,yl)..(xh,yh)`.
Translation of NJOY's `gral`.
"""
function panel_integral(xl::Real, yl::Real, xh::Real, yh::Real,
                        x1::Real, x2::Real, law::InterpolationLaw)
    x2 == x1 && return 0.0
    x2 > x1 || error("panel_integral: x2 < x1")
    brk = 0.1

    law == Histogram && return (x2 - x1) * yl
    law == LinLin && return _gral_linlin(yl, yh, xl, xh, x1, x2)

    if law == LinLog
        (xl <= 0 || xh <= 0) && return _gral_linlin(yl, yh, xl, xh, x1, x2)
        return _gral_linlog(yl, yh, xl, xh, x1, x2, brk)
    end

    if law == LogLin
        (yl < 0 || yh < 0) && return _gral_linlin(yl, yh, xl, xh, x1, x2)
        return _gral_loglin(yl, yh, xl, xh, x1, x2, brk)
    end

    if law == LogLog
        if xl <= 0 || xh <= 0
            (yl < 0 || yh < 0) && return _gral_linlin(yl, yh, xl, xh, x1, x2)
            return _gral_loglin(yl, yh, xl, xh, x1, x2, brk)
        elseif yl < 0 || yh < 0
            (xl <= 0 || xh <= 0) && return _gral_linlin(yl, yh, xl, xh, x1, x2)
            return _gral_linlog(yl, yh, xl, xh, x1, x2, brk)
        end
        b = log(yh / yl) / log(xh / xl)
        z = (b + 1) * log(x2 / x1)
        if abs(z) <= brk
            return yl * x1 * ((x1 / xl)^b) * log(x2 / x1) * (1 + z * (0.5 + z / 6))
        else
            return yl * x1 * ((x1 / xl)^b) * (((x2 / x1)^(b + 1)) - 1) / (b + 1)
        end
    end

    # CoulombPen fallback to lin-lin
    _gral_linlin(yl, yh, xl, xh, x1, x2)
end

"""
    integrate(tab::TabulatedFunction, x1, x2; coulomb_threshold=0.0) -> Float64

Integrate the tabulated function from `x1` to `x2` using analytical
panel integrals (matching NJOY's `intega`/`gral`). Zero outside range.
The `coulomb_threshold` keyword is available for future use with law 6.
"""
function integrate(tab::TabulatedFunction, x1::Real, x2::Real;
                   coulomb_threshold::Real=0.0)
    n = length(tab)
    n >= 2 || return 0.0
    x1 >= x2 && return 0.0
    a = max(x1, tab.x[1]); b = min(x2, tab.x[n])
    a >= b && return 0.0
    result = 0.0
    i, _ = _find_interval(tab, a)
    xlo = a
    while xlo < b && i < n
        xhi = min(tab.x[i+1], b)
        law = _law_for_point(tab.interp, i)
        result += panel_integral(tab.x[i], tab.y[i], tab.x[i+1], tab.y[i+1],
                                 xlo, xhi, law)
        xlo = xhi; i += 1
    end
    result
end
