# Faddeeva lookup table and fast evaluation (quickw)
#
# Direct translation of NJOY2016 reconr.f90 subroutines:
#   wtab(rw, aimw)            -- build 62x62 lookup table
#   quickw(ax, y, rew, aimw, ki, tr, ti) -- fast evaluation

# ============================================================================
# Faddeeva lookup table
# ============================================================================

"""
    FaddeevaTable

Precomputed 62x62 lookup table of the complex probability integral w(z),
used by `quickw` for fast evaluation when x^2 + y^2 < 36.

The table covers x in [-0.1, 6.0] and y in [-0.1, 6.0] at spacing 0.1.
This matches NJOY2016's `wtab` subroutine exactly.
"""
struct FaddeevaTable
    tr::NTuple{3844, Float64}   # 62*62 real parts, column-major
    ti::NTuple{3844, Float64}   # 62*62 imaginary parts, column-major
end

# Linear index into the 62x62 table (column-major, 1-based)
@inline _tindex(i, j) = i + 62 * (j - 1)

"""
    build_faddeeva_table() -> FaddeevaTable

Construct the 62x62 lookup table of w(z) values.
Direct translation of NJOY2016's `wtab` subroutine.

Grid: x[i] = -0.1 + (i-1)*0.1, y[j] = -0.1 + (j-1)*0.1 for i,j = 1..62.
"""
function build_faddeeva_table()
    dx = 0.1
    dy = 0.1
    x0 = -0.1
    y0 = -0.1

    # Build vectors of grid values
    xvals = ntuple(i -> x0 + (i - 1) * dx, Val(62))
    yvals = ntuple(j -> y0 + (j - 1) * dy, Val(62))

    # Evaluate w at every grid point
    tr_vec = Vector{Float64}(undef, 3844)
    ti_vec = Vector{Float64}(undef, 3844)

    for j in 1:62
        yj = yvals[j]
        for i in 1:62
            xi = xvals[i]
            rwt, aimt = faddeeva_w(xi, yj)
            idx = _tindex(i, j)
            tr_vec[idx] = rwt
            ti_vec[idx] = aimt
        end
    end

    return FaddeevaTable(NTuple{3844, Float64}(tr_vec),
                         NTuple{3844, Float64}(ti_vec))
end

# ============================================================================
# Fast Faddeeva evaluation (quickw)
# ============================================================================

"""
    quickw(ax, y, table::FaddeevaTable; compute_imaginary=true) -> (rew, aimw)

Fast evaluation of the complex probability integral w(z) for z = ax + iy.
Direct translation of NJOY2016's `quickw` subroutine.

Uses four evaluation regions based on `test = x^2 + y^2`:
- `test < 36`: bilinear interpolation from precomputed 62x62 table
- `36 <= test < 144`: two-pole rational approximation
- `144 <= test < 10000`: one-pole rational approximation
- `test >= 10000`: asymptotic 1/(sqrt(pi)*(x^2+y^2))

When `compute_imaginary=false`, only the real part is computed (the imaginary
part is returned as 0.0), matching `ki=0` in the Fortran code.

The function is AD-compatible: no mutation, no try-catch.
"""
function quickw(ax::Real, y::Real, table::FaddeevaTable;
                compute_imaginary::Bool=true)
    ax = Float64(ax)
    y = Float64(y)
    rpi = sqrt(Float64(pi))

    # Constants from NJOY2016
    break1 = 36.0
    break2 = 144.0
    break3 = 10000.0
    c1 = 0.2752551
    c2 = 2.724745
    c3 = 0.5124242
    c4 = 0.05176536
    c5 = 1.1283792

    aki = ax < 0.0 ? -1.0 : 1.0
    x = abs(ax)
    test = x * x + y * y

    if test < break1
        return _quickw_table(x, y, aki, table, compute_imaginary)
    elseif test < break2
        return _quickw_rational2(x, y, aki, c1, c2, c3, c4, compute_imaginary)
    elseif test < break3
        return _quickw_rational1(x, y, aki, c5, compute_imaginary)
    else
        return _quickw_asymptotic(x, y, aki, rpi, test, compute_imaginary)
    end
end

"""
Table interpolation for quickw (test < 36).
Uses 6-point bilinear interpolation on the 62x62 grid.
"""
@inline function _quickw_table(x, y, aki, table, compute_imaginary)
    ii = unsafe_trunc(Int, x * 10)
    jj = unsafe_trunc(Int, y * 10)
    i = ii + 2
    j = jj + 2
    n = j - 1
    p = 10.0 * x - ii
    q = 10.0 * y - jj
    p2 = p * p
    q2 = q * q
    pq = p * q
    hp = p / 2.0
    hq = q / 2.0
    hq2 = q2 / 2.0
    hp2 = p2 / 2.0
    a1 = hq2 - hq
    a2 = hp2 - hp
    a3 = 1.0 + pq - p2 - q2
    a4 = hp2 - pq + hp
    a5 = hq2 - pq + hq

    # Table lookups (column-major indexing)
    tr = table.tr
    rew = a1 * tr[_tindex(i, n)] + a2 * tr[_tindex(i - 1, j)] +
          a3 * tr[_tindex(i, j)] + a4 * tr[_tindex(i + 1, j)] +
          a5 * tr[_tindex(i, j + 1)] + pq * tr[_tindex(i + 1, j + 1)]

    if compute_imaginary
        ti = table.ti
        aimw = a1 * ti[_tindex(i, n)] + a2 * ti[_tindex(i - 1, j)] +
               a3 * ti[_tindex(i, j)] + a4 * ti[_tindex(i + 1, j)] +
               a5 * ti[_tindex(i, j + 1)] + pq * ti[_tindex(i + 1, j + 1)]
        aimw *= aki
        return (rew, aimw)
    else
        return (rew, 0.0)
    end
end

"""
Two-pole rational approximation for quickw (36 <= test < 144).
"""
@inline function _quickw_rational2(x, y, aki, c1, c2, c3, c4, compute_imaginary)
    a1 = x^2 - y^2
    a2 = 2.0 * x * y
    a3 = a2^2
    a4 = a1 - c1
    a5 = a1 - c2
    d1 = c3 / (a4^2 + a3)
    d2 = c4 / (a5^2 + a3)
    rew = d1 * (a2 * x - a4 * y) + d2 * (a2 * x - a5 * y)

    if compute_imaginary
        aimw = d1 * (a4 * x + a2 * y) + d2 * (a5 * x + a2 * y)
        aimw *= aki
        return (rew, aimw)
    else
        return (rew, 0.0)
    end
end

"""
One-pole rational approximation for quickw (144 <= test < 10000).
"""
@inline function _quickw_rational1(x, y, aki, c5, compute_imaginary)
    a1 = (x^2 - y^2) * 2.0
    a2 = 4.0 * x * y
    a4 = a1 - 1.0
    d1 = c5 / (a4^2 + a2^2)
    rew = d1 * (a2 * x - a4 * y)

    if compute_imaginary
        aimw = d1 * (a4 * x + a2 * y)
        aimw *= aki
        return (rew, aimw)
    else
        return (rew, 0.0)
    end
end

"""
Asymptotic approximation for quickw (test >= 10000).
"""
@inline function _quickw_asymptotic(x, y, aki, rpi, test, compute_imaginary)
    a1 = 1.0 / (rpi * test)
    rew = y * a1

    if compute_imaginary
        aimw = x * a1 * aki
        return (rew, aimw)
    else
        return (rew, 0.0)
    end
end

"""
    psi_chi(x, y, table::FaddeevaTable) -> (psi, chi)

Fast version using precomputed lookup table (via `quickw`).
"""
function psi_chi(x::Real, y::Real, table::FaddeevaTable)
    rew, aimw = quickw(x, y, table)
    return (rew, aimw)
end
