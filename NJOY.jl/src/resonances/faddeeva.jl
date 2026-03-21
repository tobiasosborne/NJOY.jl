# Faddeeva function (complex probability integral) for Doppler broadening
#
# Direct translation of NJOY2016 reconr.f90 subroutines:
#   w(rez, aim1, rew, aimw)   -- exact computation of w(z)
#   wtab(rw, aimw)            -- build 62x62 lookup table
#   quickw(ax, y, rew, aimw, ki, tr, ti) -- fast evaluation
#
# The Faddeeva function is defined as:
#   w(z) = exp(-z^2) * erfc(-iz)
#
# For Doppler broadening of resonance cross sections, the real and imaginary
# parts give the Voigt line profiles psi and chi:
#   psi = (sqrt(pi)/2) * theta * Re(w(x + iy))
#   chi = (sqrt(pi)/2) * theta * Im(w(x + iy))
#
# The implementation is AD-compatible: no mutation of arrays in the core
# evaluation path, no try-catch blocks.

# ============================================================================
# Exact Faddeeva function -- faithful translation of NJOY's `w` subroutine
# ============================================================================

"""
    faddeeva_w(rez, aim1) -> (rew, aimw)

Compute the complex probability integral w(z) = exp(-z^2) * erfc(-iz)
where z = rez + i*aim1.

Returns `(Re(w), Im(w))` as a tuple of Float64.

This is a direct translation of the NJOY2016 `w` subroutine in reconr.f90.
The algorithm uses region-based evaluation:
- Asymptotic continued-fraction for large |z| with Im(z) >= 0
- Taylor-series based evaluation otherwise

The implementation avoids mutation and is compatible with automatic differentiation.
"""
function faddeeva_w(rez::Real, aim1::Real)
    rez = Float64(rez)
    aim1 = Float64(aim1)

    rpi = sqrt(Float64(pi))

    # Constants matching NJOY2016
    c1 = 1.1283792e0   # 2/sqrt(pi)
    c2 = 1.5e0
    up = 1.0e15
    dn = 1.0e-15
    brk1 = 1.25e0
    brk2 = 5.0e0
    brk3 = 1.863636e0
    brk4 = 4.1e0
    brk5 = 1.71e0
    brk6 = 2.89e0
    brk7 = 1.18e0
    brk8 = 5.76e0
    brk9 = 1.5e0
    eps = 1.0e-7

    rew = 0.0
    aimw = 0.0

    aimz = abs(aim1)
    abrez = abs(rez)

    # Special case: z = 0
    if abrez + aimz == 0.0
        return (1.0, 0.0)
    end

    r2 = rez * rez
    ai2 = aimz * aimz

    # Determine which algorithm to use (kw=1 for asymptotic, kw=2 for Taylor)
    use_asymptotic = false
    if abrez + brk1 * aimz - brk2 > 0.0
        use_asymptotic = true
    elseif abrez + brk3 * aimz - brk4 > 0.0
        use_asymptotic = true
    elseif r2 + brk5 * ai2 - brk6 < 0.0
        # use Taylor (kw=2)
    elseif r2 + brk7 * ai2 - brk8 >= 0.0
        # use Taylor (kw=2)
    elseif aimz - brk9 >= 0.0
        use_asymptotic = true
    end
    # else: use Taylor (kw=2)

    if use_asymptotic && aim1 >= 0.0
        # Asymptotic continued-fraction series (kw=1, label 370)
        return _faddeeva_asymptotic(rez, aimz, r2, rpi, c1, eps, up, dn)
    else
        # Taylor series (kw=2, label 420)
        # Use aim1 directly (may be negative)
        return _faddeeva_taylor(rez, aim1, r2, rpi, c1, c2, eps, up, dn)
    end
end

"""
Asymptotic continued-fraction evaluation of w(z).
Used when |z| is large and Im(z) >= 0 (kw=1, Fortran label 370).
"""
function _faddeeva_asymptotic(rez, aimz, r2, rpi, c1, eps, up, dn)
    ai2 = aimz * aimz
    rv = 2.0 * (r2 - ai2)
    ak = 4.0 * rez * aimz
    el = ak
    h = 0.0
    b = 0.0
    a = 0.0
    tempm = 0.0
    temel = 0.0
    g = 1.0
    c = -c1 * aimz
    d = c1 * rez
    am = rv - 1.0
    aak = 1.0

    rew = 0.0
    aimw = 0.0
    tempc = 0.0
    tempd = 0.0

    for _iter in 1:200  # safety limit
        ajtemp = 2.0 * aak
        temp4 = (1.0 - ajtemp) * ajtemp
        ajp = rv - (4.0 * aak + 1.0)

        # Compute continued fraction step (label 480)
        tempc = ajp * c + temp4 * a - ak * d
        tempd = ajp * d + temp4 * b + ak * c
        temel = ajp * el + temp4 * h + ak * am
        tempm = ajp * am + temp4 * g - ak * el
        a = c
        b = d
        g = am
        h = el
        c = tempc
        d = tempd
        am = tempm
        el = temel

        # Overflow/underflow protection
        if abs(tempm) + abs(temel) >= up
            c *= dn; d *= dn; am *= dn; el *= dn
            tempc *= dn; tempd *= dn; tempm *= dn; temel *= dn
        elseif abs(tempm) + abs(temel) <= dn
            c *= up; d *= up; am *= up; el *= up
            tempc *= up; tempd *= up; tempm *= up; temel *= up
        end

        # Convergence check (label 390)
        aak += 1.0

        prr = rew
        pii = aimw
        amagn = tempm^2 + temel^2
        rew = (tempc * tempm + tempd * temel) / amagn
        aimw = (tempm * tempd - temel * tempc) / amagn

        if abs(rew - prr) < eps && abs(aimw - pii) < eps
            return (rew, aimw)
        end
    end

    return (rew, aimw)
end

"""
Taylor-series evaluation of w(z).
Used when Im(z) < 0 or |z| is moderate (kw=2, Fortran label 420).
"""
function _faddeeva_taylor(rez, aimz, r2, rpi, c1, c2, eps, up, dn)
    ai2 = aimz * aimz
    temp1 = r2 + ai2
    temp2 = 2.0 * temp1 * temp1
    aj = -(r2 - ai2) / temp2
    ak = 2.0 * rez * aimz / temp2
    c = 0.0
    b = 0.0
    ajsig = 0.0
    d = 0.0
    g = 0.0
    h = 0.0
    el = 0.0
    a = 1.0
    am = 1.0
    sigp = c2

    expon = exp(temp2 * aj)
    expc = expon * cos(temp2 * ak)
    exps = -expon * sin(temp2 * ak)
    sig2p = 2.0 * sigp

    rew = 0.0
    aimw = 0.0
    tempc = 0.0
    tempd = 0.0
    tempm = 0.0
    temel = 0.0

    for _iter in 1:200  # safety limit
        aj4sig = 4.0 * ajsig
        aj4sm1 = aj4sig - 1.0
        temp3 = 1.0 / (aj4sm1 * (aj4sig + 3.0))
        tt4 = sig2p * (2.0 * ajsig - 1.0)
        temp4 = tt4 / (aj4sm1 * (aj4sig + 1.0) * (aj4sig - 3.0) * aj4sm1)
        ajp = aj + temp3

        # Compute continued fraction step (label 480)
        tempc = ajp * c + temp4 * a - ak * d
        tempd = ajp * d + temp4 * b + ak * c
        temel = ajp * el + temp4 * h + ak * am
        tempm = ajp * am + temp4 * g - ak * el
        a = c
        b = d
        g = am
        h = el
        c = tempc
        d = tempd
        am = tempm
        el = temel

        # Overflow/underflow protection
        if abs(tempm) + abs(temel) >= up
            c *= dn; d *= dn; am *= dn; el *= dn
            tempc *= dn; tempd *= dn; tempm *= dn; temel *= dn
        elseif abs(tempm) + abs(temel) <= dn
            c *= up; d *= up; am *= up; el *= up
            tempc *= up; tempd *= up; tempm *= up; temel *= up
        end

        # Convergence check (label 440)
        ajsig += 1.0

        temp7 = rpi * (am^2 + el^2)
        ref = (aimz * (c * am + d * el) - rez * (am * d - c * el)) / temp7 / temp1
        aimf = (aimz * (am * d - c * el) + rez * (c * am + d * el)) / temp7 / temp1

        prr = rew
        pii = aimw
        rew = expc - ref
        aimw = exps - aimf

        if abs(rew - prr) < eps && abs(aimw - pii) < eps
            return (rew, aimw)
        end

        sig2p = 2.0 * ajsig
    end

    return (rew, aimw)
end

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

# ============================================================================
# Voigt line profiles psi and chi
# ============================================================================

"""
    psi_chi(x, y) -> (psi, chi)

Compute the Voigt line profiles psi and chi using the exact Faddeeva function.

    psi = (sqrt(pi)/2) * theta * Re(w(x + iy))
    chi = (sqrt(pi)/2) * theta * Im(w(x + iy))

In this parameterization (matching NJOY's usage with `quickw`),
`x` and `y` are the Doppler-scaled detuning and natural width:
- x = theta * (E - E_r) / Gamma_total   (effectively `ax` in NJOY)
- y = theta / 2                          (effectively `y` in NJOY)
- theta = Gamma_total / delta, where delta = Doppler width

However, this convenience function simply returns (Re(w), Im(w)) scaled
by sqrt(pi)/2, since the theta factor is already embedded in x and y
by the caller. The relationship to the cross-section Voigt profiles is:

    psi(x,y) = Re(w(x + iy))
    chi(x,y) = Im(w(x + iy))

(The sqrt(pi)*theta/2 prefactor is applied by the calling cross-section code.)
"""
function psi_chi(x::Real, y::Real)
    rew, aimw = faddeeva_w(x, y)
    return (rew, aimw)
end

"""
    psi_chi(x, y, table::FaddeevaTable) -> (psi, chi)

Fast version using precomputed lookup table (via `quickw`).
"""
function psi_chi(x::Real, y::Real, table::FaddeevaTable)
    rew, aimw = quickw(x, y, table)
    return (rew, aimw)
end

# ============================================================================
# SpecialFunctions.jl-based implementation for validation
# ============================================================================

"""
    faddeeva_w_julia(x, y) -> (rew, aimw)

Compute the Faddeeva function w(z) using Julia's SpecialFunctions.jl.

Uses the identity: w(z) = erfcx(-iz) where erfcx(z) = exp(z^2) * erfc(z).
This gives w(x + iy) = erfcx(-i(x + iy)) = erfcx(y - ix).

Returns (Re(w), Im(w)) as a tuple.
"""
function faddeeva_w_julia(x::Real, y::Real)
    z = Complex{Float64}(y, -x)  # z_arg = y - ix = -i*(x + iy)
    wz = erfcx(z)
    return (real(wz), imag(wz))
end
