# Faddeeva function (complex probability integral) -- exact evaluation
# Direct translation of NJOY2016 reconr.f90 w(rez, aim1, rew, aimw).
# w(z) = exp(-z^2) * erfc(-iz); AD-compatible (no mutation, no try-catch).

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
