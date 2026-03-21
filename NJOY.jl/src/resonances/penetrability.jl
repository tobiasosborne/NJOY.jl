# Penetrability, shift factor, and phase shift calculations
#
# Direct translation of NJOY2016 reconr.f90 subroutines:
#   facts(l, rho, se, pe)   -- penetrability P_l and shift factor S_l
#   facphi(l, rho, phi)     -- hard-sphere phase shift phi_l
#
# All formulas are analytic for l = 0, 1, 2, 3, 4.
# For l >= 5, recursive relations from Abramowitz & Stegun (10.1.19) are used.
#
# Reference: ENDF-102 Section 2.4; A&S Chapter 10.

"""
    penetrability(l::Integer, rho::Real) -> Float64

Compute the penetrability factor P_l(rho) for orbital angular momentum `l`
and channel radius parameter rho = k * a (wavenumber times channel radius).

Exact formulas for l = 0..4 matching NJOY2016 `facts` subroutine.
Uses upward recursion for l >= 5.

# Formulas (from reconr.f90 `facts`)
- l=0: P_0 = rho
- l=1: P_1 = rho^3 / (1 + rho^2)
- l=2: P_2 = rho^5 / (9 + 3*rho^2 + rho^4)
- l=3: P_3 = rho^7 / (225 + 45*rho^2 + 6*rho^4 + rho^6)
- l=4: P_4 = rho^9 / (11025 + 1575*rho^2 + 135*rho^4 + 10*rho^6 + rho^8)
"""
function penetrability(l::Integer, rho::Real)
    r2 = rho * rho

    if l == 0
        return rho
    elseif l == 1
        den = 1.0 + r2
        return r2 * rho / den
    elseif l == 2
        r4 = r2 * r2
        den = 9.0 + 3.0 * r2 + r4
        return r4 * rho / den
    elseif l == 3
        r4 = r2 * r2
        r6 = r4 * r2
        den = 225.0 + 45.0 * r2 + 6.0 * r4 + r6
        return r6 * rho / den
    elseif l == 4
        r4 = r2 * r2
        r6 = r4 * r2
        r8 = r4 * r4
        den = 11025.0 + 1575.0 * r2 + 135.0 * r4 + 10.0 * r6 + r8
        return r8 * rho / den
    else
        # Upward recursion: P_l = rho^2 * P_{l-1} / ((2l-1)^2 * (1 - S_{l-1}/(l))^2 + ... )
        # Use the recursion from spherical Bessel functions:
        # P_l = rho / ((l/P_{l-1} - S_{l-1})^2 + (l/P_{l-1})^2 ... )
        # More robustly: compute via the continued fraction representation
        # For safety, build up from l=4
        P_prev = penetrability(l - 1, rho)
        S_prev = shift_factor(l - 1, rho)
        # Recursion: P_l = rho^2 * P_{l-1} / ((l - S_{l-1})^2 + P_{l-1}^2)
        # This is the standard recursion for Coulomb/hard-sphere functions
        return r2 * P_prev / ((l - S_prev)^2 + P_prev^2)
    end
end

"""
    shift_factor(l::Integer, rho::Real) -> Float64

Compute the shift factor S_l(rho) for orbital angular momentum `l`.

Exact formulas for l = 0..4 matching NJOY2016 `facts` subroutine.
Uses upward recursion for l >= 5.

# Formulas (from reconr.f90 `facts`)
- l=0: S_0 = 0
- l=1: S_1 = -1 / (1 + rho^2)
- l=2: S_2 = -(18 + 3*rho^2) / (9 + 3*rho^2 + rho^4)
- l=3: S_3 = -(675 + 90*rho^2 + 6*rho^4) / (225 + 45*rho^2 + 6*rho^4 + rho^6)
- l=4: S_4 = -(44100 + 4725*rho^2 + 270*rho^4 + 10*rho^6) /
              (11025 + 1575*rho^2 + 135*rho^4 + 10*rho^6 + rho^8)
"""
function shift_factor(l::Integer, rho::Real)
    r2 = rho * rho

    if l == 0
        return 0.0
    elseif l == 1
        den = 1.0 + r2
        return -1.0 / den
    elseif l == 2
        r4 = r2 * r2
        den = 9.0 + 3.0 * r2 + r4
        return -(18.0 + 3.0 * r2) / den
    elseif l == 3
        r4 = r2 * r2
        r6 = r4 * r2
        den = 225.0 + 45.0 * r2 + 6.0 * r4 + r6
        return -(675.0 + 90.0 * r2 + 6.0 * r4) / den
    elseif l == 4
        r4 = r2 * r2
        r6 = r4 * r2
        r8 = r4 * r4
        den = 11025.0 + 1575.0 * r2 + 135.0 * r4 + 10.0 * r6 + r8
        return -(44100.0 + 4725.0 * r2 + 270.0 * r4 + 10.0 * r6) / den
    else
        # Recursion: S_l = rho^2 * (l - S_{l-1}) / ((l - S_{l-1})^2 + P_{l-1}^2) - l
        P_prev = penetrability(l - 1, rho)
        S_prev = shift_factor(l - 1, rho)
        return r2 * (l - S_prev) / ((l - S_prev)^2 + P_prev^2) - l
    end
end

"""
    phase_shift(l::Integer, rho::Real) -> Float64

Compute the hard-sphere phase shift phi_l(rho) for orbital angular momentum `l`.

Exact formulas for l = 0..4 matching NJOY2016 `facphi` subroutine.
Uses upward recursion for l >= 5.

# Formulas (from reconr.f90 `facphi`)
- l=0: phi_0 = rho
- l=1: phi_1 = rho - atan(rho)
- l=2: phi_2 = rho - atan(3*rho / (3 - rho^2))
- l=3: phi_3 = rho - atan((15*rho - rho^3) / (15 - 6*rho^2))
- l=4: phi_4 = rho - atan((105*rho - 10*rho^3) / (105 - 45*rho^2 + rho^4))

A small-angle test (phi/rho < 1e-6) sets the result to zero for l >= 2,
matching the Fortran behavior (reconr.f90 facphi uses `(phi/rho).lt.test`
without abs(), so negative phi/rho also triggers the cutoff).
"""
function phase_shift(l::Integer, rho::Real)
    r2 = rho * rho
    test = 1.0e-6

    if l == 0
        return rho
    elseif l == 1
        return rho - atan(rho)
    elseif l == 2
        phi = rho - atan(3.0 * rho / (3.0 - r2))
        # Matches Fortran facphi (reconr.f90:4459): (phi/rho).lt.test (no abs)
        if phi / rho < test
            phi = 0.0
        end
        return phi
    elseif l == 3
        phi = rho - atan((15.0 * rho - rho * r2) / (15.0 - 6.0 * r2))
        # Matches Fortran facphi (reconr.f90:4462): (phi/rho).lt.test (no abs)
        if phi / rho < test
            phi = 0.0
        end
        return phi
    elseif l == 4
        r4 = r2 * r2
        top = 105.0 * rho - 10.0 * r2 * rho
        bot = 105.0 - 45.0 * r2 + r4
        phi = rho - atan(top / bot)
        # Matches Fortran facphi (reconr.f90:4468): (phi/rho).lt.test (no abs)
        if phi / rho < test
            phi = 0.0
        end
        return phi
    else
        # Recursion: phi_l = phi_{l-1} + atan(P_{l-1} / (l - S_{l-1}))
        # Actually use: phi_l = rho - l*pi/2 + delta_l for Coulomb,
        # but for hard sphere, the recursion from A&S 10.1.19:
        # tan(phi_l) can be built from spherical Bessel recursion.
        # Simpler: phi_l = atan(rho * j_l(rho), -rho * n_l(rho))
        # For now, use recursion via P and S:
        P_prev = penetrability(l - 1, rho)
        S_prev = shift_factor(l - 1, rho)
        phi_prev = phase_shift(l - 1, rho)
        phi = phi_prev + atan(P_prev / (l - S_prev))
        # Matches Fortran facphi: (phi/rho).lt.test (no abs)
        if rho != 0.0 && phi / rho < test
            phi = 0.0
        end
        return phi
    end
end
