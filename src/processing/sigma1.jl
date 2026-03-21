# SIGMA1 kernel: pure mathematical functions for Doppler broadening (Proposal B)
# Correspondence: funky->f_func/f_all, hunky->h_func/h_all, hnabb->h_taylor,
#                 bsigma->sigma1_at
# All functions are pure (no mutation), AD-compatible.

using SpecialFunctions: erfc
using StaticArrays: MVector

# Exact f_n(0) values: f_0=1/2, f_1=1/(2sqrt(pi)), f_2=1/4, f_3=f_1, f_4=3/8
const _F_ZERO = (0.5, 1.0/(2.0*sqrt(pi)), 0.25, 1.0/(2.0*sqrt(pi)), 0.375)
const _H_CANCEL_TOL = 1.0e-5   # switch to Taylor when |h|/|f| < this
const _H_SMALL_ABS  = 1.0e-12  # treat as zero when |h|/|f| < this
const _POW2 = (sqrt(2.0), 2.0, 2*sqrt(2.0), 4.0, 4*sqrt(2.0))

"""
    f_func(n::Int, a) -> Real
F-function: f_n(a) = integral(a..inf) z^n exp(-z^2) dz / sqrt(pi), n=0..4.
Uses erfc for f_0 and the recursion f_{n+2} = ((n+1)*f_n + a^{n+1}*exp(-a^2)/sqrt(pi))/2.
"""
function f_func(n::Int, a)
    T = float(typeof(a)); rsqpi = inv(T(sqrt(pi))); alim = T(10)
    expo = abs(a) < alim ? exp(-a*a) : zero(T)
    f0 = abs(a) < alim ? erfc(a)/2 : zero(T)
    n == 0 && return f0
    e1 = expo*rsqpi/2;  n == 1 && return e1
    ae = a*expo*rsqpi; f2 = (f0 + ae)/2;  n == 2 && return f2
    a2e = a*ae; f3 = (2*e1 + a2e)/2;  n == 3 && return f3
    return (3*f2 + a*a2e)/2  # f4
end

"""
    f_all(a) -> NTuple{5}
Return (f_0..f_4) in a single pass (one exp + one erfc call).
"""
function f_all(a)
    T = float(typeof(a)); rsqpi = inv(T(sqrt(pi))); alim = T(10)
    expo = abs(a) < alim ? exp(-a*a) : zero(T)
    f0 = abs(a) < alim ? erfc(a)/2 : zero(T)
    ae = expo*rsqpi; f1 = ae/2
    cum = a*ae;  f2 = (f0 + cum)/2
    cum = a*cum; f3 = (2*f1 + cum)/2
    cum = a*cum; f4 = (3*f2 + cum)/2
    return (f0, f1, f2, f3, f4)
end

"""
    h_func(n::Int, a, b) -> Real
H-function: h_n(a,b) = f_n(a) - f_n(b). Falls back to Taylor series
when cancellation is detected (|diff/f| < 1e-5).
"""
function h_func(n::Int, a, b)
    fa = f_func(n, a); fb = f_func(n, b); diff = fa - fb
    a == b && return zero(diff)
    afa = abs(fa)
    afa > 0 && abs(diff) <= _H_SMALL_ABS*afa && return zero(diff)
    afa > 0 && abs(diff) < _H_CANCEL_TOL*afa && return h_taylor(n, a, b)
    return diff
end

"""
    h_all(a_old, f_old::NTuple{5}, a_new) -> (h::NTuple{5}, f_new::NTuple{5})
Compute all five H-functions and return new F-values. Stateless replacement
for NJOY's hunky (which mutates module-level arrays).
"""
function h_all(a_old, f_old::NTuple{5}, a_new)
    f_new = f_all(a_new)
    h = ntuple(5) do k
        diff = f_old[k] - f_new[k]; afa = abs(f_old[k])
        if afa > 0 && abs(diff) <= _H_SMALL_ABS*afa;  zero(diff)
        elseif afa > 0 && abs(diff) < _H_CANCEL_TOL*afa && a_old != a_new
            h_taylor(k-1, a_old, a_new)
        else; diff; end
    end
    return (h, f_new)
end

"""
    h_taylor(n::Int, aa, bb) -> Real
Direct Taylor series for h_n(a,b) when b-a is small. Matches NJOY's hnabb
with polynomial coefficient recursion and convergence to 1e-8 relative.
"""
function h_taylor(n::Int, aa, bb)
    T = float(promote_type(typeof(aa), typeof(bb)))
    rsqpi = inv(T(sqrt(pi))); aerr = T(1e30); rerr = T(1e-8); explim = T(100)
    sign_val = one(T)
    if bb < aa; a = abs(bb); b = abs(aa); sign_val = -sign_val
    else; a = abs(aa); b = abs(bb); end
    bb < zero(T) && isodd(n) && (sign_val = -sign_val)
    h_step = (b - a)*_POW2[1]; x = T(_POW2[1])*a; xx = x*x; asq = a*a
    con = asq < explim ? exp(-asq)*rsqpi/T(_POW2[n+1]) : zero(T)
    # Use MVector for stack-friendly fixed-size mutable arrays (no heap alloc)
    cm = MVector{50,T}(ntuple(_ -> zero(T), Val(50)))
    cmstar = MVector{50,T}(ntuple(_ -> zero(T), Val(50)))
    k = n; kd = 0; cm[1] = one(T)
    xk = k == 0 ? one(T) : x^k
    s = h_step*xk; fact = h_step; mflag = 0
    for m in 2:50
        fact *= h_step/m; kstar = k; kdstar = kd + 1
        for j in 1:kdstar; cmstar[j] = cm[j]; end
        k = n - m + 1
        k < 0 && (k = iseven(k) ? 0 : 1)
        kd = div(n + m - 1 - k, 2)
        cm[kd+1] = -cmstar[kdstar]; qmn = cm[kd+1]
        if kd > 0
            for j in 1:kd
                jalpha = div(2j + k - 1 - kstar, 2)
                beta = jalpha > 0 ? cmstar[jalpha] : zero(T)
                cm[j] = (2j + k - 1)*cmstar[jalpha+1] - beta
            end
            for j in 1:kd; qmn = qmn*xx + cm[kd+1-j]; end
        end
        xk = k == 0 ? one(T) : x^k
        term = fact*xk*qmn; s += term
        m < max(T(n+1), abs(h_step*x)) && continue
        # Match Fortran: test = aerr + rerr*abs(s), where aerr=1e30
        # makes convergence governed by minimum iteration count, not term size.
        test = aerr + rerr*abs(s)
        if abs(term) <= test
            mflag == 1 && return con*s*sign_val; mflag = 1
        else; mflag = 0; end
    end
    return con*s*sign_val
end

"""
    interval_contributions(h::NTuple{5}, oy, yy, xx) -> (s1, s2)
Compute s1, s2 integral contributions from one interval (NJOY's hunky).
"""
@inline function interval_contributions(h::NTuple{5}, oy, yy, xx)
    s1 = (h[3]*oy + 2*h[2])*oy + h[1]
    s2 = ((h[5] + (6*yy - xx)*h[3])*oy + (4*h[4] + (4*yy - 2*xx)*h[2]))*oy +
         (yy - xx)*h[1]
    return (s1, s2)
end

"""
    sigma1_at(E, seg_e, seg_xs, alpha) -> Float64
Evaluate Doppler-broadened cross section at energy E using SIGMA1 kernel.
`seg_e`: sorted segment energies [eV], `seg_xs`: XS values, `alpha`=AWR/(bk*T).
Two-pass structure: segments below E, then above E, with 1/v and constant
extrapolation at boundaries.
"""
function sigma1_at(E, seg_e::AbstractVector, seg_xs::AbstractVector, alpha)
    nseg = length(seg_e); atop = 4.0
    y = sqrt(alpha*E); oy_neg = -1.0/y; yy = y*y
    v = map(e -> sqrt(alpha*e), seg_e)
    k = clamp(searchsortedlast(v, y), 1, nseg - 1)
    sbt = 0.0

    # Pass 1: intervals BELOW evaluation point
    f_cur = _F_ZERO; alast = 0.0
    for l in k:-1:1
        x_lo = v[l]; x_hi = v[l+1]
        xx_hi = x_hi*x_hi; xx_lo = x_lo*x_lo
        xx_hi <= xx_lo && continue
        aa = y - x_lo
        h_vals, f_cur = h_all(alast, f_cur, aa); alast = aa
        denom = 1.0/(xx_hi - xx_lo)
        slope = (seg_xs[l+1] - seg_xs[l])*denom
        s1, s2 = interval_contributions(h_vals, oy_neg, yy, xx_hi)
        sbt += seg_xs[l+1]*s1 + slope*s2
        aa > atop && @goto below_done
    end
    # 1/v extrapolation to x=0
    aa = y; h_vals, f_cur = h_all(alast, f_cur, aa)
    sbt -= seg_xs[1]*v[1]*(oy_neg^2*h_vals[2] + oy_neg*h_vals[1])
    @label below_done

    # Pass 2: intervals ABOVE evaluation point
    f_cur = _F_ZERO; alast = 0.0; oy_pos = 1.0/y
    if k < nseg
        for l in k:(nseg-1)
            x_hi = v[l+1]; x_lo = v[l]
            xx_lo_sq = x_lo*x_lo; xx_hi_sq = x_hi*x_hi
            xx_lo_sq >= xx_hi_sq && continue
            aa = x_hi - y
            h_vals, f_cur = h_all(alast, f_cur, aa); alast = aa
            denom = 1.0/(xx_hi_sq - xx_lo_sq)
            slope = (seg_xs[l+1] - seg_xs[l])*denom
            s1, s2 = interval_contributions(h_vals, oy_pos, yy, xx_lo_sq)
            sbt += seg_xs[l]*s1 + slope*s2
            aa > atop && @goto above_done
        end
    end
    # Constant extrapolation to infinity
    factor = (f_cur[3]*oy_pos + 2*f_cur[2])*oy_pos + f_cur[1]
    sbt += seg_xs[end]*factor
    @label above_done

    # Pass 3: low-energy (negative velocity) contribution
    # Fortran broadr.f90 lines 1626-1660.
    # Computes the contribution from target nuclei with negative relative velocity.
    # Key Fortran variable trace:
    #   y = -y_orig (line 1628), aa = y_orig (line 1629)
    #   oy = -1/y_orig (line 1634: oy=-oy from oy_pos=+1/y_orig)
    #   yy remains y_orig^2 (set once at line 1558, never changed)
    aa_init = y  # = y_orig; Fortran: aa = -y where y was just negated to -y_orig
    if aa_init <= atop
        # Extend cross section to x=0.0 as 1/v
        f_cur3 = f_all(aa_init)     # funky(aa) at aa = y_orig
        alast3 = aa_init
        oy3 = -1.0/y               # Fortran oy = -1/y_orig at this point

        # 1/v extrapolation below lowest data point
        # aa = e(klow) - y = v[1] - (-y_orig) = v[1] + y_orig
        y_neg = -y  # = -y_orig, this is the "y" in pass 3
        aa = v[1] - y_neg           # = v[1] + y
        xx_zero = 0.0
        h_vals3, f_cur3 = h_all(alast3, f_cur3, aa)
        alast3 = aa
        sbt -= seg_xs[1]*v[1]*(oy3^2*h_vals3[2] + oy3*h_vals3[1])
        if aa <= atop
            # Loop over all intervals from klow to khigh-1
            for l in 1:(nseg-1)
                x = v[l+1]
                aa = x - y_neg    # = x + y_orig (always positive and increasing)
                xx = v[l]^2
                xx == x*x && continue
                h_vals3, f_cur3 = h_all(alast3, f_cur3, aa)
                alast3 = aa
                denom = 1.0/(x*x - xx)
                slope = (seg_xs[l+1] - seg_xs[l])*denom
                # yy = y_orig^2 (unchanged from pass 1/2)
                s1, s2 = interval_contributions(h_vals3, oy3, yy, xx)
                sbt -= seg_xs[l]*s1 + slope*s2
                aa > atop && @goto neg_done
            end
            # Extend cross section to infinity as constant
            factor = (f_cur3[3]*oy3 + 2*f_cur3[2])*oy3 + f_cur3[1]
            sbt -= seg_xs[end]*factor
        end
    end
    @label neg_done
    # Match Fortran sigmin=1e-15: clamp values below sigmin to zero
    return sbt < 1e-15 ? 0.0 : sbt
end
