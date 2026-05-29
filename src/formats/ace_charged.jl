# ACE charged-particle elastic angular block (LAND/AND).
#
# Port of Fortran `acecpe` (njoy-reference/src/acefc.f90:6492-6671).
# Converts MF6/MT2 LAW=5 (identical-particle Coulomb + nuclear interference,
# tabulated in LIST records per incident energy) into ACE LAW=14 tabulated
# (μ, pdf, cdf) triples + a Coulomb-corrected elastic cross section.
#
# Algorithm summary:
#   For each incident energy E:
#     - Read the LIST mu/prob pairs (LTP=12 form: prob = signi/xelas).
#     - At each μ, evaluate the Coulomb cross section sigc analytically.
#     - signi (nuclear differential) = prob * xelas at this E.
#     - sigtot(μ) = signi + sigc.
#     - Trapezoidally integrate sigtot from -1 (or first μ) to current μ
#       to get cumulative cdf-like quantity `cumm`.
#     - Adaptively insert midpoints when sigc grows >2× and dominates
#       signi (the Coulomb singularity at μ=1 needs finer sampling).
#     - Heating contribution eht += weighted (1-μ) integral.
#   After all incident energies:
#     - Log-log interpolate yys = cumm*2π onto the full ESZ grid → new
#       elastic_xs (now includes Coulomb).
#     - Update total_xs accordingly.

"""
Constants used by acecpe (acefc.f90:6510-6512). emev=1e6, fm=1e-13 (cm).
NOTE Fortran has fm=1e-12 — that's because their length unit is cm and
1 fm = 10^-13 cm; the variable named `fm` actually holds 10^-12 because
the formula uses fm² which becomes 10^-24 cm². We follow the Fortran
verbatim for FP fidelity.
"""
const _ACER_EMEV = 1.0e6
const _ACER_FM   = 1.0e-12

"""
    coulomb_sigc(μ, e, awr, awi, izai_z, target_z, spi, lidp) -> Float64

Coulomb (Mott / Rutherford) differential cross section dσ/dμ at cosine μ
and incident energy e (eV). Includes identical-particle interference
when lidp=1.

Ref: njoy-reference/src/acefc.f90:6573-6588 (acecpe sigc formula).
"""
@inline function coulomb_sigc(μ::Float64, e::Float64, awr::Float64, awi::Float64,
                                izai_z::Int, target_z::Int, spi::Float64,
                                lidp::Int)
    ai  = awi * PhysicsConstants.amassn
    at  = awr * PhysicsConstants.amassn
    cc1 = 2 * PhysicsConstants.amu * PhysicsConstants.ev * _ACER_FM^2 /
          PhysicsConstants.hbar^2
    ee  = (PhysicsConstants.ev / 1e7) * (PhysicsConstants.clight / 10)
    cc2 = ee^4 * PhysicsConstants.amu / (2 * PhysicsConstants.hbar^2 *
                                          PhysicsConstants.ev)
    wn  = at * sqrt(cc1 * e * ai) / (ai + at)
    eta = target_z * izai_z * sqrt(cc2 * ai / e)
    if lidp == 0
        return (eta^2 / wn^2) / (1 - μ)^2
    end
    # lidp == 1 (identical particles)
    i2s = round(Int, 2 * spi)
    sign_term = iseven(i2s) ? 1.0 : -1.0   # (-1)^(2*spi)
    omu2 = 1 - μ^2
    return (2 * eta^2 / wn^2) / omu2 *
           ((1 + μ^2) / omu2 +
            sign_term * cos(eta * log((1 + μ) / (1 - μ))) / (2 * spi + 1))
end

# =========================================================================
# LTP<12 Legendre reconstruction (ptlegc + coul)
# =========================================================================

"""
    _legndr_njoy!(p, x, np)

Generate Legendre polynomials `P_0..P_np` at `x` into `p[1..np+1]` by the
exact NJOY recursion. `p[l+1] = P_l(x)`.

The recursion `p(i+2) = h + g - h/(i+1)` with `g = x*p(i+1)`, `h = g - p(i)`
must be evaluated in this exact FP order to be bit-identical with the Fortran.

Ref: njoy-reference/src/mathm.f90:12-35 (legndr).
"""
@inline function _legndr_njoy!(p::Vector{Float64}, x::Float64, np::Int)
    p[1] = 1.0
    p[2] = x
    np < 2 && return p
    @inbounds for i in 1:(np - 1)
        g = x * p[i+1]
        h = g - p[i]
        p[i+2] = h + g - h / (i + 1)
    end
    return p
end

"""
    _coul(x, coeffs, nt, awp, izap, awt, izat, spi, ltp, lidp, e, p) -> Float64

Charged-particle elastic differential cross section `y` at cosine `x`,
reconstructed from the LTP=1 (nuclear-amplitude) or LTP=2 (residual-XS)
Legendre representation of File 6/LAW=5. `coeffs` is the raw LIST data array
(1-based: `coeffs[k]` corresponds to Fortran `c(6+k)`, i.e. Fortran `c(7)` =
`coeffs[1]`). `nt` = NL (Fortran `nint(c(6))`, the LIST CONT N2 field) — it is
NOT a data value and must be passed explicitly, never read from `coeffs`.

`p` is a scratch buffer of length ≥ 2*nt+1 reused across calls.

Returns `y = sigc + sigr (+ sigi)` — the analytic Coulomb part plus the
nuclear (and, for LTP=1, the Coulomb-nuclear interference) part.

NOTE on FP order: this mirrors `coul` (acefc.f90:8136-8137) where
`eta=zt*zi*sqrt(cc2*ai/e)` and `wn=at*sqrt(cc1*ai*e)/(ai+at)` — the product
`cc1*ai*e` differs in IEEE order from `acecpe`'s inline `cc1*e*ai`
(acefc.f90:6581). The subtraction `pmu=cprob-sigc` downstream uses the
`acecpe` sigc, so the two orders must be kept distinct for bit fidelity.

Ref: njoy-reference/src/acefc.f90:8103-8201 (coul).
"""
function _coul(x::Float64, coeffs::Vector{Float64}, nt::Int,
               awp::Float64, izap::Int, awt::Float64, izat::Int,
               spi::Float64, ltp::Int, lidp::Int, e::Float64,
               p::Vector{Float64})
    # acefc.f90:8124-8142 — constants.
    ai  = awp * PhysicsConstants.amassn
    at  = awt * PhysicsConstants.amassn
    zt  = Float64(izat ÷ 1000)   # int(izat/1000)
    zi  = Float64(izap ÷ 1000)   # int(izap/1000)
    i2s = round(Int, 2 * spi)
    cc1 = 2 * PhysicsConstants.amu * PhysicsConstants.ev * _ACER_FM^2 /
          PhysicsConstants.hbar^2
    ee  = (PhysicsConstants.ev / 1e7) * (PhysicsConstants.clight / 10)
    cc2 = ee^4 * PhysicsConstants.amu / (2 * PhysicsConstants.hbar^2 *
                                          PhysicsConstants.ev)
    eta = zt * zi * sqrt(cc2 * ai / e)
    wn  = at * sqrt(cc1 * ai * e) / (ai + at)

    # acefc.f90:8138-8142 — sigc (coul's own FP order).
    sigc = 0.0
    if lidp == 0
        sigc = (eta^2 / wn^2) / (1 - x)^2
    elseif lidp == 1
        sign_term = iseven(i2s) ? 1.0 : -1.0   # (-1)^i2s
        sigc = ((2 * eta^2 / wn^2) / (1 - x^2)) *
               ((1 + x^2) / (1 - x^2) +
                sign_term * cos(eta * log((1 + x) / (1 - x))) / (2 * spi + 1))
    end

    # nt = nint(c(6)) = NL. CRITICAL: c(6) is the LIST CONT field N2 (=NL),
    # NOT a data value — `coeffs` holds only the NW data values starting at
    # Fortran c(7). So nt must be passed in (= sub.nl), never read from coeffs.
    np = 2 * nt
    _legndr_njoy!(p, x, np)            # p[l+1] = P_l(x), l=0..np

    # `c(k)` (Fortran 1-based) maps to coeffs[k-6]; helper for clarity.
    @inline cc(k::Int) = coeffs[k - 6]

    if ltp == 1
        # --acefc.f90:8147-8182 — nuclear amplitude expansion.
        if lidp != 1
            sigr = cc(7) / 2
            for ip in 1:np
                sigr += (2 * ip + 1) * p[ip+1] * cc(ip + 7) / 2
            end
            cs1 = complex(cc(8 + np), cc(9 + np)) / 2
            for it in 1:nt
                cs1 += (2 * it + 1) * p[it+1] *
                       complex(cc(8 + np + 2 * it), cc(9 + np + 2 * it)) / 2
            end
            carg1 = complex(0.0, 1.0) * eta * log((1 - x) / 2)
            sigi  = (-2 * eta / (1 - x)) * real(cs1 * exp(carg1))
            return sigc + sigr + sigi
        else
            # lidp == 1 (T52): real coeffs c(7)..c(7+nt), then nt+1 complex pairs.
            sigr = cc(7) / 2
            for it in 1:nt
                sigr += (4 * it + 1) * p[2*it+1] * cc(it + 7) / 2
            end
            cs1 = complex(cc(8 + nt), cc(9 + nt)) / 2
            cs2 = cs1
            sgn = -1.0
            for it in 1:nt
                term = (2 * it + 1) * p[it+1] *
                       complex(cc(8 + nt + 2 * it), cc(9 + nt + 2 * it)) / 2
                cs1 += term
                cs2 += sgn * term
                sgn = -sgn
            end
            carg1 = complex(0.0, 1.0) * eta * log((1 - x) / 2)
            carg2 = complex(0.0, 1.0) * eta * log((1 + x) / 2)
            sigi  = (-2 * eta / (1 - x * x)) *
                    real(cs1 * (1 + x) * exp(carg1) + cs2 * (1 - x) * exp(carg2))
            return sigc + sigr + sigi
        end
    else
        # --acefc.f90:8184-8198 — ltp=2 residual cross-section expansion.
        if lidp != 1
            sigr = cc(7) / 2
            for ip in 1:np
                sigr += (2 * ip + 1) * p[ip+1] * cc(ip + 7) / 2
            end
            return sigc + sigr
        else
            sigr = cc(7) / 2
            for it in 1:nt
                sigr += (4 * it + 1) * p[2*it+1] * cc(it + 7) / 2
            end
            return sigc + sigr
        end
    end
end

"""
    ptlegc_expand(sub, awp, izai, awt, iza) -> MF6Law5Subsection

Adaptively reconstruct the (μ, cprob) angular table from the LTP=1/2 Legendre
representation, returning a new subsection with `mu`/`prob` populated and
`nl` set to the reconstructed point count (LTP preserved so the downstream
`acecpe_one_incident` still subtracts sigc via the `ltp<12` branch).

Two-phase algorithm (verbatim port):
  1. Adaptive bisection of [x(2), umin] with a stack of ≤ imax=20 points;
     a midpoint is accepted when |coul(mid) - mean| ≤ tol1*|coul(mid)|, with
     forced refinement when the interval > 0.01 or coul jumps by >5×.
  2. Linear-interpolation thinning to tol2=0.01 over windows ≤ 0.5 wide.

Ref: njoy-reference/src/acefc.f90:7980-8101 (ptlegc).
"""
function ptlegc_expand(sub::MF6Law5Subsection,
                       awp::Float64, izai::Int, awt::Float64, iza::Int)
    c    = sub.coeffs
    ltp  = sub.ltp
    lidp = sub.lidp
    spi  = sub.spi
    e    = sub.e
    nt   = sub.nl              # = NL = nint(c(6)); the LIST CONT N2 field, not data

    # Fortran parameters (acefc.f90:7993-8003).
    imax   = 20
    maxang = 4000
    tol1   = 0.001
    tol2   = 0.01
    one    = 1.0
    half   = 0.5
    hund   = 0.01
    umin   = 0.96

    # Scratch Legendre buffer (P_0..P_{2*nt}); reused by every _coul call.
    pbuf = Vector{Float64}(undef, max(2, 2 * nt + 1))
    _coul_at(xx) = _coul(xx, c, nt, awp, izai, awt, iza, spi, ltp, lidp, e, pbuf)

    x = Vector{Float64}(undef, imax)
    y = Vector{Float64}(undef, imax)
    aco   = Vector{Float64}(undef, maxang)
    cprob = Vector{Float64}(undef, maxang)

    # --adaptive reconstruction (acefc.f90:8008-8052).
    ii = 0
    i  = 2
    x[2] = -1.0
    if iza == izai
        x[2] = -umin
    end
    y[2] = _coul_at(x[2])
    x[1] = umin
    y[1] = _coul_at(x[1])

    test = 0.0
    xm = 0.0
    yt = 0.0
    while i > 0
        dy = 0.0
        if i > 1 && i < imax
            dm = x[i-1] - x[i]
            xm = half * (x[i-1] + x[i])
            xm = round_sigfig(xm, 3, 0)
            if xm > x[i] && xm < x[i-1]
                ym = half * (y[i-1] + y[i])
                yt = _coul_at(xm)
                test = tol1 * abs(yt)
                dy = abs(yt - ym)
                if dm > hund
                    dy = 2 * test
                end
                if ym != 0.0 && yt / ym > one * 5
                    dy = 2 * test
                end
                if ym != 0.0 && yt / ym < one / 5
                    dy = 2 * test
                end
            end
        end
        if dy > test
            # not converged — add midpoint to stack and continue.
            i += 1
            x[i] = x[i-1]
            y[i] = y[i-1]
            x[i-1] = xm
            y[i-1] = yt
        else
            # converged — use top point in stack.
            ii += 1
            ii > maxang && error("ptlegc_expand: too many coulomb angles " *
                                 "(maxang=$maxang) at E=$e eV, MT=2")
            aco[ii]   = x[i]
            cprob[ii] = y[i]
            i -= 1
        end
    end
    nn = ii - 1

    # --thin the distribution to the coarser tolerance (acefc.f90:8055-8091).
    i  = 1
    ii = 1
    idone = 0
    dco = 0.0
    jj = 0
    while i < nn - 1 && idone == 0
        j = i + 1
        check = 0.0
        while j < nn + 1 && check <= 0.0 && dco <= one / 2
            j += 1
            jj = j - 1
            dco = aco[j] - aco[i]
            if dco <= one / 2
                k = i
                while k < j - 1 && check <= 0.0
                    k += 1
                    f = (aco[j] - aco[k]) / dco
                    test = f * cprob[i] + (1 - f) * cprob[j]
                    diff = tol2 * abs(cprob[k])
                    check = abs(test - cprob[k]) - diff
                end
            end
        end
        if dco > one / 2 || check > 0.0
            i = jj
            ii += 1
            aco[ii]   = aco[i]
            cprob[ii] = cprob[i]
            dco = 0.0
        else
            idone = 1
        end
    end
    i = nn + 1
    ii += 1
    aco[ii]   = aco[i]
    cprob[ii] = cprob[i]

    # --load the new distribution (acefc.f90:8093-8099): c(5)=2*ii, c(6)=ii,
    #   then (aco,cprob) pairs. Downstream acecpe reads NL=c(6) and the pairs.
    mu_new   = Vector{Float64}(undef, ii)
    prob_new = Vector{Float64}(undef, ii)
    @inbounds for q in 1:ii
        mu_new[q]   = aco[q]
        prob_new[q] = cprob[q]
    end

    MF6Law5Subsection(spi, e, ltp, lidp, ii, mu_new, prob_new, c)
end

"""
    AceCpeIncident

Per-incident-energy result of `acecpe`: tabulated (μ, pdf, cdf) plus the
total integrated nuclear+Coulomb elastic xs (cumm*2π) and the heating
integral contribution (eht).
"""
struct AceCpeIncident
    e::Float64                 # incident E [eV]
    mu::Vector{Float64}
    pdf::Vector{Float64}       # normalized (sigtot/cumm)
    cdf::Vector{Float64}       # normalized cumulative (cumm/cumm) ∈ [0, 1]
    sigtot_int::Float64        # cumm*2π — total elastic xs at this E
    heating::Float64           # eht*2π (for elastic recoil heating; 0 if izai≤2004)
end

"""
    acecpe_one_incident(sub, e, xelas, awr, awi, izai_z, target_z) -> AceCpeIncident

Process one MF6/MT2 LAW=5 sub-section (one incident energy). Returns the
tabulated (μ, pdf, cdf) plus integrated elastic xs at this E.

`xelas` is the original (nuclear-only) elastic cross section at this
incident energy, linearly interpolated from the ESZ grid — the Fortran
acefc.f90:6549-6553 does exactly this lookup before the per-μ loop.
"""
function acecpe_one_incident(sub::MF6Law5Subsection, xelas::Float64,
                              awr::Float64, awi::Float64,
                              izai_z::Int, target_z::Int)
    e   = sub.e
    spi = sub.spi
    lidp = sub.lidp
    nl  = sub.nl
    ltp = sub.ltp

    # Output buffers — capacity nl + adaptive insertions.
    mu_out     = Float64[]
    sigtot_out = Float64[]
    cumm_out   = Float64[]

    cumm   = 0.0
    eht    = 0.0
    amul   = 0.0    # previous μ
    smul   = 0.0    # previous sigtot
    sigcl  = 0.0    # previous sigc
    ratrl  = 0.0    # previous ratr

    # iterp acts as a "do this midpoint, then re-loop the original endpoint"
    # state machine matching Fortran acefc.f90:6566-6624.
    jl = 1
    iterp = 0
    pending_endpoint_mu = NaN
    pending_endpoint_pmu = NaN

    while jl <= nl
        # Get the μ, pmu we'll evaluate this iteration.
        if iterp == 1
            # Midpoint subdivision: amuu = (amul + endpoint_mu)/2,
            # ratr is interpolated; signi computed from interpolated ratr.
            amuu = (amul + pending_endpoint_mu) / 2
            pmu  = pending_endpoint_pmu  # placeholder; not used in iterp branch
        else
            amuu = sub.mu[jl]
            pmu  = sub.prob[jl]
            pending_endpoint_mu  = amuu
            pending_endpoint_pmu = pmu
        end

        # Compute sigc and signi (matches Fortran acefc.f90:6573-6604).
        local sigc, signi, ratr
        retry = true
        while retry
            sigc = coulomb_sigc(amuu, e, awr, awi, izai_z, target_z, spi, lidp)
            local pmu_local = pmu
            if ltp < 12
                pmu_local = pmu - sigc
            end
            if iterp == 1
                # iterp=1 branch: signi = (ratr - 1) * sigc, where `ratr` is
                # the interpolated value carried over from the iteration that
                # triggered subdivision (set below as `ratr = (ratrl+ratr)/2`).
                # Then ratr is updated based on midpoint sigc/signi.
                signi = (ratr - 1) * sigc
            else
                signi = pmu_local * xelas
                if signi < -sigc
                    signi = -sigc
                end
            end
            ratr = (sigc + signi) / sigc
            retry = false
            # Subdivision trigger: if not first point, not already in iterp,
            # and Coulomb dominates AND grew >2× since last point.
            if jl > 1 && iterp == 0 &&
               sigc > abs(signi) && sigc > 2 * sigcl
                iterp = 1
                # Save endpoint state, switch to midpoint with interpolated ratr.
                pending_endpoint_mu  = amuu
                pending_endpoint_pmu = pmu
                amuu = (amul + amuu) / 2
                ratr = (ratrl + ratr) / 2
                # Re-evaluate at midpoint.
                retry = true
            end
        end

        # Trapezoidal cumulative.
        if !isempty(mu_out)
            cumm += (amuu - amul) * (signi + sigc + smul) / 2
            eht  += (amuu - amul) * ((1 - amuu) * (signi + sigc) +
                                      (1 - amul) * smul) / 2
        end

        push!(mu_out,     amuu)
        push!(sigtot_out, signi + sigc)
        push!(cumm_out,   cumm)

        amul  = amuu
        sigcl = sigc
        smul  = signi + sigc
        ratrl = ratr

        if iterp == 1
            # We just processed the midpoint — exit iterp and re-do the
            # original endpoint at the SAME jl in the next outer iteration.
            iterp = 0
        else
            jl += 1
        end
    end

    if cumm == 0
        return AceCpeIncident(e, mu_out, fill(0.0, length(mu_out)),
                                fill(0.0, length(mu_out)), 0.0, 0.0)
    end

    # Normalize and apply Fortran's sigfig rounding (acefc.f90:6634-6636):
    # μ and pdf at 7 sigfigs, cdf at 9 sigfigs.
    mu_rounded  = [round_sigfig(mu_out[i],     7, 0) for i in eachindex(mu_out)]
    pdf = [round_sigfig(sigtot_out[i] / cumm, 7, 0) for i in eachindex(sigtot_out)]
    cdf = [round_sigfig(cumm_out[i]   / cumm, 9, 0) for i in eachindex(cumm_out)]
    AceCpeIncident(e, mu_rounded, pdf, cdf, cumm * 2π, eht * 2π)
end

"""
    acer_charged_elastic(subs, esz_e, esz_elastic_xs, awr, awi, izai, za)
        -> (angular_block::AngularBlock,
            new_total_xs::Vector{Float64},
            new_elastic_xs::Vector{Float64},
            heating::Vector{Float64})

Build the LAND/AND-equivalent AngularBlock for charged-particle elastic
and the Coulomb-corrected total/elastic XS on the ESZ grid. Mirrors
Fortran acefc.f90:6494-6671 (acecpe + post-loop ESZ adjustment).

`esz_elastic_xs` is the original nuclear-only elastic XS column from
the PENDF (= MT=2 MF3 values interpolated onto ESZ). Output:
  - angular_block: NE incident energies + tabulated (μ, pdf, cdf) at each
  - new_elastic_xs[j] = signow at esz_e[j], log-log-interpolated from
    yys = sigtot_int across NE.
  - new_total_xs[j] is the input total minus original elastic plus new
    elastic (so MT=1 reflects the Coulomb correction).

For incident α + α (T50: izai=2004=za): heating numbers are forced to
zero (acefc.f90:6628 `if (izai.eq.nint(za))`).
"""
function acer_charged_elastic(subs::Vector{MF6Law5Subsection},
                                esz_e::Vector{Float64},
                                esz_total_xs::Vector{Float64},
                                esz_elastic_xs::Vector{Float64},
                                awr::Float64, awi::Float64,
                                izai::Int, za::Int)
    izai_z   = izai ÷ 1000
    target_z = za   ÷ 1000

    # 1. Per-incident-E processing.
    # The ACE ESZ elastic column is sigfig-7-rounded BEFORE it is read by
    # acecpe — Fortran acefc.f90:5497 (`s=sigfig(s,7,0)`) then 5499
    # (`xss(ie+j)=s`). acecpe's xelas interpolation (6549-6553) and the
    # post-loop `signi=xss(esz+3*nes+j)` (6659) both read THIS rounded
    # column, not the raw PENDF elastic. Build it once, here.
    # (Without this, xelas at incident energies that fall on/near a grid
    #  point whose raw elastic rounds differently at the 7th sigfig shifts
    #  signi=pmu*xelas, propagating into cumm and the interpolated signow —
    #  this was the T50 ESZ 7th-sigfig divergence at E=4.0e6 and 1.29e7.)
    esz_elastic_sig7 = [round_sigfig(x, 7, 0) for x in esz_elastic_xs]

    results = AceCpeIncident[]
    for sub0 in subs
        # LTP<12: reconstruct the (μ,cprob) table from the Legendre coefficients
        # BEFORE acecpe runs, exactly as Fortran calls ptlegc per incident energy
        # (acefc.f90:6541-6545). ptlegc is passed (c, awp=awi, izai, awt=awr,
        # iza=za, spi). LTP=12 is already tabulated — pass through untouched.
        sub = sub0.ltp < 12 ?
              ptlegc_expand(sub0, awi, izai, awr, za) : sub0
        # Linear-interp xelas from the sigfig-7 ESZ elastic (acefc.f90:6549-6553).
        xelas = _esz_lin_interp(esz_e, esz_elastic_sig7, sub.e)
        push!(results, acecpe_one_incident(sub, xelas, awr, awi,
                                              izai_z, target_z))
    end

    # 2. Build AngularBlock — incident energies in MeV (xss/emev convention).
    distributions = Vector{Union{TabulatedAngular, Nothing}}(undef, length(results))
    for (i, r) in enumerate(results)
        distributions[i] = TabulatedAngular(Int32(2), copy(r.mu),
                                             copy(r.pdf), copy(r.cdf))
    end
    block_energies_mev = [round_sigfig(r.e / _ACER_EMEV, 7, 0) for r in results]
    angular_block = AngularBlock(block_energies_mev, distributions)

    # 3. Log-log interpolate yys = sigtot_int onto every ESZ point and
    #    apply the Coulomb-corrected elastic + total update. Each result
    #    `r` provides (xxs[i] = r.e, yys[i] = r.sigtot_int).
    xxs = [r.e for r in results]
    yys = [r.sigtot_int for r in results]

    new_elastic = similar(esz_elastic_xs)
    new_total   = similar(esz_total_xs)
    heating     = zeros(Float64, length(esz_e))
    nes = length(esz_e)
    for j in 1:nes
        e = esz_e[j]
        # signi here is the sigfig-7 ESZ elastic column (acefc.f90:6659).
        signi_orig = esz_elastic_sig7[j]
        signow = _log_log_interp(xxs, yys, e)
        signow = round_sigfig(signow, 7, 0)
        new_elastic[j] = signow
        # total[j] - signi_orig + signow, with sigfig(., 9, 0) per acefc.f90:6659.
        new_total[j] = round_sigfig(esz_total_xs[j] - signi_orig + signow, 9, 0)
        # Heating: linearly interpolate the per-E heating contribution.
        # Suppressed for α + α (izai==za).
        if izai != za && signow > 0
            h_at_e = _terpa(results, e)
            heating[j] = round_sigfig(h_at_e / signow, 7, 0)
        end
    end

    angular_block, new_total, new_elastic, heating
end

# ----- helpers ------------------------------------------------------------

@inline function _esz_lin_interp(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    n = length(xs)
    if x <= xs[1]
        return ys[1]
    elseif x >= xs[end]
        return ys[end]
    end
    j = searchsortedlast(xs, x)
    f = (xs[j+1] - x) / (xs[j+1] - xs[j])
    ys[j] * f + ys[j+1] * (1 - f)
end

"""Log-log interpolation reproducing Fortran terp1(x1,y1,x2,y2,x,y,5)
exactly (endf.f90:1627-1632), driven by the same bracket-search the
acecpe post-loop uses (acefc.f90:6653-6657):

    jj=1; do while (xxs(jj+1) < e and jj < ne); jj=jj+1; enddo
    terp1(xxs(jj),yys(jj),xxs(jj+1),yys(jj+1),e,signow,5)

terp1's int=5 form is `y = y1*exp(log(x/x1)*log(y2/y1)/log(x2/x1))`, with
the early exits `x2==x1 → y1`, `y2==y1 → y1`, `x==x1 → y1`, `y1==0 → y1`.
We must match the ratio-then-log order (`log(x/x1)`, not `log(x)-log(x1)`)
for bit-identical FP — this is the form that produces ENDF a11 sigfig
boundaries identical to NJOY."""
@inline function _log_log_interp(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    n = length(xs)
    # Bracket search mirroring acefc.f90:6653-6657 (1-based jj over xxs).
    jj = 1
    while jj < n && xs[jj+1] < x
        jj += 1
    end
    x1 = xs[jj];   y1 = ys[jj]
    x2 = xs[jj+1]; y2 = ys[jj+1]
    # terp1 int=5 (endf.f90:1604-1632).
    x2 == x1 && return y1
    (y2 == y1 || x == x1) && return y1
    y1 == 0.0 && return y1
    y1 * exp(log(x / x1) * log(y2 / y1) / log(x2 / x1))
end

"""Linearly interpolate the heating contribution `eht` on the (x=E, y=eht)
table built by acecpe. Uses NJOY's terpa with INT=2 (lin-lin)."""
@inline function _terpa(results::Vector{AceCpeIncident}, x::Float64)
    n = length(results)
    if x <= results[1].e || n == 1
        return results[1].heating
    elseif x >= results[end].e
        return results[end].heating
    end
    j = 1
    while j < n - 1 && results[j+1].e < x
        j += 1
    end
    f = (results[j+1].e - x) / (results[j+1].e - results[j].e)
    results[j].heating * f + results[j+1].heating * (1 - f)
end
