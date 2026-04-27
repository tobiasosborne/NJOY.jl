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

        # Compute sigc and signi.
        local sigc, signi, ratr
        retry = true
        while retry
            sigc = coulomb_sigc(amuu, e, awr, awi, izai_z, target_z, spi, lidp)
            local pmu_local = pmu
            if ltp < 12
                pmu_local = pmu - sigc
            end
            if iterp == 1
                # Use interpolated ratr to back-compute signi.
                ratr_interp = (ratrl + ((sigc + pmu_local * xelas) / sigc)) / 2
                # The Fortran sets ratr = (ratrl + ratr)/2 BEFORE the inner loop
                # actually computes the new ratr; we mirror by computing the
                # endpoint's tentative ratr first.
                signi = (ratr_interp - 1) * sigc
                ratr  = ratr_interp
            else
                signi = pmu_local * xelas
                if signi < -sigc
                    signi = -sigc
                end
                ratr = (sigc + signi) / sigc
            end
            retry = false
            # Subdivision trigger: if not first point, not already in iterp,
            # and Coulomb dominates AND grew >2× since last point.
            if jl > 1 && iterp == 0 &&
               sigc > abs(signi) && sigc > 2 * sigcl
                iterp = 1
                # Insert midpoint: (amul + amuu)/2 with ratr=(ratrl+ratr)/2.
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

    # Normalize: pdf = sigtot/cumm (per acefc.f90:6634-6636), cdf = cumm/cumm.
    pdf = [sigtot_out[i] / cumm for i in eachindex(sigtot_out)]
    cdf = [cumm_out[i]   / cumm for i in eachindex(cumm_out)]
    AceCpeIncident(e, mu_out, pdf, cdf, cumm * 2π, eht * 2π)
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
    results = AceCpeIncident[]
    for sub in subs
        # Linear-interp xelas from ESZ at sub.e (acefc.f90:6549-6553).
        xelas = _esz_lin_interp(esz_e, esz_elastic_xs, sub.e)
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
        signi_orig = esz_elastic_xs[j]
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

"""Log-log interpolation matching Fortran terp1(...,...,5). For
non-positive y the routine falls back to lin-lin."""
@inline function _log_log_interp(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    n = length(xs)
    if x <= xs[1]
        return ys[1]
    elseif x >= xs[end]
        return ys[end]
    end
    j = searchsortedlast(xs, x)
    if xs[j] <= 0 || xs[j+1] <= 0 || ys[j] <= 0 || ys[j+1] <= 0
        f = (x - xs[j]) / (xs[j+1] - xs[j])
        return ys[j] + f * (ys[j+1] - ys[j])
    end
    lx = log(x); lxj = log(xs[j]); lxj1 = log(xs[j+1])
    lyj = log(ys[j]); lyj1 = log(ys[j+1])
    exp(lyj + (lx - lxj) / (lxj1 - lxj) * (lyj1 - lyj))
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
