# PURR -- Probability table generation for the unresolved resonance region
#
# Monte-Carlo ladder sampling following NJOY2016 purr.f90.
# Statistical sampling (chi-squared, Wigner) cleanly separated from physics.

# ── Chi-squared sampling via quantile table (NJOY2016 `chisq`) ───
const CHI2_QUANTILES = [
    # df=1        df=2         df=3         df=4
    1.31003e-3    0.0508548    0.206832     0.459462
    9.19501e-3    0.156167     0.470719     0.893735
    0.0250905     0.267335     0.691933     1.21753
    0.049254      0.38505      0.901674     1.50872
    0.0820892     0.510131     1.10868      1.78605
    0.124169      0.643564     1.31765      2.05854
    0.176268      0.786543     1.53193      2.33194
    0.239417      0.940541     1.75444      2.61069
    0.314977      1.1074       1.98812      2.89878
    0.404749      1.28947      2.23621      3.20032
    0.511145      1.48981      2.50257      3.51995
    0.637461      1.71249      2.79213      3.86331
    0.788315      1.96314      3.11143      4.23776
    0.970419      2.24984      3.46967      4.65345
    1.194         2.58473      3.88053      5.12533
    1.47573       2.98744      4.36586      5.67712
    1.84547       3.49278      4.96417      6.35044
    2.36522       4.17238      5.75423      7.22996
    3.20371       5.21888      6.94646      8.541
    5.58201       7.99146      10.0048      11.8359
]

"""
    chi2_sample(df, rng) -> Float64

Sample from chi-squared(df) using the NJOY 20-point quantile table.
"""
function chi2_sample(df::Int, rng)
    @assert 1 <= df <= 4
    n = clamp(ceil(Int, 20.0 * rand(rng)), 1, 20)
    return CHI2_QUANTILES[n, df]
end

"""
    wigner_spacing(D_mean, rng) -> Float64

Sample a level spacing from the Wigner distribution:
  s = D * sqrt(4/pi) * sqrt(-ln(U))
"""
function wigner_spacing(D_mean::Float64, rng)
    return D_mean * sqrt(4.0 / Float64(pi)) * sqrt(-log(rand(rng)))
end

# ── Probability Table type ────────────────────────────────────────
"""
    ProbabilityTable

Per-energy probability table with `nbins` bins of (probability, total,
elastic, fission, capture) cross sections.
"""
struct ProbabilityTable
    energies::Vector{Float64}
    nbins::Int
    prob::Matrix{Float64}       # (nbins, n_energies)
    total::Matrix{Float64}
    elastic::Matrix{Float64}
    fission::Matrix{Float64}
    capture::Matrix{Float64}
end

# ── Ladder construction ───────────────────────────────────────────
"""
    generate_ladder(seq, elow, ehigh, rng) -> (er, gnr, gfr, ggr, gxr)

Generate one resonance ladder for spin sequence `seq`, returning
resonance energies and partial width ratios (partial/total).
"""
function generate_ladder(seq::URRSpinSequence, elow::Float64,
                         ehigh::Float64, rng)
    er = Float64[]; gnr = Float64[]; gfr = Float64[]
    ggr = Float64[]; gxr = Float64[]
    dcon = seq.D * sqrt(4.0 / Float64(pi))
    e_res = elow + dcon * rand(rng)
    while e_res <= ehigh
        gn = (seq.GN0 / seq.AMUN) * chi2_sample(seq.AMUN, rng)
        gf = seq.GF
        if seq.AMUF > 0 && seq.GF > 0.0
            gf = (seq.GF / seq.AMUF) * chi2_sample(seq.AMUF, rng)
        end
        gx = seq.GX
        if seq.AMUX > 0 && seq.GX > 0.0
            gx = (seq.GX / seq.AMUX) * chi2_sample(seq.AMUX, rng)
        end
        gt = gn + gf + seq.GG + gx
        push!(er, e_res); push!(gnr, gn/gt); push!(gfr, gf/gt)
        push!(ggr, seq.GG/gt); push!(gxr, gx/gt)
        e_res += wigner_spacing(seq.D, rng)
    end
    return er, gnr, gfr, ggr, gxr
end

# ── Doppler-broadened line shape ──────────────────────────────────
@inline function line_shape(x::Float64, y::Float64)
    yy = y * y
    if abs(x) > 6.0 || y > 6.0
        denom = x * x + yy
        c1 = 0.5641895835
        return y * c1 / denom, x * c1 / denom
    end
    a1 = x * x - yy; a2 = 2.0 * x * y; a3 = a2 * a2
    t1 = a2 * x; t2 = a2 * y
    a4 = a1 - 0.2752551; a5 = a1 - 2.724745
    f1 = 0.5124242 / (a4 * a4 + a3)
    f2 = 0.05176536 / (a5 * a5 + a3)
    rew  = f1 * (t1 - a4 * y) + f2 * (t1 - a5 * y)
    aimw = f1 * (a4 * x + t2) + f2 * (a5 * x + t2)
    return rew, aimw
end

# ── Main probability table generator ─────────────────────────────
"""
    generate_ptable(model, energies; nladders=64, nbins=20, seed=12345,
                    T=300.0, bkg=(0,0,0,0)) -> ProbabilityTable

Generate probability tables by Monte Carlo ladder sampling (PURR algorithm).
"""
function generate_ptable(model::URRStatModel, energies::Vector{Float64};
                         nladders::Int=64, nbins::Int=20, seed::Int=12345,
                         T::Float64=300.0,
                         bkg::NTuple{4,Float64}=(0.0,0.0,0.0,0.0))
    rng = Random.Xoshiro(seed)
    C = NJOY.PhysicsConstants
    ne = length(energies)
    prob_out = zeros(nbins, ne); total_out = zeros(nbins, ne)
    elastic_out = zeros(nbins, ne); fission_out = zeros(nbins, ne)
    capture_out = zeros(nbins, ne)
    cwaven = sqrt(2.0 * C.amassn * C.amu * C.ev) * 1e-12 / C.hbar
    awri = model.AWRI; rat = awri / (awri + 1.0)
    aw = awri * C.amassn; aa = 0.123 * aw^(1.0/3.0) + 0.08
    Teff = max(T, 1.0); nsamp = 1000

    for (ie, E) in enumerate(energies)
        sqE = sqrt(E); k = cwaven * rat * sqE
        ab = 4.0 * Float64(pi) / k^2
        rho = k * aa; rhoc = k * model.AP
        con1 = sqrt(2901.34 * awri / (E * Teff))
        # Potential scattering
        spot = 0.0; seen_l = Set{Int}()
        for seq in model.sequences
            _, ps = urr_penetrability(seq.l, rho, rhoc)
            if !(seq.l in seen_l)
                push!(seen_l, seq.l)
                spot += ab * (2 * seq.l + 1) * sin(ps)^2
            end
        end
        dmin = minimum(s.D for s in model.sequences)
        erange = 900.0 * dmin
        elow = E - erange/2; ehigh = E + erange/2
        navoid = 100.0 / sum(1.0/s.D for s in model.sequences)
        emin = elow + navoid; emax = ehigh - navoid
        if emax <= emin; emax = emin + dmin; end
        # Accumulators
        tp = zeros(nbins); tt = zeros(nbins); te = zeros(nbins)
        tf = zeros(nbins); tc = zeros(nbins)
        total_samp = 0; bin_edges = zeros(nbins)

        for iladr in 1:nladders
            es = sort!([emin + (emax-emin)*rand(rng) for _ in 1:nsamp])
            els = fill(spot, nsamp); fis = zeros(nsamp); cap = zeros(nsamp)
            for seq in model.sequences
                Vl, ps = urr_penetrability(seq.l, rho, rhoc)
                gnx = seq.GN0 * Vl * sqE * Float64(seq.AMUN)
                gj = (2.0*seq.J+1.0) / (4.0*model.SPI+2.0)
                er, gnr_l, gfr_l, ggr_l, _ = generate_ladder(seq, elow, ehigh, rng)
                rpi = sqrt(Float64(pi))
                for (ir, Eres) in enumerate(er)
                    gt_full = gnx / gnr_l[ir]
                    y = con1 * gt_full / 2.0
                    szy = ab * gj * gnr_l[ir] * rpi * y
                    cc2 = szy * (cos(2*ps) - 1.0 + gnr_l[ir])
                    cs2 = szy * sin(2*ps)
                    ccg = szy * ggr_l[ir]; ccf = szy * gfr_l[ir]
                    for isamp in 1:nsamp
                        dx = con1 * (es[isamp] - Eres)
                        rew, aimw = if abs(dx) > 100.0 && y < 100.0
                            d2 = dx*dx + y*y
                            (y * 0.5641895835 / d2, dx * 0.5641895835 / d2)
                        else
                            line_shape(dx, y)
                        end
                        cap[isamp] += ccg * rew
                        fis[isamp] += ccf * rew
                        els[isamp] += cc2 * rew + cs2 * aimw
                    end
                end
            end
            for i in 1:nsamp
                if els[i] < -bkg[2]; els[i] = -bkg[2] + 1e-6; end
            end
            tot = [els[i]+fis[i]+cap[i]+bkg[1] for i in 1:nsamp]
            if iladr == 1
                st = sort(tot)
                for ib in 1:nbins-1
                    idx = clamp(round(Int, ib/nbins*nsamp), 1, nsamp)
                    bin_edges[ib] = st[idx]
                    if ib > 1 && bin_edges[ib] <= bin_edges[ib-1]
                        bin_edges[ib] = bin_edges[ib-1]*1.05 + 1e-30
                    end
                end
                bin_edges[nbins] = 1e6
            end
            for i in 1:nsamp
                ib = clamp(searchsortedfirst(view(bin_edges,1:nbins), tot[i]), 1, nbins)
                total_samp += 1; tp[ib] += 1.0
                tt[ib] += tot[i]; te[ib] += els[i]+bkg[2]
                tf[ib] += fis[i]+bkg[3]; tc[ib] += cap[i]+bkg[4]
            end
        end
        for ib in 1:nbins
            cnt = tp[ib]
            if cnt > 0
                total_out[ib,ie] = tt[ib]/cnt; elastic_out[ib,ie] = te[ib]/cnt
                fission_out[ib,ie] = tf[ib]/cnt; capture_out[ib,ie] = tc[ib]/cnt
            end
            prob_out[ib,ie] = cnt / total_samp
        end
    end
    return ProbabilityTable(energies, nbins, prob_out, total_out,
                            elastic_out, fission_out, capture_out)
end

# ── Bondarenko moments from probability table ─────────────────────
"""
    bondarenko_from_ptable(ptable, ie, sigma0) -> NTuple{5,Float64}

Compute Bondarenko self-shielded cross sections from probability table
at energy index `ie` for dilution `sigma0`.
Returns (total, elastic, fission, capture, transport).
"""
function bondarenko_from_ptable(ptable::ProbabilityTable, ie::Int, sigma0::Float64)
    bval = zeros(7)
    for j in 1:ptable.nbins
        p = ptable.prob[j, ie]
        if p <= 0.0; continue; end
        t = ptable.total[j, ie]
        den = sigma0 / (sigma0 + t)
        bval[1] += p * t * den
        bval[2] += p * ptable.elastic[j, ie] * den
        bval[3] += p * ptable.fission[j, ie] * den
        bval[4] += p * ptable.capture[j, ie] * den
        bval[5] += p * t * den * den
        bval[6] += p * den
        bval[7] += p * den * den
    end
    return (bval[1]/bval[6], bval[2]/bval[6], bval[3]/bval[6],
            bval[4]/bval[6], bval[5]/bval[7])
end
