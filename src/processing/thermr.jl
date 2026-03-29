# THERMR -- Thermal scattering cross sections and kernels
# PROPOSAL B: Composable, AD-compatible pure functions.
# Free-gas and S(alpha,beta) share a generic kernel interface (E,E',mu).
# Detailed balance is exactly satisfied by construction.
# Correspondence: free_gas_kernel->sig()@200, sab_kernel->sig()@2514,
#   bragg_edges->sigcoh(), incoh_elastic_xs->iel(), compute_thermal->thermr()

# ---- Data types ----

"Tabulated S(alpha,beta) thermal scattering law (ENDF MF7/MT4). Stores log(S)."
struct SABData
    alpha::Vector{Float64}; beta::Vector{Float64}
    sab::Matrix{Float64}   # log S(alpha,beta), (n_alpha, n_beta)
    sigma_b::Float64; awr::Float64; T_eff::Float64
    lasym::Int; lat::Int
end

"Precomputed Bragg edge data for coherent elastic scattering."
struct BraggData
    tau_sq::Vector{Float64}; form_factor::Vector{Float64}
    econ::Float64; scon::Float64; n_edges::Int
end

"Output of compute_thermal: thermal cross sections on an energy grid."
struct ThermalResult
    energies::Vector{Float64}; inelastic_xs::Vector{Float64}
    elastic_xs::Vector{Float64}; model::Symbol
end

# ---- Effective temperature lookup (Fortran gateff, thermr.f90:605-701) ----

# Hardcoded (mat, T, T_eff) table matching Fortran gateff
const _GATEFF_TABLE = [
    (1002, 296.0, 1396.8), (1002, 350.0, 1411.6), (1002, 400.0, 1427.4),
    (1002, 450.0, 1444.9), (1002, 500.0, 1464.1), (1002, 600.0, 1506.8),
    (1002, 800.0, 1605.8), (1002, 1000.0, 1719.8),
    (1004, 296.0, 940.91), (1004, 350.0, 961.62), (1004, 400.0, 982.93),
    (1004, 450.0, 1006.1), (1004, 500.0, 1030.9), (1004, 600.0, 1085.1),
    (1004, 800.0, 1209.0), (1004, 1000.0, 1350.0),
    (1064, 296.0, 405.64), (1064, 400.0, 484.22), (1064, 500.0, 568.53),
    (1064, 600.0, 657.66), (1064, 700.0, 749.69), (1064, 800.0, 843.63),
    (1064, 1000.0, 1035.0), (1064, 1200.0, 1229.3),
    (1065, 296.0, 713.39), (1065, 400.0, 754.68), (1065, 500.0, 806.67),
    (1065, 600.0, 868.38), (1065, 700.0, 937.64), (1065, 800.0, 1012.7),
    (1065, 1000.0, 1174.9), (1065, 1200.0, 1348.2), (1065, 1600.0, 1712.9),
    (1065, 2000.0, 2091.0),
    (1095, 296.0, 1165.9), (1095, 350.0, 1177.8), (1095, 400.0, 1191.4),
    (1095, 450.0, 1207.7), (1095, 500.0, 1226.0), (1095, 600.0, 1268.7),
    (1095, 800.0, 1373.4), (1095, 1000.0, 1497.7),
    (1096, 296.0, 317.27), (1096, 400.0, 416.29), (1096, 500.0, 513.22),
    (1096, 600.0, 611.12), (1096, 700.0, 709.60), (1096, 800.0, 808.43),
    (1096, 1000.0, 1006.8), (1096, 1200.0, 1205.7),
    (1097, 296.0, 806.79), (1097, 400.0, 829.98), (1097, 500.0, 868.44),
    (1097, 600.0, 920.08), (1097, 700.0, 981.82), (1097, 800.0, 1051.1),
    (1097, 1000.0, 1205.4), (1097, 1200.0, 1373.4),
    (1099, 296.0, 596.4), (1099, 400.0, 643.9), (1099, 500.0, 704.6),
    (1099, 600.0, 775.3), (1099, 800.0, 935.4), (1099, 1000.0, 1109.8),
    (1099, 1200.0, 1292.3),
    (1114, 296.0, 1222.0), (1114, 350.0, 1239.0),
]

"Look up T_eff from Fortran gateff hardcoded table. Returns T if not found."
function _gateff_lookup(mat::Integer, T::Real)
    for (m, t, te) in _GATEFF_TABLE
        m == mat && abs(t - T) <= 5.0 && return te
    end
    return Float64(T)  # default: T_eff = T (Fortran line 697)
end

# ---- MF7/MT4 reader ----

"""
    read_mf7_mt4(filename, mat, T) -> SABData

Read S(α,β) thermal scattering law from ENDF MF7/MT4.
Extracts data at temperature `T` [K] for material `mat`.
Returns SABData with log(S) values.

Matches Fortran thermr calcem (thermr.f90:1659-1779).
"""
function read_mf7_mt4(filename::AbstractString, mat::Integer, T::Real)
    open(filename) do io
        find_section(io, 7, 4; target_mat=mat) ||
            error("MF7/MT4 not found for MAT=$mat in $filename")
        head = read_cont(io)
        lasym = head.N1
        lat = head.L2

        # B(n) constants (LIST record)
        blist = read_list(io)
        sigma_free = blist.data[1]  # free atom XS [barn]
        dw_integral = blist.data[2]  # Debye-Waller integral
        az = blist.data[3]           # atomic mass of principal scatterer
        natom = max(1, round(Int, blist.data[4]))
        smz = sigma_free / natom
        sigma_b = smz * ((az + 1) / az)^2  # bound atom XS

        # TAB2: beta grid dimensions
        tab2 = read_tab2(io)
        nbeta = tab2.NZ

        # Read alpha grid and S(α,β) at target temperature
        alpha_grid = Float64[]
        beta_grid = Float64[]
        sab_matrix = Matrix{Float64}(undef, 0, 0)
        sabflg = -225.0
        first_beta = true
        ttol = T / 500.0  # temperature tolerance (Fortran line 1720)

        for ib in 1:nbeta
            # TAB1: temperature, beta, LT, 0, NR, NA with (alpha, S) pairs
            tab1 = read_tab1(io)
            t1 = tab1.C1
            beta_val = tab1.C2
            lt = tab1.L1  # number of additional temperatures
            na = length(tab1.x)  # number of alpha values

            push!(beta_grid, beta_val)

            if first_beta
                alpha_grid = copy(tab1.x)
                sab_matrix = Matrix{Float64}(undef, na, nbeta)
                first_beta = false
            end

            # Store S values at first temperature (tab1.y)
            if abs(t1 - T) <= ttol
                for ia in 1:na
                    v = tab1.y[ia]
                    sab_matrix[ia, ib] = v > exp(sabflg) ? log(v) : sabflg
                end
            else
                # Wrong temperature — check additional temps
                found_t = false
                for it in 1:lt
                    tlist = read_list(io)
                    if !found_t && abs(tlist.C1 - T) <= ttol
                        found_t = true
                        for ia in 1:na
                            v = tlist.data[ia]
                            sab_matrix[ia, ib] = v > exp(sabflg) ? log(v) : sabflg
                        end
                    end
                end
                found_t || error("Temperature $T not found in MF7/MT4 (available: $t1 + $lt more)")
                continue  # already consumed all LT records
            end

            # Skip remaining LT additional temperature records
            for it in 1:lt
                read_list(io)
            end
        end

        # Effective temperature for SCT fallback (Fortran thermr.f90 gateff, lines 615-697)
        # ENDF6+ files may have a TAB1 after S(α,β) data; older files use a hardcoded table.
        teff_K = _gateff_lookup(mat, T)

        SABData(alpha_grid, beta_grid, sab_matrix,
                smz, az, teff_K,
                lasym, lat)
    end
end

# ---- Free-gas model (sig() lines 2600-2608) ----

"""
    free_gas_kernel(E, E_prime, mu, A, T; sigma_b=1.0) -> Float64

Free-gas d^2sigma/(dE' dmu) [barn/eV]. S(a,b)=exp(-(a+b)^2/(4a))/sqrt(4*pi*a)
exactly satisfies detailed balance. Matches NJOY2016 sig()@200.
"""
function free_gas_kernel(E::Real, E_prime::Real, mu::Real, A::Real, T::Real;
                         sigma_b::Real=1.0)
    kT = PhysicsConstants.bk * T
    rtev = 1.0 / kT  # Fortran precomputes rtev=1/tev (line 2512)
    Ef, Epf = Float64(E), Float64(E_prime)
    (Ef <= 0 || Epf <= 0 || kT <= 0) && return 0.0
    alpha = max((Ef + Epf - 2.0*Float64(mu)*sqrt(Ef*Epf)) / (A*kT), 1e-6)  # Fortran amin=1e-6 (line 2500)
    bb = (Epf - Ef) * rtev  # Fortran bb=(ep-e)*rtev (line 2513) — multiplication, NOT division
    arg = (alpha + bb)^2 / (4.0*alpha)
    arg > 225.0 && return 0.0
    sigc = sqrt(Epf/Ef) * rtev * 0.5  # Fortran sigc=sqrt(ep/e)*rtev/2 (line 2516)
    sig = sigc * sigma_b * exp(-arg) / sqrt(4.0*PhysicsConstants.pi*alpha)
    sig < 1e-10 ? 0.0 : sig  # Fortran sigmin=1e-10 (line 2606)
end

"""
    free_gas_xs(E, A, T; sigma_b=1.0) -> Float64

Analytical free-gas total thermal XS [barn].
sigma = (sigma_b/A) * [(1+1/(2x^2))*erf(x) + exp(-x^2)/(x*sqrt(pi))],
where x = sqrt(A*E/kT).
Limits: E >> kT -> sigma_b/A; E << kT -> sigma_b*2/(A*x*sqrt(pi)) (1/v).
"""
function free_gas_xs(E::Real, A::Real, T::Real; sigma_b::Real=1.0)
    kT = PhysicsConstants.bk * T; x2 = A * Float64(E) / kT; x = sqrt(x2)
    x < 1.0e-6 && return sigma_b / A * 2.0 / (x * sqrt(PhysicsConstants.pi))
    isp = 1.0 / sqrt(PhysicsConstants.pi); ef = erf(x); ex = exp(-x2)
    sigma_b / A * ((1.0 + 1.0 / (2.0 * x2)) * ef + ex * isp / x)
end

# ---- S(alpha,beta) model (sig() lines 2514-2598) ----

"""
    read_thermal_data(alpha, beta, sab_values; sigma_b, awr, ...) -> SABData

Construct SABData. If `is_log=false`, raw S values are log-converted.
"""
function read_thermal_data(alpha::AbstractVector, beta::AbstractVector,
                           sab_values::AbstractMatrix;
                           sigma_b::Real, awr::Real, T_eff::Real=0.0,
                           lasym::Int=0, lat::Int=0, is_log::Bool=true)
    na, nb = length(alpha), length(beta)
    @assert size(sab_values) == (na, nb)
    sabflg = -225.0
    sl = is_log ? Float64.(sab_values) : [Float64(v) > exp(sabflg) ?
         log(Float64(v)) : sabflg for v in sab_values]
    sl = reshape(sl, na, nb)
    SABData(Float64.(alpha), Float64.(beta), sl,
            Float64(sigma_b), Float64(awr), Float64(T_eff), lasym, lat)
end

"Quadratic interpolation matching Fortran terpq (thermr.f90:2617-2658)."
function _terpq(x1::Real, y1::Real, x2::Real, y2::Real, x3::Real, y3::Real, x::Real)
    sabflg = -225.0; step = 2.0
    if x < x1
        # Below range: log-lin (INT=3) if increasing, else clamp
        y1 > y2 && return y1
        x1 > 0 && x > 0 || return y1
        return y1 + (y2 - y1) * (log(x) - log(x1)) / (log(x2) - log(x1))
    elseif x > x3
        # Above range: lin-lin if decreasing, else clamp
        y3 > y2 && return y3
        return y2 + (y3 - y2) * (x - x2) / (x3 - x2)
    elseif abs(y1 - y2) > step || abs(y2 - y3) > step
        # Large steps: piecewise linear
        if x < x2
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        else
            return y2 + (y3 - y2) * (x - x2) / (x3 - x2)
        end
    else
        # Quadratic (Newton forward difference)
        bq = (y2 - y1) * (x3 - x1) / ((x2 - x1) * (x3 - x2)) -
             (y3 - y1) * (x2 - x1) / ((x3 - x1) * (x3 - x2))
        cq = (y3 - y1) / ((x3 - x1) * (x3 - x2)) -
             (y2 - y1) / ((x2 - x1) * (x3 - x2))
        return y1 + bq * (x - x1) + cq * (x - x1)^2
    end
end

"""Biquadratic S(α,β) interpolation matching Fortran sig (thermr.f90:2514-2566).
Uses terpq (3-point quadratic) in alpha then beta, with cliq analytical
extrapolation for α below grid (small β), and large-step linear fallback."""
function _interp_sab(a::Real, b::Real, data::SABData)
    sabflg = -225.0; test2 = 30.0
    al, be, sab = data.alpha, data.beta, data.sab
    na, nb = length(al), length(be)
    (a <= 0 || a > al[end]) && return sabflg
    bv = data.lasym == 0 ? abs(b) : b
    (bv > be[end] || bv < be[1]) && return sabflg

    # Find grid indices (Fortran loop: ia such that a < alpha(ia+1))
    # Snap to next grid point when within FP rounding tolerance — the Fortran's
    # b*tev/tevz roundtrip can put b epsilon above a grid point while Julia's
    # equivalent puts it epsilon below, causing a one-index shift that changes
    # which sabflg corner cells are checked.
    ia = clamp(searchsortedlast(al, a), 1, na - 1)
    if ia < na - 1 && abs(a - al[ia + 1]) < 1e-10 * al[ia + 1]
        ia += 1
    end
    ib = clamp(searchsortedlast(be, bv), 1, nb - 1)
    if ib < nb - 1 && abs(bv - be[ib + 1]) < 1e-10 * be[ib + 1]
        ib += 1
    end
    bbb = data.lasym == 1 && b < 0 ? -bv : bv

    # Sabflg corner check (Fortran lines 2546-2550) — skip if small α·A and small β
    if !(a * data.awr < test2 && bv < test2)
        (sab[ia, ib] <= sabflg || sab[min(ia+1,na), ib] <= sabflg ||
         sab[ia, min(ib+1,nb)] <= sabflg || sab[min(ia+1,na), min(ib+1,nb)] <= sabflg) &&
            return sabflg
    end

    # Ensure 3 points available (Fortran lines 2552-2553)
    ia + 2 > na && (ia = max(1, na - 2))
    ib + 2 > nb && (ib = max(1, nb - 2))

    # Biquadratic: terpq in alpha for 3 beta values, then terpq in beta
    s1 = _terpq(al[ia], sab[ia, ib],   al[ia+1], sab[ia+1, ib],   al[ia+2], sab[ia+2, ib],   a)
    s2 = _terpq(al[ia], sab[ia, ib+1], al[ia+1], sab[ia+1, ib+1], al[ia+2], sab[ia+2, ib+1], a)
    s3 = _terpq(al[ia], sab[ia, ib+2], al[ia+1], sab[ia+1, ib+2], al[ia+2], sab[ia+2, ib+2], a)
    s = _terpq(be[ib], s1, be[ib+1], s2, be[ib+2], s3, bbb)
    max(s, sabflg)
end

"""
    sab_kernel(E, E_prime, mu, data::SABData, T) -> Float64

Differential kernel from tabulated S(a,b) [barn/eV]. Falls back to SCT
approximation outside table range. Detailed balance is exact.
"""
function sab_kernel(E::Real, E_prime::Real, mu::Real, data::SABData, T::Real;
                    sigma_b_override::Union{Real,Nothing}=nothing)
    kT = PhysicsConstants.bk * T; sabflg = -225.0
    Ef, Epf = Float64(E), Float64(E_prime)
    (Ef <= 0 || Epf <= 0 || kT <= 0) && return 0.0
    beta_raw = (Epf - Ef) / kT
    alpha_raw = max((Ef + Epf - 2.0*Float64(mu)*sqrt(Ef*Epf))/(data.awr*kT), 1e-6)
    tevz = 0.0253
    bt = data.lat==1 ? abs(beta_raw)*kT/tevz : abs(beta_raw)
    at = data.lat==1 ? alpha_raw*kT/tevz : alpha_raw
    sb = something(sigma_b_override, data.sigma_b) * ((data.awr+1)/data.awr)^2
    pf = sqrt(Epf/Ef) / (2.0*kT) * sb
    ls = _interp_sab(at, bt, data)
    sigmin = 1e-10  # Fortran sig line 2498/2565/2596
    if ls > sabflg
        sig = pf * exp(ls - beta_raw/2)
        return sig < sigmin ? 0.0 : sig
    end
    # SCT fallback
    Te = data.T_eff > 0 ? data.T_eff*PhysicsConstants.bk : kT
    bs, as = abs(beta_raw), alpha_raw
    arg = (as-bs)^2*kT/(4*as*Te) + (bs+beta_raw)/2
    sig = arg < 225.0 ? pf*exp(-arg)/sqrt(4*PhysicsConstants.pi*as*Te/kT) : 0.0
    sig < sigmin ? 0.0 : sig
end

"Integrated incoherent inelastic XS from S(a,b) [barn]."
function sab_xs(E::Real, data::SABData, T::Real; n_mu::Int=20, n_ep::Int=40)
    E <= 0 && return 0.0; kT = PhysicsConstants.bk*T
    Epm = max(E + data.beta[end]*(data.lat==1 ? 0.0253 : kT), 10*kT)
    dEp, dmu, total = Epm/n_ep, 2.0/n_mu, 0.0
    for iep in 1:n_ep; Ep=(iep-0.5)*dEp; Ep<=0 && continue
        ms = 0.0
        for imu in 0:n_mu; mu=-1.0+imu*dmu
            w = (imu==0||imu==n_mu) ? 0.5 : 1.0
            ms += w*sab_kernel(E, Ep, mu, data, T)
        end; total += ms*dmu*dEp
    end; total
end

# ---- Coherent elastic / Bragg edges (sigcoh(), thermr.f90:947-1207) ----

"""
    structure_factor(lat, l1, l2, l3) -> Float64

Crystal structure factor for lattice type `lat`, matching Fortran `form()`.
- lat=1: graphite (hexagonal with basis)
- lat=2: beryllium (HCP)
- lat=3: beryllium oxide
"""
function structure_factor(lat::Int, l1::Int, l2::Int, l3::Int)
    pi_v = PhysicsConstants.pi
    if lat == 1
        # graphite
        if isodd(l3)
            return sin(pi_v * (l1 - l2) / 3)^2
        else
            return (6 + 10 * cos(2 * pi_v * (l1 - l2) / 3)) / 4
        end
    elseif lat == 2
        # beryllium
        return 1 + cos(2 * pi_v * (2 * l1 + 4 * l2 + 3 * l3) / 6)
    elseif lat == 3
        # beryllium oxide
        beo1, beo2, beo3 = 7.54, 4.24, 11.31
        return (1 + cos(2 * pi_v * (2 * l1 + 4 * l2 + 3 * l3) / 6)) *
               (beo1 + beo2 + beo3 * cos(3 * pi_v * l3 / 4))
    else
        error("structure_factor: unsupported lattice type lat=$lat")
    end
end

"Coherent elastic XS [barn]: sigma(E) = scon/E * sum_{edges below E} f_i."
function bragg_edges(E::Real, bragg::BraggData)
    E <= 0 && return 0.0; tf = 0.0
    for i in 1:bragg.n_edges
        bragg.tau_sq[i] >= E*bragg.econ && break; tf += bragg.form_factor[i]
    end; bragg.scon*tf/E
end

"Sorted Bragg edge energies [eV]."
bragg_edge_energies(b::BraggData) = [b.tau_sq[i]/b.econ for i in 1:b.n_edges]

"""
    build_bragg_data(; a, c, sigma_coh, A_mass, natom=1, debye_waller, emax=5.0, lat=1)

Build Bragg data for lattice type `lat`. a,c in [cm], sigma_coh [barn].
Lattice types: 1=graphite, 2=beryllium, 3=beryllium oxide.
The crystal structure factor (Fortran `form()`) determines which reciprocal
lattice vectors actually contribute to Bragg scattering.
"""
function build_bragg_data(; a::Real, c::Real, sigma_coh::Real,
                           A_mass::Real, natom::Int=1,
                           debye_waller::Real, emax::Real=5.0,
                           lat::Int=1)
    pi_v = PhysicsConstants.pi
    amne = PhysicsConstants.amassn * PhysicsConstants.amu
    econ = PhysicsConstants.ev * 8 * amne / PhysicsConstants.hbar^2
    twopis = (2pi_v)^2
    scon = sigma_coh/natom*(4pi_v)^2 / (2*a^2*c*sqrt(3.0)*econ)
    t2 = PhysicsConstants.hbar / (2*PhysicsConstants.amu*A_mass)
    wint = 0.658173e-15 * A_mass * debye_waller
    c1, c2, ulim = 4.0/(3*a^2), 1.0/c^2, econ*emax
    tsl, ffl = Float64[], Float64[]
    tsqx = econ / 20  # Fortran line 1014: merge threshold
    phi = ulim/twopis; i1m = floor(Int, a*sqrt(phi))+1
    for i1 in 1:i1m; l1=i1-1
        disc = 3*(a^2*phi - l1^2); disc < 0 && continue
        i2m = floor(Int, 0.5*(l1+sqrt(disc)))+1
        for i2 in i1:i2m; l2=i2-1
            xv = phi - c1*(l1^2+l2^2-l1*l2)
            i3m = xv > 0 ? floor(Int, c*sqrt(xv))+1 : 1
            for i3 in 1:i3m; l3=i3-1
                w1 = l1==l2 ? 1.0 : 2.0
                w2 = (l1==0&&l2==0) ? 0.5 : (l1==0||l2==0) ? 1.0 : 2.0
                w3 = l3==0 ? 1.0 : 2.0
                for sgn in (1,-1); l2s=sgn*l2
                    tsq = (c1*(l1^2+l2s^2+l1*l2s)+l3^2*c2)*twopis
                    (tsq<=0||tsq>ulim) && continue
                    f_sf = structure_factor(lat, l1, l2s, l3)
                    wt = exp(-tsq*t2*wint)*w1*w2*w3*f_sf/sqrt(tsq)
                    # Fortran lines 1079,1101: only merge when tsq > tsqx
                    # Fortran merge is ONE-SIDED: tsl[k] <= tsq < 1.05*tsl[k]
                    found = false
                    if tsq > tsqx
                        for k in eachindex(tsl)
                            if tsq >= tsl[k] && tsq < 1.05*tsl[k]
                                ffl[k]+=wt; found=true; break
                            end
                        end
                    end
                    found || (push!(tsl,tsq); push!(ffl,wt))
    end; end; end; end
    p = sortperm(tsl)
    tsl_sorted = tsl[p]; ffl_sorted = ffl[p]
    # Fortran sigcoh lines 1134-1137: add sentinel point at ulim
    # with form factor copied from the last sorted edge
    if !isempty(ffl_sorted)
        push!(tsl_sorted, ulim)
        push!(ffl_sorted, ffl_sorted[end])
    end
    BraggData(tsl_sorted, ffl_sorted, econ, scon, length(tsl_sorted))
end

# ---- Thermal output grid (matching Fortran coh adaptive merge) ----

"""
    build_thermal_grid(bragg, input_grid, emax; tol=0.05)

Build merged thermal output grid matching Fortran thermr coh (lines 790-883).

The Fortran coh uses a LIFO convergence stack (depth 20) that processes one Bragg
edge at a time, dynamically interleaving input grid points. The input_grid should
be the elastic PENDF grid (points from iold, i.e., the broadened grid below emax).

Stack convention (matching Fortran upstk lines 924-945):
  stk[1] = highest energy (latest Bragg edge)
  stk[i] = lowest energy (oldest unprocessed point)
  upstk inserts new point at i-1, shifts old i-1 to i
  Convergence tests interval [stk[i], stk[i-1]] (lowest to next)
  Output pops stk[i] (lowest energy first)
"""
function build_thermal_grid(bragg::BraggData, input_grid::AbstractVector{Float64},
                            emax::Float64; tol::Float64=0.05)
    eps = 3.0e-5   # Fortran coh line 763: small spacing threshold
    tolmin = 1e-6   # Fortran coh line 765
    imax = 20       # Fortran coh line 774: max stack depth

    # Sorted input grid below emax (elastic PENDF from iold)
    inp = sort(filter(e -> e > 0 && e <= emax, input_grid))

    # Bragg edge energies: Fortran sigcoh (lines 1199-1200) returns each edge TWICE:
    # First call returns sigfig(elim, 7, -1) (nudged down).
    # Second call (at the nudged-down energy) returns sigfig(elim, 7, +1) (nudged up),
    # because e > sigfig(sigfig(elim,7,-1), 7, -1) triggers the +1 path.
    # So each physical Bragg edge produces TWO grid points: {nudge-down, nudge-up}.
    bragg_e = Float64[]
    for k in 1:bragg.n_edges
        e_raw = bragg.tau_sq[k] / bragg.econ
        e_raw > emax && continue
        push!(bragg_e, round_sigfig(e_raw, 7, -1))
        push!(bragg_e, round_sigfig(e_raw, 7, +1))
    end
    sort!(unique!(bragg_e))
    isempty(bragg_e) && return sort(unique(vcat(inp, [emax])))

    result = Float64[]
    iex = 1  # index into inp (next unread input point)
    nlt = 5  # sliding window size (Fortran coh line 781)
    nlt_cur = nlt  # current window size (decreases near end)

    # Sliding window of input grid energies (Fortran x[1..nlt] array)
    xwin = zeros(nlt)

    # Phase 1 (Fortran coh lines 794-809): output all input points below first Bragg edge
    # and prime the 5-point sliding window
    first_bragg = bragg_e[1]
    ix = 0  # window fill count
    while iex <= length(inp) && inp[iex] <= first_bragg * (1 + eps)
        xwin[1] = inp[iex]
        ix = 1
        push!(result, inp[iex])
        iex += 1
    end
    # Fill remaining window slots (Fortran lines 804-809: label 105 loop)
    while ix < nlt && iex <= length(inp)
        ix += 1
        xwin[ix] = inp[iex]
        iex += 1
    end

    # Helper: advance sliding window past output energy (Fortran label 170, lines 852-866)
    function advance_window!(out_e)
        while nlt_cur >= 3 && out_e > xwin[3] * (1 + eps) && iex <= length(inp)
            for k in 1:nlt_cur-1
                xwin[k] = xwin[k+1]
            end
            if iex <= length(inp)
                xwin[nlt_cur] = inp[iex]
                iex += 1
            end
            if iex > length(inp)
                nlt_cur = max(nlt_cur - 1, 1)
            end
        end
    end

    # Phase 2 (Fortran coh lines 810-817): prime stack with first Bragg edge
    stk_e = zeros(imax)
    stk_xs = zeros(imax)
    stk_e[1] = first_bragg
    stk_xs[1] = bragg_edges(first_bragg, bragg)
    i = 1
    bidx = 1  # index into bragg_e

    # Phase 3 (Fortran coh lines 819-883): main loop
    while true
        # Label 120: get next Bragg edge
        bidx += 1
        if bidx > length(bragg_e)
            break
        end
        e_next = bragg_e[bidx]
        xs_next = bragg_edges(e_next, bragg)
        # upstk: push new Bragg edge
        i += 1
        stk_e[i] = stk_e[i-1]
        stk_xs[i] = stk_xs[i-1]
        stk_e[i-1] = e_next
        stk_xs[i-1] = xs_next

        # Label 125: include input grid points using sliding window (Fortran lines 825-833)
        while true
            # Scan window x[1..nlt_cur] for first point above stk[i]
            found_input = false
            for wx in 1:nlt_cur
                if xwin[wx] > stk_e[i] * (1 + eps)
                    # Found a window point above stack bottom
                    if xwin[wx] < stk_e[i-1] * (1 - eps)
                        # It's between stk[i] and stk[i-1] → insert via upstk
                        e_in = xwin[wx]
                        xs_in = bragg_edges(e_in, bragg)
                        i += 1
                        stk_e[i] = stk_e[i-1]
                        stk_xs[i] = stk_xs[i-1]
                        stk_e[i-1] = e_in
                        stk_xs[i-1] = xs_in
                        found_input = true
                    end
                    break  # Fortran: first match only (do 130 with go to 135)
                end
            end
            found_input && continue

            # Label 140: convergence test
            if i >= imax
                # Stack full → accept top (lowest energy)
                push!(result, stk_e[i]); advance_window!(stk_e[i])
                i -= 1
                if i > 1; continue; else break; end
            end

            xm = round_sigfig(0.5 * (stk_e[i-1] + stk_e[i]), 7, 0)
            if stk_e[i-1] - stk_e[i] < eps * xm || xm <= stk_e[i] || xm >= stk_e[i-1]
                # Spacing too small → accept
                push!(result, stk_e[i]); advance_window!(stk_e[i])
                i -= 1
                if i > 1; continue; else break; end
            end

            xsm = bragg_edges(xm, bragg)
            # Linear interpolation between stack points (Fortran terp1 with INT=2)
            ym = stk_xs[i] + (xm - stk_e[i]) * (stk_xs[i-1] - stk_xs[i]) /
                 (stk_e[i-1] - stk_e[i])
            test = max(tol * abs(xsm), tolmin)

            if abs(xsm - ym) <= test
                # Converged → accept top (label 160)
                push!(result, stk_e[i]); advance_window!(stk_e[i])
                i -= 1
                if i > 1; continue; else break; end
            else
                # Subdivide → push midpoint (label 210)
                i += 1
                stk_e[i] = stk_e[i-1]
                stk_xs[i] = stk_xs[i-1]
                stk_e[i-1] = xm
                stk_xs[i-1] = xsm
                continue  # go back to label 125
            end
        end
        # i == 1: stack has one point, go get next Bragg edge (label 120)
    end

    # Output remaining stack point
    if i >= 1
        push!(result, stk_e[i])
    end

    sort!(unique!(result))
    return result
end

# Backward-compatible 4-arg signature (calcem_grid ignored, matching Fortran coh)
function build_thermal_grid(bragg::BraggData, input_grid::AbstractVector{Float64},
                            calcem_grid::AbstractVector{Float64}, emax::Float64;
                            tol::Float64=0.05)
    build_thermal_grid(bragg, input_grid, emax; tol=tol)
end

# ---- Incoherent elastic (iel(), thermr.f90:1244-1425) ----

"Incoherent elastic XS: sigma(E) = sigma_b/2 * (1 - exp(-4*E*W'))."
function incoh_elastic_xs(E::Real, sigma_b::Real, dwp::Real)
    E <= 0 && return 0.0; sigma_b/2*(1.0-exp(-4.0*E*dwp))
end

# ---- Standard thermal energy grid (calcem, thermr.f90:1587-1607) ----

const THERMR_EGRID = Float64[
    1e-5,1.78e-5,2.5e-5,3.5e-5,5e-5,7e-5,1e-4,1.26e-4,1.6e-4,2e-4,
    2.53e-4,2.97e-4,3.5e-4,4.2e-4,5.06e-4,6.15e-4,7.5e-4,8.7e-4,
    1.012e-3,1.23e-3,1.5e-3,1.8e-3,2.03e-3,2.277e-3,2.6e-3,3e-3,
    3.5e-3,4.048e-3,4.5e-3,5e-3,5.6e-3,6.325e-3,7.2e-3,8.1e-3,
    9.108e-3,1e-2,1.063e-2,1.15e-2,1.2397e-2,1.33e-2,1.417e-2,1.5e-2,
    1.6192e-2,1.82e-2,1.99e-2,2.0493e-2,2.15e-2,2.28e-2,2.53e-2,
    2.8e-2,3.0613e-2,3.38e-2,3.65e-2,3.95e-2,4.2757e-2,4.65e-2,5e-2,
    5.6925e-2,6.25e-2,6.9e-2,7.5e-2,8.1972e-2,9e-2,9.6e-2,1.035e-1,
    1.11573e-1,1.2e-1,1.28e-1,1.355e-1,1.45728e-1,1.6e-1,1.72e-1,
    1.84437e-1,2e-1,2.277e-1,2.510392e-1,2.705304e-1,2.907501e-1,
    3.011332e-1,3.206421e-1,3.576813e-1,3.9e-1,4.170351e-1,4.5e-1,
    5.032575e-1,5.6e-1,6.25e-1,7e-1,7.8e-1,8.6e-1,9.5e-1,1.05,1.16,
    1.28,1.42,1.55,1.7,1.855,2.02,2.18,2.36,2.59,2.855,3.12,3.42,
    3.75,4.07,4.46,4.9,5.35,5.85,6.4,7.0,7.65,8.4,9.15,9.85,10.0]

# ---- Driver ----

"""
    compute_thermal(pendf, T, A; model=:free_gas, emax=10.0, sigma_b=nothing,
                    sab_data=nothing, bragg=nothing, debye_waller_prime=nothing)

Compute thermal XS and update PENDF material. Replaces elastic channel below emax.
"""
function compute_thermal(pendf::PointwiseMaterial, T::Real, A::Real;
                         model::Symbol=:free_gas, emax::Real=10.0,
                         sigma_b::Union{Real,Nothing}=nothing,
                         sab_data::Union{SABData,Nothing}=nothing,
                         bragg::Union{BraggData,Nothing}=nothing,
                         debye_waller_prime::Union{Real,Nothing}=nothing)
    el_col = findfirst(==(2), pendf.mt_list)
    el_col !== nothing || error("No elastic (MT2) column in PointwiseMaterial")
    tot_col = findfirst(==(1), pendf.mt_list)
    new_e, new_el = compute_thermal_xs(
        pendf.energies, pendf.cross_sections[:,el_col], A, T;
        sigma_b=sigma_b, emax=Float64(emax), model=model, sab_data=sab_data)
    n_new, n_reac = length(new_e), length(pendf.mt_list)
    new_xs = zeros(n_new, n_reac); new_xs[:,el_col] .= new_el
    for j in 1:n_reac; j==el_col && continue
        for i in 1:n_new
            new_xs[i,j] = _linterp(pendf.energies, pendf.cross_sections[:,j], new_e[i])
        end
    end
    if tot_col !== nothing
        for i in 1:n_new
            s = sum(new_xs[i,j] for j in 1:n_reac if j!=tot_col)
            new_xs[i,tot_col] = s
        end
    end
    PointwiseMaterial(pendf.mat, new_e, new_xs, copy(pendf.mt_list))
end

function _linterp(xs::AbstractVector, ys::AbstractVector, x::Float64)
    n = length(xs)
    x <= xs[1] && return ys[1]; x >= xs[n] && return ys[n]
    i = searchsortedfirst(xs, x); i == 1 && return ys[1]
    xs[i] == x && return ys[i]
    f = (x-xs[i-1])/(xs[i]-xs[i-1]); ys[i-1]+f*(ys[i]-ys[i-1])
end

"Compute thermal elastic XS on merged grid. Supports :free_gas and :sab."
function compute_thermal_xs(energies::AbstractVector{<:Real},
                            sigma_elastic::AbstractVector{<:Real},
                            A::Real, T::Real;
                            sigma_b::Union{Real,Nothing}=nothing, emax::Real=10.0,
                            model::Symbol=:free_gas,
                            sab_data::Union{SABData,Nothing}=nothing)
    @assert length(energies)==length(sigma_elastic) && issorted(energies)
    @assert T > 0 && A > 0
    if sigma_b === nothing
        idx = clamp(searchsortedfirst(energies, emax), 1, length(energies))
        sigma_b = A * sigma_elastic[idx]
    end
    em = T > 3000.0 ? emax*T/3000.0 : emax
    te = filter(e -> e <= em, THERMR_EGRID)
    txs = if model==:free_gas
        [free_gas_xs(e, A, T; sigma_b=Float64(sigma_b)) for e in te]
    elseif model==:sab
        sab_data===nothing && error("sab_data required for :sab")
        [sab_xs(e, sab_data, T) for e in te]
    else error("Unknown model: $model") end
    me = vcat(te, energies[energies .> em])
    mxs = vcat(txs, sigma_elastic[energies .> em])
    p = sortperm(me); (me[p], mxs[p])
end

# ---- Equi-probable angular bins (Fortran sigl, thermr.f90:2660-2872) ----

"""
    sigl_equiprobable(E, E_prime, nbin, kernel, tev; tol=0.001) -> (sigma, cosines)

Compute equi-probable cosine bins for scattering from E to E_prime.
`kernel(E, E', mu)` returns d²σ/(dE'dμ). Returns total cross section
and `nbin` equi-probable cosine values (average μ within each bin).

Matches Fortran thermr.f90 sigl (lines 2660-2872) with nlin < 0 path.
"""
function sigl_equiprobable(E::Float64, E_prime::Float64, nbin::Int,
                           kernel, tev::Float64; tol::Float64=0.001,
                           awr::Float64=1.0, tev_peak::Float64=tev)
    # Constants matching Fortran (thermr.f90:2682-2694)
    imax = 20; xtol = 1e-5; ytol = 1e-3; sigmin = 1e-32
    eps_peak = 1e-3; third = 0.333333333  # Fortran truncated third (line 2687)
    half_tol = 0.5 * tol

    # Scattering peak location (Fortran sigl lines 2693-2714)
    # Fortran precomputes seep=1/sqrt(e*ep) and multiplies (line 2710)
    b = abs(E_prime - E) / tev_peak
    s1bb = sqrt(1 + b*b)
    seep = E_prime > 0 ? 1.0 / sqrt(E * E_prime) : 0.0  # Fortran line 2700
    x_peak = E_prime > 0 ? 0.5*(E + E_prime - (s1bb - 1)*awr*tev) * seep : 0.0
    abs(x_peak) > 1 - eps_peak && (x_peak = 0.99)
    x_peak = round_sigfig(x_peak, 8, 0)

    # --- Phase 1: total cross section via adaptive linearization ---
    function sig_mu(mu)
        kernel(E, E_prime, mu)
    end

    # Adaptive stack (Fortran lines 2703-2744)
    xs = Vector{Float64}(undef, imax)
    ys = Vector{Float64}(undef, imax)
    xs[3] = -1.0; ys[3] = sig_mu(-1.0)
    xs[2] = x_peak; ys[2] = sig_mu(x_peak)
    xs[1] = 1.0; ys[1] = sig_mu(1.0)
    ymax = max(ys[1], ys[2], ys[3], eps_peak)

    function adaptive_integrate()
        i = 3; total = 0.0; xl = xs[3]; yl = ys[3]
        while true
            if i == imax
                @goto accept
            end
            xm = 0.5*(xs[i-1] + xs[i])
            xm = round_sigfig(xm, 8, 0)  # Fortran sigl line 2739: sigfig(xm,8,0)
            ym_linear = 0.5*(ys[i-1] + ys[i])
            yt = sig_mu(xm)
            test = half_tol * abs(yt) + half_tol * ymax / 50
            test2 = ym_linear + ymax / 100
            if abs(yt - ym_linear) <= test && abs(ys[i-1] - ys[i]) <= test2 &&
               (xs[i-1] - xs[i]) < 0.5
                @goto accept
            end
            xs[i-1] - xs[i] < xtol && @goto accept
            i += 1; i > imax && (i = imax; @goto accept)
            xs[i] = xs[i-1]; ys[i] = ys[i-1]
            xs[i-1] = xm; ys[i-1] = yt
            continue
            @label accept
            total += 0.5 * (ys[i] + yl) * (xs[i] - xl)
            xl = xs[i]; yl = ys[i]
            i -= 1
            i > 1 && continue
            i == 1 && (total += 0.5*(ys[1] + yl)*(xs[1] - xl); break)
            break
        end
        return total
    end

    sigma_total = adaptive_integrate()
    sigma_total < sigmin && return (0.0, zeros(nbin))

    # --- Phase 2: equi-probable bins via CDF inversion ---
    fract = sigma_total / nbin
    rfract = 1.0 / fract
    cosines = zeros(nbin)

    # Re-initialize stack for second adaptive pass
    xs[3] = -1.0; ys[3] = sig_mu(-1.0)
    xs[2] = x_peak; ys[2] = sig_mu(x_peak)
    xs[1] = 1.0; ys[1] = sig_mu(1.0)

    cum_sum = 0.0; gral = 0.0; j = 0
    xl = xs[3]; yl = ys[3]
    i = 3

    while j < nbin
        # Adaptive refinement (same convergence as phase 1)
        if i > 1 && i < imax
            xm = 0.5*(xs[i-1] + xs[i])
            xm = round_sigfig(xm, 8, 0)  # Fortran sigl line 2799: sigfig(xm,8,0)
            ym_linear = 0.5*(ys[i-1] + ys[i])
            yt = sig_mu(xm)
            test = half_tol * abs(yt) + half_tol * ymax / 50
            test2 = ym_linear + ymax / 100
            if !(abs(yt - ym_linear) <= test && abs(ys[i-1] - ys[i]) <= test2 &&
                 (xs[i-1] - xs[i]) < 0.5) && (xs[i-1] - xs[i]) >= xtol
                i += 1; i > imax && (i = imax)
                xs[i] = xs[i-1]; ys[i] = ys[i-1]
                xs[i-1] = xm; ys[i-1] = yt
                continue
            end
        end

        # Panel accepted: compute contribution
        panel_area = 0.5 * (ys[i] + yl) * (xs[i] - xl)
        xs[i] == xl && (i -= 1; i >= 1 && continue; break)

        # Check if this panel completes a bin
        if j == nbin - 1 && i == 1
            # Last bin: force to μ=+1
            xn = xs[i]; j += 1
        elseif cum_sum + panel_area >= fract * 0.99999999 && j < nbin - 1
            # Bin boundary within this panel: find xn via quadratic
            j += 1
            xil = 1.0 / (xs[i] - xl)
            f_slope = (ys[i] - yl) * xil
            if yl < sigmin || abs((fract - cum_sum) * f_slope / yl^2) <= ytol
                xn = xl + (fract - cum_sum) / max(yl, sigmin)
            else
                rf = 1.0 / f_slope
                disc = (yl * rf)^2 + 2 * (fract - cum_sum) * rf
                disc < 0 && (disc = abs(disc))
                xn = f_slope > 0 ? xl - yl*rf + sqrt(disc) : xl - yl*rf - sqrt(disc)
            end
            xn = clamp(xn, xl, xs[i])
        else
            # Accumulate and move to next panel
            cum_sum += panel_area
            gral += 0.5*(yl*xs[i] - ys[i]*xl)*(xs[i]+xl) +
                    third*(ys[i]-yl)*(xs[i]^2 + xs[i]*xl + xl^2)
            xl = xs[i]; yl = ys[i]
            i -= 1
            i > 1 && continue
            i == 1 && (continue)
            break
        end

        # Compute weighted average cosine for this bin
        yn = yl + (ys[i] - yl) * (xn - xl) / (xs[i] - xl)
        gral += (xn - xl) * (yl * 0.5*(xn+xl) +
                (ys[i]-yl)/(xs[i]-xl) * (-xl*0.5*(xn+xl) +
                third*(xn^2 + xn*xl + xl^2)))
        cosines[j] = gral * rfract

        # Reset for next bin
        xl = xn; yl = yn; cum_sum = 0.0; gral = 0.0
        j >= nbin && break
        xl < xs[i] && continue
        xl = xs[i]; yl = ys[i]; i -= 1
        i > 1 && continue
        i == 1 && continue
        break
    end

    return (sigma_total, cosines)
end

# ---- MF6 record data structure ----

"""MF6 LIST record for one incident energy: secondary energies with equi-probable cosines."""
struct MF6ListRecord
    E_incident::Float64
    entries::Vector{NTuple{10, Float64}}  # (E', σ, μ₁...μ₈) per secondary energy
end

# ---- Calcem: Fortran-style thermal XS + angular computation ----

"""
    calcem_xs(sab, T, emax; nep=500, nmu=40) -> (esi, xsi)

Compute total incoherent inelastic XS on the calcem egrid.
Fast path: trapezoidal integration over (E', μ) without angular bins.
Returns incident energies and total XS [barn].
"""
function calcem_xs(sab::SABData, T::Real, emax::Real;
                   nep::Int=500, nmu::Int=40)
    kT = PhysicsConstants.bk * T
    tevz = 0.0253
    kT_eff = sab.lat == 1 ? tevz : kT

    nne = something(findfirst(e -> e > emax, THERMR_EGRID), length(THERMR_EGRID))
    egrid = THERMR_EGRID[1:nne]

    kernel = (E, Ep, mu) -> sab_kernel(E, Ep, mu, sab, T)

    esi = zeros(nne)
    xsi = zeros(nne)

    for ie in 1:nne
        E = egrid[ie]
        esi[ie] = E
        ep_max = E + sab.beta[end] * kT_eff
        dep = ep_max / nep
        total = 0.0
        sp = 0.0
        for j in 1:nep
            Ep = j * dep
            # Trapezoidal angular integration
            sig = 0.0
            for imu in 0:nmu
                mu = -1.0 + 2.0 * imu / nmu
                w = (imu == 0 || imu == nmu) ? 0.5 : 1.0
                sig += w * kernel(E, Ep, mu) * (2.0 / nmu)
            end
            total += (sig + sp) * dep * 0.5
            sp = sig
        end
        xsi[ie] = round_sigfig(total, 9, 0)
    end

    return (esi, xsi)
end

"""
    calcem(sab, T, emax, nbin; tol=0.05) -> (esi, xsi, records)

Compute thermal inelastic XS and MF6 angular distributions on the calcem
egrid, matching Fortran thermr calcem (thermr.f90:1541-2480).

Returns: (esi, xsi, records) — energies, total XS, MF6 angular data.
"""
function calcem(sab::SABData, T::Real, emax::Real, nbin::Int;
                tol::Float64=0.05)
    kT = PhysicsConstants.bk * T
    tevz = 0.0253
    kT_eff = sab.lat == 1 ? tevz : kT
    tolmin_area = 5e-7  # Fortran tolmin (line 1611)
    imax_stack = 20     # max refinement depth (Fortran line 1573)
    nl = nbin + 1       # total components: sigma + nbin cosines (Fortran nl)

    nne = something(findfirst(e -> e > emax, THERMR_EGRID), length(THERMR_EGRID))
    egrid = THERMR_EGRID[1:nne]

    kernel = (E, Ep, mu) -> sab_kernel(E, Ep, mu, sab, T)

    esi = zeros(nne)
    xsi = zeros(nne)
    records = Vector{MF6ListRecord}()

    for ie in 1:nne
        E = round_sigfig(egrid[ie], 8, 0)  # Fortran enow=sigfig(enow,8,0) line 1979
        esi[ie] = E

        # Build beta-derived E' seed points (Fortran lines 1985-2031)
        # Fortran generates seeds in beta order: jbeta=-nbeta..-1, then 1..nbeta
        # When ep == enow: nudge to sigfig(enow,8,-1) for downscatter,
        # sigfig(enow,8,+1) for upscatter (with iskip=1 for that panel)
        seeds = Float64[0.0]  # E'=0 is always the first point
        ep_down_peak = -1.0  # elastic peak lower bound (sigfig(E,8,-1))
        ep_up_peak = -1.0    # elastic peak upper bound (sigfig(E,8,+1))
        for beta in sab.beta
            ep_down = E - beta * kT_eff
            if ep_down > 0
                if abs(ep_down - E) < 1e-15 * E  # ep == enow case
                    ep_down = round_sigfig(E, 8, -1)
                    ep_down_peak = ep_down
                else
                    ep_down = round_sigfig(ep_down, 8, 0)
                end
                push!(seeds, ep_down)
            end
            ep_up = E + beta * kT_eff
            if abs(ep_up - E) < 1e-15 * E  # ep == enow case
                ep_up = round_sigfig(E, 8, +1)
                ep_up_peak = ep_up
            else
                ep_up = round_sigfig(ep_up, 8, 0)
            end
            push!(seeds, ep_up)
        end
        sort!(unique!(seeds))

        # Evaluate sigl at each seed point
        seed_data = Vector{Tuple{Float64, Float64, Vector{Float64}}}()
        for Ep in seeds
            sigma, cosines = sigl_equiprobable(E, Ep, nbin, kernel, kT; tol=tol,
                                                awr=Float64(sab.awr), tev_peak=kT_eff)
            push!(seed_data, (Ep, sigma, cosines))
        end

        # Adaptive refinement between consecutive seed points (Fortran lines 2041-2162)
        entries = NTuple{10, Float64}[]
        total_xs = 0.0
        xlast = 0.0; ylast = 0.0
        j_count = 0  # total emitted points (like Fortran j)

        function emit_point!(ep, sigma, cosines)
            j_count += 1
            if j_count > 1  # Fortran: if (j.gt.1) — accumulate from 2nd point
                total_xs += (ep - xlast) * (sigma + ylast) * 0.5
            end

            # Fortran j==3 / xsi<tolmin (lines 2095-2096): skip early low-XS region
            if j_count == 3 && total_xs < tolmin_area
                j_count = 2
                xlast = ep; ylast = sigma
                # Overwrite previous entry (Fortran overwrites scr[j=2] position)
                entry = ntuple(k -> k == 1 ? ep : k == 2 ? sigma :
                               k - 2 <= nbin ? cosines[k-2] : 0.0, 10)
                if !isempty(entries)
                    entries[end] = entry
                else
                    push!(entries, entry)
                end
                return
            end

            xlast = ep; ylast = sigma
            # Store ALL entries (Fortran stores even sig=0; jnz+1 trim handles trailing zeros)
            entry = ntuple(k -> k == 1 ? ep : k == 2 ? sigma :
                           k - 2 <= nbin ? cosines[k-2] : 0.0, 10)
            push!(entries, entry)
        end

        # Process seed points with adaptive refinement between them
        for idx in 1:length(seed_data)
            ep_hi, sig_hi, cos_hi = seed_data[idx]

            if idx > 1
                ep_lo, sig_lo, cos_lo = seed_data[idx-1]

                # Skip refinement for elastic peak panel (Fortran iskip)
                # Only skip the narrow [sigfig(E,8,-1), sigfig(E,8,+1)] panel
                is_elastic_panel = ep_down_peak > 0 && ep_up_peak > 0 &&
                    abs(ep_lo - ep_down_peak) < 1e-15 * E &&
                    abs(ep_hi - ep_up_peak) < 1e-15 * E
                if is_elastic_panel
                    emit_point!(ep_lo, sig_lo, cos_lo)
                else
                    # Convergence stack between ep_lo and ep_hi
                    stk_e = zeros(imax_stack)
                    stk_s = zeros(imax_stack)
                    stk_c = [zeros(nbin) for _ in 1:imax_stack]
                    stk_e[1] = ep_hi; stk_s[1] = sig_hi; stk_c[1] .= cos_hi
                    stk_e[2] = ep_lo; stk_s[2] = sig_lo; stk_c[2] .= cos_lo
                    depth = 2

                    while depth >= 2
                        # Check area threshold (Fortran line 2050)
                        area = 0.5 * (stk_s[depth] + stk_s[depth-1]) * (stk_e[depth-1] - stk_e[depth])
                        if depth >= imax_stack || area < tolmin_area
                            emit_point!(stk_e[depth], stk_s[depth], stk_c[depth])
                            depth -= 1
                            continue
                        end

                        xm = round_sigfig(0.5 * (stk_e[depth-1] + stk_e[depth]), 8, 0)
                        if xm <= stk_e[depth] || xm >= stk_e[depth-1]
                            emit_point!(stk_e[depth], stk_s[depth], stk_c[depth])
                            depth -= 1
                            continue
                        end

                        # Evaluate at midpoint
                        sig_m, cos_m = sigl_equiprobable(E, xm, nbin, kernel, kT; tol=tol,
                                                         awr=Float64(sab.awr), tev_peak=kT_eff)

                        # Convergence tests (Fortran lines 2057-2068)
                        # k=1 sigma test: test2 = tol*abs(yt(1)), no floor
                        ym_s = 0.5 * (stk_s[depth] + stk_s[depth-1])
                        pass = abs(sig_m - ym_s) <= tol * abs(sig_m)

                        # k>1 per-component cosine tests + accumulate for integral test
                        uu = 0.0; uum = 0.0
                        if pass
                            for k in 1:nbin
                                ym_c = 0.5 * (stk_c[depth][k] + stk_c[depth-1][k])
                                uu += cos_m[k]
                                uum += ym_c
                                if abs(cos_m[k] - ym_c) > tol
                                    pass = false; break
                                end
                            end
                        end

                        # Integral cosine test (Fortran lines 2067-2068)
                        if pass
                            test_uu = 2.0 * tol * abs(uu) + 1e-5  # uumin = 0.00001
                            if abs(uu - uum) > test_uu
                                pass = false
                            end
                        end

                        if pass
                            emit_point!(stk_e[depth], stk_s[depth], stk_c[depth])
                            depth -= 1
                        else
                            # Subdivide: push midpoint
                            depth += 1
                            if depth > imax_stack; depth = imax_stack; end
                            stk_e[depth] = stk_e[depth-1]
                            stk_s[depth] = stk_s[depth-1]
                            stk_c[depth] .= stk_c[depth-1]
                            stk_e[depth-1] = xm
                            stk_s[depth-1] = sig_m
                            stk_c[depth-1] .= cos_m
                        end
                    end
                end
            end

            # Emit the current seed point (last in pair becomes first of next)
            if idx == length(seed_data)
                emit_point!(ep_hi, sig_hi, cos_hi)
            end
        end

        # Trim trailing zero-sigma entries (Fortran lines 2133, 2206-2207)
        jnz = 0
        for k in 1:length(entries)
            entries[k][2] != 0.0 && (jnz = k)
        end
        if jnz > 0 && jnz + 1 < length(entries)
            resize!(entries, jnz + 1)
        end

        xsi[ie] = round_sigfig(total_xs, 9, 0)
        push!(records, MF6ListRecord(E, entries))
    end

    return (esi, xsi, records)
end

# Fortran free-gas beta grid (thermr.f90:1858-1911, 45 values)
const FREE_GAS_BETA = Float64[
    0, 0.1, 2, 4, 6, 8, 10, 15, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
    100, 120, 140, 160, 180, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000,
    1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3500
]

"""
    calcem_free_gas(A, T, emax, nbin; sigma_b, tol=0.05)

Compute free-gas MF6 angular distributions using the Fortran's hardcoded
45-point beta grid with adaptive refinement. Matches calcem for iinc=1.
"""
function calcem_free_gas(A::Float64, T::Float64, emax::Float64, nbin::Int;
                         sigma_b::Float64=1.0, tol::Float64=0.05)
    kT = PhysicsConstants.bk * T
    tolmin_area = 5e-7
    imax_stack = 20

    nne = something(findfirst(e -> e > emax, THERMR_EGRID), length(THERMR_EGRID))
    egrid = THERMR_EGRID[1:nne]

    kernel = (E, Ep, mu) -> free_gas_kernel(E, Ep, mu, A, T; sigma_b=sigma_b)

    records = Vector{MF6ListRecord}()
    xsi_out = zeros(nne)

    for ie in 1:nne
        E = round_sigfig(egrid[ie], 8, 0)  # Fortran enow=sigfig(enow,8,0) line 1979
        nbeta = length(FREE_GAS_BETA)

        # Match Fortran calcem sequential beta processing (lines 1985-2152)
        entries = NTuple{10, Float64}[]
        total_xs = 0.0
        xlast = 0.0; ylast = 0.0
        j_count = 0

        function emit_point!(ep, sig, cosv)
            j_count += 1
            if j_count > 1
                total_xs += (ep - xlast) * (sig + ylast) * 0.5
            end
            if j_count == 3 && total_xs < tolmin_area
                j_count = 2
                xlast = ep; ylast = sig
                entry = ntuple(k -> k == 1 ? ep : k == 2 ? sig :
                               k - 2 <= nbin ? cosv[k-2] : 0.0, 10)
                !isempty(entries) ? (entries[end] = entry) : push!(entries, entry)
                return
            end
            xlast = ep; ylast = sig
            push!(entries, ntuple(k -> k == 1 ? ep : k == 2 ? sig :
                           k - 2 <= nbin ? cosv[k-2] : 0.0, 10))
        end

        # Convergence stack subroutine (Fortran labels 330-410)
        function run_stack!(x1_e, x1_s, x1_c, x2_e, x2_s, x2_c, do_iskip)
            if do_iskip
                emit_point!(x2_e, x2_s, x2_c)
                return
            end
            stk_e = zeros(imax_stack); stk_s = zeros(imax_stack)
            stk_c = [zeros(nbin) for _ in 1:imax_stack]
            stk_e[1] = x1_e; stk_s[1] = x1_s; stk_c[1] .= x1_c
            stk_e[2] = x2_e; stk_s[2] = x2_s; stk_c[2] .= x2_c
            depth = 2
            while depth >= 2
                area = 0.5 * (stk_s[depth] + stk_s[depth-1]) * (stk_e[depth-1] - stk_e[depth])
                if depth >= imax_stack || area < tolmin_area
                    emit_point!(stk_e[depth], stk_s[depth], stk_c[depth])
                    depth -= 1; continue
                end
                xm = round_sigfig(0.5 * (stk_e[depth-1] + stk_e[depth]), 8, 0)
                if xm <= stk_e[depth] || xm >= stk_e[depth-1]
                    emit_point!(stk_e[depth], stk_s[depth], stk_c[depth])
                    depth -= 1; continue
                end
                sig_m, cos_m = sigl_equiprobable(E, xm, nbin, kernel, kT; tol=tol, awr=A)
                ym_s = 0.5 * (stk_s[depth] + stk_s[depth-1])
                pass = abs(sig_m - ym_s) <= tol * abs(sig_m)
                uu = 0.0; uum = 0.0
                if pass
                    for k in 1:nbin
                        ym_c = 0.5 * (stk_c[depth][k] + stk_c[depth-1][k])
                        uu += cos_m[k]; uum += ym_c
                        if abs(cos_m[k] - ym_c) > tol
                            pass = false; break
                        end
                    end
                end
                if pass
                    test_uu = 2.0 * tol * abs(uu) + 1e-5
                    abs(uu - uum) > test_uu && (pass = false)
                end
                if pass
                    emit_point!(stk_e[depth], stk_s[depth], stk_c[depth])
                    depth -= 1
                else
                    depth += 1
                    depth > imax_stack && (depth = imax_stack)
                    stk_e[depth] = stk_e[depth-1]; stk_s[depth] = stk_s[depth-1]
                    stk_c[depth] .= stk_c[depth-1]
                    stk_e[depth-1] = xm; stk_s[depth-1] = sig_m; stk_c[depth-1] .= cos_m
                end
            end
        end

        # Initialize with E'=0 (Fortran line 1985-1988)
        x1_e = 0.0
        x1_s, x1_c = sigl_equiprobable(E, 0.0, nbin, kernel, kT; tol=tol, awr=A)
        iskip = false

        # Sequential beta loop matching Fortran (lines 1992-2148)
        jbeta = -nbeta
        while true
            jbeta == 0 && (jbeta = 1)

            # Compute next E' from beta (Fortran lines 2005-2027)
            bi = abs(jbeta)
            beta_val = FREE_GAS_BETA[bi]
            if jbeta < 0
                ep = E - beta_val * kT
                ep = round_sigfig(ep, 8, 0) == round_sigfig(E, 8, 0) ?
                     round_sigfig(E, 8, -1) : round_sigfig(ep, 8, 0)
            else
                ep = E + beta_val * kT
                if round_sigfig(ep, 8, 0) == round_sigfig(E, 8, 0)
                    ep = round_sigfig(E, 8, +1)
                    iskip = true
                else
                    ep = round_sigfig(ep, 8, 0)
                end
            end

            # Skip if ep <= previous right endpoint (Fortran line 2029)
            if ep <= x1_e
                jbeta += 1
                jbeta == 0 && (jbeta = 1)
                jbeta > nbeta && break
                continue
            end
            ep = round_sigfig(ep, 8, 0)  # Fortran line 2033

            # Set up panel: x2=previous, x1=new seed (Fortran lines 1999-2038)
            x2_e = x1_e; x2_s = x1_s; x2_c = copy(x1_c)
            x1_e = ep
            x1_s, x1_c = sigl_equiprobable(E, ep, nbin, kernel, kT; tol=tol, awr=A)

            # Run convergence stack (Fortran labels 330-410)
            run_stack!(x1_e, x1_s, x1_c, x2_e, x2_s, x2_c, iskip)
            iskip = false

            # Advance beta (Fortran line 2147)
            jbeta += 1
            jbeta == 0 && (jbeta = 1)
            jbeta > nbeta && break
        end

        # Emit final zero point (Fortran lines 2149-2152 → label 430)
        emit_point!(x1_e, 0.0, zeros(nbin))

        # Trim trailing zeros (Fortran lines 2204-2205)
        jnz = 0
        for k in 1:length(entries)
            entries[k][2] != 0.0 && (jnz = k)
        end
        if jnz > 0 && jnz + 1 < length(entries)
            resize!(entries, jnz + 1)
        end

        xsi_out[ie] = round_sigfig(total_xs, 9, 0)
        push!(records, MF6ListRecord(E, entries))
    end
    return (egrid, xsi_out, records)
end

"""
    compute_mf6_thermal(pendf_energies, pendf_elastic, A, T, nbin;
                        model=:free_gas, sab_data=nothing, mtref=221)

Compute MF6 angular distribution data for thermal inelastic scattering.
Returns Vector{MF6ListRecord}, one per incident energy.
"""
function compute_mf6_thermal(pendf_energies::Vector{Float64},
                             pendf_elastic::Vector{Float64},
                             A::Float64, T::Float64, nbin::Int;
                             model::Symbol=:free_gas,
                             sab_data::Union{SABData, Nothing}=nothing,
                             sigma_b::Union{Float64, Nothing}=nothing,
                             emax::Float64=0.0)
    kT = PhysicsConstants.bk * T
    tev = kT
    emax = emax > 0 ? emax : 10.0 * kT  # default: 10*kT for backward compat

    # Bound scattering cross section
    if sigma_b === nothing
        idx = clamp(searchsortedfirst(pendf_energies, 10.0), 1, length(pendf_energies))
        sigma_b = A * pendf_elastic[idx]
    end

    # Kernel function
    kernel = if model == :free_gas
        (E, Ep, mu) -> free_gas_kernel(E, Ep, mu, A, T; sigma_b=sigma_b)
    elseif model == :sab && sab_data !== nothing
        (E, Ep, mu) -> sab_kernel(E, Ep, mu, sab_data, T)
    else
        error("Invalid model=$model or missing sab_data")
    end

    # Use the standard thermr energy grid for incident energies (up to emax)
    nne = something(findfirst(e -> e > emax, THERMR_EGRID), length(THERMR_EGRID))
    egrid = THERMR_EGRID[1:nne]

    records = MF6ListRecord[]
    for E in egrid
        # Secondary energy grid from beta values (Fortran calcem approach)
        ep_grid = Float64[]
        beta_grid = model == :free_gas ? FREE_GAS_BETA :
                    (sab_data !== nothing ? sab_data.beta : FREE_GAS_BETA)
        kT_eff = (model == :sab && sab_data !== nothing && sab_data.lat == 1) ? 0.0253 : kT
        for beta in beta_grid
            ep_down = E - beta * kT_eff
            ep_down > 0 && push!(ep_grid, round_sigfig(ep_down, 8, 0))
            ep_up = E + beta * kT_eff
            push!(ep_grid, round_sigfig(ep_up, 8, 0))
        end
        sort!(unique!(ep_grid))
        filter!(ep -> ep > 0, ep_grid)

        entries = NTuple{10, Float64}[]
        for Ep in ep_grid
            sigma, cosines = sigl_equiprobable(E, Ep, nbin, kernel, tev; tol=0.001)
            entry = ntuple(k -> k == 1 ? Ep : k == 2 ? sigma :
                           k-2 <= nbin ? cosines[k-2] : 0.0, 10)
            push!(entries, entry)
        end
        push!(records, MF6ListRecord(E, entries))
    end
    return records
end
