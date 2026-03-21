# Standard neutron group structures for GROUPR processing
# Correspondence to NJOY2016 groupr.f90 gengpn: ign -> NTuple of bounds [eV].
# Bounds are in ascending energy order; bounds[g] = lower edge of group g.

"LANL 30-group neutron energy structure [eV], 31 bounds."
const LANL_30 = NTuple{31,Float64}((
    1.39e-4,  1.52e-1,  4.14e-1,  1.13e0,   3.06e0,   8.32e0,
    2.26e1,   6.14e1,   1.67e2,   4.54e2,   1.235e3,  3.35e3,
    9.12e3,   2.48e4,   6.76e4,   1.84e5,   3.03e5,   5.00e5,
    8.23e5,   1.353e6,  1.738e6,  2.232e6,  2.865e6,  3.68e6,
    6.07e6,   7.79e6,   1.00e7,   1.20e7,   1.35e7,   1.50e7,
    1.70e7,
))

"EPRI-CPM/WIMS 69-group neutron energy structure [eV], 70 bounds."
const WIMS_69 = NTuple{70,Float64}((
    1.0e-5,    0.005e0,   0.01e0,    0.015e0,   0.02e0,    0.025e0,
    0.03e0,    0.035e0,   0.042e0,   0.05e0,    0.058e0,   0.067e0,
    0.08e0,    0.1e0,     0.14e0,    0.18e0,    0.22e0,    0.25e0,
    0.28e0,    0.3e0,     0.32e0,    0.35e0,    0.4e0,     0.5e0,
    0.625e0,   0.78e0,    0.85e0,    0.91e0,    0.95e0,    0.972e0,
    0.996e0,   1.02e0,    1.045e0,   1.071e0,   1.097e0,   1.123e0,
    1.15e0,    1.3e0,     1.5e0,     2.1e0,     2.6e0,     3.3e0,
    4.0e0,     9.877e0,   15.968e0,  27.7e0,    48.052e0,  75.501e0,
    148.728e0, 367.262e0, 906.898e0, 1425.1e0,  2239.45e0, 3519.1e0,
    5530.0e0,  9118.0e0,  1.503e4,   2.478e4,   4.085e4,   6.734e4,
    1.11e5,    1.83e5,    3.025e5,   5.0e5,     8.21e5,    1.353e6,
    2.231e6,   3.679e6,   6.0655e6,  1.0e7,
))

# VITAMIN-J 175-group (ign=17): eg15a[1:84] + eg15b[1:81] + e175 + eg15b[82:91]
"VITAMIN-J 175-group neutron energy structure [eV] (ORNL-5510), 176 bounds."
const VITAMINJ_175 = NTuple{176,Float64}((
    # eg15a[1:84]
    1.0e-5,     1.0e-1,     4.1399e-1,  5.3158e-1,  6.8256e-1,
    8.7642e-1,  1.1254e0,   1.4450e0,   1.8554e0,   2.3824e0,
    3.0590e0,   3.9279e0,   5.0435e0,   6.4760e0,   8.3153e0,
    1.0677e1,   1.3710e1,   1.7603e1,   2.2603e1,   2.9023e1,
    3.7267e1,   4.7851e1,   6.1442e1,   7.8893e1,   1.0130e2,
    1.3007e2,   1.6702e2,   2.1445e2,   2.7536e2,   3.5358e2,
    4.5400e2,   5.8295e2,   7.4852e2,   9.6112e2,   1.2341e3,
    1.5846e3,   2.0347e3,   2.2487e3,   2.4852e3,   2.6126e3,
    2.7465e3,   3.0354e3,   3.3546e3,   3.7074e3,   4.3074e3,
    5.5308e3,   7.1017e3,   9.1188e3,   1.0595e4,   1.1709e4,
    1.5034e4,   1.9305e4,   2.1875e4,   2.3579e4,   2.4176e4,
    2.4788e4,   2.6058e4,   2.7000e4,   2.8500e4,   3.1828e4,
    3.4307e4,   4.0868e4,   4.6309e4,   5.2475e4,   5.6562e4,
    6.7379e4,   7.2000e4,   7.9500e4,   8.2500e4,   8.6517e4,
    9.8037e4,   1.1109e5,   1.1679e5,   1.2277e5,   1.2907e5,
    1.3569e5,   1.4264e5,   1.4996e5,   1.5764e5,   1.6573e5,
    1.7422e5,   1.8316e5,   1.9255e5,   2.0242e5,
    # eg15b[1:81]
    2.1280e5,   2.2371e5,   2.3518e5,   2.4724e5,   2.7324e5,
    2.8725e5,   2.9452e5,   2.9720e5,   2.9850e5,   3.0197e5,
    3.3373e5,   3.6883e5,   3.8774e5,   4.0762e5,   4.5049e5,
    4.9787e5,   5.2340e5,   5.5023e5,   5.7844e5,   6.0810e5,
    6.3928e5,   6.7206e5,   7.0651e5,   7.4274e5,   7.8082e5,
    8.2085e5,   8.6294e5,   9.0718e5,   9.6164e5,   1.0026e6,
    1.1080e6,   1.1648e6,   1.2246e6,   1.2873e6,   1.3534e6,
    1.4227e6,   1.4957e6,   1.5724e6,   1.6530e6,   1.7377e6,
    1.8268e6,   1.9205e6,   2.0190e6,   2.1225e6,   2.2313e6,
    2.3069e6,   2.3457e6,   2.3653e6,   2.3852e6,   2.4660e6,
    2.5924e6,   2.7253e6,   2.8650e6,   3.0119e6,   3.1664e6,
    3.3287e6,   3.6788e6,   4.0657e6,   4.4933e6,   4.7237e6,
    4.9659e6,   5.2205e6,   5.4881e6,   5.7695e6,   6.0653e6,
    6.3763e6,   6.5924e6,   6.7032e6,   7.0469e6,   7.4082e6,
    7.7880e6,   8.1873e6,   8.6071e6,   9.0484e6,   9.5123e6,
    1.0000e7,   1.0513e7,   1.1052e7,   1.1618e7,   1.2214e7,
    1.2523e7,
    # e175 inserted at position 166
    1.284e7,
    # eg15b[82:91]  (shifted down by 1 from eg15b indexing)
    1.3499e7,   1.3840e7,   1.4191e7,   1.4550e7,
    1.4918e7,   1.5683e7,   1.6487e7,   1.6905e7,   1.7333e7,
    1.9640e7,
))

# SAND-II 620-group (ign=12): built algorithmically matching gengpn logic.
function _build_sandii_620()
    deltl = (5.0, 7.5, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0)
    ndelta = (2, 6, 10, 19, 23, 28, 36, 40, 46)
    sandb = 1.0e-6
    sanda = 1.0e-4
    sandc = 2.8e-4
    sandd = 1.0e6
    sande = 1.0e5
    ngp = 621
    egn = Vector{Float64}(undef, ngp)
    egn[1] = sanda
    for ig in 1:8
        delta = deltl[ig] * sandb
        n1 = ndelta[ig]
        n2 = ndelta[ig + 1] - 1
        for n in n1:n2
            egn[n] = egn[n - 1] + delta
        end
    end
    egn[21] = sandc
    for ig in 46:450
        egn[ig] = egn[ig - 45] * 10.0
    end
    egn[451] = sandd
    for ig in 452:ngp
        egn[ig] = egn[ig - 1] + sande
    end
    return NTuple{ngp, Float64}(Tuple(egn))
end

"SAND-II 620-group neutron energy structure [eV], 621 bounds."
const SANDII_620 = _build_sandii_620()

# XMAS 172-group (ign=18): eg18 reversed to ascending order.
"XMAS 172-group neutron energy structure [eV], 173 bounds."
const XMAS_172 = NTuple{173,Float64}((
    1.00001e-5, 3.00000e-3, 5.00000e-3, 6.90000e-3,
    1.00000e-2, 1.50000e-2, 2.00000e-2, 2.50000e-2,
    3.00000e-2, 3.50000e-2, 4.20000e-2, 5.00000e-2,
    5.80000e-2, 6.70000e-2, 7.70000e-2, 8.00000e-2,
    9.50000e-2, 1.00001e-1, 1.15000e-1, 1.34000e-1,
    1.40000e-1, 1.60000e-1, 1.80000e-1, 1.89000e-1,
    2.20000e-1, 2.48000e-1, 2.80000e-1, 3.00000e-1,
    3.14500e-1, 3.20000e-1, 3.50000e-1, 3.91000e-1,
    4.00000e-1, 4.33000e-1, 4.85000e-1, 5.00000e-1,
    5.40000e-1, 6.25000e-1, 7.05000e-1, 7.80000e-1,
    7.90000e-1, 8.50000e-1, 8.60000e-1, 9.10000e-1,
    9.30000e-1, 9.50000e-1, 9.72000e-1, 9.86000e-1,
    9.96000e-1, 1.02000e0, 1.03500e0, 1.04500e0,
    1.07100e0, 1.09700e0, 1.11000e0, 1.12535e0,
    1.15000e0, 1.17000e0, 1.23500e0, 1.30000e0,
    1.33750e0, 1.37000e0, 1.44498e0, 1.47500e0,
    1.50000e0, 1.59000e0, 1.67000e0, 1.75500e0,
    1.84000e0, 1.93000e0, 2.02000e0, 2.10000e0,
    2.13000e0, 2.36000e0, 2.55000e0, 2.60000e0,
    2.72000e0, 2.76792e0, 3.30000e0, 3.38075e0,
    4.00000e0, 4.12925e0, 5.04348e0, 5.34643e0,
    6.16012e0, 7.52398e0, 8.31529e0, 9.18981e0,
    9.90555e0, 1.12245e1, 1.37096e1, 1.59283e1,
    1.94548e1, 2.26033e1, 2.49805e1, 2.76077e1,
    3.05113e1, 3.37201e1, 3.72665e1, 4.01690e1,
    4.55174e1, 4.82516e1, 5.15780e1, 5.55951e1,
    6.79041e1, 7.56736e1, 9.16609e1, 1.36742e2,
    1.48625e2, 2.03995e2, 3.04325e2, 3.71703e2,
    4.53999e2, 6.77287e2, 7.48518e2, 9.14242e2,
    1.01039e3, 1.23410e3, 1.43382e3, 1.50733e3,
    2.03468e3, 2.24867e3, 3.35463e3, 3.52662e3,
    5.00451e3, 5.53084e3, 7.46586e3, 9.11882e3,
    1.11378e4, 1.50344e4, 1.66156e4, 2.47875e4,
    2.73944e4, 2.92830e4, 3.69786e4, 4.08677e4,
    5.51656e4, 6.73795e4, 8.22975e4, 1.11090e5,
    1.22773e5, 1.83156e5, 2.47235e5, 2.73237e5,
    3.01974e5, 4.07622e5, 4.50492e5, 4.97871e5,
    5.50232e5, 6.08101e5, 8.20850e5, 9.07180e5,
    1.00259e6, 1.10803e6, 1.22456e6, 1.35335e6,
    1.65299e6, 2.01897e6, 2.23130e6, 2.46597e6,
    3.01194e6, 3.67879e6, 4.49329e6, 5.48812e6,
    6.06531e6, 6.70320e6, 8.18731e6, 1.00000e7,
    1.16183e7, 1.38403e7, 1.49182e7, 1.73325e7,
    1.96403e7,
))

"ECCO 33-group neutron energy structure [eV], 34 bounds."
const ECCO_33 = NTuple{34,Float64}((
    1.000010e-05, 1.000000e-01, 5.400000e-01, 4.000000e+00,
    8.315287e+00, 1.370959e+01, 2.260329e+01, 4.016900e+01,
    6.790405e+01, 9.166088e+01, 1.486254e+02, 3.043248e+02,
    4.539993e+02, 7.485183e+02, 1.234098e+03, 2.034684e+03,
    3.354626e+03, 5.530844e+03, 9.118820e+03, 1.503439e+04,
    2.478752e+04, 4.086771e+04, 6.737947e+04, 1.110900e+05,
    1.831564e+05, 3.019738e+05, 4.978707e+05, 8.208500e+05,
    1.353353e+06, 2.231302e+06, 3.678794e+06, 6.065307e+06,
    1.000000e+07, 1.964033e+07,
))

"Named group structure identifier for the standard built-in structures."
@enum GroupStructureId begin
    IGN_LANL30     = 3
    IGN_WIMS69     = 9
    IGN_SANDII620  = 12
    IGN_VITAMINJ   = 17
    IGN_XMAS172    = 18
    IGN_ECCO33     = 19
end

"Return group boundary tuple for built-in structure. Ascending energy [eV]."
function get_group_structure(ign::GroupStructureId)
    ign == IGN_LANL30    && return LANL_30
    ign == IGN_WIMS69    && return WIMS_69
    ign == IGN_SANDII620 && return SANDII_620
    ign == IGN_VITAMINJ  && return VITAMINJ_175
    ign == IGN_XMAS172   && return XMAS_172
    ign == IGN_ECCO33    && return ECCO_33
    error("Unknown group structure: $ign")
end

get_group_structure(ign::Integer) = get_group_structure(GroupStructureId(ign))

"Number of energy groups from a bounds collection (length - 1)."
num_groups(bounds) = length(bounds) - 1

"Check that group bounds are strictly ascending and positive."
function validate_group_bounds(bounds)
    n = length(bounds)
    n >= 2 || throw(ArgumentError("need at least 2 bounds (1 group)"))
    for i in 1:n-1
        bounds[i] < bounds[i+1] ||
            throw(ArgumentError("bounds not strictly ascending at index $i: " *
                                "$(bounds[i]) >= $(bounds[i+1])"))
    end
    bounds[1] > 0 || throw(ArgumentError("bounds must be positive; got $(bounds[1])"))
    nothing
end

"Find group index g such that bounds[g] <= E < bounds[g+1]; 0 if out of range."
function find_group(bounds, E::Real)
    ng = length(bounds) - 1
    (E < bounds[1] || E >= bounds[end]) && return 0
    lo, hi = 1, ng
    while lo < hi
        mid = (lo + hi + 1) >> 1
        bounds[mid] <= E ? (lo = mid) : (hi = mid - 1)
    end
    lo
end
