# REVIEW: Wave 3 THERMR -- Skeptical Review

**File:** `src/processing/thermr.jl` (292 lines)
**Tests:** `test/runtests.jl` lines 2167-2447 (281 lines)
**Fortran reference:** `njoy-reference/src/thermr.f90` (3387 lines)
**All 6690 tests PASS.**

---

## 1. Diff Against Fortran

### Free gas cross section formula (free_gas_xs)
**PASS.** Analytical formula `sigma = (sigma_b/A)*[(1+1/(2x^2))*erf(x) + exp(-x^2)/(x*sqrt(pi))]`
where `x = sqrt(A*E/kT)` is the standard free-gas total thermal XS. Not directly in
NJOY2016's `sig()` (which computes the differential kernel), but is the correct analytical
integral. The Fortran code computes this implicitly via numerical integration of `sig()`.

### Free gas kernel (free_gas_kernel vs sig()@200)
**PASS.** Line-by-line match:
- `alpha = (E+Ep-2*mu*sqrt(E*Ep))/(A*kT)` matches Fortran `a=(e+ep-2*u*sqrt(e*ep))/(az*tev)`
- `beta = (Ep-E)/kT` matches Fortran `bb=(ep-e)*rtev` (signed, not absolute)
- `arg = (alpha+beta)^2/(4*alpha)` matches Fortran `arg=(a+bb)**2/(4*a)`
- `sigma_b*sqrt(Ep/E)/(2*kT)*exp(-arg)/sqrt(4*pi*alpha)` matches Fortran `sigc*sb*exp(-arg)/(c*sqrt(a))`
- Minimum alpha clamped to `1e-10` (Julia) vs `1e-6` (Fortran `amin`). Minor difference, acceptable.
- Cutoff at `arg > 225` matches Fortran `sabflg = -225`.

### S(alpha,beta) kernel (sab_kernel vs sig()@2514)
**PASS with caveats.**
- Prefactor `sqrt(Ep/E)/(2*kT)*sb` matches Fortran `sigc*sb`.
- `sb = sigma_b*((awr+1)/awr)^2` matches Fortran `sb=smz*((az+1)/az)**2`.
- Symmetrization `exp(log_S_sym - beta_raw/2)` matches Fortran `exp(s-bb/2)`.
- `lat==1` scaling by `kT/tevz` (where `tevz=0.0253`) matches Fortran.
- SCT fallback `arg = (|beta|-alpha)^2*kT/(4*alpha*teff) + (|beta|+beta_raw)/2` matches Fortran label 170.

**Caveats:**
- **Interpolation:** Julia uses bilinear; Fortran uses quadratic (`terpq`). Lower accuracy on coarse SAB tables.
- **Missing cliq extrapolation:** Fortran small-alpha extrapolation (lines 2539-2544) not implemented.
- **Missing sb2/az2 second-atom SCT term:** Fortran lines 2588-2594 handle materials with two SCT atoms. Julia only supports one.

### Bragg edges (build_bragg_data + bragg_edges vs sigcoh())
**FAIL -- Missing structure factor.**
- `econ`, `scon`, `wint`, `t2` formulas all match Fortran exactly.
- `tausq` formula matches Fortran's `tausq()` contained function.
- Loop structure over `(l1, l2, l3)` with `sgn` for `+/-l2` matches Fortran.
- Weight multiplicities `w1, w2, w3` match Fortran.
- Grouping of nearby shells (5% tolerance) matches Fortran.
- **BUT:** Julia does NOT multiply by `form(lat,l1,l2,l3)`. The Fortran computes `f = w*form(lat,l1,l2,l3)` where `form()` gives crystal structure factors for graphite, Be, and BeO. Julia stores just `w` (the Debye-Waller weighted multiplicity) without the structure factor. This produces incorrect Bragg edge intensities for any real material. The `bragg_edges()` evaluation function itself correctly implements `scon*sum(f)/E`, matching Fortran label 210.

### Incoherent elastic (incoh_elastic_xs vs iel())
**PASS.** Formula `sigma(E) = sigma_b/2*(1-exp(-4*E*W'))` matches the ENDF-6 manual MF7/MT2 standard formula. The Fortran `iel()` computes both the total XS and angular distributions; Julia correctly implements the total XS portion.

### Standard thermal energy grid (THERMR_EGRID)
**PASS.** 118 entries, all values match Fortran `egrid(ngrid)` to full precision.

---

## 2. Physical Invariants

### Detailed balance: k(E->E')/k(E'->E) = (E'/E)*exp(-(E'-E)/kT)
**PASS.** Exact by construction for both `free_gas_kernel` and `sab_kernel`.
The ratio cancels to `(Ep/E)*exp(-beta)` algebraically. Tested with `rtol=1e-10`.

### Free gas high-E limit: sigma -> sigma_b/A
**PASS.** Tested at `E = 1000*kT` and `E = 10000*kT` with `rtol=1e-4`.

### Free gas low-E limit: 1/v behavior
**PASS.** For `x = sqrt(A*E/kT) < 1e-6`, returns `sigma_b/A*2/(x*sqrt(pi))`.
Leading-order expansion of the full formula confirms this is correct.

### Bragg: sigma = 0 below first edge
**PASS.** Loop breaks before accumulating any form factors when all `tau_sq[i] >= E*econ`.

---

## 3. Edge Cases, AD Compatibility, File Size

### Edge cases
- **E=0, Ep=0, kT=0:** All kernels return 0.0. Correct.
- **Large arg cutoff:** `arg > 225` returns 0.0, matching Fortran `sabflg = -225`.
- **Alpha clamping:** `max(alpha, 1e-10)` in free gas, `max(alpha, 1e-6)` in SAB. Prevents division by zero.
- **SAB table boundary:** Returns `sabflg` when alpha or beta exceed table range. Correct.

### AD compatibility
- `free_gas_xs`: **Compatible** (ForwardDiff). Uses `erf`, `exp`, `sqrt` -- all differentiable.
- `free_gas_kernel`: **Partially compatible.** `Float64()` casts will strip ForwardDiff `Dual` numbers. Would need to remove explicit Float64 conversion for full AD support.
- `sab_kernel`: **Partially compatible.** Same `Float64()` cast issue. `searchsortedlast` is integer-valued (OK for ForwardDiff on continuous parameters like E).
- `bragg_edges`: **Compatible** (ForwardDiff). Simple arithmetic.
- `incoh_elastic_xs`: **Fully compatible.** Pure arithmetic.
- `build_bragg_data`: **Not AD-compatible** (precomputation step; push!, sortperm). This is acceptable since it is called once at setup, not during evaluation.

### File size
292 lines for thermr.jl. Compact and well-organized. The Fortran is 3387 lines, so this is roughly 8.6% of the original, which is appropriate given that Julia omits I/O, tape management, and angular distribution binning.

---

## 4. Test Coverage

| Component | Tests | Verdict |
|-----------|-------|---------|
| `free_gas_xs` high-E limit | Yes (rtol=1e-4) | PASS |
| `free_gas_xs` low-E limit | Yes (rtol=1e-2) | PASS |
| `free_gas_xs` positivity | Yes (4 A x 4 T x 5 E) | PASS |
| `free_gas_xs` monotonicity | Yes (7 energies) | PASS |
| `free_gas_kernel` detailed balance | Yes (4 points, rtol=1e-10) | PASS |
| `free_gas_kernel` non-negativity | Yes (3x3x5 grid) | PASS |
| `THERMR_EGRID` | Yes (sorted, count, bounds) | PASS |
| `compute_thermal_xs` | Yes (basic functionality) | PASS |
| `sab_kernel` vs free gas | Yes (1 point, rtol=0.05) | PASS |
| `sab_kernel` detailed balance | Yes (2 points, rtol=1e-10) | PASS |
| `read_thermal_data` | Yes (3x3 table) | PASS |
| `bragg_edges` zero below first edge | Yes | PASS |
| `bragg_edges` positive above first edge | Yes | PASS |
| `bragg_edges` 1/E behavior | Yes (rtol=1e-10) | PASS |
| `bragg_edges` step increase | Yes (5 edges) | PASS |
| `incoh_elastic_xs` E=0 | Yes | PASS |
| `incoh_elastic_xs` saturation | Yes (rtol=1e-6) | PASS |
| `incoh_elastic_xs` monotonicity | Yes (6 energies) | PASS |
| `incoh_elastic_xs` small-E linear | Yes (rtol=1e-4) | PASS |

### Missing test coverage
- **SCT fallback path:** No test exercises the short-collision-time approximation in `sab_kernel`.
- **lat==1 scaling:** No test with `lat=1` to verify `kT/tevz` scaling.
- **lasym==1:** No test for asymmetric S(alpha,beta) tables.
- **compute_thermal driver:** No test for the full `compute_thermal(pendf, ...)` path with PointwiseMaterial.
- **Bragg edge intensities:** Tests verify structure (edges exist, 1/E behavior, step increases) but do not validate absolute magnitudes against known reference values.
- **AD pass-through:** No test with ForwardDiff `Dual` numbers.

---

## Summary

| Check | Verdict |
|-------|---------|
| Free gas kernel vs sig()@200 | **PASS** |
| Free gas XS formula | **PASS** |
| SAB kernel vs sig()@2514 | **PASS** (bilinear vs quadratic interp noted) |
| SAB symmetrization | **PASS** |
| SCT fallback | **PASS** (missing sb2 second-atom term) |
| Bragg edges vs sigcoh() | **FAIL** -- missing `form()` structure factor |
| Incoherent elastic vs iel() | **PASS** |
| Energy grid | **PASS** (118 entries, exact match) |
| Detailed balance | **PASS** (exact by construction) |
| High-E / Low-E limits | **PASS** |
| AD compatibility | **PARTIAL** -- Float64 casts in kernels |
| Test coverage | **PASS** (adequate, some gaps noted) |

### Critical finding
**`build_bragg_data` omits the crystal structure factor `form(lat,l1,l2,l3)`** from Fortran's `sigcoh()`. This means Bragg edge weights/intensities will be physically incorrect for graphite, beryllium, and beryllium oxide. The existing tests do not catch this because they only check structural properties (edges exist, 1/E scaling, step behavior) rather than absolute magnitudes against reference data.

### Recommendations
1. Implement `form()` structure factor in `build_bragg_data` for at least graphite (lat=1).
2. Add a reference-value test for graphite Bragg edge XS at a known energy.
3. Add a test that exercises the SCT fallback path.
4. Remove `Float64()` casts in `free_gas_kernel` and `sab_kernel` for full ForwardDiff compatibility (use `float()` or let type promotion handle it).
