# Wave 3 Review: BROADR (Doppler Broadening)

**Reviewer:** Skeptical Reviewer (automated)
**Date:** 2026-03-21
**Files reviewed:**
- `src/processing/sigma1.jl` (229 lines)
- `src/processing/broadr.jl` (156 lines)
- `test/runtests.jl` (BROADR sections, lines ~1733-1953)

**Reference:** `njoy-reference/src/broadr.f90` (subroutines funky, hunky, hnabb, bsigma, broadn, thinb)

---

## Verdict: CONDITIONAL PASS

The implementation is solid and demonstrates careful translation of the SIGMA1 kernel. The F-functions, H-functions, and sigma1_at kernel are faithfully rendered. However, there are several issues ranging from a genuine convergence-test discrepancy to missing edge-case tests. The code earns a conditional pass: the two MUST-FIX items must be resolved before merging.

---

## 1. Diff Against Fortran

### 1.1 F-functions (funky) -- PASS

**sigma1.jl `f_func` vs broadr.f90 lines 1806-1839**

| Fortran | Julia | Match? |
|---------|-------|--------|
| `f(1) = half*erfc(a)` | `f0 = erfc(a)/2` | Yes |
| `expo = expo*resqpi; f(2) = half*expo` | `e1 = expo*rsqpi/2` | Yes |
| `expo = expo*a; f(3) = half*(f(1)+expo)` | `ae = a*expo*rsqpi; f2 = (f0 + ae)/2` | Yes |
| `expo = expo*a; f(4) = half*(two*f(2)+expo)` | `a2e = a*ae; f3 = (2*e1 + a2e)/2` | Yes |
| `expo = expo*a; f(5) = half*(three*f(3)+expo)` | `return (3*f2 + a*a2e)/2` | Yes |

The Fortran accumulates `expo` by multiplying by `a` at each step: `expo` goes through `exp(-a^2)*resqpi`, `a*exp(-a^2)*resqpi`, `a^2*exp(-a^2)*resqpi`, `a^3*exp(-a^2)*resqpi`. The Julia computes `ae = a*expo*rsqpi`, then `a2e = a*ae`, then `a*a2e`. These produce the same values. Verified correct.

The `alim = 10` guard is consistent between both implementations. The `_F_ZERO` constants match `fzero` in bsigma.

The `f_all` function correctly mirrors `f_func` in a single pass.

### 1.2 H-functions (hunky) -- PASS

**sigma1.jl `h_func`/`h_all` vs broadr.f90 lines 1762-1804**

| Property | Fortran | Julia | Match? |
|----------|---------|-------|--------|
| Small threshold | `small = 1.e-12` | `_H_SMALL_ABS = 1.0e-12` | Yes |
| Cancel threshold | `toler = 1.e-5` | `_H_CANCEL_TOL = 1.0e-5` | Yes |
| Cancellation test | `abs(h(k)) <= toler*abs(f(k))` | `abs(diff) < _H_CANCEL_TOL*afa` | Yes |
| Taylor fallback | `h(k) = hnabb(n, alast, aa)` | `h_taylor(k-1, a_old, a_new)` | Yes |
| s1/s2 computation | `s1 = (h(3)*oy+2*h(2))*oy+h(1)` | `s1 = (h[3]*oy+2*h[2])*oy+h[1]` | Yes |
| s2 computation | Matches term by term | `interval_contributions` | Yes |

Key design difference: Fortran hunky mutates module-level `h`, `f`, `s1`, `s2` arrays; Julia `h_all` returns tuples and `interval_contributions` is a separate pure function. This is a clean refactoring.

The `h_all` function correctly passes `a_old` (the previous `alast`) and `a_new` (the current `aa`) to `h_taylor` when cancellation is detected. The index mapping `k-1` in `h_taylor(k-1, ...)` correctly converts from 1-based tuple index to the 0-based n value.

### 1.3 Taylor Series (hnabb) -- PASS with ISSUE

**sigma1.jl `h_taylor` vs broadr.f90 lines 1841-1944**

The coefficient recursion logic is a line-by-line translation:

| Property | Fortran | Julia | Match? |
|----------|---------|-------|--------|
| pow2 constants | `pow2 = (/1.4142..., 2.0, 2.828..., 4.0, 5.656.../)` | `_POW2 = (sqrt(2.0), 2.0, ...)` | Yes |
| explim | `100.0` | `100` | Yes |
| rerr | `1.e-8` | `1e-8` | Yes |
| sign logic | `sign = -sign` if bb < aa; negate for odd n if bb < 0 | Same logic | Yes |
| k < 0 parity | `kk=mod(k,2); k=0; if(kk.ne.0) k=1` | `k = iseven(k) ? 0 : 1` | Yes |
| jalpha calc | `jalpha = (2*j+k-1-kstar)/2` | `jalpha = div(2j+k-1-kstar, 2)` | Yes |
| min iteration | `m < max(xn1, xn2)` | `m < max(T(n+1), abs(h_step*x))` | Yes |

**ISSUE (MUST-FIX): Convergence criterion mismatch.**

Fortran line 1856,1933:
```fortran
real(kr),parameter::aerr=1.e30_kr
test = aerr + rerr*abs(s)       ! test = 1e30 + 1e-8*|s|
if (abs(term).gt.test) go to 200  ! almost never true
```

Julia line 112:
```julia
if abs(term) <= rerr*abs(s)      ! purely relative, no absolute term
```

The Fortran's `aerr = 1e30` is intentional: it makes the absolute error criterion essentially vacuous, meaning the convergence test always passes as soon as `m >= max(n+1, |h*x|)`. The true convergence governor is the minimum-iteration check, not the term-size check.

The Julia version imposes a real relative convergence test (`abs(term) <= 1e-8*abs(s)`), which is STRICTER than the Fortran. This means:
1. The Julia code may iterate more than the Fortran for the same inputs.
2. When `s` is near zero (as can happen for near-cancelling integrals), the Julia code may never converge and exhaust all 50 iterations, returning an unconverged result.

In practice, for typical nuclear data this likely produces acceptable results (and may even be more accurate), but it is a semantic deviation from the Fortran. The mflag double-convergence check is correctly implemented in both.

**Recommendation:** Add `aerr = T(1e30)` to match Fortran exactly, or document the intentional deviation with justification. At minimum, add a test that exercises the Taylor series with near-zero `s` to confirm the 50-iteration fallback produces acceptable accuracy.

### 1.4 bsigma -> sigma1_at -- PASS with NOTE

**sigma1.jl `sigma1_at` vs broadr.f90 lines 1510-1675**

| Feature | Fortran | Julia | Match? |
|---------|---------|-------|--------|
| Velocity transform | `e(n) = sqrt(alpha*e(n))` (in broadn) | `v = map(e -> sqrt(alpha*e), seg_e)` | Yes |
| y computation | `y = en` (already velocity) | `y = sqrt(alpha*E)` | Yes |
| Panel search | `do while (en >= e(k+1)); k=k+1` | `searchsortedlast(v, y)` | Equivalent |
| f-function init | `f(ll) = fzero(ll)` | `f_cur = _F_ZERO` | Yes |
| Pass 1 loop direction | `l = k+1; do ll=klow,k; l=l-1` (backward) | `for l in k:-1:1` (backward) | Yes |
| Pass 1 aa | `aa = y - x` where x=e(l) | `aa = y - x_lo` | Yes |
| Pass 1 slope index | `s(i,l+1)*s1 + slope*s2` | `seg_xs[l+1]*s1 + slope*s2` | Yes |
| Pass 1 xx | `xx = xp*xp` (upper endpoint) | `xx_hi = x_hi*x_hi` | Yes |
| 1/v extrapolation | `sbt -= s(klow)*e(klow)*(oy^2*h(2)+oy*h(1))` | Same formula | Yes |
| Pass 2 oy sign | `oy = -oy` (becomes +1/y) | `oy_pos = 1.0/y` | Yes |
| Pass 2 slope index | `s(i,l)*s1 + slope*s2` | `seg_xs[l]*s1 + slope*s2` | Yes |
| Pass 2 xx | `xx = xm*xm` (lower endpoint) | `xx_lo_sq = x_lo*x_lo` | Yes |
| Constant extrapolation | `factor = (f(3)*oy+2*f(2))*oy+f(1)` | Same | Yes |
| Pass 3 (neg velocity) | Lines 1626-1660 | Lines 191-226 | See below |
| Output clamping | `if (sbt(i) < sigmin) sbt(i) = 0` (sigmin=1e-15) | `max(sbt, 0.0)` | See below |

**Pass 3 (negative velocity contribution):**

The Pass 3 translation is the most complex part. Tracing through:

Fortran (lines 1626-1660):
- `y = -y` (negate y); `aa = -y` (= original y)
- `call funky` to compute f at aa = y_orig
- `oy = -oy` (at this point oy was +1/y, so becomes -1/y)
- `aa = e(klow) - y` (y is now negative, so aa = e(klow) + y_orig)
- 1/v extrapolation: `sbt -= s(klow)*e(klow)*(oy^2*h(2)+oy*h(1))`
- Loop: `aa = x - y` = `x + y_orig`
- `sbt -= s(l)*s1 + slope*s2` (subtract contributions)
- Constant extrapolation: `sbt -= s(khigh)*factor`

Julia (lines 191-226):
- `aa_init = y` (= y_orig)
- `f_cur3 = f_all(aa_init)` (f at y_orig)
- `oy3 = -1.0/y` (= -1/y_orig)
- `y_neg = -y` (= -y_orig)
- `aa = v[1] - y_neg` (= v[1] + y_orig)
- 1/v: `sbt -= seg_xs[1]*v[1]*(oy3^2*h_vals3[2] + oy3*h_vals3[1])`
- Loop: `aa = x - y_neg` = `x + y_orig`
- `sbt -= seg_xs[l]*s1 + slope*s2`
- Constant: `sbt -= seg_xs[end]*factor`

This matches the Fortran trace faithfully. The variable naming with `y_neg` and `oy3` is clearer than the Fortran's mutation of `y` and `oy`.

**Output clamping difference:** Fortran uses `sigmin = 1e-15` (clamps values below 1e-15 to zero), Julia uses `max(sbt, 0.0)` (clamps negative to zero). This is a minor behavioral difference: Fortran preserves very small positive values (1e-16), Julia preserves them too; Fortran zeros values in [0, 1e-15], Julia does not. In practice this is inconsequential since 1e-15 barns is physically meaningless, but for exact reproducibility it differs.

### 1.5 broadn -> doppler_broaden -- PASS

The Fortran `broadn` implements its own adaptive refinement stack (lines 1256-1508). The Julia `doppler_broaden` delegates to `adaptive_reconstruct` from `adaptive_grid.jl`, which is the generic infrastructure reviewed in Wave 2.

This is a valid architectural decision. The key parameters map correctly:
- `tol` -> `errthn`
- `errmax` -> `errmax` (default `10*tol`)
- `errint` -> `errint` (default `tol/20000`)

The `_prepare_grid` and `_enrich_broadr_grid` functions add midpoints near sharp slope changes, which partially replicates the Fortran's node selection logic in broadn.

### 1.6 thinb -> thin_xs -- PASS with NOTE

| Property | Fortran | Julia | Match? |
|----------|---------|-------|--------|
| step_max | `stpmax = 1.24` | `step_max = 1.24` | Yes |
| Error criterion | `abs(sp-s(i,k)) > errthn*s(i,k)` | `abs(sp-xs[j,r]) > tol*abs(xs[j,r])` | Equivalent |
| Multi-reaction | Tests all reactions, union grid | `_compute_thin_mask` tests all columns | Yes |
| First/last point | Always kept | `keep[1]=true; keep[n]=true` | Yes |

**Note:** Fortran thinb also checks `e(k)*(1+eps) > thnmax` (line 1724), stopping thinning above the max broadening energy. The Julia thin_xs does not have this thnmax guard -- it thins the entire grid. This is acceptable because `doppler_broaden` handles the thnmax split in the `PointwiseMaterial` overload (lines 67-69 of broadr.jl), separating high-energy points before calling `doppler_broaden_multi`.

---

## 2. Numerical Stability

### 2.1 F-functions at extreme arguments -- PASS

- `a = 0`: Tested explicitly for all n=0..4 (lines 1738-1746 of runtests.jl). Values match analytical.
- `a` very large (15.0): Tested to return 0.0 exactly (line 1761). The `alim = 10` guard prevents overflow in `exp(-a^2)`.
- Intermediate values (0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0): Tested against analytical formulas.

### 2.2 H-functions cancellation -- PASS

- `a = b`: Returns exactly 0.0 (line 1779).
- Well-separated a, b: Direct subtraction matches (line 1784-1785).
- Near-cancellation (a=1.5, b=1.5+1e-8): Taylor series produces finite, physically meaningful result (lines 1790-1804). The derivative check for n=0 gives confidence the Taylor series is working.
- h_taylor direct verification (lines 1807-1812): Matches h_func at rtol=1e-6 for small intervals.

### 2.3 sigma1_at numerical behavior -- PASS

- Constant XS at various energies returns the constant within 1% (lines 1940-1951).
- The `@goto` early exit when `aa > atop` prevents wasted computation without introducing numerical issues.
- The clamping `max(sbt, 0.0)` prevents negative broadened cross sections.

---

## 3. Edge Cases

### 3.1 Single-point cross section -- CONCERN (SHOULD-FIX)

The `doppler_broaden` function requires `n >= 2` (line 19: `@assert n >= 2`). Single-point cross sections will throw an assertion error. This is documented behavior but there is no test for it.

### 3.2 Two-point cross section -- PASS

With n=2, `sigma1_at` will have a single segment. The code handles this correctly: the loop `for l in k:-1:1` with k=1 iterates once (or not at all), and boundary extrapolations cover the rest. No explicit test exists but the constant-XS test with 4 points exercises similar logic.

### 3.3 Temperature T=0 -- CONCERN (MUST-FIX)

`doppler_broaden` asserts `T > 0` (line 19), which is correct since T=0 makes `alpha -> inf`. However, the `PointwiseMaterial` overload computes `T_eff = T_new - T_old` and errors if `T_eff <= 0`. This means:
- T_new = T_old (no broadening needed): throws error instead of returning input unchanged.
- T_new < T_old: throws error (correct, can't un-broaden).

**The T_new == T_old case should return the input unchanged rather than throwing an error.** This is a usability bug.

### 3.4 Very high temperature -- NOT TESTED

No test exercises extreme broadening (T > 10000K). At very high T, `alpha` becomes small, `y = sqrt(alpha*E)` becomes small, and `1/y` becomes large. The `atop = 4.0` cutoff means the integration window expands. No overflow risk is apparent (the `alim = 10` guard in f_func handles large arguments), but a test would be prudent.

### 3.5 Negative cross sections -- NOT TESTED

MLBW interference terms can produce negative cross sections. The `max(sbt, 0.0)` clamping in sigma1_at will zero them out. This may or may not be the desired behavior for interference-broadened data. No test covers this case.

---

## 4. AD Compatibility

### 4.1 Mutation in sigma1.jl -- CONCERN

`h_taylor` allocates and mutates two local arrays: `cm = zeros(T, 50)` and `cmstar = zeros(T, 50)` (line 90). These are local temporaries, so they do not violate purity in the mathematical sense (the function has no side effects), but:

1. **Array mutation blocks reverse-mode AD.** If ForwardDiff is the only AD target, `zeros` + mutation works. For Zygote or Enzyme in reverse mode, the `cm[j] = ...` mutations inside the loop will fail or produce incorrect gradients.
2. **Heap allocation per call.** Every call to `h_taylor` allocates two 50-element arrays. In a tight broadening loop, this generates GC pressure.

The rest of sigma1.jl is mutation-free and AD-clean: `f_func`, `f_all`, `h_func`, `h_all`, `interval_contributions`, and `sigma1_at` use only local scalars and tuples (with the exception of `sigma1_at`'s `v = map(...)` which allocates a vector, but is not mutated after construction).

**Recommendation:** Consider using `SVector` or tuple-based accumulation in `h_taylor` to avoid mutation and allocation, or at minimum document that `h_taylor` is ForwardDiff-only compatible.

### 4.2 broadr.jl -- MINOR CONCERN

- `_enrich_broadr_grid` uses `push!` and `sort!`/`unique!` (lines 150-154). These are called only during grid setup, not in the hot evaluation path. AD would not typically differentiate through grid construction.
- `doppler_broaden_multi` uses `sort!` and `unique!` on the initial grid (line 49). Same note applies.
- `_compute_thin_mask` mutates `keep` (a `BitVector`). Not in the AD-sensitive path.

### 4.3 No try-catch -- PASS

No try-catch blocks in either file.

---

## 5. Physical Invariants

### 5.1 Constant XS preservation -- PASS

Tested at line 1826-1842. A constant sigma = 10.0 returns within 10% after broadening. The tolerance is generous (rtol=0.1) but appropriate given boundary effects on a finite energy grid.

### 5.2 1/v XS invariance -- PASS

Tested at line 1844-1863. The product sigma(E)*sqrt(E) is checked to be constant (within 5%) for energies away from boundaries. This is the correct invariant: if sigma(E) = C/sqrt(E), then sigma_broadened(E)*sqrt(E) should also equal C.

### 5.3 Step smoothing -- PASS

Tested at line 1865-1882. A step function is confirmed to be smoothed (all values positive, plateau values approached far from the step).

### 5.4 Integral preservation (Barker's theorem) -- NOT TESTED

Doppler broadening preserves the integral of sigma(E)*dE over all energies. No test verifies this. This would be a strong validation of the kernel correctness.

---

## 6. Test Coverage

### 6.1 F-function orders n=0..4 -- PASS

All five orders tested individually at a=0 (lines 1738-1746) and a>0 (lines 1750-1757 for n=0,1). Large-a tested for all orders (line 1760-1762). f_all consistency checked for all orders at multiple a values (lines 1764-1770).

**Gap:** Only n=0 and n=1 are tested against analytical formulas at non-zero a. The recursion for n=2,3,4 is only verified indirectly through f_all consistency. Adding explicit analytical checks for f_2, f_3, f_4 at specific a values would strengthen coverage.

### 6.2 H-function cancellation threshold -- PASS

The threshold is tested implicitly: the near-cancellation test (eps_val = 1e-8) forces the cancellation path. The derivative check for n=0 confirms the Taylor series produces physically correct results. However, no test explicitly verifies that the threshold of 1e-5 triggers the Taylor path (vs direct subtraction).

### 6.3 Broadening against analytical results -- PARTIAL

The 1/v and constant-XS invariant tests provide analytical validation. However, no test compares against a known exact broadened cross section (e.g., a single Breit-Wigner resonance with known Doppler-broadened psi/chi functions). This would be the gold standard validation.

### 6.4 Missing tests -- SHOULD ADD

1. **T=0 or T_new == T_old handling** -- no test.
2. **Two-point (single segment) cross section** -- no explicit test.
3. **Very high temperature** -- no test.
4. **Negative input cross sections** -- no test.
5. **Integral preservation** -- no test.
6. **sigma1_at at very low energy (< 0.001 eV)** -- boundary extrapolation regime not explicitly tested.
7. **sigma1_at with a single-element segment** -- not tested.

---

## Summary of Findings

### MUST-FIX (2 items)

1. **h_taylor convergence criterion mismatch** (Section 1.3): The Fortran uses `aerr = 1e30` in the convergence test, making it effectively a minimum-iteration-only check. The Julia uses a strict relative test `abs(term) <= 1e-8*abs(s)` which may fail to converge when `s` is near zero. Either match the Fortran exactly or add explicit documentation and a regression test.

2. **T_new == T_old should not error** (Section 3.3): The `PointwiseMaterial` overload of `doppler_broaden` throws an error when T_new == T_old. It should return the input unchanged.

### SHOULD-FIX (3 items)

3. **h_taylor allocates two arrays per call** (Section 4.1): In a tight loop, this creates GC pressure. Consider stack-allocated alternatives (MVector or NTuple with manual indexing).

4. **Output clamping uses 0.0 instead of Fortran's 1e-15** (Section 1.4): Minor behavioral difference. Consider using `sigmin = 1e-15` to match Fortran exactly, or document the intentional deviation.

5. **Add analytical F-function tests for n=2,3,4** (Section 6.1): Currently only n=0 and n=1 have explicit analytical verification at non-zero arguments.

### NICE-TO-HAVE (4 items)

6. Add integral preservation test (Barker's theorem).
7. Add single Breit-Wigner resonance broadening comparison test.
8. Add extreme temperature (T > 10000K) test.
9. Add negative cross section input test.

---

## Code Quality Notes

- The separation of sigma1.jl (pure math kernel) from broadr.jl (pipeline orchestration) is clean and well-motivated.
- The correspondence comments at the top of both files mapping Julia functions to Fortran subroutines are excellent for maintainability.
- The use of `@goto`/`@label` for early exit mirrors the Fortran `go to` structure faithfully without introducing control-flow complexity.
- The `interval_contributions` function is correctly marked `@inline`.
- The `_prepare_grid` / `_enrich_broadr_grid` enrichment is a thoughtful addition that is not present in the Fortran (which relies on its paging scheme to have sufficient surrounding data).
