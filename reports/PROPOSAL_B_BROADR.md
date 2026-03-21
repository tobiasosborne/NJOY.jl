# Proposal B: BROADR Implementation -- Numerical Robustness + AD Compatibility

## Design Philosophy

This proposal separates the SIGMA1 kernel (pure math) from the broadening pipeline
(orchestration), producing two files that are each under 300 lines:

1. **`sigma1.jl` (229 lines)** -- Pure mathematical functions with no mutation
2. **`broadr.jl` (156 lines)** -- Pipeline orchestration reusing `adaptive_reconstruct`

Key design principles:
- **No mutation in hot paths**: Every function returns values, never modifies arguments
- **AD-compatible**: No `try/catch`, no mutation, no global state in the kernel
- **Clean reuse**: `adaptive_reconstruct` from `adaptive_grid.jl` handles all bisection
  logic -- we do NOT reimplement the broadn loop
- **Cancellation avoidance**: Taylor series fallback in h_func when |f_n(a)-f_n(b)|/|f_n(a)| < 1e-5

## Architecture

```
User API
  doppler_broaden(energies, xs, T, awr)     [broadr.jl]
    |
    v
  Build closure: E -> sigma1_at(E, seg_e, seg_xs, alpha)    [sigma1.jl]
    |
    v
  adaptive_reconstruct(closure, grid, config)    [adaptive_grid.jl -- REUSED]
    |
    v
  thin_xs(out_e, out_xs)    [broadr.jl]
    |
    v
  (broadened_energies, broadened_xs)
```

## File: sigma1.jl -- SIGMA1 Kernel Functions

### Functions

| Function | NJOY Equivalent | Purpose |
|----------|----------------|---------|
| `f_func(n, a)` | `funky` | F-function: integral(a..inf) z^n exp(-z^2) dz / sqrt(pi), n=0..4 |
| `f_all(a)` | `funky` (all 5) | Return (f_0..f_4) in one exp+erfc call |
| `h_func(n, a, b)` | `hunky` (one) | H-function with cancellation detection |
| `h_all(a_old, f_old, a_new)` | `hunky` (all 5) | Stateless replacement for NJOY's module-level arrays |
| `h_taylor(n, aa, bb)` | `hnabb` | Taylor series for small |b-a| |
| `interval_contributions(h, oy, yy, xx)` | `hunky` s1/s2 | Compute integral terms from one segment |
| `sigma1_at(E, seg_e, seg_xs, alpha)` | `bsigma` | Full broadened XS at a single energy |

### Numerical Robustness

The F-functions use the exact recursion from NJOY:
```
f_0(a) = erfc(a)/2
f_1(a) = exp(-a^2)/(2*sqrt(pi))
f_{n+2}(a) = ((n+1)*f_n(a) + a^{n+1}*exp(-a^2)/sqrt(pi)) / 2
```

The H-function avoidance strategy (matching NJOY's hunky):
1. Compute `diff = f_n(a) - f_n(b)` directly
2. If `|diff| < 1e-12 * |f_n(a)|`: return zero (below noise floor)
3. If `|diff| < 1e-5 * |f_n(a)|`: switch to `h_taylor` (Taylor series)
4. Otherwise: use direct subtraction

The Taylor series (`h_taylor`) is a faithful translation of NJOY's `hnabb` with
coefficient recursion, converging to 1e-8 relative precision. It requires two
consecutive converged terms to declare convergence (the `mflag` protocol).

### Two-Pass + Pass 3 Structure (bsigma)

The `sigma1_at` function implements NJOY's bsigma with three passes:

1. **Pass 1 (below)**: Loop from the panel containing E downward. Each segment
   contributes `sigma_hi * s1 + slope * s2` where s1/s2 come from the H-functions.
   Extends with 1/v extrapolation to E=0.

2. **Pass 2 (above)**: Loop from the panel containing E upward. Each segment
   contributes `sigma_lo * s1 + slope * s2`. Extends with constant extrapolation
   to infinity.

3. **Pass 3 (negative velocity)**: Accounts for target nuclei moving opposite
   to the neutron. This is the Fortran lines 1626-1660 (y=-y section) and is
   important at low energies where E is comparable to the Doppler width.

Each pass early-exits when `|a| > 4` (the Gaussian tail is negligible).

### AD Compatibility

All functions are pure:
- `f_func`, `f_all`: Take a scalar, return a scalar or tuple
- `h_func`: Takes two scalars, returns a scalar
- `h_all`: Takes (scalar, tuple, scalar), returns (tuple, tuple) -- no mutation
- `sigma1_at`: Takes vectors by reference but does not modify them;
  creates the velocity-space `v` array via `map` (allocation, not mutation)

The only imperative construct is `@goto` for early exit from the below/above
loops, which has no side effects.

## File: broadr.jl -- Pipeline Orchestration

### Functions

| Function | Purpose |
|----------|---------|
| `doppler_broaden(energies, xs, T, awr)` | Single-reaction broadening |
| `doppler_broaden_multi(energies, xs_matrix, T, awr)` | Multi-reaction broadening |
| `doppler_broaden(pendf::PointwiseMaterial, ...)` | PointwiseMaterial wrapper |
| `thin_xs(energies, xs)` | Thinning (vector or matrix) |

### Key Design Decision: Reuse adaptive_reconstruct

Instead of reimplementing NJOY's `broadn` bisection loop (which Proposal A does),
we wrap `sigma1_at` in a closure and pass it to the existing `adaptive_reconstruct`.
This means:

- **Zero duplicated adaptive logic** -- convergence tests, stack management,
  thermal tightening, step-ratio guards are all reused from Wave 1/2
- **Automatically benefits from future improvements** to the adaptive algorithm
- **Trivially composable**: the same pattern works for THERMR or any other module
  that needs adaptive linearization of a computed function

The closure pattern:
```julia
broadened_eval = E -> (sigma1_at(E, seg_e, seg_xs, alpha),)
out_e, out_v = adaptive_reconstruct(broadened_eval, initial_grid, config)
```

### Thinning

The `thin_xs` function matches NJOY's `thinb`:
- Tests all reactions simultaneously (union grid preservation)
- Step ratio guard (stpmax=1.24) prevents skipping sharp features
- First and last points always preserved

The implementation uses a `_compute_thin_mask` function that returns a BitVector,
keeping the thin logic separate from the indexing logic.

### Grid Enrichment

The `_enrich_broadr_grid` helper adds midpoints near sharp slope changes in the
input cross section, giving the adaptive algorithm better starting points near
resonance peaks.

## Tests Added

All tests pass (6450 total, 0 failures).

### F-function tests
- Analytical values at a=0 for all n=0..4
- Comparison with erfc(a)/2 for f_0
- Comparison with exp(-a^2)/(2*sqrt(pi)) for f_1
- Large-a limit (all functions -> 0 at a=15)
- f_all consistency with individual f_func calls

### H-function cancellation tests
- h(n, a, a) = 0 exactly
- h(n, a, b) matches f(a)-f(b) for well-separated arguments
- h(n, a, a+1e-8) is finite and matches expected derivative value (CRITICAL)
- h_taylor matches h_func for small intervals
- h_all consistency

### Broadening invariant tests
- Constant XS remains constant (Fe-56, T=300K, energies 1-1000 eV)
- 1/v XS: sigma*sqrt(E) remains constant (Westcott theorem)
- Step function is smoothed (all values remain positive)

### Thinning tests
- Constant function thins dramatically (500 -> <50 points)
- Linear ramp thins well (<100 points from 200)
- Interpolation accuracy preserved within tolerance

### Kernel test
- sigma1_at returns correct value for constant cross section

## Comparison with Proposal A (anticipated)

| Aspect | Proposal B | Proposal A (anticipated) |
|--------|-----------|-------------------------|
| Kernel purity | Pure functions, no mutation | May use mutable state |
| AD compatibility | Full (no try/catch, no mutation) | May have barriers |
| Adaptive loop | Reuses adaptive_grid.jl | Custom bisection (broadn translation) |
| Code volume | 385 lines total | Likely larger |
| Taylor series | Direct hnabb translation | Unknown |
| Testing | Invariant-based + analytical | Unknown |
