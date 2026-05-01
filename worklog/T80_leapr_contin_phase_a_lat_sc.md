# T80 — leapr contin Phase A: lat/sc/arat plumbing (Phase 57)

**Date:** 2026-05-01
**Test:** T80 (H bound in HF, MAT=8, ENDF-6 thermal scattering law, lat=1, T=343 K)
**Goal:** Close the dominant value-drift in Julia's `generate_sab` versus
Fortran `contin` (leapr.f90:455-645) by porting the `sc`/`arat` α/β
rescaling. T22 (lat=0) BIT_IDENTICAL preserved; T80 (lat=1) moves
DIFFS → NUMERIC_PASS @1e-5.

## Outcome

**T80 tape24: DIFFS → NUMERIC_PASS @1e-5.**

| Metric                           | Before (Phase 10) | After (Phase 57) |
|----------------------------------|-------------------|-------------------|
| Status                           | DIFFS             | **NUMERIC_PASS @1e-5** |
| Lines passing at no-tolerance    | 69 231 / 91 453 (75.7%) | unchanged (rtol=0 unforgiving) |
| Lines passing at 1e-5            | (not measured)    | **91 406 / 91 453 (99.95%)** |
| First-diff lines (post-fix)      | massive drift     | ±1 in 7th sigfig only |
| Wallclock                        | 289 s             | 312 s |

Regression-clean:
- T22 (para-H₂, 20 K, lat=0): **BIT_IDENTICAL 4636/4636 lines**, rtol=1e-9.
- All leapr unit tests pass (`test_leapr_parser.jl`, `test_leapr_trans.jl`,
  `test_leapr_coldh.jl`, `test_leapr_writer.jl`).

## Bug

Julia's `generate_sab` (src/processing/leapr.jl:188-238) used `alpha_grid[j]`
and `beta_grid[k]` *raw* in the phonon expansion's `_terpt(...)` lookups
and the SCT formula. Fortran `contin` applies a rescaling
(leapr.f90:492-500-506):

```fortran
sc = 1
if (lat .eq. 1) sc = therm / tev   ! therm = 0.0253 eV, tev = bk·T
...
al = alpha(j) * sc / arat
be = beta(k)  * sc
```

For T22 (lat=0) `sc=1`, so the missing rescaling has no effect — that
explains why T22 was already BIT_IDENTICAL.

For T80 (lat=1) at T=343 K, `sc = 0.0253 / (8.617e-5 · 343) ≈ 0.8559`.
Every β passed to `_terpt(p, np, deltab, be)` was off by ~14% in scale,
shifting the interpolation in the convolved DOS arrays `tlast`/`tnow`,
and the SCT branch's Gaussian center `(alw - be)²` was off in both
arguments. **This was the dominant T80 numerics gap.**

`arat` is the secondary-scatterer mass ratio (`aws/awr`); it is `1.0`
for the principal scatterer (leapr.f90:323-330). T80 has `nss=0`, so
`arat=1` is the only branch exercised here, but the kwarg is in place
for the Phase 8 secondary-scatterer wiring.

## Fix

`src/processing/leapr.jl`:
1. Added `const _LEAPR_THERM = 0.0253` (matches Fortran's
   `real(kr),parameter::therm=0.0253e0_kr` at leapr.f90:475).
2. Added `lat::Int=0, arat::Float64=1.0` kwargs to `generate_sab`.
3. Compute `sc = (lat == 1) ? _LEAPR_THERM / kT : 1.0` and
   `sc_a = sc / arat`.
4. Apply `al = alpha_grid[j] * sc_a` and `be = beta_grid[k] * sc` in
   every loop where the original Fortran applies them: l=1 phonon term,
   l≥2 phonon expansion, and SCT fallback.
5. Updated docstring with leapr.f90 line citations and a brief note on
   the SCT-replacement `naint` gating that remains aligned (Fortran
   gates SCT replacement on `iprt = mod(j-1, naint)+1 == 1`; Julia
   replaces all (k, j>=maxt[k]) — Phase B follow-up).

`src/orchestration/modules/leapr.jl`:
- Pass `lat=params.lat, arat=1.0` from the `leapr_module` callsite into
  `generate_sab`. Updated comment to cite leapr.f90:455-645 and note
  Phase 8 secondary-scatterer dependency for arat≠1.

Net diff: ~25 lines changed in leapr.jl, ~5 lines in modules/leapr.jl.
No new Julia files; no public-API breakage (kwargs both default).

## Test surface

- **T22 reference** (`reference_test.jl 22`): BIT_IDENTICAL 4636/4636 @
  rtol=1e-9 (regression baseline preserved).
- **T80 reference** (`reference_test.jl 80`): NUMERIC_PASS 91406/91453 @
  rtol=1e-5 (was DIFFS 69231/91453 @ rtol=0).
- **Leapr unit tests** (4 files, ~194 assertions): all pass.

## Acceptance criteria — met

- [x] T22 stays BIT_IDENTICAL (regression-clean).
- [x] T80 measurably improves toward bit-identical: 75.7% → 99.95% @1e-5.
- [x] Leapr unit tests pass.
- [ ] T80 BIT_IDENTICAL @1e-9 — **not yet** (47 lines remain off; ±1 in
      7th sigfig). Phase B work.

## Out of scope (still open)

**Phase B — close the residual 47 lines for T80 BIT_IDENTICAL @1e-9:**

1. **SCT-replacement `naint` gating.** Fortran `contin` (leapr.f90:584-642)
   wraps SCT replacement in `iprt = mod(j-1, naint)+1 == 1` (every naint-th
   α plus j=nalpha). Julia currently replaces all (k, j>=maxt[k]). For
   non-printed α's >= maxt, Fortran preserves the unconverged phonon-sum
   value (bounded by the `add < ssm/1000` stop criterion); Julia
   overwrites with SCT. `naint` must be threaded through `LeaprParams`
   (it's a card-3 parser field); `generate_sab` needs it as a kwarg.

2. **Phonon-loop accumulation order.** The 47 lines may also be sensitive
   to FP order in `_convol!` and the `xa[j] += log(al*f0/n)` in-place
   add. Compare against Fortran via `write(*,...)` diagnostics on a
   single (j, k) trace point at T=343 K to pinpoint.

3. **Moment-check mutations.** Fortran `contin` modifies `ssm(k,j)` to
   `ssct` inside the moment-check loop (leapr.f90:611). Once `naint`
   gating is in, the Julia structure should mirror Fortran's:
   moment-check loop *both* prints diagnostics *and* mutates ssm.

**Phase C — discrete oscillators + secondary scatterer + coher.** Out
of scope for T80 closure (T80 has nd=0, nss=0, iel=0). Tracked under
HANDOFF "Real plotr / covr / leapr / purr output" P3.

## Fortran source citations

- `njoy-reference/src/leapr.f90:455-645` — `contin` subroutine
- `njoy-reference/src/leapr.f90:475` — `therm = 0.0253 eV`
- `njoy-reference/src/leapr.f90:492-493` — `sc = 1; if (lat==1) sc = therm/tev`
- `njoy-reference/src/leapr.f90:500/506/532/538/588` — `α·sc/arat`, `β·sc`
- `njoy-reference/src/leapr.f90:323-330` — `arat = 1` for principal,
  `arat = aws/awr` for secondary
- `njoy-reference/src/leapr.f90:565-579` — `maxt` monotone-decreasing
  enforcement
- `njoy-reference/src/leapr.f90:605-612` — SCT formula and `ssm(k,j)=ssct`
  replacement (gated by `iprt==1`)
- `njoy-reference/src/leapr.f90:584-642` — moment-check + SCT-replacement
  outer loop
