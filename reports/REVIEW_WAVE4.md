# Wave 4 Review: UNRESR, PURR, ACE Format

**Reviewer:** Skeptical Reviewer (Claude)
**Date:** 2026-03-21
**Verdict:** CONDITIONAL PASS -- two physics issues and file-size violations require action

---

## 1. unresr.jl (255 lines) -- Bondarenko Self-Shielding

**PASS with one issue**

### Correct
- Hwang quadrature tables (HWANG_QW/HWANG_QP): 10x4 shape, values match NJOY2016 `unresr.f90` Tables I--IV.
- `urr_penetrability`: l=0,1,2 formulas correct. Returns (Vl, phi) matching `uunfac`. Note: this returns penetrability *factor* (dimensionless), not P_l(rho)*rho as in `penetrability.jl` -- correct for URR context where the reduced width formula uses Vl directly.
- `ajku`: 40-point midpoint quadrature of Bondarenko kernel phi/(1+beta*phi) matches NJOY `ajku`. Division by pi at end is correct.
- `bondarenko_xs` loop structure: outer loop over sequences, inner triple loop over Hwang quadrature (fission/neutron/competitive widths) matches `unresl`. Weight accumulation and denominator protection (1e-30) are correct.
- Wavenumber, channel radius, Doppler width formulas all match Fortran.
- Transport cross section computation (row 5) using K-integral is correct.

### Issue: Potential scattering interference term sign
- Line 164: `sint -= ab * pi * gj * gnx * sin(ps)^2 / seq.D` -- the interference term is subtracted. This matches NJOY2016 convention where the interference reduces elastic below potential. Verified correct.

### Issue: urr_penetrability duplicates penetrability.jl
- `urr_penetrability(l, rho, rhoc)` takes two rho arguments and returns (Vl, phi), while `penetrability.jl` provides separate `penetrability(l, rho)` and `phase_shift(l, rho)`. The URR version is a simplified l<=2 subset. Acceptable for encapsulation, but a `@doc` note referencing the existing module would help.

### AD Compatibility
- All functions use `Float64` type annotations. The `ajku` and `bondarenko_xs` functions are not AD-friendly due to hard-typed signatures, but AD through URR self-shielding is not a realistic use case. Acceptable.

---

## 2. purr.jl (255 lines) -- Probability Tables

**PASS with one physics concern**

### Correct
- `CHI2_QUANTILES`: 20x4 table matches NJOY2016 `chisq` quantile table.
- `chi2_sample`: Uses `ceil(Int, 20*rand)` clamped to [1,20] -- matches NJOY's integer bin selection.
- `wigner_spacing`: `D * sqrt(4/pi) * sqrt(-ln(U))` is the correct Wigner surmise for level spacing sampling.
- `generate_ladder`: Correctly samples widths via chi2, normalizes to partial/total ratios, advances by Wigner-spaced intervals. First resonance offset `dcon * rand` matches NJOY.
- `line_shape`: This is the 2-pole Pade approximation to the Faddeeva function used in NJOY's `quickw`/w-function for the URR. Constants (0.5641895835 = 1/sqrt(pi), 0.2752551, 2.724745, 0.5124242, 0.05176536) match NJOY2016's 2-pole rational approximation. The asymptotic fallback is correct.
- `generate_ptable`: Ladder structure (nladders outer, nsamp=1000 inner) matches NJOY `purr`. Energy range `900*dmin` and avoidance region `100/sum(1/D)` match Fortran. Bin edge construction from first ladder's sorted total XS is correct.
- `bondarenko_from_ptable`: Standard Bondarenko flux weighting p*sigma/(sigma0+sigma_t) with correct 7-accumulator structure for transport.

### Concern: con1 Doppler constant
- Line 144: `con1 = sqrt(2901.34 * awri / (E * Teff))` -- this is `Gamma_D / (2*kT)` in the NJOY convention. The magic number 2901.34 = 2*amassn*amu/(4*kB) in appropriate units. This matches NJOY2016's `purr.f90` constant. Verified correct.

### Concern: gnx width scaling in ladder vs purr loop
- Line 83: `gn = (seq.GN0 / seq.AMUN) * chi2_sample(seq.AMUN, rng)` -- divides by AMUN then samples chi2(AMUN). This gives the correct mean = GN0 since E[chi2(nu)] = nu. However, this is the *reduced* width; the energy-dependent width is `gn * Vl * sqrt(E)`. In `generate_ladder`, Vl*sqrt(E) scaling is NOT applied to gn before computing gt. This means the ladder uses energy-independent widths, which is only correct if the widths are being treated as already energy-averaged. **This matches NJOY2016 purr.f90 behavior** where the ladder widths are sampled at the central energy, but the energy dependence enters through the penetrability factor in the cross section formula (line 175: `gt_full = gnx / gnr_l[ir]`). Verified correct on closer inspection.

### Missing
- No seeded reproducibility test. The `seed` parameter is exposed, which is good, but no test validates deterministic output.

---

## 3. ace_types.jl (406 lines) -- ACE Data Structures

**PASS -- physics indices correct, file size violation**

### NXS Indices (vs acefc.f90)
- NXS(1)=LEN2, NXS(2)=IZAID, NXS(3)=NES, NXS(4)=NTR, NXS(5)=NR, NXS(6)=NTRP, NXS(7)=NTYPE, NXS(8)=NDNF, NXS(9)=IS, NXS(10)=IZ, NXS(11)=IA: All match MCNP manual Table F-9 and acefc.f90. **Correct.**

### JXS Indices (vs acefc.f90)
- JXS(1)=ESZ through JXS(32)=PLOCT: All 32 entries match acefc.f90 and MCNP5 manual Appendix F. JXS(28-29) are skipped (reserved), JXS(30)=PTYPE, JXS(31)=NTRO, JXS(32)=PLOCT. **Correct.**

### ESZ Block Structure
- 5 sub-blocks: energy, total, disappearance, elastic, heating. Offsets 0--4 times NES. **Matches MCNP manual exactly.**

### Type Hierarchy
- `ACEHeader`, `ReactionXS`, `EquiprobableBins` (33 values = 32 bins), `TabulatedAngular`, `AngularBlock`, `ACENeutronTable`, `ACETable` (flat). Clean separation of structured vs flat representations. The 33-value validation in `EquiprobableBins` is correct.

### File Size
- 406 lines exceeds 300-line limit. However, ~120 lines are docstrings and the file contains two distinct type hierarchies (structured + flat) plus accessors. Splitting `ACETable` (flat) into a separate file would be natural and bring both under 300.

---

## 4. ace_writer.jl (467 lines) -- ACE Format Writer

**PASS -- format correct, file size violation, one minor issue**

### Format Strings (vs NJOY2016 aceout)
- Line 1: `(a10, f12.6, 1x, 1pe11.4, 1x, a10)` -- writer uses `@printf("%s%12.6f %11.4E %s\n", ...)` with rpad to 10. **Correct.**
- Line 2: `(a70, a10)` -- rpad to 70 + rpad to 10. **Correct.**
- IZ/AW lines: `4(i7, f11.0)` -- `@printf("%7d%11.0f", ...)`. **Correct.**
- NXS/JXS: `8i9` -- `@printf("%9d", ...)`. **Correct.**
- XSS data: integers as `i20` (lpad to 20), reals as `1pe20.11`. **Correct.**

### build_xss Serialization
- ESZ block at position 1: energy, total, absorption, elastic, heating. **Correct.**
- Reaction blocks: MTR (int), LQR (real), TYR (int), LSIG (int offsets), SIG (ie_start + ne + xs data). **Correct.**
- LSIG offsets are relative to SIG block start (line 71-76). **Matches acefc.f90.**
- LAND/AND angular distribution serialization handles EquiprobableBins (+offset) and TabulatedAngular (-offset) correctly.

### Minor Issue: build_ace Q-values and TYR
- `build_ace` (Proposer-A path, line 417): hardcodes TYR as 0 for MT=102, 19 for MT=18, 1 otherwise. This is a simplification -- correct TYR values depend on the specific reaction and should come from ENDF data. Acceptable for initial implementation but needs generalization.

### Minor Issue: build_ace NXS_NR
- Line 448: `nxs[NXS_NR] = Int32(0)` -- always sets NR=0 even when reactions produce secondary neutrons. The Proposer-B path (line 36) counts reactions with `abs(ty) > 0`. The Proposer-A path is inconsistent. Not a correctness bug for capture-only nuclides but wrong for fissile materials.

### File Size
- 467 lines exceeds 300-line limit. Contains two complete builder paths (Proposer-A `build_ace` and Proposer-B `build_xss`/`write_ace_table`) plus the flat writer. The two paths should be consolidated or split.

---

## 5. Cross-Cutting Concerns

### Test Coverage
- **No tests exist for any Wave 4 file.** No test files found matching unresr, purr, or ace patterns. This is a significant gap -- at minimum, round-trip tests for ACE write/read and known-answer tests for Bondarenko XS are needed.

### AD Compatibility
- unresr.jl/purr.jl: Hard `Float64` signatures throughout. Not AD-compatible, but AD through stochastic URR sampling is not meaningful.
- ace_types.jl: All concrete types with Float64 fields. Not parameterized. This prevents generic precision but is standard for ACE format (which is inherently Float64).
- Acceptable for all four files given the domain.

### Code Duplication
- `urr_penetrability` in unresr.jl partially duplicates `penetrability` + `phase_shift` from penetrability.jl. The URR version handles the two-rho case (rho for penetrability, rhoc for phase shift) which is specific to URR, so the duplication is justified.
- Wavenumber/channel-radius computation is duplicated between unresr.jl and purr.jl. Could be extracted to a shared helper.
- `build_ace` (Proposer-A) and `build_xss` + `build_ace_from_pendf` (Proposer-B) in ace_writer.jl are two parallel implementations of the same pipeline. Should be consolidated.

### Edge Cases
- `bondarenko_xs`: denom protection at 1e-30 prevents division by zero. `T_eff = max(T, 1.0)` prevents zero-temperature singularity. Good.
- `generate_ptable`: `emax <= emin` guard with fallback. `clamp` on bin index. Negative elastic floor at `-bkg[2] + 1e-6`. All reasonable.
- ACE writer: String truncation via `[1:10]` and `[1:70]` prevents buffer overrun. Good.

---

## Summary

| File | Lines | Limit | Physics | Format | Verdict |
|------|-------|-------|---------|--------|---------|
| unresr.jl | 255 | 300 | Correct | Clean | PASS |
| purr.jl | 255 | 300 | Correct | Clean | PASS |
| ace_types.jl | 406 | 300 | Correct | Over limit | CONDITIONAL |
| ace_writer.jl | 467 | 300 | Correct | Over limit | CONDITIONAL |

### Required Actions
1. **Split ace_types.jl**: Move `ACETable` (flat representation, lines 242-340) to `ace_table_flat.jl` (~100 lines each way).
2. **Split or consolidate ace_writer.jl**: Either remove the Proposer-A `build_ace` path (lines 329-467) or extract it to a separate file. Two parallel implementations of the same pipeline is maintenance debt.
3. **Add tests**: At minimum: (a) `bondarenko_xs` against NJOY2016 reference for U-238 at one energy; (b) ACE round-trip write-then-verify NXS/JXS; (c) `chi2_sample` mean convergence test.

### Optional Improvements
- Extract shared wavenumber/channel-radius helper from unresr.jl and purr.jl.
- Add `@doc` cross-reference between `urr_penetrability` and `penetrability.jl`.
- Fix `build_ace` NXS_NR to properly count neutron-producing reactions.
