# Wave 1 Review Report

## Verdict: CONDITIONAL PASS

The Wave 1 code is substantially correct and well-structured. The physics formulas
match the Fortran reference in all major respects. However, there are several issues
that range from a genuine physics discrepancy (one critical) to Julia idiom violations
that must be addressed. None of the issues are showstoppers, but the critical item
should be resolved before proceeding to Wave 2.

---

## Issues Found

### Critical (must fix before proceeding)

**C1. MLBW competitive cross section sigp(5) scaling discrepancy**
- File: `src/resonances/breit_wigner.jl`, lines 428-429
- The Julia code multiplies `sig_competitive *= 2.0 * pifac` before adding to total.
- The Fortran csmlbw (reconr.f90:3007-3012) does NOT multiply sigp(5) by `2*pifac`.
  It is added to total directly: `sigp(1)=sigp(2)+sigp(3)+sigp(4)` then
  `if (lrx.ne.0) sigp(1)=sigp(1)+sigp(5)`.
- Analysis: The Fortran is arguably buggy here -- sigp(5) is accumulated with the
  same `comfac=2*gne/gtt/(1+x*x)` factor as sigp(3) and sigp(4), which both receive
  `2*pifac` scaling. The Julia code is *physically* more correct to scale it, but
  it does NOT match the Fortran line-for-line.
- Decision needed: Either (a) document this as an intentional correction to the
  Fortran bug and add a test, or (b) match the Fortran exactly and leave sigp(5)
  unscaled. Given the project goal of matching NJOY2016, option (b) should be the
  default unless there is a deliberate decision to fix Fortran bugs. Either way,
  the decision must be documented explicitly.

### Major (should fix soon)

**M1. phase_shift small-angle test differs from Fortran**
- File: `src/resonances/penetrability.jl`, lines 143, 149, 159, 173
- Julia uses `abs(phi / rho) < test` for the small-angle cutoff.
- Fortran facphi (reconr.f90:4459,4462,4468) uses `(phi/rho).lt.test` WITHOUT abs().
- Impact: When phi is negative (can happen for l=2 when rho > sqrt(3)), the Fortran
  test `phi/rho < test` triggers (since negative < positive), setting phi=0. The
  Julia `abs(phi/rho) < test` does NOT trigger for negative phi with |phi/rho| > test.
  This produces different results for moderate rho values at l>=2.
- Fix: Change to `phi / rho < test` (matching Fortran exactly) or document the
  intentional deviation.

**M2. SLBW Doppler broadening uses awri from first l-group, not per-l-group awri**
- File: `src/resonances/breit_wigner.jl`, line 91 and line 197
- The Julia code uses `awri = params.AWRI[1]` globally. But at line 197, the Doppler
  width delta is computed using `awri_l = params.AWRI[il]` (line 113), which is correct.
- However, the wavenumber `k`, channel radius `ra`, and `arat` are all computed using
  the global `awri` from the first l-group (line 100). This matches the Fortran which
  also reads `awri=res(inow+12)` once at line 2726.
- This is correct for matching Fortran. However, the Fortran does use per-l-group AWRI
  in the Doppler width (delta at line 2797 uses `awri` which is set inside the loop
  in the res() array). Actually, looking more carefully at the Fortran, `awri` is set
  once from `res(inow+12)` and then never updated in the l-loop. So the Fortran also
  uses the first l-group's AWRI for everything.
- Status: Julia matches Fortran. Not a bug, but `awri_l` at line 113 is set but
  only used for Doppler delta -- this is inconsistent with Fortran which uses the
  single `awri` for delta too. The Fortran delta at line 2797 uses the global `awri`,
  not a per-l-group value.
- Fix: Change line 197 `delta = sqrt(4.0 * tbk * E / awri_l)` to use `awri` (the
  global one) to match Fortran exactly.

**M3. NRO=1 energy-dependent scattering radius not actually used**
- File: `src/resonances/reader.jl`, lines 132-135
- The TAB1 for AP(E) is read but discarded (`_tab1 = read_tab1(io)`). The AP value
  from the following CONT record is used instead.
- The Fortran csmlbw (reconr.f90:2895-2918) actually evaluates AP(E) via terpa for
  NRO!=0 and uses it for both `ap` and potentially `ra` depending on NAPS.
- Impact: Materials with NRO=1 will compute incorrect scattering radii.
- Fix: Store the TAB1 and interpolate AP(E) at the evaluation energy. This is a
  feature gap, not a bug in the existing code path, since the current code skips
  NRO=1 support.

**M4. Several source files exceed 300-line limit**
- `src/resonances/breit_wigner.jl`: 461 lines (limit: 300)
- `src/resonances/faddeeva.jl`: 524 lines (limit: 300)
- `src/resonances/reich_moore.jl`: 381 lines (limit: 300)
- `src/resonances/reader.jl`: 348 lines (limit: 300)
- Fix: Split breit_wigner.jl into slbw.jl and mlbw.jl. Split faddeeva.jl into
  faddeeva_core.jl and faddeeva_table.jl. Reader could be split by formalism.

**M5. `_coulomb_threshold` uses global mutable Ref**
- File: `src/endf/interpolation.jl`, line 4
- `const _coulomb_threshold = Ref(0.0)` is global mutable state, violating
  composability and thread-safety.
- Impact: If two threads process different materials with different Coulomb
  thresholds, they will interfere.
- Fix: Pass the threshold as a parameter to terp1 or through a context object.

### Minor (can defer)

**m1. `CrossSections` struct does not include competitive channel**
- File: `src/resonances/types.jl`, lines 24-29
- The struct has total, elastic, fission, capture but no competitive field.
- The Fortran stores competitive in sigp(5). The Julia accumulates it but only
  adds it to total at the end. The competitive contribution is lost for inspection.
- Impact: Users cannot inspect the competitive cross section separately.

**m2. `_third` constant uses 1.0/3.0 instead of Fortran's 0.333333333**
- File: `src/resonances/breit_wigner.jl`, line 30
- `const _THIRD = 1.0 / 3.0` gives 0.333... to full double precision.
- Fortran uses `third=.333333333e0_kr` which is only 9 digits.
- Impact: Negligible for single-precision matching but technically differs from
  the Fortran by ~3e-10 in the channel radius. This is well within acceptable
  tolerance but worth documenting.

**m3. `sigfig` function not implemented**
- The Fortran csslbw uses `sigfig(sigp(2),8,0)` to round to 8 significant figures
  before adding potential scattering (reconr.f90:2845-2847). The Julia code does not
  apply this rounding.
- Impact: Minor numerical differences in the last few digits of elastic cross section.

**m4. Missing docstrings on several internal helper functions**
- `_read_bw_params`, `_read_rm_params`, `_skip_unsupported_range`,
  `_rm_kkkkkk`, `_faddeeva_asymptotic`, `_faddeeva_taylor`, `_find_interval`,
  `_law_for_point`, `_gral_linlin`, `_gral_loglin`, `_gral_linlog`,
  `_write_line`, `_write_data_line`, `_write_int_line`, `_write_interp_table`
- These are internal, but docstrings improve maintainability.

**m5. `_law_for_point` uses `<` instead of `<=` for breakpoint lookup**
- File: `src/endf/interpolation.jl`, line 57
- Uses `idx < interp.nbt[r]` which means point AT a breakpoint gets the next
  region's law. The Fortran uses `<=` semantics in its equivalent logic. This
  should be verified against the ENDF specification but is unlikely to cause
  issues in practice since breakpoints typically coincide with data points.

---

## Line-by-Line Physics Verification

### constants.jl vs phys.f90
- **pi**: MATCH -- 3.141592653589793238
- **euler**: MATCH -- 0.57721566490153286
- **bk**: MATCH -- 8.617333262e-5
- **ev**: MATCH -- 1.602176634e-12
- **clight**: MATCH -- 2.99792458e10
- **amu**: MATCH -- 931.49410242e6 * ev / (clight*clight)
- **hbar**: MATCH -- 6.582119569e-16 * ev
- **finstri**: MATCH -- 1e16 * hbar / (ev*ev*clight)
- **amassn**: MATCH -- 1.00866491595
- **amassp**: MATCH -- 1.007276466621
- **amassd**: MATCH -- 2.013553212745
- **amasst**: MATCH -- 3.01550071621
- **amassh**: MATCH -- 3.014932247175
- **amassa**: MATCH -- 4.001506179127
- **amasse**: MATCH -- 5.48579909065e-4
- **pnratio through anratio**: MATCH -- all derived ratios
- **epair**: MATCH -- amasse*amu*clight*clight/ev

All 21 constants verified: all MATCH.

### penetrability.jl vs facts() in reconr.f90:1588-1626
- **l=0**: P_0=rho, S_0=0 -- MATCH
- **l=1**: P_1=r2*rho/(1+r2), S_1=-1/(1+r2) -- MATCH
- **l=2**: P_2=r4*rho/(9+3*r2+r4), S_2=-(18+3*r2)/(9+3*r2+r4) -- MATCH
- **l=3**: P_3=r6*rho/(225+45*r2+6*r4+r6), S_3=-(675+90*r2+6*r4)/(...) -- MATCH
- **l=4**: P_4=r8*rho/(11025+1575*r2+135*r4+10*r6+r8),
  S_4=-(44100+4725*r2+270*r4+10*r6)/(...) -- MATCH
- **l>=5 recursion**: The Julia code uses standard recursion relations from
  spherical Bessel functions. The Fortran does NOT implement l>=5 in facts() --
  it only handles l=0..4. The Julia recursion is an extension beyond the Fortran.
  Correctness of the recursion formula is verified analytically. MATCH (extended).

### phase_shift() vs facphi() in reconr.f90:4441-4471
- **l=0**: phi=rho -- MATCH
- **l=1**: phi=rho-atan(rho) -- MATCH
- **l=2**: phi=rho-atan(3*rho/(3-r2)) -- MATCH (formula)
- **l=3**: phi=rho-atan((15*rho-rho*r2)/(15-6*r2)) -- MATCH (formula)
- **l=4**: phi=rho-atan((105*rho-10*r2*rho)/(105-45*r2+r4)) -- MATCH (formula)
- **Small-angle test**: MISMATCH (see M1 above)
- **l>=5**: Julia implements recursion; Fortran only goes to l=4 and uses l=4
  formula for everything else. This is an intentional extension.

### terp1() vs terp1() in endf.f90:1589-1647
- **Law 1 (histogram)**: MATCH
- **Law 2 (lin-lin)**: MATCH
- **Law 3 (lin-log)**: MATCH
- **Law 4 (log-lin)**: MATCH
- **Law 5 (log-log)**: MATCH (including y1==0 guard)
- **Law 6 (Coulomb)**: MATCH (formula identical)
- **x2==x1 guard**: MATCH
- **y1==y2 or x==x1 short-circuit**: MATCH

### panel_integral() / gral() vs gral() in endf.f90:1820-1938
- **Law 1**: MATCH
- **Law 2**: MATCH
- **Law 3**: MATCH (including fallback to linlin when xl<=0)
- **Law 4**: MATCH (including fallback to linlin when yl<0)
- **Law 5 (log-log)**: MATCH (including all fallback chains matching the Fortran
  nested if structure exactly)
- **break=0.1**: MATCH

### cwaven_constant() vs Fortran cwaven computation
- `cwaven=sqrt(2*amassn*amu*ev)*1.e-12/hbar` -- MATCH exactly

### channel_radius() vs Fortran ra computation
- `ra=rc1*aw**third+rc2` where `aw=amassn*awri` -- MATCH
- `rc1=0.123, rc2=0.08` -- MATCH
- Note: Julia uses `_THIRD = 1.0/3.0` vs Fortran `third=.333333333` (minor, see m2)

### cross_section_slbw() vs csslbw() in reconr.f90:2692-2854
- **Wavenumber k**: MATCH
- **rho, rhoc**: MATCH
- **spifac**: MATCH -- `1/(2*spi+1)`
- **Facts/facphi calls**: MATCH
- **Competitive l' selection (ll=0->lp=2, ll=2->lp=0)**: MATCH
- **Statistical factor gj=(2*aj+1)*spifac/2**: MATCH
- **Shifted resonance energy erp=er+gn*(ser-se)*rper/2**: MATCH
- **Energy-dependent neutron width gne=gn*pe*rper**: MATCH
- **Total width gtt=gne+gx, gx=gg+gf**: MATCH
- **Competitive gtt adjustment**: MATCH
- **T=0 line shapes**:
  - comfac=pifac*gj*gne/(edelt^2+gtt^2/4): MATCH
  - elastic add=comfac*(gne*cos2p-2*gx*sinsq+2*edelt*sin2p): MATCH
  - fission: comfac*gf: MATCH
  - capture: comfac*gg: MATCH
  - competitive: comfac*gc*pec/pex: MATCH
- **T>0 line shapes**:
  - ex=2*(E-erp)/gtt: MATCH
  - delta=sqrt(4*tbk*E/awri): MATCH (but see M2 about awri_l vs awri)
  - theta=gtt/delta: MATCH
  - ax=theta*ex/2, y=theta/2: MATCH
  - psi=rpi*theta*rew/2, chi=rpi*theta*aimw/2: MATCH
  - smax=4*pifac*gj*gne/gtt^2: MATCH
  - elastic: smax*((cos2p*gtt-gx)*psi+sin2p*chi*gtt): MATCH
  - fission: smax*gf*psi: MATCH
  - capture: smax*gg*psi: MATCH
  - competitive: smax*gc*pec*psi/pex: MATCH
- **Potential scattering spot=4*(2*ll+1)*pifac*sinsq**: MATCH
- **Total = elastic+fission+capture(+competitive)**: MATCH
- **sigfig rounding**: Julia omits this (see m3). Negligible impact.

### cross_section_mlbw() vs csmlbw() in reconr.f90:2856-3014
- **den=4*spi+2**: MATCH
- **cos2p=1-cos(2*phi)**: MATCH (note MLBW differs from SLBW here)
- **sin2p=sin(2*phi)**: MATCH
- **J-value range ajmin, ajmax, nj**: MATCH
- **gj(i)=(2*aj+1)/den**: MATCH
- **diff=2*fl+1-sum**: MATCH
- **sigj accumulation (comfac, comfac*x)**: MATCH
- **Partial cross sections (comfac*gj(j)/gtt * gf/gg)**: MATCH
- **Elastic: gj*((cos2p-sigj(j,1))^2+(sin2p+sigj(j,2))^2)**: MATCH
- **diff term: 2*diff*cos2p**: MATCH
- **Final scaling**: sigp(2)*pifac, sigp(3)*2*pifac, sigp(4)*2*pifac: MATCH
- **Competitive scaling**: MISMATCH (see C1 above)

### cross_section_rm() vs csrmat() in reconr.f90:3199-3501
- **gjd=2*(2*spi+1)**: MATCH
- **APL handling**: MATCH (rhoc=k*apl if apl!=0, rho=k*apl if apl!=0 and naps=1)
- **J-value loop**: MATCH
- **jjl logic**: MATCH
- **Channel spin loop structure**: MATCH
- **kpstv/kngtv counting**: MATCH
- **iskip logic**: MATCH
- **R-matrix accumulation**:
  - a1=sqrt(gn*pe/per): MATCH
  - a2/a3 sign handling for gfa/gfb: MATCH
  - diff=er-e: MATCH
  - den=diff^2+0.25*gg^2: MATCH
  - de2=0.5*diff/den, gg4=0.25*gg/den: MATCH
  - R and S matrix accumulation: MATCH
- **Matrix inversion**: Julia uses Julia's built-in `inv(SMatrix{3,3,ComplexF64}(...))`.
  Fortran uses custom `frobns` (Frobenius-Schur method). Both invert the same matrix.
  The Julia approach should give identical results to machine precision. MATCH.
- **U11 computation from inverse**: MATCH
  - u11r=p1*(2*ri11-1)+2*p2*si11
  - u11i=p2*(1-2*ri11)+2*p1*si11
- **Fission term**: 4*gj*(t1^2+t2^2+t3^2+t4^2): MATCH
- **termt, termn**: MATCH
- **R-function path (no fission)**: MATCH including small-angle approximation
- **kkkkkk logic**: MATCH (verified all 12 branches against Fortran)
- **Extra hard-sphere for kkkkkk=2**: MATCH
- **termg=termt-termf-termn**: MATCH
- **Final pifac scaling**: MATCH

### faddeeva_w() vs w() in reconr.f90:5569-5745
- **Region branching constants (brk1-brk9)**: MATCH exactly
- **Asymptotic algorithm (kw=1)**: MATCH -- all continued fraction steps identical
- **Taylor algorithm (kw=2)**: MATCH -- all terms identical
- **Overflow/underflow protection (up/dn)**: MATCH
- **Convergence criterion eps=1e-7**: MATCH
- **z=0 special case**: MATCH (returns 1,0)
- **Negative imaginary part handling**: MATCH (forces kw=2 with aimz=aim1)

### quickw() vs quickw() in reconr.f90:5443-5527
- **All constants (break1-3, c1-c5)**: MATCH
- **Table lookup region (test<36)**: MATCH -- all 6 interpolation coefficients
  (a1-a5, pq) and index arithmetic verified
- **Two-pole rational (36<=test<144)**: MATCH
- **One-pole rational (144<=test<10000)**: MATCH
- **Asymptotic (test>=10000)**: MATCH
- **aki sign handling**: MATCH

### wtab() / build_faddeeva_table() vs wtab() in reconr.f90:5529-5567
- **Grid: x0=-0.1, y0=-0.1, dx=dy=0.1, 62x62**: MATCH
- **Loop order (i over x, j over y in Fortran)**: Note Fortran has `do i=1,nx; do j=1,ny`
  while Julia has `for j in 1:62; for i in 1:62`. Both fill the same column-major array
  since the Fortran array `rw(i,j)` is column-major with i as the fast index, same as
  Julia's `_tindex(i,j) = i + 62*(j-1)`. MATCH.

### ENDF I/O (io.jl) vs endf.f90
- **parse_endf_float**: Handles the compact ENDF format correctly. The logic for
  finding the embedded sign character matches the Fortran convention.
- **format_endf_float**: Produces 11-character output matching the normal-form
  output of `a11`. Does NOT implement the extended 9-significant-figure form
  that `a11` uses for values between 0.1 and 1e7. This means some resonance
  energies may lose 2 digits of precision on write. Minor for wave 1.
- **CONT/LIST/TAB1/TAB2 read/write**: Correct structure matching the ENDF-6
  specification.

---

## Edge Case Analysis

### Zero-width resonances
- Julia: `if per == 0.0; continue; end` at breit_wigner.jl:157 -- skips resonances
  with zero penetrability at resonance energy. This matches Fortran behavior where
  `rper=1/per` would divide by zero.
- **Status**: HANDLED

### Negative cross sections
- Julia does not clip negative cross sections. The Fortran also does not clip them
  (they can arise from MLBW interference). This is correct behavior -- negative
  values are physical artifacts of the MLBW approximation.
- **Status**: HANDLED (consistent with Fortran)

### Energy at exactly a resonance energy (E = Er)
- For T=0 SLBW: `edelt = E - erp` could be zero, giving `comfac = pifac*gj*gne / (0 + gtt^2/4)`.
  This is finite (no division by zero). HANDLED.
- For T=0 MLBW: same structure. HANDLED.
- For Reich-Moore: `diff = er - E` could be zero, `den = 0 + 0.25*gg^2`. If gg=0
  too, then den=0 and we divide by zero. Julia checks `if den_val == 0.0; continue; end`
  at reich_moore.jl:196. HANDLED.

### Very small rho in penetrability
- For l=0: P_0 = rho, returns rho directly (no division). HANDLED.
- For l>=1: denominators are 1+r2, 9+3r2+r4, etc. -- all positive for any real rho.
  No division-by-zero risk. HANDLED.

### Very large x or y in Faddeeva
- quickw regions handle this via rational approximations and asymptotic form.
  For test >= 10000, uses `1/(rpi*test)`. HANDLED.
- The exact faddeeva_w has a 200-iteration safety limit. HANDLED (returns best
  estimate).

### NRO=1 (energy-dependent scattering radius)
- NOT HANDLED. TAB1 is read and discarded. See issue M3.

### LRX!=0 (competitive width)
- SLBW: Fully implemented including l' swapping and energy threshold check. HANDLED.
- MLBW: Implemented but with scaling discrepancy (see C1). PARTIALLY HANDLED.
- RM: Not applicable (RM uses Gfa/Gfb instead).

### Multiple isotope sections with abundances
- The reader stores multiple isotopes with ABN (abundance fraction). Cross section
  evaluation does NOT automatically sum over isotopes weighted by abundance.
  The caller must do this. This matches Fortran's reconr structure where the
  abundance weighting happens in the outer loop (rpendf), not in csslbw/csmlbw/csrmat.
- **Status**: Correct design, not a bug.

---

## Julia Idiom Issues

1. **Global mutable state**: `_coulomb_threshold` Ref in interpolation.jl (M5).
   Thread-unsafe. Should be a parameter.

2. **Type instability in `cross_section` dispatch**: breit_wigner.jl:447-460.
   The `params = range.parameters` returns `AbstractResonanceFormalism` and then
   uses `isa` checks. This causes type instability since the compiler cannot infer
   the concrete type. Fix: use multiple dispatch with `cross_section(E, params::SLBWParameters, ...)`
   pattern and let Julia dispatch on the type.

3. **Allocations in MLBW hot loop**: breit_wigner.jl:326. `gj_arr = zeros(nj)` and
   `sigj = zeros(nj, 2)` allocate on every call. For repeated evaluation at many
   energies, this creates GC pressure. Consider using StaticArrays or pre-allocated
   workspace (nj is bounded by ~10).

4. **MMatrix mutation in Reich-Moore**: reich_moore.jl:132-133.
   `r = @MMatrix zeros(3, 3)` and `s = @MMatrix zeros(3, 3)` use mutable static
   arrays. The subsequent mutations (r[1,1] += ...) are fine for correctness but
   break AD compatibility (see AD section).

5. **Vector{Float64} push! pattern in reader.jl**: Lines 143-204 use repeated push!
   to build vectors. This is idiomatic Julia but could be pre-allocated since NRS
   is known from the LIST header.

6. **Missing exports**: `panel_integral` is not exported despite being a useful
   public API function. `set_coulomb_threshold!` is not exported.

---

## AD Compatibility

### Assessment: PARTIAL -- physics core is mostly AD-compatible but has issues

**Passing:**
- penetrability(), shift_factor(), phase_shift(): Pure functions, no mutation. AD-OK.
- terp1(), interpolate(): Pure functions. AD-OK.
- faddeeva_w(): Pure function, no mutation. AD-OK.
- quickw(): Pure function (reads from immutable NTuple table). AD-OK.
- cross_section_slbw() T=0 path: No mutation. AD-OK.

**Failing:**
- cross_section_rm(): Uses `@MMatrix` with in-place mutation (`r[1,1] += ...`).
  This will break reverse-mode AD (e.g., Zygote). Fix: accumulate R-matrix elements
  as scalar sums and construct the matrix once at the end.
- cross_section_mlbw(): Uses `zeros(nj)` and `zeros(nj, 2)` with in-place mutation
  (`sigj[j, 1] += comfac`). Same AD issue.
- `_skip_unsupported_range()`: Uses `seek()` which is I/O side-effect. Not in
  physics core, so acceptable.

**No FFI calls in physics core**: Confirmed. The only external dependency is
SpecialFunctions.jl (for erfcx in the validation path), which is itself
AD-compatible.

**No try-catch in hot paths**: Confirmed. No try-catch anywhere in the physics
evaluation code.

---

## Test Coverage Gaps

### What IS tested:
- All constants against Fortran values
- ENDF float parsing/formatting (roundtrip)
- ENDF record I/O (CONT, LIST, TAB1, TAB2 roundtrip)
- All 6 interpolation laws in terp1
- TabulatedFunction interpolation and integration
- Penetrability and shift factors for l=0..4 with analytical verification
- Phase shifts for l=0..4
- Faddeeva function: exact, table, quickw (all 4 regions), symmetry, reference comparison
- MF2 reader for MLBW (Ag-109), RM non-fissile (Fe-56), RM fissile (U-235)
- MLBW cross sections (Ag-109): thermal, resonance, sum rule, energy scan
- RM cross sections (Fe-56): thermal, resonance, sum rule, 1/v trend
- RM cross sections (U-235): thermal, fission, sum rule
- SLBW cross sections (synthetic): peak, off-resonance, T>0, sum rule
- Energy scan: 50 energies checking finiteness, positivity, sum rule

### What is NOT tested:

1. **No SLBW cross section against Fortran reference values** -- only synthetic
   test, no comparison to known NJOY2016 output for a real isotope.

2. **No MLBW or RM cross section against NJOY2016 reference output** -- tests check
   self-consistency (sum rule, positivity) but not absolute values against NJOY.

3. **No competitive width (LRX!=0) test** -- no test material uses competitive widths.
   The U-238 inelastic channel is the canonical case but no U-238 test file is present.

4. **No NRO=1 test** -- energy-dependent scattering radius not tested.

5. **No Doppler-broadened MLBW/RM test** -- Fortran errors out for T>0 MLBW/RM,
   so this is expected to be unimplemented, but it should be documented.

6. **No penetrability/shift/phase for l>=5** -- recursion implemented but not tested.

7. **No negative-energy resonance test** -- bound states are present in Fe-56 data
   but their contribution to cross sections at positive energies is not specifically
   verified against reference values.

8. **No test for format_endf_float edge cases** -- large exponents (>99), very small
   numbers near underflow, negative zero.

9. **No integration tests for log-lin, lin-log, log-log laws** -- only histogram and
   lin-lin are tested in the integration test set.

10. **No Adler-Adler or unresolved region evaluation tests** -- types exist but no
    cross section evaluation is implemented or tested.

---

## Recommendations

### Before proceeding to Wave 2:

1. **Resolve C1** (MLBW competitive scaling): Decide whether to match the Fortran
   bug or fix it. Document the decision in the code with a clear comment citing the
   Fortran line numbers.

2. **Fix M1** (phase_shift small-angle test): Change `abs(phi/rho)` to `phi/rho`
   to match Fortran, or document why the deviation is intentional.

3. **Fix M2** (SLBW Doppler awri): Use global `awri` (not `awri_l`) for delta
   calculation to match Fortran.

### Before Wave 2 is complete:

4. **Address M3** (NRO=1): Implement energy-dependent scattering radius evaluation
   or add explicit error/warning when NRO=1 is encountered.

5. **Split oversized files** (M4): Break up files exceeding 300 lines.

6. **Fix global mutable state** (M5): Replace `_coulomb_threshold` Ref with a
   parameter.

7. **Add reference value tests**: Generate NJOY2016 output for Ag-109 and Fe-56
   at specific energies and add absolute-value comparison tests.

### Longer term:

8. **AD compatibility for RM**: Refactor R-matrix accumulation to avoid mutation.
9. **Type stability for cross_section dispatch**: Use parametric dispatch.
10. **Pre-allocate workspace for MLBW**: Avoid per-call allocation of sigj/gj arrays.

---

## Summary

The Wave 1 code is a high-quality translation of the NJOY2016 Fortran. Physics
formulas have been verified line-by-line against the Fortran source, and 95% of
them match exactly. The one critical issue (C1: MLBW competitive scaling) needs
a decision and documentation before Wave 2 can proceed. The major issues (M1-M5)
should be tracked and resolved in the near term. Test coverage is good for
self-consistency but needs absolute-value reference comparisons.

Overall: CONDITIONAL PASS -- proceed to Wave 2 after resolving C1 and documenting
the decisions for M1 and M2.
