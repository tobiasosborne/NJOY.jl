# Wave 2 Review Report

## Verdict: CONDITIONAL PASS

Wave 2 delivers a functional RECONR pipeline that correctly handles the H-2 (LRU=0)
test case and provides a well-structured adaptive grid algorithm. However, several
deviations from the Fortran reference, edge-case gaps, file-length violations, and
a critical silent-failure bug must be addressed before this can be considered
production-ready.

---

## Issues Found

### Critical

**C1. Silent exception swallowing in `sigma_mf2` and `read_mf3_sections`**

- File: `src/processing/reconr.jl:629-636` (`sigma_mf2`) and `reconr.jl:483-490`
  (`read_mf3_sections`)
- Both functions use bare `try ... catch` blocks with no logging, no rethrow, and
  no error message. If `cross_section()` throws for a valid resonance (e.g., due to
  a parameter-reading bug or an untested formalism path), the exception is silently
  swallowed and that resonance's contribution is simply omitted from the result.
  This can produce silently wrong cross sections with no diagnostic.
- In `read_mf3_sections`, a parsing failure silently drops an entire MF3 section.
  The user will never know data was lost.
- **Fix**: At minimum, log a warning. Better: catch only expected exceptions
  (e.g., `UnsupportedFormalismError`) and let genuine bugs propagate.

**C2. Step-ratio guard missing `xm < eresr` condition**

- File: `src/processing/adaptive_grid.jl:282-286`
- Fortran (reconr.f90:2419): `if (in.gt.3.and.dx.gt.est.and.xm.lt.eresr) go to 175`
- The Fortran applies the step-ratio guard ONLY in the resolved resonance region
  (`xm < eresr`). The Julia code omits this third condition entirely. This means
  the step-ratio guard fires unnecessarily in the unresolved region or above
  resonances, causing excessive grid refinement and potentially the stack-overflow
  error for wide energy ranges.
- The `adaptive_reconstruct` function is designed to be generic (knows nothing about
  resonance ranges), so this guard condition must be injected externally -- e.g.,
  via an optional upper-energy parameter in `AdaptiveConfig`, or by restricting the
  grid passed to `adaptive_reconstruct` to the resolved range only.
- Currently `reconr.jl:437-445` does filter to `[eresl, eresh]` but this includes
  the unresolved region. The Fortran specifically uses `eresr` (resolved upper bound),
  not `eresh` (overall upper bound).

### Major

**M1. Deviation in convergence comparison: `abs(mid_val[j])` vs `sig(j+1)`**

- File: `src/processing/adaptive_grid.jl:272,275`
- Julia: `dm > err_eff * abs(mid_val[j])` -- compares against the TRUE midpoint value.
- Fortran (reconr.f90:2395): `dm(j).gt.errn*sig(j+1)` -- compares against `sig(j+1)`,
  which is the TRUE midpoint cross section for that component.
- These are actually equivalent because `sig` is populated by `call sigma(xm,sig,sigd)`
  at line 2373, so `sig(j+1)` IS the true value at the midpoint. The Julia
  `abs(mid_val[j])` is the same. **No issue here on closer inspection.** However,
  the Julia uses `abs()` while the Fortran does not -- this matters when cross
  sections could be negative (which they shouldn't be after clamping, but before
  clamping in some intermediate states they could be). Minor behavioral difference.

**M2. `merge_background!` omits non-primary MT channels**

- File: `src/processing/reconr.jl:300-324`
- The function only handles MT=2 (elastic), MT=18/19 (fission), and MT=102
  (capture). NJOY's `emerge` (reconr.f90:4757-4770) also handles additional
  resonance MTs via `mmtres` (e.g., competitive reactions, MT=51 for inelastic
  when LRX != 0, and extended MT lists for LRF=7 evaluations).
- For materials with only standard resonance channels this is fine, but for
  materials with competitive reactions (LRX flag) or R-Matrix Limited (LRF=7)
  evaluations, backgrounds will be missing.

**M3. `merge_background!` missing overlap logic for resolved-unresolved boundary**

- File: `src/processing/reconr.jl:313`
- Julia: `if e >= eresr && e < eresh` -- suppresses backgrounds in the
  resolved-unresolved overlap region for primary channels.
- Fortran (reconr.f90:4800): `if (eg.ge.eresr.and.eg.lt.eresh.and.itype.gt.0.and.itype.lt.5) sn=0`
- The Fortran additionally requires `itype > 0` (i.e., this is a resonance reaction)
  AND `itype < 5` (only primary channels). The Julia check is slightly broader
  because it does not verify `itype > 0` -- it will suppress backgrounds for
  non-resonance reactions that happen to be MT=2/18/19/102 in the overlap region.
  For standard materials this is equivalent, but could differ for exotic evaluations.

**M4. `build_grid` does not include the `lunion`-style 1/v grid tightening**

- File: `src/processing/reconr.jl:142-170`
- Fortran's `lunion` (reconr.f90:2112-2174) performs extensive linearization of
  MF3 data during grid construction: it adds midpoints to ensure 1/v shapes are
  represented within tolerance, forces 1-2-5-10 decade points, and linearizes
  nonlinear interpolation regions. The Julia `build_grid` simply collects MF3
  breakpoints without any linearization.
- For H-2 (mostly linear MF3 data) this works, but for materials with log-log
  or lin-log interpolated MF3 data, the initial grid will be too coarse and the
  adaptive algorithm may struggle or produce incorrect results since `merge_background!`
  only does linear interpolation at the final grid points.

**M5. `reconr.jl` is 679 lines and `pendf_writer.jl` is 380 lines**

- File: `src/processing/reconr.jl` (679 lines), `src/processing/pendf_writer.jl`
  (380 lines)
- The project's 300-line-per-file limit is exceeded by both files. `reconr.jl` is
  more than double the limit.
- **Recommendation**: Split `reconr.jl` into at least three files:
  - `reconr_evaluator.jl` (build_evaluator, sigma_mf2, ~100 lines)
  - `reconr_grid.jl` (build_grid, _add_mf2_nodes!, peak nodes, ~130 lines)
  - `reconr_merge.jl` (merge_background!, merge_background_legacy, ~100 lines)
  - `reconr.jl` (reconstruct, reconr, MF3 reader, ~200 lines)
- Similarly, `pendf_writer.jl` could separate the legacy interface into its own file.

### Minor

**m1. `_write_fend` produces duplicate SEND for MF3**

- File: `src/processing/pendf_writer.jl:125-128` and `pendf_writer.jl:188`
- `_write_mf3_from_matrix` writes a SEND record (mat, 3, 0) at line 188 for each
  MT section. Then `_write_fend` writes another (mat, 3, 0) at line 126, followed
  by FEND (mat, 0, 0) at line 127.
- Comparing with the reference tape (referenceTape100): each MF3 section ends with
  its own SEND record (e.g., line 282: `128 3  099999` for MT=1, line 543:
  `128 3  0  261` for MT=2), and then a single FEND (line 851: `128 0  0  569`)
  appears after all MF3 sections.
- The Julia code writes an EXTRA SEND record between the last MT's SEND and the
  FEND. This is a format deviation. The `_write_fend` function should only write
  the FEND record `(mat, 0, 0)`, not the additional `(mat, 3, 0)`.

**m2. TPID record MAT field**

- File: `src/processing/pendf_writer.jl:116-117`
- The TPID record writes `mat` in the MAT field position. The reference tape
  (line 1) shows `1 0  0    0` -- i.e., MAT=1, MF=0, MT=0 for the TPID line.
  The Julia code writes the actual MAT number (e.g., 128) which differs from
  NJOY's convention of using MAT=1 for the TPID.

**m3. Sequence number handling in PENDF output**

- File: `src/processing/pendf_writer.jl:122`
- `_write_record_sep` always writes `99999` as the sequence number. The reference
  tape shows that SEND records use `99999` only for the first occurrence (line 14,
  282), while subsequent SEND/FEND records use sequential numbering (line 543:
  `128 3  0  261`, line 851: `128 0  0  569`). This means the Julia PENDF output
  uses non-standard sequence numbering for delimiters.

**m4. `_write_mf3_from_matrix` writes ZA=0, AWR=0 in HEAD**

- File: `src/processing/pendf_writer.jl:169`
- The HEAD record for each MF3 section writes ZA=0.0 and AWR=0.0. The reference
  tape (line 283, 544) shows `1.002000+3 1.996800+0` (i.e., ZA and AWR from the
  material). While this may be accepted by downstream codes, it is a deviation
  from the reference.

**m5. `round_sigfig` uses `floor(Int, aa)` while Fortran uses `int(aa)` with conditional**

- File: `src/processing/adaptive_grid.jl:75`
- Julia: `ipwr = floor(Int, aa)` -- always rounds toward negative infinity.
- Fortran (util.f90:379-380): `ipwr=int(aa)` then `if (aa.lt.zero) ipwr=ipwr-1`.
  `int()` in Fortran truncates toward zero, then the conditional adjusts for
  negative values. For negative `aa` (values < 1.0), `floor` and `int-1` give
  the same result. For positive `aa`, `floor` and `int` give the same result.
  **These are actually equivalent.** No action needed.

**m6. `_pendf_parse_int` duplicates `_parse_int`**

- File: `src/processing/pendf_writer.jl:376-380` vs `src/endf/io.jl:62`
- These are nearly identical functions. The only difference is `_parse_int` returns
  `Int32` and `_pendf_parse_int` returns `Int`. This duplication should be
  eliminated by using `_parse_int` and converting as needed.

**m7. Allocation in `_write_data_values` inner loop**

- File: `src/processing/pendf_writer.jl:139-154`
- The `" "` string literal (11 spaces) is allocated on each call for padding.
  This is minor but could be a `const` for clarity and to avoid repeated allocation.

---

## Fortran Comparison

### adaptive_grid.jl vs resxs (reconr.f90:2240-2569)

| Feature | Fortran | Julia | Match? |
|---------|---------|-------|--------|
| 3-tier convergence | primary/relaxed+integral/step-guard | Same structure | YES |
| Thermal tightening threshold | `trange=0.4999` | `thermal_threshold=0.4999` | YES |
| Thermal factor | `/5` for both err and errmax | `/5.0` for both | YES |
| Thermal check location | `x(i-1).lt.trange` (upper stack point) | `sx[i-1] < config.thermal_threshold` | YES |
| sigfig ndig logic | `ndig=9; if (xm>0.1 and xm<1) ndig=8` | Same | YES |
| sigfig convergence guard | `xm>sigfig(x(i),ndig,+1) and xm<sigfig(x(i-1),ndig,-1)` | Same | YES |
| sigfig rounding | 7-digit or ndig fallback | Same | YES |
| Step-ratio guard | `estp=4.1`, `in>3 and dx>est and xm<eresr` | `_STEP_RATIO=4.1`, missing `xm<eresr` | **NO** (C2) |
| Stack management | `i` grows upward, `x(1)=eg, x(2)=egl` | Same convention | YES |
| Stack overflow | `if (i.gt.ndim) call error(...)` | `error(...)` if `i > max_depth` | YES |
| ndim (max_depth) | 30 | 30 | YES |

### reconr.jl vs sigma (reconr.f90:2571-2667)

| Feature | Fortran | Julia | Match? |
|---------|---------|-------|--------|
| Energy range check | `e.ge.elt(i).and.e.lt.eht(i)` | `E >= rng.EL && E < rng.EH` | YES |
| LRU=0 skip | via `mode=0 -> csnorp` | `Int(rng.LRU) > 0` skips LRU=0 | Equivalent |
| Abundance weighting | `sig(j)=sig(j)+sigp(j)*abn` | `sig_t += max(0.0, sigp.total) * abn` | YES |
| Negative clamping | `if (sigp(j).lt.zero) sigp(j)=0` | `max(0.0, ...)` | YES |
| Unresolved contribution | `call sigunr(e,sigp)` added for `e>=eresu and e<eresh` | Not implemented | **MISSING** |
| Formalism modes | 0-7 (SLBW, MLBW, RM, AA, hybrid, RML) | SLBW, MLBW, RM only | Partial |
| mode 11/12 skip | `if (mode.ne.11.and.mode.ne.12)` | Not needed (no URR modes) | OK for now |

### reconr.jl vs lunion (reconr.f90:1771-2238)

| Feature | Fortran | Julia | Match? |
|---------|---------|-------|--------|
| MF3 breakpoints | Yes, with linearization | Yes, without linearization | **PARTIAL** (M4) |
| 1-2-5-10 decade points | Yes (lines 2116-2128) | No | **MISSING** |
| Thermal 0.0253 eV | Yes (line 2127) | Yes (line 160) | YES |
| 1/v grid tightening | Yes (lines 2131-2141, `stpmax`) | No | **MISSING** |
| Resonance nodes | Er, Er+/-hw with sigfig rounding | Same approach | YES |
| Range boundary shading | sigfig(el,7,-1), sigfig(el,7,+1) | Same | YES |
| Deduplication | Yes, with resonance-boundary filtering | sort/unique/filter | Simpler but OK |
| Nonlinear linearization | Yes (lines 2143-2153) | No | **MISSING** |

### reconr.jl vs emerge (reconr.f90:4646-4982)

| Feature | Fortran | Julia | Match? |
|---------|---------|-------|--------|
| Primary channels (MT=2,18,102) | Yes | Yes | YES |
| Resolved-unresolved overlap | Suppresses backgrounds for itype 1-4 | Similar check | Close (M3) |
| Elastic floor | `sn=small` (1e-8) | `elastic=1e-8` | YES |
| sigfig rounding | `sn=sigfig(sn,7,0)` | `round_sigfig(values[i,j],7,0)` | YES |
| Redundant reaction accumulation | Accumulates totals for MT=1,3,4,18,101 | Recomputes total=e+f+c | Different approach |
| Additional resonance MTs | mmtres array for LRF=7 | Not handled | **MISSING** (M2) |
| Threshold handling | Complex threshold logic with sigfig | Not implemented | **MISSING** |

### pendf_writer.jl vs recout (reconr.f90:4984-5441)

| Feature | Fortran | Julia | Match? |
|---------|---------|-------|--------|
| TPID record | MAT=1 convention | Uses actual MAT | **DIFFERS** (m2) |
| MF1/MT451 HEAD | ZA, AWR, LRP, LFI, etc. | Simplified | Functional |
| Dictionary | Complex with redundant MTs | Basic version | Simpler but OK |
| MF2/MT151 | ZA, AWR, spin, scattering radius | Minimal placeholder or copy | OK |
| MF3 HEAD | ZA, AWR, 0, 99 | ZA=0, AWR=0, 0, 99 | **DIFFERS** (m4) |
| SEND/FEND/MEND/TEND | Standard delimiters | Extra SEND before FEND | **DIFFERS** (m1) |
| Sequence numbers | Sequential | 99999 for all delimiters | **DIFFERS** (m3) |

---

## Edge Case Analysis

### E exactly at range boundary (EL or EH)

- `build_evaluator` uses `E >= rng.EL && E < rng.EH` (line 101). The upper bound
  is exclusive, matching Fortran's `e.ge.elt(i).and.e.lt.eht(i)`. This is correct.
- However, resonance nodes include `round_sigfig(eh, 7, +1)` (line 197), which
  places a grid point slightly ABOVE EH. At this point, the evaluator returns zero
  for that range. This is correct behavior -- the adaptive algorithm will capture
  the transition.

### Material with no resonances (LRU=0 only, e.g., H-2)

- Handled correctly in both `reconstruct` (line 424-434) and `reconr` (line 551-578).
- The `_resonance_bounds` function returns `eresl=Inf, eresh=0.0` when all ranges
  have LRU=0, which triggers the no-resonance path.
- `reader.jl` (line 77-91) correctly handles LRU=0 by reading the SPI/AP CONT
  record and creating a minimal SLBWParameters with empty resonance lists.
- **This is the LRU=0 fix mentioned in the review scope and it works correctly.**

### Material with overlapping resonance ranges

- Multiple ranges with overlapping [EL, EH] intervals: the evaluator sums
  contributions from ALL ranges where E falls within bounds (line 100-107).
  This matches Fortran behavior.
- `_resonance_bounds` correctly computes the overall min/max across all ranges.

### Very narrow resonances (high-Q reactions)

- The `_push_peak_triplet!` function (line 257-271) computes `ndig` based on
  `log10(er / max(hw/10, 1e-30))`, clamped to [5, 9]. This ensures narrow
  resonances get more significant figures in their grid nodes, preventing
  the triplet from collapsing. This matches the Fortran approach.
- However, the adaptive algorithm could still struggle if the initial grid
  misses a very narrow resonance entirely (no grid node near Er). This is
  mitigated by the peak triplet nodes, but only if the resonance parameters
  are read correctly.

### MF3 section with zero cross section (below threshold)

- `merge_background!` (line 308) checks `bg == 0.0 && continue`, so zero
  backgrounds are correctly skipped.
- However, there is no threshold logic in `merge_background!` -- if MF3
  data exists below threshold (e.g., fission below threshold energy), it
  will still be interpolated and added. The Fortran `emerge` has explicit
  threshold handling (lines 4792-4794) that sets `sn=0` below threshold.
  This could produce small nonzero fission cross sections below threshold.

### Empty MF3 sections

- If `read_mf3_sections` encounters an MF3 section that fails to parse, it
  silently catches the exception and moves on (line 483-490). Empty sections
  are not explicitly handled but would be caught by the `tab.x` being empty,
  in which case `interpolate` returns 0.0 for all energies.

---

## Numerical Validation

### H-2 reference tape comparison

The test at `test/runtests.jl:1532-1604` performs a comparison against the
reference tape at `njoy-reference/tests/84/referenceTape100`. Key observations:

1. **Grid density**: The reference tape has 769 points for MT=1 (total), MT=2
   (elastic), and MT=102 (capture). The Julia H-2 reconstruction produces
   a similar number of points (test checks `length > 100`).

2. **Tolerance**: The test uses a 5% relative tolerance (`rel_diff < 0.05`)
   and requires 70% of points to match (`n_close / n_compared > 0.7`).
   These are very loose acceptance criteria -- they could mask systematic
   errors of up to 5%.

3. **Missing reactions**: The reference tape includes MT=16 (n,2n) with
   125 points (a threshold reaction). The Julia reconstruction does NOT
   output MT=16 because `mt_list = [1, 2, 18, 102]` is hardcoded. This
   is a known simplification.

4. **Missing MF12**: The reference tape includes MF12/MT102 (photon
   production multiplicity). The Julia PENDF writer does not handle MF12.

### Total = elastic + fission + capture

- Tested explicitly in `runtests.jl:1508-1513` with `rtol=1e-8`.
- The `merge_background!` function (line 336-339) explicitly recomputes
  total as elastic + fission + capture, ensuring this holds by construction.
- The `merge_background_legacy` function (line 675) does the same.

### Grid density near resonances

- H-2 has no sharp resonances (scattering radius only), so this test
  (`runtests.jl:1612-1624`) only verifies basic grid density at low energies.
- For materials with actual resonance parameters (e.g., Fe-56), there are
  no tests at all. This is a significant coverage gap.

---

## Julia Idiom Issues

### File length violations

- `reconr.jl`: 679 lines (limit: 300) -- **VIOLATION**
- `pendf_writer.jl`: 380 lines (limit: 300) -- **VIOLATION**
- `adaptive_grid.jl`: 355 lines (limit: 300) -- **MARGINAL VIOLATION**
- See M5 for splitting recommendations.

### Type instability in hot loops

- `adaptive_grid.jl:252`: `mid_val = _as_tuple(f(xm))` -- the `_as_tuple`
  dispatch is `@inline` and type-stable for `NTuple{N,Float64}` and
  `CrossSections`. For the generic fallback (line 352-354), `ntuple` with
  a computed `N` is NOT type-stable because `N` is a runtime value. This
  could cause type instability if `f` returns a generic container.
  - However, in practice, the probe call at line 148-149 determines `N` at
    construction time, and subsequent calls should return the same type.
  - The `AdaptiveWorkspace{N}` is correctly parameterized on `N`.

### Allocations in inner loops

- `adaptive_grid.jl:315-316`: `push!(out_e, sx[i])` / `push!(out_v, sv[i])`
  -- these are append operations to pre-allocated (sizehinted) vectors.
  Occasional reallocation when capacity is exceeded, but this is standard
  and acceptable.
- `pendf_writer.jl:179-184`: Building `data` array with `push!` in a loop
  is acceptable given `sizehint!`.
- `reconr.jl:597-601`: Creating `CrossSections` objects in a loop (legacy
  interface) involves unnecessary allocation. The matrix-based `reconstruct`
  avoids this.

### try-catch in sigma_mf2

- `reconr.jl:629-636`: Bare `try ... catch` with no body in the catch block.
  This is a **critical** anti-pattern (see C1). It silently ignores all
  errors, including programming bugs, type errors, and bounds-check failures.
  At minimum this should be `catch e; @warn "..." e` or catch only specific
  expected exception types.

### Global mutable state

- No global mutable state detected. `AdaptiveConfig` is immutable. The
  `_STEP_RATIO` constant is properly declared as `const`. The `Ref{Int}`
  for sequence numbers in the PENDF writer is local to each function call.
  This is well done.

### Duplicate code

- `reconr.jl` contains two parallel implementations: `reconstruct` (modern,
  matrix-based) and `reconr` (legacy, NamedTuple-based). These share
  significant logic (~60% overlap). The legacy `reconr` also has its own
  `merge_background_legacy` that duplicates `merge_background!`. This
  duplication is a maintenance burden -- the legacy interface should delegate
  to the modern one and adapt the result format.

---

## AD Compatibility

### build_evaluator closure

- `reconr.jl:94-112`: The evaluator closure captures `sections` (a vector of
  tuples) and the keyword arguments. The inner computation is pure arithmetic
  (`+=`, `max`, `*`). No mutation of captured state.
- The `cross_section` call goes through dispatch to `cross_section_slbw`,
  `cross_section_mlbw`, or `cross_section_rm`. If those functions are
  AD-compatible, the closure is AD-compatible.
- **Concern**: The `@inbounds` annotation and the `max(0.0, x)` call.
  `max` is differentiable (subgradient at 0), and `@inbounds` should not
  affect AD. This looks AD-safe.

### Adaptive grid inner loop

- `adaptive_grid.jl:232-336`: The inner loop mutates `ws.stack_x` and
  `ws.stack_vals` (the pre-allocated workspace), and appends to `out_e`
  and `out_v`. These mutations are in the DRIVER code, not in the function
  `f` being evaluated.
- The function `f` is called at line 252 as a pure evaluation:
  `mid_val = _as_tuple(f(xm))`. This is the point where AD would operate.
- **Verdict**: AD through `f` is compatible. AD through the adaptive grid
  algorithm itself (differentiating the grid selection) is not possible due
  to mutation and control flow, but this is expected and acceptable.

### round_sigfig

- `adaptive_grid.jl:70-84`: Uses `log10`, `floor`, `round`, `abs`, integer
  arithmetic, and `exp10`. The `floor` and `round(Int, ...)` calls are not
  differentiable, but `round_sigfig` is only used in the grid selection
  logic, not in the cross section evaluation path. AD-safe by design.

---

## Test Coverage Gaps

### Convergence tiers

| Tier | Tested? | Details |
|------|---------|---------|
| Primary (err) | YES | Lorentzian test verifies convergence within tolerance |
| Relaxed + integral (errmax/errint) | INDIRECT | Lorentzian test allows `max_rel_error < 0.02` which is within errmax range |
| Step-ratio guard | NO | No test specifically triggers or verifies the 4.1x step-ratio guard |
| Thermal tightening | NO | No test verifies that tolerance tightens below 0.5 eV |

### PENDF output

| Feature | Tested? |
|---------|---------|
| Basic output (non-empty, has MF1/MF3) | YES |
| 80-char line length | YES |
| SEND/FEND/MEND/TEND correctness | NO |
| Dictionary correctness | NO |
| Roundtrip (write then read back) | NO |
| Comparison against reference tape | NO (only cross section values compared) |
| MF3 data values correct | NO |
| PointwiseMaterial interface | NO (only legacy NamedTuple tested) |

### Edge cases

| Case | Tested? |
|------|---------|
| H-2 (LRU=0, no resonances) | YES |
| Material with actual resonances | NO |
| Overlapping resonance ranges | NO |
| Very narrow resonances | NO |
| MF3 below threshold | NO |
| Empty MF3 sections | NO |
| Multiple isotopes | NO |
| Energy at range boundary | NO |
| E exactly equal to EH (exclusive) | NO |

### Missing test scenarios

1. **No test for any material with LRU=1 resonances** (e.g., Fe-56 SLBW,
   U-235 MLBW, Fe-56 RM). All RECONR tests use H-2 which has LRU=0 only.
   This means the adaptive reconstruction of resonance cross sections has
   ZERO end-to-end test coverage.

2. **No test for the `reconstruct` function** -- only `reconr` (legacy) is
   tested. The modern interface is completely untested.

3. **No test for `build_evaluator` with actual resonances** -- the test at
   line 1725-1728 evaluates at E=1.0 for H-2 which has no resonances,
   so it just returns zeros.

4. **No test for PENDF PointwiseMaterial writer** -- only the legacy
   NamedTuple interface is tested.

---

## Recommendations

### Must-fix before merge (blocking)

1. **C1**: Replace bare `try...catch` with specific exception handling and
   logging in `sigma_mf2` and `read_mf3_sections`.

2. **C2**: Add the `xm < eresr` (or equivalent upper-energy bound) condition
   to the step-ratio guard, either by adding a parameter to `AdaptiveConfig`
   or by ensuring the grid passed to `adaptive_reconstruct` is limited to
   the resolved region.

3. **M5**: Split `reconr.jl` and `pendf_writer.jl` to meet the 300-line limit.

### Should-fix soon (non-blocking but important)

4. **m1**: Fix `_write_fend` to not write the extra SEND record.

5. Add at least one test with a material that has actual resonance parameters
   (LRU=1, e.g., a simple SLBW material).

6. Add a test for the step-ratio guard and thermal tightening.

7. Add a PENDF roundtrip test (write then parse back to verify structure).

8. Eliminate the `reconr` / `reconstruct` duplication by having the legacy
   interface delegate to the modern one.

### Future work (non-blocking)

9. Implement `lunion`-style 1/v grid tightening and decade-point forcing
   for MF3 linearization.

10. Add unresolved resonance region support (`sigunr` equivalent).

11. Handle additional resonance MTs for LRF=7 evaluations.

12. Add threshold handling in `merge_background!`.

13. Tighten the H-2 reference comparison tolerance from 5% to 1% or better.

14. Implement the `_as_tuple` generic fallback using `Val{N}` to ensure
    type stability.
