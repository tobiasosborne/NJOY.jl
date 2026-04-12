# T03 Completion + T04 Progress — Phase 7 Handoff

## Date: 2026-04-12

## Summary

**T03 DONE.** tape34, tape36, tape37 are **100% BIT-IDENTICAL**. Official
`execute.py` test passes at `rel_tol=1e-12` (much tighter than the 1e-9 default).

**T04 substantial progress.** tape23 passes at `rel_tol=1e-7`, tape24 passes
at `rel_tol=1e-5`, tape25 has 11 residual value diffs at 1e-5.

Commits: `73c21e8` (T03), `6f50e09` (T04).

---

## T03 — what was fixed

5 root causes diagnosed via `write(99, ...)` diagnostics in Fortran
`gpanel` / `gtff` and per-panel comparison to Julia.

### 1. Last-panel linear-rate extrapolation (gaminr.f90:980-982)

Fortran's `rr = b + (a-b)*(eq-elo)/(ehi_local-elo)` uses the **original
group upper bound** on the last panel (`ehigh = delta*ehi < ehi`), not
the delta-nudged `ehigh`. Julia was doing standard trap with σW evaluated
at the nudged boundary. Applied the extrapolation formula in:

- `_gaminr_group_average`
- `_gaminr_heating_average`
- `_gaminr_coherent_matrix!`
- `_gaminr_incoherent_matrix!`
- `_gaminr_pair_production_matrix!`

```julia
if is_last_of_group
    t1_pb = dp / (ehi - pa)
    rr_pb = b + (a - b) * t1_pb
    num += (b + rr_pb) * dp / 2.0
else
    num += (b + a) * dp / 2.0   # standard trap, unchanged for middle panels
end
```

### 2. MF26 coherent: GL-6 → GL-2

`gtff MT=502` returns `nq=0` → `gpanel` uses `nq=2` (Lobatto-2 =
trapezoidal). Julia was using GL-6 for "accuracy"; that's not faithful.

### 3. MF26 incoherent: GL-6 → Lobatto-10

Discovered via diagnostic: `gtff MT=504` returns `nq=8` (from
`nq=nq+2` at gaminr.f90:1495 where the inner `nq=6` for `nl<=6` gets
bumped). `gpanel` then adds another 2 → 10 (capped at 10). So
**always Lobatto-10** for MF26/MT=504.

### 4. Kinematic Compton breakpoints (gaminr.f90:1477-1495)

`gtff MT=504` returns `enext = egg[i]/(1-2·c3·egg[i])` for `egg[i]<1/c3`.
These are Compton-edge energies that split panels inside groups where
`c3·elo < 0.5`. Added `kin_breaks` list to `_gaminr_incoherent_matrix!`
with **discontinuity handling** (delta on approach, rndoff after —
matching Fortran's `idisc=1` path at gaminr.f90:935).

### 5. Per-MT flux for MT=504

MT=504's flux integral differs from MT=501's base flux because of the
kinematic panel splits. `_gaminr_scatter_matrix` now returns
`(ans, scat_flux)`; caller uses per-MT flux as the MF26 flux column and
toth denominator for MT=504.

### Results

| Tape | Before | After |
|---|---|---|
| tape33 (GENDF) | 48.7% | 97.0% (only ULP diffs remain) |
| tape34 (DTF) | 86.5% | **100% BIT-IDENTICAL** |
| tape36 (plot) | 98.4% | **100% BIT-IDENTICAL** |
| tape37 (PS) | 32.4% | **100% BIT-IDENTICAL** |

---

## T04 — what was fixed

5 fixes against referenceTape23/24/25.

### 1. errorr MF1/MT451 N1 field

Fortran `errorr.f90:5927` writes `b(5) = -11` — a hardcoded marker, not
a computed `-NW`. Julia was writing `-length(egn)` which gave `-10`.

### 2. errorr 1st-call sequence resets

Standard ENDF: `seq` resets to 1 at each SEND (new section within same
MF) and FEND (new MF). Julia was writing continuous `seq` through the
entire file.

### 3. errorr 2nd-call tape25 hybrid convention

`tape25` uses a weird hybrid: TPID `seq=0` + **continuous** numbering
through Phase A (copied nin content + MEND separator + 2nd material's
MF1/MT451), then **switches to standard ENDF** (`SEND=99999`, `FEND=0`,
`seq` resets per section) for the 2nd material's MF3/MF33 sections.

Fortran's `errorr.f90` 2nd-call writer produces this pattern; Julia now
matches it.

### 4. groupr tape24 structure quirks

- MF1/MT451 goes **straight to FEND** (no SEND).
- MF3 sections use `SEND=99999` but **SKIP FEND** — MEND follows
  directly.
- TEND at the very end (Julia was missing it).

### 5. MF31 LB=1/2 covariance expansion

Julia's 2nd errorr call was **skipping LB=1 and LB=2 blocks** entirely
(only LB=5 was processed), producing covariance values ~20× smaller
than Fortran's. Fix:

- `expand_lb1` (LB=0/1): diagonal, unchanged.
- **New `expand_lb2`**: fully-correlated outer product
  `C[i,j] = F_k × F_l`, using group **lower bound** for the interval
  lookup (matches Fortran covcal `un(jg)`/`un(jh)` at
  errorr.f90:2087-2090).
- `expand_covariance_block` dispatches LB=2 to `expand_lb2`.
- 2nd call refactored to build a `CovarianceBlock` list and call
  `multigroup_covariance(blocks, egn2)`.

### Results

| Tape | Before | After byte-match | execute.py pass |
|---|---|---|---|
| tape23 | 17.1% | 97.6% | **rel_tol=1e-7** |
| tape24 | 73.0% | 75.7% | **rel_tol=1e-5** |
| tape25 | 15.1% | 84.9% | 11 bad at 1e-5 |

---

## T04 residual — tape25

**11 lines still differ at `rel_tol=1e-5`**, all in MF33/MT452
covariance matrix entries (rows 1-7, cols 6-7). Diffs are 3-6% relative.

**Root cause:** Fortran errorr's `covcal` works on a **union grid**
(all block breakpoints ∪ group boundaries), computing the covariance
per fine interval and then collapsing to the target group grid with a
`T * C_fine * T'` projection. Julia's `multigroup_covariance` function
does this pipeline, but only when `total_e <= 10000`; for U-235 MF31
the blocks have breakpoints at fine energies (4e6, 6e6, 8e6, 1.2e7,
1.6e7, etc.) that straddle group 7 `[1e6, 1e7]`. The collapse math
doesn't quite match Fortran's per-entry convention for outer-product
(LB=2) blocks whose F intervals are finer than the target groups.

**Not yet diagnosed:** whether Fortran applies a flux-weighted collapse
(`sig(jg)*sig1(jh)` weighting at errorr.f90:2283/2294) that Julia's
direct expansion misses, or whether there's a sigfig rounding at a
specific step.

**Impact:** Low — tape25 passes structurally (same line count, same
section layout) and byte-matches 85%. The 11 bad lines cluster in the
nubar-nubar covariance corner of the matrix.

---

## Files changed (T03 + T04)

| File | Changes |
|---|---|
| `src/orchestration/modules/gaminr.jl` | +162/-85 — 5 last-panel fixes, GL-2 coherent, Lobatto-10 incoherent, kin_breaks, per-MT flux |
| `src/orchestration/modules/errorr.jl` | +38/-45 — MF1 N1=-11, seq resets, 2nd-call hybrid convention |
| `src/orchestration/modules/groupr.jl` | +4/-8 — remove stray SEND, add TEND |
| `src/processing/errorr.jl` | +25/-5 — `expand_lb2`, dispatch LB=2 properly |

---

## How to verify

```bash
# T03 — all 4 tapes at 100% (or ULP for tape33)
rm -rf /tmp/t03_test ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/03/input"; work_dir="/tmp/t03_test")'
# Compare: tape34/36/37 should be 100% byte-identical with /tmp/t03_fortran equivalents

# T04 — tape23 passes at 1e-7, tape24 at 1e-5
rm -rf /tmp/t04_test ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/04/input"; work_dir="/tmp/t04_test")'
# referenceTape23/24 comparison via execute.py-style float matching

# Regression: T01, T02 unchanged
```

---

## Next steps

### Priority 1 — T04 tape25 covariance residual (MEDIUM IMPACT)

Diagnose why the union-grid collapse for MF31 LB=2 blocks doesn't match
Fortran's covcal output for cols 6-7 of the 7×7 matrix. Likely fix is in
`multigroup_covariance` or its `_collapse_matrix` weighting. Compare
Fortran's `covcal` output directly via a diagnostic printing
`cov(jg, jh)` before the final sigfig.

### Priority 2 — T03 tape33 ULP diffs (LOW IMPACT)

18 remaining diffs in tape33 are ULP-level (1-ULP in last sigfig) for
MT=621 heating and some MF26 entries. Doesn't affect tape34/36/37.
Would need FP accumulation order tweaking in
`_gaminr_incoherent_matrix!` heating accumulation.

### Priority 3 — Other tests (T05+)

T03 and T04 both now run end-to-end. T01 (32750/32962) and T02
(86.8% at 1e-5) are unchanged. Next candidates: T05, T06, T08+ —
generate oracles if missing, compare per-module.
