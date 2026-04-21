# T15 covcal MT=77 diagnosis (NJOY.jl-f8k canary)

## Date: 2026-04-21

## Status

**Diagnosis only**, no implementation yet. Root cause pinned down via
Fortran `write(*,...)` diagnostic trace + Julia code reading.

## Grind Method walkthrough

### 1. Pick canary MT — MT=77 (n,n7' inelastic, U-238 JENDL-3.3)

Per HANDOFF Open Work P1 "Covcal content drift": Julia 30 lines vs ref 20,
values differ. Smallest fully-divergent MT in T15 tape26.

### 2. Generate oracle — `njoy-reference/tests/15/referenceTape26` (already exists)

### 3. Run Julia, save tape26 to `/tmp/t15_grind/tape26`

Phase-49 errorr_module run produced tape26 with MAT=9237 MF=33 MT=77.

### 4. Find first diff

`/tmp/t15_grind/ref_mt77.txt` vs `/tmp/t15_grind/jul_mt77.txt`:

| | Reference | Julia |
|---|---|---|
| HEAD record | `0 0 0 0 0 3` (NK=3 sub-sections) | `0 0 0 0 0 1` (NK=1) |
| (mt2=77) self LB-block | LB=5 LT=1, 5 nonzero rows + 1 zero stub = 14 lines | LB=5 LT=1, 11 nonzero rows = 28 lines |
| (mt2=91) cross | all-zero stub, 3 lines | **missing** |
| (mt2=102) cross | all-zero stub, 3 lines | **missing** |
| **First differing data line** | row 20: `2.987998-2 3.150980-2 2.611777-2 2.012270-2 7.056241-3` | `4.000000-2 3.600000-2 3.200000-2 2.600000-2 4.000000-3` |

Two distinct bugs:

- **Bug A — missing zero-stub cross-pairs.** Julia's writer emits only
  pairs in `cov_matrices` ∪ `listed_pairs` (plus `nc_derived_mts`-style
  synthesis); MT=77 has neither (91, 102) in `listed_pairs` nor
  zero-cov entries in `cov_matrices`. Reference Fortran covout writes a
  zero-stub for every reaction MT pair declared in `mfcov`'s reaction
  list with `mt2 > mt`.

- **Bug B (the canary itself) — wrong values + extra rows.** Julia's
  `expand_lb5_symmetric` (src/processing/errorr.jl:70-82) does
  midpoint-sampling: for output group (g_i, g_j) it copies the single
  input bin value at the midpoint with no group-averaging. This is why
  the values are "round" (`4.0e-2`, `3.6e-2`, etc.) — they're verbatim
  from MT=77's LB=5 block. The Fortran does proper XS·flux-weighted
  expansion on a UNION grid + collapse to output groups.

### 5. Classify diff

- Bug A: line-count + sub-section count diff. Easy fix in writer (auto-
  synthesize zero-stub pairs for cross-MTs in `reaction_mts`).
- Bug B: structural / value diff. Requires architectural change to
  errorr_module covariance pipeline (union-grid + XS + flux weighting).

### 6. Read Fortran covcal — confirm algorithm

`njoy-reference/src/errorr.f90` `covcal` subroutine (lines 1770-2417):

LB=5 branch (lines 2208-2235):
```fortran
! Find input bin k containing union-group jg lower bound `eg`
350 k = k+1
   if (eg.ge.scr(locip6+k-1).and.eg.lt.scr(locip6+k)) go to 360
   if (k.lt.nk1) go to 350
   go to 510
! Find input bin l containing union-group jh lower bound `ehr`
360 l = 0
370 l = l+1
   ...
390 ifloc = locip6 + nk1*np/2 - (np-l+1)*(np-l)/2 + k - l   ! upper-tri index
   if (l.ge.k) ifloc = locip6 + nk1*np/2 - (np-k+1)*(np-k)/2 + l - k
400 cov(jh) = cov(jh) + scr(ifloc+np)*sig(jg)*sig1(jh)      ! ← XS weighted
```

Per-row write (lines 2336-2353):
```fortran
do ij=1,jgend
   ic=ic+1
   ...
   b(ic) = cov(ij)*flx(jg)*flx(ij)                           ! ← flux weighted
   ...
enddo
```

Net per-(jg, jh) on union grid:
```
b[jg, jh] = LB5_fvals[k(jg), k(jh)] · sig(jg) · sig(jh) · flx(jg) · flx(jh)
```

Then covout collapses to output groups by summation over (jg→ig, jh→igp).
Final relative cov:
```
C_out_rel[ig, igp] = Σ b[jg, jh] / (sig_out[ig]·flx_out[ig] · sig_out[igp]·flx_out[igp])
```

For piecewise-constant XS within an output group, this collapses to
flux-weighted averaging of the input cov values within the (output_group_ig
× output_group_igp) rectangle.

### 7. gdb-equivalent: Fortran `write(*,...)` trace

Patched `errorr.f90` covcal LB=5 branch (line 2230) with:
```fortran
if (mts(nmts).eq.77 .and. un(jg).ge.1.35e6 .and. un(jg).lt.1.74e6) then
   write(*,'...') jg, jh, un(jg), un(jh), k, l, scr(ifloc+np), &
      sig(jg), sig1(jh), scr(ifloc+np)*sig(jg)*sig1(jh)
endif
```

Recompiled `njoy-reference/build/njoy`, ran T15 input deck on `/tmp/t15_fortran_diag`.
Yielded 990 trace lines for union groups 56-64 covering LANL-30 group 20.

#### Key trace insights

For LANL group 20 = [1.353e6, 1.738e6], the union grid contains 9
sub-groups (jg ∈ {56..64}) covering this range — **much finer than
LANL alone**. Union grid edges include LANL boundaries + input MT=77
bin edges + extras from other reaction grids and the GENDF processing.

XS varies dramatically across union groups within LANL group 20:

| un(jg) | sig(jg) | input bin k |
|---|---|---|
| 1.350e+06 | 7.461067e-05 | 1 |
| 1.400e+06 | 5.255324e-03 | 2 |
| 1.450e+06 | 5.255324e-03 | 2 |
| 1.500e+06 | 5.255324e-03 | 3 |
| 1.550e+06 | 5.255324e-03 | 3 |
| 1.650e+06 | 5.255324e-03 | 3 |
| 1.700e+06 | 5.255324e-03 | 4 |
| 1.738e+06 | 8.133546e-03 | 4 |

**Bin 1 is sub-threshold** (XS = 7.46e-5, near zero) and contributes
zero to the LB=5 fvals (also zero in input). The other bins have:
- bin 2 (`[1.4e6, 1.5e6]`): fvals (k=2, l=2) = 4.0e-2
- bin 3 (`[1.5e6, 1.7e6]`): fvals (k=3, l=3) = 3.8e-2
- bin 4 (`[1.7e6, 1.9e6]`): fvals (k=4, l=4) = 3.6e-2

**Reference C_out[20, 20] = 2.987998e-2** is the XS·flux-weighted
average of these (4.0e-2, 3.8e-2, 3.6e-2) over union groups in LANL
group 20, with bin 1's zero contribution dragging the average down.

**Julia C_out[20, 20] = 4.0e-2** = midpoint-bin sample (group 20
midpoint 1.5455e6 falls in input bin 3? actually it should be 3.8e-2
from bin 3 — but Julia outputs 4.0e-2 = bin 2's value. Possible bug in
`_find_bin` boundary handling, but moot until the union-grid
reformulation lands).

Trace file: `/tmp/t15_fortran_diag/stdout.log` (990 DIAG lines).

### 8. Fortran patch reverted

`njoy-reference/src/errorr.f90` reverted via `git checkout`; binary
rebuilt clean. No artifact left in `njoy-reference/`.

## Required fix (next session)

Substantial refactor of `src/orchestration/modules/errorr.jl` and/or
`src/processing/errorr.jl`:

1. **Union-grid expansion**: build `ugrid = sort(unique(egn ∪ all_block.energies))`.
2. **Per-union-group XS**: derive from `group_xs` (output XS, available)
   by piecewise-constant assignment per output group containing each
   union group. Closer-to-Fortran would interpolate from PENDF, but
   piecewise-constant should be a strong improvement.
3. **Per-union-group flux**: Fortran covcal reads union-group fluxes
   from groupr's GENDF tape (`flx` array). For Julia first-pass use
   lethargy weights `flx_u[g] = log(ugrid[g+1] / ugrid[g])` or bin
   width `ugrid[g+1] - ugrid[g]`. Should add real flux read in a
   later iteration.
4. **Sandwich**:
   ```julia
   b[g, h] = LB5_fvals[k(g), k(h)] * sig_u[g] * sig_u[h] * flx_u[g] * flx_u[h]
   ```
5. **Collapse to output**:
   ```julia
   for ig in 1:ngn, igp in 1:ngn
     num = sum(b[g, h] for g→ig, h→igp)
     den = (Σ sig_u[g]·flx_u[g] for g→ig) · (Σ sig_u[h]·flx_u[h] for h→igp)
     C_out[ig, igp] = num / den
   end
   ```
6. **Bug A fix (writer)**: when iterating sub_keys for an MT in
   `reaction_mts`, auto-synthesize zero-stub entries for every
   `mt2 > mt in reaction_mts` whose pair was declared in the input
   ENDF (`listed_pairs`) — already partially handled for NC-derived
   MTs; extend to non-NC MTs whose input has multiple sub-sections in
   `listed_pairs`. Need to verify against reference how Fortran
   chooses which pairs to emit.

### Test design (RED→GREEN)

Add to `test/validation/test_errorr_mf33_sparse.jl`:
- Element-wise check `|jul_mt77[20, 20] − 2.987998e-2| < 1e-6`.
- `jul_pair_lines[(77, 91)] >= 2` (zero-stub presence).
- `jul_pair_lines[(77, 102)] >= 2`.
- After fix expect MT=77 line count to drop from 30 → ~20.

## Files referenced

- `src/processing/errorr.jl:70-82` — `expand_lb5_symmetric` (the bug)
- `src/orchestration/modules/errorr.jl:168` — caller passes `egn`
  directly (no union grid)
- `src/orchestration/modules/errorr.jl:1002-1011` — writer pair
  iteration (Bug A surface)
- `njoy-reference/src/errorr.f90:1770-2417` — `covcal` subroutine
- `njoy-reference/src/errorr.f90:2208-2235` — LB=5 branch
- `njoy-reference/src/errorr.f90:2336-2353` — per-row write with flux
  multiplication
- `njoy-reference/tests/resources/J33U238` — input ENDF (MT=77 LB=5
  starts at line 537 of MF=33 MT=77 in the per-section count)
- `njoy-reference/tests/15/referenceTape26` — reference output
- `/tmp/t15_fortran_diag/stdout.log` — Fortran trace (990 DIAG lines)
- `/tmp/t15_grind/{ref,jul}_mt77.txt` — extracted MT=77 sub-sections
