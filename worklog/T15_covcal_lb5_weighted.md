# T15 covcal LB=5 σ·flx-weighted collapse (NJOY.jl-f8k Bug B)

## Date: 2026-04-22

## Outcome

MT=77 self-cov C[20, 20] now matches reference **exactly** (0.02987998
vs 0.02987998, max diff on the 5×5 nonzero block [20:24, 20:24] = 0.0).

T15 tape26 lines: 4240 → 4275 (MT=77 content changed, +35 net across
all MTs). T15 MF33 test still at 36/36 pass; T04 MF31 unchanged
(107/119 lines match on tape25, identical to pre-change baseline).

## RED → GREEN

### RED test

`test/validation/test_errorr_covcal_lb5.jl`: runs T15 groupr + second
errorr (MF33 path) into a temp dir, parses the MT=77 self-cov
sub-section back from `tape26` into a 30×30 matrix using
`NJOY.find_section` + `read_cont` + `read_list`, and asserts

```
|C[20, 20] − 2.987998e-2| < 1e-6      # canary from Phase-50 trace
max |jul − ref| on [20:24, 20:24] < 1e-5
```

Pre-fix reading: jul = 0.04 — direct midpoint sample of a single
LB=5 input bin by `expand_lb5_symmetric`, bypassing the union-grid +
σ·flx-weighted collapse path the Fortran uses.

### GREEN — four pieces

1. **`src/processing/errorr.jl` — `_collapse_matrix`** gained a
   `weights::Union{Nothing,AbstractVector{<:Real}}=nothing` kwarg. When
   given, each row `T[g, k]` is `weights[k]·overlap/u_w` divided by
   `Σ weights[k']·overlap/u_w` over union cells overlapping group `g`.
   The legacy `flux` kwarg is kept (falls back to `weights` semantics
   when only `flux` is provided).

2. **`src/processing/errorr.jl` — `multigroup_covariance`** gained
   `xs_row`, `xs_col`, and `ugrid` kwargs. When either xs is given,
   `T_row = _collapse_matrix(..., weights=xs_row·flux)` and
   `T_col = _collapse_matrix(..., weights=xs_col·flux)` are built
   independently, and `C_out = T_row · C_fine · T_col'`. Mirrors
   Fortran covcal `b[jg, jh] = fvals · sig(jg)·sig1(jh) · flx(jg)·flx(jh)`
   (errorr.f90:2208-2235 LB=5 accumulate, 2336-2353 per-row write)
   followed by covout's union→output collapse (errorr.f90:7431-7438).

3. **`src/orchestration/modules/errorr.jl` — block collection**.
   Replaced the inline `C = expand_covariance_block(block, egn)` +
   accumulate pattern with a two-pass scheme: first pass builds
   `pair_blocks :: Dict{(mt, mt2) => Vector{CovarianceBlock}}`; second
   pass (`_collapse_pair_blocks!`) routes each pair's LB=5/6 blocks
   through `multigroup_covariance` with:

   - `ugrid` = `egn ∪ all block.energies`, clipped to `egn`'s range.
     Each union cell lies entirely within one output (LANL) group and
     one LB=5 input bin, so midpoint-sampling produces piecewise-
     constant fvals per cell — the correct input to a weighted
     collapse.
   - `xs_row[k]`, `xs_col[k]` = output-group-averaged σ keyed by the
     LANL group containing `ugrid[k]` (the union cell lower bound).
     Matches Fortran `sig(jg)` semantics: `rdsig`/`rdgout` populates
     σ per GENDF group and covcal uses `eg = un(jg)` (lower bound)
     for lookups.
   - `flx_u[k]` = `ugrid[k+1] − ugrid[k]` (bin width). For T15's
     `iwt=2` (flat weighting), the GENDF group flux is proportional
     to the group's energy width, so bin width per union cell matches
     Fortran `flx(jg)` up to an overall constant that cancels in the
     normalised collapse.

   LB=0/1/2 blocks stay on the legacy direct-expansion path on `egn`
   (no change in behaviour for T04 MF31, other nubar/angle paths).

4. **RED test file** added under `test/validation/test_errorr_covcal_lb5.jl`.

## Ground-truth derivation (why bin width, not lethargy)

The canary was the 3% gap between Julia-with-lethargy (0.02891) and
Fortran (0.02988) after the first union-grid + σ·flx pass landed.

Raw LB=5 data for MT=77 (from `njoy-reference/tests/resources/J33U238`,
MAT=9237 MF=33 MT=77):

- `ek = [1e-5, 1.4e6, 1.5e6, 1.7e6, 1.9e6, ...]` (25 breakpoints, 24 bins)
- `fvals(1, *) = 0` (sub-threshold row)
- Diagonals `fvals(k, k) = 0.04` for `k ≥ 2`
- Off-diagonals: `fvals(k, k+1) = 0.038`, `fvals(k, k+2) = 0.036`, …

LANL-30 group 20 = `[1.353e6, 1.738e6]`. Four MT=77-only union cells:

| union cell          | bin | fvals(k, k) | width  |
|---------------------|-----|-------------|--------|
| [1.353e6, 1.4e6]    | 1   | 0           | 0.047  |
| [1.4e6, 1.5e6]      | 2   | 0.04        | 0.100  |
| [1.5e6, 1.7e6]      | 3   | 0.04        | 0.200  |
| [1.7e6, 1.738e6]    | 4   | 0.04        | 0.038  |

With weight `w_k = u_width` (σ constant within LANL 20 cancels), the
4×4 σ·flx-weighted σum over `(jg, jh)` with fvals:

```
Σ fvals(k_jg, k_jh) · w_jg · w_jh
= fvals(2,2)·w²_2 + 2·fvals(2,3)·w_2·w_3 + 2·fvals(2,4)·w_2·w_4
+ fvals(3,3)·w²_3 + 2·fvals(3,4)·w_3·w_4 + fvals(4,4)·w²_4
(row 1 contributes 0)
= 0.04·0.01 + 2·0.038·0.02 + 2·0.036·0.00376
+ 0.04·0.04 + 2·0.038·0.0076 + 0.04·0.001444
= 4e-4 + 1.52e-3 + 2.736e-4 + 1.6e-3 + 5.776e-4 + 5.776e-5
= 4.4296e-3
```

Normalised by `(Σ w)² = 0.385² = 0.148225`:
`C_out[20, 20] = 4.4296e-3 / 0.148225 = 0.02989` — matches 0.02987998
to within the 5th significant figure.

Lethargy weights `log(u_hi/u_lo)` give a different distribution across
cells (0.034, 0.069, 0.125, 0.022) → weighted σum ≈ 0.0334 / (Σ²) ≈
different result. The GENDF flux for `iwt=2` is width-proportional (the
T15 GENDF MT=3 MT=77 flux values 5.30e5, 3.85e5, 4.94e5, 6.33e5 for
groups 19, 20, 21, 22 precisely match the LANL group widths 0.530e6,
0.385e6, 0.494e6, 0.633e6). So bin width is correct for iwt=2.

For iwt=3 (1/E) the analogous relation would pick out lethargy. A
follow-up should read real per-union-cell flux from the input GENDF
(MT=1 MF=3 records; Fortran rdgout at errorr.f90:2788-2927) for full
iwt-agnostic fidelity. Not blocking any current test.

## What this does NOT fix (still open)

- **Bug A — writer zero-stub cross-pairs**: T15 MT=77 still emits NK=1
  (self only) while reference has NK=3 (self + zero-stub (77, 91) +
  zero-stub (77, 102)). Values are right; structure count differs.
  Worklog/HANDOFF P1 has the fix plan; separate phase.
- **T15 tape26 total lines**: 4275 vs ref 5958. Closing the remaining
  gap is largely the NK=3 zero-stub count (Bug A) plus a handful of
  smaller cross-MT cases not yet pinned.
- **LB=1/2 blocks**: still go through `expand_lb1`/`expand_lb2` directly
  on `egn`. Same class of bug in theory (midpoint sampling loses
  within-group XS structure), but no canary has surfaced for them.
  T04 tape25 MF31 residual 11-of-119 diffs (HANDOFF P2) may be this.

## Files touched

- `src/processing/errorr.jl` — `multigroup_covariance` +
  `_collapse_matrix` kwargs.
- `src/orchestration/modules/errorr.jl` — block collection +
  `_collapse_pair_blocks!` + `_find_output_group` helpers.
- `test/validation/test_errorr_covcal_lb5.jl` — RED→GREEN test (new).

## Regression check

- `test/validation/test_errorr_mf33_sparse.jl`: 36/36 pass.
- API unit test (xs_row/xs_col/flux backward-compat): 8/8 pass.
- Reference T04 (MF31, LB=2): tape24 NUMERIC_PASS 56/74, tape25 DIFFS
  107/119 — identical to pre-change baseline.

## Reference

- Phase-50 diagnosis: `worklog/T15_covcal_mt77_diagnosis.md`
- Fortran covcal: `njoy-reference/src/errorr.f90:1770-2417`
- Fortran covout pair iteration: `errorr.f90:7393-7440`
- LB=5 branch: `errorr.f90:2208-2235`
- Per-row σum-and-write: `errorr.f90:2336-2353`
