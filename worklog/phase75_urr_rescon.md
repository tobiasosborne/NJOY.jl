# Phase 75 — URR rescon port (T15 tape26 EXACT line-count match)

**Date**: 2026-05-19
**Branch**: master
**Status**: GREEN — T15 tape26 line count Julia **5958 = ref 5958** (gap −6 → **0**). Every MF=33 sub-section row layout matches reference byte-for-byte at the LIST-header level. Phase 72c/73/74 canaries preserved. T22 BIT_IDENTICAL 4636/4636 preserved.

## Summary

Phase 74 left T15 tape26 at −6 lines vs reference. The Phase 74 worklog hypothesised this was a `_resonance_group_window` boundary bug ("MT=102 row 14 / MT=2 rows 13-15 — group-window edge"). Diagnosis this session disproved that: `_resonance_group_window` faithfully mirrors Fortran resprx (errorr.f90:3093-3108). The real root cause was that **URR (LRU=2) rescon was not ported** — Phase 72 deferred it with a `@info ... skipping` and the gap had been waiting on it ever since.

U-238 JENDL MF=32 has 10 LRU=1 ranges (max EH = 10 keV) plus one LRU=2 range [10 keV, 150 keV]. The URR range spans LANL groups 13..15. Without the URR sandwich, those cells in (1,1), (2,2), (2,102), (102,102) sub-sections were below the writer's eps=1e-20 threshold and the rows were skipped. With the URR sandwich those cells rise above eps and the rows appear — closing the gap precisely.

## Fortran reference

The URR rescon chain in Fortran covout:

| Subroutine | errorr.f90 lines | Purpose |
|---|---|---|
| `rdumrd2` | 5091-5225 | Walk MF=2 to populate `amur(3, NLRU2)` = (AMUN, AMUF, AMUX) per (L, J). LRU=2/LRF=1/LFW=0 walks but leaves amur untouched; LRU=2/LRF=1/LFW=1 walks via LIST records; LRU=2/LRF=2 fills amur from LIST words 9, 10, 12 (lines 5206-5208). |
| `resprx` | 3066-3247 | Top-level rescon orchestrator. Computes `iest, ieed` per range (3093-3108) — the [first group above elr, last group ≥ ehr] window. Dispatches to rpxlc12 / rpxlc2 / rpxlc0 / rpxsamm / **rpxunr** by LRU/LRF/LCOMP. |
| `rpxunr` | 4785-4985 | URR sensitivity build + sandwich. Walks NLS L-state LISTs (each holding NJS J-states × 6 words = D, AJ, GNO, GG, GF, GX), then one final cov LIST (MPAR × NPAR, NPAR = MPAR · ΣNJS). Loops `loop = 1..NPAR+1`: loop=1 unperturbed, loop>1 perturbs one (J, param) by 1.01×. Steps E from `elg = egn(iest)` by ratio 1.015 with forced breakpoints at elr / ehr / ehg. Calls ggunr1, group-averages via rpxgrp, builds `sens(ch, p, ig) = 100·(σ_p - σ_0)/orig · cflx · abn` for ig ∈ [iest, ieed]. Sandwiches into uff/ugg/uee/utt/uef/ueg/ufg with absolute cov `cov(i,j) = cov_REL(i,j) · b(i) · b(j)`. |
| `ggunr1` | 6800-6905 | URR XS at one E. Walks L-states, per-J calls `eunfac` for penetrability + phase shift, accumulates potential scattering once per L, fluctuation integrals via `egnrl` for elastic/capture/fission. |

## Fix

Two files, ~310 LOC added net (zero LOC subtracted from existing paths).

### `src/processing/mf32_reader.jl`

Replaced the LRU=2 skip-stub (Phase 72's `@info "URR cov port deferred to next session"`) with a real reader. Added three structs (`MF32URRJState`, `MF32URRLState`, `MF32UnresolvedRange`) and extended `MF32IsotopeData` with `unresolved_ranges::Vector{MF32UnresolvedRange}`. The reader parses:

1. URR head CONT: `[SPI, AP, LSSF, 0, NLS, ISR]`.
2. If ISR=1: one CONT for global ΔAP.
3. NLS LIST records: each holds NJS J-states × 6 body words (D, AJ, GNO, GG, GF, GX).
4. Final cov LIST: L1=MPAR, N1=NW=NPAR(NPAR+1)/2, N2=NPAR. Body is the upper-triangular packed relative-relative cov matrix, unpacked into a full symmetric `cov_RP`.

Errors loudly on `length(body) ≠ 6·NJS` or `NW ≠ NPAR(NPAR+1)/2`. NPAR is checked against MPAR·ΣNJS by the rescon caller (one layer up — the reader doesn't know MPAR until it reads the cov LIST itself).

### `src/processing/rescon.jl`

New section (≈260 LOC) before `apply_rescon!`:

- `_urr_param_offsets(mpar, lfw)` — table of which offsets into `[D, AJ, GNO, GG, GF, GX]` are perturbed for each MPAR value. MPAR=1: D. MPAR=2: D, GNO. MPAR=3: D, GNO, GG (U-238). MPAR=4: depends on LFW. MPAR=5: D, GNO, GG, GF, GX. Mirrors Fortran rpxunr offsets at errorr.f90:4865-4897.
- `_flatten_urr_params(urr)` — flatten MF=32 L-states into per-J 6-vectors, plus L and AWRI tracking arrays.
- `_urr_dofs_from_mf2(rng)` — extract (AMUN, AMUF, AMUX) → (mu, nu, lamda) from URR2Data or URRData. Returns Vector{NTuple{3,Int}} in (L, J) order — must align with `_flatten_urr_params` (it does, per the (L, J) walk of MF=2 reader).
- `_ggunr1(E, params, l_per_j, awri_per_j, dofs, spi, ap)` — port of Fortran ggunr1. Single-E URR XS via SLBW + fluctuation integrals (`_unfac`, `_gnrl`, both pre-existing in `unresolved.jl`). Spot scattering accumulated once per new L (matches `_csunr1`'s `ll != prev_l` pattern). Returns (tot, el, fis, cap).
- `_build_urr_xs_grid(urr, params, ..., elg, ehg)` — dense E grid via ×1.015 stepping with elr/ehr/ehg breakpoints. Outside [elr, ehr]: XS=0. Mirrors errorr.f90:4904-4930.
- `_apply_rescon_urr_range!(cov_matrices, urr, mf2_urr_range, egn, group_xs, weight_fn)` — full sensitivity build + sandwich for one URR range. Returns the count of (mt, mt2) pair contributions applied.

`apply_rescon!` extended with a parallel loop after the existing LRU=1 loop, walking `iso.unresolved_ranges` and matching MF=2 URR ranges by [EL, EH] (MF=32 LRU=2 reports LRF=1 in its range head regardless of the MF=2 LRF — see ENDF-6 §32.2 — so matching by LRF is wrong).

### `src/NJOY.jl`

Added exports: `MF32UnresolvedRange`, `MF32URRLState`, `MF32URRJState`.

### `test/validation/test_mf32_reader.jl`

Added a 25-assertion testset `read_mf32 — U-238 URR range (10 keV..150 keV)` covering the head CONT (EL, EH, LRU, LRF, SPI, AP, LSSF, NLS), per-L-state shape (1+2+2 J-states for L=0/1/2), first-J body values (D=11.014, AJ=0.5, GNO=1.0194e-3, GG=2.3994e-2, GF=0, GX=0), and cov matrix shape (15×15) + symmetry + first diagonal value (1.5585e-2). Plus one assertion in the top-level testset for `length(iso.unresolved_ranges) == 1`.

## Before / after

| Quantity                                | Before    | After     | Reference |
|-----------------------------------------|-----------|-----------|-----------|
| T15 tape26 line count                   | 5952      | **5958**  | 5958      |
| Net gap vs ref                          | **−6**    | **0**     | 0         |
| T15 tape26 MF=33 line count             | 5649      | **5655**  | 5655      |
| (1, 1) row layout                       | 29 rows   | **30**    | 30        |
| (2, 2) row layout (cols 13-15 short)    | drift     | **exact** | exact     |
| (2, 102) row 14 row 14                  | absent    | **present**| present  |
| (102, 102) row 14                       | absent    | **present**| present  |
| RP-cov sandwich contributions logged    | 70        | **77**    | —         |
| MT=102 C[1,1] (Phase 72c canary)        | 2.658408e-4 | preserved | 2.658914e-4 |
| T22 BIT_IDENTICAL                       | 4636/4636 | 4636/4636 | —         |

## Test results

| Test file | PASS | Notes |
|---|---|---|
| `test_mf32_reader.jl` | 83/83 (+26 new) | 8 + 25 + 21 + 20 + 2 + 7 |
| `test_errorr_mf33_sparse.jl` | 75/75 | T15 tape26 = 5958 = ref |
| `test_errorr_covcal_lb5.jl` | 50/50 | Phase 72c MT=102 row-1 canary GREEN (10/10 nonzero cols) |
| `test_errorr_writer_mf_dispatch.jl` | 51/51 | (45 + 4 + 2) |
| `test_errorr_nc_expansion.jl` | 9/9 | Phase 73/74 NC-block canaries |
| `test_errorr_gendf_readback.jl` | 38/38 | MF31 (tape25) readback unchanged |
| Reference test **T22 (leapr light-water)** | BIT_IDENTICAL 4636/4636 | @rtol=1e-9 |

Total: **306/306** unit assertions PASS, zero regressions.

## Surprises

- The Phase 74 worklog framed this as a "rescon group-window edge" issue. It was not — `_resonance_group_window` was always correct. The real diagnosis required dumping per-(MT, MT2) per-row layouts of Julia's tape26 vs reference; the column-13/14 absence pattern across rows 13/14/15 in exactly the four pairs that share an URR contribution made the missing port the only consistent explanation. CLAUDE.md Rule 2 ("verify HANDOFF claims against current Fortran source and current Julia code — both drift") fired hard here: trusting the worklog's label would have wasted a session chasing a non-existent boundary bug.
- The MF=32 LRU=2 range header reports `LRF=1` in its CONT regardless of MF=2's LRF value. ENDF-6 §32.2 specifies a single URR cov format. So the wiring loop must match by [EL, EH] alone, not by LRF — matching by LRF fails for U-238 (MF=32 says LRF=1; MF=2 says LRF=2).
- Setting `sens` only inside `[iest, ieed]` and not multiplying by `cflx`/`abn` (matching the LRU=1 Julia convention, not the Fortran rpxunr convention which carries them through and then divides them back out in covout) gave bit-identical line counts on first compile. The two formulations are algebraically equivalent post-relative-cov-conversion, so we kept the simpler Julia LRU=1 convention.

## Remaining work for T15 to data-level BIT_IDENTICAL

Line counts now match exactly. Sub-line-count residue:

1. **MT=2/mt2=2 rows 10-12 sub-ULP FP precision** (Phase 72c follow-up #3): a `~5e-8` relative difference on the MT=102 C[1,1] canary remains the dominant data-level residue, plausibly Julia's tanh-stretched perturbation grid vs Fortran's adaptive rpendf eskip mesh + a 6-digit vs 5-digit `wt6b` constant in `weight_functions.jl` (3e-6 fractional thermal effect). Same FP class likely applies to a handful of URR-affected cells; sub-line-count so doesn't affect tape26 line totals. Defer until a wider FP-grind pass.
2. **LCOMP=0 / LCOMP=2 / LRF=7 (RML/SAMMY) MF=32 readers**: rpxlc0, rpxlc2, rpxsamm Fortran subroutines are not yet ported. These are deferred; T15 U-238 JENDL uses LCOMP=1 LRF=3 (resolved) + LRU=2 LRF=2 (URR, now handled) and is fully covered.

## Follow-up

- (P3) MT=2/mt2=2 rows 10-12 sub-ULP FP precision — defer to a future FP-grind.
- (P3) Multi-isotope MF=32: Julia URR port currently ignores `abn` (single-isotope materials have abn=1, so no effect on T15). If a future test hits a multi-isotope URR cov, factor abn back into the sandwich.
- (P2) Generalize `_urr_dofs_from_mf2` to handle URR2Data with non-uniform J-state counts vs MF=32 — the current 1:1 (L, J) walk assumes MF=2 and MF=32 list the same J-states in the same order. Verified for U-238 JENDL; document the invariant or add a defensive check.
