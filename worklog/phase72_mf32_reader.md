# Phase 72 — MF=32 reader + apply_rescon! foundation

**Date:** 2026-05-08
**Test:** T15 (U-238 JENDL, MAT=9237) — RED canary `MT=102 C[1,1] ≈ 2.658914e-4`
**Goal:** Land the foundation for the multi-session rescon port — MF=32
reader (LRU=1/LRF=1,2,3/LCOMP=1/NRO=0/NLRS=0), `apply_rescon!` skeleton,
and tests against the actual U-238 JENDL ENDF data.

## Outcome

**Foundation landed; canary stays RED (by design).** `read_mf32`
parses the U-238 MF=32 tape into typed Julia structs (10 resolved
ranges, 317 resonances, 951 uncertain parameters); URR (LRU=2) range
is skipped cleanly with a clear info log. `apply_rescon!` is wired
into errorr_module replacing the Phase-71 `@warn` — it parses + logs
+ returns without modifying `cov_matrices`. The Phase-71 RED canary
in `test_errorr_covcal_lb5.jl` stays `@test_broken` (rows 1..14 of
MT=102 self-cov are still zero in Julia output) until the
next-session sensitivity-builder port lands.

This is intentional senior-quality scoping: the sensitivity builder
requires faithfully porting `rpendf` (adaptive RM grid with
eskip1/2/3/4 step factors), `rpxgrp` (lethargy-weighted trapezoidal
integration with `egtwtf` flux weights), the MF=2 cross-resonance
match, and the 78×2-call ±0.0001/±0.01 perturbation loop. Splitting
that off avoids landing partial/buggy sensitivity code that could
silently corrupt MF=33 output.

## What landed

### 1. `src/processing/mf32_reader.jl` (NEW, 213 LOC)

Reader for MF=32 LRU=1 covariance, LCOMP=1 only:
- **Top-level**: `MF32Data` (za, awr, isotopes).
- **Per-isotope**: `MF32IsotopeData` (zai, abn, lfw, resolved_ranges).
- **Per-range**: `MF32ResolvedRange` (elr, ehr, lrf, naps, spi, ap,
  nls, isr, dap, subsections).
- **Per-subsection**: `MF32ResolvedSubsection` (awri, mpar, nrb,
  params [nrb, 6], cov [npar, npar]).
- **Helpers**: `mf32_resonance_count`, `mf32_param_count` for
  sensitivity-builder array sizing.

Mirrors Fortran `resprx` (errorr.f90:3011-3250) range walk + `rpxlc12`
(errorr.f90:4252-4274) subsection LIST parse. URR (LRU=2) is skipped
cleanly: per-format CONT + optional ISR=1 CONT + NLS L-state LISTs +
one covariance LIST, mirroring the rpxunr record sequence
(errorr.f90:4824-4837). Emits `@info` log when skipping; orthogonal to
the resolved-range path.

Errors loudly on unsupported branches (LCOMP=0/2, LRF=4/7, NRO≠0,
NLRS≠0, LRU=2 except via skip path, ISR=1 with LRF=3) per Rule 6.

### 2. `src/processing/rescon.jl` (NEW, 86 LOC)

- **`RESCON_PAIRS`**: const tuple of the seven (mt, mt2) pairs that
  receive an RP-cov contribution, in Fortran rescon `itp` order
  (errorr.f90:8531-8539).
- **`rescon_pair_index(mt, mt2)`**: maps a pair to itp 1..7 or 0.
  Canonical order mt ≤ mt2.
- **`rescon_supports_pair(mt, mt2)`**: bool predicate.
- **`apply_rescon!(cov_matrices, endf_path, mat, egn, group_xs)`**:
  parses MF=32, validates, logs the per-MAT count summary, returns
  without mutating cov. Embedded comment block specifies the
  next-session sensitivity-builder contract verbatim (channel mapping
  for sens[1..4], output array dispatch for itp=1..7).

### 3. Wired into `src/orchestration/modules/errorr.jl`

Replaced the Phase-71 `@warn` block (lines 222-245) with a direct
`apply_rescon!(cov_matrices, endf_path, params.mat, egn, group_xs)`
call. Behaviour unchanged for the seven (mt, mt2) pairs — the
function logs and returns without mutation. Comment block updated
to reference Phase 72 status + next-session sensitivity port.

### 4. `test/validation/test_mf32_reader.jl` (NEW, 134 LOC, 57 assertions)

Five testsets, all GREEN against the real `J33U238` ENDF tape:
- **top-level header**: ZA=92238, AWR=236.006, NIS=1, NER=11
  (10 resolved + 1 URR), LFW=0, ABN=1.0.
- **range 1 (1e-5..1000 eV, RM/LCOMP=1)**: EL/EH, LRF=3, SPI=0,
  AP=0.942848, LCOMP=1, NLS=0, ISR=0; subsection AWRI=236.006,
  MPAR=3, NRB=26; first resonance (ER, AJ, GN, GG, GFA, GFB) =
  (6.674, 0.5, 1.493e-3, 2.300e-2, 0.0, 9.99e-9).
- **covariance symmetry**: cov is symmetric, diagonal ≥ 0 across
  all 10 resolved ranges.
- **aggregate counts**: ≥26 resonances, NPAR = 3·NRB (MPAR
  consistency).
- **range 2 (1000..2000 eV)**: NRB=28, first resonance ER=1054.65 —
  proves multi-range parsing.

### 5. Exports added to `src/NJOY.jl`

```julia
export MF32Data, MF32IsotopeData, MF32ResolvedRange, MF32ResolvedSubsection
export read_mf32, mf32_resonance_count, mf32_param_count
export apply_rescon!, rescon_pair_index, rescon_supports_pair, RESCON_PAIRS
```

## Regression verification

| Test | Pre-Phase-72 | Post-Phase-72 |
|------|--------------|---------------|
| `test_errorr_covcal_lb5.jl` Bug A NK suite | 5 PASS | 5 PASS |
| `test_errorr_covcal_lb5.jl` Phase 71 RED canary | 18 PASS + 2 BROKEN | 18 PASS + 2 BROKEN |
| `test_mf32_reader.jl` (NEW) | n/a | 57 PASS |
| T15 errorr tape26 line count | 5964 | 5964 |
| MT=102 row-1 non-zero col count | jul=0 ref=10 | jul=0 ref=10 |

T15 errorr behavioural output unchanged — `apply_rescon!` parses +
logs + returns. Bug A geometry preserved. Phase 71 canary stays
@test_broken (next session flips it).

## What's next (Phase 73)

Sensitivity builder port — the bulk of the remaining work. Order:

1. **Port `rpendf`** (errorr.f90:5015-5089, ~75 LOC): adaptive
   pointwise σ(E) generator over a single MF=2 resonance. Uses
   eskip1/2/3/4 step-factor table (proximity-weighted) for the
   adaptive grid, reuses existing `cross_section_rm` (or analogous
   LRF=1/2 evaluator) at each E. Output: `Vector{NTuple{5,Float64}}`
   of (sig_total, sig_elastic, sig_fission, sig_capture, energy).
   Two `rpendf` calls per perturbation pair.

2. **Port `rpxgrp`** (errorr.f90:5227-5367, ~140 LOC): lethargy-
   weighted trapezoidal integration of pointwise σ over the output
   group grid. Uses `egtwtf` (per-output-group flux weight). Crosses
   group boundaries by linear interpolation. Returns
   `Matrix{Float64}(4, ngn)`.

3. **Port the rpxlc12 perturbation loop** (errorr.f90:4399-4523,
   ~125 LOC, the inner core): for each subsection, for each
   `(loopm, loopn)` ∈ NRB × MPAR: identify the (resonance,
   parameter) being perturbed, ±0.01% for ER (loopn=1) or ±1% for
   widths (loopn≥2), call rpendf+rpxgrp twice, central-difference
   into `sens[channel, p, ig]`.

4. **Port the sandwich fill** (errorr.f90:4541-4593, ~55 LOC):
   triangular `cgg`/`cee`/`cff`/`ctt` for self-cov pairs; full-matrix
   `cef`/`ceg`/`cfg` for cross pairs. Cross-term doubling when
   `i ≠ j`.

5. **Wire into `apply_rescon!`** (rescon.jl): per range, per
   subsection, dispatch to `_build_sensitivities_lcomp1` then
   accumulate into `cov_matrices[(mt, mt2)]` for the seven pairs.

6. **GREEN the canary**: `MT=102 C[1, 1] ≈ 2.658914e-4` within 1e-7;
   `MT=102 C[1, 9] < 0` (negative-sandwich signature).

Estimated cost: 1–2 sessions for rpendf+rpxgrp+perturbation +
sandwich; then a focused FP-grind session if values are close but
not within 1e-7. The Fortran perturbation factors (×0.9999/×1.0001
for ER, ×0.99/×1.01 for widths) and central-difference
denominator (`gwidth*2`) are bit-faithful — no FP-order ambiguity.

## Open follow-ups (smaller scope, not blocking)

- **URR cov port** (rpxunr, errorr.f90:4785-5013): fills uee/uff/ugg/
  utt/uef/ueg/ufg. Different perturbation strategy (×1.01 only,
  forward-difference). For T15 the URR range [10 keV, 150 keV] is
  orthogonal to the MT=102 row 1..14 canary which lives entirely in
  resolved [1e-5, 1000] eV.
- **LCOMP=0 reader** (rpxlc0, errorr.f90:3734): legacy format; not
  used by U-238 JENDL but needed for some older evaluations.
- **LCOMP=2 reader** (rpxlc2, errorr.f90:4634): compact format with
  INTG correlation rows + NDIGIT rounding.
- **LRF=7 reader** (rpxsamm, errorr.f90:3252): SAMMY/RML, fills the
  4-D `crr` array directly. Required for some modern evaluations
  (ENDF/B-VIII U-235, etc.).

## Files

NEW:
- `src/processing/mf32_reader.jl`
- `src/processing/rescon.jl`
- `test/validation/test_mf32_reader.jl`
- `worklog/phase72_mf32_reader.md` (this file)

MODIFIED:
- `src/NJOY.jl` — include order + exports.
- `src/orchestration/modules/errorr.jl` — `@warn` → `apply_rescon!`.

## Fortran source citations

- `errorr.f90:3011-3250` — `resprx` (MF=32 dispatch).
- `errorr.f90:3088-3231` — per-range CONT + DAP handling.
- `errorr.f90:4108-4632` — `rpxlc12` (LCOMP=1/2 reader, the bulk).
- `errorr.f90:4252-4274` — LIST parse + MPAR/NRB extraction.
- `errorr.f90:4528-4535` — upper-triangular packed cov read.
- `errorr.f90:4824-4837` — rpxunr URR record sequence (drives the
  Phase-72 URR-skip path).
- `errorr.f90:8513-8819` — `rescon` itself (the in-place add into
  `cova` — the part that's structurally trivial once sensitivities
  exist).
