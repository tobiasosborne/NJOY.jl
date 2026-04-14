# Phase 10 — Batch dispatch + scope/reader fixes

## Date: 2026-04-14

## Summary

Landed five independent fixes from the priority-ordered list in
`worklog/T09_acer_dispatch.md` §Recommendations, then re-swept:

1. **Dispatch `:purr`** in the pipeline (stub: copy npendf_in → npendf_out).
2. **Dispatch `:leapr`** in the pipeline (stub: write empty output tape).
3. **Fix `UndefVarError: h` in heatr** — uninitialised local when `mt`
   doesn't match any branch in the kerma elseif chain.
4. **Gate Bragg lookup by `icoh`** in `_thermr_sab!` / `_recompute_thermr_mf6!`
   and graceful degrade when lattice parameters are unknown.
5. **Accept `INT=0` in MF TAB1 reader** — coerce to INT=2 (linear) in the
   `InterpolationTable` constructor to match Fortran NJOY silent acceptance.

All five were implemented **RED → GREEN**: the test framework's
`run_reference_test(N)` served as the failing assertion (pre-fix it
surfaced the exact crash string from the baseline sweep), and the fixes
were verified test-by-test before running the full 84-test sweep.

## Sweep result — post-Phase-10

`reports/REFERENCE_SWEEP.md`, 9784.8s (163 min):

| Status | Pre-Phase-10 | Post-Phase-10 | Δ |
|--------|--------------|---------------|----|
| BIT_IDENTICAL    | 1  | 1  | = |
| NUMERIC_PASS     | 1  | 1  | = |
| DIFFS            | 12 | 48 | **+36** |
| MISSING_TAPE     | 14 | 17 | +3 |
| NO_REFERENCE     | 1  | 1  | = |
| **CRASH**        | **55** | **16** | **−39** |

**39 tests transitioned CRASH → DIFFS.** Pipeline now runs end-to-end on
68 of 84 tests.

## Files changed

| File | Change |
|------|--------|
| `src/NJOY.jl` | Include `orchestration/modules/{purr,leapr}.jl` |
| `src/orchestration/pipeline.jl` | `:purr` + `:leapr` dispatch branches; Bragg-stub `try/catch` wrapper |
| `src/orchestration/modules/purr.jl` | **NEW** stub `purr_module(tapes, params)` |
| `src/orchestration/modules/leapr.jl` | **NEW** stub `leapr_module(tapes, params)` |
| `src/orchestration/modules/thermr.jl` | Gate Bragg by icoh > 0; try/catch missing lattice with graceful inelastic-only fallback |
| `src/orchestration/input_parser.jl` | Added `LeaprParams` struct + `parse_leapr(mc)` |
| `src/endf/types.jl` | `InterpolationTable` constructor coerces INT=0 → INT=2 |
| `src/processing/heatr.jl` | Initialise `h = 0.0` before kerma elseif chain (line ~1180) |
| `test/validation/test_phase10_dispatches.jl` | **NEW** RED→GREEN assertion suite |
| `reports/REFERENCE_SWEEP.md` | Regenerated |

## Per-item detail

### 1. `:purr` dispatch — STUB

- `parse_purr` already existed (from an earlier port) but no orchestration
  wrapper. New file `src/orchestration/modules/purr.jl` implements a
  stub: `cp(npendf_in, npendf_out; force=true)` so downstream acer/viewr
  see a readable tape.
- Full impl (generating MT152 ENDF sections from
  `generate_ptable`) is deferred — tracked in the "Missing for
  bit-identical" section below.
- Verified T35–T42, T28, T63 transitioned CRASH → DIFFS.

### 2. `:leapr` dispatch — STUB

- No `parse_leapr` existed; added `LeaprParams(nout::Int)` +
  `parse_leapr(mc)` (reads only card 1's first integer). The leapr deck
  has many additional cards (oscillators, DOS, descriptions) — the stub
  doesn't parse them.
- `leapr_module` calls `touch(tape{nout})` — empty output.
- T22, T23, T80 transitioned CRASH → DIFFS (structural, not
  numeric).
- **Exposed new downstream failure**: T09 now crashes with `MF7/MT4 not
  found for MAT=101 in tape24` because thermr tries to read real S(α,β)
  data from the empty leapr output. Was already broken upstream (no
  leapr dispatch); now fails earlier with a clearer, more actionable
  message.

### 3. heatr `h` scope bug — ONE-LINE FIX

At `src/processing/heatr.jl:~1180`, the per-energy kerma loop:

```julia
for ie in 1:ne
    sigma = ... ; sigma <= 0.0 && continue
    if mt == 2
        h = elastic_heating_aniso(...) * sigma
        elastic[ie] += h
    elseif mt == 102 ... h = ...
    elseif mt == 18 || ... ; h = fission_heating(...) * sigma
    elseif 51 <= mt <= 90 ; h = (E + q0 - ebar) * sigma
    elseif mt == 91      ; h = ...
    elseif mt == 16 || mt == 17 || mt == 37 ; h = ...
    elseif mt > 100      ; h = (E + Q) * sigma
    end
    total[ie] += h   # ← UndefVarError if no branch matched
end
```

MTs that fall through **all** branches (MT=22–29, 41–45 `(n,γ+charged
particle)`, MT=152+, plus delayed-fission markers) had `h` uninitialised
when the `total[ie] += h` ran. Fix: initialise `h = 0.0` at the top of
the inner loop, matching the implicit contribution-zero semantics
(these MTs have no heating model yet).

Resolved T08, T13, T21, T26, T49, T79.

### 4. Bragg lookup gating + graceful missing-entry handling

Two bugs in one:

**(a) `icoh == 0` was not gated.** For T32, T68 (thermr with no coherent
elastic), `_thermr_sab!` unconditionally called `lookup_bragg_params`,
erroring on MAT=58 / MAT=1 even though the user explicitly turned
coherent elastic off. Fixed by wrapping the Bragg construction in
`if icoh > 0 ... end`.

**(b) `icoh > 0` with unknown MAT hard-crashed.** For T25/T67/T69/T70/T74
the user *does* want coherent elastic, but `BRAGG_LATTICE_PARAMS` has no
entries for MAT 1/7/15/53/58. **These materials use ENDF-6 format**
where NJOY reads the lattice info directly from MF7/MT2 (see
`thermr.f90:406-415`) — the lattice *table lookup* is a legacy ENDF-4/5
codepath. The proper fix is to port `coh()`'s MF7/MT2 reader; tracked as
a follow-up.

For now: wrap the lookup in `try/catch`, emit a `@warn` naming the MAT,
skip the MT=mtref+1 (coherent elastic) section, continue with inelastic
only. The five tests now produce a partial PENDF (MT=mtref only) and
run to completion. This matches the no-approximations rule — we emit
no fake lattice data.

**Decision deliberately taken**: I did **not** add hardcoded
`BRAGG_LATTICE_PARAMS` entries for MAT 1/7/15/53/58 even though the
worklog item said "add Bragg lattice entries". Those entries would
either be wrong (matde=1 = H-in-H2O is liquid, has no lattice) or
approximate (ZrH for matde=7/58 has no canonical NJOY built-in, would
need to be fabricated). The graceful-degrade path is the
no-approximations-compliant option and leaves the correct fix
(MF7/MT2 reader port) clearly scoped for a future session.

### 5. INT=0 in MF TAB1 reader

`src/endf/types.jl`:

```julia
function InterpolationTable(nbt, law)
    normed = [li == 0 ? Int32(2) : Int32(li) for li in law]
    new(Int32.(nbt), InterpolationLaw.(normed))
end
```

JENDL-3.3 U-238 evaluations (T15, T17 — MAT=9237) emit INT=0 in some
MF3/MF4 TAB1 records. This is **technically illegal per ENDF-6
standard**, but Fortran NJOY silently accepts it as linear-linear
(INT=2). Julia's `@enum InterpolationLaw` starting at 1 rejected the 0
with `ArgumentError`. Fix: coerce 0 → 2 at the struct constructor.

Verified T15/T17 got past the reader: T15 now crashes much later with a
`BoundsError: attempt to access 0-element Vector{Float64} at index
[1]` during broadr. This is a *different* JENDL-specific bug that the
INT=0 fix revealed. Not a regression — the pre-fix crash was masking
the second issue.

## Missing for bit-identical output

Tracked here so the next session can pick up directly.

### Real purr output (MT152)
`generate_ptable` already runs the MC ladder algorithm. Needs:
- ENDF MF2/MT152 section writer (table format: LIST of probability bins
  per energy).
- Integration with `read_pendf` → `copy_with_modifications` to splice
  MT152 into the output PENDF.
- Only iopt=1 / lssf=0 path matters for the 13 purr-gated tests.

### Real leapr output (MF7)
This is the biggest remaining gap:
- Parse the full leapr deck: nout, title, ntempr, mat/za, awr/spr/npr/iel/ncold,
  energy grid, oscillators, T-effective, Bragg (card 18), descriptions.
- Wire to `generate_sab` (already ported in `src/processing/leapr.jl`).
- Emit MF7/MT2 (elastic) + MF7/MT4 (inelastic) ENDF sections.
- Unblocks T22/T23/T33/T80 (leapr-only) and T09, T67–70, T74 (thermr
  chains that consume leapr output).

### MF7/MT2 lattice reader (for Bragg in thermr)
Currently the Bragg lookup uses a 7-entry hardcoded
`BRAGG_LATTICE_PARAMS` table keyed by matde. Modern ENDF-6 evaluations
(including most of the failing tests) put the lattice reciprocal-cell
constants directly in MF7/MT2. The Fortran routine is
`thermr.f90::coh()` / `rdelas()`; port would be ~150 lines.

### T15/T17 JENDL-3.3 BoundsError
Separate failure mode surfaced by the INT=0 fix. Need to trace the
`0-element Vector{Float64}` — likely an empty interpolation table or
MF3 section getting indexed. The Fortran reference runs cleanly, so
this is purely a Julia bug (probably in the MF3 reader or broadr
thinning).

### Other crashes (post-Phase-10)
- **T05, T06**: covr dispatch (covr is still not wired in pipeline).
- **T12**: broadr tape23 — broadr writes to wrong unit or in wrong mode.
- **T16**: covr/errorr-related tape36.
- **T18**: groupr MF6 transfer matrix missing for U-233.
- **T20**: moder card with tape unit 0 (input parsing edge case).
- **T24, T27, T34**: unresr or groupr output tapes — need dispatch or
  tape-unit plumbing.
- **T43**: broadr at T=0 → NaN → Int conversion.
- **T47**: chained unresr output missing.
- **T60**: Fe-nat IRDFF-II (MF10-only, no MF3) — broadr MT=1 assumption.
- **T65**: unresr tape51 missing.

## Recommendations (priority order)

### Immediate (small, high-yield)
1. **Dispatch `covr`** — unblocks T05, T06, T16. `covr_module` already
   exists in `src/processing/covr.jl`.
2. **Investigate T15/T17 BoundsError** — JENDL U-238. Now that the INT=0
   mask is off, run T15 with verbose, locate the empty-vector access,
   fix. Most likely a one-liner.
3. **T20 moder tape-unit-0** — input-parser edge case; the `moder` card
   uses `0` as a "no tape" sentinel. Handle that in `moder_module`.

### Medium-term
4. **Real purr MT152 writer** — deferred but straightforward, unblocks
   finer DIFFS reduction for T28, T35–T42, T63.
5. **Full `leapr` port** — biggest single remaining module. Will
   probably need its own multi-cycle task.
6. **MF7/MT2 lattice reader** — unblocks real Bragg output for T25,
   T67–70, T74.

### Long-running
7. Grind the 48 DIFFS cases to bit-identical. Most are small
   per-tape deltas once the pipeline is complete (e.g. T46: only MT=1
   has ±1 diffs; T34, T20: known R-matrix FP class).

## Framework health

The Phase 8 test framework is continuing to pay off. Adding five
independent dispatches + fixes and measuring the exact impact in one
re-sweep took ~3 hours total (implementation + verification +
163-minute sweep). Without the framework this would have been at least
a full day of per-test manual orchestration.

## Commit

- `src/orchestration/modules/purr.jl`, `src/orchestration/modules/leapr.jl`: new.
- `src/NJOY.jl`, `src/orchestration/pipeline.jl`: include + dispatch.
- `src/orchestration/input_parser.jl`: LeaprParams/parse_leapr.
- `src/endf/types.jl`: INT=0 → INT=2 coercion.
- `src/orchestration/modules/thermr.jl`: Bragg gating + graceful degrade.
- `src/processing/heatr.jl`: h=0.0 initialization.
- `test/validation/test_phase10_dispatches.jl`: new RED→GREEN suite.
- `reports/REFERENCE_SWEEP.md`: regenerated.
- `HANDOFF.md`, this worklog.
