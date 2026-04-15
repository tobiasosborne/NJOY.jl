# Phase 12 — Small-fix batch (plotr + T20 errorr + T43 broadr)

## Date: 2026-04-15

## Summary

Three independent CRASH → DIFFS transitions in the style of Phase 10,
each with a dedicated RED → GREEN test. The fourth planned item
(T15/T17 BoundsError) turned out to have a different failure mode on
this machine — broadr completes instead of crashing, but takes 471 s on
U-238 JENDL-3.3 and trips the harness hard-timeout. Filed as a
performance bead rather than folded into this batch.

## Items

### 1. `plotr` dispatch stub

T06 was classified as CRASH in the Phase-10 sweep. The T10 worklog
listed it under "covr dispatch", but T06's input deck uses `plotr`
(which was `@info "skipped (visualization)"` in `pipeline.jl:223`) not
covr. Stubbed the same way as covr/purr/leapr.

- `src/orchestration/modules/plotr.jl` — **NEW** `plotr_module` touches
  the output tape (`nplt`).
- `src/orchestration/input_parser.jl` — `PlotrParams(nplt)` + `parse_plotr`.
- `src/orchestration/pipeline.jl` — replaced `plotr: skipped` with
  proper dispatch.
- `src/NJOY.jl` — include plotr.jl.

RED: `SystemError: opening file ".../tape31"` in viewr.
GREEN: T06 runs end-to-end, 0/1 tapes pass (38 s). Output still 0 vs
22,422 lines — stub produces no plot commands.

### 2. errorr "999 option" (dummy MF33 insertion) stub

T20 uses Fortran errorr's special "999" mode — card 1 is just the
integer 999, card 2 has `(nin, nout)`, and subsequent cards list MTs
for which to synthesize placeholder MF33 covariance sections. NJOY
rarely uses it but the reference test relies on it to produce the
reconr input tape.

`parse_errorr` couldn't handle the card layout (it treats card 1 as
six tape-unit integers). Added a dedicated detector in `pipeline.jl`'s
`:errorr` branch:

```julia
if !isempty(cards) && length(cards[1]) == 1 &&
   _parse_int_token(cards[1][1]) == 999 && length(cards) >= 2 &&
   length(cards[2]) >= 2
    nin  = abs(_parse_int_token(cards[2][1]))
    nout = abs(_parse_int_token(cards[2][2]))
    errorr_dummy_mf33_stub!(tapes, nin, nout)
else
    # normal errorr
end
```

`errorr_dummy_mf33_stub!` (added in `src/orchestration/modules/errorr.jl`)
copies `nin → nout`. Dummy MF33 synthesis is deferred. Also added
`params.nout <= 0` graceful-skip guard in the normal path.

RED: `tape unit 0 is invalid` (`resolve(tapes, 0)`).
GREEN: T20 runs full deck (7 module calls) in 168 s, 0/2 tapes pass.
tape23 reaches STRUCTURAL_FAIL (line-count diff) rather than CRASH.

### 3. broadr T=0 K pass-through

T43's input deck broadens at T=0 K. In the loop over temperatures,
`t_eff = temp - t_old` comes out to 0.0, `alpha = awr / (bk * 0)` is
`Inf`, and downstream `Int(alpha)` throws `InexactError: Int64(NaN)`.
Fortran broadr short-circuits this case — at T=0, broadening is a
delta function (no-op on the pointwise grid). Added matching guard in
`broadr_module`:

```julia
if t_eff <= 0.0
    @info "broadr: T=$(temp)K ΔT=$(t_eff) ≤ 0 — pass-through"
    push!(all_temp_results, (energies=cur_energies, total=cur_total,
                             partials=cur_xs, partial_mts=partial_mts,
                             temperature=temp))
    t_old = temp
    continue
end
```

RED: `InexactError: Int64(NaN)`.
GREEN: T43 runs end-to-end in 38 s, 0/1 tapes pass. The pass-through
emits the reconr grid verbatim as the "broadened at T=0" result.

### 4. (deferred) T15/T17 BoundsError

Sweep report classified these as `BoundsError: attempt to access
0-element Vector{Float64}`. On this machine (post-Phase 11), T15 no
longer crashes:

```
reconr: 246931 points written
broadr: MAT=9237 thnmax=10000.0 ntemp=1
broadr: T=300.0K alpha=9129.1
broadr: T=300.0K → 102946 points
broadr: wrote .../tape23
```

— broadr completes in 471 s. The harness hard-timeout (300 s per test)
fires during the *next* module, so T15 is now classified as
**TIMEOUT** rather than CRASH/BoundsError. That is a *move in the
taxonomy* but not a proper win.

Fortran broadr processes this same evaluation in seconds. The 471 s
is almost certainly a Julia-side performance regression — most likely
an O(N²) hot path in `broadn_grid` or `sigma1_at` for large point
counts (~250k points at 1e-3 reconstruction tolerance on U-238).

Tracked as bead **NJOY.jl-326** (broadr performance on JENDL U-238).
Not folded into this batch because fixing it is a performance grind,
not a crash-class fix.

## Sweep impact — estimated

| Test | Pre-P12 | Expected post-P12          |
|------|---------|----------------------------|
| T06  | CRASH   | DIFFS (0/1)                |
| T20  | CRASH   | DIFFS (0/2)                |
| T43  | CRASH   | DIFFS (0/1)                |
| T15  | CRASH   | TIMEOUT (broadr slow path) |
| T17  | CRASH   | TIMEOUT (assumed — not re-run; same evaluation)  |

So **3 tests move CRASH → DIFFS**, 2 tests move CRASH → TIMEOUT.
Full sweep not re-run.

Updated crash bucket estimate post-P12: **13 CRASHes** from the
pre-P11 16, minus T05 (fixed P11), minus T06/T20/T43 (this phase),
plus caveat that T15/T17 shift class rather than resolve. Phase
framework health unchanged: per-fix cost stays under an hour each.

## Tests

`test/validation/test_phase12_small_batch.jl` — 8/8 assertions pass in
2m 43s (T20 is the dominant contributor at ~168 s).

## Recommendations (priority order, post-P12)

### Immediate

1. **Plotr real output** — biggest remaining dispatch-shaped module.
   T06, and several other plotr/viewr tests (to survey) need real plot
   commands emitted.
2. **errorr `mt_list` parse + real dummy MF33 synthesis** — makes T20
   viable bit-identical. The MT list is already in the deck (cards
   3+); just needs a `synthesize_dummy_mf33(mt_list)` that writes
   Fortran-matching placeholder sections.
3. **T15/T17 broadr performance** (NJOY.jl-326) — profile
   `broadn_grid` / `sigma1_at` on the U-238 reconr output. Likely a
   searchsortedfirst or linear scan in an inner loop.

### Medium-term

4. **Real covr / leapr / purr output** (from Phase 10 / 11 leftovers).
5. **MF7/MT2 lattice reader** for Bragg.
6. **Grind the ~50 DIFFS** tests toward NUMERIC_PASS / BIT_IDENTICAL.

### Other crashes remaining (post-P12)

Per the pre-P11 taxonomy, still crashing (estimated):

- T09 MF7/MT4 not found (leapr stub too thin).
- T12/T18/T24/T27/T34/T47/T65 tape-unit plumbing (need inspection).
- T60 Fe-nat IRDFF-II MT=1 missing in broadr.

## Framework health

Still one hour per small fix. Filing a bead for the one item that
wasn't actually a crash this session (T15/T17) avoids conflating
"CRASH bucket shrinkage" with "TIMEOUT bucket growth".

## Commit

- `src/orchestration/modules/plotr.jl`: new.
- `src/orchestration/modules/errorr.jl`: `errorr_dummy_mf33_stub!` + guard.
- `src/orchestration/modules/broadr.jl`: T=0 pass-through.
- `src/orchestration/input_parser.jl`: `PlotrParams` / `parse_plotr`.
- `src/orchestration/pipeline.jl`: plotr dispatch, errorr 999 detection.
- `src/NJOY.jl`: include plotr.jl.
- `test/validation/test_phase12_small_batch.jl`: new RED→GREEN suite.
- `HANDOFF.md`, this worklog.
