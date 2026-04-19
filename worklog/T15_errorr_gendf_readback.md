# T15/T16/T17 errorr — GENDF MF3 readback (colaps same-ign path)

## Date: 2026-04-19

## Summary

T15 tape26 MF=3 sections: **0 → 36 MTs** (auto-expanded `3 /` from the
previous phase now visible in the covariance output). tape25 MF=3: **0
→ 3** nubar MTs. tape26 line count: 7953 → 8205 (closer to reference
5958; remaining ~2.2k delta is the HANDOFF-known MF33 over-expansion,
separate bug). Second immediate-next-step of HANDOFF Phase 45 done.

## Root cause

`errorr_module` (`src/orchestration/modules/errorr.jl`:90-94) gated
`group_xs` population on `params.npend > 0`:

```julia
group_xs = Dict{Int, Vector{Float64}}()
if params.npend > 0 && params.iread == 0
    pendf_path = resolve(tapes, params.npend)
    group_xs = _errorr_group_average(pendf_path, params.mat, egn, params.iwt)
end
```

T15/T16/T17 errorr decks run with `npend=0, ngout=91` — they rely on
the groupr GENDF for group-averaged cross sections, since the deck
skips a per-temperature PENDF pass. Julia left `group_xs` empty → the
MF3 writer at `_write_errorr_tape` line 716 skipped every MT via
`haskey(group_xs, mt) || continue`. For T15, all 36 MF33-covariance
MTs were missing their MF3 cross-section output.

Fortran handles this branch via `colaps` (errorr.f90:9097-9532): walk
MF=3 on ngout, emit one LIST record per output group with
`(flux, [nubar,] sigma)`.

## Fix

New helper `_errorr_read_gendf_xs` (`src/orchestration/modules/
errorr.jl`) that reads per-group sigma from the input GENDF. Wired
into the main dispatch via a second arm of the `group_xs` guard:

```julia
elseif params.npend == 0 && abs(params.ngout) > 0 && params.iread == 0
    gendf_path = resolve(tapes, abs(params.ngout))
    group_xs = _errorr_read_gendf_xs(gendf_path, params.mat, egn)
end
```

**Scope**: same-ign only. For T15/T16/T17, groupr uses `ign=3`
(LANL-30, 30 groups) and errorr uses `ign=3` → group structures
match and direct readback is correct. Cross-ign flux-weighted
collapse (the full Fortran `colaps` algorithm at errorr.f90:9255-
9283) is NOT implemented — the helper logs a warning and returns
empty when group counts disagree. Filed as `NJOY.jl-<next>` for the
cross-ign case once a test exercises it.

### Reader layout

GENDF MF=3/MT sections (matches Julia groupr output at
`groupr.jl:317`, mirroring Fortran `grpav`):

- HEAD:   `ZA, 0, NL=1, NZ=1, 0, NGN` (seq=1 sentinel)
- Per g:  CONT `(temp, 0, NG2, NL, NW, IG=g)` + body line with NW
          floats; **position 2 of the body is always the relevant
          value** — sigma for standard MTs, nubar for MT=452/455/456.

Reader parses line-by-line using column-fixed slicing (matches
existing `_read_gendf_nubar` / `_read_gendf_group_bounds` style).

## Verification

### RED → GREEN test

`test/validation/test_errorr_gendf_readback.jl`:

| | Pre-fix | Post-fix |
|-|---------|----------|
| `@test` pass | 2 | 38 |
| `@test` fail | 36 | 0 |

Exercises both errorr calls in T15 (mfcov=33 XS cov → tape26;
mfcov=31 nubar cov → tape25). Skips the slow broadening by driving
`errorr_module` directly against the cached broadr PENDF + a
freshly-generated Julia GENDF.

### Tape line-count summary (T15)

| Tape | Pre | Post | Reference | Notes |
|------|-----|------|-----------|-------|
| tape25 (MF31 nubar cov) | 392 | 413 | 392 | +21 lines (MF3 section overhead) |
| tape26 (MF33 XS cov) | 7953 | 8205 | 5958 | +252 MF3 lines; remaining +2247 is MF33 over-expansion (HANDOFF Phase-44 known issue) |
| tape27 (MF34 angular cov) | 112 | unchanged | 15 | not in scope |

### Regression (T04, full pipeline)

Zero regression vs Phase-45 baseline:
- tape23 BIT_IDENTICAL 82/82
- tape24 NUMERIC_PASS 56/74 @ 1e-7 (groupr, unrelated)
- tape25 DIFFS 108/119

T04's first errorr has `npend=-22 > 0` → existing PENDF path.
T04's second errorr has `iread=1` (grid inheritance) → different
branch; my guard `params.iread == 0` keeps it untouched.

## Files changed

- `src/orchestration/modules/errorr.jl` — wired `elseif` for GENDF
  readback; added `_errorr_read_gendf_xs` (~55 LOC including doc).
- `test/validation/test_errorr_gendf_readback.jl` — new RED→GREEN
  test.

## Trap (NEW)

**Trap (errorr MF3 source by tape)**: errorr's group-averaged MF3
cross sections come from **either** the PENDF (`npend>0`, via
`grpav`-style panel averaging) **or** the input GENDF (`npend=0`,
`ngout>0`, via `colaps` re-collapse or direct readback when group
structures match). Julia must dispatch on `(npend, ngout)` just like
Fortran does at errorr.f90:1257-1262 (top of main) and call either
`_errorr_group_average` or `_errorr_read_gendf_xs`. Gating only on
`npend>0` silently produces a covariance tape with empty MF3
sections — MF33 blocks write fine but the companion cross-section
data the covariance *refers to* is missing. 216+ missing output
lines for T15/T17 tape26.

## Follow-ups

- **Cross-ign colaps collapse**: port the flux-weighted re-mapping at
  errorr.f90:9255-9283 for cases where GENDF ign ≠ errorr ign. Not
  blocking T15/T16/T17 (same ign). File a bead if a test exercises it.
- **T15 tape26 MF33 over-expansion** (+2247 lines): HANDOFF Phase-44
  known issue. Julia writes full NGN×NGN row blocks vs Fortran's LB=5
  sparse triangle. Next T15 blocker.
