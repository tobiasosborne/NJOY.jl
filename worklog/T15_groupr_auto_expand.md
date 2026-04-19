# T15/T16/T17 groupr — `3 /` auto-expand sentinel fix

## Date: 2026-04-19

## Summary

Fixed the Fortran auto-expand sentinel for bare groupr reaction cards
(`3 /` — mfd=3 with no mtd). T15 GENDF tape91 MF=3 MT count went from
**3 → 39** (Fortran reference: 40). Top HANDOFF Phase-44 next-step
blocker cleared; unblocks downstream errorr MT coverage for T15/T16/T17.

## Root cause

Fortran `groupr.f90:610-634` initialises `mtdp=-1000` before the
list-directed `read(nsysi,*) mfd,mtdp,strng`. A deck card like `3 /`
(tokens: only `mfd`) leaves `mtdp` at the sentinel. Line 622 branches
to label 382, sets `iauto=1`, and calls `nextr` (groupr.f90:1087-1123)
to walk MF=3 on the PENDF tape. `nextr` yields every MT passing the
inclusion filter (line 1115-1117):

```fortran
if (mtd.le.200) idone=1
if (mtd.ge.203.and.mtd.le.207) idone=1
if (mtd.gt.300) idone=1
```

Thermal (201-202) and engineering/derived (208-300) are skipped — they
must be named explicitly (which is why T15's deck has `3 251 'mubar'`
and `3 252 'xi'` *after* the `3 /` card).

### Julia pre-fix

`parse_groupr` (`input_parser.jl:439-446`) read missing mtd as `0` via
`_fint(..., default=0)` — indistinguishable from an explicit `3 0`
card. `groupr_module` fell into the `else` branch at line 57, ran
`haskey(mf3, 0)` which was false, and dropped the entry silently. T15
tape91 contained only the 5 explicit MTs (of which Julia emits 3: 452,
455, 456; MT=251/252 need MF=4/6 derivation, separate bead `NJOY.jl-cdy`).

## Fix

**Two places, matching the Fortran split between deck parse + runtime
reaction walk:**

1. `src/orchestration/input_parser.jl::parse_groupr` — distinguish
   "missing mtd" from "explicit mtd=0" by token count:
   ```julia
   mtd = length(cards[ci]) >= 2 ? _parse_int_token(cards[ci][2]) : -1000
   ```
   The `-1000` mirrors Fortran's exact sentinel value (groupr.f90:610).

2. `src/orchestration/modules/groupr.jl::_groupr_expand_auto` (new) —
   expands sentinel entries against `keys(mf3)` after the PENDF read,
   applying `_nextr_filter` and sorting ascending (Fortran reads the
   tape sequentially, MTs are always ENDF-ordered):
   ```julia
   _nextr_filter(mt) = mt > 0 && (mt <= 200 || (203 <= mt <= 207) || mt > 300)
   ```

Scoped to `mfd=3`. Other auto-expand mfd values (6/8/10/16-36) depend on
the `conver` subroutine's mf4/mf6/mf10 reaction lists, not yet ported.

## Verification

### RED→GREEN test

`test/validation/test_groupr_auto_expand.jl` — runs Julia groupr against
the cached `oracle_cache/test15/after_broadr.pendf` (skips the ~500s
U-238 broadening) and asserts the expected MT set.

| | Pre-fix | Post-fix |
|-|---------|----------|
| `@test` pass | 5 | 39 |
| `@test` fail | 37 | 0 |
| `@test_broken` | 0 | 3 |

The 3 `@test_broken` entries flag orthogonal pre-existing bugs
discovered during the grind:
- `NJOY.jl-cdy`: MT=251 (mubar) / MT=252 (xi) — Fortran derives from
  MF=4/6; Julia only does MF=3 lookup.
- `NJOY.jl-5oi`: MT=37 emitted with all-zero groups — threshold 17.82
  MeV above LANL-30 top 17.0 MeV. Fortran silently drops such MTs.

### Regression (T04, full pipeline)

Identical to HANDOFF Phase-44 baseline, zero regression:

| Tape | Status | Detail |
|------|--------|--------|
| tape23 | BIT_IDENTICAL | 82/82 lines match |
| tape24 | NUMERIC_PASS | 56/74 @ 1e-7, 74/74 @ 1e-5 |
| tape25 | DIFFS | 108/119 pass (known MF31 LB=2 cov residual) |

T01/T02 groupr decks use only explicit-mtd cards (length ≥ 2 tokens on
every reaction line) — untouched by the parser change.

## Files changed

- `src/orchestration/input_parser.jl` — sentinel detection in
  `parse_groupr` (4 lines: 2 docstring-comment + 1 logic).
- `src/orchestration/modules/groupr.jl` — added `_nextr_filter` and
  `_groupr_expand_auto` helpers (~35 LOC); wired `mt_list = _expand_...`
  before the dispatch loop.
- `test/validation/test_groupr_auto_expand.jl` — new RED→GREEN test.

## Trap (NEW)

**Trap (groupr `3 /` auto-expand)**: Fortran groupr init'd `mtdp=-1000`
before each `read(nsysi,*) mfd,mtdp,strng`. A card with only `mfd` leaves
the sentinel; goto 382 → `iauto=1` → `nextr` walks MF=3 on the PENDF and
yields MTs passing the filter `mt<=200 OR 203<=mt<=207 OR mt>300`.
Thermal (201-202) and derived (208-300) are excluded and must be named
explicitly in the deck. Any Julia deck parser reading missing mtd as
`default=0` collapses the sentinel into a dropped `(mfd, 0)` entry and
silently omits 30+ MTs from T15/T16/T17 tape91, starving downstream
errorr of group-averaged XS data.

## Follow-ups

- `NJOY.jl-cdy` (P2): compute MT=251/252 from MF=4/6 angular data.
- `NJOY.jl-5oi` (P2): skip MTs with all-zero group-averaged XS at write
  time (matches Fortran's silent-drop behavior for threshold-above-grid
  MTs like U-238 MT=37 under LANL-30).
- Downstream T15/T16/T17 errorr tape26 line-count parity — now that the
  GENDF has the 35 missing MTs, errorr's `colaps`-style readback (HANDOFF
  Phase-44 next-step #2) becomes the active blocker.
