# T15/T17 errorr output size — wrong group structure fix

## Date: 2026-04-18

## Summary

Julia's errorr module was producing **2.5 GB / 30 M-line tape26** for
T15/T17 (U-238 JENDL). Fortran reference is **5958 lines** (~500 KB).
4000× blowup.

After fix: tape26 = **7953 lines** (~640 KB), ~34 % over reference.
Remaining delta is a missing-MF3 + MF33 overcounting issue (covered
below, deferred to a follow-up).

## Root cause

`errorr_module` (`src/orchestration/modules/errorr.jl`) used
`_build_errorr_grid_from_endf` — the union of user_egn + all MFcov LIST-
record breakpoints — as the **output** group structure for every `ign`.

For T15/T17 which request `ign=3` (LANL-30, 30 groups), this union
contained **2305 points** (every MF33 breakpoint in U-238). Each
covariance matrix was then an expensive `2305 × 2305` row-major write →
30 M lines × 36 reactions → gigabytes.

Fortran `egngpn` (errorr.f90:9716) does this correctly:

```fortran
call gengpn(ign,ngn,egn)        ! library structure for ign >= 2
...
call uniong(nendf)              ! union of user + MFcov breakpoints
if (ign.eq.-1) then             ! REPLACE egn only when ign == -1
   allocate(egn(nunion+1))
   egn(1:nunion+1) = un(1:nunion+1)
endif
```

The union grid is always computed — but it only replaces the output
grid for `ign == -1` ("arbitrary structure, read in"). For `ign == 1`
(user-supplied with explicit ngn+1 boundaries) and `ign >= 2` (library
structures), the output grid stays at the user / library structure.
Internally, covcal may use the union grid for flux-weighted collapse,
but that is not the output.

HANDOFF's "T15/T17 INT=0 fix exposes downstream empty-vector" analysis
was incomplete — the downstream bug (groupr MT=455 LIST-skip, fixed in
`60f5bdc`) unmasked the crash, and the crash unmasked this
output-size bug.

## Fix

`src/orchestration/modules/errorr.jl`:

```julia
# Output group structure — matches Fortran egngpn (errorr.f90:9716).
# ign=-1: union of user_egn + MFcov breakpoints (replaces egn).
# ign=1:  user_egn exactly.
# ign>=2: library structure (LANL-30 for ign=3, etc.).
egn = if params.ign == -1
    _build_errorr_grid_from_endf(endf_path, params.mat, mfcov, params)
else
    _errorr_output_grid(params)
end
```

New helper `_errorr_output_grid` returns `user_egn` for `ign <= 1` and
the library structure (via existing `get_group_structure`) for
`ign >= 2`, falling back to LANL-30 with a warning on unknown ign.

## Verification

### T04 (ign=-1 first errorr, ign=1 second errorr) — no regression

```
tape23: NUMERIC_PASS 81/82 lines (passes at 1e-07, 1e-05, 1e-03)
tape24: NUMERIC_PASS 56/74 lines (passes at 1e-05, 1e-03)
tape25: DIFFS 107/119 lines (residual MF31 LB=2 covcal, from T03_phase7)
```

Identical to pre-fix baseline on all tapes.

### T15 (ign=3 errorr ×3, U-238 JENDL)

Pre-fix: `CRASH` (BoundsError until MT=455 fix), then 30M lines tape26.
Post-fix: all 10 tapes classified `STRUCTURAL_FAIL` (line counts differ
but file is sane). Key numbers:

| Tape | Pre-fix lines | Post-fix lines | Ref lines |
|------|---------------|----------------|-----------|
| tape25 (MF31 cov) | — | 950 | 392 |
| tape26 (MF33 cov) | 30 253 210 | **7953** | 5958 |
| tape27 (MF34 cov) | — | 112 | 15 |

### T17 (U-238 + U-235 + Pu-239 JENDL, larger deck)

Pre-fix: `CRASH`. Post-fix: completes end-to-end, no crash.

## Deferred follow-up: remaining 2k-line overshoot

The ~2k line excess on tape26 has two contributors:

1. **MF3 sections missing** (-216 lines). Julia's MF3 write loop at
   `_write_errorr_tape` line 714 is gated on
   `haskey(group_xs, mt)`. For T15/T17, `npend == 0` so the existing
   `_errorr_group_average` path (gated on `params.npend > 0`) never
   populates `group_xs`. Fortran in this case uses `colaps` to read
   group-averaged cross sections from the GENDF input tape (`ngout`).

2. **MF33 over-expansion** (+2.2k lines). Julia writes a full `ngn × ngn`
   row block for each MT pair, even when the Fortran collapses to a
   sparse block format (LB=5 triangle). Visible in the per-MT counts:
   MT=1 self: 271 (REF) vs 287 (JUL) — close. MT=2 self: 1400 vs 103
   — Julia has *fewer* here because many off-diagonal cov entries are
   being dropped. The structure is different enough that per-MT
   comparison is not directly meaningful; needs a closer look at
   `_write_errorr_tape` MFcov loop and Fortran `covout`.

Both items require wider work (read Fortran `colaps` / `covout`, wire up
GENDF reading for `group_xs` fallback, restructure the covariance
writer). Not tackled in this patch to keep the fix surface minimal.

There is also a **Julia groupr `3 /` auto-expand** gap — T15's groupr
deck uses the Fortran sentinel `3 /` (mfd=3 with no mtd) which Fortran
expands to "process all MF=3 MTs automatically" (groupr.f90:622 →
label 382, iauto=1). Julia's parser reads this as (mfd=3, mtd=0) and
produces a GENDF with only the 6 explicitly-listed MTs, so the
downstream errorr wouldn't find the 36 covariance MTs even with the
colaps read path fixed. This is the top blocker for getting T15 to
line-count parity.

## Files changed

- `src/orchestration/modules/errorr.jl` — dispatch output grid
  construction on `params.ign`; new `_errorr_output_grid` helper.

## Traps (NEW)

**Trap (errorr output grid)**: `ign` semantics for errorr's **output**
group structure are:

- `ign == -1`: union of user_egn + MFcov breakpoints (Fortran
  `errorr.f90:9787`).
- `ign == 1`: user_egn exactly (nothing more).
- `ign >= 2`: library structure (LANL-30 etc., via `gengpn`).

The MFcov breakpoint union is NEVER the output grid for `ign >= 1`.
Using it as output for `ign == 3` inflates LANL-30's 30-group write
to a 2305-group write, blowing the file size up ~4000×.

**Trap (groupr `3 /` auto-expand)** — unfixed: a groupr input-deck
card of the form `<mfd> /` with no mtd is a Fortran sentinel meaning
"auto-process all MTs for that MF". Handled at `groupr.f90:622` with
`iauto=1` / `nextr`. Julia's `parse_groupr` reads the missing mtd as
default 0 and produces a single `(mfd=3, mtd=0)` entry, missing the
auto-expansion. This matters for T15/T17/T16 decks that use `3 /` to
request all cross sections for errorr covariance grouping.
