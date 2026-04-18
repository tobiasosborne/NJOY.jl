# T15/T17 JENDL U-238 — MT=455 LIST-skip crash fix

## Date: 2026-04-18

## Summary

Last 2 CRASHes in the 84-test sweep (T15, T17 — both JENDL U-238, MAT=9237)
drained to DIFFS. `BoundsError: attempt to access 0-element Vector{Float64}
at index [1]`.

Also ran a fresh full sweep for current ground truth — **93.9 min** on this
machine (prev 163 min, ~42% faster). Pre-fix state:

| Status | Post-P10 (2026-04-14) | Pre-fix (2026-04-18) | Δ |
|--------|-----------------------|----------------------|---|
| BIT_IDENTICAL | 1 | 1 | = |
| NUMERIC_PASS  | 1 | 1 | = |
| DIFFS         | 48 | 62 | +14 |
| MISSING_TAPE  | 17 | 17 | = |
| NO_REFERENCE  | 1 | 1 | = |
| CRASH         | 16 | **2** | **-14** |

The -14 CRASH reflects T11-T14 phases (covr / plotr / broadr T=0 / errorr
999 / moder extract / gaspr / heatr plot-stub / broadr MT=1 fallback /
thermr SAB fallback) cascading through the sweep for the first time.

## Root cause

HANDOFF §Sweep and §3 described T15/T17 as an "INT=0 fix exposing a
downstream empty-vector indexing bug" — imprecise. The real location is
**groupr**, specifically `_groupr_nubar_records` calling `_interp_linlin`
with an empty `nu_e` (delayed nubar energies).

Traceback:

```
_interp_linlin  (groupr.jl:244)
_groupr_nubar_records  (groupr.jl:165)
groupr_module          (groupr.jl:53)
```

`_read_nubar` for MT=455 (delayed ν̄) with LNU=2 was reading the TAB1
immediately after the HEAD, but per ENDF-6 §1.5 — and Fortran `getyld`
at `groupr.f90:6472-6483` (label 110) — MT=455 has a **LIST of delayed-
neutron precursor decay constants** between the HEAD and the TAB1. JENDL
U-238 has NPL=6 (6 decay constants = 6 precursor families, standard).

Reading the TAB1 at the wrong file offset returned NP=0 → empty `tab.x`
and `tab.y` → empty `nubar_data.energies` → BoundsError.

MT=452 (total) and MT=456 (prompt) have TAB1 immediately after HEAD for
LNU=2; they were never affected.

## Fix

`src/orchestration/modules/groupr.jl::_read_nubar`, LNU=2 branch: skip a
LIST record before reading the TAB1 when `mt == 455`.

```julia
elseif lnu == 2
    # Tabulated. For MT=455 (delayed), a LIST of NNF decay
    # constants precedes the TAB1 (ENDF-6 §1.5; Fortran groupr.f90
    # getyld label 110 at line 6472).
    mt == 455 && read_list(io)
    tab = read_tab1(io)
    return (energies=collect(Float64, tab.x), values=collect(Float64, tab.y))
end
```

Gated on `mt == 455` — MT=452/456 paths unchanged, zero regression surface.

## Verification

**Pre-fix probe** (U-238 JENDL `J33U238`):

```
MT=452 → 7 energies, 7 values    ← LNU=2, TAB1 direct, OK
MT=455 → 0 energies, 0 values    ← LNU=2, TAB1 misread past LIST
MT=456 → 5 energies, 5 values    ← LNU=2, TAB1 direct, OK
```

**Post-fix probe**:

```
MT=455 → 4 energies, 4 values
  energies = [1.0e-5, 4.5e6, 9.0e6, 2.0e7]
  values   = [0.04634, 0.04634, 0.02853, 0.02853]
```

**Full-run T15 and T17**: both complete end-to-end, no crash.

**Regression subset** (T01/T02/T04, via `sweep_reference_tests.jl`):

| Test | Tape | Status | Δ vs full sweep |
|------|------|--------|-----------------|
| T01 | tape25 | NUMERIC_PASS 32812/32962 @ 1e-5 | = |
| T02 | tape28 | NUMERIC_PASS 12519/13873 @ 1e-5 | = |
| T04 | tape23 | NUMERIC_PASS 81/82 @ 1e-7 | = |
| T04 | tape24 | NUMERIC_PASS 56/74 @ 1e-5 | = |
| T04 | tape25 | DIFFS 107/119 | = |

Zero regression.

## Known follow-up (out of scope)

T15/T17 errorr MF33 covariance writer produces a **2.5 GB tape26** (30 M
lines). Dense full-matrix emit where Fortran writes sparse per-group
blocks. This is a correctness/perf issue in errorr — separate from the
crash fix. Does not affect pass/fail classification (DIFFS either way),
but does mean T15/T17 take ~10 min each in the sweep instead of ~500s of
broadr + a few seconds of errorr.

## Files changed

- `src/orchestration/modules/groupr.jl` — 4-line additive fix in
  `_read_nubar` LNU=2 branch.

## Trap (NEW)

**Trap N (MT=455 LIST)**: MF1/MT=455 (delayed ν̄) with LNU=2 has a LIST
of NNF precursor-family decay constants *between* the HEAD and the TAB1.
MT=452 (total) and MT=456 (prompt) have the TAB1 directly after the HEAD.
Any reader handling all three MTs must special-case MT=455 to skip the
LIST first. See Fortran `groupr.f90:6472-6483`, ENDF-6 §1.5.
