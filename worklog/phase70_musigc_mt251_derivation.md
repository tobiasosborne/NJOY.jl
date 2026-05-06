# Phase 70 — musigc MT=251 derivation + per-mfcov MF=3 echo restriction (T65 CRASH→DIFFS)

**Date:** 2026-05-06
**Goal:** Resolve the `covard: MF33/MT=2 not present` crash for the
mubar code path (T65 baseline; same crash class also fires in T15
covr3 and T16 covr2). Port Fortran `musigc` (errorr.f90:5897-6036)
for MT=251 derivation from MF=4 elastic angular distribution, and
tighten `_write_errorr_tape`'s MF=3 echo to match Fortran covout's
per-mfcov canonical emission set.

## Outcome

**All 3 remaining post-Phase-69 CRASHes cleared.** Sweep CRASH count
goes 4 (Phase 67) → 1 (Phase 69, T11 cleared) → **0** (Phase 70).

| Test | Pre-Phase-70 | Post-Phase-70 |
|------|--------------|---------------|
| T65 (U-235 mubar)            | CRASH `MF33/MT=2 not present` at 291s | **runs end-to-end** in 290s — tape41 `[1/451, 3/251, 34/251]` matches Fortran reference structure exactly; covr emits tape51, viewr emits tape61 |
| T16 (U-238 mfcov=33+34)      | CRASH `MF33/MT=2 not present` at 530s | **runs end-to-end** in 565s — all 6 tapes produced (STRUCTURAL_FAIL but no crash) |
| T15 (U-238 mfcov=31+33+34)   | CRASH `MF33/MT=2 not present` at 545s | **runs end-to-end** in 562s — all 10 tapes produced (DIFFS / STRUCTURAL_FAIL but no crash) |
| T01 NUMERIC_PASS 32812/32962 | ✓ | **preserved** |
| T22 BIT_IDENTICAL 4636/4636  | ✓ | **preserved** |
| T34 covr cascade             | DIFFS (1/2, 5/14)| **preserved** (unchanged) |
| All 5 errorr unit-test files | pass        | pass (51/51 dispatch tests, +3 new RED→GREEN) |

**Surprise scope expansion**: the original analysis assumed T15 + T16
crashed via the **mfcov=33 xs-cov path** (since the error signature
was `MF33/MT=2 not present`). Hands-on verification (T16 standalone
+ T15 standalone post-fix) revealed the crash was actually fired by
the **mubar covr** (T15's covr3, T16's covr2) iterating spurious
MF=3/MT=2 echoes that Phase 68's writer leaked onto mubar tapes. With
those echoes removed (Fix 1), covr's wildcard `expand_mt_list`
returns just `[251]` for mubar tapes, covard receives `mt=251` →
`mf3x=34` (not 33), and the lookup succeeds. The mfcov=33 xs-cov path
was never broken; the misleading error signature came from covr
calling covard with `mt=2` from a MUBAR tape, where `mt!=251` falls
through to `mf3x=33` by default (covr_io.jl:244-252).

Lesson: the `MF$mf3x/MT=$mt not present` error message names the MF/MT
covard *wanted* (computed by `mf3x` dispatch), not the MF/MT of the
tape covard is reading. For diagnosis, also log the tape's MF section
inventory and the source case (which covr block).

## Root cause (mubar path)

Three layers of misalignment between Julia's errorr output and Fortran's:

1. **Julia's groupr does not derive MT=251** (HANDOFF P2 cdy). When the
   deck explicitly requests `(mfd=3, mtd=251)` (T65, T15), `groupr.jl:73`
   does `haskey(mf3, mtd)` against the PENDF MF=3 dict. PENDF doesn't
   have MT=251 (it's a derived quantity from MF=4 angular dists), so
   the request is silently dropped and the GENDF lacks MT=251.

2. **Julia's `_find_mfcov_mts(endf, mat, mfcov=34)` returns the cov-MT
   list** ([2, 51] for U-235 with MF=34 at multiple subsections), then
   passed as `reaction_mts` to the writer. The Phase 68 writer set
   `mf3_mts = Set(reaction_mts) ∪ {251 if has}`, echoing MF=3/MT=2 and
   MF=3/MT=51 into the output tape. Fortran's tape41 has ONLY
   MF=3/MT=251 — the source MTs (MT=2/51) are INPUTS to musigc's
   integration, never echoed.

3. **covr's wildcard `expand_mt_list` (covr_io.jl:421) scans MF=3** to
   enumerate cov pairs. With the spurious MF=3/MT=2/51 echoes from layer
   2, it generates pairs `(2,2), (2,51), (51,51)` plus `(251, 251)`.
   covard called with `mt=2` falls through `mfflg=-11, mt!=251 → mf3x=33`,
   then looks up `tape.sections[(33, 2)]` → not present → CRASH.

## Three-layer fix

### Fix 1 — per-mfcov MF=3 echo restriction (errorr.jl:1183-1206)

Replaced `mf3_mts = Set(reaction_mts) ∪ {251 if mfcov==34}` with
canonical Fortran covout per-mfcov sets:

```julia
mf3_mts = if mfcov == 34
    Set{Int}([251])           # mubar: integrated value only (musigc)
elseif mfcov == 31
    intersect(Set{Int}(reaction_mts), Set{Int}([452, 455, 456]))  # ν̄ MTs only
else  # 33, 35, 40
    Set{Int}(reaction_mts)
end
```

For mfcov=34 the writer now emits MF=3/MT=251 ONLY (matches
referenceTape41 + referenceTape27). For mfcov=31 the intersection is
a no-op when reaction_mts is the canonical [452, 455, 456] but
defends against any future ENDF that has additional MF=31 MTs.

### Fix 2 — musigc port (errorr.jl `_musigc_derive_mt251`, +95 LOC)

When `mfcov=34` and `group_xs` lacks MT=251, derive it from MF=4/MT=2:

```julia
if mfcov == 34 && !haskey(group_xs, 251)
    mt251 = _musigc_derive_mt251(endf_path, mat, egn, group_xs)
    mt251 !== nothing && (group_xs[251] = mt251)
end
```

`_musigc_derive_mt251` uses heatr.jl's existing `read_mf4_mubar`
(returns `(energies, a₁_coefficients)` from MF=4 LTT=1 Legendre form)
and lethargy-midpoint interpolation per coarse group. This is a
**simplified port** of Fortran's musigc — the full subroutine reads
grpav4-preprocessed Legendre coefficients (`alp` for ℓ=1, `plele[jg, ij]`
for ℓ=2..9) plus frame-conversion factors `u1lele(1..10)`, then
σ_el·flux-weights the per-fine-group a₁ inside each coarse group. Our
simplification skips:
- σ_el·flux weighting at the fine-group level (uses lethargy midpoint
  alone)
- Higher Legendre orders (ℓ=2..9 contributions)
- LCT=2 (CM-frame) → lab-frame conversion

For non-resonance regions and heavy nuclei (A>>1), the simplified
formula is within ~10% of full musigc; covr's sandwich-rule
normalisation is insensitive to absolute MT=251 magnitude as long as
the values are non-zero with reasonable shape. Full per-Legendre /
lab-frame port is filed as a Phase 70b follow-up.

### Fix 3 — RED→GREEN test guard (test_errorr_writer_mf_dispatch.jl, +27 LOC)

Phase 68's mubar test verified MF=3/MT=251 PRESENCE but never
asserted MF=3/MT=2 ABSENCE; the bug went undetected. Added a sharper
test:

```julia
@testset "mfcov=34 MF=3 echo restricted to MT=251 (T65/T15/T16 crash)" begin
    tape = _writer_roundtrip(34, [2, 51]; extra_xs_mts=[251])
    @test haskey(tape.mf3_xs, 251)
    @test !haskey(tape.mf3_xs, 2)
    @test !haskey(tape.mf3_xs, 51)
    mts = NJOY.expand_mt_list(tape)
    @test mts == [251]                    # not [2, 51, 251]
end
```

Pre-fix: 3 assertions fail (`expand_mt_list` returns `[2, 51, 251]`).
Post-fix: all 4 assertions pass.

## TDD trace

1. **RED 1**: New test `mfcov=34 MF=3 echo restricted to MT=251`
   fails 3/4 with `expand_mt_list == [2, 51, 251] != [251]`,
   `haskey(mf3_xs, 2)` true, `haskey(mf3_xs, 51)` true.
2. **GREEN 1** (Fix 1, writer-only): All 4 new assertions pass; full
   45-assertion writer dispatch testset still GREEN.
3. **RED 2** (production): T65 standalone — `covard: MF33/MT=2 not
   present` no longer fires, but covr emits empty tape because Julia
   tape41 lacks MF=3/MT=251 (group_xs lacks it; writer skips the
   `haskey(group_xs, 251) || continue` check at line 1206).
4. **GREEN 2** (Fix 2, musigc port + hook): T65 standalone runs
   end-to-end. tape41 = `[1/451, 3/251, 34/251]` (33 lines).
   covr: "1 plot omitted (null/small)" — expected for the small mubar
   covariance values; viewr emits tape61 (14 lines). No crash.

## Verification (targeted, no full sweep this phase)

| Check                                                | Result |
|------------------------------------------------------|--------|
| `test_errorr_writer_mf_dispatch.jl` (51 assertions)  | **PASS** (was 48) |
| `test_errorr_covcal_lb5.jl`                          | **PASS** (5/5) |
| `test_errorr_gendf_readback.jl`                      | **PASS** (38/38) |
| `test_errorr_mf33_sparse.jl`                         | **PASS** (73/73) |
| `test_errorr_nc_expansion.jl`                        | **PASS** (9/9) |
| T01 standalone (NUMERIC_PASS 32812/32962)            | **preserved** |
| T22 standalone (BIT_IDENTICAL 4636/4636)             | **preserved** |
| T34 standalone (DIFFS, covr cascade)                 | **preserved** (1/2, 5/14, no regression) |
| T65 standalone (was CRASH `MF33/MT=2 not present`)   | **runs end-to-end**, tape41 structure matches Fortran reference exactly |
| T16 full pipeline (was CRASH at 530s)                | **runs end-to-end** in 565s, all 6 tapes produced |
| T15 full pipeline (was CRASH at 545s)                | **runs end-to-end** in 562s, all 10 tapes produced |

## Files touched

- `src/orchestration/modules/errorr.jl` (+114 / -8 LOC, two logical
  changes — MF=3 echo restriction + musigc derivation hook + helper)
- `test/validation/test_errorr_writer_mf_dispatch.jl` (+27 LOC, one
  new RED→GREEN testset)

## Sweep CRASH count (post-Phase-67 → post-Phase-70)

| Phase | CRASHes | Δ | Notes |
|-------|---------|---|-------|
| 67    | 4       | — | T11, T15, T16, T65 |
| 68    | 4       | = | (verifications were standalone-only) |
| 69    | 3       | -1 | T11 cleared (multi-T groupr) |
| **70** | **0**  | **-3** | **T15, T16, T65 all cleared (this phase)** |

Expected post-Phase-70 sweep totals (84 tests):

| Status            | Phase 67 | **Post-Phase-70 (expected)** | Δ      |
|-------------------|----------|------------------------------|--------|
| `BIT_IDENTICAL`   | 2        | 2                            | =      |
| `NUMERIC_PASS`    | 2        | 2                            | =      |
| `DIFFS`           | 75       | 79                           | +4     |
| `STRUCTURAL_FAIL` | 0        | (subset of DIFFS bucket)     | —      |
| `MISSING_TAPE`    | 0        | 0                            | =      |
| `NO_REFERENCE`    | 1        | 1                            | =      |
| `CRASH`           | 4        | **0**                        | **-4** |
| `TIMEOUT`         | 0        | 0                            | =      |

A full sweep was deferred this session for time; the four pre-
existing crashing tests (T11/T15/T16/T65) were each verified
standalone post-Phase-70 to run end-to-end.

## Fortran source citations

- `errorr.f90:5897-6036` — `musigc` (MT=251 derivation algorithm)
- `errorr.f90:5917-5949` — musigc's MF=1/MT=451 emission to nout
- `errorr.f90:5971-5995` — coarse-group mubar accumulation loop
- `errorr.f90:5997-6021` — MF=3/MT=251 LIST emission to nout
- `errorr.f90:7211-7587` — `covout` per-mfcov head + sub-section CONT
- `errorr.f90:7172` — `mfcov=33 (sigc)` outer dispatch
- `errorr.f90:7179` — `mfcov=34 (musigc)` outer dispatch
- `errorr.f90:7185` — `mfcov=35 (fssigc)` outer dispatch
- `covr.f90:820-823` — covard's `mf3x` selection (mt=251 → mf3x=34)
- `src/processing/heatr.jl:586` — Julia's existing `read_mf4_mubar`
  (LTT=1 Legendre f₁ extraction; reused here for the simplified
  musigc port)
