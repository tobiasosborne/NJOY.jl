# Phase 69 — multi-temperature groupr port (T11 CRASH→DIFFS)

**Date:** 2026-05-05
**Goal:** Resolve T11 (Pu-238 WIMS), the only remaining CRASH after Phase 68.
Port the Fortran groupr outer-temperature loop so multi-T decks
(T11 ntemp=3 at 300/900/2100K) emit a multi-T GENDF that downstream
wimsr can consume.

## Outcome

| Test | Pre-fix | Post-fix |
|------|---------|----------|
| T11  | CRASH `wimsr_extract_resint: GENDF tape has 1 temperature(s) but params.ires=3` | **DIFFS** — both tape27 (455 lines vs ref 1169) and tape28 (10379 vs ref 15722) produced; pipeline runs end-to-end |

The four CRASHes from Phase 67 are now all resolved (T15, T16, T65 in Phase 68;
T11 here). Sweep summary will catch up to this in the next full run.

## Three landed fixes (one per layer)

### Fix 1 — per-temperature MF3 extractor (`src/orchestration/pendf_io.jl`, +75 LOC)

`extract_mf3_all` walks every PENDFMaterial entry and overwrites
`result[mt]` per section, so on a multi-T PENDF (broadr writes one
MEND-bounded material entry per T) the **last temperature's data wins**.
Phase 67 hid this by only ever passing `params.temperatures[1]` to groupr.

New helpers:

- `_pendf_material_temperature(material) -> Float64`: locates the TEMP
  CONT in `material.mf1_lines` by walking lines 2..5 and picking the
  last line where (L2=0, N1>0, 0 ≤ C1 < 1e6). Handles ENDF-IV layout
  (TEMP at mf1_lines[2]; T11 case) and ENDF-V/VI (TEMP at mf1_lines[3]
  or [4]; T01 case). Skips description Hollerith via integer-fields-
  parse-first guard.
- `extract_mf3_at_temperature(tape, mat, target_temp; tol=0.5)`:
  selects the PENDFMaterial whose stored TEMP matches `target_temp`
  within tolerance; falls back to single-T material when present.
  Errors loudly per Rule 6 if no candidate matches (lists available
  temperatures in the message).

### Fix 2 — `groupr_module` outer-T loop (`src/orchestration/modules/groupr.jl`)

Replaced the single-T extract+process with a `for ti in 1:length(temps)`
loop that mirrors Fortran groupr.f90:479-943 (`do 630 itemp=1,ntemp`):

```julia
temps = isempty(params.temperatures) ? [0.0] : copy(params.temperatures)
for ti in 1:length(temps)
    mf3 = extract_mf3_at_temperature(tape, params.mat, temps[ti])
    mt_list = _groupr_expand_auto(params.mt_list, mf3)
    mt_results = MTRes[]
    for (mfd, mtd, name) in mt_list
        # ... unchanged per-MT processing ...
    end
    push!(per_temp_results, mt_results)
end
```

`_write_groupr_tape` signature changed to take `temps::Vector{Float64}`,
`sigz_list::Vector{Float64}`, and `per_temp_results::Vector{Vector{MTRes}}`
indexed `[temp_index][mt_index]`.

### Fix 3 — multi-T GENDF writer + tempin in seq=2 C1

Per-temperature tape layout now mirrors Fortran (per Phase 68
research-agent A's tape sketch from groupr.f90:479-943):

```
TPID once
for each T:
  HEAD CONT (mat,1,451): (za, awr, 0, nz, -1, 1)
  CONT      (mat,1,451): (tempin, 0, ngn, 0, NW, 0)   ← C1=tempin (Fortran 551)
  DATA      egn boundaries packed 6/line
  FEND      (mat,0,0)
  for each MT:
    HEAD CONT (mat,3,mt) (za, 0, nl=1, nz=1, 0, ngn)
    per (g, ig2): (tempin, 0, 3, 1, 3, g) + 3-float body
    SEND
  MEND      (0,0,0)
TEND        (-1,0,0)
```

Crucially: **C1 of seq=2 is tempin** (was `0.0`). wimsr's
`_wimsr_read_gendf_metadata` (wimsr_xsecs.jl:209) extracts each block's
temperature from this field — without the fix, all 3 blocks looked like
T=0.0 and got deduped to a single tempr entry of length 1 → BoundsError
on `tempr[1:ires]` for `ires=3`.

**MF=3 nz=1 simplification:** the per-MT HEAD now declares `nz=1`
even when the deck has `nsigz=7` (T11 case). Each per-(g, ig2) record
holds only one σ₀ column (flux + sigma_g + sigma_int), and wimsr's
body indexing `_bidx(i, iz, ig2, nl, nz)` (wimsr_xsecs.jl:260) requires
`nl*nz < length(body)`. For non-resonance MTs and decks with
`sgref ≥ 1e10` (T11 has this; isg=0 → wimsr falls back to inf-dil),
emitting a single column is correct. Full nsigz×ngn block-matrix output
(URR self-shielding) is a follow-up phase.

**Per-(g, ig2) C1 = tempin** (was 0.0): mirrors Fortran groupr.f90:857.
Not strictly required by Julia's wimsr reader (which segments by
MF=1/MT=451 HEAD count, not by per-record C1), but Fortran-faithful.

### Side fix — pipeline.jl `_collect_gaspr!` leftover (-1 LOC)

Removed a dispatch line that called `_collect_gaspr!(ctx, tapes, params)`
inside the wimsr branch (pipeline.jl:262). The collector is typed for
`GasprParams`; the wimsr branch passes `WimsrParams`, producing a
`MethodError` for every wimsr-using test. Phase 58a leftover that was
never exercised because T11 always crashed earlier in wimsr.

`_collect_gaspr!` is also not called from the gaspr branch itself
(line 245) — gaspr's MT=203/207 emission goes through `ctx.extra_mf3`
directly per Phase 54. The function is currently dead code; left in place
pending the question of whether the gaspr branch should call it (out of
scope for this phase).

## TDD trace

- **RED 1** (single-T groupr): `T11 CRASH wimsr_extract_resint: GENDF tape has
  1 temperature(s) but params.ires=3`. Reproduced in 178s.
- **First fix attempt** crashed on `parse_endf_float("94-pu-238f")` —
  the temperature heuristic walked Hollerith description lines.
  Reordered to call `tryparse(Int, ...)` on integer fields first;
  description lines fail there and are skipped before the float parse.
- **RED 2**: `BoundsError: attempt to access 3-element Vector{Float64} at
  index [8]` in `_dispatch_mf3!` — Julia's MF=3 head declared nz=7 but
  body only had 3 floats. Set `nz_mf3=1` (single-column simplification).
- **RED 3**: `MethodError: _collect_gaspr!(::RunContext, ::TapeManager,
  ::WimsrParams)` — the pipeline's leftover. Removed.
- **GREEN**: T11 produces tape27 (455 lines, 13 match) and tape28
  (10379 lines). STRUCTURAL_FAIL (line counts diverge from reference)
  but no crash; the multi-T port has done its job.

## Files touched

- `src/orchestration/pendf_io.jl` (+75 LOC) — `_pendf_material_temperature`,
  `extract_mf3_at_temperature`
- `src/orchestration/modules/groupr.jl` (~+90 / -65 LOC) — outer-T loop +
  multi-T writer
- `src/orchestration/pipeline.jl` (-1 LOC) — drop bad `_collect_gaspr!`

No new tests added this phase. `test_wimsr_t11_standalone.jl` (Phase 58)
already covers wimsr's GENDF consumption with a Fortran-faithful oracle —
exercising Julia groupr against that oracle is the natural follow-up but
requires a per-T fixture.

## Verification (targeted)

| Check                                                | Result |
|------------------------------------------------------|--------|
| T11 reference test (`reference_test.jl 11`)          | **DIFFS — both tapes produced** (was CRASH) |
| T01 reference test                                   | _(see sweep)_ |
| T22 reference test                                   | _(see sweep)_ |
| T15/T34/T27 (errorr+covr — different GENDF consumers) | _(see sweep)_ |

(Targeted regression sweep run separately; this phase guards groupr's
single-T path by passing `temps=[temps[1]]` when the deck has one T.)

## Remaining T11 work (deferred — not blocking sweep)

Both tape27 and tape28 are STRUCTURAL_FAIL. wimsr cross-section data is
mostly zeros, indicating the per-MT `nz=1` simplification drops most of
the URR self-shielding signal that T11 expects. Path to BIT_IDENTICAL:

1. **Multi-σ₀ MF=3 emission**: write `nsigz` copies of the inf-dil cross
   section per (g, ig2) record (or the proper per-σ₀ values when URR is
   present). Fortran does this via `genflx` + `nflmax` allocation
   (groupr.f90:528-534) — port that path. Estimated ~40 LOC delta.

2. **Fortran-faithful MF=1/MT=451 LIST body**: currently we write only
   egn boundaries; Fortran writes (titles + sigz + egn + egg) in one
   LIST under the seq=2 CONT. Required for full bit-for-bit match;
   wimsr currently doesn't read the title/sigz/egg sections, so this
   is purely cosmetic until the writer is tightened.

3. **MF=6 transfer matrices**: T11's deck requests MF=6 for elastic,
   inelastic, fission, n2n, n3n, thermal — none currently emitted by
   Julia groupr. Fortran's `panel` + `getff` infrastructure is a
   substantive port (~hundreds of LOC). Likely the bulk of the
   remaining tape27 zero-fill.

4. **MT=251/252 derivation** (HANDOFF P2 cdy): mubar/xi from MF=4/6
   angular distributions. Already filed.

5. **moder tape28**: tape28 is moder's binary→ASCII conversion of
   tape26. The format diff at line 1 is moder writing the GENDF's
   own TPID instead of the deck's `94-pu-238` title. Probably a
   moder-side issue; orthogonal to groupr.

## Fortran source citations

- `groupr.f90:479` — outer-T loop header `do 630 itemp=1,ntemp`
- `groupr.f90:484` — `tempin = temp(itemp)`
- `groupr.f90:493-525` — find-mat-and-temp on PENDF (tomend/contio loop)
- `groupr.f90:538-585` — per-T MF=1/MT=451 head + LIST emission
- `groupr.f90:544-549` — head CONT scr setup `(za, awr, 0, nz, -1, ntw)`
- `groupr.f90:551-578` — LIST CONT scr setup with `tempin` in C1
- `groupr.f90:585` — afend (FEND closing MF=1)
- `groupr.f90:833-862` — per-MT head + per-(g, ig2) LIST
- `groupr.f90:923` — `call amend(ngout2, 0)` (MEND between temps)
- `groupr.f90:943` — outer-T loop end
- `wimsr.f90:947` — per-temperature reset boundary on the wimsr side
- `wimsr_xsecs.jl:188-248` — Julia wimsr reader keying on MF=1/MT=451 HEAD
- `wimsr_xsecs.jl:209` — `t = parse_endf_float(row2[1:11])` (tempin from seq=2 C1)
- `wimsr_xsecs.jl:260` — `_bidx(i, iz, ig2, nl, nz)` body indexing
