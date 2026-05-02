# Phase 59 — mixr + resxsr + powr scaffold (23/23 modules wired)

**Date:** 2026-05-02 (continuation of T11/wimsr Phase 58 session)
**Goal:** Close the "unwired modules" gap. After Phase 58, 20/23 NJOY modules
were dispatched in `pipeline.jl`; this phase brings the count to 23/23.

## Outcome

| Module  | Status this phase                                                   | LOC (Julia) |
|---------|---------------------------------------------------------------------|-------------|
| mixr    | **BIT-IDENTICAL** standalone test (carbon weight=1.0 + self-mix 0.5+0.5) | ~280 |
| resxsr  | **BIT-IDENTICAL** standalone test (carbon @296K, eps=0.001, first try) | ~340 |
| powr    | **Phase A scaffold** — parser + dispatch wired; module errors loudly per CLAUDE.md Rule 6 ("powr lib=N not yet ported"). Phases B-D follow. | ~50 |

T01 regression-clean: tape25 stays at NUMERIC_PASS 32812/32962 @ rtol=1e-5.

All three new dispatches are listed in `AVAILABLE_MODULES` (input_parser.jl)
and have an `elseif mc.name == :…` branch in `pipeline.jl`.

## Files added / modified

- `src/orchestration/modules/mixr.jl` (NEW, ~280 LOC): full port of
  `njoy-reference/src/mixr.f90:17-390`. `mixr_module`, MF=1/MT=451 builder
  (iverf-aware: 4/5/6 layouts, NSUB lookup from AWI), MF=3 builder (TAB1
  with INT=2 lin-lin, sigfig(7,0) on XS), custom tape writer (Fortran's
  `dictio + afend` skips MT=451 SEND so we can't compose with the shared
  `write_pendf_tape`), lin-lin interpolator that returns 0 outside src range.
- `src/orchestration/modules/resxsr.jl` (NEW, ~340 LOC): full port of
  `njoy-reference/src/resxsr.f90:10-502`. `resxsr_module`, per-material
  union grid construction across MT=2/18/102 + multi-temp interpolation,
  3-point linear-deviation thinner mirroring Fortran's exact restart
  pattern (including the end-of-grid stack-discard quirk), CCCC-IV binary
  writer reusing `_write_record` / `_record_buf` from `src/formats/ccccr.jl`.
- `src/orchestration/modules/powr.jl` (NEW, ~50 LOC): Phase A scaffold —
  recognises lib ∈ {1, 2, 3}, errors loudly with mode-named TODO. Header
  documents the multi-phase plan (B=fast/lib=1, C=therm/lib=2, D=cpm/lib=3).
- `src/orchestration/input_parser.jl`: `MixrParams` + `parse_mixr` (6 cards),
  `ResxsrParams` + `parse_resxsr` (variable cards depending on nholl/nmat),
  `PowrParams` + `parse_powr` (cards 1-2; mode-specific cards retained as
  raw_cards for future phases). All three names added to AVAILABLE_MODULES.
- `src/orchestration/pipeline.jl`: dispatch branches for `:mixr`, `:resxsr`,
  `:powr`.
- `src/NJOY.jl`: includes + exports for the three new modules.
- `test/validation/test_mixr_standalone.jl` (NEW): 14 assertions across two
  oracles (single-input weight=1.0 and two-input self-mix 0.5+0.5).
- `test/validation/test_resxsr_standalone.jl` (NEW): 16 assertions, full
  byte-for-byte equality against the Fortran-generated CCCC binary.
- `test/validation/oracle_cache/mixr_standalone/` and
  `…/resxsr_standalone/` (gitignored, regen recipes inside the test files):
  Fortran NJOY reference fixtures.

## Fortran-faithful quirks worth pinning

### mixr
1. **TPID is 66 blanks**, not the description text. Fortran writes
   `math=1; afend(nout, 0)` (mixr.f90:198-199) — `afend` with mat=1 emits
   a blank-data marker line, the description never reaches the TPID.
2. **No SEND between MT=451 dictionary and MF=1 FEND.** Fortran's `dictio`
   does not auto-emit SEND, and mixr.f90:255 calls `afend` directly. The
   shared `write_pendf_tape` always emits SEND after MF=1 → had to write
   the tape inline in `_mixr_write_tape`.
3. **Dictionary entries write `(22x, 4i11)`**, not floats: cols 1-22 are
   blank (not "          0          0"). Caught immediately by L6 diff
   in the first oracle run.

### resxsr
4. **Description goes nowhere.** The user's `holl` cards are read but
   never written. Fortran resxsr.f90:460-462 forces the set-hollerith
   record to all blanks. The user description shows up only in the run
   log, not on the binary tape.
5. **Material name read with `(a6)` then stored in `k8` slot.** For
   "cnat" (4 chars) the on-disk pattern is `cnat    ` (4 chars + 4 trailing
   blanks). `rpad(name, 8)` matches.
6. **Thinning intentionally drops uncommitted stack contents at end of
   grid.** Fortran resxsr.f90:366 (`if (ie.eq.ne) go to 350`) emits only
   the very last input point as the final thinned point — anything still
   in the un-flushed stack from the previous break is silently lost. We
   mirror exactly so the byte-for-byte total nener matches.

### powr (scaffold)
- Parser validates `lib ∈ {1, 2, 3}` and errors otherwise. The cards 3+
  layouts diverge wildly per lib mode — keeping them as `raw_cards` until
  Phases B-D port the per-mode drivers (`fast`, `therm`, `cpm`).

## Test commands

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/test_mixr_standalone.jl
julia --project=. test/validation/test_resxsr_standalone.jl
julia --project=. test/validation/reference_test.jl 1     # T01 regression
```

## Phase 60+ — powr port plan

1. **Phase B (lib=1, fast/GAMTAP)** — port `fast`, `pgam`, `pgamff`, `mergt`,
   `gamll`, `gamxs`, `gamff` (~1000 LOC Fortran → ~600 LOC Julia).
   Build a self-contained oracle: GENDF(carbon @296K) + minimal powr deck
   with `lib=1`. Capture Fortran tape50 binary as reference. Standalone
   test in `test/validation/test_powr_lib1_standalone.jl`.
2. **Phase C (lib=2, thermal/LIBRAR)** — port `therm`, `tprnt`, `packa`
   (~500 LOC).
3. **Phase D (lib=3, cpm/CLIB)** — port `cpm`, `pinit`, `sfile2`, `sfile3`,
   `sfile4`, `sfile6` and the `sf*out` writers (~2000 LOC; biggest by far).
4. Each phase: red bar via standalone oracle test → port → bit-identical →
   commit + worklog entry.
