# Phase 66 — powr lib=1 delayed-neutron spectra (MF=5/MT=455)

**Date:** 2026-05-03 (continuation of Phase 60-65 powr lib=1 port)
**Goal:** Port the delayed-neutron-spectrum path of powr lib=1: gamll's
`nfs += nl` accumulator, gamxs's MT=455 dispatch (both mfh=3 and mfh=5
branches), and fast()'s emission of `nfs` delayed-chi blocks after the
prompt chi.

## Outcome

| Phase | Test                                     | Status                                |
|-------|------------------------------------------|----------------------------------------|
| 66    | U-235 nfs=6 delayed neutrons (10611 B)   | **BIT-IDENTICAL on first run**         |

All 6 powr lib=1 standalone tests now pass:
- carbon abs-only (Phase 60), carbon + MF6 elastic (Phase 61), carbon +
  inelastic (Phase 62), U-235 fission (Phase 63), Fe-56 multi-T multi-σ₀
  (Phases 64-65), **U-235 delayed neutrons (Phase 66, this commit)**.

T01 regression-clean at NUMERIC_PASS 32812/32962 @ 1e-5 (unchanged).

## What landed

### `src/orchestration/modules/powr.jl` (~140 LOC delta)

1. **Reader extended** (`_powr_read_gendf_for_fast`):
   - MF=5/MT=455 records collected at the reference temperature into
     `mf5_455::Vector{NamedTuple}`. Each record carries `(ig, ig2lo, ng2,
     nl, nz, data)`. `nfs += nl_section` per record (gamll powr.f90:932-934).
   - MF=3/MT=455 records routed to a separate `mf3_455::Vector{NamedTuple}`
     instead of the regular `mf3` dict (data layout differs — MF=3/MT=455
     LIST data is `(d0, d1, d2)` per record, not the standard `(flux, xs)`
     pair). Empty in our oracle but populated defensively for future-proofing.

2. **Pointer layout corrected**:
   ```
   locchi = locsf0 + ngnd            (size ngnd, prompt only)
   locdla = locchi + ngnd            (size ngnd*nfs, delayed chi blocks)
   locnus = locdla + ngnd*nfs
   ```
   Previously `locnus = locchi + ngnd` (which works only for nfs=0).

3. **Delayed-chi accumulation** (mirror of gamxs.f90:1226-1262):
   - MF=3/MT=455 branch (mfh!=5): for each record, set
     `dnorm += d0*d1*d2`, latch `dnu` and `jgdnu` on first `d2 >= 0.1`,
     accumulate `a[locnus + jg - 1] += d1*d2*d0*cflux[jg]`.
   - MF=5/MT=455 branch (mfh==5): per (ig, k=2..ng2, jgt=1..nl):
     - `a[locchi + jgc - 1] += val * dnorm`
     - `a[locdla + ngnd*(jgt-1) + jgc - 1] += val * dnu`
     - **OVERWRITE** at `jgc == ngnd`: `a[locdla + ngnd*(jgt-1) + ngnd - 1]
       = a[jgt + locb_offset]` where `locb_offset = 6` reads from the
       absorption region (Fortran's `a[jgt+locb]` quirk — see "Quirk
       worth pinning" below).
     - `sumd[jgt] += val * dnu` when `jgc < ngnd-1`.
     - `cnorm += val * dnorm`.

4. **Final-pass normalization** (mirror of gamxs.f90:1366-1393):
   - Prompt chi divided by `cnorm` (1/Σ chi).
   - If `dnu != 0`: delayed chi scaled by `rnorm = 1/a[locnus + jgdnu - 1]`
     for ig=1..ngnd-2; replaced with `sumd[il]*rnorm` at ig=ngnd-1; cell
     at ig=ngnd untouched (stays at locb-replacement value).

5. **Writer extended** (`_powr_fast_one_material` chi block emission):
   - After the prompt chi (`(i6,i2,10a4)` header with `i2=0` + chi data),
     emit `nfs` delayed-chi blocks each with `i2 = ifs ∈ {1..nfs}` and
     chi data drawn from `a[locdla + ngnd*(ifs-1) .. locdla + ngnd*ifs - 1]`.

6. **Loud-TODO removed**: the previous `nfs > 0 → error` guard at
   `_powr_fast_one_material` is replaced by the proper implementation. The
   guard for unsupported MF=6 MTs remains.

### `test/validation/test_powr_lib1_delayed.jl` (new, 130 LOC)

6-assertion bit-identical test against the U-235 + delayed-neutrons oracle.
Includes the regen recipe in the footer.

### `test/validation/oracle_cache/powr_lib1_delayed/` (gitignored)

U-235 ENDF/B-VIII oracle with `5 455/` explicitly added to groupr's card 7
sequence (groupr's auto-finder excludes MF=5; explicit request mandatory):
- `tape20`: `n-092_U_235-ENDF8.0.endf` (37 MB).
- `tape21..tape23`: PENDF chain (reconr, broadr, groupr).
- `tape50`: Fortran NJOY reference (131 lines, 10611 bytes — 1 prompt chi
  + 6 delayed chi blocks + the standard nscr block).

## One Fortran quirk worth pinning

**`a[jgt+locb]` reads from the OUTPUT array, not from `scr`.** At
gamxs.f90:1255:
```fortran
if (jgc.eq.ngnd) a(ngnd*(jgt-1)+locd) = a(jgt+locb)
```
With `locb = l + lz - 1 = 6` (l=1, lz=6), this is `a[7..12]` — sitting
INSIDE the absorption region (locab0=1..68 for ngnd=68). For our U-235
oracle, those cells happen to hold the capture XS at jg=7..12 after the
MF=3 absorption loop has run. The Fortran is essentially copying capture
XS values into the last cell of each delayed-chi block. The agent's
deep-read framing as "time constants" was incorrect — these values come
straight from the absorption accumulator.

For bit-identical reproduction, our Julia port reads from `a` (not `data`)
at the same offset. We confirmed against the oracle: each delayed-chi
block is all zeros except the last cell, which holds 1.342, 1.327, 1.301,
1.276, 1.248, 1.273 respectively (= U-235 capture XS at jg=7..12).

## Why it landed first try

- LAW 2 (Fortran first): direct read of powr.f90:1056-1058, 1061-1066,
  1135, 1224-1264, 1365-1393 caught the locdla pointer correction (Phase
  64 worklog had this wrong) and the `a[jgt+locb]` quirk (initial agent
  framing was misleading) BEFORE any code was written.
- LAW 1 (oracle-driven TDD): the oracle was built and the failing test
  was confirmed before the implementation started. Since dnu=dnorm=0 in
  our oracle (no MF=3/MT=455), most of the gamxs MT=455 logic is a no-op
  — only the locb-replacement contributes. The implementation faithfully
  preserves the no-op paths so that future oracles with MF=3/MT=455 will
  exercise them correctly.

## Test commands

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/test_powr_lib1_delayed.jl     # Phase 66 (NEW)
julia --project=. test/validation/test_powr_lib1_standalone.jl  # Phase 60
julia --project=. test/validation/test_powr_lib1_elastic.jl     # Phase 61
julia --project=. test/validation/test_powr_lib1_matrices.jl    # Phase 62
julia --project=. test/validation/test_powr_lib1_fission.jl     # Phase 63
julia --project=. test/validation/test_powr_lib1_selfshielding.jl  # Phases 64-65
julia --project=. test/validation/reference_test.jl 1           # T01 regression
```

## Status: powr lib=1 ~98% covered

After Phase 66, powr lib=1 handles all observed code paths bit-identical:
absorption, MF=6 matrices (elastic / inelastic / n2n), fission (MF=3 +
MF=6 chi/nu), multi-temperature + multi-σ₀ self-shielding, and now
delayed-neutron spectra. Remaining: the `matd<0` (iread=1)
"absorption-only direct read" corner — the worklog inherited from Phase
60 notes this is dead code in practice (the outer `if (matd.gt.0)` filter
in fast() prevents iread=1 from being set), so it can be left as a loud
TODO indefinitely. lib=2 (LIBRAR thermal) and lib=3 (CLIB cpm) remain
Phase A scaffolds.
