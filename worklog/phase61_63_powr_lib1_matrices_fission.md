# Phases 61-63 — powr lib=1 MF=6 matrices + fission

**Date:** 2026-05-02 (continuation of Phase 60 carbon-only port)
**Goal:** Round out powr lib=1 with the full MF=6 matrix paths (elastic +
inelastic + n2n) and the fission code path (MF=3/MT=18 sigma_f, MF=6/MT=18
chi + nu*sigf). All bit-identical to Fortran NJOY's tape50.

## Outcome

| Phase | Test                                       | Status                                        |
|-------|--------------------------------------------|------------------------------------------------|
| 60    | carbon abs-only (16 lines)                 | BIT-IDENTICAL (commit 8077d21)                |
| 61    | carbon + MF=6 elastic (84 lines)           | BIT-IDENTICAL (commit 96cb744)                |
| 62    | carbon + MF=6 elastic+inelastic (113 lines)| BIT-IDENTICAL (commit 96cb744)                |
| 63    | U-235 fission + chi + nu (53 lines)        | BIT-IDENTICAL (commit ea6da24)                |

T01 + mixr + resxsr standalone tests all unchanged across all four
incremental landings.

## Phase 61 — MF=6 elastic matrix

Carbon GENDF with `6 2/` added to the groupr deck. Output gains 68 lines
(34 P0 + 34 P1 elastic transfer matrix in packed band-diagonal storage).

Fortran-faithful behaviors:
- gamll detection: `lap = min(igc)`, `ldp = max(igc - ig2lo)`. With
  ngn=ngnd, icgrp(ig)=ig so the formulas use raw group indices.
- `xla(3) = ngnd - lap + 1`, `xld(3) = ldp`, `xlol(3) = xla(3) * (xld(3)+1)
  - lx*(lx-1)/2` where `lx = max(0, xla + xld - ngnd - 1)`.
- gamxs storage: per-row offset `lc = (ldp+1)*(jgp-1)` plus a triangular
  correction when `jg > kdp = ngnd+1-ldp`. Cell `loc = loce0 + lc + (jg2c - jg)`.
- P0 contribution: `data[offset + 1] * data[1] * cflux[jg]` with
  `offset = nl * ((izref - 1) + nz * (k - 1))` (sigma-zero index izref).
- P1 contribution: `3 * data[offset + 2] * data[2] * cfluxp1[jg]`
  (or cflux[jg] fallback when cfluxp1 zero), per gamxs.f90:1318.

Implementation: `_powr_pack_matrix(recs; kind=:elastic, ...)` — see
src/orchestration/modules/powr.jl.

## Phase 62 — MF=6 inelastic + n2n matrices

Same packing scheme as elastic, simpler accumulation:
- Inelastic (MT∈[51,91]): only P0 at first sigma-zero; `offset = nl*nz*(k-1)`.
- N2n (MT∈{6-9, 16, 46-49}): same as inelastic but divide accumulated value
  by 2 (gamxs.f90:1348 — n2n produces 2 neutrons).

Both kinds share the same `_powr_pack_matrix` body; the kind kwarg toggles
the index formula and the divide-by-2 factor. Output emission order matches
fast() lines 572-613: elastic P0 → elastic P1 → inelastic → n2n.

The carbon oracle exercises elastic + inelastic (MT=51, MT=91); n2n is
exercised by the same code path but no current oracle hits it (carbon has
no MF=6/MT=16 in the GENDF — would require a heavier nuclide).

## Phase 63 — Fission path (the most subtle)

U-235 GENDF with `6 18/` added to groupr. Adds three blocks to the output:
- chi (fission spectrum) — emitted FIRST, with a `(i6,i2,10a4)` header
  line (mergt-equivalent: kscr block prepended to nscr block).
- sigma_f — after absorption, before matrices.
- nu-bar — after sigma_f.

### The sigf-corruption quirk (worth pinning)

gamxs.f90:1206-1208 in the spectrum branch:
```fortran
jg2c = ngnd - ig2c - 1     ! sign mismatch — matrix branch uses +1
locc = locchi + jg2c - 1
a(locc) = a(locc) + cspc(i) * scr(loca) * scr(loca-1)
```

For ngn=ngnd and `i ∈ {ngnd-1, ngnd}`:
- i = ngnd-1 → jg2c = 0 → locc = locchi - 1 = locsf0 + ngnd - 1 = sigf[ngnd]
- i = ngnd     → jg2c = -1 → locc = locchi - 2 = locsf0 + ngnd - 2 = sigf[ngnd-1]

Fortran silently writes to the last two sigf cells. The corruption is small
(cspc[i] for high i is small) but visible — and reproduced, because
NJOY's reference output bakes it in. To match exactly, we use a single
unified `a[]` array and Fortran-style pointer arithmetic so the same
out-of-chi writes land in the same sigf cells.

### Other fission semantics

- MF=3 fission: MT=18 latch (`i318`); if MT=18 present, skip
  MT ∈ {19, 20, 21, 38} (Fortran gamxs.f90:1167-1168). Sigma_f goes to
  BOTH `locsf0` AND `locab0` (fission contributes to absorption).
- MF=6/MT=18 ig=0 record: cspc setup. Each `data[1 + nl*nz*(k-1)]` →
  `cspc[ig2lo + k - 1]` for k=1..ng2.
- MF=6/MT=18 ig>0, ig2lo=0: spectrum branch. For each i in 1..ngn,
  accumulate `cspc[i] * nu_sf * flux * cflux[jg]` into nus and
  `cspc[i] * nu_sf * flux` into chi[jg2c=ngnd-ig2c-1]. Cnorm gets
  `nu_sf * flux` per i (so cnorm = ngn × per-record sum).
- Post-loop: `nu = nus / sigf` per group; `chi *= 1/cnorm` per group.

### Bug fix landed alongside

`_powr_read_gendf_for_fast` was filtering MF=6 records on `ig > 0`,
discarding the cspc setup record (which has ig=0). Relaxed to `ig ≥ 0`
for MF=6 only — MF=3 still requires `ig > 0` (constant-XS records that
don't apply to powr).

## Phase B remaining items (each fail-loud at its guard)

- **Multi-temperature loop**: gamll currently only honors the reference
  temperature; ntp > 1 would need to walk all temps in the GENDF and
  emit per-temp arrays. (Phase 64 candidate.)
- **gamff self-shielding factors**: gated on `iwr = 1` (requires
  `nsigz > 1` AND `iff = 1`). Adds the locabs / locsf blocks and a
  separate sigma-zero × temperature × group output. (Phase 65.)
- **Delayed-neutron spectra**: `nfs > 0` with MF=5/MT=455 → adds
  delayed chi blocks via locdla. (Phase 66.)
- **iread=1 absorption-direct-read**: matd<0 in card 3. Fortran code
  for this branch is dead anyway (the outer `if (matd.gt.0)` filters
  it out before iread=1 can be set), so probably skip.

## Test commands

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/test_powr_lib1_standalone.jl     # Phase 60
julia --project=. test/validation/test_powr_lib1_elastic.jl        # Phase 61
julia --project=. test/validation/test_powr_lib1_matrices.jl       # Phase 62
julia --project=. test/validation/test_powr_lib1_fission.jl        # Phase 63
```

## Status: 23/23 modules dispatched + powr lib=1 ~70% covered

After Phase 63, powr lib=1 handles carbon (abs-only and abs+matrices) and
U-235 (fission). The remaining single-temperature, single-sigma-zero code
paths are all fenced off with loud errors; the multi-temperature and
shielding-factor paths are the next items. lib=2 (LIBRAR thermal) and
lib=3 (CLIB cpm) remain as Phase A scaffolds.
