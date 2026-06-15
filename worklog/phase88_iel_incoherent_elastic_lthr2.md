# Phase 88 — thermr incoherent-elastic (iel, LTHR=2) wired end-to-end (bead NJOY_jl-69j Stage 2)

Date: 2026-06-16

## Goal
Stage 2 of `NJOY_jl-69j`: wire the incoherent-elastic (`iel`, LTHR=2) path end-to-end
so the thermr PENDF (`tape60`) gains MF3/MT222 + MF6/MT222 for **T74** (H-in-ZrH) and
**T69** (Zr-in-ZrH). LTHR=2 only (LTHR=3 / T67 MT223 is a separate later stage).
Stage 1 (the pure-function layer: `incoh_elastic_xs`, `read_mf7_mt2` LTHR=2 returning
`IncoherentElasticData`, `incoh_dwp_at`, `build_incoh_elastic_records`) was already done.

## Fortran ground truth (LAW 2)
- `iel` (thermr.f90:1244-1425): MF6 written with `ltt=6` flag, product TAB1 LIP=`scr(3)=-2`,
  LAW=1; TAB2 LANG=3 (equiprobable cosines), NE=`nne` (calcem grid count, =92 here);
  σ(E) computed as `xie` on the calcem grid `esi` then interpolated onto the MT221 grid
  via `ej(nj)=terp(esi,xie,nne,ex(1),3)` (line 1412) — `il1=3` is the **3-point Lagrange
  ORDER** (linear x/y), NOT a lin-log mode. Last grid point zeroed (line 1413).
  `ncdse = 5 + ne*((nbin+11)/6)` (line 1423) — `ne` = ELASTIC (MT221) grid count.
- dispatch (thermr.f90:404-435): MF7/MT2 present ⇒ `lthr=l1h; icoh=10*lthr`; LTHR=2 ⇒
  icoh=20 ⇒ `iel(20,…,mtref+1)` (line 431). mtref+1 = MT222.
- save-elastic cutoff (thermr.f90:383-395, `up=1.00001` line 161): for the incoherent
  (non-coherent) grid the cutoff is `…, emax (held), up*emax (σ=0)`, then tpend appends
  `etop=2e7` (line 3212). NO `sigfig(emax,7,-1)` point (that is the COHERENT/sigcoh
  grid's last point, only present on the Bragg path).
- tpend ltt=6 (thermr.f90:3293-3301): copies the iel MF6 LIST records verbatim, setting
  only `scr(8)=1` — NO xsi normalization and NO sigfig of the cosines.

## Changes
1. **`src/orchestration/modules/thermr.jl`** `_thermr_sab!`:
   - capture `incoh_endf` from `read_mf7_mt2` (declared `lthr`/`incoh_endf` before the
     `icoh>0` block); LTHR=2 logs the iel branch and leaves `bragg=nothing`.
   - **grid construction split** by path: `bragg !== nothing` → `build_thermal_grid` +
     `_append_emax_sentinels!` (UNCHANGED coherent path, T70-preserving). `bragg === nothing`
     (LTHR=2) → `{el_e ≤ emax} ∪ {emax, up*emax=1.00001, 2e7}` (save-elastic cutoff;
     NO `sigfig(emax,7,-1)`, and `up*emax` not `sigfig(emax,7,+1)`).
   - new iel branch (after the bragg MF3 block): `dwp = incoh_dwp_at(incoh_endf, temp)`;
     `xie,_ = build_incoh_elastic_records(esi_sab, incoh_endf.sigma_b, dwp, nbin, natom)`;
     `iel_xs = [_terp_lagrange(esi_sab, xie, e, 3) for e in thermal_e]`; zero where
     `thermal_e[i] > emax`; `added_mf3[mtref+1]=(thermal_e,iel_xs)`; `push!(thermr_mts,mtref+1)`.
2. **`src/orchestration/pipeline.jl`** `_recompute_thermr_mf6!` (iinc==2, icoh>0):
   - LTHR=2 branch: `xie,iel_records = build_incoh_elastic_records(esi,…)`;
     `ctx.mf6_records[mtref+1]=iel_records`; `ctx.mf6_emax[mtref+1]=emax`; **does NOT set
     `ctx.mf6_xsi`** (iel cosines/weights are unnormalized). Computes the directory-NC quirk
     `ctx.mf6_iel_nc[mtref+1] = 5 + (NP_MF3−1)*((nbin+11)÷6)` (NP read back from the -50 tape).
   - new field `mf6_iel_nc::Dict{Int,Int}` on `RunContext` (+ constructor).
   - `_collect_thermr!` already reads back `extra_mf3[mtref+1]` and pushes `mtref+1`
     when `icoh>0` — confirmed it fires for T74/T69 (card icoh=1).
3. **`src/processing/pendf_writer.jl`** `write_full_pendf` + `_write_mf6_section`:
   - new kwargs `lip::Int=-1` and `sigfig_cosines::Bool=true` on `_write_mf6_section`.
   - iel MTs (keys of `mf6_iel_nc`) emit `lip=-2` and `sigfig_cosines=false` (raw cosines).
   - directory-NC quirk: `section_lines[(6,mt)] = mf6_iel_nc[mt]` overrides the generic
     line-count formula so the MF1/MT451 directory shows the quirk value (1730/1640) while
     the section BODY stays the real 466 lines.
   - new kwarg `mf6_iel_nc::Dict{Int,Int}` threaded write_full_pendf ← final_assembly! ←
     pipeline (ctx.mf6_iel_nc).

## terp helper
`terp(esi,xie,nne,ex(1),3)` is **3-point Lagrange (linear x/y)** — `il1=3` is the order,
not a lin-log INT mode. The existing `_terp_lagrange(xs,ys,arg,il)` (modules/thermr.jl,
direct port of thermr.f90:1427-1538) IS that function; used as `_terp_lagrange(esi_sab,xie,e,3)`.
(The bead's "lin-log mode 3" wording was a misreading of `il1`.)

## MF6 NC quirk override
Threaded a `mf6_iel_nc::Dict{Int,Int}` (MT → quirk NC). Computed in `_recompute_thermr_mf6!`
as `5 + (NP_MF3−1)*div(nbin+11,6)` (the iel `ncdse`, thermr.f90:1423; `ne`=MT221 grid count
=NP_MF3−1 via tpend `scr(6)=ne+1`). The writer uses it ONLY for the MF1/MT451 directory
entry; the emitted MF6 body stays the actual 466 lines.

## Results (fixed-column ENDF parse — reference-test regex is unreliable, Rule 2)
| test | section | lines (J/ref) | NP/dir-NC | diffs | maxreldev |
|---|---|---|---|---|---|
| T74 | MF3/MT221 | 119/119 | 346 / 118 | 0 (byte-identical) | 0 |
| T74 | MF3/MT222 | 119/119 | 346 / 118 | 4 | 3.7e-7 |
| T74 | MF6/MT222 | 466/466 | — / **1730** | **0 (byte-identical)** | 0 |
| T69 | MF3/MT221 | 113/113 | 328 / 112 | 12 | 4.85e-5 (pre-existing inelastic) |
| T69 | MF3/MT222 | 113/113 | 328 / 112 | 5 | 1.69e-7 |
| T69 | MF6/MT222 | 466/466 | — / **1640** | **0 (byte-identical)** | 0 |

- σ(1e-5): T74 81.96615 (exact), T69 6.337533 (exact). First MF3/MT222 data lines and the
  cutoff tails (`…1.000000+0 2.413888+0 1.000010+0 0.000000+0 / 2.000000+7 0.000000+0`)
  byte-match the targets, including the `up*emax=1.00001` cutoff.
- MF6/MT222 first 6 lines byte-match (HEAD L2=1/N1=1; TAB1 LIP=-2 LAW=1 over [1e-5,emax];
  TAB2 LANG=3 LEP=1 NE=92).

## Residual diffs (follow-up bead candidate)
MF3/MT222 has ≤5 interior lines off by ≤4e-7 (7th-sigfig). Root-cause hypothesis: FP-order
grind in `_terp_lagrange` (3-point Lagrange) of `xie` onto the thermal grid — same sub-ULP
class as the pre-existing MT221 inelastic-calcem grind (T70 MF3/MT221 14 diffs ≤9.4e-7).
First point is exact; geometry exact. Not chased (per task scope).

## Regression (each cache-nuked, serial, Rule 9)
- **T70** (coherent, LTHR=1 — refactor-sensitive): MF3/MT222 **0 diffs (byte-identical)**,
  MF6/MT222 stub **0 diffs**, MF3/MT221 14 diffs ≤9.4e-7 — **identical to Phase 87's
  established T70 state** (no regression; coherent path preserved).
- **T09**: **BIT_IDENTICAL 1830/1830**, ALL PASS @ rtol=1e-9.
- **T01**: (see HANDOFF) — NUMERIC_PASS.
- `test_thermr_incoherent_elastic.jl`: **72/72**.

## LTHR=3 (still needed, not in scope)
Dispatch (thermr.f90:432-434, icoh>20): `coh(10,…,mtref+1)` THEN `iel(20,…,mtref+2)` →
coherent MT222 + incoherent **MT223**. Needs: read_mf7_mt2 already returns both bragg+incoh
for LTHR=3; wire the iel branch to fire alongside the coherent branch with `mt_iel=mtref+2`
(not mtref+1); the readback (`_collect_thermr!`) and `_recompute_thermr_mf6!` must register
mtref+2 as well. Driver test T67 (D-in-LiD, natom=2). Tracked under NJOY_jl-69j.
