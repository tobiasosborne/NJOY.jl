# Phase 88b — thermr LTHR=3 mixed elastic (coherent MT222 + incoherent MT223) wired (bead NJOY_jl-yms)

Date: 2026-06-16

## Goal
Generalize the Phase-88 incoherent-elastic (`iel`, LTHR=2) machinery to **LTHR=3 mixed
mode** (coherent + incoherent elastic), driver test **T67** (D-in-7LiD). Dispatch
(thermr.f90:432-434, icoh=30): `coh(10,…,mtref+1)` THEN `iel(20,…,mtref+2)` → coherent
elastic at **MT222** (mtref+1), incoherent elastic at **MT223** (mtref+2), separate sections.

## RED baseline (truncated T67, current Julia BEFORE the fix)
T67 emitted ONLY inelastic, and that grid was STARVED: MF3/MT221 = 103 lines (ref 2357),
MF3/MT222 = MF3/MT223 = 0 (ref 2357 each), MF6/MT222 = MF6/MT223 = 0 (ref 4 / 596).
**Root:** `_thermr_sab!`'s coherent gate was `lthr == 1` only → for LTHR=3 it fell to the
hardcoded-lattice fallback (MAT=15 has no params) → bragg=nothing → the 3374 ENDF Bragg
edges were never merged into the thermal grid (MT221 grid = 299 pts), and no elastic emitted.
So wiring the coherent Bragg for LTHR=3 is the crux — it feeds the merged grid AND emits MT222
AND enables MT223.

## Fortran ground truth (LAW 2)
- dispatch thermr.f90:428-434: icoh=10*lthr; LTHR=3 ⇒ icoh=30 ⇒ coh(10,…,mtref+1) + iel(20,…,mtref+2).
- LTHR=3 directory-NC QUIRK (tpend thermr.f90:3120 `if(mfd.eq.6) nc=ncdse`): the MF6 loop sets
  nc=ncds for MT221, then nc=ncdse for BOTH MT222 and MT223. Because `iel` runs AFTER `coh`
  and overwrites the shared `ncdse = 5+ne*((nbin+11)/6)` (thermr.f90:1423, ne=merged-grid count),
  the coherent MF6 STUB (4 lines) ALSO gets dir NC = the iel value (35300 for T67). For LTHR=1
  (T70) no iel runs → coherent stub keeps NC=3.
- MF6/MT222 coherent stub: LIP = -nbragg = -3374, LAW=0 (coh). MF6/MT223 incoherent: the iel
  LANG=3 equiprobable-cosine records, LIP=-2, raw (un-sigfig'd) cosines (tpend ltt=6).

## Changes (on top of Phase-88 LTHR=2)
1. `src/orchestration/modules/thermr.jl` `_thermr_sab!`: coherent gate `lthr==1` → `(lthr==1||lthr==3)`
   (ENDF lat=10 Bragg fires for LTHR=3 → merged into `thermal_e` → MT222 emitted); iel MF3 gate
   `lthr==2` → `(lthr==2||lthr==3)` with `mt_iel = lthr==3 ? mtref+2 : mtref+1`.
2. `src/orchestration/pipeline.jl` `_collect_thermr!`: also reads back `extra_mf3[mtref+2]` and
   pushes `mtref+2` to `thermr_mts` (LTHR=3); `_recompute_thermr_mf6!`: two independent gates —
   iel records (`lthr==2||lthr==3` → mt_iel) and coherent stub (`lthr==1||lthr==3`); populates
   `mf6_iel_nc` for the iel MT AND (LTHR=3) the coherent stub mtref+1, each `5+(NP_MF3-1)*div(nbin+11,6)`.
3. `src/processing/pendf_writer.jl` `write_full_pendf`: moved the `mf6_iel_nc` directory-NC override
   to a FINAL pass so it also overrides the coherent STUB's NC for LTHR=3. LIP=-2/raw-cosines apply
   ONLY to the `mf6_records` (iel/MT223) path; the coherent stub keeps LIP=-nbragg via
   `_write_mf6_coherent_stub`.

## Results (truncated T67, independently re-verified by the parent — fixed-column diff)
| section | julia | ref | status |
|---|---|---|---|
| MF6/MT222 (coherent stub) | 4 | 4 | **BYTE-IDENTICAL** (LIP=-3374, LAW=0) |
| MF6/MT223 (incoherent) | 596 | 596 | **BYTE-IDENTICAL** (σ(1e-5)=1.026886 exact) |
| MF3/MT221/222/223 | 2351 | 2357 | −6 lines (grid undershoot → km2) |
| MF6 222/223 dir NC | 35215 | 35300 | tracks the undershot grid (quirk mechanism correct) |

The MF6 elastic sections (the deliverable) are byte-identical. The MF3 line-count gap is a
**pre-existing `build_thermal_grid` limitation** on the dense 3374-edge Bragg structure (merged
grid = 7042 pts vs ref 7060; drops input points 1.5, 5.0, and a high-E 5.0…9.6875 sequence;
matches exactly through index 2312) — bead **NJOY_jl-km2**. Also heavy: ~5.1e9 allocations,
1.4 GB peak RSS, 13.5 min wall (a 600s-timeout attempt was killed; true-background run completes).

## Regression (independently re-verified; serial Julia, cache-nuked)
- **T70 LTHR=1**: MF3/MT221 498 ln (14 diffs ≤9.4e-7 = Phase-87), MF3/MT222 + MF6/MT222 stub
  BYTE-IDENTICAL, **MF6/MT222 dir NC=3 (NOT 35300) — the LTHR=3 quirk did NOT leak into LTHR=1.**
- **T74 LTHR=2**: MF3/MT221 byte-identical, MF6/MT222 466 ln byte-identical, dir NC=1730 (agent + Phase-88).
- **T09**: BIT_IDENTICAL 1830/1830. Unit test `test_thermr_incoherent_elastic.jl` 72/72.

## Follow-up
- `NJOY_jl-km2` (P2): `build_thermal_grid` under-includes input points + O(n²)-ish cost on dense
  Bragg → full T67 structural pass (MF3 2357 lines, MF6 dir NC 35300) + reasonable runtime.
- `NJOY_jl-c9f` (P3): MF3/MT222 (LTHR=2) ≤4e-7 3-pt-Lagrange FP grind.
