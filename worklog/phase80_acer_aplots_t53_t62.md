# Phase 80 — ACER aplots: heating + threshold + aploxp (T53 + T62 FULL BIT_IDENTICAL)

**Date:** 2026-05-31
**Bead:** NJOY_jl-zsh (ACER aplots: port heating/threshold/particle-production blocks)
**Outcome:** **2 NEW FULL BIT_IDENTICAL: T62, T53.** Continues Phase 79 (which made
T50/T52 tape33 BIT_IDENTICAL). All edits in `src/formats/ace_aplots.jl` (+587/−89).

Orchestrated session (Opus coders, parallel Sonnet researchers for read-only scoping;
serial Julia; verified-before-trust per Rule 2). The Phase-79 scaffold had already placed
`_aplots_unported(...)` guards at the **correct positions** in the page-emission sequence,
so this work was replacing guards with byte-faithful emission code at existing slots — not
restructuring the page order.

## What landed

Ported the four remaining tape33 blocks that T53/T62 hit (Fortran `acefc.f90`):

1. **aplotr threshold-reaction page** (`acefc.f90:16970-17093`) — T62's sole missing block,
   also T53 page 5. Replaced the guard in `_aplotr`. One curve per qualifying threshold
   reaction (MT filter: excludes mt<5, fission 18-21/38, inelastic 50-90 for izai=1,
   207<mt<600, sub-threshold xss(iaa)<1e-6); `1 0 2 1/` literal axes (NO xtag/ytag, unlike
   the inelastic page); `_ascle(4,...)` linear axes; thinned `1p,2e14.6` data; the
   `0 0 0 N/` color records (iwcol=3) and 10-char MT name labels.

2. **log-log heating page** (`acefc.f90:16191-16226`) — T53 page 2. ESZ heating column
   `xss[esz+4*nes-1+i]`, floor below `hmin=1e-10`; fires if `ymax>ymin`; `ymin<ymax/scale`
   floor (scale=1e6); both axes `_ascll` (log); `4 0 2 1/`; log thinning.

3. **lin-lin heating page** (`acefc.f90:16557-16595`) — T53 page 4. Fires if heating
   nonzero; `_ascle(4,...)` both axes; `1 0 2 1/`; energy threshold e≥0.2 MeV; linear
   thinning; always-write-last (i==nes); no hmin floor.

4. **aploxp particle-production** (new `_aploxp`, `acefc.f90:18787-19460`) — T53 pages 7-13,
   gated on `nxs[NXS_NTYPE]>0`. Sub-blocks:
   - C1 particle heating contributions (`18816-18934`): per-type curves, heating
     `xss[hpd+1+naa+i]`, `1 0 2 1`+xtag/ytag (`_e12_4_no1p`), y-label `<m>e<v>/collision`,
     inter-curve separator `2/ / 0 0 0 N/`.
   - C2 recoil heating (`18936-19004`): recoil = esz_heat − Σ particle heating; single curve;
     `ymin<0 && ymax<-ymin/2` clamp.
   - C3 production XS (`19006-19119`): always fires; xs `xss[hpd+1+i]`; photons/5 with label
     `photons/5`; xtag `test=30` branch.
   - C5 3D angular pages (`19121-19458`): per type/MT, fires when `na=xss[landh+imt-1]>0`.
     **Reuses the existing `_aplof4` perspective machinery** — the C5 angular page is
     byte-identical to `_aplof4`'s `k<0` tabulated path apart from subtitle and
     `andb`→`andh`; the perspective body was factored into a shared `_aplof4_perspective!`.

## Fortran surprises / subtleties pinned

- **`1p,e14.6` drops the `E` for 3-digit exponents.** `1.155100e-117` prints `  1.155100-117`
  (mantissa, sign, 3 digits, no `E`), still width 14. C `%14.6E` keeps the `E`
  (` 1.155100E-117`). Fixed `_e14_6` to strip `E` only when the exponent has ≥3 digits;
  2-digit exponents (the form T50/T52's BI pages rely on) are untouched. This case was
  **never exercised** by the principal-XS pages T50/T52 made BI in Phase 79 — it first
  appeared in the threshold-page sub-threshold underflow values.
- **C5 angular subtitle** is `angular distribution for <len_trim(name)> <particle-word>`
  (name trailing-blanks STRIPPED), e.g. `(d,p*0) triton`. This differs from `_aplof4`'s
  elastic subtitle which keeps the full 10-char padded name (`elastic   `). Particle words
  (`19387-19410`): 1→neutron, 9→proton, 31→deuteron, 32→triton, 33→`3he`, 34→alpha.
- **Incident-particle name swap** `name(2:2)='d'` for izai=1002 (`19312-19326`) is caller-side
  (a helper `_aplots_name_izai`), not in `mtname` — so `(n,p*0)` → `(d,p*0)`.
- **aploxp locator offsets** (per type i, from `19124-19133`): `hpd=ploct+10*(i-1)`,
  `landh=+5`, `andh=+6`, `lsigh=+3`, `sigh=+4`; `ipt=ptype+i-1`; energy `e=xss[esz+iaa-2+i]`.
- The **mtname table** was extended from a 50-entry stub to the full 500 entries
  (transcribed from `acecm.f90:31-138`) so MT≥600 charged partials resolve (MT=600→(d,p*0)).
- **law 4/44/61 emission-spectrum sub-block** (aploxp, `19144-19305`) is NOT ported — guarded
  with a loud `_aplots_unported` (Rule 6); no T50-54/62 reaction triggers it.

## Verification (cache-nuked, serial)

- **T62**: tape33 BIT_IDENTICAL 5460/5460 (`cmp`-clean), tape34 7221/7221, tape35 — FULL BI.
- **T53**: tape33 BIT_IDENTICAL 17236/17236 (`cmp`-clean), tape34 12030/12030, tape35 — FULL BI.
- **Regressions — none**: T50 tape33 432/432, T52 tape33 6042/6042, T61 tape71 49211/49211 —
  all still BIT_IDENTICAL.

## Other aplots-using tests (no regression, blocked elsewhere)

- **T54** (p+H-3): tape33 now structurally complete **11340/11340 lines** (was empty stub);
  DIFFS 9917/11340 trace 1-ULP to T54's tape34 ACE values (e.g. total XS `1.427856E+08` vs ref
  `1.427855E+08`). aplots port is correct — blocked by the pre-existing tape34 FP residual
  (bead **NJOY_jl-53h**, the triton recoil-heating word). Flips to FULL BI once that lands.
- **T51** (1002.10h) and **T71** (78184.10o): aplots generalizes (no crash), but tape33 is
  wrong-content because their **ACE (tape34) is structurally wrong** (T51 unported MF6 LAW=6,
  T71 ACE over-produces 72594 vs 3052 lines). Pre-existing, ACE-level, separate from aplots.

## Net

ACER FULL BIT_IDENTICAL: T50, T52, T61 → **+T53, +T62**. The aplots-porting scope of bead
NJOY_jl-zsh is complete for every test whose ACE is byte-identical; the remaining aplots
sub-paths (MT=444 damage, URR, aplonu nubar, aplodn delayed-neutron, aplopf higher-fission,
aplof4 equiprobable contours, aploxp law-4/44/61 emission) stay loud-guarded and ungated by
any current test.
