# Phase 71 — Covcal MT=102 row-1 diagnosis: rescon (MF=32 → MF=33) is the missing piece

**Date:** 2026-05-07
**Test:** T15 (U-238 JENDL, MAT=9237, second errorr call → tape26)
**Goal:** Diagnose T15 tape26 MT=102 self-cov line drift (Julia 41 vs ref 79 lines, –38) and identify the smallest tractable canary for the next session's port.

## Outcome

**Diagnosis only — no production code changes.** Discovered that the
HANDOFF P1 "Bug-B siblings (LB=1/LB=2 midpoint sampling)" hypothesis is
**not** the root cause of the T15 covcal content drift.

The dominant cause is missing resonance-parameter uncertainty
propagation: Fortran's `covout` (errorr.f90:7465) calls `rescon`
(errorr.f90:8513-8819) to add MF=32 RP-cov contributions to MF=33
covariances for seven specific (mt, mt2) pairs:

| itp | (mt, mt2) | Branch        |
|-----|-----------|---------------|
| 1   | (18, 18)  | fission/fission |
| 2   | (18, 102) | fission/capture |
| 3   | (102, 102)| capture/capture |
| 4   | (2, 2)    | elastic/elastic |
| 5   | (2, 18)   | elastic/fission |
| 6   | (2, 102)  | elastic/capture |
| 7   | (1, 1)    | total/total |

U-238 JENDL has an 8092-line MF=32 section (single MT=151 RP
covariance, RML or SAMMY). Julia has zero MF=32 → MF=33 propagation
today. This explains the per-MT line drift pattern in T15 tape26
(Phase 56 baseline):

| MT  | Julia | Ref  | Δ    | rescon contributes? |
|-----|-------|------|------|---------------------|
| 1   | 171   | 271  | -100 | (1, 1) yes          |
| 2   | 1506  | 1400 | +106 | (2, 2)/(2, 18)/(2, 102) yes — but Julia *over*, separate issue (NC) |
| 18  | 266   | 228  | +38  | (18, 18)/(18, 102) yes — Julia over, also NC-related |
| 102 | 41    | 79   | -38  | (102, 102) yes — purely missing rescon |
| 4   | 1106  | 1106 | 0    | no — but cross-pairs drift from NC |

MT=102 is the cleanest single-cause canary: the entire –38 line drift
is rescon (102, 102). MT=1's –100 lines is rescon (1, 1). MT=2 / MT=18
*overflow* (Julia over ref) is NC-expansion-related, orthogonal.

## Evidence

### Matrix-level dump for (102, 102) self-cov — the canary

Extracting Julia's MT=102 30×30 cov matrix from `tape26` and reference
`referenceTape26` shows a sharp boundary at LANL group 15 (energy
~67380..183156 eV — straddles the 150 keV input-bin-1 boundary):

```
ig | jul cols (count) | ref cols (count)
 1 | jul: —           | ref: 1..10 (10)   ← all zero in jul
 2 | jul: —           | ref: 1..10 (10)   ← all zero in jul
...
14 | jul: —           | ref: 13..15 (3)   ← all zero in jul
15 | jul: 15..18 (4)  | ref: 13..18 (6)   ← jul missing cols 13..14
16 | jul: 15..18 (4)  | ref: 15..18 (4)   ← exact match
17 | jul: 15..18 (4)  | ref: 15..18 (4)   ← exact match
...
30 | jul: 24..30 (7)  | ref: 24..30 (7)   ← exact match
```

**Rows 16..30 match exactly.** This proves the LB=5 weighted-collapse
path (Phase 51) IS correct above the resonance region. **Rows 1..14
are entirely missing in Julia.** Row 15 is partial (Julia has the
high-E half, missing the low-E sub-bin contributions).

The reference values in row 1 columns 1..10:
```
2.658914e-4  3.299130e-4  3.496066e-4  4.091654e-4  2.493738e-4
9.086889e-5  1.684544e-5  8.319187e-8 -4.685423e-6 -3.113298e-7
```

These are pure resonance-parameter contributions:
- **Negative values** in cols 9..10 cannot come from MF=33 LB=5 (whose
  fvals are all non-negative for U-238 capture); they're sandwich-rule
  outputs from RP cov × sensitivity matrices.
- **Diagonal magnitude ~3e-4** corresponds to ~1.6% relative SD —
  consistent with U-238 thermal capture xs RP uncertainty propagation.

### Why rows 1..14 are zero in Julia

The input MF=33/MT=102 has a single LB=5 LT=1 block, NE=22 energies
starting at E_1 = 1e-5 eV, with the entire first row of fvals zero
(`fvals[1, *] = 0.0`):

```
input MF=33/MT=102 fvals matrix (21×21 upper triangular):
  fvals(1, *) = 0  ← bin 1 = [1e-5, 150 keV] sub-threshold
  fvals(2, 2) = 3.6e-4
  fvals(3, 3) = 5.5e-4
  ...
  fvals(13, 13) = 0.249
  fvals(14, 14) = 0.063
  ...
```

LANL-30 groups 1..14 fall entirely below 150 keV (input bin 1). Both
Julia (`expand_lb5_symmetric` + xs/flx-weighted union-grid collapse,
Phase 51) and Fortran's covcal (errorr.f90:2208-2235) correctly
compute zero LB=5 contribution for any (jg, jh) with both endpoints in
input bin 1. **The missing piece is rescon's additive RP-cov
contribution at line 7465** — Fortran's covcal-only output for MT=102
would also be zero in rows 1..14; rescon is what fills them in.

### Where in Fortran

```
errorr.f90:7465   (covout)   call rescon(ix, ixp, csig, cova, ...)  if mf32 != 0
errorr.f90:8513   (rescon entry) — 306 lines
errorr.f90:3011   (resprx)   resonance-parameter cov reader (RML/LC0/LC1/LC2)
errorr.f90:3252   (rpxsamm)  SAMMY-format MF=32 reader
errorr.f90:3734   (rpxlc0)   LC0 (compatible) reader
errorr.f90:4108   (rpxlc12)  LC1/LC2 reader
errorr.f90:4634   (rpxlc2)   LC2 (long-range) cov
errorr.f90:5411   (grpav4)   group-averaging for MF=4 + RP sensitivity build
```

The data flow:
- `resprx` reads MF=32 → builds RP covariance + parameter covariances
- `rpxlc0`/`rpxlc12`/`rpxsamm` per format
- Sensitivity computation: dσ_X/dRP at each group energy (this is the
  big one — needs MLBW/RM/SAMMY evaluation + analytic derivatives)
- `rescon` does the final sandwich: cova[ig, ig2] += J_X[ig, :] · Cov_RP · J_X[ig2, :]'
  via precomputed `cff`, `cfc`, `cee`, `cec`, `ccc`, `cef`, `cfx` arrays

This is a multi-day port (~1500-2500 LOC of Fortran across 8+
subroutines). NOT a single-session effort.

## Why the HANDOFF was wrong

HANDOFF P1 said:
> Bug B siblings: `expand_lb1` / `expand_lb2` likely have the
> midpoint-sampling class of bug Bug B fixed for LB=5 (smaller blast
> radius — diagonal only).

But MT=102's input is a single **LB=5** block (not LB=1/2), and that
block already routes through Phase 51's σ·flx-weighted union-grid
collapse via `_collapse_pair_blocks!` (verified `have_xs = true` for
mt=102). The matrix dump confirms the LB=5 path produces exactly the
right values for rows 16..30. The drift is pure rescon-missing.

Lesson: the "Bug-B siblings" framing was a hypothesis from
LB-extrapolation; the matrix-level diff against reference data was
the diagnostic that broke the wrong frame. Per CLAUDE.md Rule 1
(skepticism) and Rule 0 (oracle-driven): always pull the actual
reference matrix into a Julia REPL and compare element-wise before
trusting a HANDOFF root-cause label.

## What's worth landing this session

1. **Worklog** (this file) — captures the diagnosis so the next
   session doesn't repeat the investigation.
2. **HANDOFF update** — P1 Covcal section: rescon port is the
   actual remaining work, not Bug-B siblings.
3. **RED test** for the MT=102 row-1 canary so progress can be
   measured.
4. **Scaffold** for the future port: MF=32 detection + loud-TODO at
   the would-be rescon call site, returning early when no MF=32 is
   present (so existing tests stay green). Keeps the door open for
   incremental landing of the rescon branches one at a time.

## Acceptance for the eventual port (multi-session)

- [ ] T15 tape26 MT=102 lines: 41 → 79 (matches ref).
- [ ] T15 tape26 MT=1 lines: 171 → 271.
- [ ] T15 MT=102 self-cov C[1, 1] ≈ 2.658914e-4 (within 1e-7).
- [ ] T15 MT=102 self-cov C[1, 9] ≈ -4.685423e-6 (negative — proves RP
      sandwich, not MF=33 LB=5).
- [ ] T15 tape26 total: 5964 → 5958 (or closer to ref).
- [ ] T22 BIT_IDENTICAL preserved.
- [ ] T04 NUMERIC_PASS preserved (tape23 81/82).

## Files

NEW:
- `worklog/phase71_rescon_diagnosis.md` (this file)

MODIFIED (this phase):
- `HANDOFF.md` — P1 Covcal content drift section updated to point at
  rescon as the missing piece.
- `test/validation/test_errorr_covcal_lb5.jl` — added new RED testset
  asserting MT=102 row-1 canary values; expected to fail until rescon
  lands.
- `src/orchestration/modules/errorr.jl` — MF=32 detection + loud-TODO
  scaffold at the would-be rescon call site (no functional change;
  early-returns when MF=32 absent).

## Fortran source citations

- `njoy-reference/src/errorr.f90:7465` — `call rescon(ix, ixp, csig, cova, ...)`
- `njoy-reference/src/errorr.f90:8513-8819` — `rescon` subroutine
- `njoy-reference/src/errorr.f90:3011-3250` — `resprx` (MF=32 reader)
- `njoy-reference/src/errorr.f90:3252-3732` — `rpxsamm` (SAMMY format)
- `njoy-reference/src/errorr.f90:3734-4106` — `rpxlc0` (LC0 format)
