# T50 — ACER module promotion, Phases 1-5

**Date**: 2026-04-24 (Phases 1-3) + 2026-04-27 (Phase 4 + Phase 5)
**Test**: T50 (α-particle + He-4, MAT=228, ENDF/B-VIII.0, charged-particle ACE)
**Goal**: Promote `acer_module` from iopt=1 MF3-only stub to real ACE generator targeting T50's 163-line `referenceTape34`.

## Score at session close

| Artefact | Before | Phase 3 | Phase 4 | **Phase 5** | Ref | Status |
|----------|--------|---------|---------|-------------|-----|--------|
| tape34 exact lines | — | 11 / 163 | 16 / 163 | **35 / 59 produced** | 163 | ESZ values bit-identical (37/37) |
| tape34 total lines | — | 49 | 53 | **59** | 163 | 104 lines short (angular block) |
| tape34 NES | — | 29 | 32 | **37** | 37 | EXACT match ✓ |
| T01 tape25 (regression) | 32750/32962 | 32750/32962 | 32750/32962 | **32812/32962** | 32962 | NUMERIC_PASS @1e-5; unchanged from Phase 4 baseline ✓ |
| T02 tape28 (regression) | 12037/13873 | 12037/13873 | 12037/13873 | **12519/13873** | 13873 | NUMERIC_PASS @1e-5; unchanged from Phase 4 baseline ✓ |
| T50 run | crash | RAN_OK | RAN_OK | RAN_OK | runs | ✓ |

(T01/T02 numbers in earlier rows were stale at the time of writing; the
Phase-4 baseline is what Phase 5 measures and confirms unchanged.)

Six commits pushed: `bd43751` (P1), `fe5aeec` (P2), `c061919` (P3),
`820d64e` (HANDOFF + worklog P1-3), `b952b33` (P4), `e543e56` (P5).

## Phases

### Phase 1 — Parser + iopt=7 passthrough

Commit `bd43751`.

**Bug**: `AcerParams` card-1 parser read only 4 slots (`nendf npendf nace ndir`) but
Fortran `acer.f90` has 5 (`nendf npend ngend nace ndir`). This shifted EVERY acer
test that has a nonzero ngend — T14, T48, T50-54, T62, T71. T50's `-21 -22 0 31 32`
was being parsed as `nendf=21, npendf=22, nace=0, ndir=31`, dropping the ACE
output unit entirely.

**Fix**: `AcerParams` now has 13 fields including `ngend`, `iprint`, `itype`,
`title`, `nplot`. Card 2 position 4 is parsed as suffix string (`.XX` →
two-digit tag) or `nplot` for iopt=7.

**iopt=7 module**: `_acer_iopt7` copies the input ACE tape (unit `npend`)
verbatim to the output unit (`nace`), touches `ngend` as a stub viewr plot
tape, and writes a 1-line xsdir summary to `ndir`. Matches Fortran iopt=7
check-mode behaviour for un-thinned input.

### Phase 2 — Header byte-identical modulo date

Commit `fe5aeec`.

Six fixes bringing T50 header from 0/6 to 6/6 lines bit-identical (execute.py
masks the date):

1. `ACEHeader` constructor: `hz` and `hd` use `lpad` not `rpad`. Fortran's
   `a10` RIGHT-justifies strings shorter than the field width. Ref
   `"  2004.10a"` (2 leading spaces); Julia was writing `"2004.10a  "`.

2. `mat_string` no longer `strip()`'d in `ACEHeader`. The three leading
   spaces in `"   mat 228"` ARE the alignment and must survive into the
   10-char `hm` buffer. Ref line 2 cols 71-80; Julia was flipping to
   `"mat 228   "`.

3. IZ/AW pair format: Fortran `f11.0` emits `"         0."` (trailing
   decimal point); C/Julia `%11.0f` drops the `.`. Changed writer
   format to `"%7d%10.0f."`.

4. `parse_acer` strips surrounding quotes from the title token. NJOY
   titles are `'…'`-delimited in the input deck; the tokenizer preserved
   the quotes, which ended up inside the ACE comment field.

5. `build_ace_from_pendf` gains `za=`, `awr=`, `date=` kwargs. Julia
   was defaulting `awr = Float64(za % 1000)` (A, mass number, not AWR).
   For He-4: ref `3.968219` vs trial `4.0` — 0.8% mass error.

6. New `read_mf1_mt451_header` in `src/endf/readers.jl` parses the three
   MF1/MT451 CONT records and returns ZA, AWR, LRP, NLIB, ELIS, LIS,
   LISO, NFOR, AWI, EMAX, LREL, **NSUB**, NVER. NSUB is the incident-
   particle indicator (ENDF-6 §1.1):

   | NSUB | Incident | ACE letter |
   |------|----------|-----------|
   | 10    | neutron   | c |
   | 0     | photon    | p |
   | 10010 | proton    | h |
   | 10020 | deuteron  | o |
   | 10030 | triton    | r |
   | 20030 | He-3      | s |
   | 20040 | **alpha** | **a** |
   | 113   | elec-nuc  | n |
   | 3     | electron  | e |

   New `acer_incident_letter(nsub)` maps NSUB → letter. `acer_module` reads
   MF1/MT451 from the ENDF (nendf, not npendf — matching Fortran
   `acefc.f90 first()` at lines 280-360), derives the letter, and replaces
   the parser's provisional letter with the NSUB-driven one. T50 →
   ZAID `2004.10a`.

### Phase 3 — JXS layout + total-XS reconstruction

Commit `c061919`.

**JXS null-pointer convention**: Fortran `aceout` (acefc.f90:5367-5371)
treats MTR/LQR/TYR/LSIG/SIG as "structurally present, possibly zero-
length" — their JXS pointers always equal `end-of-ESZ+1` even when
NTR=0. Julia was leaving JXS[3-7]=0 when `n_react==0`. Ref shows
JXS[3..8] all = 186 for T50 (end of 185-word ESZ).

Similarly LDLW/DLW use the end-of-XSS sentinel, not 0, when NR=0.
JXS[22] (END) = `length(xss)` (last valid word), not `length+1`.

**Total-XS reconstruction**: `build_ace_from_pendf` previously wrote
zeros for the ESZ total XS column when MT=1 was absent. Charged-
particle evaluations commonly have only MT=2 (α+He-4 has literally
MT=2 and nothing else). Now reconstructs total as sum of non-redundant
partials, matching Fortran `unionx`.

### Phase 4 — MF6 incident-energy union

Commit `b952b33`.

**Bug**: ESZ grid was the raw 29-point MF3 grid. Fortran's `unionx`
guarantees every MF6 tabulation point lands on an ESZ grid point —
otherwise the ACE angular/energy block pointers (LAND/AND, LDLW/DLW)
have nowhere to cross-reference at lookup time.

**New reader** `read_mf6_incident_energies(filename, mat, mt)` in
`src/endf/readers.jl`. Walks first MF6/MT subsection: HEAD →
yield TAB1 (whose L2 carries LAW) → if LAW ∈ {1,5,6,7} read TAB2 +
NE record CONT headers, extracting C2 (incident E) from each. NJOY
standard assumes all NK subsections share the same incident grid
(verified by `topfil` at acefc.f90:3500+).

Subtle bug fixed during port: ENDF MF6 subsections start with the
**yield TAB1's HEAD record** — there is NOT a separate sub-section
CONT before the yield TAB1. Initial implementation read the sub-head
as a CONT then called `_discard_tab1` which tried to read another
HEAD that wasn't there. Fix: `read_tab1` consumes the whole yield
including its HEAD; LAW comes from `yield.L2`.

For T50 MF6/MT=2 LAW=5: returns `[0.295, 3.0, 4.0, 12.9, 16.6, 20.0]`
MeV. Three of these (0.295, 3.0, 20.0) are already in MF3; three are
new (4.0, 12.9, 16.6).

**Wiring**: `acer_module` iterates over every MT in the PENDF MF3
dict, attempts `read_mf6_incident_energies`, and unions the result
into `master_e`. Existing per-MT linear interpolation code re-samples
each MT's XS onto the unioned grid.

**Impact**: T50 NES 29 → 32; XSS length 146 → 161; tape34 exact
11 → 16. T01 / T02 unchanged (verified — neither uses iopt=1 + MF6).

### Phase 5 — `unionx` step-1.2 + dedup → ESZ exact

Commit `e543e56`.

**Pre-flight diagnostics**: patched `unionx` in
`njoy-reference/src/acefc.f90` with `write(*,...)` between MF6-union
and step-pass, and inside `aceout`'s ESZ loop. Re-compiled and re-ran
Fortran on T50 to capture the exact intermediate grids:

  | Stage | Points |
  |-------|--------|
  | MF3 (post-stage-1, pre-stage-2) | **32** (29 MF3 + 3 unique MF6 anchors, but with quirks below) |
  | Post-step-pass | **39** |
  | Final ESZ via `gety1` | **37** |

The Phase-4 worklog's claim that the Fortran does "thin(3)=0.001
midpoint-on-line drop" was wrong — `unionx`'s `thin` parameter array
is all-zero by default (T50 deck has no card-7 thinning) and ithopt=0
disables thinning entirely. The actual 3-point drop comes from a
**c-buffer overwrite quirk** in the MF6-union loop, not adaptive
linearization.

**Three Fortran quirks reproduced verbatim**:

1. **MF6-union c-overwrite drop** (acefc.f90:1608-1656). The inner
   `do while (eg.lt.ee)` writes `c` from the *previous* `finda`, then
   advances `k` and re-reads `c`. After the inner loop exits, the
   explicit `c(1)=ee` overwrites `c` before the next `iee` iteration's
   first `loada` runs — so the MF3 point that was just read into `c`
   is lost, AND the next iter's first `loada` writes the anchor again,
   producing a duplicate. Effect on T50: drops MF3 4.4416 (between
   MF6 anchor 4.0 and next MF3 4.7004), drops 13.578 (between 12.9
   and 14.416), drops 17.906 (between 16.6 and 20.0). Adds duplicates
   at 4.0, 12.9, 16.6.

2. **Step-1.2 off-by-one drop of second-to-last** (acefc.f90:1658-1686).
   The outer `do while (k.lt.nold)` advances k inside the body; the
   inner `if (k.lt.nold)` skips processing for the last k. Net effect:
   processes panels [old[1],old[2]]…[old[nold-2],old[nold-1]] (writing
   each left-endpoint + step pads), then the post-loop `loada(jt=-j)`
   writes old[nold] — old[nold-1] is **never written**. For T50 this
   silently drops the second 16.6 duplicate, so post-step has only
   one 16.6 instead of two. Step pads (8 of them) inserted at 0.72,
   0.88, 1.30, 1.80, 2.40, 2.90, 3.60, 8.40 MeV.

3. **gety1 dedup at aceout ESZ readback** (acefc.f90:5341-5355).
   When MT=1 is read back from the scratch tape via `gety1` to
   populate the ACE ESZ block, adjacent equal energies (TAB1 step
   discontinuities) collapse to single entries. Removes the remaining
   duplicates at 4.0 and 12.9 → final 37 ESZ points.

**Julia port**: new `_acer_unionx_charged(mf3_grid, mf6_anchors)` in
`src/orchestration/modules/acer.jl` reproduces all three stages with
matching FP semantics (`round_sigfig(1.2*eg, 2, 0)` mirrors
`sigfig(step*eg, 2, 0)`). Gated on `MF1/MT451 NSUB ≠ 10` so the
neutron path (T01/T02 etc.) is untouched.

**Impact**: T50 NES 32 → **37** (matches reference exactly). All 37
ESZ energies bit-identical to reference. tape34 exact lines 16 → **35
of 59 produced** (35/163 of reference). T01 32812/32962 unchanged
(NUMERIC_PASS @1e-5). T02 tape28 12519/13873 unchanged.

## Outstanding work

### LAND/AND elastic angular block (418 XSS words)

Completely missing. Ref structure:
- LAND[0] = 1 at XSS[186] (offset into AND block)
- AND[1] = 6 (NE = number of incident energies with angular dist)
- AND[2..7] = 6 incident energies (0.295, 3, 4, 12.9, 16.6, 20 MeV)
- AND[8..13] = 6 LC locators (all negative → tabulated CDF form)
- AND[14..604] = 6 tabulated (μ, pdf, cdf) tables

Source: MF6/MT=2 with LAW=5 (identical particles, Coulomb+nuclear).
Conversion from ENDF MF6 LAW=5 → ACE LAW=12/14 (tabulated CDF)
is the core work.

**Next step**: port Fortran `acelod` LAW handling (acefc.f90:4890-6320),
specifically the MF6 LAW=5 → ACE conversion path. This is the main
angular-distribution converter and will be reusable across T51-T54,
T62, T71.

### Known-untested bugs

- `acer_incident_letter` NSUB mappings for deuteron/triton/He-3 are
  speculative — verify against T51 (proton+H-2, letter `h`? or `o`?),
  T53 (deuteron+H-2), T54 (proton+H-3), T62 (deuteron+He-3).
- Suffix numeric portion always 2 digits — Fortran MCNPX can use 3
  digits (`.001a`). Not exercised by current tests.
- iopt=7 nplot=-1 vs nplot=other not differentiated. Fortran uses
  nplot to drive aplots; Julia always writes an empty viewr tape.
- `read_mf6_incident_energies` only walks the FIRST subsection of NK.
  All subsequent subsections are skipped with no validation that they
  share the same incident grid. Fortran `topfil` does verify this for
  ENDF-6 evaluations; if a test ever surfaces with mismatched grids
  the reader will need extension.
- LAW=2 path in the MF6 reader is a no-op. Standard MF6 elastic
  evaluations don't use LAW=2; angular-only distributions do. Port
  when a test surfaces it.

## References consulted

- Fortran source: `njoy-reference/src/acer.f90` (dispatcher), `acefc.f90`
  (iopt=1 `acetop`, `unionx`, `acelod`, `aceout`), `acecm.f90` (XSS
  format helpers).
- ENDF-6 manual §1.1 (MF1/MT451 format + NSUB encoding), §6 (MF6
  LAW=5 identical-particle elastic).
- T50 input deck: `njoy-reference/tests/50/input`.
- T50 reference: `njoy-reference/tests/50/referenceTape34` (163 lines,
  XSS=604 words).
