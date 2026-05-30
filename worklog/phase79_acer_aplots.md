# Phase 79 — ACER aplots: viewr plot tape (tape33) (2026-05-30)

Orchestrated, subagent-driven session. Three parallel read-only Sonnet researchers
(Fortran aplots spec / reference-tape byte-format reverse-engineering / Julia
integration + gap analysis) synthesised their reports into a single spec; one Opus
coder then worked serial (Rule 9, cache-cleared before every Julia run), targeting
T50 first (smallest reference at 432 lines) before widening to T52. Orchestrator
re-verified all outcomes independently before trusting the report (Rule 2).
Commit `90b1754` (already pushed).

## Headline

**T50 and T52 `tape33` BIT_IDENTICAL, completing both tests' full ALL-PASS
@1e-9 status.** `tape33` was the *sole remaining failing tape* for T50/T52/T53/T62
after the Phase-78 run; porting `aplots` flips T50 and T52 immediately. T61 (thermal
iopt=7) retains its BIT_IDENTICAL status — the class-letter gate correctly emits no
plot for class `'t'`. T53 and T62 preserve their existing tape34/tape35 BIT_IDENTICAL
results with no crash, and tape33 is deferred gracefully via the
`AplotsNotPortedError` fallback path.

Running BIT_IDENTICAL count: T03, T22, T50, T52, T61, T53, T62 (T50/T52 now fully
ALL PASS @1e-9).

## Why aplots was the chosen next step

The fresh post-Phase-78 sweep showed T50/T52/T53/T62 each had tape34 (ACE) and
tape35 (xsdir) BIT_IDENTICAL. `tape33` (ngend, the viewr plot tape) was the only
remaining failing tape for all four. Because tape34 is byte-identical, the ACE
arrays that `aplots` consumes are provably correct; `aplots` is therefore a pure
function of verified data with zero upstream risk. T50 (bare H-1 elastic, 432-line
reference) was the smallest and cleanest possible bite.

## What landed (commit `90b1754`, 1178 insertions)

### `src/formats/ace_reader.jl` (116 lines, new)

`read_ace_ascii(lines) -> ACETable` — inverts the Type-1 ASCII ACE writer
(Fortran `acefix` itype=1 read path, acefc.f90:13992-14013). Reads the fixed-column
header (a10/f12.6/e11.4/a10, then a70/a10), the 16 (IZ,AW) pairs at `4(i7,f11.0)`,
the 48 NXS+JXS integers at `8i9`, and the XSS body at `4e20.0`. Targets the
existing `ACETable` struct and named JXS/NXS constants in `ace_types.jl`. Because
the file being read is the same byte-identical tape34 just written, the reader
functions as a verified inverse of an already-verified serialisation. Helper
`_parse_ace_real` handles Fortran's drop of the `E` character for 3-digit exponents
(e.g. `1.15510000000-117`), which `Julia parse` rejects without it.

### `src/formats/ace_aplots.jl` (1031 lines, new)

Port of Fortran `aplots` (acefc.f90:15212-16830) and its delegates:

- **`_acer_aplots`** — main entry: global header card, log-log principal XS page,
  resonance/URR/heating/damage guards (raise `AplotsNotPortedError` for any
  not-yet-ported block), lin-lin principal XS page, heating/damage/non-threshold
  lin-log guards, delegated `_aplotr`/`_aplopf`/`_aplof4`, `aplonu`/`aploxp`/
  `aplodn` guards, and the `99/` + `stop` terminator.
- **`_aplotr`** (acefc.f90:16832-17094) — inelastic levels + threshold reactions.
  Both inner loops are fully ported; for `ntr=0` they emit nothing. For
  not-yet-ported cases (inelastic-level pages for izai=1, threshold-reaction pages)
  raises `AplotsNotPortedError`.
- **`_aplopf`** (acefc.f90:17096-17203) — higher fission reactions (MT=20/21/38).
  Loop ported; emits nothing when `ntr=0`.
- **`_aplof4`** (acefc.f90:17205-17428) — angular-distribution 3D perspective
  pages. Full tabulated (`k<0`) path ported; equi-probable-contour (`k>=0`) path
  raises `AplotsNotPortedError`. This is the page that fires for T50's elastic
  angular distribution (`-1 2/` header, 3D perspective).
- **`_ascll` / `_ascle`** (acefc.f90:19684-19724, 19541-19682) — log and linear
  axis scaling, ported to match Fortran's exact branching and rounding.
- **Three Fortran float formatters**: `_e14_6` (`1p,e14.6`, data fields),
  `_e12_3` (`1p,e12.3`, axis limit fields), `_e12_4_no1p` (scale/legend-tag coords
  — see Fortran Surprises below).
- **`_aplots_mtname`** (acecm.f90:21-224) — MT → 10-char reaction name with
  charged-particle alternates, verbatim trailing-blank padding for byte identity.
- **`AplotsNotPortedError`** — dedicated exception type; `_acer_iopt7` catches it
  specifically and falls back to an empty stub, leaving all other output tapes
  intact. Any other exception propagates and aborts loudly (Rule 6).
- **`_aplots_nonthreshold!`** and **`_aplots_nonthreshold_linlog!`** — log-log and
  lin-log non-threshold reaction pages (acefc.f90:16290-16398, 16670-16793). For
  `ntr=0` the reaction loop never iterates and nothing is emitted.

### `src/orchestration/modules/acer.jl` (30-line diff)

`_acer_iopt7` stub (empty `touch`) replaced with a real `_acer_aplots` call.
Gated on the ZAID class letter: `ht ∈ {'c','h','o','r','s','a'}` calls
`read_ace_ascii` then `_acer_aplots` via a render-into-buffer pattern. Class `'t'`
(thermal) and others do `touch` (Fortran writes nothing to the plot unit for those
classes). The buffer pattern is the key safety mechanism: a not-yet-ported block
raises `AplotsNotPortedError`, the catch leaves the on-disk file as the prior empty
stub, and execution continues to tape35 — so T53/T62 keep their tape34/35 BI
status without a crash. Any other exception re-raises immediately (Rule 6). Ref:
acefc.f90:14198-14206 (acefix `aplots` gate).

### `src/NJOY.jl`

Two `include(...)` lines for `ace_reader.jl` and `ace_aplots.jl`.

## Fortran surprises

**1. Two distinct `e12.4` forms: legend-tag coords use the non-`1p` form.**
Every data and axis field in the plot tape uses the `1p` edit descriptor (Fortran
`1p,e14.6` or `1p,e12.3`), which produces a mantissa in `[1.0, 10.0)` (e.g.
`  2.950000E-01` → `2.95` × 10⁻¹). The scale/legend-tag coordinate card
(acefc.f90:15312, `write(nout,'(''4 0 2 1'',2e12.4,''/'')') xtag,ytag`) uses bare
`e12.4` WITHOUT `1p`, producing a mantissa in `[0.1, 1.0)` (e.g. `  0.9999E+00`).
Julia's `@sprintf("%12.4E", v)` gives the `1p` form and was wrong; a custom
`_e12_4_no1p` function was required. The reference line `4 0 2 1  0.9999E+00  0.1661E-01/`
confirmed the form. This was the only diff on the first run after everything else
matched.

**2. `aplof4`'s elastic curve name is `"elastic   "` (3 trailing blanks).**
Fortran declares `character(10) :: name='elastic'` (acefc.f90:17242). The string
literal has 7 characters; the declaration pads to 10 with spaces, and all 10
characters are written into the quoted label as a Fortran `a` edit descriptor.
The Julia port must use the string `"elastic   "` (3 trailing spaces) verbatim.
Omitting the padding produced a 3-character diff on that label line; it was the
first regression caught against the reference bytes.

**3. T50's ACE has `ntr=0, gpd=0, iurpt=0` — only 3 plot pages fire.**
An early research hypothesis (based on the Phase 78 HANDOFF) suggested T50 might
also emit heating/non-threshold pages. Reading the actual reference bytes disproved
this: the reference is 432 lines covering exactly 3 pages (log-log principal XS,
lin-lin principal XS, elastic angular-distribution perspective). Every guarded block
correctly skips: no reactions, no gamma production, no URR, all-zero heating array.
The 3-page structure was the bit-identity spec.

## Verification (cache-nuked, serial)

| Test | Tape | Lines | Result |
|------|------|-------|--------|
| T50  | tape33 | 432/432 | BIT_IDENTICAL — ALL PASS @1e-9 |
| T52  | tape33 | 6042/6042 | BIT_IDENTICAL — ALL PASS @1e-9 (bonus full pass) |
| T61  | tape33 | (no plot, class 't') | BIT_IDENTICAL preserved |
| T62  | tape34 | 7221/7221 | BIT_IDENTICAL preserved |
| T62  | tape35 | 1/1 | BIT_IDENTICAL preserved |
| T53  | tape34 | 12030/12030 | BIT_IDENTICAL preserved |
| T53  | tape35 | 1/1 | BIT_IDENTICAL preserved |

T53 and T62 tape33 are deferred: the log-log/lin-lin heating block (T53,
acefc.f90:16177-16226 and 16544-16595) and the aplotr threshold-reaction page (T62,
T51, T71, acefc.f90:16850-17094) are not yet ported. Both raise
`AplotsNotPortedError` and leave tape33 as the prior empty stub. No crash. No
tape34/tape35 regression.

## Beads

Bead `NJOY_jl-cnh.4` closed (Phase 79). Follow-up bead `NJOY_jl-zsh` (P2) tracks
the remaining aplots blocks: log-log + lin-lin heating (T53/T54,
acefc.f90:16177-16226, 16544-16595), aplotr threshold-reaction page (T62/T51/T71,
acefc.f90:16850-17094), aploxp particle-production (T53).

## Recommended next steps

1. **Port aplots heating + aplotr threshold-reaction block (bead NJOY_jl-zsh)** —
   the only remaining barrier to full ALL-PASS on T62 (tape33) and T53 (tape33 +
   the pre-existing T54 FP residual from Phase 78). Four further full-test flips
   from a single feature cluster.
2. **Fresh full sweep** — `reports/REFERENCE_SWEEP.md` predates the Phase-78 four
   BI flips *and* the Phase-79 T50/T52 completions; a sweep is the authoritative
   current state.
3. **T54 FP residual (bead NJOY_jl-53h)** — 1 triton recoil heating word at
   >1e-5; flips tape34 to NUMERIC_PASS@1e-5 once addressed.
4. **acelcp unported LAWs (bead NJOY_jl-5tm)** — LAW=1/6/7/Kalbach/phase-space/MF4
   for T51/T14/T71; T51 has a separate value bug (6si).
