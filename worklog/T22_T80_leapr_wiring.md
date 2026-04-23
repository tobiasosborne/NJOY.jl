# T22/T80 — leapr end-to-end wiring

## Date: 2026-04-23

## Outcome

**T22 (para-H₂ at 20K) → BIT_IDENTICAL.** 4636/4636 lines match the
Fortran reference tape20 at 1e-9 rtol. First leapr test to reach
bit-identical. Joins 8 existing RECONR-level bit-identical tests
(T01, T02, T08, T18, T27, T45, T47, T84).

**T80 (H_HF, 343K, 200×1325 α/β) → DIFFS, structural match.**
91453/91453 lines, correct MF1/MF7 record layouts. Numerical
agreement: 76.8% of meaningful S(α,β) values bit-identical with
Fortran, 77.9% pass 1e-5, 81.3% pass 1e-1; ~19% differ by >10%.
Tape goes from MISSING_TAPE (stub) → DIFFS (structural + mostly-right
numerics). Path to bit-identical filed as Phase 11 (port contin).

## Phase-by-phase

| Phase | Bead ID | Scope | Result |
|-------|---------|-------|--------|
| 1 | `NJOY_jl-g82` | 4-agent Fortran deep-read (trans, coldh, skold, coher, endout, parser) | Closed |
| 2 | `NJOY_jl-ixz` | `LeaprParams` 1→37 fields + full parser for T22/T23/T80 | 119/119 tests green |
| 3 | `NJOY_jl-b51` | `trans!` + stable/terps/sbfill/besk1 (translational diffusion) | 54/54 tests green |
| 4 | `NJOY_jl-n67` | `coldh!` + 8 helpers (cold H/D rotational convolution) | 114/114 tests green |
| 5 | `NJOY_jl-1xj` | `skold!` (Sköld intermolecular coherence) | Ported, wired |
| 6 | `NJOY_jl-z1t` | `coher` (Bragg coherent elastic) | Deferred (T23+) |
| 7 | `NJOY_jl-2ot` | `endout` MF1/MT451 + MF7/MT4 ENDF writer | 9/9 structural tests green |
| 8 | `NJOY_jl-8ij` | Secondary scatterer two-pass | Deferred (T23+, T09) |
| 9 | `NJOY_jl-43o` | `leapr_module` end-to-end wire + T22 green | **T22 BIT_IDENTICAL** |
| 10 | `NJOY_jl-n81` | T23/T80/T33/T09 green | In progress (T80 DIFFS) |
| 11 | `NJOY_jl-63f` | Port Fortran `contin` for T80 bit-identical | Filed (new) |

**Totals**: ~1700 LOC Julia added, 296 unit tests, 8 commits, all
pushed to master. Five leapr files authored: `leapr.jl` expansions
(trans/coldh/skold + 13 helpers), new `leapr_writer.jl` (350 LOC
MF1+MF7 emitter), extended `input_parser.jl` (parser + `*…*` string
delimiter).

## RED → GREEN (T22)

### RED test
`test/validation/reference_test.jl 22` — runs the T22 input deck
through `run_njoy` and compares the produced `tape20` against
`njoy-reference/tests/22/referenceTape20` at rtol ∈ {1e-9, 1e-7,
1e-5, 1e-3}. Pre-session: STUB wrote an empty file. After each
phase, the reference test was rerun to check progress.

### GREEN — Three sub-passes, each fixing specific diffs

**Sub-pass 1** (after Phases 2-4, 7 lands): 4617/4636 lines match.
Diffs confined to hollerith comment lines (missing leading space).
Root cause: `_leapr_strip_title` was stripping the leading space
that NJOY preserves inside `' … '` quotes. Fix: swap final
`String(strip(t))` for direct return after quote removal.

**Sub-pass 2**: 4631/4636 lines match. Three concrete gaps:
- MF1/MT451 DICT NC entry: Julia 24 vs ref 23 (off-by-one).
  Fortran's NC excludes SEND. Fix: `_leapr_mf1_nc` drops the `+ 1`.
- MF7/MT4 DICT NC entry: Julia 4607 vs ref 2315 (wildly wrong).
  Root cause: Fortran's `ncards = 5 + nbeta·(2 + (2·nalpha+4)÷6)`
  formula uses `nbeta` (not `2·nbeta-1`) and specific integer-
  division order. This is a deliberate Fortran undercount for
  isym ∈ {1,3} — even though the tape writes 2·nbeta-1 β-rows,
  the DICT entry stays at the nbeta-based count. Fix:
  `_leapr_mf7_mt4_nc` ported Fortran's exact integer arithmetic.
- LIST B-coefs position 2: Julia emitted β_max·therm (7.59); ref
  has raw β_max (300). Fortran layout: `scr(8) = beta(nbeta)`
  (raw), `scr(10) = sigfig(therm·beta(nbeta), 7, 0)` (scaled).

**Sub-pass 3** (final): 4634/4636 lines at 1e-5, two remaining
format-only diffs. Ref: `3.969624+1`; Julia: `39.6962390`. My
writer called `format_endf_float(x)` with default `extended=true`,
which enables 9-sigfig fixed-point for values in (0.1, 1e7). MF7
uses 7-sigfig scientific throughout. Fix: `_leapr_a11` now passes
`extended=false`.

Rerun: `T22 BIT_IDENTICAL 4636/4636 lines (passes at 1e-9, 1e-7,
1e-5)`.

## T80 analysis

T80 differs structurally from T22 in ways that expose `generate_sab`:

- T22: ncold=2, twt>0. The cold-H `coldh!` rotational convolution
  runs on top of `generate_sab`'s output; coldh dominates. Small
  residual generate_sab errors are masked.
- T80: ncold=0, twt>0 (minimal trans), nsk=0, nss=0. `generate_sab`
  output flows straight to the MF7 tape with no coldh to smooth it.

Result (from `/tmp/t80_diff/tape24` vs `referenceTape24`):

| Metric | Count | % of meaningful |
|--------|-------|------------------|
| All float values | 540,668 | — |
| Meaningful (|R| ≥ 1e-10) | 299,378 | 100% |
| Bit-identical | 229,996 | 76.8% |
| Pass 1e-5 | 233,112 | 77.9% |
| Pass 1e-3 | 235,014 | 78.5% |
| Pass 1e-1 | 243,405 | 81.3% |

Worst single-value diff: rel=1056 at line 55276 (R=1.19659e-13,
J=1.26556e-10 — both near `smin`, clamp boundary).

`src/processing/leapr.jl` header explicitly says "Proposal B…
AD-compatible" — it was designed as an alternative phonon-expansion,
not a transliteration of Fortran contin. Byte-exact T80 requires
porting contin directly (Phase 11).

## What this does NOT fix

- **T23 (BeO, 8 temps, secondary scatterer, iel=3)**: still blocked on
  Phase 6 (coher) + Phase 8 (secondary scatterer). Each is ~150-250
  LOC; T23 is the largest remaining leapr-only test.
- **T09 (H in H₂O, moder→reconr→broadr→leapr→thermr chain)**: leapr
  card itself uses `nss=1` (primary=H, secondary=O). Needs Phase 8.
  Also has full reconr+broadr prerequisites.
- **T33 (D in D₂O, nsk=2)**: Phase 5 (skold) is ported and wired but
  not yet exercised end-to-end.
- **T80 BIT_IDENTICAL**: Phase 11 (port contin) remaining.

## Files touched

- `src/orchestration/input_parser.jl` — LeaprParams 1→37 fields,
  parse_leapr, `*…*` string delimiter, title preserves whitespace
- `src/processing/leapr.jl` — appended ~600 LOC:
  trans!/stable/terps/sbfill/besk1 (Phase 3), coldh!/bfill/exts/sint/
  bt/sumh/cn/sjbes/terpk (Phase 4), skold! (Phase 5)
- `src/processing/leapr_writer.jl` (new, 350 LOC) — TPID + MF1/MT451
  + MF7/MT4 emitter
- `src/orchestration/modules/leapr.jl` — full replacement
  (stub → real); parse → _start → generate_sab → trans → coldh →
  [skold] → write_leapr_tape
- `src/NJOY.jl` — added `include("processing/leapr_writer.jl")` after
  `input_parser.jl` (LeaprParams dependency order)
- `test/validation/test_leapr_parser.jl` (new) — 119 tests for
  T22/T23/T80 parse fixtures
- `test/validation/test_leapr_trans.jl` (new) — 54 tests (besk1
  vs SpecialFunctions, stable vs analytic Gaussian, trans! smoke)
- `test/validation/test_leapr_coldh.jl` (new) — 114 tests (helpers +
  smoke)
- `test/validation/test_leapr_writer.jl` (new) — 9 structural tests
- `HANDOFF.md` — T22 + T80 status rows

## Regression check

- T22 reference test: BIT_IDENTICAL (4636/4636 @ 1e-9) — verified
  holding after each phase commit, including after skold addition
- All 296 new unit tests green
- Existing RECONR/BROADR/ERRORR etc. unchanged (no edits to those
  files); tokenizer change for `*` scoped to mean "alternative
  string delimiter in module blocks, comment rule preserved for
  single-`*` lines" so other module decks unaffected

## Reference

- Fortran leapr: `njoy-reference/src/leapr.f90:1-3625`
- Particular subroutines:
  - leapr main entry: `:55-453`
  - input deck parse: `:230-400`
  - trans + helpers: `:844-1318`
  - bfill/exts/sint: `:1798-1934`
  - coldh + helpers: `:1936-2466`
  - skold: `:2816-2922`
  - coher: `:2489-2814` (deferred)
  - endout: `:2972-3623`
- Reference tapes: `njoy-reference/tests/22/referenceTape20`,
  `njoy-reference/tests/80/referenceTape24`
- Beads memory: `bd memories leapr` surfaces the MF7-NC-quirk entry
  from this session
