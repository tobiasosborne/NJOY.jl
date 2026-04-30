# T05 / T16 — Full COVR Port (Phase 55)

**Date:** 2026-04-30
**Goal:** Replace the empty-tape `covr_module` stub with a complete
Fortran-faithful port of `njoy-reference/src/covr.f90` (2249 lines, 18
subroutines), bit-identical to the reference plot tapes.

## Outcome

**3/3 covr-isolation tapes BIT-IDENTICAL:**

| Tape                                        | Bytes (J/Ref)         | Status         |
|---------------------------------------------|-----------------------|----------------|
| T05 referenceTape34 (C-12 MF33 plot)        | 1 839 463 / 1 839 463 | BIT-IDENTICAL  |
| T16 referenceTape36 (U-238 MF33 plot)       |   270 922 /   270 922 | BIT-IDENTICAL  |
| T16 referenceTape37 (U-238 MF34 mubar plot) |    23 768 /    23 768 | BIT-IDENTICAL  |

`test/validation/test_covr_isolated.jl` locks all three in. T01
(NUMERIC_PASS 32812/32962) and T22 (BIT-IDENTICAL 4636/4636) regression-
checked: unchanged from prior baseline.

The full-pipeline T05/T16 reference tests still diff because Julia's
`errorr` (HANDOFF P1, "Covcal content drift") emits a sparse subset of
the reference covariance tape; covr's input is therefore wrong. When
errorr's NK-stub + LB=5 union-grid path lands, T05/T16 will pass end-to-
end without further covr changes.

## Files

NEW (port subroutines):
- `src/processing/covr_labels.jl`  — `elem`, `mtno`, `matmes`, `smilab`
- `src/processing/covr_io.jl`      — `covard`-equivalent (errorr tape reader)
- `src/processing/covr_corr.jl`    — `corr`, `truncg`
- `src/processing/covr_plot.jl`    — `plotit`, `matshd`, `level`, `patlev`
- `src/processing/covr_press.jl`   — `press`, `setfor` (boxer-format library)

REWRITTEN:
- `src/orchestration/modules/covr.jl` — main `covr_module` (mirrors
  Fortran `covr` outer loop), library + plot mode dispatch
- `src/orchestration/input_parser.jl` — full `CovrParams` + `parse_covr`
  reading every Fortran card (cards 1, 2, 2', 2a, 3a, 2b, 3b, 3c, 4),
  including the `//` multi-card-per-line shorthand for
  consecutive-default-cards
- `src/NJOY.jl` — added 5 includes after `pendf_io.jl`

NEW (test):
- `test/validation/test_covr_isolated.jl` — three @testset blocks
  feeding `referenceTapeNN` directly into `covr_module` and asserting
  byte-for-byte equality with `referenceTape{plot_unit}`.

## Bugs landed (each oracle-driven)

1. **Title padding** (covr.f90:1025 declares `character(80)::strng`).
   Fortran's `write '(a)'` emits the full 80 chars including trailing
   spaces. Julia originally emitted only the active text. Fix: rpad
   each plot title to 80 before writing. Diff dropped from 3088 lines
   to ~1700.

2. **MT name table row 13 typo** — entered `)g<)` instead of `]g<)`.
   The right-bracket-gamma form is the (n,γ) reaction marker. Single-
   char typo collapsed the entire MT=102 frame title. Fix in
   `_MTNO_HIRA[13]`.

3. **`(mt NNN)` fallback width** — Fortran writes the digits into cols
   4-6 of the fixed-width `(mt     ` template (i3 right-justified):
   for MT=251 → `(mt251  `. My Julia did `"(mt " * lpad(...,3)` which
   shifts the digits right, producing `(mt 251 ` and the truncated
   inamel=6 view becomes `(mt 25` — visible in T16's mubar title.
   Fix: `"(mt" * lpad(string(mt), 3) * "  "`.

4. **Truncated-grid indexing in `plotit`** — Fortran calls
   `plotit(x(ixmin:),y(ixmin:),...,rsdx(ixmin:),rsdy(ixmin:),...)` so
   inside the routine, local index 1 maps to global index `ixmin`. My
   Julia passed the full ixmax-length arrays and indexed with the
   local `i`, producing wrong x-axis bounds when `ixmin > 1`. Fix:
   slice `xs_x = view(cr.x, ixmin:ixmax+1)`, `ys_y = view(cr.y, ixmin:ixmax+1)`,
   and `rsdx = collect(view(cr.rsdx, ixmin:ixmax))`. (Same fix for
   rsdy.) Surfaced as the entire (n,inel) frame plotting from the
   wrong starting energy on T05.

5. **Spurious NC sub-subsection skip** in `_read_mf33_subsection`. My
   reader interpreted the sub-section CONT's N1 field as `NC` and
   consumed N1 LIST records before the matched-data read. For MF34,
   errorr stores `NSS=1` in N1 (the number of sub-subsections per
   (mat1,mt1) pair), so my code ate the first row of MF34/MT=251's
   covariance, producing rsdy[13]=0 instead of `sqrt(1.293e-3)≈0.036`
   and starting the rotated frame one group too late. Fortran covard
   never reads N1 — it only uses N2 (=NI=`kgp`) and reads kgp many
   LIST records terminating when `kl ≥ kgp`. Fix: removed the
   `for _ in 1:nc: skip_one_list` loop. Caught by T16 tape37 frame 2.

6. **Empty-default cards eaten by tokeniser** — `parse_njoy_input`
   discards cards with no tokens, but covr decks rely on `/` (and `//`)
   as explicit "use defaults" markers. Without the raw line trail,
   `parse_covr` couldn't tell card 2a from card 3a when both were
   defaulted. Fix: add `raw_lines::Vector{String}` to `ModuleCall`,
   and have `parse_covr` walk those, splitting each line on every `/`
   so every terminator yields one card (empty or not). Tokeniser
   behaviour for non-covr modules is unchanged (they still read
   `mc.raw_cards`).

## Acceptance criteria — met

- T05 referenceTape34 (78 796 lines) bit-identical via covr_module
- T16 referenceTape36 (14 082 lines) bit-identical
- T16 referenceTape37 (1 234 lines) bit-identical
- Library/boxer-format mode (`press` + `setfor`) ported but no test
  coverage yet (no T## currently exercises `nout > 0`)
- T01 NUMERIC_PASS preserved
- T22 BIT-IDENTICAL preserved

## Next steps (out of scope for this phase)

- End-to-end T05/T16 will pass when errorr's "Covcal content drift"
  (HANDOFF P1) lands — no further covr work needed.
- Library-mode (boxer) regression test once a candidate test surfaces.

## Fortran source citations

All ported routines cite line ranges of `njoy-reference/src/covr.f90`:
`covr` 48-505, `expndo` 508-576, `corr` 578-718, `covard` 720-937,
`truncg` 939-1011, `plotit` 1014-1308, `matshd` 1310-1599, `level`
1601-1619, `patlev` 1621-1647, `smilab` 1649-1703, `matmes` 1705-1751,
`elem` 1753-1779, `mtno` 1781-1890, `press` 1991-2218, `setfor`
2220-2247.
