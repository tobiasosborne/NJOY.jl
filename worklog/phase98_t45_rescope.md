# Phase 98 — T45 current-state re-scope

Date: 2026-08-22

Diagnosis and handoff only. No production or test code changed in this phase.
The user requested a clean stopping point after updating the handoff, commit,
and push.

## Outcome

The stale `NJOY_jl-1kf` three-record metadata premise was replaced with a
fixed-column audit of the current Phase 97 tree.

- Oracle and Julia tape40 both contain **7,188 records**.
- All section boundaries, grids, trailers, MF1/MT451 records 2-69, MF2, and
  MF12/MF13 photon bodies are exact.
- Raw comparison has **412 differing records**; the official Phase 97 T45 run
  reports 6,777/7,188 strict matches because one raw residual is hidden by the
  comparator tolerance.
- No implementation was started, so the repository stops at a reproducible
  diagnostic boundary.

## Exact residual map

| Scope | Records | First-last line | Classification |
|---|---:|---:|---|
| TPID | 1 | 1 | wrong default text |
| MF3/MT1 | 1 | 80 | upstream BROADR ULP |
| MF3/MT102 | 2 | 2371-2563 | upstream BROADR ULP |
| MF3/MT113 | 1 | 3713 | upstream BROADR ULP |
| MF3/MT800 | 1 | 6015 | upstream BROADR ULP |
| MF3/MT203 | 42 | 4302-4343 | missing GASPR yields |
| MF3/MT204 | 49 | 4365-4413 | missing GASPR yields |
| MF3/MT205 | 1 | 4419 | inherited MT113 ULP |
| MF3/MT207 | 314 | 4736-5049 | missing GASPR yields |

All MF3 energy fields agree. The 405 substantive GASPR residuals are value
errors, not grids or serialization.

## TPID ownership

T45's RECONR input card supplies the two-blank label `'  '`. Fortran writes a
blank TPID in `reconr.f90:165-192`; BROADR, GASPR, and final MODER copy it
unchanged (`broadr.f90:334-340`, `gaspr.f90:229-232`,
`moder.f90:142-158`).

Julia strips the label, maps the empty title to `nothing`, and then substitutes
`"reconstructed data"` in the PENDF writer. `NJOY_jl-1kf` is re-scoped to
this single exact line-1 metadata defect. Its old structural acceptance is
already satisfied.

## GASPR source gaps

Fortran scans every eligible MF3 reaction on the MT1 grid, in tape order, and
rounds only when writing (`gaspr.f90:280-869,1055-1109`). Two absent semantics
explain the large clusters:

1. **MT51-91 TAB1 LR yields.** LR=28 adds one proton; LR=35 adds one deuteron
   and two alphas (`gaspr.f90:330-346,567-611`). This starts T45 MT203 at
   8.9 MeV and MT204 at 7.2 MeV and contributes to high-energy MT207.
2. **Residual-nuclide gas.** Fortran adds gas carried by the residual identity
   after the primary reaction mapping (`gaspr.f90:821-826`). In T45, MT105's
   residual O-8 adds two alphas beyond the explicit triton. Julia's static
   MT-only multiplicity table handles the explicit MT113 mapping but omits
   this residual addition, leaving MT207 low across all 941 points.

Tracked as `NJOY_jl-9h4`. Its core scope is MT203/204/207; the single MT205
record follows the upstream MT113 ULP.

## BROADR ULP lane

GASPR copies pre-gas MF3 sections unchanged, so the five isolated records are
an independent BROADR grind. Source-backed candidates for instrumentation:

- Fortran `mathm.erfc` Chebyshev implementation versus `SpecialFunctions.erfc`;
- Fortran's truncated sqrt-two constants in `hnabb` versus full-precision
  `sqrt(2.0)`;
- Fortran `hunky`'s independent `if` statements and `<=` boundary versus
  Julia's `elseif` and `<`.

No candidate is accepted without a raw `tt` trace. Tracked as
`NJOY_jl-bdu`, with the exact five tape lines in its acceptance criteria.

## Beads after re-scope

- `NJOY_jl-1kf` — blank TPID line 1 only; open.
- `NJOY_jl-9h4` — GASPR LR/residual-nuclide yields; open, P2.
- `NJOY_jl-bdu` — five isolated BROADR ULP records; open, P3.

## Next

Resume with `NJOY_jl-1kf` as the smallest exact red/green change, then take
`NJOY_jl-9h4` through the core 3+1 workflow. Keep `NJOY_jl-bdu` separate and
instrument the Fortran raw broadening state before changing any constants or
branch boundaries.
