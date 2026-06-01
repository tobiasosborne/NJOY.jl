# T54 ACER triton recoil-heating[1] — root-cause diagnosis (`_terpa_scr2` below-range branch)

**Date:** 2026-06-01
**Bead:** `NJOY_jl-53h` (ACER: T54 tape34 FP residual to BIT_IDENTICAL, 6383/7327)
**Status:** DIAGNOSIS ONLY — no code landed (session stopped before the fix/verify phases at user request).

---

## What this session was

An orchestrated diagnosis of the single >1e-5 word blocking T54 (ACER, p+H-3
triton ACE). Per the bead, T54's `acelcp` block structure is already
byte-perfect (7327/7327 lines); the residual is ~2000 words at ~7e-7 (Coulomb
1-ULP FP-floor class, same as T01/T34 — the BIT_IDENTICAL stretch) plus **one
word > 1e-5**: triton MT=2 recoil heating `[1]` at E=100 eV, where
`ref = 0.0 exactly` but Julia emits `6.008e-6`. Fixing that one word is what
flips T54 DIFFS → NUMERIC_PASS (the 7e-7 words pass at the 1e-5 bar).

Orchestration: one read-only Sonnet agent mapped the Julia `acelcp` path; a
second read-only Sonnet agent (Fortran `acefc.f90` recoil heating) was launched
in parallel but **stopped mid-run** when the session was wrapped up — its
corroboration was not needed because the Julia-side agent already cited the
governing Fortran (`terpa`, endf.f90) directly.

## Root cause (single read-only agent; NOT yet oracle-verified — see caveat)

The MT=2 recoil heating is built in `_build_andh!`
(`src/formats/ace_lcp_build.jl:499-551`): a scratch heating-kernel TAB1 `scr`
(INT=5, log-log) is filled with `(E_i, 2·amass·E_i·(1+ubar)/(1+amass)^2)`, then
interpolated onto the ESZ grid via `_terpa_scr2` and accumulated as
`heating[ie] += h · prodxs[ie]`.

The first ESZ energy (~1e-5 eV ≈ 1e-11 MeV) lies **below** the first tabulated
heating-kernel energy `Ei(1)`. The Julia interpolator does the wrong thing in
that below-range region:

```julia
# src/formats/ace_lcp_build.jl:601  (_terpa_scr2)
e <= Ei(1) && return Vi(1)     # BUG: extrapolates below-range to the first value
```

Fortran `terpa` instead returns **0** for any x below the first tabulated point
(endf.f90:1812-1816, label 170: `y=0; xnext=a(jp-2); idis=1; return`). So at the
first point Fortran yields `heating[1] = 0`, while Julia yields
`Vi(1)·prodxs[1] ≈ 6.008e-6`.

The existing `delt = _LCP_DELT = 1e-10` clamp at `ace_lcp_build.jl:240`
(`hv < delt && (hv = 0.0)`, after dividing by the total XS) does **not** catch
this — `6.008e-6 ≫ 1e-10`. The zero must come from the interpolator's
below-range semantics, not the post-divide clamp.

## Proposed fix (NOT applied)

In `_terpa_scr2` (`src/formats/ace_lcp_build.jl:601`), match Fortran `terpa`
label 170: return `0.0` for `e` strictly below the first tabulated point, and
keep the exact-match `e == Ei(1)` path returning `Vi(1)` via the normal lin/
log interp:

```julia
e < Ei(1) && return 0.0        # Fortran terpa label 170: below range → 0
```

Then (serial Julia, cache-nuked):
1. `julia --project=. test/validation/reference_test.jl 54` — confirm the word
   is now 0 and T54 status moves DIFFS → NUMERIC_PASS.
2. Regression-check **T53** (d+H-2, also exercises the shared `acelcp` /
   `_terpa_scr2` path and is currently FULL BIT_IDENTICAL) — must stay
   4636-style BI, i.e. tape34 12030/12030 unchanged.

## Caveat / what's left (honest state)

- **Not oracle-verified.** Phase 2 of the plan (reproduce the failing diff with
  `reference_test.jl 54` to pin ref=0 vs jul=6.008e-6 at the exact tape34 line)
  was **not run** this session — Law 1 (start from a red bar) is still owed
  before the fix lands. The root cause above is from one read-only Sonnet
  research agent reading source; treat as a strong hypothesis, not a confirmed
  fix. Reproduce first, then patch, then verify.
- The Fortran `terpa` semantics should be re-confirmed against
  `njoy-reference/src/endf.f90` (the below-range `y=0` branch) and against how
  Fortran `acelcp` fills/queries the heating `scr` TAB1 (`acefc.f90` ~9566-9574;
  `scr(llht+7)=5` log-log code) — the stopped Fortran-reader agent was partway
  into exactly this.
- The broader **~2000-word Coulomb 7e-7 residual** (the BIT_IDENTICAL stretch
  for T54, same FP-floor class as T01/T34) is **untouched**. Once the recoil
  word lands and T54 reaches NUMERIC_PASS, that residual is the remaining gap to
  FULL BI — likely a separate FP-order grind worth its own child bead.

## Files referenced (read-only this session)

- `src/formats/ace_lcp_build.jl` — `_build_andh!` (499-551), `_terpa_scr2`
  (593-608, bug at 601), `delt` clamp (240), `_LCP_DELT` alias.
- `src/formats/ace_lcp.jl:60` — `_LCP_DELT = 1.0e-10`.
- `njoy-reference/src/endf.f90:1812-1816` — `terpa` below-range `y=0` branch.
- `njoy-reference/src/acefc.f90` — `acelcp` recoil-heating (~9566-9574).
