# Phase 81 — T54 ACER tape34 DIFFS → NUMERIC_PASS (two real bugs); tape33 FP-floor reframe

**Date:** 2026-06-02
**Bead:** `NJOY_jl-53h` (ACER: T54 tape34 FP residual to BIT_IDENTICAL)
**Commits:** `648c6cc` (recoil below-range), `8c622c7` (AND-locator integer format)
**Outcome:** T54 **tape34 DIFFS → NUMERIC_PASS** (7326/7327 @1e-5; 1e-9 6383→6398). Zero regressions. T54 the *test* still DIFFS (tape33 FP-floor, reframed below).

Orchestrated: 2 Opus coding subagents (one per bug, each running the full Law-1 reproduce→fix→verify→regress loop, serial Julia), orchestrator verified every result independently (Rule 2) and committed.

---

## Bug 1 — `_terpa_scr2` below-range branch (recoil heating[1])

`_terpa_scr2` (`src/formats/ace_lcp_build.jl:601`) is a port of Fortran `terpa`
(endf.f90). It did `e <= Ei(1) && return Vi(1)`, extrapolating the first
tabulated value into the below-range region. Fortran `terpa` returns **y=0** for
x strictly below the first point (endf.f90:1812-1817, label 170) and the first
value only *at* the point (label 140). The first ESZ energy (~1e-11 MeV) is
below the first heating-kernel energy Ei(1), so the triton MT=2 recoil
heating[1] came out `6.008442e-6` where ref = 0 exactly.

**Fix:** `e < Ei(1) && return 0.0` (e == Ei(1) now falls through to `_terp1` →
Vi(1), matching label 140). Upper-bound line left untouched (T53's shared LAW=2/4
path relies on it). tape34 lines 153 + 4733 → 0. tape34 1e-9 6383→6385.

## Bug 2 — AND-block LC locators written as ENDF floats (the bead premise was incomplete)

The first Opus agent honestly flagged that the bead's "byte-perfect structure,
one word > 1e-5" premise was wrong: tape34 had **15** lines failing at 1e-5 — 2
numeric (recoil, Bug 1) + **13 residue** (lines 4785–4797, 52 words). The
reference writes the AND-block tabulated-distribution LC locators `-106..-4543`
as right-justified i20 **integers**; Julia emitted ENDF floats.

Root cause: `_build_andh!` (`ace_lcp_build.jl:502`) did
`_set!(x, nb + ie, -(next - andh + 1))` — **missing `isint=true`**, unlike every
sibling locator write (the direct analog `_andh_law2!:675`, and `landh` writes at
483/623, all `isint=true`). Fortran: `acefc.f90:9520` `xss(nb+ie)=-(next-andh+1)`;
ANDH write path `acefc.f90:13724` `write_integer_list` → `acecm.f90 typen`
iflag=1 → `write(...,'(i20)')`.

**Why the FULL-BI cohort (T53/T50/T52/T62) was unaffected:** their angular data
goes through `_andh_law2!`/other (isint) paths; only T54's data reaches the
`_build_andh!` line-502 path (`_build_andh!` is called only from
`ace_lcp_build.jl:204`, the LCP builder). Confirmed: T53/T50/T52/T62 stayed FULL
BIT_IDENTICAL, T02/T08 unchanged.

**Fix:** add `isint=true` at line 502. tape34 1e-9 6385→6398, **NUMERIC_PASS**.

---

## tape33 reframe — T54 the TEST is still DIFFS, and the HANDOFF claim was wrong

The prior HANDOFF asserted "T54 tape33 flips to FULL BI once the tape34 recoil
word lands." **False.** After both fixes, `reference_test.jl 54` (orchestrator's
own run) still reports:

```
tape33  DIFFS         9919/11340
tape34  NUMERIC_PASS  6398/7327   (passes at 1e-5)
tape35  BIT_IDENTICAL 1/1
```

Robust paired comparison of tape33 (ref vs jul): **zero** token-count or
non-numeric residue mismatches; max numeric reldiff **1.84e-5**; only **4 lines**
exceed 1e-5 (3109 1.32e-5, 3115 1.14e-5, 3130 1.35e-5, 3132 1.84e-5), all in the
**'recoil heating'** plot curve (title at line 3040). aplots builds that curve as
`esz_heat − Σ particle_heating` (a subtraction); the underlying ~7e-7
charged-heating FP residual (same floor class as T01/T34, the acer_charged_elastic
Coulomb total/elastic noise) is **cancellation-amplified** to ~1.8e-5 in the small
difference, pushing those 4 points just over the 1e-5 bar → tape33 DIFFS instead
of NUMERIC_PASS. The identical residual is tape34's 949-line gap at 1e-9.

So a single root cause — the ~7e-7 charged heating/elastic FP residual — blocks
both tape33 (NUMERIC_PASS) and tape34 (BIT_IDENTICAL). The recoil + locator
portion of bead 53h is **done**; the FP-residual portion remains (deep grind,
deferred).

### Remaining work to flip T54 → NUMERIC_PASS (bead 53h)
- Easiest lever: match Fortran aplots' exact `esz_heat − Σ particle_heating`
  accumulation order in `src/formats/ace_aplots.jl` (the recoil-heating curve) so
  the 4 outliers drop below 1e-5 — without needing to kill the ACE-level residual.
- Harder / full BI: reduce the ~7e-7 charged heating/elastic FP residual at the
  source (acer_charged_elastic Coulomb total/elastic; T01/T34 class).

---

## Lessons
- **The HANDOFF was wrong twice over** (Rule 2): it under-counted the tape34
  failures (missed the 13 locator-format lines) and mis-attributed tape33 to "the
  recoil word." A subagent's reproduce step + the orchestrator's own paired diff
  caught both. Trust the run, not the worklog.
- **A subtractive plot curve amplifies FP floor.** tape34's stored heating passes
  1e-5, but aplots' `esz_heat − Σparticle` difference of those values spikes to
  1.84e-5. A NUMERIC_PASS ACE can still yield a DIFFS plot tape.
- The "one word" framing of FP-residual beads is a trap — always re-measure the
  full failing-line set against the live reference before scoping.
