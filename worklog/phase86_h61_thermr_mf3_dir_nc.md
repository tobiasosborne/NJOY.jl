# Phase 86 — T70 thermr MF1/MT451 directory-NC off-by-one (bead `NJOY_jl-h61`)

**Date:** 2026-06-14
**Outcome:** T70 tape60 MF3/MT221+222 directory NC **497 → 498** (== ref); first
structural diff moved **line 69 → line 223** (now purely the c3q value-grind);
total 213800 == ref; matching lines **128064 → 128066** (+2 = exactly the two
fixed directory lines). T01 NUMERIC_PASS 32812/32962 + T68 directory NC=129
preserved. No regression.

## Orchestration

- **Phase A — research (parallel, read-only, NO Julia; Rule 9):** 3 sonnet agents,
  non-overlapping — (1) Fortran `tpend` NC ground truth, (2) Julia data-flow trace
  + git-blame, (3) reference-tape verification (T70/T68/T01 NP + directory NC).
- **Phase B–D — TDD (serial Julia, one opus agent):** RED (`reference_test.jl 70`
  first-diff at line 69, 497 vs 498) → GREEN (formula fix) → unit test → regress.
- Orchestrator independently re-verified the agent's report (Rule 2): re-ran the
  unit test (10/10), `reference_test.jl 70` (directory fixed, first-diff→223), and
  `reference_test.jl 1` (NUMERIC_PASS 32812/32962).

## Root cause

The PENDF writer's MF1/MT451 directory-NC loop (`pendf_writer.jl:276-284`) had a
thermr-specific branch that substituted a **stale, smaller `thermr_coh_ne`**
(~1482, read from a different pipeline stage in `_recompute_thermr_mf6!`,
`pipeline.jl:580` = `length(mf3[mtref][1]) - 2`) in place of the **actual emitted
np** (1485) for MF3/MT221+222:

```julia
if imt in thermr_mts && thermr_coh_ne >= 0
    ne_param = np > thermr_coh_ne ? thermr_coh_ne : np
    section_lines[(3, imt)] = 3 + div(ne_param + 2, 3)   # → 497
else
    section_lines[(3, imt)] = 3 + cld(np, 3)
end
```

The branch was introduced in commit `775092c` ("T01 pipeline grind … MF1 NC …")
and the `thermr_coh_ne` setter added later in `2bc4cb8`; T70 was STRUCTURAL_FAIL at
the time and never a validation target. T01 was only *accidentally* correct — its
`thermr_coh_ne` got set large by the graphite second-thermr (mtref=229), so
`np > thermr_coh_ne` was false and MT221 silently fell through to the right value.

## Fortran ground truth (LAW 2) — overturned the bead's own suggested fix

The bead (and the Julia data-flow agent) proposed `3 + cld(np, 3)`. **Rejected per
Rule 1.** Fortran thermr `tpend` computes:

```
nc = 3 + (ne+2)/3            ! thermr.f90:3087  (integer division)
scr(6) = ne+1                ! thermr.f90:3198  (written NP = ne+1)
if (ib.ge.ne+1) ex(1)=etop   ! thermr.f90:3211-3213  (one etop=20e6 sentinel)
```

So `ne = NP_written - 1`: tpend appends **one** `etop=20e6` sentinel that is counted
in the written NP field but **not** in `ne`. In terms of the emitted np:

```
nc = 3 + div((np - 1) + 2, 3) = 3 + div(np + 1, 3)
```

This **diverges from `3 + cld(np,3) = 3 + div(np+2,3)` exactly when `np ≡ 1 (mod 3)`**,
where Fortran *undercounts the true record count by 1*. None of T70 (1485≡0),
T68 (378≡0), T01 (146≡2) hit that case, so the two formulas agree on every current
test — but the Fortran-faithful `div(np+1,3)` is canonical, and the undercount is a
quirk we **reproduce, not "fix"** (THE GROUND-TRUTH PRINCIPLE). reconr's own tpend
(`reconr.f90:5115-5116`, `3+int((np+2)/3)`) has no appended sentinel → the else
branch correctly keeps `3 + cld(np, 3)`.

## The fix (`src/processing/pendf_writer.jl`)

Two literate, Fortran-cited helpers + a one-line branch replacement:

```julia
_thermr_mf3_dir_nc(np::Int) = 3 + div(np + 1, 3)   # thermr.f90:3087/3198
_reconr_mf3_dir_nc(np::Int) = 3 + cld(np, 3)        # reconr.f90:5115-5116
...
section_lines[(3, imt)] = imt in thermr_mts ?
    _thermr_mf3_dir_nc(np) : _reconr_mf3_dir_nc(np)
```

`thermr_coh_ne` dropped from both the condition and the formula (kwarg left in the
signature with a deprecation note → removal tracked by `NJOY_jl-czw`).

## Tests

- New `test/validation/test_thermr_mf3_nc.jl` (10/10): asserts `_thermr_mf3_dir_nc`
  gives 498/52/129 for np = 1485/146/378, **498 (not 499) for np=1486** (the
  `np≡1 mod 3` undercount), the `_reconr_mf3_dir_nc(1486)==499` contrast, and
  agreement of both forms for np ≢ 1 mod 3.
- `reference_test.jl 70`: tape60 first-diff line 69 → 223; total 213800; match
  128064 → 128066. (tape55/71/76 STRUCTURAL_FAIL are pre-existing acer-thermal
  stubs under c3q/p9q, unrelated.)
- `reference_test.jl 1`: NUMERIC_PASS 32812/32962 @1e-5 (unchanged).
- T68 end-to-end TIMEOUTs in a downstream thermr stage (pre-existing, per
  `reports/REFERENCE_SWEEP.md`); its directory NC=129 confirmed via the reference
  tape (`tests/68/referenceTape60:13`) + the unit test instead.

## Follow-ups

- `NJOY_jl-czw` (P3, depends on h61): remove the now-vestigial `thermr_coh_ne`
  plumbing (kwarg, RunContext field/default, setter, moder.jl call sites, doc).
- T70 tape60 first-diff at line 223 onward is the c3q σ/cosine value-grind (±1 ulp);
  remains under `NJOY_jl-c3q`.
