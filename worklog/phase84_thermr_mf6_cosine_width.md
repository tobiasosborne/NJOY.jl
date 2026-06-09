# Phase 84 — Oracle realign + thermr MF6/MT221 cosine-width fix

**Date:** 2026-06-09
**Beads:** `NJOY_jl-c3q` (premise corrected, width fix landed, stays in_progress), `NJOY_jl-rby` (setup_reference.sh fetch bug — closed)
**Outcome:** T70 tape60 **112818 → 213826 lines** (ref 213800), match **25562 → 125483 / 213800** (22% → 59%). T01 NUMERIC_PASS + T09 BIT_IDENTICAL preserved. One mechanical fix corrects the thermal scattering matrix for **all 10 `nbin=20` tests**.

Orchestrated session (parallel read-only Sonnet researchers, Opus coder, serial Julia run by the orchestrator; verified-before-trust).

---

## 0. Blocker found first: the Fortran oracle was a ghost commit

`njoy-reference` HEAD was `ac5adf5` — the rebased-away 2026-05-30 baseline (the
`njoy-reference-oracle-pin` memory flags it as unreachable in upstream history).
The canonical pin is `2c64dfb`. `bash scripts/setup_reference.sh` **failed**:
`git checkout 2c64dfb` aborted with `fatal: reference is not a tree`.

Root cause: the script did a plain `git fetch --quiet origin` (branch tips only).
Upstream `develop` rebases, so the pin (a March commit) is on **no** branch tip
and a plain fetch never downloads it. A direct `git fetch origin <SHA>` **does**
retrieve it (GitHub serves any reachable-by-SHA object that way).

**Fix (`scripts/setup_reference.sh`):** after the plain `git fetch origin`, also
`git fetch origin "${SHA}"` (both guarded `|| true`). Manually proven: direct
fetch of `2c64dfb` succeeded, checkout aligned the oracle. Bead `NJOY_jl-rby`
closed. **Lesson: always align the oracle FIRST — comparisons against a ghost
oracle are unsound.** Note T68/T70 reference tapes are byte-identical between
`ac5adf5` and the pin (thermr is pin-insensitive; the only known pin diff is the
leapr EMAX of PR #372), but the alignment is still mandatory for rigor.

## 1. Rule-2 catch: the c3q premise was inverted

`c3q` claimed "MF3 grid over-produced (Julia 213800 vs ref 112817)". Against the
**corrected** oracle the numbers are **swapped**:

- T70: ref tape60 = **213800**, Julia = **112818** → Julia **UNDER**-produces.
- The bulk is the **MF6/MT221 thermal scattering MATRIX** (LAW=1, NE=101), not an
  MF3 grid: ref 185295 lines, Julia 84312 (~46%). MF3/MT221 is only 498 lines.
- "112817/112818" was Julia's own count; "213800" is the reference. The Phase-83
  triage "correction" (claiming ref=112818 at the pin) was itself wrong.

So `c3q`'s "lat=10 coherent / MF3 grid over-production" framing was a wrong-oracle
artifact. HANDOFF/bead disagreed with the source — trust the (correctly pinned)
Fortran (CLAUDE.md Rule 2).

## 2. Root cause (verified + empirically exact)

The MF6/MT221 LIST records store a per-secondary-energy row of **10** values
(E′, σ, μ₁…μ₈) — hardcoded — instead of **nbin+2**:

- `src/processing/thermr.jl:1004` — `entries::Vector{NTuple{10, Float64}}`
- `src/processing/thermr.jl` packing — `ntuple(k -> …, 10)` (5 sites across
  `calcem`, `calcem_free_gas`, legacy `compute_mf6_thermal`)
- `src/processing/pendf_writer.jl:598` — `nl = 10`; and the directory-NC at ~300

Fortran ground truth (`thermr.f90`): `nl = nbin+1` (1630), `scr(6) = nl+1`
(2230), `NW = (nl+1)·j` (2225-2232) → **nbin+2** words per secondary energy.

For `nbin=20`: 10/22 = **0.4545** of correct size. Observed deficit
84312/185295 = **0.455** — exact. Empirical confirmation of the first
incident-energy LIST: ref `NW=10120 stride=22 → NEP=460`; Julia
`NW=4600 stride=10 → NEP=460`. **Same 101 incident × 460 secondary energies —
only the per-row width was wrong.** The matrix was never missing rows; each row
was missing 12 of its 22 values.

## 3. Fix

`src/processing/thermr.jl`
- New helper `_mf6_entry(ep, σ, cosines, nbin) -> Vector{Float64}` = single source
  of truth; builds `[E′, σ, μ₁…μ_nbin]` of length nbin+2; zero-pads if
  `length(cosines) < nbin` (mirrors the old `k-2<=nbin ? cosines[k-2] : 0.0`).
- `MF6ListRecord.entries`: `Vector{NTuple{10,Float64}}` → `Vector{Vector{Float64}}`.
- All 5 packing sites + 3 `entries =` initializers converted; the Fortran
  `j==3 / total_xs<tolmin_area` overwrite branch preserved verbatim.

`src/processing/pendf_writer.jl`
- LIST writer + directory-NC: `nl = n_ep>0 ? length(rec.entries[1]) : 2` (= nbin+2);
  `nw = n_ep*nl`. Inner write loop (sigma 7/6-sigfig norm, cosines 9-sigfig) and
  `n_lines += 1 + cld(nw,6)` unchanged.

**Faithfulness:** for `nbin=8`, `_mf6_entry` yields exactly `[E′,σ,μ₁…μ₈]`
(length 10) — byte-identical to the old truncation. T01/T09/T11 (all `nbin=8`)
are protected by construction.

## 4. Verification (serial, cache-nuked, oracle at pin)

- **T70** (Al-27, nbin=20): tape60 **112818 → 213826** lines (ref 213800),
  match **25562 → 125483 / 213800**. tape55/71/76 still 0 (acer iopt=2 stub).
- **T01** (nbin=8): tape25 **NUMERIC_PASS 32812/32962** — unchanged.
- **T09** (nbin=8): tape24 **BIT_IDENTICAL 1830/1830** — preserved.

## 5. Impact + remaining

One fix corrects the thermal matrix structure for all `nbin=20` thermr tests:
**T24, T25, T32, T49, T67, T68, T69, T70, T72, T74** (the sweep "TIMEOUTs"
T67/68/74 were the 300s cap + contention — T68 completes in 265s, T70 in 143s —
not hangs).

Remaining on T70 tape60 (folded into `c3q`, still in_progress):
1. MF3/MT221+MT222 grid off-by-one (Julia 497 vs ref 498; T68 MF3 213 vs 214) —
   separate pre-existing grid issue.
2. MF6/MT221 +26-line overshoot (185320 vs 185294) — minor, introduced-adjacent;
   investigate the trailing-point/dedup count.
3. Matrix VALUE grind to 1e-7/1e-9 — the ~41% non-matching lines are σ/cosine FP
   (T34-class).
4. tape55/71/76 EMPTY → `acer iopt=2` thermal not ported (bead `NJOY_jl-p9q`).
