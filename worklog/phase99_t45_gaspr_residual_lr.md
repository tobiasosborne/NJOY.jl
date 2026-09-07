# Phase 99 — T45 GASPR residual/LR gas yields + blank TPID

**Date:** 2026-09-07
**Beads:** `NJOY_jl-9h4` (closed), `NJOY_jl-1kf` (closed), `NJOY_jl-xpf` (confirmed, notes updated),
`NJOY_jl-jbi` + `NJOY_jl-<egas>` (filed)
**Result:** T45 tape40 residual records **412 → 6**; T45 `DIFFS → NUMERIC_PASS`.

## Starting point

Phase 98 stopped at diagnosis. Re-running T45 on the pushed Phase 98 tree
reproduced its audit exactly: tape40 is structurally exact (7,188/7,188 lines,
582,228 bytes — the same size as the oracle) with 412 differing records:

| Lane | Records | Bead |
|---|---:|---|
| MF3/MT207 (alpha) | 314 | `9h4` |
| MF3/MT204 (deuteron) | 49 | `9h4` |
| MF3/MT203 (proton) | 42 | `9h4` |
| TPID line 1 | 1 | `1kf` |
| MF3 MT1/MT102×2/MT113/MT800 | 5 | `bdu` |
| MF3/MT205 (mirrors the MT113 ULP) | 1 | `bdu` |

## What the oracle said before any code was read

T45 is B-10 (MAT 525, ZA 5010). Decomposing the oracle's MT207 against the
other MF3 sections on the same 941-point grid pinned the defect to seven
significant digits *before* opening the Fortran:

```
ref MT207 = MT107 + 2*MT113 + 2*MT105      (agrees to <2.2e-7, i.e. tape quantisation)
jl  MT207 = MT107 + 2*MT113                (exactly what Julia produced)
```

So Julia was short by `2 * sigma(MT105)` at every energy. MT203/MT204 diverged
only above 8.806 MeV and 7.155 MeV, which are precisely the thresholds of MT65
(LR=28) and MT62 (LR=35).

## Ground truth (Law 2)

`njoy-reference/src/gaspr.f90` does **not** carry an MT → multiplicity table.
For every MF3 reaction it tracks two things at once (gaspr.f90:495, 507-818):

* `izr` — the residual-nucleus ZA, initialised to `za + zain` and decremented
  by the ZA of each emitted particle;
* `y203..y207` — the particles the reaction *explicitly* emits.

and only then converts a light residual into gas (gaspr.f90:820-825):

```fortran
if (izr.eq.1001) y203=y203+1     if (izr.eq.2003) y206=y206+1
if (izr.eq.1002) y204=y204+1     if (izr.eq.2004) y207=y207+1
if (izr.eq.1003) y205=y205+1     if (izr.eq.4008) y207=y207+2
```

The last line is the whole bug. For B-10 MT105 (n,t):
`izr = 5010 + 1 - 1003 = 4008` — Be-8, which is unbound and breaks into two
alphas. Hence 2 alphas on top of the explicit triton, exactly the missing
`2*sigma(MT105)`.

For MT51-91 the gas comes entirely from the LR breakup flag in the MF3 TAB1
L2 field, applied after the emitted neutron is removed (gaspr.f90:565-611).
B-10's evaluation carries LR=22 on 14 levels, LR=35 on 14, LR=28 on 2.

Julia's `gas_multiplicity` was a static MT → yield table with neither the
residual arithmetic nor any LR branch, so all 30 LR-bearing levels contributed
nothing and every residual-nucleus alpha was lost.

Two further findings from reading the filter (gaspr.f90:475-491):

* Levels MT600-849 are skipped so the MT103-107 totals are not double counted.
  This is why B-10's MT700 (n,t0), whose values equal MT105's, is *not* also
  summed — the tape alone cannot distinguish the two paths, only the source can.
* Julia's ad-hoc `_GASPR_SKIP_MTS` wrongly skipped MT16/17/37 and MT102, all of
  which the Fortran keeps. They contribute nothing for B-10, but MT16 on Be-9
  leaves Be-8 and MT102 on H-1 leaves a deuteron, so the set was latently wrong.

## Changes

**`src/processing/gaspr.jl`** — replaced the static table with a faithful port:
`_GASPR_MT_CHANNEL` and `_GASPR_LR_CHANNEL` pair each reaction's ΔZA with its
explicit emissions in one table (so the two halves cannot drift),
`_GASPR_RESIDUAL_GAS` holds the six residual rules, and `gaspr_skips_mt`
implements the Fortran's real filter including the iverf-dependent level bands.
New `gas_channel(mt, lr, za, zain)` returns what GASPR actually accumulates;
`gas_residual_za` and `gas_threshold_candidate` expose the other two Fortran
predicates. `gas_multiplicity`/`gas_yield` are kept, now documented as *only*
the explicitly-emitted half.

**`src/orchestration/modules/gaspr.jl`** — reads `za`/`zain`/`iverf` from the
original ENDF MF1/MT451 via `read_mf1_header_info` (Fortran does the same at
gaspr.f90:82-102, and uses that same `za` for the output section header at
1067-1068), parses each section's LR from TAB1 cols 34-44, and switches the
threshold and accumulation loops to the new predicates. Also **sorts the
accumulation MTs ascending**: it was iterating a `Dict`, so the floating-point
summation order into `sgas` was nondeterministic — invisible while only three
reactions contributed, a live hazard once the LR terms were added.

**`src/orchestration/modules/reconr.jl`** — dropped `isempty(title) ? nothing :
title`. RECONR's label card is mandatory (reconr.f90:168) and its content is
written verbatim to the TPID (reconr.f90:192), then passed through unchanged by
BROADR/GASPR/MODER. An empty parsed title means "the deck asked for a blank
label", never "substitute a default"; T45's deck line is literally `'  '/`.
Surveyed all 86 decks first: every one supplies a title card and T45 is the only
blank, so nothing depended on the `"reconstructed data"` fallback. (That string
is a real title in the T83-T86 decks, which pass it explicitly and are
unaffected.)

## Result

```
T45 tape40 differing records: 412 -> 6
```

The 6 remaining are exactly the pre-existing BROADR ULP cohort owned by
`NJOY_jl-bdu`: MT1 line 80, MT102 lines 2371/2563, MT113 line 3713, MT800
line 6015, and MT205 line 4419 (the triton channel is
`sigma(MT105)+sigma(MT113)`, so it inherits the MT113 ULP verbatim).

T45 moves `DIFFS → NUMERIC_PASS` (7,183/7,188 at 1e-5).

New test `test/validation/test_t45_gaspr_gas_yields.jl`: 15 unit assertions on
`gas_channel` hand-derived from the Fortran (including a non-B-10 case — Li-6
(n,t) leaves an alpha), plus tape-level assertions that MT203/204/206/207 are
byte-exact, that MT205's only differing record is the documented `bdu` one, and
that TPID line 1 is 66 blanks.

## Filed, not fixed

* **`egas(1)` must be the synthetic `thrg`**, not the first MT1 knot ≥ `thrg`
  (gaspr.f90:449-458 sets `en=enext=thrg` on the first iteration). Julia filters
  `e >= thrg` instead. Latent: B-10's `thrg` is 1e-5, which *is* MT1's first
  point, so T45 is unaffected.
* **`NJOY_jl-jbi`** — MF6/MT5 energy-dependent yields (the `y=111` sentinel,
  gaspr.f90:501-506/839-868) remain unported. Pre-existing, not a regression.
* **`NJOY_jl-xpf` CONFIRMED, not stale.** `_collect_gaspr!` (pipeline.jl:471-480)
  has zero call sites — dead code — while broadr/heatr/thermr all call their
  `_collect_*!` sibling. When heatr or thermr populate `ctx.extra_mf3`,
  `final_assembly!` rebuilds the tape from the RunContext and never reads
  gaspr's output, dropping MT203-207. T13's tape28 shows NXC 30 vs 23,
  consistent with this. Left for its own oracle-gated change.

## Method note

The GASPR mechanism was established from the oracle tape *first* (pure
arithmetic on the reference MF3 sections) and only then confirmed against
`gaspr.f90`. Seven read-only research lanes ran in parallel, each followed by an
adversarial citation-checker; six came back SOUND and one **PARTIALLY_REFUTED** —
the ENDF-inventory lane had reflowed its "verbatim" ENDF quotes into a
prettified form that appears nowhere in the file. Its substantive counts held
(LR=22 ×14, LR=35 ×14, matching an independent parse), but the refutation is a
concrete instance of Law 2's "never trust a subagent's paraphrase". No subagent
ran Julia (Rule 9); every Julia invocation this phase was serial on the main
thread.
