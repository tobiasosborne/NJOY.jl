# NJOY.jl Session Handoff

> Live state only. The complete pre-cleanup Phase 91 handoff is preserved
> verbatim in `handoff_archive/handoff_phase91_pre_hygiene.md`; older history is
> indexed by `handoff_archive/INDEX.md`.

## What this project is

NJOY.jl is a faithful, tape-driven Julia replacement for NJOY2016. It covers
the 23 NJOY modules exercised by the reference suite and targets the same
PENDF, ACE, GENDF, covariance, and auxiliary outputs as the canonical Fortran.

The scope is deliberately frozen: no new physics, no reinterpretation, and no
"cleaning up" NJOY quirks that changes output. Julia may express the algorithms
idiomatically, but the Fortran controls constants, branch decisions, rounding,
accumulation order, sentinels, records, and tape plumbing.

Primary references:

- `AGENTS.md` / `CLAUDE.md` — governing workflow and repository rules.
- `njoy-reference/src/*.f90` — canonical implementation truth.
- `REFERENCE_PIN` — frozen NJOY2016 oracle revision.
- `reports/REFERENCE_SWEEP.md` — latest full-suite result.
- `worklog/INDEX.md` — chronological phase/session record.
- `reports/ACCEPTANCE_CRITERIA.md` — numeric and structural acceptance rules.

## The two laws

1. **Oracle-driven TDD.** Every non-trivial fix begins with a failing comparison
   against the Fortran oracle or a reference tape. Name the exact tape and line
   range in the regression assertion before changing implementation code.
2. **Fortran before Julia.** Read the exact Fortran routine and relevant ENDF-6
   section before writing processing, evaluator, reader, or writer code. Cite
   the routine and line range in the Julia implementation.

Structural agreement is mandatory at every tolerance. Bit-identical output at
`1e-9` is the stretch target; `1e-7` is the first-round last-digit floor. Never
accept an apparent floating-point floor without tracing it: every supposedly
irreducible floor investigated so far has exposed a real bug.

## Reference oracle pin

The local `njoy-reference/` tree is git-ignored and must be aligned on every
machine:

```bash
bash scripts/setup_reference.sh
```

Canonical pin:

- SHA: `2c64dfb3d7cd57bea9fa8e36b995eee4a9d88d58`
- Upstream branch at pin time: `develop`
- Commit date: 2026-03-18
- Pin recorded: 2026-06-01

The old `ac5adf5` baseline is a rebased-away, unreachable commit and must not be
used. `reference_test.jl` and `sweep_reference_tests.jl` run a fail-soft pin
preflight; fix any warning before trusting a comparison.

## Current state — Phase 93 (2026-07-15)

The post-Phase-93 full sweep covers all 86 reference tests at the pinned
Fortran oracle and the documented 300-second per-test timeout. It completed in
56.0 minutes:

| Status | Count |
|---|---:|
| `BIT_IDENTICAL` | 11 |
| `NUMERIC_PASS` | 5 |
| `DIFFS` | 68 |
| `STRUCTURAL_FAIL` | 0 |
| `MISSING_TAPE` | 0 |
| `NO_REFERENCE` | 1 |
| `CRASH` | 0 |
| `TIMEOUT` | 1 |

Bit-identical full tests are T03, T09, T22, T50, T52, T53, T61, T62,
**T81, T84**, and T86. Numeric-pass tests are T01, T33, T54, T80, and
**T85**; inspect the per-tape tolerance columns because T85 currently passes
at `1e-5`, not the first-round `1e-7` floor.

Phase 93 closed two direct-RECONR defects:

- `NJOY_jl-gkw`: direct PENDF output now carries the full source-derived
  ENDF-6 MF1 header, canonical blank description, per-section directory NC,
  source ZA/AWR, faithful redundant L2/source LR, and signed-`npend` sequence
  behavior. T81 moved from 2/45,369 aligned records to **45,372/45,372 BI**;
  T84 moved from 2/1,112 to **1,115/1,115 BI**.
- `NJOY_jl-rhm`: unresolved MT152 totals remain independent after background
  addition and the sixth payload value copies total rather than elastic.
  T85 is structurally exact and now NUMERIC_PASS with seven residual lines.

T83's remaining `csunr2`/downstream-grid cluster is `NJOY_jl-fod`. T85's
six-node last-digit grind is `NJOY_jl-c33`. T17 remains the sole timeout and
again reached ERRORR before the 300-second cap. See
`worklog/phase93_reconr_direct_pendf.md`.

### Prior Phase 92 context

The latest full sweep is the clean post-Phase-92 run of all 86 reference tests
at `a7b9335`, using the documented default 300-second per-test timeout. It
completed in 55.2 minutes:

| Status | Count |
|---|---:|
| `BIT_IDENTICAL` | 9 |
| `NUMERIC_PASS` | 4 |
| `DIFFS` | 71 |
| `STRUCTURAL_FAIL` | 0 |
| `MISSING_TAPE` | 0 |
| `NO_REFERENCE` | 1 |
| `CRASH` | 0 |
| `TIMEOUT` | 1 |

Bit-identical full tests: T03, T09, T22, T50, T52, T53, T61, T62, and T86.
Numeric-pass tests under the sweep's configured tolerance ladder: T01, T33,
T54, and T80. Read the per-tape tolerance columns in
`reports/REFERENCE_SWEEP.md`; `NUMERIC_PASS` does not by itself mean `1e-7`.

Phase 92 closed two orchestration milestones:

- Live documentation and Beads state were reconciled. The pre-cleanup HANDOFF
  and design document are preserved verbatim in their archive directories;
  README, HANDOFF, architecture, design, tutorial, and API docs now describe
  the tape-first implementation and current oracle pin.
- `NJOY_jl-2k4`: RECONR now unionizes the first total MF12 `LO=1` and MF13
  TAB1 records even on the LRU=0 path, writes `NK=1` lin-lin sections, and
  BROADR preserves non-MF2/MF3 sections. T45 photon bodies match the oracle
  exactly in columns 1-66: MF12/102 NP=1014, MF13/4 NP=360, MF13/103 NP=228;
  directory NC values are 341/123/79 and MF14 remains intentionally absent.
  The focused regression passes 7/7, T84 preserves its nonzero MF12 TAB1
  metadata exactly, the negative-Q/eligibility unit regression passes 6/6,
  and T50/T52/T53/T61/T62 remain bit-identical.

The full T45 tape is still structurally short by three MF1 records
(7185 versus 7188); that separate metadata problem is `NJOY_jl-1kf`. The sweep
confirms the Phase-92 photon gain: T45 moved from 6634 to 7185 generated lines,
and T84 from 849 to 1112 against 1115 oracle lines. The 9 bit-identical and 4
numeric-pass tests were preserved.

T17 is the only timeout. It reached ERRORR and was interrupted at the next
yield after the 300-second limit (311.8 seconds observed). Its partial tape
comparisons are unchanged from Phase 91, when the same deck completed in 612.9
seconds under a longer allowance, so this is a default-sweep performance gate
rather than a fidelity regression. Track it as `NJOY_jl-5tu`.

Phase 91 landed three milestones:

- `NJOY_jl-czw` (`f3c68c4`): removed dead `thermr_coh_ne` plumbing.
- `NJOY_jl-4k2` (`895795e`): changed the resolved-resonance channel-radius
  exponent from `1/3` to Fortran's truncated `0.333333333`. T49's MLBW elastic
  value at 110487.7 eV became bit-identical.
- `NJOY_jl-k1v` (`abc1e93`): lunion now shades only breakpoints whose local
  interpolation law is histogram (`inta==1`). T42 MT1 shrank from 52,307 to the
  exact Fortran grid of 51,211 points; tape44 dropped 44,372 spurious lines.

All 9 bit-identical and 4 numeric-pass tests were preserved. The six prior
contention-related timeouts completed as ordinary `DIFFS`; there were no status
regressions. See `worklog/phase91_reconr_third_lunion_shading.md`.

## Beads status and live work

Beads runs in embedded Dolt mode. `.beads/issues.jsonl` is the git-tracked,
cross-machine record; there is no Dolt remote and no `bd dolt push` target.

```bash
bd prime
bd blocked
bd ready
bd show <id>
bd update <id> --claim
```

Every non-trivial change starts with a bead. Some broad umbrella issues predate
later phase closures, so confirm an issue's current premise against the latest
worklog, Fortran, and Julia before claiming it. Do not revive a completed task
merely because an old HANDOFF snapshot calls it open.

Immediate Phase 93 follow-ups:

- **`NJOY_jl-5tu` (open):** make T17 complete under the documented default
  300-second sweep limit; the current bottleneck is ERRORR, not BROADR.
- **`NJOY_jl-fod` (open):** match T83's `csunr2` accumulation/state-machine
  and downstream grid propagation. Direct serialization and MT152 packing are
  already faithful; start at tape50 line 103.
- **`NJOY_jl-c33` (open):** eliminate T85's seven-line MT152 numerical
  residual and reach at least the `1e-7` first-round floor.
- **`NJOY_jl-1kf` (open):** reproduce T45's final TPID/MF1/MT451 ownership and
  three missing metadata records. Photon sections and their directory tuples
  are already exact; keep this follow-up isolated from `NJOY_jl-2k4`.
- **`NJOY_jl-6lg` (open):** BROADR's new passthrough directory entries preserve
  section bodies but currently default MOD to zero; retain incoming nonzero MOD
  values from MF1/MT451 in a separately oracle-driven change.
- **`NJOY_jl-c9f` (open, researched):** trace the 4–5 interior MF3/MT222
  incoherent-elastic values off by at most `4e-7` in T69/T74. The three-point
  Lagrange loop matches Fortran; instrument T69 first and inspect upstream
  Debye-Waller/`exp` intermediates before blaming libm.
- **T49 metadata:** the residual HEAD-record `0 0` versus `0 99` difference is
  separate from the resolved `4k2` channel-radius bug.

Larger active families remain tracked in Beads: ACER distributions and modes,
ERRORR covariance fidelity, GROUPR derived reactions, LEAPR secondary/coherent
paths, PURR probability tables, WIMSR, CCCCR, and PLOTR.

## Recent phases

| Phase | Date | Outcome | Worklog |
|---:|---|---|---|
| 93 | 2026-07-15 | T81/T84 bit-identical; T85 numeric; direct RECONR serialization and MT152 layout fixed | `phase93_reconr_direct_pendf.md` |
| 92 | 2026-07-14 | Stale docs/Beads reconciled; T45 MF12/MF13 totals unionized and preserved exactly | `phase92_staleness_t45_photons.md` |
| 91 | 2026-06-26 | `_THIRD` and mixed-law lunion shading fixed; dead THERMR field removed; sweep refreshed | `phase91_reconr_third_lunion_shading.md` |
| 90 | 2026-06-20 | BROADR thermal MT1 rebuilt from widened, deduplicated broadened partials | `phase90_e5n_broadr_mt1_partials.md` |
| 89 | 2026-06-19 | Dense Bragg grid fixed and accelerated; free-gas cutoff corrected | `phase89_thermr_grid_freegas_e5n.md` |
| 88b | 2026-06-16 | THERMR LTHR=3 mixed coherent/incoherent elastic wired | `phase88b_iel_lthr3_mixed.md` |
| 88 | 2026-06-16 | THERMR LTHR=2 incoherent elastic wired end-to-end | `phase88_iel_incoherent_elastic_lthr2.md` |
| 87 | 2026-06-15 | THERMR coherent cutoff held-point sentinel restored | `phase87_thermr_emax_sentinel.md` |
| 86 | 2026-06-14 | THERMR MF3 directory-NC quirk reproduced | `phase86_h61_thermr_mf3_dir_nc.md` |
| 85 | 2026-06-14 | `calcem` convergence now uses Fortran's `terp1` chord | `phase85_thermr_calcem_terp1_convergence.md` |
| 84 | 2026-06-09 | THERMR MF6 cosine width corrected; oracle setup made direct-SHA safe | `phase84_thermr_mf6_cosine_width.md` |
| 83 | 2026-06-04 | T09 BI, T33 numeric, T72 crash fixed, T70 header fixed | `phase83_o2l_lho_orchestrated.md` |

All paths are relative to `worklog/`. Older phases are indexed in
`worklog/INDEX.md`; historical HANDOFF prose is under `handoff_archive/`.

## Architecture

NJOY.jl follows NJOY's tape architecture:

```text
input deck
   ↓
run_njoy(input; work_dir)
   ↓
module runner reads input tape(s)
   ↓
module computes and writes output tape(s)
   ↓
next module re-reads those tapes
```

Do not add shared mutable cross-module state or pass in-memory processed data
around the tape boundary. `RunContext` contains transitional assembly plumbing;
new work should reduce rather than expand that exception.

Key source areas:

| Path | Responsibility |
|---|---|
| `src/orchestration/pipeline.jl` | deck execution, tape context, module dispatch |
| `src/orchestration/modules/` | per-module tape runners |
| `src/processing/` | processing algorithms and PENDF assembly |
| `src/resonances/` | SLBW, MLBW, Reich-Moore, SAMMY/RML, URR |
| `src/endf/` | ENDF record parsing/writing and section readers |
| `src/formats/` | ACE, GENDF, MATXS, WIMS, CCCC, DTF, POWR formats |
| `test/validation/` | reference runner, sweep, oracle diagnostics |

The Fortran module files in `njoy-reference/src/` are the architecture map:
`reconr.f90`, `broadr.f90`, `heatr.f90`, `thermr.f90`, `groupr.f90`,
`errorr.f90`, `acer.f90`/`acefc.f90`, and the remaining module peers.

## Critical traps

1. Clear `~/.julia/compiled/v1.12/NJOY*` before every Julia run. Never run two
   Julia processes concurrently.
2. `round_sigfig` includes Fortran's bias. Use `_dedup_tol!`, not `unique!`, for
   energy grids.
3. Preserve Fortran constants exactly: `ehigh=20e6`, `elow=1e-5`,
   `elim=min(0.99e6,eresr)`, and `third=0.333333333`.
4. MF3 section AWR may differ from MF2 AWR; thresholds use the section AWR.
5. Multi-material reads must pass `target_mat`.
6. Physical constants intentionally use the CGS values in `phys.f90`.
7. MF2 can contain multiple resolved and unresolved ranges.
8. MF1/MT455 with `LNU=2` has a LIST preamble before TAB1.
9. ERRORR `ign=-1` uses the union grid; `ign=1` uses `user_egn`; `ign>=2` uses
   the library structure.
10. `emerge` reads MF3 from lunion's scratch tape, not directly from ENDF.
11. Initial duplicate shading uses `sigfig(e,7,0)`; mid-data shading uses
    `sigfig(e,7,-1)`.
12. MT18 is redundant when MT19 exists (`mtr18`); MT103–107 are redundant only
    when their partial ranges exist.
13. Reich-Moore inversion must use `_frobns`/`_thrinv!`/`_abcmat`, not `inv`.
14. The reference-test float regex can misread dense fixed-column ENDF lines;
    use a six-by-eleven-character parse when a structural claim looks absurd.
15. Read `err` from each reference input deck; historical tables have been
    wrong.
16. Positive RECONR `npend` carries NS across copied sections/files, but a
    computed redundant section still writes `SEND=99999` and restarts NS.
17. MF2/MT152 stores total independently from elastic+fission+capture; its
    sixth per-energy value is a copy of total.

## Test workflow

One Julia process only:

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/reference_test.jl 7
```

Full sweep, at session end or in the sole Julia slot:

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/sweep_reference_tests.jl
```

Generate a focused Fortran oracle:

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/diagnose_harness.jl <N>
```

Grind method:

1. Reproduce the failing reference tape and exact line/section.
2. Compare structural slices first: section presence, directory counts, line
   counts, grids, then values. For MF3, start with columns 1–66.
3. Classify the first difference: grid, threshold, missing feature, formatting,
   or value accumulation.
4. Read the exact Fortran routine and ENDF section.
5. Add the failing assertion, implement a small fix, and rerun the driver test.
6. Rerun the relevant BI canaries and end with a passing reference result.

Canonical regression set for resonance-sensitive changes:
T01, T02, T08, T27, T34, T45, and T46. Add the module-specific BI cohort.

## T01 and diagnostic how-tos

Full T01 reference chain:

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/reference_test.jl 1
```

Direct pipeline execution:

```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/01/input"; work_dir="/tmp/t01")'
```

For a deep FP trace:

1. Add tightly filtered `write(*,*)` diagnostics to the exact Fortran routine.
2. Rebuild: `cmake --build njoy-reference/build --target njoy`.
3. Run the smallest truncated oracle deck and capture the tagged output.
4. Print the same Julia intermediates in the same loop order.
5. Restore the Fortran checkout and rebuild before running the oracle again.

Never trust a diagnostic built from a dirty Fortran tree. Do not use formatted
Fortran diagnostics when exact intermediate precision matters; list-directed
`write(*,*)` has repeatedly exposed the first divergent value.

Useful T01 files:

| File | Focus |
|---|---|
| `src/processing/broadr.jl` | Doppler broadening and partial selection |
| `src/processing/sigma1.jl` | Sigma1 Doppler kernel |
| `src/processing/heatr.jl` | KERMA and damage |
| `src/processing/thermr.jl` | thermal kernels, grids, MF7 readers |
| `src/processing/pendf_writer.jl` | PENDF records and directory metadata |
| `src/processing/reconr.jl` | resonance reconstruction orchestration |
| `test/validation/reference_test.jl` | end-to-end reference comparator |

## Session close

Follow the complete protocol in `AGENTS.md`: inspect status, stage specific
files, commit, push, update the worklog/HANDOFF if state changed, and record a
Beads memory for genuinely surprising cross-session lessons. Work is not
complete until `git push` succeeds.
