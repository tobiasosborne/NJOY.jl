# Phase 82 — T86 BIT_IDENTICAL (groupr ign=1 + card9a extended format + MF10 activation)

**Date:** 2026-06-02
**Beads:** `NJOY_jl-57z` (ign=1) + `NJOY_jl-ow6` (card9a/MF10) — both CLOSED
**Commits:** `675fdb5` (Cycle 1: parsing + panel quadrature), `2375874` (Cycle 2: reconr MF10 write + groupr MF10 averaging)
**Outcome:** **T86 tape25 STRUCTURAL_FAIL 2/12 → BIT_IDENTICAL 52/52.** A NEW full bit-identical test (Hf-177 TENDL2021, MAT=7234). Regression-clean.

Orchestrated: 2 Opus coding cycles (serial Julia), preceded by an oracle-alignment foundation step + a 4-agent read-only triage. Orchestrator verified every result independently (Rule 2).

---

## Foundation — oracle was stale (had to align to the pin first)
The local `njoy-reference/` clone was at `3c9873d` (2025-12-15 laptop oracle), NOT the pin `2c64dfb` (2026-03-18). Consequences caught this session:
- T22's local `referenceTape20` EMAX = `0.000000+0` (old), but the committed leapr fix emits the nonzero pin value → **T22 was actually DIFFS against the stale local oracle** despite the HANDOFF claiming BI. Aligning to the pin recovered T22 → BIT_IDENTICAL 4636/4636.
- tests/86 + tests/87 were **not checked out** at 3c9873d, so T86/T87 literally could not run.
The pin was already present in the local clone; `git -C njoy-reference checkout 2c64dfb` is a safe forward move (local HEAD is an ancestor of the pin). Only 09/22/23/33/80 references change (all the EMAX fix Julia already emits) + 86/87 added; every other test's references are byte-identical between the two commits. `njoy-reference/` is git-ignored → per-machine step, no NJOY.jl commit.

## Triage (4 read-only Sonnet agents) — all candidates PARTIAL; T86 chosen as the tractable milestone
- `3rp` (T14 ACER suffix `80c→00c`): clean 2-line fix but T14 needs ~300-500 LOC more (MF6 LAW=1 Kalbach, photon prod). Not a milestone.
- `65i` (T87 LUNR): only ~3% of T87's gap (needs photon prod + secondary-neutron dist + purr fidelity). Deep.
- `c3q` (T70 lat=10): **bead premise stale** — Julia 213800 == local ref 213800 (not over-production); real bug is a corrupt MF1/MT451 header. Re-scoped on the bead.
- `57z` (T86): ~140 LOC, no algorithmic unknowns, small (52-line) target → chosen.

## Cycle 1 (commit 675fdb5) — parsing + panel quadrature
- `parse_groupr`: read card6a `ngn` + card6b `ngn+1` boundaries when `abs(ign)==1` (gengpn, groupr.f90:4156-4165) into new `GrouprParams.user_egn`; parse card9a (`mfd==-1` → next card `file sec izar lfs`) and direct compound mfd (≥1e7) into `GrouprActivationEntry` (groupr.f90:625-713). `_groupr_group_structure` ign==1 branch returns user_egn.
- **Surprise (honest agent disclosure):** the brief assumed iwt=3 averaging "already worked." It did not — the pre-existing `group_integrate` linearizes 1/E over the PENDF grid and gave g1 flux ~10× low (2.46 vs 25.61). It had **never been validated for bit-identical GENDF** (other groupr tests check downstream errorr tapes, not the GENDF directly). Required a faithful panel Lobatto-2 quadrature port (`_groupr_getwtf`/`_groupr_getsig`/`_groupr_panel_xs`, groupr.f90:5115-6091) for nz=nl=1/MF3/iwt∈{2,3}/lord=0. `getwtf` `sigfig(enext,6,1)` step boundary (5305) fixed the last-ULP g2 flux.
- Writer: per-record ng2=2 for normal MTs (3 only for 251/252/253/452/455/456); group constants via `format_endf_float(...; extended=false)` (7-sigfig a-format, matches Fortran listio). Bonus: this fixed T15's tape91 MF3 group constants (were mis-formatted as `11.5383566`).
- Result: T86 tape25 lines 1-16 (dict + MF3/MT4 C2=0, all 4 groups) **bit-identical**; first-diff → line 17. Regression-clean (T01/T04/T11 unchanged).

## Cycle 2 (commit 2375874) — MF10 data path
- **reconr** (`reconr.jl`/`reconr_types.jl`/`pendf_writer.jl`/`pendf_io.jl`): new `_reconstruct_mf10` reconstructs every MF10 sub-section onto the union grid (INT=2 lin interp, `sigfig(·,7,0)`, gety1 `thresh` + `ith` trim, reconr.f90:4734/4792/4834/4918) and writes them to the PENDF with an updated MF1/MT451 directory. **Gated on MF10 presence** — returns empty for materials without MF10, so the path is a no-op for every non-activation evaluation (T01/T02/T08 untouched).
- **groupr** (`groupr.jl`): per `GrouprActivationEntry` with file index 4 (MF10) → `extract_mf10_subsections` match by target ZA (izar) + final state (lfs), group-average with the EXACT iwt=3 panel quadrature, emit MF3/MTd with the entry's `fzam` in C2. `igzero` zero-group skip (groupr.f90:893) drops below-threshold groups so the lfs=29/30 excited-state sections correctly emit only groups 2-4 / 3-4.
- Result: **T86 BIT_IDENTICAL 52/52** — all 5 sections (MF3/MT4 C2=0 + four MF10 activation sections, izam C2 = 72177.0 / 721770.0 / 72177.029 / 72177.03).

## Regression (orchestrator-verified, post-Cycle-2)
T01 NUMERIC_PASS 32812/32962; T02 tape28 NUMERIC_PASS 12519/13873; T08 STRUCTURAL_FAIL 467/14680 + 3/16878; T04 81/82 + 56/74 + 107/119 — **all exactly at baseline**.

## Lessons
- **A "reuse the existing path" premise is a trap.** The iwt=3 group averaging was assumed working; it was 10× wrong at g1 and had simply never been checked against a bit-identical GENDF. The agent's reproduce step exposed it. Always reproduce, never assume a path is validated.
- **Validate the oracle before validating the code.** The local oracle was 3 months stale; T22's "BI" was only true against the pin. One `checkout` recovered a test and unblocked two more. Check `reference oracle at pin` in the preflight.
- A subagent rate-limited mid-finish (Cycle 2) left a complete, applied diff; the orchestrator verified it end-to-end (BI + regressions + no hardcoding + gate) without the agent's report — the run is the source of truth, not the report.
