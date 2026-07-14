# NJOY.jl

A Julia reimplementation of [NJOY2016](https://github.com/njoy/NJOY2016) — the
standard nuclear data processing system used worldwide for reactor physics,
criticality safety, and radiation transport (ENDF → PENDF / ACE / GENDF /
covariance). The Fortran reference is ~120 k lines of Fortran 90 in
`njoy-reference/src/` (39 files, ~100 k code-only); this port is ~33 k lines of
idiomatic Julia in `src/` (107 files, code-only per `cloc`).

## Project goal

Reproduce NJOY2016 faithfully in Julia: the same tape-driven architecture, the
same sections and grids, and ultimately the same bytes on NJOY's own reference
suite. This is not a new nuclear-data system and does not reinterpret or
"improve" NJOY's physics. The implementation may be idiomatic Julia; its
observable behaviour is frozen by the canonical Fortran.

The Fortran is canonical truth. Where NJOY rounds to 8 sigfigs before adding
potential scattering, or accumulates in a particular IEEE-754 order, or carries
a quirk, we reproduce it exactly — expressed idiomatically (multiple dispatch,
parametric types, scoped state), never blindly transliterated. See `AGENTS.md`
or `CLAUDE.md` for the working agreement and `HANDOFF.md` for the living state.

## Status at a glance

*Verified against the frozen NJOY2016 oracle at `2c64dfb3` (see
`REFERENCE_PIN`) via the committed Phase-91 86-test sweep; 0 crashes and 0
timeouts.*

Sweep summary: **9 `BIT_IDENTICAL`, 4 `NUMERIC_PASS`, 72 `DIFFS`, 1
`NO_REFERENCE`**. There are no structural failures or missing tapes.

Full test bit-identical (`rtol = 1e-9`, **every produced tape**):

| Test | Chain | Tapes | Result | Note |
|------|-------|-------|--------|------|
| T03 | moder → reconr | tape37 | 9274 / 9274 | photoatomic |
| T09 | leapr | tape24 | 1830 / 1830 | discrete oscillators |
| T22 | leapr | tape20 | 4636 / 4636 | S(α,β) |
| T50 | acer (α + He-4) | tape33/34/35 | 432 + 163 + 1 | aplots plot tape ported (Phase 79) |
| T52 | acer (p + H-1) | tape33/34/35 | 6042 + 3986 + 1 | aplots ported (Phase 79) |
| T53 | acer (d + H-2) | tape33/34/35 | 17236 + 12030 + 1 | aplots heating + aploxp ported (Phase 80) |
| T61 | acer (iopt=7 thermal) | tape71/72 | 49211 + 1 | thermal ACE passthrough |
| T62 | acer (d + He-3) | tape33/34/35 | 5460 + 7221 + 1 | aplots aplotr threshold ported (Phase 80) |
| T86 | groupr | tape25 | 52 / 52 | user groups + MF10 activation |

`NUMERIC_PASS` means every produced tape passes at least one tolerance in the
sweep's configured ladder (`1e-7`, `1e-5`, `1e-3` after the strict `1e-9`
check). It is not synonymous with the project's `1e-7` first-round target.

| Test | Chain | Loosest tolerance required |
|------|-------|---------------------------:|
| T01 | reconr → broadr → heatr → thermr → groupr | `1e-5` |
| T33 | leapr | `1e-5` |
| T54 | charged-particle acer + aplots | `1e-3` |
| T80 | leapr | `1e-5` |

Additional ground covered:

- **RECONR is the most mature module** — ~19 reference materials are
  bit-identical on the RECONR-produced PENDF across every resonance formalism
  (LRU=0, SLBW, MLBW, Reich-Moore, SAMMY/RML, URR). The *full pipelines* built
  on top of those then diverge in the downstream module being worked on; see the
  per-test table in `reports/REFERENCE_SWEEP.md`.
- The remaining tests run end-to-end and land in `DIFFS` (last-digit / FP grind,
  or a not-yet-ported sub-path) — these are the active grind, not breakage.

## Module maturity

The 23 NJOY modules, as the **code** stands today (not aspirations):

| Maturity | Modules |
|----------|---------|
| **Full** (real port, no stub paths in the common case) | `reconr`, `broadr`, `unresr`, `heatr`¹, `gaspr`, `covr`, `mixr`, `resxsr`, `matxsr`, `dtfr`, `viewr` |
| **Partial** (core paths work; some options `@warn`/stubbed) | `thermr`², `leapr`³, `groupr`⁴, `errorr`⁵, `gaminr`, `acer`⁶, `powr`⁷, `moder`⁸ |
| **Stub** (structure only / not yet wired) | `purr`⁹, `wimsr`¹⁰, `ccccr`¹¹, `plotr`¹² |

¹ heatr KERMA/damage full; photon plot tape stubbed.
² thermr free-gas + S(α,β) + coherent and incoherent elastic (LTHR=1/2/3);
  small MF3/MT222 interpolation residuals and downstream thermal ACER remain.
³ leapr phonon/translational/cold-H/Sköld + discrete oscillators; coherent
  elastic (`iel>0`) and the `b7<=0` secondary-scatterer path remain open.
⁴ groupr MF3 + nubar + multi-T/multi-σ₀; MF6/8/10/16+ transfer matrices and
  MT=251/252 derivation incomplete.
⁵ errorr MF31/33 covariance + LB=0/1/2/5/6 + NC + rescon (LCOMP=1); LTY=1/2/3
  standards, MF34 multi-section, LCOMP=0/2 still open.
⁶ acer iopt=1 fast-neutron + charged-particle paths + the tested iopt=7 path;
  iopt=2/3/4/5 are stubs. Aplots emits the Phase-79/80-tested pages and fails
  loudly on remaining unported blocks.
⁷ powr lib=1 (GAMTAP) full; lib=2 (EPRI-CPM) / lib=3 error out.
⁸ moder tape copy/convert; material-filter extraction + multi-tape merge missing.
⁹ purr emits correct tape *structure* with placeholder (zero) probability tables.
¹⁰ wimsr parses decks + writes a structurally-valid but zero-content WIMS tape.
¹¹ ccccr writes ISOTXS/BRKOXS/DLAYXS marker files only — no binary content.
¹² plotr touches an empty output tape.

## Repository layout

```
src/
  endf/           ENDF-6 tape reader/writer, record types, MF/MT registry
  resonances/     SLBW, MLBW, Reich-Moore, SAMMY/RML, URR reconstruction
  processing/     reconr, broadr, heatr, thermr, unresr, errorr, ... physics
  formats/        ACE, GENDF, MATXS, CCCC, WIMS, DTF, POWR writers
  orchestration/  run_njoy, input-deck parser, module wrappers, tape manager
  viewr/          viewr PostScript engine (graph.f90 port)
  visualization/  plotr/visualization backends
  NJOY.jl         top-level module
test/
  validation/     reference_test.jl (execute.py port), sweep, diagnose harness
docs/             design.md, architecture.md, api.md, tutorial.md
njoy-reference/   NJOY2016 Fortran source + tests (git-ignored; clone yourself)
reports/          acceptance criteria, sweep results, review waves
worklog/          per-session debug journals (T*.md, phase*.md)
```

## Quick start

The Fortran reference tree is **git-ignored** and must be cloned locally — it
provides both the canonical Fortran source and the reference test tapes the
sweep compares against.

```bash
# 1. Clone/fetch NJOY2016 and check out the frozen commit from REFERENCE_PIN
bash scripts/setup_reference.sh

# 2. Build the Fortran oracle (optional — only to regenerate reference tapes)
cmake -S njoy-reference -B njoy-reference/build
cmake --build njoy-reference/build --target njoy

# 3. Run the Julia tests against the reference tapes.
#    ALWAYS clear the precompile cache first (concurrent Julia corrupts it):
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/reference_test.jl 1        # single test
julia --project=. test/validation/sweep_reference_tests.jl   # full sweep (~90 min)
```

Issue tracking uses **Beads** in embedded Dolt mode. Durable cross-machine state
is exported to the tracked `.beads/issues.jsonl`; there is no Dolt remote. Run
`bd prime`, then `bd blocked` and `bd ready`, before selecting work.

## Acceptance tolerances

See `reports/ACCEPTANCE_CRITERIA.md` for the full hierarchy with maintainer
quotes. In brief:

| Tolerance | Meaning |
|-----------|---------|
| **1e-9** (bit-identical) | byte-for-byte after date wildcarding — the stretch goal |
| **1e-7** (numeric pass) | ±1 in the 7th sigfig — the cross-compiler floor (even Fortran↔Fortran fails 1e-9 across architectures) |
| **structural** | same line counts, same sections, same grid sizes — non-negotiable at any tolerance |

## Contributing

Read **`CLAUDE.md`** (the working agreement — the two laws: oracle-driven TDD,
and Fortran-before-Julia) and **`HANDOFF.md`** (the living state: current phase,
open work, traps) before touching code. License: GPL-3.0.
