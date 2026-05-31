# NJOY.jl

A Julia reimplementation of [NJOY2016](https://github.com/njoy/NJOY2016) — the
standard nuclear data processing system used worldwide for reactor physics,
criticality safety, and radiation transport (ENDF → PENDF / ACE / GENDF /
covariance). The Fortran reference is ~120 k lines of Fortran 90 in
`njoy-reference/src/` (39 files, ~100 k code-only); this port is ~33 k lines of
idiomatic Julia in `src/` (105 files, code-only per `cloc`).

## Project goal (north star)

A full Julia reimplementation of NJOY that is **bit-identical** with the Fortran
reference, and that is also **maintainable, extensible, faster, and more
pleasant to work with** than the original.

We get there in two phases, in this order:

1. **Mini north star — bit-identical first.** Reproduce NJOY2016's output
   byte-for-byte (or within the published cross-compiler tolerance) on NJOY's
   own reference test suite. This is the trust we build everything else on.
2. **North star — then make it better.** Once an area is bit-identical and we
   have confidence in it, pursue the maintainability / extensibility / speed /
   ergonomics gains that a modern Julia design allows — *without* losing the
   bit-identical guarantee.

The Fortran is canonical truth. Where NJOY rounds to 8 sigfigs before adding
potential scattering, or accumulates in a particular IEEE-754 order, or carries
a quirk, we reproduce it exactly — expressed idiomatically (multiple dispatch,
parametric types, no globals), never transliterated. See `CLAUDE.md` for the
working agreement and `HANDOFF.md` for the living state.

## Status at a glance

*Verified against a fresh clone of NJOY2016 `develop @ ac5adf5` (2026-04-06) via
the full 86-test sweep (`reports/REFERENCE_SWEEP.md`); 0 crashes, 0 timeouts.*

Full test bit-identical (`rtol = 1e-9`, **every produced tape**):

| Test | Chain | Tapes | Result | Note |
|------|-------|-------|--------|------|
| T03 | moder → reconr | tape37 | 9274 / 9274 | photoatomic |
| T61 | acer (iopt=7 thermal) | tape71/72 | 49211 / 49211 | thermal ACE passthrough |
| T50 | acer (α + He-4) | tape33/34/35 | 432 + 163 / all | aplots plot tape ported (Phase 79) |
| T52 | acer (p + H-1) | tape33/34/35 | 6042 + 3986 / all | aplots ported (Phase 79) |
| T53 | acer (d + H-2) | tape33/34/35 | 17236 + 12030 / all | aplots heating + aploxp ported (Phase 80) |
| T62 | acer (d + He-3) | tape33/34/35 | 5460 + 7221 / all | aplots aplotr threshold ported (Phase 80) |

Numeric pass (`rtol = 1e-5`, the cross-compiler floor):

| Test | Chain | Tape | Result |
|------|-------|------|--------|
| T01 | reconr → broadr → heatr → thermr → groupr | tape25 | 32812 / 32962 |
| T80 | leapr (S(α,β)) | tape24 | 91405 / 91453 |

Additional ground covered:

- **RECONR is the most mature module** — ~19 reference materials are
  bit-identical on the RECONR-produced PENDF across every resonance formalism
  (LRU=0, SLBW, MLBW, Reich-Moore, SAMMY/RML, URR). The *full pipelines* built
  on top of those then diverge in the downstream module being worked on; see the
  per-test table in `reports/REFERENCE_SWEEP.md`.
- The remaining tests run end-to-end and land in `DIFFS` (last-digit / FP grind,
  or a not-yet-ported sub-path) — these are the active grind, not breakage.

> **Note: the committed sweep is stale.** `reports/REFERENCE_SWEEP.md` predates
> Phases 78–80 and still shows only 2 `BIT_IDENTICAL` and the four ACER tests as
> `DIFFS`. The aplots plot tape (`tape33`) has since been ported — T50/T52
> (Phase 79) and T53/T62 (Phase 80) are now byte-perfect on **all** their tapes
> (verified via targeted cache-nuked runs). The classifier marks a *test*
> bit-identical only when *every* reference tape passes, so always read the
> per-tape detail; a fresh full sweep will refresh the headline count.

> **Reference drift.** This port is validated against a specific NJOY2016
> baseline. The upstream `develop` branch moves; e.g. T22 (leapr) is currently
> 4635/4636 — a single MF1/MT451 header field the newer reference emits
> differently. When resuming on a new machine, expect to pin or reconcile the
> `njoy-reference` checkout. See the open issues.

## Module maturity

The 23 NJOY modules, as the **code** stands today (not aspirations):

| Maturity | Modules |
|----------|---------|
| **Full** (real port, no stub paths in the common case) | `reconr`, `broadr`, `unresr`, `heatr`¹, `gaspr`, `covr`, `mixr`, `resxsr`, `matxsr`, `dtfr`, `viewr` |
| **Partial** (core paths work; some options `@warn`/stubbed) | `thermr`², `leapr`³, `groupr`⁴, `errorr`⁵, `gaminr`, `acer`⁶, `powr`⁷, `moder`⁸ |
| **Stub** (structure only / not yet wired) | `purr`⁹, `wimsr`¹⁰, `ccccr`¹¹, `plotr`¹² |

¹ heatr KERMA/damage full; photon plot tape stubbed.
² thermr free-gas + S(α,β) + coherent-elastic; incoherent-elastic (LTHR=2) not
  yet wired; Debye-Waller defaults to graphite.
³ leapr phonon/translational/cold-H/Sköld full; coherent-elastic (iel>0) and
  secondary scatterer (nss>0) stubbed.
⁴ groupr MF3 + nubar + multi-T/multi-σ₀; MF6/8/10/16+ transfer matrices and
  MT=251/252 derivation incomplete.
⁵ errorr MF31/33 covariance + LB=0/1/2/5/6 + NC + rescon (LCOMP=1); LTY=1/2/3
  standards, MF34 multi-section, LCOMP=0/2 still open.
⁶ acer iopt=1 fast-neutron + full charged-particle path + iopt=7; iopt=2/3/4/5
  (thermal/dosimetry/photoatomic/photonuclear) are stubs; aplots tape is a stub.
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
# 1. Clone the NJOY2016 reference (Fortran source + reference tapes)
git clone https://github.com/njoy/NJOY2016.git njoy-reference
# (to reproduce the validated baseline exactly, pin to the commit the port
#  was last verified against — see the open issues / HANDOFF.md)

# 2. Build the Fortran oracle (optional — only to regenerate reference tapes)
cd njoy-reference && cmake -S . -B build && cmake --build build --target njoy
cd ..

# 3. Run the Julia tests against the reference tapes.
#    ALWAYS clear the precompile cache first (concurrent Julia corrupts it):
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. test/validation/reference_test.jl 1        # single test
julia --project=. test/validation/sweep_reference_tests.jl   # full sweep (~50 min)
```

Issue tracking uses **beads** (`bd`); the durable state lives in
`.beads/issues.jsonl` and is loaded with `bd import` on a fresh clone. Run
`bd ready` to see available work.

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
