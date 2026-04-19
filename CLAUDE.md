# CLAUDE.md — NJOY.jl

## What is this?
NJOY.jl is a **100% faithful drop-in Julia replacement** for NJOY2016 — the standard nuclear data processing system (ENDF → PENDF / ACE / GENDF / covariance). 23 Fortran modules, 84 reference tests, ~15k lines of Julia replacing ~120k lines of Fortran 90. No GUI. No new features. No reinterpretation. Only bit-for-bit reproduction of the canonical output, expressed in idiomatic Julia.

See `HANDOFF.md` for the living state (current phase, sweep results, per-test status).
See `worklog/T*.md` for per-session debug journals.
See `reports/ACCEPTANCE_CRITERIA.md` for the tolerance hierarchy (1e-9 bit-identical → 1e-7 first-round → structural match non-negotiable).
See `njoy-reference/src/*.f90` for **canonical truth**.

---

## THE TWO LAWS (non-negotiable, read first, every session)

**LAW 1 — ORACLE-DRIVEN TDD.** Every non-trivial change starts with a *failing* oracle comparison. Generate the Fortran oracle (`diagnose_harness.jl` or the binary in `njoy-reference/build/njoy`), run Julia, diff MF3 cols 1-66 (or the equivalent structural slice), confirm the diff reproduces the bug, then fix. Bug fixes begin with a failing reference-test assertion that names the exact tape and line range. No code before the red bar. "Obvious" fixes are where oracles find surprises ("nothing is irreducible" — every claimed FP precision floor in this repo has turned out to be a real bug).

**LAW 2 — FORTRAN BEFORE JULIA.** Before writing any processing, evaluator, reader, or writer code, open and read the ground-truth source: the exact Fortran subroutine in `njoy-reference/src/*.f90`. Ground truth = (a) the specific routine, (b) the specific line range, (c) the ENDF-6 manual section when applicable. Never code from memory. Never code from the HANDOFF alone (it has been wrong multiple times). Never trust a subagent's paraphrase. Cite the source in code:
```julia
# Ref: njoy-reference/src/reconr.f90:2845-2847 (csslbw)
# sigfig(sigp(2),8,0) and sigfig(spot,8,0) BEFORE adding potential scattering
```

---

## THE RULES (numbered, follow to the letter)

0. **LAWS 1 & 2 APPLY.** Oracle-driven TDD. Fortran before Julia.

1. **FORTRAN IS GROUND TRUTH.** Not unit tests (fix the test if it disagrees). Not HANDOFF claims. Not previous sessions. Not LLM recall. The `njoy-reference/src/*.f90` tree is the oracle. Reference tapes in `njoy-reference/tests/NN/` are the acceptance bar. Bit-identical stretch goal; ±1 ULP at 7 sigfigs is the first-round floor; structural match (line counts, section presence) is non-negotiable at any tolerance.

2. **SKEPTICISM.** Verify every subagent report twice. Verify HANDOFF claims against current Fortran source and current Julia code (both drift). Verify previous "irreducible" FP precision floors — every one so far has been a real bug (see Phase 34, T18, T27, T49). Verify your own assumptions against a REPL check or a `write(*,...)` in Fortran + recompile.

3. **ALL BUGS ARE DEEP.** No bandaids. No hardcoded constants to silence a diff. No "good enough for this test" shims. Root-cause every failure. The coincidence-shading bug that cost T34 five MTs looked like FP noise; it was a one-line comparison against the wrong shaded value. The MT=18 fission doubling looked like a URR evaluator bug; it was a missing `mtr18` flag. Assume every bug is that deep.

4. **TIERED WORKFLOW.** Scale effort to change size:
   - **Trivial** (<5 LOC, typo/rename/comment): direct fix, run one reference test, no subagents.
   - **Small** (<30 LOC, one function, one MT fix): 1 research subagent (read the Fortran subroutine) + TDD + self-review.
   - **Core** (new module / new reader / >30 LOC / cross-module): **3+1 pattern** — 3 parallel read-only research subagents (Fortran reader + ENDF inspector + Julia-state comparator) + 1 Julia runner. Proposers must not see each other's output. Orchestrator synthesizes. This is the pattern that found the MF=12 histogram breakthrough, the Frobenius-Schur inversion, and the XQ indexing bug.

5. **GET FEEDBACK FAST.** Never code 500 lines before running a reference test. Every ~50 lines targeting a known test: run it. Always:
   ```bash
   rm -rf ~/.julia/compiled/v1.12/NJOY*   # non-negotiable before every run
   julia --project=. test/validation/reference_test.jl <N>
   ```
   Full 84-test sweep only in background or at session end (~90 min).

6. **FAIL FAST, FAIL LOUD.** `error("MF7/MT4 not found for MAT=$mat at offset $pos")`, not `return nothing`. Silent fallbacks corrupt downstream tapes. When the Fortran would abort (`call error`), Julia must abort too — and with enough context (MAT, MF, MT, file offset, numeric values) to diagnose without a rerun. Assertions over defensive defaults.

7. **LITERATE CODING.** Every non-trivial function docstring: WHAT it computes, WHICH Fortran subroutine it mirrors, WHICH ENDF-6 manual section defines the format. Cite line ranges. Comments explain the *why* — the subtle Fortran-specific behavior being matched (sigfig bias, accumulation order, scratch-tape vs MF3 read path, per-section AWR). The code shows *what*; the comment shows *why this particular Fortran quirk*.

8. **JULIA IDIOMATIC, FORTRAN SEMANTIC.** No transliteration. Parametric types over tagged unions. Multiple dispatch over `isa`-cascades. ScopedValues / function args over globals. But when Fortran rounds to 8 sigfigs before adding potential scattering, or accumulates in a specific IEEE 754 order, or uses a sentinel value, **match it exactly** — just express it idiomatically. A Vector comprehension is fine; a Fortran-style `do k=1,n` loop with 1-based mutation of a workspace is not. See existing `reich_moore.jl` for the pattern: Fortran's `frobns`/`thrinv`/`abcmat` ported as clean Julia helpers with matching FP semantics.

9. **NO PARALLEL JULIA AGENTS.** Precompilation cache corruption is real and costs hours. Read-only research subagents (Fortran reading, ENDF parsing, Grep/Glob exploration) **can and should** run in parallel. Anything invoking `julia` — tests, REPL checks, oracle generation, sweep runs — **serial only, one process at a time**. If a subagent's description mentions running Julia, it does not run in parallel with any other Julia-running subagent.

10. **RESEARCH FIRST.** Before touching an unfamiliar module: 15 minutes in `HANDOFF.md` (Architecture + Traps + relevant Phase history), then the Fortran subroutine map (HANDOFF §"Fortran subroutines you'll need to read"), then the specific `.f90` file. Before a new ENDF record type: read ENDF-6 manual for that MF/MT. Before a new covariance format (LB=0/1/2/5/6): read both the manual and the Fortran reader.

11. **LOC LIMIT (SOFT).** New source files target ~300 lines. This is a *guideline, not enforced* — splitting a cohesive Fortran port across arbitrary file boundaries hurts more than a 500-line file that tracks one subroutine. Split when a file covers multiple unrelated concerns or a future reader would struggle; don't split for its own sake.

12. **DEMAND ELEGANCE (BALANCED).** On non-trivial changes, pause: "is there a more Julia way?" Multiple dispatch, broadcasting, clean types — yes. Exotic metaprogramming to save 3 lines — no. Fortran-style state machines when a clean Julia alternative exists — no. Skip this step for one-line fixes.

13. **END WITH A PASSING REFERENCE TEST.** Every logical feature group ends with at least one T## moving BIT_IDENTICAL (or NUMERIC_PASS at 1e-7, or from CRASH → DIFFS). Record the delta in the session's `worklog/TNN_*.md` entry and in the next HANDOFF update. No silent "this should work now" — prove it with the sweep.

14. **TAPE-BASED ARCHITECTURE.** Modules talk via tapes, not shared mutable state. `run_njoy(input)` dispatches each deck line to `<module>_module(...)`. Each module reads its input tape(s), writes its output tape(s), period. No globals, no module-level caches of PENDF data, no side channels. Mirrors Fortran's tape-unit plumbing. If a new module needs data from upstream, it reads the tape — not a Julia struct passed in memory.

15. **SESSION CLOSE PROTOCOL** (mandatory):
    1. `git status` — see what changed
    2. `git add <specific files>` — stage (never `-A`; avoids .env, /tmp leftovers)
    3. `git commit -m "..."` — descriptive, references T## and phase
    4. `git push` — to remote
    5. `bd dolt push` — beads to Dolt
    6. Write `worklog/TNN_*.md` if a phase closed (mirrors `HANDOFF.md` entries)
    7. Update `HANDOFF.md` "Current State" + "Immediate Next Steps" if status changed
    8. `bd remember "<surprising insight>"` for cross-session lessons
    Work is NOT complete until `git push` succeeds.

16. **REREAD THIS FILE** at session start and after any context compression.

---

## THE GROUND-TRUTH PRINCIPLE

**The Fortran is the architecture. If a Julia module produces output the Fortran wouldn't, the Julia is wrong, even if it "looks more correct".**

We are not writing a modern nuclear data code. We are not fixing NJOY's quirks. Sigfig bias of `1.0000000000001`; `third = 0.333333333` (truncated, not `1/3`); dead code `scr(5)=1` that is overwritten two lines later; MT=18 redundancy flag `mtr18`; the scratch-tape data path in `emerge`; cross-compiler ±1 ULP at 7 sigfigs — all of it ships verbatim. Every "cleaner" formulation that changes the output is wrong. Every "obvious mathematical equivalent" whose FP accumulation order differs is wrong. The reference is the reference.

The only freedom is *how* the Julia expresses it: dispatch, broadcasting, types, module structure, names. The *what* is frozen.

---

## Build & test

```bash
# Before every run — non-negotiable
rm -rf ~/.julia/compiled/v1.12/NJOY*

# Single reference test (fast, preferred during dev)
julia --project=. test/validation/reference_test.jl 7

# Full 84-test sweep (~90 min on this machine, write reports/REFERENCE_SWEEP.md)
julia --project=. test/validation/sweep_reference_tests.jl

# Unit tests + reference sweep via Pkg.test
julia --project=. -e 'using Pkg; Pkg.test()'

# Single module pipeline (e.g. T02 end-to-end)
julia --project=. -e 'using NJOY; run_njoy("njoy-reference/tests/02/input"; work_dir="/tmp/t02")'
```

**Never** launch multiple `julia` processes concurrently. Subagents that run Julia are serial; readers are parallel (see Rule 9).

---

## Grind Method cheat sheet

1. **Pick a test** (ordered: LRU=0 → SLBW → MLBW → Reich-Moore → SAMMY → URR; then broadr/thermr/heatr/groupr/errorr chains).
2. **Generate oracle**: `julia --project=. test/validation/diagnose_harness.jl <N>`
3. **Run Julia reconr → PENDF** (or full chain via `run_njoy`).
4. **Diff MF3 cols 1-66** (ignore sequence numbers in 76-80). First differing byte → root cause.
5. **Classify the diff**: grid (energies differ), XS ±1 (sigfig boundary), large XS (missing feature), format (a11 variant), line count (extra/missing grid points from shading/thresholds). See HANDOFF §Grind Method for common causes.
6. **Read the relevant Fortran subroutine**. Patch with `write(*,...)` if needed, `cd njoy-reference && make` to recompile.
7. **Fix, rerun, repeat until PERFECT**. Then rerun T01/T02/T08/T27/T34/T45/T46 to check for regression.

3+1 agent pattern for the hard diffs: 3 parallel read-only agents (Fortran reader / ENDF inspector / Julia grid comparator) + 1 Julia runner.

---

## Issue tracking: Beads

```bash
bd ready                         # available work
bd show <id>                     # details
bd update <id> --claim           # claim atomically
bd close <id> --reason "..."     # complete
bd create --title="..." --description="..." -t bug|feature|task -p 0-4
bd dep add <issue> <depends-on>  # issue depends on depends-on
bd blocked                       # dependency view
bd remember "<insight>"          # cross-session memory (use for surprises only)
bd memories <keyword>            # search
bd dolt push                     # MANDATORY at session close
```

Rules:
- Every non-trivial change starts with a bead. File the bead BEFORE writing the failing test.
- Do NOT use TodoWrite, TaskCreate, or markdown TODOs.
- Respect dependencies. `bd blocked` before picking up work.
- Session end: `bd dolt push` is part of the close protocol.

---

## Julia Idiom Cheatsheet

### DO
```julia
# Parametric types for dispatch
struct MF3Section{T<:AbstractFloat}
    mt::Int; mf::Int; awr::T; tab::TAB1{T}
end

# Multiple dispatch on resonance formalism
_add_peak_nodes!(nodes, p::SLBWParameters, el, eh) = ...
_add_peak_nodes!(nodes, p::MLBWParameters, el, eh) = ...
_add_peak_nodes!(nodes, p::ReichMooreParameters, el, eh) = ...
_add_peak_nodes!(nodes, p::SAMMYParameters, el, eh) = ...  # was missing — Bug 7

# Concrete small Unions for type stability
const ResonanceParameters = Union{SLBWParameters, MLBWParameters,
                                  ReichMooreParameters, SAMMYParameters}

# Cite the Fortran source
# Ref: reconr.f90:2845-2847 (csslbw)
elastic = round_sigfig(sigp2, 8, 0) + round_sigfig(spot, 8, 0)
```

### DO NOT
```julia
# ✗ isa cascades instead of dispatch
if p isa SLBWParameters ... elseif p isa ReichMooreParameters ... end

# ✗ Fortran-style state in globals
CURRENT_MAT = 0; CURRENT_TEMP = 0.0     # ✗ pass as args or ScopedValue

# ✗ "Mathematically equivalent" rearrangements that change FP order
total = sum(sigfig_partials)             # ✗ if Fortran does tot += sn in loop order

# ✗ Silent fallbacks
haskey(mf3, 1) || return pendf_in        # ✗ at least warn; prefer error unless verified

# ✗ Hardcoded "test-specific" constants
thnmax = 4.81207e6                       # ✗ — compute from reaction thresholds
```

---

## NJOY.jl-specific lessons (from HANDOFF §Traps + worklog sessions)

1. **Precompile cache corrupts** under concurrent Julia. Always nuke before a run.
2. **`sigfig` bias** — `round_sigfig` multiplies by `1.0000000000001`. Use `_dedup_tol!`, not `unique!`, to dedupe.
3. **Constants are hardcoded in Fortran** — `ehigh=20e6`, `elow=1e-5`, `elim=min(0.99e6, eresr)`, `third=0.333333333` (truncated). Never read from MF2.
4. **Per-section AWR** — MF3 HEAD record's AWR may differ from MF2's. Use per-section.
5. **Multi-material tapes** — always filter `find_section(io, 2, 151; target_mat=MAT)`.
6. **CGS units** in `constants.jl`, matching `phys.f90`. Don't "fix" to SI.
7. **MF2 has multiple ranges** (resolved + URR). Reader must parse both.
8. **MT=455 LIST preamble** (LNU=2, delayed ν̄) — between HEAD and TAB1. MT=452/456 don't.
9. **errorr output grid by `ign`** — `-1` → union; `1` → user_egn; `≥2` → library structure. Union is NEVER the output grid for `ign ≥ 1`.
10. **Fortran emerge reads MF3 from lunion's scratch tape**, not the original ENDF. Julia's direct MF3 interpolation differs subtly for INT≠2 sections.
11. **Initial vs mid-data duplicate shading differs** — `sigfig(e,7,0)` for initial (label 207), `sigfig(e,7,-1)` for mid-data (label 270).
12. **MT=18 redundant when MT=19 exists** — `mtr18` flag; lunion skips MT=18, recout computes as sum.
13. **MT=103-107 redundant when partials (600-849) exist** — sums only; otherwise output normally.
14. **T01 full-chain thermal diffs** at 1e-5 are sigma1 Doppler broadening FP accumulation class — sub-ULP at errt boundary, same class as T34.
15. **err values: always read from `njoy-reference/tests/NN/input`** — HANDOFF tables have had wrong values (T19 was 0.02 not 0.001; T04 is 0.10 not 0.005).
16. **Frobenius-Schur for Reich-Moore** — use `_frobns`/`_thrinv!`/`_abcmat`, not `inv(SMatrix{3,3,ComplexF64})`. Different FP rounding.
17. **3+1 agent pattern found every subtle bug** — MF=12 histogram, coincidence shading, XQ indexing, ajku Faddeeva.
18. **"Irreducible FP precision" has always been a real bug.** Never accept ±1 as unfixable without gdb tracing the Fortran.

---

## Session Close Protocol — MANDATORY

Before saying "done":

```bash
[ ] git status                         # what changed?
[ ] git add <specific files>           # stage (never -A)
[ ] git commit -m "T## phase summary"  # descriptive, reference phase/test
[ ] git push                           # to remote
[ ] bd dolt push                       # beads to Dolt
[ ] worklog/TNN_*.md written?          # if phase closed
[ ] HANDOFF.md updated?                # Current State + Next Steps
[ ] bd remember "<surprise>"?          # if a non-obvious lesson landed
```

Work is NOT complete until `git push` succeeds. If push fails, resolve and retry. Never leave work stranded locally.
