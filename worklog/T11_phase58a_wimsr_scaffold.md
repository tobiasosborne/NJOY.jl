# T11 Phase 58a — wimsr orchestration scaffold

**Date:** 2026-05-02
**Test:** T11 (Pu-238, ENDF/B-IV t404, WIMS library generation)
**Goal (Phase 58 overall):** Wire `:wimsr` into NJOY.jl's pipeline so Test 11
moves from `MISSING_TAPE` (no dispatcher) to `DIFFS` (real attempt) to
`BIT_IDENTICAL` (full Fortran-faithful port). This phase (58a) lands the
scaffold; xsecs/resint data extraction is Phase 58b.

## Why this module

Per HANDOFF "Open Work", four NJOY modules are unwired: **wimsr**, mixr,
powr, resxsr. Of these, **wimsr is the only one with a built-in reference
test** — T11 (`njoy-reference/tests/11`) calls the moder→reconr→broadr→unresr→
thermr→groupr→**wimsr**→moder chain and ships `referenceTape27` /
`referenceTape28` as the bit-identical acceptance bar. None of the 84
reference tests calls mixr, powr, or resxsr. So wimsr is the natural next
target.

Existing prep already in place:

- `src/formats/wimsr.jl` (516 LOC) — `WIMSMaterial`, `WIMSResonanceTable`,
  `WIMSP1Block`, validators, `write_wims(io, ::WIMSMaterial)` writer.
- All upstream T11 modules already produce real output (Pu-238 reconr is
  bit-identical via T02/T10/T11 reconr-only shared baseline).

## Outcome

**Phase 58a — scaffold landed; standalone test promoted from
`wimsr_module not defined` (precompile error) to `2/3 assertions pass,
1/3 fails on bit-identical (stub data)`.**

Concrete deltas:

| Metric                                  | Before   | After (Phase 58a)        |
|-----------------------------------------|----------|--------------------------|
| `:wimsr` in `AVAILABLE_MODULES`         | NO       | YES                      |
| Pipeline dispatcher case for `:wimsr`   | NO       | YES (calls `wimsr_module`)|
| `parse_wimsr` (8-card parser)           | NO       | YES                      |
| `WimsrParams` struct                    | NO       | YES (28 fields)          |
| `wimsr_module(tapes, params)`           | NO       | YES (stub-data builder)  |
| Standalone test fixture                 | NO       | YES (`test_wimsr_t11_standalone.jl`) |
| Oracle pair (gendf input + WIMS output) | NO       | YES (under `oracle_cache/test11/wimsr_standalone/`) |
| Standalone test status                  | (none)   | **2/3 pass** (red bar on bit-identical, expected) |

Regression-clean: T01 unchanged at NUMERIC_PASS 32812/32962 @1e-5
(confirmed via `reference_test.jl 1`). T22 / leapr / reconr untouched —
no `:wimsr` in those decks, no dispatch path change.

## What landed

### 1. Card parser (`src/orchestration/input_parser.jl`)

Added `WimsrParams` struct (28 fields) + `parse_wimsr(::ModuleCall)` per
the wimsr.f90 driver card layout:

| Card | Fields | Conditional |
|------|--------|-------------|
| 1    | `ngendf nout`                                | always |
| 2    | `iprint iverw igroup`                        | always (defaults 0/4/0) |
| 2a   | `ngnd nfg nrg igref`                         | only if `igroup==9` |
| 3    | `mat nfid rdfid iburn`                       | always |
| 4    | `ntemp nsigz sgref ires sigp mti mtc ip1opt inorf isof ifprod jp1` | always (12 fields, all optional with defaults) |
| 5    | `ntis efiss`                                 | only if `iburn>0` |
| 6a-c | burnup chain `(ident, yield)` pairs          | only if `iburn>0` |
| 7    | `glam(1:nrg)`                                | always |
| 8    | `p1flx(1:jp1)`                               | only if `jp1>0` |

Card-3 `nfid`/`rdfid` cross-fill: Fortran uses `nfid=nint(rdfid)`. T11's
deck `1050 1 1050.` provides both — we pick whichever was supplied if the
other defaulted to 0.

Cited Fortran lines: wimsr.f90:48-307 (driver), 170-172 (card-2 defaults),
186-189 (card-2a gating), 202-205 (card-3), 217-218 (card-4), 242-255
(burnup data).

### 2. Pipeline dispatch (`src/orchestration/pipeline.jl`)

```julia
elseif mc.name == :wimsr
    params = parse_wimsr(mc)
    wimsr_module(tapes, params)
```

`AVAILABLE_MODULES` gained `:wimsr`. `NJOY_MODULES` already contained it —
so previously T11's `wimsr` keyword was recognized but produced no output
(silent miss in the dispatch loop).

### 3. Module skeleton (`src/orchestration/modules/wimsr.jl`)

`wimsr_module(tapes, params)`:
- Resolves `params.ngendf` to a real GENDF input path; bails to a stub
  empty tape if missing (matches the existing `gaspr_module` /
  `acer_module` graceful-degrade pattern).
- Calls `_wimsr_build_material_stub(in_path, params)` (Phase 58a stub —
  produces a minimally-valid `WIMSMaterial` with all numeric fields zeroed,
  resonance/fission flags off so `validate(::WIMSMaterial)` passes).
- Writes via `write_wims(io, wmat)` from `src/formats/wimsr.jl`.

The stub builder is **explicitly temporary**: when xsecs / resint port
lands (Phase 58b), it's deleted and the orchestration calls those
extractors directly.

### 4. Failing standalone test (`test/validation/test_wimsr_t11_standalone.jl`)

Drives `wimsr_module` against an oracle pair:
- **Input**: `gendf_input.tape` — Fortran's full-precision ASCII GENDF
  (Pu-238 at 300/900/2100K, 69-group structure, MT={1, 2, 16, 17, 18, 102,
  221} + transfer matrices). 10 379 lines.
- **Target**: `wims_reference.tape` — Fortran's `referenceTape27`
  (byte-identical). 1 169 lines.

The oracle pair is regenerable; the deck used is documented at the bottom
of the test file (inject `moder \n -26 29/` immediately before the wimsr
call to dump groupr's binary tape26 as ASCII tape29 with full precision
preserved). Approach validated:
`cmp /tmp/t11_oracle_v2/tape27 njoy-reference/tests/11/referenceTape27`
returns 0 (Fortran reproduces reference exactly).

Current test state: **2/3 pass**:
- ✅ `wimsr_module` is defined
- ✅ Output file is written
- ❌ Output is bit-identical (84 differing lines vs ref 1169 — stub data
  is mostly zeros; first-diff at line 1 because Fortran emits the burnup
  sentinel `999999999 3` while the stub has no burnup).

## Open: Phase 58b

The remaining work to flip the third assertion:

### 1. Port `xsecs` (njoy-reference/src/wimsr.f90:863-1422)

Loops `(temperature × sigma_zero × MT)`. Routes per:

- MF=3 MT=1 → flux/p1flx ; MT=2 (skip, used in RI only) ;
  MT=16/24/875-891 → sn2n ; MT=17/25 → sn2n×2 ;
  MT=18-21,38 → sfi/sf0/ab0 ; MT=102-150 → abs1/abs2 ;
  MT=mti(221) → nth boundary ; MT=252 → xi
- MF=5 MT=452 → snu, MT=455 → snus
- MF=6 MT=2 → elastic matrix ; MT=18-21,38 → fission spectrum matrix ;
  MT=51-58,91,221 → inelastic/thermal matrix

Critical: when `sgref ≥ 1e10` ⇒ `isg=0` ⇒ infinite-dilution mode (read
last sig0 column; T11 has `sgref=1e10`). Otherwise `isg=iz` selects the
sig0 list index matching the reference dilution.

### 2. Port `resint` (wimsr.f90:425-676)

Builds the 2D resonance integral tables (sigma_zero × temperature) for
groups `nfg+1..nfg+nrg`. T11 has `ires=3` (3 temps × 7 sig0 = 21 integrals
per resonance group), 13 resonance groups, `ifis=3` ⇒ both absorption RI
and fission RI tables.

Goldstein-Cohen normalization: `σ_abs,RI = σ_b·σ_a / (σ_b + σ_a)` where
`σ_b = sigp + glam[jg]·sigp` (line 657: `siglam=spot(nfg+jg)*glam(jg)`).

### 3. Audit `src/formats/wimsr.jl` writer

Per the research punchlist (research-agent-B), 6 verify-or-fix items:
1. Identifier block emission for WIMS-D vs WIMS-E
2. Burnup chain emission with `iburn=0` (WIMS-E sets a leading 0)
3. Temperature-independent XS block ordering: `4*nrg + 2*nnt = 106`
   values for T11 — verify nscr2 layout matches Julia's separate-block
   approach
4. Resonance table flat-data ordering `[ntemp temps, nsigz sig0, ntemps×nsigz integrals]`
5. Sentinel format `      999999999` — 15-char right-aligned (Fortran
   `(i15)`)
6. Engineering notation `1p,e15.8` — Julia's `@printf("%15.8E")` should
   match exactly (1.23456780E+03)

### 4. Iterate to bit-identical

Grind Method on T11 standalone: diff Julia tape27 vs Fortran tape27, fix
root cause, re-run. Then full-pipeline `reference_test.jl 11` (Julia
groupr → Julia wimsr) — that's where cross-module integration bugs
surface.

### 5. ip1opt path (T11 has ip1opt=1)

T11's deck card 4 omits the ip1opt field, defaulting to `1` (skip P1
matrices). So Phase 58b can land WITHOUT porting `p1scat` (wimsr.f90:1685-1935)
— filed as future work for tests that actually exercise `ip1opt=0`.

## Fortran source citations

- `njoy-reference/src/wimsr.f90:48-307` — `wimsr` driver
- `njoy-reference/src/wimsr.f90:170-172` — card-2 defaults
- `njoy-reference/src/wimsr.f90:186-189` — card-2a gating on `igroup==9`
- `njoy-reference/src/wimsr.f90:202-205` — card-3 `mat/nfid/rdfid/iburn`
- `njoy-reference/src/wimsr.f90:217-218` — card-4 12-field free read
- `njoy-reference/src/wimsr.f90:242-255` — card-5/6 burnup conditional
- `njoy-reference/src/wimsr.f90:863-1422` — `xsecs`
- `njoy-reference/src/wimsr.f90:425-676` — `resint`
- `njoy-reference/src/wimsr.f90:1989-2147` — `wimout` (writer)
- `njoy-reference/tests/11/input` — T11 deck (lines 68-73 are the wimsr
  call)
