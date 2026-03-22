# NJOY.jl Session Handoff

## What is this project

NJOY.jl is a Julia port of NJOY2016 — the standard nuclear data processing system used worldwide for reactor physics, criticality safety, and radiation transport. The original is 119,613 lines of Fortran 90. Our Julia version is ~15,000 lines.

The goal: produce **bit-identical** PENDF/ACE output matching all 85 of NJOY's own reference test problems. The port must be idiomatic Julia — composable, differentiable, no global state — not a transliteration.

**Repo:** https://github.com/tobiasosborne/NJOY.jl (GPL-3.0)

---

## MANDATORY RULES — READ THESE FIRST

### Rule 1: 100% Bit Agreement with Fortran

The **only** acceptable outcome is byte-for-byte identical MF3 output (columns 1-66 of each data line) with the Fortran NJOY2016 reference. "Close" (0.01% error, last-digit rounding) is still a bug. There is NO tolerance. One digit wrong in one number = failure.

### Rule 2: Fortran is Canonical Truth

If a unit test breaks after a change that makes output match Fortran, **fix the test** — the Fortran is ground truth. Never "fix" code to match a test if the Fortran disagrees.

### Rule 3: Read the Fortran Before Writing Julia

The Fortran source in `njoy-reference/src/` is the authoritative reference. Every formula, constant, rounding step, and accumulation order matters. Read the specific subroutine before implementing or fixing Julia code. Copy nothing blindly from the HANDOFF — verify against the source.

### Rule 4: No Parallel Julia Processes

Precompilation cache corruption is real and wastes hours. **Always** run:
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
```
before every test. One Julia process at a time. Weird runtime errors after agent work = cache corruption, not a code bug.

### Rule 5: Use Bundled ENDF Files Only

All test data lives in `njoy-reference/tests/resources/`. Never download or use external ENDF files.

### Rule 6: Be Skeptical of Everything

Verify claims from this HANDOFF, from agents, from comments in the code. The codebase has had multiple agents working on it. Previous "fixes" have been wrong. Read the Fortran, test the claim, then proceed.

### Rule 7: Idiomatic Julia

No Fortran transliterations. Use multiple dispatch, broadcasting, clean types. But when the Fortran does something specific (like rounding to 8 sigfigs before adding potential scattering), match it exactly — just express it idiomatically.

---

## The Grind Method

This is how we achieve bit-identical output. It works. Don't skip steps.

### Per-test workflow

1. **Pick a test** (ordered by complexity: LRU=0 → SLBW → MLBW → Reich-Moore → SAMMY)
2. **Generate the Fortran oracle**: `julia --project=. test/validation/diagnose_harness.jl <test_number>`. This runs the Fortran binary and caches per-module PENDF output in `test/validation/oracle_cache/test<NN>/`
3. **Run Julia RECONR** and write PENDF:
```bash
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e '
using NJOY
result = reconr("njoy-reference/tests/resources/<file>"; mat=<MAT>, err=<ERR>)
write_pendf_file("/tmp/julia_test.pendf", result; mat=<MAT>, err=<ERR>)
'
```
4. **Compare MF3 data** (columns 1-66, ignoring sequence numbers in columns 76-80):
```bash
julia --project=. -e '
using NJOY
function parse_all_mf3(fn)
    result = Dict{Int, Vector{String}}()
    lines = readlines(fn); idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = NJOY._parse_int(p[71:72]); mt = NJOY._parse_int(p[73:75])
        mat = NJOY._parse_int(p[67:70])
        if mf == 3 && mt > 0 && mat > 0
            idx += 3; data = String[]
            while idx <= length(lines)
                length(lines[idx]) < 75 && (idx += 1; continue)
                pd = rpad(lines[idx], 80)
                mfc = NJOY._parse_int(pd[71:72]); mtc = NJOY._parse_int(pd[73:75])
                (mfc != 3 || mtc != mt) && break
                push!(data, pd[1:66]); idx += 1
            end; result[mt] = data
        else; idx += 1; end
    end; return result
end
function compare(name, jf, ff)
    j = parse_all_mf3(jf); f = parse_all_mf3(ff)
    mts = sort(collect(union(keys(j), keys(f)))); np = 0
    for mt in mts
        jd = get(j, mt, String[]); fd = get(f, mt, String[])
        if jd == fd; np += 1
        elseif isempty(jd); println("$name MT=$(lpad(mt,3)): MISSING")
        else
            d = findfirst(i -> i > length(jd) || i > length(fd) || jd[i] != fd[i],
                          1:max(length(jd),length(fd)))
            println("$name MT=$(lpad(mt,3)): DIFF $(length(jd))v$(length(fd)) @$d")
            d !== nothing && d <= min(length(jd),length(fd)) &&
                println("  J: $(jd[d])\n  F: $(fd[d])")
        end
    end
    println("$name: $np/$(length(mts)) PERFECT")
end
compare("T02", "/tmp/julia_test.pendf",
        "test/validation/oracle_cache/test02/after_reconr.pendf")
'
```
5. **Find the first byte that differs**. Trace to root cause. Read the Fortran. Fix.
6. **Repeat** until all MTs show PERFECT.

### Debugging a diff

When you see a diff, classify it:
- **Grid diff** (different energies): issue is in `lunion_grid`, `adaptive_reconstruct`, or the initial grid construction
- **XS diff at same energy** (±1 in last digit): issue is in sigma evaluation, sigfig rounding, float accumulation order, or threshold interpolation
- **Large XS diff**: missing feature (e.g., unresolved resonance evaluation was completely missing)
- **Format diff** (same value, different string): issue is in `format_endf_float` (the `a11` equivalent)

For ±1 diffs, common root causes (all encountered and fixed in this project):
- Missing intermediate `round_sigfig(x, 8)` call before addition (SLBW)
- Wrong float accumulation order (non-associativity of IEEE 754)
- Testing total instead of partials in adaptive convergence
- Hardcoded linear interpolation instead of using MF3's actual law (LinLog, LogLog)
- Missing boundary-crossing forced convergence in adaptive reconstruction

---

## Current State

### What WORKS — bit-identical

| Test | Material | Formalism | MTs | Grid | Status |
|------|----------|-----------|-----|------|--------|
| Test 84 | H-2 (deuterium) | LRU=0 | 4/4 | 769 pts exact | **BIT-IDENTICAL** |
| Test 01 | C-nat (carbon) | LRU=0 | 29/29 | 1033 pts exact | **BIT-IDENTICAL** |
| Test 02 | Pu-238 | SLBW + URR | 14/17 | 3567 pts exact | **14/17 PERFECT** |

Test 02 PERFECT MTs: 4, 16, 17, 18, 51, 52, 53, 54, 55, 56, 57, 58, 59, 91

### Test 02 remaining diffs (3 MTs)

**MT=1 (total): 753 XS diffs (±1 in 7th sigfig)**

Root cause: Julia recomputes total as `other_bg + elastic + fission + capture` (sum of individually rounded components). Fortran emerge processes MT=1 using MF3/MT=1 background + resonance total (`res(1)`) directly. The sum of rounded components differs by ±1 from the pre-summed MF3/MT=1 value.

The fix is tricky: simply using MF3/MT=1 drops contributions from non-primary MTs (MT=51-91 etc.) that go through `other_bg`. The Fortran handles this because MT=1 in the ENDF file already includes all partial contributions. But switching to MF3/MT=1 requires ensuring the `other_bg` contributions are redundant with MT=1, which they are in principle but might not be numerically.

An attempted fix (using `total = res_total + total_bg` where total_bg = MF3/MT=1) broke Test 01 because it dropped `other_bg`. The correct approach may be: for MT=1 output in `_get_legacy_section`, use MF3/MT=1 background interpolated at each energy + resonance total, instead of the recomputed sum. This separates the output computation from the merge computation.

**MT=2 (elastic): 1 XS diff at 200.0001 eV**

44.83975 vs 44.83976. This is at the resolved/unresolved boundary. The URR table gives 44.83976 at E=200.0, but log-log interpolation at 200.0001 gives 44.83975. The Fortran might compute this differently (e.g., evaluating csunr1 directly at 200.0001 in resxs rather than interpolating from the table).

**MT=102 (capture): 1 XS diff at 191.58 eV**

3.227129 vs 3.227130. SLBW evaluation at a specific energy rounds to one side of a boundary. May require additional intermediate `sigfig` rounding (like the `sigfig(sigp(2),8,0)` fix for elastic, but for capture or total).

---

## Architecture

### Key source files

| File | What it does | Fortran equivalent |
|------|-------------|-------------------|
| `src/processing/reconr.jl` | Top-level RECONR pipeline. `reconr()` returns NamedTuple, `reconstruct()` returns PointwiseMaterial. | reconr.f90 main |
| `src/processing/reconr_grid.jl` | Grid construction: `lunion_grid()` (union grid builder matching Fortran lunion) | reconr.f90 lunion |
| `src/processing/reconr_evaluator.jl` | `merge_background_legacy()` adds MF3 backgrounds. `sigma_mf2()` evaluates resonance XS. | reconr.f90 emerge, sigma |
| `src/processing/pendf_writer.jl` | PENDF output. `_get_legacy_section()` handles thresholds, redundant sums. | reconr.f90 emerge/recout |
| `src/processing/adaptive_grid.jl` | Generic adaptive linearization: `adaptive_reconstruct()`, `round_sigfig()` | reconr.f90 resxs/panel |
| `src/resonances/slbw.jl` | Single-Level Breit-Wigner cross section evaluation | reconr.f90 csslbw |
| `src/resonances/unresolved.jl` | **NEW** — Unresolved resonance XS: `_csunr1`, `_gnrl`, table builder + interpolator | reconr.f90 csunr1, gnrl, genunr, sigunr |
| `src/resonances/reader.jl` | MF2 reader: SLBW, MLBW, RM, SAMMY, **URR (LRU=2)** | reconr.f90 rdfil2, rdf2u1 |
| `src/endf/io.jl` | ENDF I/O: `format_endf_float` (with 9-sigfig `a11` extension), `parse_endf_float` | endf.f90 a11, lineio |
| `src/constants.jl` | Physics constants in CGS matching `phys.f90` | phys.f90 |

### Fortran subroutines you'll need to read

| Subroutine | reconr.f90 lines | Purpose |
|-----------|-----------------|---------|
| `lunion` | 1771-2238 | Union energy grid from MF3 sections |
| `resxs`/`panel` | 2240-2569 | Adaptive reconstruction in resonance range |
| `sigma` | 2571-2667 | Dispatch to SLBW/MLBW/RM + unresolved |
| `csslbw` | 2669-2856 | SLBW cross section formula |
| `emerge` | 4646-4982 | Merge grids + evaluate + write output |
| `genunr` | 1628-1735 | Build unresolved XS table |
| `sigunr` | 1737-1769 | Interpolate from unresolved table |
| `csunr1` | 3826-4077 | Unresolved average XS (mode=11) |
| `gnrl` | 4498-4644 | Fluctuation integrals (Hwang quadrature) |
| `unfac` | 4473-4496 | Penetrability with AMUN |
| `rdf2u1` | 1312-1425 | Read LRU=2 unresolved MF2 data |
| `a11` | endf.f90:882-981 | Float → 11-char ENDF format |
| `sigfig` | util.f90:361-393 | Round to N sigfigs with shading |

### Pipeline flow (reconr.jl for LRU=1 + URR materials)

```
1. read_mf2 → MF2Data (includes resolved + unresolved ranges)
2. read_mf3_sections → Vector{MF3Section}
3. build_unresolved_table → URRTable (csunr1 + MF3 bkg at egridu nodes)
4. _add_mf2_nodes! → MF2 peak/width nodes + URR table energies
5. lunion_grid → bg_grid (union of MF3 panels + MF2 nodes + URR nodes)
6. filter to [eresl, eresh) → res_grid
7. adaptive_reconstruct(xs_partials, res_grid) → res_energies
   - xs_partials returns (elastic, fission, capture) — NOT total
   - force_boundaries=[eresr] prevents subdivision across resolved/unresolved boundary
8. merge bg_outside + res_energies → all_energies
9. for each energy: sigma_mf2 (resolved) + eval_unresolved (URR table) → res_xs
10. merge_background_legacy(all_energies, res_xs, mf3_sections) → final XS
11. write_pendf_file → PENDF output
```

---

## Traps and Lessons (from 4 debugging sessions)

### Critical traps

1. **Precomp cache corruption** — `rm -rf ~/.julia/compiled/v1.12/NJOY*` before EVERY run. Non-negotiable.

2. **`sigfig` bias** — `round_sigfig` multiplies by 1.0000000000001. This creates near-duplicates that `unique!` misses. Use `_dedup_tol!` instead.

3. **Constants are hardcoded in Fortran** — `ehigh = 20e6`, `elow = 1e-5`, `elim = 0.99e6`, `emax = 19e6`, `third = 0.333333333` (truncated, NOT 1/3). Never read these from MF2.

4. **Multi-material tapes** — `find_section(io, 2, 151; target_mat=MAT)` must filter by MAT. Without it, you get another material's data.

5. **CGS units** — `constants.jl` uses ergs, grams, cm/s matching `phys.f90`. A previous session changed to SI and broke everything. Don't touch.

6. **MF2 has MULTIPLE ranges** — Pu-238 has LRU=1 (resolved, 1-200 eV) AND LRU=2 (unresolved, 200-10000 eV). The reader must parse BOTH. `eresh = max(EH_resolved, EH_unresolved)`, `eresr = EH_resolved_only`.

### Subtle Fortran behaviors that must be matched

7. **SLBW sigfig rounding** — `sigfig(sigp(2),8,0)` and `sigfig(spot,8,0)` BEFORE adding potential scattering to elastic. Without this, values at 7th-sigfig boundaries round wrong (±1 in last digit). See reconr.f90:2845-2847.

8. **Adaptive convergence tests partials only** — The Fortran tests j=1..nsig-1 = (elastic, fission, capture), NOT total. Testing total makes convergence stricter, producing extra grid points. See reconr.f90:2394.

9. **Boundary-crossing forced convergence** — Panels straddling eresr (resolved/unresolved boundary) are forced to converge without subdivision. See reconr.f90:2353-2355. Implemented via `force_boundaries` field in AdaptiveConfig.

10. **Threshold interpolation uses MF3's own law** — Near thresholds, the first breakpoint is modified to (thrxx, 0.0) and the section's own interpolation law (LinLog, LogLog, etc.) is used — NOT hardcoded linear. Implemented via `_threshold_interp()`.

11. **Float accumulation order** — `total = other_bg + elastic + fission + capture` (matching Fortran: `sig(1) = sig(1) + sig(2) + sig(3) + sig(4)` where sig(1) already holds other_bg). Float non-associativity causes 8th-sigfig differences visible in the extended format.

12. **ENDF float format has two modes** — The Fortran `a11` uses 9-sigfig fixed-point for values in (0.1, 1e7) with genuine precision, otherwise 7-sigfig scientific. Falls back to scientific when trailing zeros indicate only 7 sigfigs. Only applied to energy values (odd positions), not XS values (even positions in data pairs). See `pair_data` parameter in `_write_data_values`.

13. **URR table needs egridu intermediate points** — The Fortran `rdf2u1` adds points from a standard 78-point energy grid between fission-width energy nodes when the step ratio > 1.26. It also adds `sigfig(EH, 7, -1)` as the last node. Without these, log-log interpolation of the URR table gives wrong values.

14. **lunion skips MT=251-300** — Line 1892: `if ((mth.ge.251.and.mth.le.300).and.mth.ne.261) go to 150`. These are non-cross-section quantities (mubar, xi, gamma). They do NOT contribute grid points.

15. **Duplicate MF3 breakpoints** — MF3 data often has duplicate x-values at resonance range boundaries (e.g., x=1.0 appears twice with y=14.35 and y=0.0). This encodes a discontinuity. The interpolation function handles this correctly.

---

## Immediate Next Steps

### Priority 1: Fix Test 02 MT=1 (753 ±1 diffs)

This is the largest remaining diff. The root cause is understood: Julia recomputes total from components, Fortran uses MF3/MT=1 directly.

**Approach**: In `_get_legacy_section` (pendf_writer.jl), for MT=1, instead of returning `result.total`, compute: interpolate MF3/MT=1 at each energy + add resonance total from `res_xs`. This separates the PENDF output from the internal merge computation, matching how Fortran emerge processes MT=1 (itype=0). The `_write_legacy_mf3` function calls `_get_legacy_section` which returns (energies, xs_values) — you can override the MT=1 case.

### Priority 2: Fix Test 02 MT=2 and MT=102 (1 diff each)

These might require more intermediate sigfig rounding in the SLBW or URR evaluation. Compare the Julia and Fortran values at the exact energy of the diff to find the divergence point.

### Priority 3: Grind Test 03 or another SLBW test

After Test 02, try another material with a different formalism (MLBW, Reich-Moore). Look at tests in `njoy-reference/tests/` and check their input decks.

### Priority 4: Grind BROADR

The oracle cache has `after_broadr.pendf` for Tests 01 and 02. This is the next module in the processing chain.

---

## Test Details

| Test | MAT | ENDF file | err | Formalism | Chain | Status |
|------|-----|-----------|-----|-----------|-------|--------|
| 84 | 128 | n-001_H_002-ENDF8.0.endf | 0.001 | LRU=0 | RECONR only | **BIT-IDENTICAL** |
| 01 | 1306 | t511 | 0.005 | LRU=0 | RECONR→BROADR→... | RECONR **BIT-IDENTICAL** |
| 02 | 1050 | t404 | 0.005 | SLBW+URR | RECONR→BROADR→UNRESR→... | RECONR **14/17 PERFECT** |

---

## Session History (what was done across 4 grind phases)

### Phase 1-2: LRU=0 materials (Tests 84, 01) → 33/33 PERFECT

Built `lunion_grid`, fixed 17 bugs in grid construction, XS evaluation, threshold handling, PENDF writing.

### Phase 3: SLBW resolved (Test 02) → 13/17 PERFECT

- Restructured LRU=1 pipeline to use lunion_grid
- Fixed threshold interpolation to use correct MF3 law
- Added 9-sigfig ENDF float format
- Fixed float accumulation order

### Phase 4: Unresolved resonances + remaining fixes → 14/17 PERFECT

- Implemented complete LRU=2 unresolved evaluation: `_csunr1`, `_gnrl`, `build_unresolved_table`, `eval_unresolved`
- Parsed LRU=2 MF2 data (`_read_urr_lfw1`)
- Added egridu dense grid + EH-shaded boundary node to URR table
- SLBW `sigfig(x,8,0)` rounding before potential scattering
- Partials-only adaptive convergence
- Boundary-crossing forced convergence (`force_boundaries`)
- Expanded adaptive range to [eresl, eresh) for full resonance + unresolved coverage
