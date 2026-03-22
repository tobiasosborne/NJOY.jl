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
| Test 02 | Pu-238 | SLBW + URR | 17/17 | 3567 pts exact | **BIT-IDENTICAL** |

### In Progress — Test 08 (Ni-61, Reich-Moore)

| Test | Material | Formalism | MTs | Status |
|------|----------|-----------|-----|--------|
| Test 08 | Ni-61 | Reich-Moore (LRF=3) | 3/18 | **3/18 PERFECT** |

MAT=2834, err=0.01, ENDF file `eni61`. Single resolved range [1e-5, 70000 eV].

**Root causes of 15 failing MTs (deeply investigated):**

1. **Grid difference: 41 extra Fortran points** (affects MT=1, 2, 102 — 1582 diffs each). Split: 17 extra below 1 eV, 24 extra above 1 eV. The Julia `adaptive_reconstruct` has the correct thermal tightening (err/5 below 0.4999 eV), step guard (estp=4.1), and sigfig midpoint rounding. The convergence test at the critical panel [0.0253, 0.028387] was verified: capture error 0.0030 < threshold 0.00485 → correctly converges. Yet the Fortran subdivides further. Likely cause: the Fortran's cumulative lunion grid (each MF3 section processed sequentially, merging into the previous grid) produces slightly different grid points than the Julia's all-at-once approach. The difference propagates through the step guard (est = 4.1 × (current_E - last_accepted_E)), causing cascading differences.

   **Key diagnostic**: Grid diverges at index 82 (E=0.0253). Julia: [0.0253, 0.028387], Fortran: [0.02529999, 0.02530001, 0.02530006, ..., 0.028387]. XS values match at common grid points.

2. **Pseudo-threshold skip** (affects MT=4, 51-58, 91, 103, 107 — 192-864 diffs each). The Fortran `emerge` tracks `ith` — the first grid point where sn > 0 — and only outputs from that point. For MT=4, this means output starts at E=69999.99 (resonance boundary, where the first nonzero inelastic background appears). The Julia `_get_legacy_section` for MT=4 starts at the computed physical threshold of the first inelastic level (E=68136.5), which includes leading zeros that Fortran skips.

3. **Threshold boundary shading**: The Fortran starts threshold sections with sigfig(thrxx, 7, -1) (XS=0) and sigfig(thrxx, 7, +1). The Julia uses the unshaded threshold.

**What's verified working:**
- Reich-Moore reader (`_read_rm_params`) parses correctly
- Reich-Moore evaluator (`cross_section_rm`) produces correct XS — bit-identical at 50+ compared grid points near 70000 eV boundary
- MF3 interpolation at discontinuities (duplicate breakpoints) works correctly
- Thermal tightening, step guard, sigfig rounding all implemented correctly
- Oracle cache generated at `test/validation/oracle_cache/test08/`

### Not Yet Attempted

| Test | Material | Formalism | Notes |
|------|----------|-----------|-------|
| Test 07 | U-235 | SLBW + URR(LRF=2) | Needs `_read_urr_lrf2` + `_csunr2` (mode=12) |
| Test 15-17 | U-238 JENDL | Likely RM | ENDF file not yet identified |

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

11. **Total = sum of sigfig'd sections, NOT component sum** — Fortran `lunion` skips MT=1 (line 1882), so emerge never processes MT=1. The total is accumulated from all non-redundant sections, each `sigfig(sn,7,0)`'d. Then `recout` applies a final `sigfig(total,7,0)`. In Julia: `round_sigfig` each primary channel before summing, then `round_sigfig` the total. MF3 backgrounds for primary channels (MT=2,18,102) must be added UNROUNDED to the resonance value — the combined value is then rounded. Non-primary backgrounds are rounded individually.

12. **ENDF float format has two modes** — The Fortran `a11` uses 9-sigfig fixed-point for values in (0.1, 1e7) with genuine precision, otherwise 7-sigfig scientific. Falls back to scientific when trailing zeros indicate only 7 sigfigs. Only applied to energy values (odd positions), not XS values (even positions in data pairs). See `pair_data` parameter in `_write_data_values`.

13. **URR table needs rdfil2 boundary nodes** — The Fortran `rdfil2` adds `sigfig(EL,7,-1)` (negative, overlap marker) and `sigfig(EL,7,+1)` to `eunr`, then sorts and deduplicates. This means the URR table has a node at `sigfig(EL,7,+1)` (e.g., 200.0001) NOT at `sigfig(EL,7,0)` (200.0). The `sigfig(EL,7,0)` entry is removed by dedup (too close to previous). Without this, sigunr interpolates at boundary energies instead of returning the exact table value, causing ±1 diffs. Also adds egridu intermediate points when step ratio > 1.26.

14. **lunion skips MT=251-300** — Line 1892: `if ((mth.ge.251.and.mth.le.300).and.mth.ne.261) go to 150`. These are non-cross-section quantities (mubar, xi, gamma). They do NOT contribute grid points.

15. **Duplicate MF3 breakpoints** — MF3 data often has duplicate x-values at resonance range boundaries (e.g., x=1.0 appears twice with y=14.35 and y=0.0). This encodes a discontinuity. The interpolation function handles this correctly.

---

## Immediate Next Steps

### Priority 1: Fix Test 08 grid near thermal

The adaptive reconstruction grid differs near 0.0253 eV. The Fortran resxs (reconr.f90:2240-2569) adds dense points near the thermal enode because the capture 1/v cross section changes within panels. Compare the Fortran resxs panel bisection logic against Julia's `adaptive_reconstruct` to find the exact divergence point. Key files: `src/processing/adaptive_grid.jl`, Fortran `resxs`/`panel`.

### Priority 2: Fix Test 08 threshold shading

The Fortran emerge outputs threshold reactions (MT=4, 51-91, etc.) starting with `sigfig(thrxx, 7, -1)` and `sigfig(thrxx, 7, +1)` as the first two energies. The Julia `_get_legacy_section` doesn't generate these shaded boundaries. Fix in `pendf_writer.jl`.

### Priority 3: Implement LRU=2/LRF=2 reader (mode=12)

Test 07 (U-235) needs `_read_urr_lrf2` (matching Fortran `rdf2u2`, reconr.f90:828-829) and `_csunr2` evaluator. This unlocks U-235 and Mn-55 tests.

### Priority 4: Grind BROADR

The oracle cache has `after_broadr.pendf` for Tests 01 and 02. This is the next module in the processing chain (Doppler broadening).

---

## Test Details

| Test | MAT | ENDF file | err | Formalism | Chain | Status |
|------|-----|-----------|-----|-----------|-------|--------|
| 84 | 128 | n-001_H_002-ENDF8.0.endf | 0.001 | LRU=0 | RECONR only | **BIT-IDENTICAL** |
| 01 | 1306 | t511 | 0.005 | LRU=0 | RECONR→BROADR→... | RECONR **BIT-IDENTICAL** |
| 02 | 1050 | t404 | 0.005 | SLBW+URR | RECONR→BROADR→UNRESR→... | RECONR **BIT-IDENTICAL** |
| 08 | 2834 | eni61 | 0.01 | Reich-Moore | RECONR→BROADR→HEATR→... | RECONR **3/18 PERFECT** |
| 07 | 1395 | t511 | 0.005 | SLBW+URR(LRF=2) | RECONR→BROADR→... | **Needs mode=12 reader** |

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

### Phase 5: Final 3 diffs → 17/17 BIT-IDENTICAL

Three fixes to achieve 100% bit agreement on Test 02:

1. **MT=1 total (753 diffs → 0)**: Fortran `lunion` skips MT=1 — the total is accumulated from all non-redundant sections, each `sigfig(sn,7,0)`'d, with a final `sigfig(total,7,0)` in recout. Fixed `merge_background_legacy` to `round_sigfig` each primary channel before summing, then `round_sigfig` the total. MF3 backgrounds for primary channels are now added UNROUNDED to the resonance value (matching Fortran where `gety1` bg is added to resonance before `sigfig`).

2. **MT=102 capture (1 diff → 0)**: At E=191.58 eV, the combined (bg + resonance) value was at a `sigfig(7)` rounding boundary. `format_endf_float` uses `@sprintf` (round-half-to-even) while Fortran's `sigfig` has a tiny upward bias (`10^(ndig-11)`). Fixed by applying `round_sigfig(elastic/fission/capture, 7)` explicitly in `merge_background_legacy` before storing, matching Fortran emerge line 4832.

3. **MT=2 elastic (1 diff → 0)**: At E=200.0001 eV (resolved/unresolved boundary), the URR table's first node was at `sigfig(200,7,0)=200.00000000002` but the Fortran table had `sigfig(200,7,+1)=200.0001` (from rdfil2 boundary nodes). The `sigfig(200,7,0)` entry was removed by Fortran's dedup. Fixed `build_unresolved_table` to include rdfil2 boundary nodes (`sigfig(EL,7,±1)` with overlap marking, `sigfig(EH,7,±1)`), then sort and deduplicate matching Fortran rdfil2 lines 856-869.
