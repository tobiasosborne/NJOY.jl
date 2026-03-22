# NJOY.jl Session Handoff

## What is this project

NJOY.jl is a Julia port of NJOY2016 — the standard nuclear data processing system used worldwide for reactor physics, criticality safety, and radiation transport. The original is 119,613 lines of Fortran 90. Our Julia version is ~14,000 lines.

The goal: produce **bit-identical** PENDF/ACE output matching all 85 of NJOY's own reference test problems. The port must be idiomatic Julia — composable, differentiable, no global state — not a transliteration.

**Repo:** https://github.com/tobiasosborne/NJOY.jl (GPL-3.0)

## MANDATORY RULES — READ FIRST

### The Standard: 100% Bit Agreement

The **only** acceptable outcome is byte-for-byte identical MF3 output with the Fortran NJOY2016 reference. "Close" (0.01% error, last-digit rounding) is still a bug. If unit tests fail after a change that makes output match Fortran, **fix the tests** — the Fortran is canonical ground truth.

### The Grind Phase Method

We work **module by module, test by test**:

1. Pick a test (start simple: LRU=0 → SLBW → MLBW → Reich-Moore → SAMMY)
2. Run Fortran oracle to get ground-truth PENDF for each module in the chain
3. Run Julia, write PENDF, compare MF3 data (first 66 chars of each line, ignoring sequence numbers)
4. Find first byte that differs. Trace to root cause. Read the Fortran. Fix.
5. Repeat until `diff` shows zero lines
6. Move to next test

### Operational Rules

1. **100% bit agreement with Fortran** — anything else = failure
2. **Fortran is canonical truth** — if unit tests break after matching Fortran, fix the tests
3. **No full test suite runs** — targeted tests only (full suite takes ages)
4. **No parallel Julia processes** — precomp cache corruption. One process at a time. Always `rm -rf ~/.julia/compiled/v1.12/NJOY*` before running.
5. **Module by module** — compare Fortran output vs Julia output for each processing step
6. **Be skeptical of everything** — verify claims from HANDOFF/agents, don't trust without checking
7. **Read the Fortran before writing Julia** — `njoy-reference/src/` is the authority
8. **Use bundled ENDF files** from `njoy-reference/tests/resources/` only
9. **Idiomatic Julia** — no Fortran transliterations. Multiple dispatch, broadcasting, clean types.
10. **3+1 agent workflow** for significant architectural changes (dual proposals, orchestrator, reviewer)

## Current State (after Grind Phases 3+4)

### What WORKS — bit-identical

| Test | Material | Formalism | MTs | Grid | XS Values |
|------|----------|-----------|-----|------|-----------|
| Test 84 | H-2 (deuterium) | LRU=0 | 4/4 PERFECT | 769 pts exact | All identical |
| Test 01 | C-nat (carbon) | LRU=0 | 29/29 PERFECT | 1033 pts exact | All identical |
| Test 02 | Pu-238 | SLBW+URR | 14/17 PERFECT | 3567 pts exact | 14 MTs identical |

Test 02 is the first resonance material (SLBW resolved + unresolved). The PERFECT MTs:
- All threshold reactions: MT=4, 16, 17, 51-59, 91
- Fission: MT=18 (including unresolved range 200-10000 eV)

### What does NOT yet work for Test 02

3 MTs remain with ±1 differences in the 7th significant digit:
- **MT=1 (total)**: 753 XS diffs — the total is recomputed from components in Julia (elastic+fission+capture+other_bg) but Fortran uses MF3/MT=1 background directly. Sigfig rounding of individual components produces a sum that differs by ±1 from the MF3/MT=1 total. Fix: use MF3/MT=1 + resonance_total for MT=1 output (but careful — non-primary MTs must still be included).
- **MT=2 (elastic)**: 1 XS diff at 200.0001 eV (44.83975 vs 44.83976) — URR table edge rounding
- **MT=102 (capture)**: 1 XS diff at 191.58 eV (3.227129 vs 3.227130) — SLBW sigfig rounding edge case

### Other items NOT yet done

- **PENDF headers**: MF1/MT451 format, MF2 EH value, MF12/MT=102 photon production not output
- **BROADR, HEATR, THERMR**: Not validated against per-module oracle
- **Unit tests**: May need updating — some test expectations were written against old (incorrect) behavior

## What was done in Grind Phase 3 (this session)

### 5 bugs fixed for SLBW resonance material support

**Pipeline restructure (`src/processing/reconr.jl`):**

1. **LRU=1 path now uses `lunion_grid`** — The old code used `build_grid` + `adaptive_reconstruct` + manual MF3 extension with `linearize_one_over_v!`. This produced too few grid points outside the resonance range. The new code matches the Fortran flow (reconr.f90:358): lunion → resxs → emerge. Specifically: `lunion_grid(mf3_sections, err; nodes=mf2_nodes)` builds the full union grid FIRST, then `adaptive_reconstruct` refines the resonance range, then results are merged. MF2 nodes (peaks, widths) seed the initial grid.

2. **Don't add exact eresl/eresh to resonance grid** — Fortran resxs starts from the FIRST lunion grid point ≥ eresl (the shaded version from `sigfig(EL,7,+1)`), not the exact boundary. Adding exact boundaries created duplicate near-boundary points.

**Threshold interpolation (`src/processing/pendf_writer.jl`, `reconr_evaluator.jl`):**

3. **Use correct interpolation law for threshold adjustment** — The old code hardcoded linear interpolation (`frac * tab.y[2]`) for near-threshold cross sections. But MF3 sections can use LinLog, LogLog, etc. The fix temporarily modifies the first breakpoint to `(thrxx, 0.0)` and calls `interpolate(tab, e)` which respects the actual law. This matches Fortran emerge line 1936 + `gety1`. New `_threshold_interp()` helper in pendf_writer.jl and `_threshold_interp_legacy()` in reconr_evaluator.jl.

**ENDF float formatting (`src/endf/io.jl`):**

4. **Extended 9-sigfig format matching Fortran `a11`** — The Fortran `a11` subroutine (endf.f90:882-981) uses fixed-point format for values in (0.1, 1e7) when they have genuine 9-sigfig precision. For example, adaptive grid energies like 18.2814725 format as `" 18.2814725"` (9 sigfigs in 11 chars) instead of `" 1.828147+1"` (7 sigfigs). Falls back to scientific when trailing zeros indicate only 7 significant digits. Only applied to energy values (odd positions in data pairs) via the `pair_data` parameter.

**Float accumulation order (`src/processing/reconr_evaluator.jl`):**

5. **Match Fortran total XS accumulation** — Changed `total = elastic + fission + capture + other_bg` to `total = other_bg + elastic + fission + capture` to match Fortran emerge's `sig(1) = sig(1) + sig(2) + sig(3) + sig(4)`. This prevents 8th-sigfig differences from float non-associativity that the extended format would expose.

## What was done in Grind Phase 2

### 16 RECONR bugs fixed

**Grid construction (`src/processing/reconr_grid.jl`):**

1. **New `lunion_grid()` function** — unified DFS matching Fortran `lunion` (reconr.f90:1771-2238). Single bisection loop that simultaneously forces decade points (1,2,5 × 10^n), checks energy ratios for linear sections, and checks interpolation error for nonlinear sections. This replaces the old separated approach (`_linearize_mf3!` + `linearize_one_over_v!`).

2. **Arithmetic midpoints** — was `sqrt(lo*hi)` (geometric), now `(lo+hi)/2` (arithmetic) matching Fortran line 2138. Rounded to 7 sigfigs via `round_sigfig(mid, 7)`.

3. **Single stpmax for linear sections** — was applying 5x tighter tolerance for ALL sections below 0.5 eV. Fortran only does this for NONLINEAR sections (line 2150).

4. **elim = 0.99e6** — was 1.0e6. The Fortran constant is `0.99e6_kr` (line 1797), then `elim = min(0.99e6, eresr)`.

5. **ehigh = 20e6 constant** — was using MF2 EH. The Fortran `ehigh` is a parameter `20.e6_kr` (line 144), NOT read from MF2.

6. **Panel splitting at existing grid nodes** — `_bisect_panel!` now splits MF3 panels at existing grid points before bisection, matching Fortran lunion labels 210-280 which merge the existing grid into panel processing.

7. **Tolerance-based deduplication** — `_dedup_tol!` removes near-duplicate energies (within 1e-9 relative) caused by the `sigfig` bias. Exact `unique!` missed these.

8. **Boundary filtering** — post-processing removes interior grid points at eresl/eresr/eresh, matching Fortran lunion lines 2196-2226. Not applied for LRU=0 materials.

9. **Threshold adjustment** — computes physical threshold `thrx = -QI * (AWR+1)/AWR`, shades up with `sigfig(thrx,7,+1)`, adds to grid. Only when first MF3 breakpoint is below threshold AND not in a zero-XS pseudo-threshold skip region.

10. **Pseudo-threshold skip** — doesn't add threshold when first TWO MF3 breakpoints both have zero XS, matching Fortran lines 1973-1976.

**XS evaluation (`src/processing/reconr_evaluator.jl`):**

11. **sigfig rounding** — each interpolated MF3 background value is rounded to 7 significant figures via `round_sigfig(bg, 7)`, matching Fortran emerge line 4836.

12. **Threshold suppression** — XS set to 0 at the threshold energy, matching Fortran emerge line 4795.

13. **Threshold-adjusted interpolation** — near the threshold, interpolate from the adjusted threshold energy `thrxx` instead of the original MF3 breakpoint, matching Fortran's modified first data point (line 1936).

**PENDF writer (`src/processing/pendf_writer.jl`):**

14. **Section ordering** — outputs MTs in ENDF file order (matching Fortran emerge), not hardcoded primary-first order.

15. **Pseudo-threshold skip in output** — skips leading zero-XS panels for threshold reactions, matching Fortran emerge.

16. **MT=4 redundant sum** — computed as sum of MT=51-91 with threshold-adjusted interpolation and sigfig rounding.

**Reader fix (`src/processing/reconr_types.jl`):**

17. **MF3Section QM** — was reading HEAD.C1 (=ZA) as QM. Fixed to read TAB1.C1 (the actual Q value).

## How to continue the grind

### Immediate next: Test 02 (U-234, SLBW resonances)

```bash
# 1. Check oracle cache exists
ls test/validation/oracle_cache/test02/

# 2. Check test input
cat njoy-reference/tests/02/input
# MAT=1050, err=0.001, ENDF: t404 (multi-material tape)

# 3. Run Julia RECONR and write PENDF
rm -rf ~/.julia/compiled/v1.12/NJOY*
julia --project=. -e '
using NJOY
result = reconr("njoy-reference/tests/resources/t404"; mat=1050, err=0.001)
write_pendf_file("/tmp/julia_test02.pendf", result; mat=1050, err=0.001)
println("Points: $(length(result.energies))")
'

# 4. Compare MF3 data (columns 1-66, ignoring sequence numbers)
# Use this pattern for each MT section:
julia --project=. -e '
using NJOY, Printf

function parse_all_mf3(filename)
    result = Dict{Int, Vector{String}}()
    lines = readlines(filename)
    local idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = NJOY._parse_int(p[71:72]); mt = NJOY._parse_int(p[73:75]); mat = NJOY._parse_int(p[67:70])
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

j = parse_all_mf3("/tmp/julia_test02.pendf")
f = parse_all_mf3("test/validation/oracle_cache/test02/after_reconr.pendf")
all_mts = sort(collect(union(keys(j), keys(f))))
for mt in all_mts
    jd = get(j, mt, String[]); fd = get(f, mt, String[])
    if isempty(jd) && !isempty(fd); @printf("MT=%3d: MISSING\n", mt)
    elseif jd == fd; @printf("MT=%3d: PERFECT ✓\n", mt)
    else; @printf("MT=%3d: %d vs %d lines, DIFF\n", mt, length(jd), length(fd))
    end
end
'
```

### Key insight: Test 02 uses the RESONANCE code path

Test 84 and Test 01 are both LRU=0 (no resonance parameters). Test 02 has actual SLBW resonances. The code path diverges at `reconr.jl` line ~250: instead of `lunion_grid`, it uses `build_grid` → `adaptive_reconstruct` → grid extension → `merge_background_legacy`. This path has NOT been grind-tested yet and likely has bugs.

The `reconr()` function (legacy interface) and `reconstruct()` function have duplicated but divergent code. Both need to produce identical output. Focus on whichever the test executor calls.

### After RECONR: validate BROADR

The oracle cache has `after_broadr.pendf` for Test 01. Compare Julia BROADR output against it. The test executor calls `execute_broadr` which calls `doppler_broaden`.

## Architecture of key files

| File | What it does | Fortran equivalent |
|------|-------------|-------------------|
| `src/processing/reconr.jl` | Top-level RECONR pipeline. Two interfaces: `reconr()` (legacy NamedTuple) and `reconstruct()` (PointwiseMaterial). | reconr.f90 main |
| `src/processing/reconr_grid.jl` | Grid construction: `lunion_grid()` (LRU=0 path), `build_grid()` (resonance path), `linearize_one_over_v!()`, `_bisect_panel!()` | reconr.f90 lunion |
| `src/processing/reconr_evaluator.jl` | XS evaluation: `merge_background_legacy()` adds MF3 backgrounds to resonance XS. Has threshold logic. | reconr.f90 emerge |
| `src/processing/reconr_types.jl` | Types: `MF3Section`, `CrossSections`, `PointwiseMaterial`. Reader: `read_mf3_sections()` | — |
| `src/processing/pendf_writer.jl` | PENDF output: `write_pendf()`, `_collect_reactions()`, `_get_legacy_section()`. Handles thresholds, redundant sums, pseudo-threshold skip. | reconr.f90 emerge/recout |
| `src/processing/adaptive_grid.jl` | Generic adaptive refinement: `adaptive_reconstruct()`, `round_sigfig()` | reconr.f90 panel loop |
| `src/processing/broadr.jl` | Doppler broadening: Sigma1 kernel | broadr.f90 |
| `test/validation/reference_comparator.jl` | Parses PENDF files, compares XS values | — |
| `test/validation/fortran_oracle.jl` | Runs Fortran NJOY binary with truncated input decks | — |

## Fortran reference: key subroutines to read

| Subroutine | File:Lines | What it does |
|-----------|-----------|-------------|
| `lunion` | reconr.f90:1771-2238 | Builds the union energy grid from MF3 sections. Bisects panels with forced decade points + ratio/interpolation tests. |
| `emerge` | reconr.f90:4646-4982 | Walks the union grid, evaluates each MF3 section, writes output. Handles thresholds, pseudo-thresholds, redundant sums. |
| `anlyzd` | reconr.f90:5442-5747 | Analyzes the MF3 dictionary to identify redundant reactions (MT=1,3,4,101, gas production). |
| `sigma` | reconr.f90:2571-2667 | Evaluates resonance cross sections at a single energy. |
| `sigfig` | util.f90:361-393 | Rounds to N significant figures with optional shading (+/-1 in last digit). Has `bias=1.0000000000001` multiplier. |
| `panel` | reconr.f90:2256-2569 | Adaptive reconstruction within resonance range. |

## Traps and lessons (from 3 sessions of debugging)

1. **Precomp cache corruption** — `rm -rf ~/.julia/compiled/v1.12/NJOY*` before EVERY test run. Weird errors after agent work = cache corruption.

2. **`sigfig` bias** — `round_sigfig` multiplies by 1.0000000000001 (matching Fortran). This creates near-duplicates that exact `unique!` misses. Use `_dedup_tol!` (tolerance-based).

3. **Fortran constants vs MF2 values** — `ehigh = 20e6` is a Fortran parameter, NOT read from MF2 EH. Similarly `elow = 1e-5`, `elim = 0.99e6`, `emax = 19e6` are all hardcoded.

4. **Threshold formula** — `thrx = -QI * (AWR+1)/AWR` where AWR is from MF2 (NOT AWR/amassn). `awin = 1` for neutron-induced reactions in both ENDF-6 and older formats.

5. **Pseudo-threshold skip** — Fortran skips MF3 panels where BOTH current AND next XS are zero (lines 1973-1976). This means threshold reactions with zero-XS padding at the start are trimmed.

6. **Panel splitting** — Fortran's lunion merges the existing grid INTO each section's panel processing. A panel [A,B] is split by any existing grid point G into [A,G] and [G,B]. This prevents over-bisection.

7. **Decade points** — forced within bisection (not injected upfront). Range: ipwr=-4 to 5 giving 1, 0.5, 0.2 × 10^ipwr. Thermal 0.0253 eV. Only below `elim`.

8. **ENDF format columns** — MAT=cols 67-70, MF=71-72, MT=73-75, NS=76-80. Data in cols 1-66. Compare cols 1-66 only (ignore sequence numbers).

9. **Julia 1.12 scoping** — `for` loops inside top-level scripts can't assign to globals. Wrap in functions or use `local`.

10. **CGS units** — constants.jl uses CGS (ergs, grams, cm/s) matching Fortran phys.f90. The previous session's CGS→SI "fix" was WRONG and was reverted.

11. **Multi-material tapes** — `find_section(io, 2, 151; target_mat=MAT)` must filter by MAT. Without this, multi-material tapes return wrong material's resonance data.

12. **MF3Section.QM was ZA** — the reader was storing HEAD.C1 (ZA) as QM instead of TAB1.C1 (actual Q value). Fixed.

## How to run things

```bash
# Clear precomp (ALWAYS do this before running Julia)
rm -rf ~/.julia/compiled/v1.12/NJOY*

# Quick targeted test for a specific material
julia --project=. -e '
using NJOY
result = reconr("njoy-reference/tests/resources/t511"; mat=1306, err=0.005)
println("C-nat: $(length(result.energies)) points")
'

# Write PENDF and compare against Fortran oracle
julia --project=. -e '
using NJOY
result = reconr("path/to/endf"; mat=MAT, err=ERR)
write_pendf_file("/tmp/julia_output.pendf", result; mat=MAT, err=ERR)
'
# Then compare MF3 data columns 1-66 as shown above

# Generate Fortran oracle for a test
julia --project=. test/validation/diagnose_harness.jl 2

# Check beads issues
bd ready
bd list
```

## Test details

| Test | MAT | ENDF file | err | Formalism | Chain | Status |
|------|-----|-----------|-----|-----------|-------|--------|
| 84 | 128 | n-001_H_002-ENDF8.0.endf | 0.001 | LRU=0 | RECONR only | **BIT-IDENTICAL** |
| 01 | 1306 | t511 | 0.005 | LRU=0 | RECONR→BROADR→HEATR→THERMR→GROUPR | RECONR **BIT-IDENTICAL**, rest untested |
| 02 | 1050 | t404 | 0.005 | SLBW+URR | RECONR→BROADR→UNRESR→GROUPR | RECONR **14/17 PERFECT** (MT=1: 753 ±1 diffs, MT=2,102: 1 diff each) |

## Open beads issues

Run `bd ready` to see current work. Key issues:
- `NJOY.jl-qkv` — Grind Test 02: Pu-238/U-234 RECONR bit-identical (P0)
- `NJOY.jl-gun` — Grind BROADR: validate against per-module oracle (P1)
- `NJOY.jl-dkm` — PENDF format: MF1 header, MF2, MF12/MT=102 (P2)
