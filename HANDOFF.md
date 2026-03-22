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
before every test. One Julia process at a time. Weird runtime errors after agent work = cache corruption, not a code bug. **Use subagents for research (Fortran reading, ENDF parsing, code exploration) but NEVER for running Julia tests in parallel.**

### Rule 5: Use Bundled ENDF Files Only

All test data lives in `njoy-reference/tests/resources/`. Never download or use external ENDF files.

### Rule 6: Be Skeptical of Everything

Verify claims from this HANDOFF, from agents, from comments in the code. The codebase has had multiple agents working on it. Previous "fixes" have been wrong. Previous HANDOFF analysis has been wrong (e.g., "coincident breakpoint shading" was misidentified — the real cause was MF=12 histogram shading). Read the Fortran, test the claim, then proceed.

### Rule 7: Idiomatic Julia

No Fortran transliterations. Use multiple dispatch, broadcasting, clean types. But when the Fortran does something specific (like rounding to 8 sigfigs before adding potential scattering), match it exactly — just express it idiomatically.

### Rule 8: Use Subagents Wisely

Launch multiple **research** subagents in parallel for maximum efficiency:
- **Explore agents** to read Fortran source, search for patterns, understand data structures
- **Research agents** to analyze ENDF files, extract breakpoints, check data formats
- **Comparison agents** to run a SINGLE Julia diagnostic (grid extraction, MF3 comparison)

**NEVER** launch multiple agents that run Julia processes — Rule 4 applies. Only ONE agent should execute Julia at a time. Others should do read-only research.

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
compare("T07", "/tmp/julia_test.pendf",
        "test/validation/oracle_cache/test07/after_reconr.pendf")
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
- **Line count diff** (same data, shifted lines): extra or missing grid points from shading or threshold handling

For ±1 diffs, common root causes (all encountered and fixed in this project):
- Missing intermediate `round_sigfig(x, 8)` call before addition (SLBW elastic)
- Wrong float accumulation order (non-associativity of IEEE 754)
- Testing total instead of partials in adaptive convergence
- Hardcoded linear interpolation instead of using MF3's actual law (LinLog, LogLog)
- Missing boundary-crossing forced convergence in adaptive reconstruction
- MF3 backgrounds added rounded instead of unrounded to resonance values
- Missing `round_sigfig(x, 7)` on primary channel output values
- URR table missing rdfil2 boundary nodes (sigfig(EL, 7, ±1) shading)
- MF=12 histogram sections not processed (creates shaded pairs at shared breakpoints)

### Effective subagent strategy (proven in Phase 6)

For complex debugging, launch 3+ research subagents **in parallel**:
1. **Fortran source reader** — have an Explore agent deeply read the relevant Fortran subroutine
2. **ENDF data inspector** — have an agent examine the raw ENDF file (grep for specific energies, list breakpoints, check interpolation laws, find MF=12/13/23 sections)
3. **Grid comparator** — have ONE agent run Julia and compare grids/MF3 output

This 3+1 pattern (3 read-only researchers + 1 Julia runner) was how the MF=12 breakthrough was found. The data inspector discovered that MF=12 MT=102 had histogram interpolation at the exact mystery energies — something that manual Fortran code tracing missed for hours.

---

## Current State

### What WORKS — bit-identical

| Test | Material | Formalism | MTs | Grid | Status |
|------|----------|-----------|-----|------|--------|
| Test 84 | H-2 (deuterium) | LRU=0 | 4/4 | 769 pts exact | **BIT-IDENTICAL** |
| Test 01 | C-nat (carbon) | LRU=0 | 29/29 | 1033 pts exact | **BIT-IDENTICAL** |
| Test 02 | Pu-238 | SLBW + URR(mode=11) | 17/17 | 3567 pts exact | **BIT-IDENTICAL** |
| Test 08 | Ni-61 | Reich-Moore (LRF=3) | 18/18 | 4825 pts exact | **BIT-IDENTICAL** |

### In Progress — Test 07 (U-235, SLBW + URR mode=12)

**Test details**: MAT=1395, err=0.005, ENDF file `t511`, resolved range [1, 82] eV (SLBW, LRF=1), unresolved range [82, 25000] eV (LRF=2, mode=12). 27 MTs total.

**Status: 0/27 PERFECT** — but mode=12 reader (`_read_urr_lrf2`) and evaluator (`_csunr2`) are implemented. The code compiles, runs, and produces output. Oracle cache exists at `test/validation/oracle_cache/test07/`.

**Grid analysis**: Julia has 6924 points, Fortran has 6944 (20 fewer). 114 total grid diffs. The diffs start at index 693 (~6.3 eV) — **inside the SLBW resolved range (1-82 eV)**, NOT in the URR region. The energies are very close but slightly different (e.g., Julia 6.2953 vs Fortran 6.295313), suggesting cascading differences from the initial grid.

**Root cause hypothesis**: The initial MF2 peak nodes for U-235's dense resonance structure (many narrow resonances in 1-82 eV) differ slightly between Julia and Fortran. This cascades through adaptive reconstruction. Test 02 (Pu-238) doesn't show this because its resonances are fewer and wider.

**What's been done for Test 07:**
1. Implemented `_read_urr_lrf2` in `reader.jl` — reads LRU=2/LRF=2 LIST records with 6-word stride per energy
2. Implemented `_csunr2` in `unresolved.jl` — interpolates ALL parameters (D, GX, GNO, GG, GF) at each energy, computes fluctuation integrals
3. Added `URR2Sequence` and `URR2Data` types for mode=12 data
4. Added dispatch in `read_mf2` for LRU=2 && LRF==2
5. Added `build_unresolved_table(::URR2Data, ...)` with mode=12 grid generation (1/1000 gap threshold)
6. Updated `eval_unresolved` to support linear interpolation (mode=12 INT=2) alongside log-log (mode=11)
7. Fixed `elim = min(0.99e6, eresr)` — was using default 0.99e6 for all materials
8. Checked: U-235 MF=12 sections (MT=4, 18, 102) all use linear interpolation (INT=2), NOT histogram — so MF=12 histogram shading is NOT a factor for Test 07

**What needs debugging:**
1. **Initial grid nodes** — Compare `_add_bw_peaks!` output for U-235 vs what Fortran's rdf2bw produces. U-235 has 130+ resonances in 1-82 eV; the `ndig` calculation for peak triplets might differ.
2. **URR XS accuracy** — After fixing the grid, verify the mode=12 cross sections match. The `_csunr2` evaluator may need tuning (the interpolation bracketing loop, small-value clamping, etc.).
3. **elim effects** — With elim=82 for U-235, lunion only forces decade points below 82 eV. Verify this is handled correctly.

### Not Yet Attempted

| Test | Material | Formalism | Notes |
|------|----------|-----------|-------|
| Test 15-17 | U-238 JENDL | Likely Reich-Moore | ENDF file not yet located in resources/ |

### Formalisms status

| Formalism | Reader | Evaluator | Tested |
|-----------|--------|-----------|--------|
| LRU=0 (no resonances) | N/A | N/A | Tests 84, 01 BIT-IDENTICAL |
| SLBW (LRF=1) | `_read_bw_params` | `cross_section_slbw` | Test 02 BIT-IDENTICAL |
| MLBW (LRF=2) | `_read_bw_params` | `cross_section_mlbw` | Not tested |
| Reich-Moore (LRF=3) | `_read_rm_params` | `cross_section_rm` | Test 08 BIT-IDENTICAL |
| SAMMY/RML (LRF=7) | `_read_sammy_params` | `build_rml_evaluator` | Not tested |
| URR mode=11 (LRU=2, LRF≤1, LFW=1) | `_read_urr_lfw1` | `_csunr1`/`_gnrl` | Test 02 BIT-IDENTICAL |
| URR mode=12 (LRU=2, LRF=2) | `_read_urr_lrf2` | `_csunr2`/`_gnrl` | Test 07 (runs, 0/27, grid diffs) |

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
| `src/resonances/mlbw.jl` | Multi-Level Breit-Wigner + Reich-Moore dispatch | reconr.f90 csmlbw, csrmat |
| `src/resonances/reich_moore.jl` | Reich-Moore R-matrix cross section evaluation | reconr.f90 csrmat |
| `src/resonances/unresolved.jl` | Unresolved resonance XS: `_csunr1`, `_csunr2`, `_gnrl`, table builder + interpolator | reconr.f90 csunr1, csunr2, gnrl, genunr, sigunr |
| `src/resonances/reader.jl` | MF2 reader: SLBW, MLBW, RM, SAMMY, URR (LRU=2 modes 11+12) | reconr.f90 rdfil2, rdf2u1, rdf2u2 |
| `src/processing/reconr_types.jl` | MF3 + MF12 readers: `read_mf3_sections`, `read_mf12_lo1_sections` | reconr.f90 lunion (MF processing) |
| `src/endf/io.jl` | ENDF I/O: `format_endf_float` (with 9-sigfig `a11` extension), `parse_endf_float` | endf.f90 a11, lineio |
| `src/constants.jl` | Physics constants in CGS matching `phys.f90` | phys.f90 |

### Fortran subroutines you'll need to read

| Subroutine | reconr.f90 lines | Purpose |
|-----------|-----------------|---------|
| `lunion` | 1771-2238 | Union energy grid from MF3+MF12 sections. Key: decade points (2112-2127), histogram shading (2050-2064), singularity shading (1930-1937, 2024-2033, 2093-2096), step ratio (2136-2141) |
| `resxs` | 2240-2569 | Adaptive reconstruction in resonance range. Key: thermal tightening (2390), step guard estp=4.1 (2419), sigfig midpoint ndig=8/9 rounding (2360-2372) |
| `sigma` | 2571-2667 | Dispatch to SLBW/MLBW/RM + unresolved (sigunr at 2659-2665) |
| `csslbw` | 2669-2856 | SLBW cross section formula. Key: sigfig(sigp(2),8,0) at line 2845 |
| `csrmat` | 2858-3023 | Reich-Moore R-matrix cross section formula |
| `emerge` | 4646-4982 | Merge grids + evaluate + write output. Key: itype dispatch (4756-4770), MF3 bg suppression in URR range (4800), sigfig(sn,7,0) at line 4832, ith pseudo-threshold tracking (4834) |
| `recout` | 4984-5441 | Write PENDF. Key: redundant reactions at line 5308 |
| `genunr` | 1628-1735 | Build unresolved XS table. Dispatches: mode=11→csunr1, mode=12→csunr2 (1682-1683) |
| `sigunr` | 1737-1769 | Interpolate from unresolved table |
| `rdfil2` | 700-881 | Read all MF2 ranges. Key: eunr boundary nodes (759-775), eunr sort+dedup (856-869) |
| `rdf2u1` | 1312-1425 | Read LRU=2 unresolved data (mode=11, energy-dependent fission widths) |
| `rdf2u2` | 1428-1540 | Read LRU=2 unresolved data (mode=12, all widths energy-dependent). 6-word stride per energy, gap filling with egridu at 1/1000 threshold |
| `csunr2` | 4079-4317 | Unresolved XS for mode=12. ALL parameters interpolated via terp1 at each energy. Clamps small values <1e-8 to 0. Final XS interpolation uses per-J INT (typically 2=linear) |
| `rdf2bw` | 884-1310 | Read SLBW/MLBW resonance parameters. Adds peak nodes to enode with sigfig rounding |
| `a11` | endf.f90:882-981 | Float → 11-char ENDF format |
| `sigfig` | util.f90:361-393 | Round to N sigfigs with shading. Key: bias=1.0000000000001 at line 390 |

### Pipeline flow (reconr.jl for LRU=1 + URR materials)

```
1. read_mf2 → MF2Data (includes resolved + unresolved ranges)
2. read_mf3_sections → Vector{MF3Section}
3. read_mf12_lo1_sections → Vector{MF3Section} (photon multiplicities)
4. build_unresolved_table → URRTable (csunr1/csunr2 + MF3 bkg at egridu nodes)
   - Includes rdfil2 boundary nodes: sigfig(EL,7,±1), sigfig(EH,7,±1)
   - Sort, dedup, overlap marking (negative energies for E < eresr)
5. _add_mf2_nodes! → MF2 peak/width nodes + URR table energies
6. lunion_grid → bg_grid (union of MF3+MF12 panels + MF2 nodes + URR nodes)
   - Shades duplicate MF3 breakpoints to sigfig(x,7,±1)
   - Shades histogram interior breakpoints to sigfig(x,7,±1) pairs
   - elim = min(0.99e6, eresr) — decade/ratio checks only below elim
   - Processes sections cumulatively (each sees previous grid)
7. filter to [eresl, eresh) → res_grid
8. adaptive_reconstruct(xs_partials, res_grid) → res_energies
   - xs_partials returns (elastic, fission, capture) — NOT total
   - force_boundaries=[eresr] prevents subdivision across resolved/unresolved boundary
   - Thermal tightening: err/5 below 0.4999 eV
   - Step guard: estp=4.1
   - Midpoint rounding: ndig=9 (or 8 for 0.1<xm<1), sigfig comparison check
9. merge bg_outside + res_energies → all_energies
10. for each energy: sigma_mf2 (resolved) + eval_unresolved (URR table) → res_xs
11. merge_background_legacy(all_energies, res_xs, mf3_sections) → final XS
    - Primary channels (MT=2,18,102): add UNROUNDED MF3 bg to resonance
    - Non-primary channels: round MF3 bg to 7 sigfigs before accumulating
    - Round each primary channel to 7 sigfigs (matching emerge line 4832)
    - Total = round_sigfig(other_bg + sigfig(el) + sigfig(fi) + sigfig(ca), 7)
12. write_pendf_file → PENDF output
```

---

## Traps and Lessons (from 6 debugging sessions)

### Critical traps

1. **Precomp cache corruption** — `rm -rf ~/.julia/compiled/v1.12/NJOY*` before EVERY run. Non-negotiable.

2. **`sigfig` bias** — `round_sigfig` multiplies by 1.0000000000001. This creates near-duplicates that `unique!` misses. Use `_dedup_tol!` instead.

3. **Constants are hardcoded in Fortran** — `ehigh = 20e6`, `elow = 1e-5`, `elim = 0.99e6`, `emax = 19e6`, `third = 0.333333333` (truncated, NOT 1/3). Never read these from MF2.

4. **Multi-material tapes** — `find_section(io, 2, 151; target_mat=MAT)` must filter by MAT. Without it, you get another material's data.

5. **CGS units** — `constants.jl` uses ergs, grams, cm/s matching `phys.f90`. A previous session changed to SI and broke everything. Don't touch.

6. **MF2 has MULTIPLE ranges** — Pu-238 has LRU=1 (resolved, 1-200 eV) AND LRU=2 (unresolved, 200-10000 eV). U-235 has SLBW (1-82 eV) + URR/LRF=2 (82-25000 eV). The reader must parse BOTH. `eresh = max(EH_resolved, EH_unresolved)`, `eresr = EH_resolved_only`.

### Subtle Fortran behaviors that must be matched

7. **SLBW sigfig rounding** — `sigfig(sigp(2),8,0)` and `sigfig(spot,8,0)` BEFORE adding potential scattering to elastic. Without this, values at 7th-sigfig boundaries round wrong (±1 in last digit). See reconr.f90:2845-2847.

8. **Adaptive convergence tests partials only** — The Fortran tests j=1..nsig-1 = (elastic, fission, capture), NOT total. Testing total makes convergence stricter, producing extra grid points. See reconr.f90:2394.

9. **Boundary-crossing forced convergence** — Panels straddling eresr (resolved/unresolved boundary) are forced to converge without subdivision. See reconr.f90:2353-2355. Also applies at eresu and eresm boundaries. Implemented via `force_boundaries` field in AdaptiveConfig.

10. **Threshold interpolation uses MF3's own law** — Near thresholds, the first breakpoint is modified to (thrxx, 0.0) and the section's own interpolation law (LinLog, LogLog, etc.) is used — NOT hardcoded linear. Implemented via `_threshold_interp()`.

11. **Total = sum of sigfig'd sections, NOT component sum** — Fortran `lunion` skips MT=1 (line 1882), so emerge never processes MT=1. The total is accumulated from all non-redundant sections, each `sigfig(sn,7,0)`'d. Then `recout` applies a final `sigfig(total,7,0)`.

12. **ENDF float format has two modes** — The Fortran `a11` uses 9-sigfig fixed-point for values in (0.1, 1e7) with genuine precision, otherwise 7-sigfig scientific. Falls back to scientific when trailing zeros indicate only 7 sigfigs. Only applied to energy values (odd positions), not XS values.

13. **URR table needs rdfil2 boundary nodes** — sigfig(EL,7,±1) and sigfig(EH,7,±1) must be in the URR energy grid. Without this, sigunr interpolates at boundary energies instead of returning exact table values.

14. **lunion skips MT=251-300** (except 261) and MT=1, MT=3, MT=101. Also conditionally skips MT=4 (if mtr4>0), MT=103-107 (if redundant), MT=18 (if mtr18>0).

15. **Duplicate MF3 breakpoints** — Encode discontinuities. Shaded to sigfig(x,7,-1) and sigfig(x,7,+1). Implemented in lunion_grid.

16. **elim = min(0.99e6, eresr)** — NOT a fixed constant. For U-235 (eresr=82), elim=82. For Pu-238 (eresr=200), elim=200. For Ni-61 (eresr=70000), elim=70000. Decade forcing and ratio checks only below elim.

17. **MF=12 histogram shading in lunion** — The Fortran lunion processes MF=12 LO=1 sections alongside MF=3. Histogram interpolation (INT=1) triggers the two-pass shading mechanism (reconr.f90:2050-2064): each interior breakpoint E → {sigfig(E,7,-1), sigfig(E,7,+1)}. For Ni-61, MF=12 MT=102 has INT=1; for U-235, all MF=12 sections use INT=2 (linear, no shading). **Check each material's MF=12 interpolation law!**

18. **Redundant MT=4 pseudo-threshold** — Fortran recout tracks first nonzero index (`ith`) and starts MT=4 from the energy before. Without this, output includes leading zeros. Implemented as pseudo-threshold skip in `_get_legacy_section`.

19. **Mode=12 csunr2 clamps small widths** — GX and GF values below 1e-8 are clamped to exactly 0 (reconr.f90:4243-4244). Mode=11 doesn't do this.

20. **Mode=12 csunr2 potential scattering condition** — Uses `if (j.le.1)` (Fortran line 4249) vs mode=11's `if (j.eq.1)`. The `.le.` means potential scattering is added for j=0 AND j=1, not just j=1.

21. **Mode=12 final XS interpolation uses per-J INT** — The INT variable in csunr2 is overwritten by each J-state (line 4219). The final terp1 call at line 4312 uses the LAST J-state's INT, typically 2 (linear). Mode=11 hardcodes log-log (intlaw=5).

---

## Immediate Next Steps

### Priority 1: Debug Test 07 grid diffs (SLBW resolved range)

The 114 grid diffs start at 6.3 eV (inside SLBW range 1-82 eV). The energies differ by tiny amounts (6.2953 vs 6.295313), suggesting the initial MF2 peak nodes differ.

**Debugging approach:**
1. Use a subagent to read Fortran `rdf2bw` (lines 884-1310) — specifically how it computes `ndig` for peak triplets and adds nodes to enode
2. Use another subagent to extract Julia's `_add_bw_peaks!` output for U-235 (print all nodes added)
3. Compare node-by-node: are the resonance peak triplets at the same energies?
4. The Fortran rdf2bw uses `ndig = clamp(2 + log10(Er/hw), 5, 9)` — verify the Julia matches

### Priority 2: Verify mode=12 URR XS accuracy

Once the grid matches, check if the URR cross sections are correct. Compare the Julia URR table output (from `build_unresolved_table`) with the Fortran's at a few specific energies.

### Priority 3: Grind more tests

After Test 07, try other tests with MLBW, different materials. Each new test exercises different code paths.

### Priority 4: Grind BROADR

The oracle cache has `after_broadr.pendf` for Tests 01, 02, 07. This is the next module in the processing chain (Doppler broadening).

---

## Test Details

| Test | MAT | ENDF file | err | Formalism | Chain | Status |
|------|-----|-----------|-----|-----------|-------|--------|
| 84 | 128 | n-001_H_002-ENDF8.0.endf | 0.001 | LRU=0 | RECONR only | **BIT-IDENTICAL** |
| 01 | 1306 | t511 | 0.005 | LRU=0 | RECONR→BROADR→... | RECONR **BIT-IDENTICAL** |
| 02 | 1050 | t404 | 0.005 | SLBW+URR(mode=11) | RECONR→BROADR→UNRESR→... | RECONR **BIT-IDENTICAL** |
| 08 | 2834 | eni61 | 0.01 | Reich-Moore(LRF=3) | RECONR→BROADR→HEATR→... | RECONR **BIT-IDENTICAL** |
| 07 | 1395 | t511 | 0.005 | SLBW+URR(mode=12) | RECONR→BROADR→... | **0/27** (runs, grid diffs in resolved range) |

---

## Session History

### Phase 1-2: LRU=0 materials (Tests 84, 01) → 33/33 PERFECT

Built `lunion_grid`, fixed 17 bugs in grid construction, XS evaluation, threshold handling, PENDF writing.

### Phase 3: SLBW resolved (Test 02) → 13/17 PERFECT

Restructured LRU=1 pipeline to use lunion_grid. Fixed threshold interpolation, 9-sigfig ENDF float format, float accumulation order.

### Phase 4: Unresolved resonances + remaining fixes → 14/17 PERFECT

Implemented complete LRU=2 mode=11 evaluation: `_csunr1`, `_gnrl`, `build_unresolved_table`, `eval_unresolved`. SLBW sigfig rounding, partials-only convergence, boundary-crossing forced convergence.

### Phase 5: Test 02 final 3 diffs → 17/17 BIT-IDENTICAL

Fixed MT=1 total accumulation (round each channel before summing), MT=102 rounding boundary, MT=2 URR boundary node.

### Phase 6: Test 08 complete → 18/18 BIT-IDENTICAL

**Key breakthrough**: discovered that Fortran lunion processes MF=12 LO=1 sections (photon multiplicities) alongside MF=3. MF=12 MT=102 for Ni-61 uses histogram interpolation (INT=1), triggering the two-pass shading mechanism. This created {sigfig(E,7,-1), sigfig(E,7,+1)} pairs at 8 breakpoints (0.0253, 4000, 100000, 250000, 500000, 750000, 2000000, 5000000 eV) — explaining all 41 missing grid points. Previous analysis incorrectly attributed this to the coincidence check at line 1996 (which never fires). Also added pseudo-threshold skip for redundant MT=4.

### Phase 7: Test 07 scaffolding → mode=12 reader + evaluator implemented

Implemented `_read_urr_lrf2`, `_csunr2`, `URR2Sequence`, `URR2Data`. Fixed `elim = min(0.99e6, eresr)`. Test runs but 0/27 — grid diffs cascade from initial MF2 peak nodes in the dense 1-82 eV SLBW resolved range.
