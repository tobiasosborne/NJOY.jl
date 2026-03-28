# NJOY.jl Acceptance Criteria & Validation Standards

## Our Standards

**Stretch goal**: 1e-9 relative tolerance (bit-identical with Fortran NJOY2016 at ENDF precision).
This is achievable for RECONR (19 tests already pass) and is pursued relentlessly for all modules.

**First-round acceptance**: 1e-7 relative tolerance on all float values, plus structural match
(same line counts, same sections, same energy grid sizes).

**Minimum viable**: 1e-4 relative tolerance (0.01%), matching the NJOY21 C++ reimplementation
standard. Anything worse than this is a bug.

**Per-test tolerance**: Each test's RECONR `err` parameter (typically 0.001 = 0.1%) defines the
physical accuracy target. Values agreeing within `err` are physically identical; disagreement
beyond `err` indicates an algorithm error.

---

## What the NJOY Maintainers Say

### The 1e-9 tolerance is compiler-specific

**GitHub Issue #229** (Jean-Christophe Sublet, 2022-02-03):
https://github.com/njoy/NJOY2016/issues/229

> "Traditionally NJOY87, NJOY99, NJOY12 have always been supported on multiple OS & Fortran
> compiler. If the free Linux as OS and GNU as Fortran certainly took over one must also
> recognise that the Intel Ifort compiler played an important role. The quality and thoroughness
> of the regression test suite of NJOY2016 makes it compiler specific. **NJOY2016.65 when compiled
> with ifort version 19 and above delivers the same physics as its GNU counterpart but 25% of
> the test fail on a too strict numerical precision criteria** imposed when compared with GNU
> Fortran reference output. The opposite would also be true. It would make sense in my view for
> the regression test suite reference outputs to be compiler specific and/or the numerical
> precision relaxed. Historically, Bob always relied on 64-bit SUN then Intel compilers to
> program NJOY capabilities."

**Wim Haeck** (NJOY maintainer, LANL), responding on 2022-10-05:
https://github.com/njoy/NJOY2016/issues/229

> "arm64 architecture using gcc11.3 (which does not have the ERRORR compiler error issue for
> some reason) **also fail on multiple tests - due to precision settings**. I got a new laptop at
> work and it is an M1 (i.e. arm64 architecture). I'll play around with the precision settings
> in the comparisons to see if I can alleviate the issue a bit."

### Same Fortran, different arch: 30% failure rate

**GitHub Discussion #322** (YaqiWang, 2023-11-28):
https://github.com/njoy/NJOY2016/discussions/322

On Apple M1 Ultra with gfortran 11.0.1: **24 of 78 tests failed** (69% pass rate).
YaqiWang loosened tolerance to 1e-5 and most passed, but Tests 48 and 59 still failed
due to **different line counts** (different number of energy grid points generated).

### Haeck's definitive tolerance analysis

**GitHub PR #140** (Wim Haeck, 2021-01-06):
https://github.com/njoy/NJOY2016/pull/140

When updating physical constants to CODATA2018, Haeck performed the most detailed tolerance
analysis in the NJOY project history:

> "To perform these tests, I updated the execute.py for these tests, now using Python's own
> isclose function to compare two numbers. By default, the relative difference must be greater
> than 1e-9 for two numbers to be flagged as different, as long as one of those is not zero.
> In that case, the absolute difference may be up to 1e-10 for these numbers to be considered
> equal. **At those settings, 47 out of 58 tests are flagged as failed.**"

> "At a relative difference of 1e-9, **ENDF formatted numbers (with 6 digits after the decimal
> point) that differ even in their last digit get flagged as different.** Setting the relative
> difference to 1e-5 and 5e-5 will allow us to filter out the cases where the difference is in
> the last significant digit in standard ENDF notation. **At 1e-5 relative difference, only 14
> out of 58 tests still continue to be flagged as 'failed'.**"

Haeck classified those 14 remaining failures:

- Tests 12, 14 (VIEWR PostScript): *"only differences in VIEWR produced post-script files --
  **acceptable**"*
- Test 1 (format difference): *"A single reference tape has ~500 different lines due to a
  difference in fixed and scientific notation (the numbers are actually the same) --
  **acceptable**"*
- Tests 2, 11, 16 (scattering matrices): *"small differences in scattering matrices --
  **acceptable**. The differences are detected in GENDF files... These constitute differences
  in the 4th or 5th digit"*
- Tests 19, 37, 39 (probability tables): *"probability tables are slightly different --
  **acceptable**"*
- Tests 21, 32 (line counts): *"output files appear to have different number of lines --
  **might require further investigation**. Most likely some tolerance may have lead to the
  addition of an additional point or something."*
- Tests 22, 25, 49 (thermal scattering): *"These are mainly in the 4th and 5th digit after
  the decimal point"*

### Last-digit differences are "not physically meaningful"

**GitHub Issue #309** (Sublet, 2023-09-22):
https://github.com/njoy/NJOY2016/issues/309

With NJOY2016.71 and Intel ifort: 76% tests passed, 19 failed. Sublet's assessment:

> "Comments quickly brushed aside if one takes the step to look at the Artificial Intelligence
> style reason(s) the tests failed by browsing tests/\*/\*\_diff: `-0.00` or `0.00` in a PDF
> file!! ... or comparing the file production dates !! or claiming that **the numerical
> difference between 0.94739773 and 0.94739774 is physically meaningful, scientifically
> important.**"

> "**Bottom line: it is safe and robust to compile NJOY2016 with Intel oneAPI ifort**, even if
> some of the test cases are said to failed because the reference output tapes have been
> produced using NJOY2016 compiled with gcc on another OS brand."

### Line count differences are the real concern

**GitHub PR #99** (Wim Haeck and Jeremy Conlin, 2018-09-11):
https://github.com/njoy/NJOY2016/pull/99

> **jlconlin**: "I tested with 1E-5 relative tolerance. This didn't change the test results.
> When the comparison is done, the first thing that is done is to count the number of lines in
> the generated file and the reference file. **If the number of lines are not the same, then no
> other checking is done.**"

> **whaeck**: "Well, we do test non regression, so **if the number of lines is different there is
> definitely something going on.** In those cases a by hand verification is warranted."

---

## Published Standards for Processing Code Reimplementations

### NJOY21 V&V (the C++ reimplementation)

**LA-UR-20-27056** (Haeck, Conlin, Parsons, 2020), cited in our `VALIDATION_PLAN.md`:

| Metric | Acceptance |
|--------|------------|
| PENDF pointwise cross sections | relative difference < 0.01% (100 ppm) |
| ACE total cross section at thermal | relative difference < 0.001% (10 ppm) |
| Multigroup cross sections | relative difference < 0.1% (1000 ppm) |
| k_eff impact | < 10 pcm for all ICSBEP benchmarks tested |

This establishes precedent that a reimplementation achieves < 0.01% on PENDF cross sections.

### ACE file comparison (LANL standard practice)

**LA-UR-13-21822** (Conlin, 2013):

| Metric | Tolerance |
|--------|-----------|
| Cross section data | relative tolerance 0.01% |
| Angular/energy distributions | relative tolerance 0.1% |
| MT 1 total at thermal (0.0253 eV) | must agree within 0.001% |
| MT 2 elastic at thermal | must agree within 0.001% |
| MT 18 fission at thermal | must agree within 0.01% |
| MT 102 capture at thermal | must agree within 0.01% |

### Inter-code agreement (NJOY vs PREPRO vs AMPX)

**D.E. Cullen, "Accuracy of processed nuclear data"**:
https://www.researchgate.net/publication/236539541

> "AMPX, NJOY and PREPRO independently each define their own energy grid to produce results
> that **agree with each other to within 0.2%.**"

> "NJOY 9 digit results show results similar to the PREPRO 9 digit results, with **maximum
> differences of less about 0.01%**, an order of magnitude below the target 0.1%."

> "These differences are not due to any changes WITHIN these codes, but rather solely due to
> changes in their output to the ENDF format."

### FRENDY vs NJOY (JAEA)

**Journal of Nuclear Science and Technology**, 2021:
https://www.tandfonline.com/doi/full/10.1080/00223131.2021.1921631

> "The relative differences between FRENDY/MG and NJOY are generally small (**10^-3 -- 10^-4**)
> for LWR, FR, and 1/E spectra."

### AXSP vs NJOY2016 (Xi'an Jiaotong University)

**Frontiers in Energy Research**, 2022:
https://www.frontiersin.org/journals/energy-research/articles/10.3389/fenrg.2022.1009515

> "the relative error of multigroup flux and multigroup cross sections of various reaction
> types are **less than 0.01%**."

> "The multiplication factor difference between AXSP and NJOY2016 is **less than 20 pcm**."

---

## The Official Test Infrastructure

### execute.py (the sole test runner)

**File**: `njoy-reference/tests/execute.py`

Hardcoded tolerances (no per-test override):
```python
identical = identicalLines(reference_lines, trial_lines, diff_file, 1E-9, 1e-10)
```

Pass/fail logic:
1. Run NJOY, produce output tapes
2. For each `referenceTapeNN`, compare against `tapeNN`
3. First check: byte-identical (`filecmp.cmp`)
4. If not identical: parse floats, compare with `math.isclose(rel_tol=1e-9, abs_tol=1e-10)`
5. Date patterns (`MM/DD/YY`) masked before comparison
6. **Different line counts = immediate failure** (no value comparison attempted)

### ENDF format precision

From the ENDF-102 manual (CSEWG Document ENDF-102, BNL-90365-2009):

Standard ENDF notation provides **6-7 significant figures** (e.g., `1.234567+3`).
The 9-digit fixed-point extension (e.g., `12345.6789`) provides up to 9 significant figures.

The ENDF-102 specification itself defines NO numerical tolerance for processing codes.
It defines format rules, sum rules (MT=1 = MT=2 + MT=3 + ...), energy range requirements
(1e-5 to 20 MeV), and interpolation laws. Validation is delegated to CSEWG testing procedures.

### Compiler flags affecting reproducibility

From `njoy-reference/CMakeLists.txt`:
```
set( njoy_GNU_Linux_RELEASE_flags "-O3" "-DNDEBUG" )
set( njoy_GNU_Linux_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-frounding-math" "-fsignaling-nans" )
```

The `-O3` flag allows the compiler to reorder FP operations, use FMA instructions, and
vectorize loops. This means the **reference tapes are tied to a specific compiler version
and architecture**. The README warns that gcc-11.3 specifically fails to compile NJOY2016.

---

## Summary: Tolerance Hierarchy

| Level | Tolerance | Who uses it | Status for NJOY.jl |
|-------|-----------|-------------|---------------------|
| **Stretch goal** | **1e-9** | Fortran same-binary regression | Achieved for 19 RECONR tests. Pursued relentlessly. |
| **First-round acceptance** | **1e-7** | Last digit of 7-sigfig ENDF format | Achievable for most modules after structural fix. |
| **NJOY21 reimplementation** | 1e-4 (0.01%) | LANL C++ port (LA-UR-20-27056) | Already achieved for RECONR+BROADR. |
| **Inter-code agreement** | 1e-3 (0.1%) | NJOY/PREPRO/AMPX comparison | Already achieved for full pipeline. |
| **Physical accuracy** | `err` (0.001-0.005) | RECONR reconstruction tolerance | Already achieved. |
| **Integral validation** | 10-50 pcm | ICSBEP k_eff benchmarks | Not yet tested. |

**The structural requirement (same line counts) is non-negotiable at ANY tolerance level.**
The maintainers explicitly flag line count differences as requiring investigation, even when
value differences are accepted as trivial.
