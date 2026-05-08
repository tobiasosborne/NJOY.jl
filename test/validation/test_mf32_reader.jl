# =========================================================================
# MF=32 (resonance-parameter covariance) reader tests
#
# Phase 72 (2026-05-08): Foundation tests for `read_mf32`. Validates the
# reader against the actual U-238 JENDL-3.3 ENDF tape used by T15 (the
# canonical first target for the rescon port).
#
# Expected values pulled directly from the ENDF text at
# `njoy-reference/tests/resources/J33U238` lines 8345-8359 (the file
# header + first MF=32/MT=151 isotope/range/subsection records). Cross-
# checked against the Fortran port: rpxlc12 (errorr.f90:4252-4274) reads
# the same LIST and slices `mpar = nint(a(lptr+2))`, `nrb = nint(a(lptr+5))`,
# with body length `6*nrb + npar*(npar+1)/2` where `npar = mpar*nrb`.
# =========================================================================

using NJOY: read_mf32, mf32_resonance_count, mf32_param_count,
            MF32Data, MF32IsotopeData, MF32ResolvedRange,
            MF32ResolvedSubsection
using Test
using LinearAlgebra: diag

const T15_ENDF = joinpath(@__DIR__, "..", "..",
                          "njoy-reference", "tests", "resources", "J33U238")
const T15_MAT = 9237  # U-238 in JENDL-3.3

# -------------------------------------------------------------------------
# Top-level: file-level header + isotope count
# -------------------------------------------------------------------------

@testset "read_mf32 — U-238 JENDL-3.3 top-level header" begin
    isfile(T15_ENDF) || (@warn "skipping: $T15_ENDF missing"; return)
    data = read_mf32(T15_ENDF, T15_MAT)

    # First MF=32 line: ZA=92238.0, AWR=236.006, NIS=1
    @test data.za ≈ 92238.0
    @test data.awr ≈ 236.006
    @test length(data.isotopes) == 1

    iso = data.isotopes[1]
    # Second MF=32 line: ZAI=92238.0, ABN=1.0, LFW=0, NER=11.
    # NER counts ALL ranges (10 resolved + 1 URR); resolved_ranges
    # holds only LRU=1 (URR is skipped — port deferred).
    @test iso.zai ≈ 92238.0
    @test iso.abn ≈ 1.0
    @test iso.lfw == 0
    @test length(iso.resolved_ranges) == 10
end

# -------------------------------------------------------------------------
# Range 1: [1e-5, 1000] eV — the resolved range that drives MT=102 group 1
# -------------------------------------------------------------------------

@testset "read_mf32 — U-238 range 1 (1e-5..1000 eV, RM/LCOMP=1)" begin
    isfile(T15_ENDF) || (@warn "skipping: $T15_ENDF missing"; return)
    data = read_mf32(T15_ENDF, T15_MAT)
    rng = data.isotopes[1].resolved_ranges[1]

    # Range header (line 8347): ELR=1e-5, EHR=1000, LRU=1, LRF=3, NRO=0, NAPS=0
    @test rng.elr ≈ 1.0e-5
    @test rng.ehr ≈ 1000.0
    @test rng.lrf == 3      # Reich-Moore
    @test rng.naps == 0

    # Per-format CONT (line 8348): SPI=0.0, AP=0.942848, LCOMP=1, NLS=0, ISR=0
    @test rng.spi ≈ 0.0
    @test rng.ap ≈ 0.942848
    @test rng.nls == 0
    @test rng.isr == 0
    @test rng.dap ≈ 0.0

    # NSRS=1 sub-section (line 8349 CONT: AWRI=236.006, NSRS=1, NLRS=0).
    @test length(rng.subsections) == 1
    sub = rng.subsections[1]
    @test sub.awri ≈ 236.006
    @test sub.mpar == 3       # ER, GN, GG (U-238 RM has no fission)
    @test sub.nrb == 26       # 26 resonances in [1e-5, 1000] eV
    @test size(sub.params) == (26, 6)
    @test size(sub.cov) == (78, 78)

    # First resonance (line 8351): ER=6.674, AJ=0.5, GN=1.493e-3,
    # GG=2.300e-2, GFA=0.0, GFB=9.99e-9. ENDF stores RM as
    # (ER, AJ, GN, GG, GFA, GFB) at indices 1..6.
    @test sub.params[1, 1] ≈ 6.674       atol=1e-6
    @test sub.params[1, 2] ≈ 0.5         atol=1e-6
    @test sub.params[1, 3] ≈ 1.493e-3    atol=1e-9
    @test sub.params[1, 4] ≈ 2.300e-2    atol=1e-9
    @test sub.params[1, 5] ≈ 0.0         atol=0
    @test sub.params[1, 6] ≈ 9.99e-9     atol=1e-12
end

# -------------------------------------------------------------------------
# Covariance-matrix shape and basic invariants
# -------------------------------------------------------------------------

@testset "read_mf32 — covariance symmetry and structure" begin
    isfile(T15_ENDF) || (@warn "skipping: $T15_ENDF missing"; return)
    data = read_mf32(T15_ENDF, T15_MAT)

    for (ir, rng) in enumerate(data.isotopes[1].resolved_ranges)
        for (is, sub) in enumerate(rng.subsections)
            # Symmetric (we filled both triangles from upper packed).
            @test sub.cov ≈ sub.cov'  rtol=0
            # Diagonal must be ≥ 0 — variances cannot be negative.
            @test all(>=(0.0), diag(sub.cov))  # imported below via `using LinearAlgebra`
        end
    end
end

# -------------------------------------------------------------------------
# Aggregate counters used by the sensitivity builder for sizing
# -------------------------------------------------------------------------

@testset "read_mf32 — aggregate counts" begin
    isfile(T15_ENDF) || (@warn "skipping: $T15_ENDF missing"; return)
    data = read_mf32(T15_ENDF, T15_MAT)

    # Sum of NRB over the 11 ranges. We verify the first range (NRB=26)
    # exactly; the totals just need to be self-consistent.
    @test mf32_resonance_count(data) >= 26
    # Each resonance contributes MPAR=3 uncertain parameters.
    @test mf32_param_count(data) == 3 * mf32_resonance_count(data)
end

# -------------------------------------------------------------------------
# Range 2 sanity (1000..2000 eV, NRB=28) — proves multi-range parsing
# -------------------------------------------------------------------------

@testset "read_mf32 — range 2 (1000..2000 eV)" begin
    isfile(T15_ENDF) || (@warn "skipping: $T15_ENDF missing"; return)
    data = read_mf32(T15_ENDF, T15_MAT)
    @test length(data.isotopes[1].resolved_ranges) >= 2
    rng2 = data.isotopes[1].resolved_ranges[2]
    @test rng2.elr ≈ 1000.0
    @test rng2.ehr ≈ 2000.0
    @test rng2.lrf == 3
    @test length(rng2.subsections) == 1
    @test rng2.subsections[1].nrb == 28
    # Range 2 first resonance: ER=1054.65, AJ=0.5 (per ENDF line 8895).
    @test rng2.subsections[1].params[1, 1] ≈ 1054.65  atol=1e-3
end
