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
            MF32ResolvedSubsection, MF32UnresolvedRange,
            MF32URRLState, MF32URRJState
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
    # NER counts ALL ranges (10 resolved + 1 URR). Phase 75 added URR
    # capture — both buckets must be populated.
    @test iso.zai ≈ 92238.0
    @test iso.abn ≈ 1.0
    @test iso.lfw == 0
    @test length(iso.resolved_ranges) == 10
    @test length(iso.unresolved_ranges) == 1
end

# -------------------------------------------------------------------------
# URR range 11 ([10 keV, 150 keV]) — Phase 75 LRU=2 capture
# -------------------------------------------------------------------------

@testset "read_mf32 — U-238 URR range (10 keV..150 keV)" begin
    isfile(T15_ENDF) || (@warn "skipping: $T15_ENDF missing"; return)
    data = read_mf32(T15_ENDF, T15_MAT)
    @test length(data.isotopes[1].unresolved_ranges) == 1
    urr = data.isotopes[1].unresolved_ranges[1]

    # Range head CONT (J33U238 line 16405): EL=1e4, EH=1.5e5, LRU=2, LRF=1
    @test urr.elr ≈ 1.0e4
    @test urr.ehr ≈ 1.5e5
    @test urr.lrf == 1

    # URR per-format CONT (line 16406): SPI=0, AP=0.94285, LSSF=0, NLS=3
    @test urr.spi ≈ 0.0
    @test urr.ap ≈ 0.94285  atol=1e-6
    @test urr.lssf == 0
    @test urr.nls == 3
    @test urr.isr == 0
    @test urr.dap ≈ 0.0

    # 3 L-states: L=0 (1 J), L=1 (2 J), L=2 (2 J).
    @test length(urr.l_states) == 3
    @test urr.l_states[1].l == 0 && length(urr.l_states[1].j_states) == 1
    @test urr.l_states[2].l == 1 && length(urr.l_states[2].j_states) == 2
    @test urr.l_states[3].l == 2 && length(urr.l_states[3].j_states) == 2

    # First J-state of L=0 (line 16408): D=11.014, AJ=0.5, GNO=1.0194e-3,
    # GG=2.3994e-2, GF=0, GX=0
    j00 = urr.l_states[1].j_states[1]
    @test j00.D ≈ 11.014       atol=1e-6
    @test j00.AJ ≈ 0.5
    @test j00.GNO ≈ 1.0194e-3  atol=1e-9
    @test j00.GG ≈ 2.3994e-2   atol=1e-9
    @test j00.GF ≈ 0.0
    @test j00.GX ≈ 0.0

    # Total J-states = 1+2+2 = 5; MPAR=3 (per cov LIST head 16415:
    # L1=MPAR=3, N1=NW=120, N2=NPAR=15). NPAR = MPAR · ΣNJS = 3·5 = 15.
    @test urr.mpar == 3
    @test size(urr.cov_RP) == (15, 15)
    # Cov is symmetric and diagonal must be ≥ 0.
    @test urr.cov_RP ≈ urr.cov_RP'  rtol=0
    @test all(>=(0.0), diag(urr.cov_RP))
    # First diagonal entry (ER variance) — first cov body word at
    # line 16416 col 1: 1.558500e-2.
    @test urr.cov_RP[1, 1] ≈ 1.5585e-2  atol=1e-6
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
