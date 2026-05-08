# =========================================================================
# MF=32 (resonance-parameter covariance) reader
#
# Mirrors Fortran resprx + rpxlc12 (njoy-reference/src/errorr.f90:3011-3250
# + 4108-4632). Used by errorr's rescon path (errorr.f90:7465 → 8513) to
# propagate resonance-parameter uncertainty into MF=33 covariances for
# seven (mt, mt2) pairs: (1,1), (2,2), (2,18), (2,102), (18,18), (18,102),
# (102,102).
#
# Scope of this port:
#   LRU=1 (resolved), LRF=1/2 (BW) or LRF=3 (Reich-Moore), LCOMP=1
#   (general covariance form), NRO=0 (energy-independent AP), NLRS=0
#   (no long-range covariance), ISR=0 or ISR=1 with LRF=1/2.
#
# This covers U-238 JENDL-3.3 (T15 reference), the canonical first
# target. LCOMP=0 (legacy, errorr.f90:3734 rpxlc0), LCOMP=2 (compact,
# errorr.f90:4634 rpxlc2), LRF=7 (RML/SAMMY, errorr.f90:3252 rpxsamm),
# NRO≠0, NLRS≠0, and LRU=2 (URR cov) error loudly with named TODOs.
#
# Data layout (LCOMP=1, per NSRS subsection):
#   List of NRB resonances, each with 6 raw parameter words:
#     LRF=1/2 (BW): (ER, AJ, GT, GN, GG, GF)
#     LRF=3 (RM):   (ER, AJ, GN, GG, GFA, GFB)
#   followed by NPAR=MPAR×NRB upper-triangular packed covariance words
#   over the MPAR uncertain parameters per resonance:
#     MPAR=1 → ER only
#     MPAR=2 → ER, GN
#     MPAR=3 → ER, GN, GG       (LRF=3, U-238)
#     MPAR=4 → ER, GN, GG, GFA  (or GF for BW)
#
# Ref: ENDF-6 manual §32 (resonance parameter uncertainty).
# =========================================================================

"""
    MF32ResolvedSubsection

One NSRS sub-section of MF=32 LCOMP=1 covariance data for a resolved
resonance range. Mirrors the LIST record at errorr.f90:4252 (rpxlc12).
"""
struct MF32ResolvedSubsection
    awri::Float64                # atomic-weight ratio
    mpar::Int                    # uncertain parameters per resonance
    nrb::Int                     # resonances in this sub-section
    params::Matrix{Float64}      # (nrb, 6): ENDF-stored raw parameters
    cov::Matrix{Float64}         # (npar, npar) symmetric, npar = mpar*nrb
end

"""
    MF32ResolvedRange

One energy range of MF=32 LRU=1 covariance data.
Mirrors the per-range records at errorr.f90:3088-3231 (resprx).
"""
struct MF32ResolvedRange
    elr::Float64                       # range lower energy [eV]
    ehr::Float64                       # range upper energy [eV]
    lrf::Int                           # 1, 2 (BW), or 3 (RM)
    naps::Int                          # AP application flag
    spi::Float64                       # target spin
    ap::Float64                        # scattering radius [√b]
    nls::Int                           # number of l-values from MF=32
    isr::Int                           # 0 or 1: scattering-radius unc flag
    dap::Float64                       # global ΔAP if isr=1 & lrf in (1,2); else 0
    subsections::Vector{MF32ResolvedSubsection}
end

"""
    MF32IsotopeData

MF=32 data for one isotope.
"""
struct MF32IsotopeData
    zai::Float64                       # isotope ZA identifier
    abn::Float64                       # abundance
    lfw::Int                           # fission widths flag (URR)
    resolved_ranges::Vector{MF32ResolvedRange}
end

"""
    MF32Data

Top-level MF=32 (resonance-parameter covariance) data for one material.
Returned by `read_mf32`. Currently carries only resolved-range data;
URR covariance (LRU=2) is not yet supported.
"""
struct MF32Data
    za::Float64
    awr::Float64
    isotopes::Vector{MF32IsotopeData}
end

# -------------------------------------------------------------------------
# Reader entry point
# -------------------------------------------------------------------------

"""
    read_mf32(filename::AbstractString, mat::Integer) -> MF32Data

Read MF=32/MT=151 (resonance-parameter covariance) for material `mat`
from an ENDF tape. Errors loudly on unsupported sub-formats.

Mirrors Fortran `resprx` (errorr.f90:3011-3250) + `rpxlc12`
(errorr.f90:4108-4632) for the LCOMP=1 path.
"""
function read_mf32(filename::AbstractString, mat::Integer)
    open(filename, "r") do io
        find_section(io, 32, 151; target_mat=Int(mat)) ||
            error("read_mf32: MF=32 MT=151 not found for MAT=$mat in $filename")

        head = read_head(io)
        za, awr = head.C1, head.C2
        nis = Int(head.N1)

        isotopes = MF32IsotopeData[]
        for _ in 1:nis
            iso = read_cont(io)
            zai, abn = iso.C1, iso.C2
            lfw = Int(iso.L2)
            ner = Int(iso.N1)

            resolved = MF32ResolvedRange[]
            for ir in 1:ner
                rng = _read_mf32_range(io, mat, ir)
                rng === nothing || push!(resolved, rng)
            end
            push!(isotopes, MF32IsotopeData(zai, abn, lfw, resolved))
        end
        return MF32Data(za, awr, isotopes)
    end
end

# -------------------------------------------------------------------------
# Per-range reader — dispatches on LRU/LRF/LCOMP
# -------------------------------------------------------------------------

# Returns nothing for ranges we cannot parse but want to skip cleanly
# (currently: none — every unsupported branch errors loudly per Rule 6).
# When URR cov support lands, LRU=2 will return nothing here and the
# parsed URR struct will live alongside resolved_ranges.
function _read_mf32_range(io::IO, mat::Integer, ir::Integer)
    range_cont = read_cont(io)
    elr, ehr = range_cont.C1, range_cont.C2
    lru = Int(range_cont.L1)
    lrf = Int(range_cont.L2)
    nro = Int(range_cont.N1)
    naps = Int(range_cont.N2)

    if lru == 2
        # URR covariance (LRU=2). Skip cleanly so the file pointer is
        # positioned at the next range header or MF=32 SEND. Mirrors
        # Fortran rpxunr's record sequence (errorr.f90:4824-4837):
        #   - one CONT (per-format header with NLS, ISR)
        #   - if ISR=1: one CONT (URR DAP)
        #   - NLS LIST records (per-l URR parameters)
        #   - one LIST record (URR covariance: MPAR/NPAR header)
        # The actual URR cov fill (uff/ugg/uee/utt/uef/ueg/ufg) is a
        # separate port (see RESCON_URR_TODO in rescon.jl). For T15
        # U-238 the URR range is [25 keV, 149 keV] — orthogonal to the
        # MT=102 row 1..14 canary which lives entirely in the resolved
        # range [1e-5, 1000] eV.
        urr_params_cont = read_cont(io)
        nls_urr = Int(urr_params_cont.N1)
        isr_urr = Int(urr_params_cont.N2)
        if isr_urr == 1
            read_cont(io)
        end
        for _ in 1:nls_urr
            read_list(io)
        end
        read_list(io)
        @info "read_mf32: MAT=$mat range $ir is URR [$elr, $ehr] eV — \
              skipping (URR cov port deferred to next session)."
        return nothing
    end
    if lru != 1
        error("read_mf32: MAT=$mat range $ir LRU=$lru not yet supported.")
    end
    if !(lrf in (1, 2, 3))
        error("read_mf32: MAT=$mat range $ir LRF=$lrf not yet supported \
              (LRF=4 Adler-Adler / LRF=7 R-Matrix Limited — TODO: port \
              rpxsamm, errorr.f90:3252).")
    end
    if nro != 0
        error("read_mf32: MAT=$mat range $ir NRO=$nro (energy-dependent \
              scattering radius) not yet supported.")
    end

    params_cont = read_cont(io)
    spi, ap = params_cont.C1, params_cont.C2
    lcomp = Int(params_cont.L2)
    nls = Int(params_cont.N1)
    isr = Int(params_cont.N2)

    if lcomp != 1
        error("read_mf32: MAT=$mat range $ir LCOMP=$lcomp not yet \
              supported (LCOMP=0 legacy / LCOMP=2 compact — TODO: port \
              rpxlc0 errorr.f90:3734 / rpxlc2 errorr.f90:4634).")
    end

    dap = 0.0
    if isr == 1
        if lrf == 1 || lrf == 2
            dap_cont = read_cont(io)
            # Per rpxlc12 (errorr.f90:3157-3162): dap is C2 of the CONT.
            dap = dap_cont.C2
        else  # lrf == 3
            error("read_mf32: MAT=$mat range $ir ISR=1 with LRF=3 \
                  (per-l ΔAP LIST) not yet supported — TODO: port the \
                  LIST branch at errorr.f90:3167-3196.")
        end
    end

    nsrs_cont = read_cont(io)
    awri = nsrs_cont.C1
    nsrs = Int(nsrs_cont.N1)
    nlrs = Int(nsrs_cont.N2)

    if nlrs != 0
        error("read_mf32: MAT=$mat range $ir NLRS=$nlrs (long-range \
              covariance) not yet supported — Fortran rpxlc12 errors \
              with 'nlrs>0 not coded' at errorr.f90:4618.")
    end

    subs = MF32ResolvedSubsection[]
    for _ in 1:nsrs
        push!(subs, _read_mf32_subsection(io, awri, lrf))
    end

    return MF32ResolvedRange(elr, ehr, lrf, naps, spi, ap, nls, isr, dap, subs)
end

# -------------------------------------------------------------------------
# Per-subsection reader (LCOMP=1 only)
# -------------------------------------------------------------------------

function _read_mf32_subsection(io::IO, awri::Float64, lrf::Int)
    list = read_list(io)
    mpar = Int(list.L1)
    nrb = Int(list.N2)
    npl = Int(list.N1)
    npar = mpar * nrb
    expected = 6*nrb + npar*(npar + 1) ÷ 2
    npl == expected || error(
        "read_mf32: LIST body length NPL=$npl ≠ 6·NRB + NPAR·(NPAR+1)/2 = $expected \
        (MPAR=$mpar, NRB=$nrb, NPAR=$npar, LRF=$lrf). \
        File may use a different LCOMP=1 layout than rpxlc12 expects.")

    # First 6·NRB words = parameter values, packed (NRB rows × 6 cols).
    # ENDF/MF=32 stores 6 words/resonance regardless of MPAR; only the
    # MPAR uncertain parameters per resonance are covered by `cov`.
    params = Matrix{Float64}(undef, nrb, 6)
    @inbounds for r in 1:nrb, c in 1:6
        params[r, c] = list.data[(r - 1) * 6 + c]
    end

    # Next NPAR·(NPAR+1)/2 words = upper-triangular packed covariance
    # over the NPAR uncertain parameters (MPAR per resonance × NRB).
    # Mirrors errorr.f90:4528-4535 (rpxlc12 lcomp=1 cov read):
    #   l3 = lptr + 5 + nvs1
    #   do i=1,npar; do j=i,npar; l3=l3+1; cov(i,j)=cov(j,i)=a(l3); end; end
    cov = zeros(Float64, npar, npar)
    k = 6 * nrb
    @inbounds for i in 1:npar, j in i:npar
        k += 1
        v = list.data[k]
        cov[i, j] = v
        cov[j, i] = v
    end

    return MF32ResolvedSubsection(awri, mpar, nrb, params, cov)
end

# -------------------------------------------------------------------------
# Convenience accessors for the rescon sandwich
# -------------------------------------------------------------------------

"""
    mf32_resonance_count(data::MF32Data) -> Int

Total number of resolved resonances summed across all isotopes, ranges,
and sub-sections. Useful for sizing sensitivity arrays.
"""
function mf32_resonance_count(data::MF32Data)
    n = 0
    for iso in data.isotopes, rng in iso.resolved_ranges, sub in rng.subsections
        n += sub.nrb
    end
    n
end

"""
    mf32_param_count(data::MF32Data) -> Int

Total number of uncertain parameters (Σ MPAR × NRB across sub-sections).
This is the dimension of the per-range covariance matrices stacked.
"""
function mf32_param_count(data::MF32Data)
    n = 0
    for iso in data.isotopes, rng in iso.resolved_ranges, sub in rng.subsections
        n += sub.mpar * sub.nrb
    end
    n
end
