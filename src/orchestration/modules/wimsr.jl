# wimsr module runner — WIMS-D / WIMS-E reactor library generator.
#
# Reads a GENDF tape produced by groupr, extracts cross sections and
# scattering matrices for a single material, optionally builds resonance
# self-shielding tables and P1 components, then emits the WIMS-format
# library on the output unit.
#
# Ref: njoy-reference/src/wimsr.f90 — driver subroutine `wimsr` (lines 48-307);
# core processors `wminit` (309-423), `xsecs` (863-1422), `resint` (425-676),
# `p1scat` (1685-1935), `wimout` (1989-2147).
#
# Phase 58 (2026-05-02) — Phase A scaffolding only. Wires the dispatcher,
# parses the 8-card deck, and chains into the existing `src/formats/wimsr.jl`
# writer. The xsecs / resint / p1scat data extraction passes (tasks #3, #4)
# fill in the WIMSMaterial fields. Without them, this writes a stub WIMS
# tape that will not match referenceTape27 — that is the red bar covered by
# `test/validation/test_wimsr_t11_standalone.jl`.

"""
    wimsr_module(tapes::TapeManager, params::WimsrParams)

Build a WIMS library from a GENDF input tape.

Inputs:
  - `params.ngendf`: GENDF input unit (registered on `tapes`)
  - `params.nout`:   WIMS output unit (created in `tapes.work_dir`)

Drives the data extraction (xsecs / resint / p1scat) and serializes the
result via `write_wims(::IO, ::WIMSMaterial)`.

Skips silently when `nout == 0` (Fortran wmout no-op path, wimsr.f90:284).
"""
function wimsr_module(tapes::TapeManager, params::WimsrParams)
    @info "wimsr: ngendf=$(params.ngendf) nout=$(params.nout) " *
          "mat=$(params.mat) iverw=$(params.iverw) ngnd=$(params.ngnd) " *
          "nfg=$(params.nfg) nrg=$(params.nrg) " *
          "ntemp=$(params.ntemp) nsigz=$(params.nsigz) " *
          "ires=$(params.ires) ip1opt=$(params.ip1opt) iburn=$(params.iburn)"

    if params.nout == 0
        @warn "wimsr: nout=0 — nothing to write"
        return nothing
    end

    out_path = resolve(tapes, params.nout)

    if params.ngendf == 0
        @warn "wimsr: no GENDF input — writing empty stub tape"
        touch(out_path); register!(tapes, params.nout, out_path)
        return nothing
    end

    in_path = resolve(tapes, params.ngendf)
    if !isfile(in_path)
        @warn "wimsr: GENDF input (unit $(params.ngendf)) not found at $in_path " *
              "— writing empty stub tape"
        touch(out_path); register!(tapes, params.nout, out_path)
        return nothing
    end

    # Tasks #3-#5: build a real WIMSMaterial from the GENDF tape.
    # Phase 58 scaffolding stubs the data extraction — output will diff
    # referenceTape27 by ~all lines until the xsecs / resint port lands.
    wmat = _wimsr_build_material_stub(in_path, params)

    open(out_path, "w") do io
        write_wims(io, wmat)
    end
    register!(tapes, params.nout, out_path)
    nothing
end

"""
Stub WIMSMaterial builder — Phase 58 scaffolding only.

Produces a minimally-valid `WIMSMaterial` so the writer doesn't crash on
validate(). Real data extraction (Fortran `xsecs` + `resint`) replaces this
in tasks #3 and #4. Until then the output tape will not match Fortran's
referenceTape27.

When the real builders land, this stub is deleted and the orchestration
calls them directly.
"""
function _wimsr_build_material_stub(in_path::AbstractString, p::WimsrParams)
    nnt = p.nfg + p.nrg
    nthermal = p.ngnd - nnt
    ntemp = max(p.ntemp, 1)

    # Goldstein lambdas — already parsed; pad to nrg if deck supplied fewer.
    lam = if length(p.glam) >= p.nrg
        p.glam[1:p.nrg]
    else
        vcat(p.glam, ones(p.nrg - length(p.glam)))
    end

    WIMSMaterial(
        p.iverw, p.nfid, p.rdfid,
        0.0,         # atomic_weight — populated by xsecs
        0,           # z_number — populated by xsecs
        p.ngnd, p.nfg, p.nrg,
        0,           # fission_flag — populated by xsecs (1, 2, or 3)
        ntemp,
        false,       # has_resonance_tables — stub: data not yet extracted
        false,       # has_fission_spectrum — stub
        nothing,     # burnup — TODO when iburn>0
        lam,
        zeros(nnt),  # potential_xs
        zeros(nnt),  # scatter_xs
        zeros(nnt),  # transport_xs
        zeros(nnt),  # absorption_xs
        Float64[],   # nu_fission (filled if fission_flag>1)
        Float64[],   # fission_xs
        Float64[],   # nonthermal_matrix
        zeros(ntemp),         # temperatures
        [zeros(nthermal) for _ in 1:ntemp],  # thermal_transport
        [zeros(nthermal) for _ in 1:ntemp],  # thermal_absorption
        [Float64[] for _ in 1:ntemp],        # thermal_nu_fission
        [Float64[] for _ in 1:ntemp],        # thermal_fission
        [Float64[] for _ in 1:ntemp],        # thermal_matrices
        WIMSResonanceTable[],                # resonance_tables (filled by resint)
        WIMSResonanceTable[],                # fission_resonance_tables
        Float64[],                           # fission_spectrum
        Vector{WIMSP1Block}[],               # p1_data
    )
end
