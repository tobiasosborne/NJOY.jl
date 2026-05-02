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

    # Read MF=1/MT=451 from the GENDF tape to extract awr, ZA, fission flag,
    # and group/sigma_zero/temperature lists. Fortran wmsr's `wminit` reads
    # the same record (wimsr.f90:309-423).
    awr_endf, za_endf, fission_present, ntemp_tape = _wimsr_read_mf1_mt451(in_path, p.mat, p.mti)

    # Atomic weight per Fortran wimsr.f90:2042: awt = awr*amassn (CODATA neutron mass).
    amassn = 1.00866491595
    awt = awr_endf * amassn

    z_number = div(round(Int, za_endf), 1000)

    # Fortran fission flag (`ifis` in wimsr.f90):
    #   0 = no fission; 1 = nubar present, no resonance; 2 = MF=3/MT=18 fission;
    #   3 = MT=18 fission with resonance fission tables (T11 path).
    # For T11 (Pu-238) we have ires>0 + MT=18 → ifis=3. The proper value is
    # decided by xsecs once it scans the tape; here we infer a reasonable
    # default from the tape contents + deck flags.
    ifis = if fission_present && p.ires > 0
        3
    elseif fission_present
        2
    else
        0
    end

    nrestb = p.ires > 0  # has_resonance_tables
    nsigz_eff = p.nsigz == 0 ? 7 : p.nsigz  # default 7 for T11; xsecs scans tape later

    # Stub resonance tables — one per resonance group, all zeros. Replaced by
    # the resint port (Phase 58b/c). Keeps validate() happy so the upstream
    # WIMS sections (lines 1-921 of T11) are testable.
    res_stubs = nrestb ? [WIMSResonanceTable(p.rdfid, p.ires, nsigz_eff,
                                              zeros(p.ires), zeros(nsigz_eff),
                                              zeros(p.ires * nsigz_eff))
                          for _ in 1:p.nrg] : WIMSResonanceTable[]
    fres_stubs = (nrestb && ifis == 3) ? [WIMSResonanceTable(p.rdfid, p.ires, nsigz_eff,
                                                              zeros(p.ires), zeros(nsigz_eff),
                                                              zeros(p.ires * nsigz_eff))
                                          for _ in 1:p.nrg] : WIMSResonanceTable[]

    # Use the actual ntemp from the GENDF if the deck card-4 ntemp was 0
    # ("use all temperatures on tape" — wimsr.f90:407-408).
    final_ntemp = (p.ntemp == 0) ? ntemp_tape : p.ntemp
    final_ntemp = max(final_ntemp, 1)

    # Populate spot[ngr0..ngr1] = sigp (Fortran wimsr.f90:950: spot(i)=sigp).
    spot = fill(p.sigp, p.nrg)

    # Run xsecs port to extract MF=3 + MF=6 (P0) data.
    xs_data = wimsr_extract_xsecs(in_path, p)
    spot_w = xs_data.spot[p.nfg+1:p.nfg+p.nrg]
    sdp_w  = xs_data.sdp[p.nfg+1:p.nfg+p.nrg]
    xtr_w  = xs_data.xtr[1:nnt]
    ab0_w  = xs_data.ab0[1:nnt]

    # Build per-temperature thermal arrays from xsecs output. Each thermal
    # array slices the [nt0..ngnd] portion of the per-temp xtr_t / ab0_t / etc.
    # Fortran wimsr.f90:1595-1611 emits in this layout.
    nthermal = p.ngnd - nnt
    nt0 = nnt + 1
    n_temps_eff = min(length(xs_data.tempr), final_ntemp)
    therm_tr  = Vector{Float64}[]
    therm_abs = Vector{Float64}[]
    therm_nuf = Vector{Float64}[]
    therm_fis = Vector{Float64}[]
    therm_mat = Vector{Float64}[]
    for ti in 1:final_ntemp
        si = min(ti, n_temps_eff, length(xs_data.xtr_t))
        if si == 0
            push!(therm_tr,  zeros(nthermal))
            push!(therm_abs, zeros(nthermal))
            push!(therm_nuf, ifis > 1 ? zeros(nthermal) : Float64[])
            push!(therm_fis, ifis > 1 ? zeros(nthermal) : Float64[])
            push!(therm_mat, Float64[])
            continue
        end
        push!(therm_tr,  xs_data.xtr_t[si][nt0:p.ngnd])
        push!(therm_abs, xs_data.ab0_t[si][nt0:p.ngnd])
        if ifis > 1
            push!(therm_nuf, xs_data.snus_t[si][nt0:p.ngnd])
            push!(therm_fis, xs_data.sf0_t[si][nt0:p.ngnd])
        else
            push!(therm_nuf, Float64[])
            push!(therm_fis, Float64[])
        end
        # Thermal scatter matrix sparse flat (same format as nonthermal)
        flat = Float64[]
        for ia in nt0:p.ngnd
            lone = xs_data.l1_t[si][ia]
            ltwo = xs_data.l2_t[si][ia]
            if lone != p.ngnd || ltwo != 1
                push!(flat, Float64(ia - lone + 1))
                push!(flat, Float64(ltwo - lone + 1))
                if lone != ltwo
                    for j in lone:ltwo
                        push!(flat, xs_data.xs_t[si][ia, j])
                    end
                end
            else
                push!(flat, 1.0)
                push!(flat, 0.0)
            end
        end
        push!(therm_mat, flat)
    end

    # Build the nonthermal scattering matrix flat representation
    # (Fortran wimsr.f90:1533-1574 in xseco):
    #   for each nonthermal group `ia` (1..nnt):
    #     if l1[ia] != ngnd or l2[ia] != 1:  # has scattering
    #       k += (offset = ia - lone + 1)
    #       k += (nbands = ltwo - lone + 1)
    #       k += xs[ia, lone..ltwo]
    #     else:  # zero row
    #       k += (1, 0)
    nonthermal = Float64[]
    for ia in 1:nnt
        lone = xs_data.l1[ia]
        ltwo = xs_data.l2[ia]
        if lone != p.ngnd || ltwo != 1
            push!(nonthermal, Float64(ia - lone + 1))
            push!(nonthermal, Float64(ltwo - lone + 1))
            if lone != ltwo
                for j in lone:ltwo
                    push!(nonthermal, xs_data.xs[ia, j])
                end
            end
        else
            push!(nonthermal, 1.0)
            push!(nonthermal, 0.0)
        end
    end

    WIMSMaterial(
        p.iverw, p.nfid, p.rdfid,
        awt, z_number,
        p.ngnd, p.nfg, p.nrg,
        ifis,
        final_ntemp,
        nrestb,
        false,            # has_fission_spectrum — stub (isof=0 for T11 anyway)
        nothing,          # burnup — iburn=0 default; writer emits the (0.0, ident) pair
        lam,
        spot_w,           # potential_xs (length nrg) — from xsecs
        sdp_w,            # scatter_xs (sdp; length nrg) — from xsecs (=0 until MF=6 ports)
        xtr_w,            # transport_xs (xtr; length nnt) — from xsecs
        ab0_w,            # absorption_xs (ab0; length nnt) — from xsecs
        ifis > 1 ? xs_data.snu[1:nnt] .* xs_data.sf0[1:nnt] : Float64[],  # nu*fission
        ifis > 1 ? xs_data.sf0[1:nnt] : Float64[],  # fission_xs (inf-dil)
        nonthermal,       # nonthermal_matrix — packed sparse from xs[]/l1[]/l2[]
        xs_data.tempr[1:final_ntemp],
        therm_tr,
        therm_abs,
        therm_nuf,
        therm_fis,
        therm_mat,
        res_stubs,
        fres_stubs,
        Float64[],
        Vector{WIMSP1Block}[],
    )
end

"""
Read MF=1/MT=451 from a GENDF tape, returning `(awr, za, fission_present, ntemp_count)`.

Fortran wimsr's `wminit` (wimsr.f90:309-423) extracts equivalent metadata.
This is a minimal version; the full sigma-zero list / group bounds /
temperature list extraction belongs to the xsecs port.
"""
function _wimsr_read_mf1_mt451(in_path::AbstractString, target_mat::Int, mti::Int)
    awr = 0.0
    za = 0.0
    fission_present = false
    temps = Set{Int}()  # round to int, count distinct

    open(in_path, "r") do io
        for line in eachline(io)
            length(line) < 75 && continue
            p = rpad(line, 80)
            mat_val = tryparse(Int, strip(p[67:70]))
            mat_val === nothing && continue
            (target_mat != 0 && mat_val != target_mat && mat_val != 0) && continue
            mf = tryparse(Int, strip(p[71:72]))
            mt = tryparse(Int, strip(p[73:75]))
            seq = tryparse(Int, strip(p[76:80]))
            (mf === nothing || mt === nothing || seq === nothing) && continue

            # MF=1/MT=451 HEAD record (seq=1): ZA at cols 1-11, AWR at 12-22.
            if mat_val == target_mat && mf == 1 && mt == 451 && seq == 1 && awr == 0.0
                za  = parse_endf_float(p[1:11])
                awr = parse_endf_float(p[12:22])
            end

            # MF=1/MT=451 record 2 (seq=2): C1 = temperature.
            if mat_val == target_mat && mf == 1 && mt == 451 && seq == 2
                t = parse_endf_float(p[1:11])
                push!(temps, round(Int, t))
            end

            # MF=3/MT=18 presence flags fission.
            if mat_val == target_mat && mf == 3 && mt == 18 && seq == 1
                fission_present = true
            end
        end
    end
    return (awr, za, fission_present, length(temps))
end
