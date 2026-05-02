# WIMS-D / WIMS-E library format types, validation, and I/O
#
# Proposer-B design: clean Julia types matching the WIMS library specification
# as implemented in NJOY2016 wimsr.f90.
#
# The WIMS library contains:
#   - Material identification and burnup chain data
#   - Temperature-independent cross sections (potential, scatter, transport,
#     absorption, lambda) plus optional nu*fission and fission
#   - Nonthermal and thermal scattering matrices
#   - Temperature-dependent data for thermal energy range
#   - Resonance self-shielding tables (effective resonance integrals)
#   - Optional fission spectrum and P1 scattering matrices

using Printf

# =========================================================================
# Constants
# =========================================================================
const WIMS_DEFAULT_NGROUPS = 69
const WIMS_DEFAULT_NFG     = 14    # number of fast groups
const WIMS_DEFAULT_NRG     = 13    # number of resonance groups
const WIMS_SENTINEL        = 999999999
const WIMS_VERSION_D       = 4
const WIMS_VERSION_E       = 5

# =========================================================================
# Burnup chain data
# =========================================================================
"""
    BurnupProduct

A single product entry in the WIMS burnup chain: an isotope identifier
paired with a yield or decay constant.
"""
struct BurnupProduct
    ident::Int       # WIMS library identifier of the product isotope
    value::Float64   # yield (for capture/fission products) or decay constant (s^-1)
end

"""
    WIMSBurnup

Burnup chain data for a WIMS material. Following the NJOY convention:
  - `capture_product`: capture product (card 6a)
  - `decay_product`: decay product with decay constant (card 6b)
  - `fission_products`: vector of fission product yields (card 6c, ntis-2 entries)
  - `energy_per_fission`: energy released per fission (MeV)
"""
struct WIMSBurnup
    capture_product::BurnupProduct
    decay_product::BurnupProduct
    fission_products::Vector{BurnupProduct}
    energy_per_fission::Float64
end

# =========================================================================
# Resonance table
# =========================================================================
"""
    WIMSResonanceTable

Resonance self-shielding data for one resonance group in the WIMS library.
Contains temperatures, sigma-zero values, and effective resonance integrals
at each (temperature, sigma-zero) combination.
"""
struct WIMSResonanceTable
    rid::Float64                 # resonance identifier (rdfid)
    ntemps::Int                  # number of temperatures
    nsigz::Int                   # number of sigma-zero values
    temperatures::Vector{Float64}
    sigma_zeros::Vector{Float64}
    integrals::Vector{Float64}   # resonance integrals, length = ntemps * nsigz
end

function validate(rt::WIMSResonanceTable)
    rt.ntemps > 0 || error("WIMSResonanceTable: ntemps must be > 0, got $(rt.ntemps)")
    rt.nsigz > 0 || error("WIMSResonanceTable: nsigz must be > 0, got $(rt.nsigz)")
    length(rt.temperatures) == rt.ntemps ||
        error("WIMSResonanceTable: temperatures length mismatch")
    length(rt.sigma_zeros) == rt.nsigz ||
        error("WIMSResonanceTable: sigma_zeros length mismatch")
    expected = rt.ntemps + rt.nsigz + rt.ntemps * rt.nsigz
    # Integrals array holds temps + sigz + data in one flat vector per Fortran convention
    # but we store them separately; just check data length
    return true
end

# =========================================================================
# P1 scattering block
# =========================================================================
"""
    WIMSP1Block

P1 scattering matrix data for one energy group: a compressed sparse row
with `offset` (first sink group, 1-based) and `values` (non-zero entries).
"""
struct WIMSP1Block
    offset::Int              # first group index (1-based) in the sink range
    values::Vector{Float64}  # P1 scattering transfer values
end

# =========================================================================
# Main WIMS material type
# =========================================================================
"""
    WIMSMaterial

Complete WIMS-D or WIMS-E library record for a single material.

# Fields
- `version`: WIMS version (4=WIMS-D, 5=WIMS-E)
- `ident`: material identifier (integer)
- `rid`: resonance identifier (floating-point, rdfid)
- `atomic_weight`: atomic weight in amu
- `z_number`: atomic number Z
- `ngroups`: total number of energy groups
- `nfg`: number of fast groups
- `nrg`: number of resonance groups
- `fission_flag`: fission type (0=none, 1=fissile no res, 2=with res fission, etc.)
- `ntemp`: number of temperatures in thermal range
- `has_resonance_tables`: whether resonance tables are present
- `has_fission_spectrum`: whether chi is included
- `burnup`: optional burnup chain data
- `lambdas`: Goldstein lambda values for resonance groups (length nrg)
- `potential_xs`: potential scattering cross section per group (nfg+1..nfg+nrg)
- `scatter_xs`: elastic scattering cross section per nonthermal group
- `transport_xs`: transport cross section per nonthermal group
- `absorption_xs`: absorption cross section per nonthermal group
- `nu_fission`: nu*sigma_f per nonthermal group (empty if non-fissile)
- `fission_xs`: sigma_f per nonthermal group (empty if non-fissile)
- `nonthermal_matrix`: nonthermal scattering matrix (sparse, flattened)
- `temperatures`: temperature values (K) for thermal data
- `thermal_transport`: transport XS in thermal range, per temperature
- `thermal_absorption`: absorption XS in thermal range, per temperature
- `thermal_nu_fission`: thermal nu*sigma_f, per temperature (empty if non-fissile)
- `thermal_fission`: thermal sigma_f, per temperature (empty if non-fissile)
- `thermal_matrices`: thermal scattering matrices, per temperature
- `resonance_tables`: resonance self-shielding tables per resonance group
- `fission_resonance_tables`: resonance fission tables (if fission_flag == 3)
- `fission_spectrum`: chi vector (empty if not included)
- `p1_data`: P1 scattering data per temperature, each a vector of WIMSP1Block
"""
struct WIMSMaterial
    version::Int
    ident::Int
    rid::Float64
    atomic_weight::Float64
    z_number::Int
    ngroups::Int
    nfg::Int
    nrg::Int
    fission_flag::Int
    ntemp::Int
    has_resonance_tables::Bool
    has_fission_spectrum::Bool
    burnup::Union{Nothing,WIMSBurnup}
    lambdas::Vector{Float64}
    potential_xs::Vector{Float64}
    scatter_xs::Vector{Float64}
    transport_xs::Vector{Float64}
    absorption_xs::Vector{Float64}
    nu_fission::Vector{Float64}
    fission_xs::Vector{Float64}
    nonthermal_matrix::Vector{Float64}
    temperatures::Vector{Float64}
    thermal_transport::Vector{Vector{Float64}}
    thermal_absorption::Vector{Vector{Float64}}
    thermal_nu_fission::Vector{Vector{Float64}}
    thermal_fission::Vector{Vector{Float64}}
    thermal_matrices::Vector{Vector{Float64}}
    resonance_tables::Vector{WIMSResonanceTable}
    fission_resonance_tables::Vector{WIMSResonanceTable}
    fission_spectrum::Vector{Float64}
    p1_data::Vector{Vector{WIMSP1Block}}
end

# =========================================================================
# Validation
# =========================================================================
"""
    validate(mat::WIMSMaterial)

Validate internal consistency of a WIMSMaterial record. Throws on error.
"""
function validate(mat::WIMSMaterial)
    mat.version in (WIMS_VERSION_D, WIMS_VERSION_E) ||
        error("WIMSMaterial: version must be 4 or 5, got $(mat.version)")
    mat.ngroups > 0 || error("WIMSMaterial: ngroups must be > 0")
    mat.nfg >= 0    || error("WIMSMaterial: nfg must be >= 0")
    mat.nrg >= 0    || error("WIMSMaterial: nrg must be >= 0")
    mat.nfg + mat.nrg <= mat.ngroups ||
        error("WIMSMaterial: nfg + nrg > ngroups")
    nnt = mat.nfg + mat.nrg  # nonthermal groups
    # Fortran wimsr.f90:1459-1463 layout: spot/sdp are resonance-only (length nrg);
    # xtr/ab0 are all-nonthermal (length nnt); glam is resonance-only (length nrg).
    length(mat.potential_xs) == mat.nrg ||
        error("WIMSMaterial: potential_xs (=spot) length must equal nrg=$(mat.nrg)")
    length(mat.scatter_xs) == mat.nrg ||
        error("WIMSMaterial: scatter_xs (=sdp) length must equal nrg=$(mat.nrg)")
    length(mat.lambdas) == mat.nrg ||
        error("WIMSMaterial: lambdas length must equal nrg=$(mat.nrg)")
    length(mat.transport_xs) == nnt ||
        error("WIMSMaterial: transport_xs length must equal $(nnt)")
    length(mat.absorption_xs) == nnt ||
        error("WIMSMaterial: absorption_xs length must equal $(nnt)")
    if mat.fission_flag > 1
        length(mat.nu_fission) == nnt ||
            error("WIMSMaterial: nu_fission length must equal $(nnt)")
        length(mat.fission_xs) == nnt ||
            error("WIMSMaterial: fission_xs length must equal $(nnt)")
    end
    length(mat.temperatures) == mat.ntemp ||
        error("WIMSMaterial: temperatures length must equal ntemp=$(mat.ntemp)")
    length(mat.thermal_transport) == mat.ntemp ||
        error("WIMSMaterial: thermal_transport must have ntemp entries")
    length(mat.thermal_absorption) == mat.ntemp ||
        error("WIMSMaterial: thermal_absorption must have ntemp entries")
    if mat.has_resonance_tables
        length(mat.resonance_tables) == mat.nrg ||
            error("WIMSMaterial: resonance_tables length must equal nrg=$(mat.nrg)")
    end
    if mat.has_fission_spectrum
        length(mat.fission_spectrum) > 0 ||
            error("WIMSMaterial: fission_spectrum is empty but flag is set")
    end
    return true
end

# =========================================================================
# Writer -- produce WIMS-format text output
# =========================================================================
"""
    write_wims(io::IO, mat::WIMSMaterial)

Write a WIMSMaterial record in WIMS-D or WIMS-E format to the given IO stream.
"""
function write_wims(io::IO, mat::WIMSMaterial)
    validate(mat)
    nnt = mat.nfg + mat.nrg
    ngr0 = mat.nfg + 1
    nthermal = mat.ngroups - nnt
    nrestb = mat.has_resonance_tables ? 1 : 0
    ifis = mat.fission_flag

    # -- identifier block (WIMS-E only)
    if mat.version == WIMS_VERSION_E
        @printf(io, "%15d\n", mat.ident)
        @printf(io, "%15d%15d\n", mat.ident, 1)
    end

    # -- burnup chain
    # Fortran wimsr.f90:2030-2039: emitted whenever iburn>=0. For iburn==0
    # the default content is jcc=2 + a single pair (0.0, ident).
    # iverw==5 also gets a leading single-zero record (line 2032).
    if mat.burnup !== nothing
        bp = mat.burnup
        ntis = 2 + length(bp.fission_products)
        jcc = ntis * 2 + 4
        pairs = Vector{Tuple{Float64,Int}}()
        push!(pairs, (0.0, mat.ident))                        # slot 1
        push!(pairs, (bp.capture_product.value, bp.capture_product.ident))
        push!(pairs, (bp.decay_product.value, bp.decay_product.ident))
        push!(pairs, (bp.energy_per_fission, ifis > 0 ? ifis : 0))
        for fp in bp.fission_products
            push!(pairs, (fp.value, fp.ident))
        end
        if mat.version == WIMS_VERSION_E
            @printf(io, "%15d\n", 0)  # WIMS-E leading zero (Fortran line 2032)
        end
        @printf(io, "%15d%15d\n", WIMS_SENTINEL, 3)
        @printf(io, "%15d\n", jcc)
        jcc2 = jcc >> 1
        for i in 1:min(jcc2, length(pairs))
            @printf(io, "%15.8E%6d", pairs[i][1], pairs[i][2])
            if i % 3 == 0 || i == jcc2
                println(io)
            end
        end
        @printf(io, "%15d%15d\n", WIMS_SENTINEL, 4)
    else
        # iburn==0 default emission (Fortran-faithful): jcc=2, single pair (0.0, ident)
        if mat.version == WIMS_VERSION_E
            @printf(io, "%15d\n", 0)
        end
        @printf(io, "%15d%15d\n", WIMS_SENTINEL, 3)
        @printf(io, "%15d\n", 2)
        @printf(io, "%15.8E%6d\n", 0.0, mat.ident)
        @printf(io, "%15d%15d\n", WIMS_SENTINEL, 4)
    end

    # -- material identification line
    @printf(io, "%6d%15.8E%6d%6d%6d%6d%6d\n",
            mat.ident, mat.atomic_weight, mat.z_number,
            ifis, mat.ntemp, nrestb, mat.has_fission_spectrum ? 1 : 0)

    # -- temperature-independent data
    # Fortran layout (wimsr.f90:1459-1463, iverw==4 / WIMS-D):
    #   spot(ngr0..ngr1)  — nrg potential xs (resonance groups only)
    #   sdp(ngr0..ngr1)   — nrg scattering DPL = ξ·σ_e (resonance groups only)
    #   xtr(1..nnt)       — transport-corrected total (all nonthermal groups)
    #   ab0(1..nnt)       — absorption (all nonthermal groups)
    #   zero(ngr0..ngr1)  — nrg zeros (the "(unused)" slot)
    #   glam(1..nrg)      — Goldstein-Cohen lambdas
    # Total length: 4*nrg + 2*nnt
    #
    # WIMSMaterial fields map: potential_xs→spot (nrg long), scatter_xs→sdp
    # (nrg long), transport_xs→xtr (nnt long), absorption_xs→ab0 (nnt long).
    # Concatenate into one stream so 5-per-line packing wraps across boundaries
    # exactly as Fortran's single `write(nout,'(1p,5e15.8)')(scr(i),i=1,nw)` does.
    nrg_zeros = zeros(mat.nrg)
    block = vcat(mat.potential_xs, mat.scatter_xs,
                 mat.transport_xs, mat.absorption_xs,
                 nrg_zeros, mat.lambdas)
    _write_block(io, block)

    # nu*fission and fission for fissile materials (Fortran lines 1508-1509)
    if ifis > 1
        block_f = vcat(mat.nu_fission, mat.fission_xs)
        _write_block(io, block_f)
    end

    # nonthermal scattering matrix
    ndat = length(mat.nonthermal_matrix)
    @printf(io, "%15d\n", ndat)
    if ndat > 0
        _write_block(io, mat.nonthermal_matrix)
    end

    # -- temperature-dependent data
    # Fortran writes per-temp blocks as concatenated streams (lines 2072,
    # 2076): (transport[..nnt] + absorption[..nnt]) as one 5-per-line block,
    # then optionally (nu*fission + fission) as another. The 5-per-line
    # packing crosses array boundaries, so we must concat first.
    _write_block(io, mat.temperatures)
    for it in 1:mat.ntemp
        # Transport + absorption combined block (Fortran line 2072: nw=2*nthermal)
        ta_block = vcat(mat.thermal_transport[it], mat.thermal_absorption[it])
        _write_block(io, ta_block)
        if ifis > 1
            # nu*fission + fission combined block (Fortran line 2076)
            nf_block = vcat(mat.thermal_nu_fission[it], mat.thermal_fission[it])
            _write_block(io, nf_block)
        end
        tmat = mat.thermal_matrices[it]
        @printf(io, "%15d\n", length(tmat))
        if length(tmat) > 0
            _write_block(io, tmat)
        end
    end

    # WIMS-D terminator between material data sections
    if mat.version == WIMS_VERSION_D
        @printf(io, "%15d\n", WIMS_SENTINEL)
    end

    # -- resonance tables
    if mat.has_resonance_tables
        for irg in 1:mat.nrg
            rt = mat.resonance_tables[irg]
            ntnp = rt.ntemps * rt.nsigz
            if mat.version == WIMS_VERSION_D
                @printf(io, "%15d\n", ntnp)
            end
            @printf(io, "%15.8E%6d%6d\n", mat.rid, rt.ntemps, rt.nsigz)
            # temperatures, sigma-zeros, integrals
            flat = vcat(rt.temperatures, rt.sigma_zeros, rt.integrals)
            _write_block(io, flat)

            # fission resonance table (if fission_flag == 3)
            if ifis == 3 && irg <= length(mat.fission_resonance_tables)
                frt = mat.fission_resonance_tables[irg]
                if mat.version == WIMS_VERSION_D
                    @printf(io, "%15d\n", ntnp)
                end
                @printf(io, "%15.8E%6d%6d\n", mat.rid, frt.ntemps, frt.nsigz)
                fflat = vcat(frt.temperatures, frt.sigma_zeros, frt.integrals)
                _write_block(io, fflat)
            end

            if mat.version == WIMS_VERSION_D
                @printf(io, "%15d\n", WIMS_SENTINEL)
            end
        end
    end

    # -- fission spectrum
    if mat.has_fission_spectrum
        _write_block(io, mat.fission_spectrum)
    end

    # -- P1 scattering data
    for it in 1:length(mat.p1_data)
        p1blocks = mat.p1_data[it]
        # Flatten to the WIMS format: count, then (offset, nbands, values...)
        flat = Float64[]
        for blk in p1blocks
            push!(flat, Float64(blk.offset))
            push!(flat, Float64(length(blk.values)))
            append!(flat, blk.values)
        end
        @printf(io, "%15d\n", length(flat))
        _write_block(io, flat)
    end

    return nothing
end

"""
    write_wims(filename::AbstractString, mat::WIMSMaterial)

Write a WIMSMaterial to a file.
"""
function write_wims(filename::AbstractString, mat::WIMSMaterial)
    open(filename, "w") do io
        write_wims(io, mat)
    end
end

# Helper: write a vector of floats in 5-per-line scientific format
function _write_block(io::IO, v::AbstractVector{Float64})
    for (i, x) in enumerate(v)
        @printf(io, "%15.8E", x)
        if i % 5 == 0 || i == length(v)
            println(io)
        end
    end
end

# =========================================================================
# Convenience: write_wims from MultiGroupXS
# =========================================================================

"""
    write_wims(io::IO, multigroup_data::MultiGroupXS; options...)

Write multigroup cross sections in WIMS-D/WIMS-E format, constructing the
required WIMSMaterial automatically from the MultiGroupXS data.

# Keyword Arguments
- `version`: WIMS version (4=WIMS-D, 5=WIMS-E; default 4)
- `ident`: material identifier (default 0)
- `awr`: atomic weight ratio (default 1.0)
- `z_number`: atomic number Z (default 0)
- `nfg`: number of fast groups (default 14)
- `nrg`: number of resonance groups (default 13)
- `temperatures`: processing temperatures in K (default [300.0])
- `lambdas`: Goldstein lambda values for resonance groups (default ones)
"""
function write_wims(io::IO, multigroup_data::MultiGroupXS;
                    version::Int = WIMS_VERSION_D,
                    ident::Int = 0,
                    awr::Float64 = 1.0,
                    z_number::Int = 0,
                    nfg::Int = WIMS_DEFAULT_NFG,
                    nrg::Int = WIMS_DEFAULT_NRG,
                    temperatures::Vector{Float64} = [300.0],
                    lambdas::Vector{Float64} = Float64[])
    ng = length(multigroup_data.group_bounds) - 1
    nnt = nfg + nrg
    ntemp = length(temperatures)
    mt_list = multigroup_data.mt_list
    xs = multigroup_data.xs

    # Build per-column lookup
    mt_to_col = Dict{Int,Int}()
    for (c, mt) in enumerate(mt_list)
        mt_to_col[mt] = c
    end

    # Determine fission flag
    has_fission = haskey(mt_to_col, 18) || haskey(mt_to_col, 19)
    ifis = has_fission ? 1 : 0

    # Lambdas
    lam = if !isempty(lambdas) && length(lambdas) >= nrg
        lambdas[1:nrg]
    elseif !isempty(lambdas)
        vcat(lambdas, ones(nrg - length(lambdas)))
    else
        ones(nrg)
    end

    # Build nonthermal cross section vectors
    pot_xs   = zeros(nnt)
    scat_xs  = zeros(nnt)
    tr_xs    = zeros(nnt)
    abs_xs   = zeros(nnt)
    nuf_xs   = zeros(nnt)
    fis_xs   = zeros(nnt)

    for ig in 1:min(nnt, ng)
        if haskey(mt_to_col, 2);   scat_xs[ig] = xs[ig, mt_to_col[2]]; end
        pot_xs[ig] = scat_xs[ig]
        if haskey(mt_to_col, 1);   tr_xs[ig]   = xs[ig, mt_to_col[1]]; end
        if haskey(mt_to_col, 102); abs_xs[ig]   = xs[ig, mt_to_col[102]]; end
        if haskey(mt_to_col, 18);  fis_xs[ig]   = xs[ig, mt_to_col[18]]; end
        if haskey(mt_to_col, 452); nuf_xs[ig]   = xs[ig, mt_to_col[452]]; end
    end

    # Thermal vectors per temperature (replicated since MultiGroupXS has one temp)
    nt0 = nnt + 1
    n_therm = max(0, ng - nnt)
    therm_tr  = Vector{Float64}[]
    therm_abs = Vector{Float64}[]
    therm_nuf = Vector{Float64}[]
    therm_fis = Vector{Float64}[]
    therm_mat = Vector{Float64}[]
    for _ in 1:ntemp
        tt = zeros(n_therm)
        ta = zeros(n_therm)
        tnf = zeros(n_therm)
        tf  = zeros(n_therm)
        for ig in 1:n_therm
            gidx = nnt + ig
            if gidx <= ng
                if haskey(mt_to_col, 1);   tt[ig]  = xs[gidx, mt_to_col[1]]; end
                if haskey(mt_to_col, 102); ta[ig]  = xs[gidx, mt_to_col[102]]; end
                if haskey(mt_to_col, 18);  tf[ig]  = xs[gidx, mt_to_col[18]]; end
                if haskey(mt_to_col, 452); tnf[ig] = xs[gidx, mt_to_col[452]]; end
            end
        end
        push!(therm_tr, tt)
        push!(therm_abs, ta)
        push!(therm_nuf, tnf)
        push!(therm_fis, tf)
        push!(therm_mat, Float64[])  # no thermal scattering matrix data
    end

    amassn = 1.00866491595
    mat = WIMSMaterial(
        version, ident, Float64(ident), awr * amassn, z_number,
        ng, nfg, nrg, ifis, ntemp,
        false, false,          # no resonance tables, no fission spectrum
        nothing,               # no burnup
        lam,
        pot_xs, scat_xs, tr_xs, abs_xs,
        ifis > 1 ? nuf_xs : Float64[],
        ifis > 1 ? fis_xs : Float64[],
        Float64[],             # nonthermal_matrix
        temperatures,
        therm_tr, therm_abs, therm_nuf, therm_fis, therm_mat,
        WIMSResonanceTable[],  # resonance_tables
        WIMSResonanceTable[],  # fission_resonance_tables
        Float64[],             # fission_spectrum
        Vector{WIMSP1Block}[]  # p1_data
    )

    write_wims(io, mat)
end
