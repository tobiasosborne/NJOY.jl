# DTF-IV and ANISN format types, validation, and I/O
#
# Proposer-B design: clean Julia types for the DTF/ANISN cross section table
# format as implemented in NJOY2016 dtfr.f90.
#
# The DTF format contains:
#   - Table header: isotope name, material, temperature, sigma-zero selection
#   - Neutron tables: ngroup rows x table_length columns with positions for
#     total, absorption, nu*sigma_f, in-group scattering, and band scattering
#   - Edit cross sections (optional, in-table or separate CLAW format)
#   - Photon production tables: ngamma_groups x ngroups
#   - Multiple Legendre orders (P0, P1, ..., P_{nlmax-1})

using Printf

# =========================================================================
# Edit specification
# =========================================================================
"""
    DTFEdit

Specification for one edit cross section column in a DTF table.

# Fields
- `name`: 6-character name (e.g. "nusigf", "absorp")
- `position`: column position in the table
- `mt`: ENDF reaction number (or 300 for total scatter, 1000 for computed)
- `multiplicity`: multiplier when summing this MT into the edit position
"""
struct DTFEdit
    name::String
    position::Int
    mt::Int
    multiplicity::Int
end

function validate(e::DTFEdit)
    e.position > 0 || error("DTFEdit: position must be > 0, got $(e.position)")
    return true
end

# =========================================================================
# Table layout parameters
# =========================================================================
"""
    DTFLayout

Describes the column layout of a DTF neutron table.

# Fields
- `nlmax`: number of Legendre orders (neutron tables)
- `ngroups`: number of neutron groups
- `iptotl`: column position of total cross section
- `ipingp`: column position of in-group scattering
- `itabl`: total table length (columns per group)
- `ntherm`: number of thermal groups (0 = none)
- `mti`: MT for thermal incoherent scattering (0 = none)
- `mtc`: MT for thermal coherent scattering (0 = none)
- `nlc`: number of coherent Legendre orders
- `edits`: vector of edit specifications
- `nptabl`: number of photon tables (0 = none)
- `ngamma`: number of photon groups (0 = none)
"""
struct DTFLayout
    nlmax::Int
    ngroups::Int
    iptotl::Int
    ipingp::Int
    itabl::Int
    ntherm::Int
    mti::Int
    mtc::Int
    nlc::Int
    edits::Vector{DTFEdit}
    nptabl::Int
    ngamma::Int
end

function validate(lay::DTFLayout)
    lay.nlmax >= 0  || error("DTFLayout: nlmax must be >= 0")
    lay.ngroups > 0 || error("DTFLayout: ngroups must be > 0")
    if lay.nlmax > 0
        lay.iptotl > 0 || error("DTFLayout: iptotl must be > 0")
        lay.ipingp > lay.iptotl ||
            error("DTFLayout: ipingp must be > iptotl")
        lay.itabl >= lay.ipingp ||
            error("DTFLayout: itabl must be >= ipingp")
    end
    lay.ntherm >= 0 && lay.ntherm <= lay.ngroups ||
        error("DTFLayout: ntherm out of range")
    return true
end

"""
    claw_layout(; nlmax=5, ngroups=30) -> DTFLayout

Create a standard CLAW (iedit=1) layout with 48 edits and the predefined
edit map from NJOY dtfr.f90.
"""
function claw_layout(; nlmax::Int=5, ngroups::Int=30)
    iptotl = 43 + 3
    ipingp = iptotl + 1
    itabl = iptotl + ngroups

    # Standard CLAW edit map (abbreviated to first 8 most common)
    claw_edits = [
        DTFEdit("  els", 1, 2, 1),
        DTFEdit("  ins", 2, 4, 1),
        DTFEdit("  n2n", 3, 16, 1),
        DTFEdit("  ngm", 5, 102, 1),
        DTFEdit("  nal", 6, 107, 1),
        DTFEdit("   np", 7, 103, 1),
        DTFEdit("  nnf", 8, 19, 1),
        DTFEdit(" ftot", 11, 18, 1),
    ]

    DTFLayout(nlmax, ngroups, iptotl, ipingp, itabl, 0, 0, 0, 0,
              claw_edits, 0, 0)
end

# =========================================================================
# Neutron cross section table
# =========================================================================
"""
    DTFNeutronTable

A single Legendre-order neutron cross section table in DTF format.

# Fields
- `legendre_order`: Legendre order l (0-based; P0=0, P1=1, ...)
- `data`: matrix of size (ngroups, itabl) -- each row is one group,
          columns hold total, absorption, nu*sigma_f, in-group, band scatter
"""
struct DTFNeutronTable
    legendre_order::Int
    data::Matrix{Float64}
end

function validate(t::DTFNeutronTable, lay::DTFLayout)
    size(t.data, 1) == lay.ngroups ||
        error("DTFNeutronTable: row count must equal ngroups=$(lay.ngroups)")
    size(t.data, 2) == lay.itabl ||
        error("DTFNeutronTable: column count must equal itabl=$(lay.itabl)")
    return true
end

# =========================================================================
# Photon production table
# =========================================================================
"""
    DTFPhotonTable

Photon production cross section table.

# Fields
- `legendre_order`: Legendre order (0-based)
- `data`: matrix of size (ngroups, ngamma) -- row = source neutron group,
          column = sink photon group
"""
struct DTFPhotonTable
    legendre_order::Int
    data::Matrix{Float64}
end

# =========================================================================
# Complete DTF material
# =========================================================================
"""
    DTFMaterial

Complete DTF-format cross section data for one isotope at one temperature
and sigma-zero selection.

# Fields
- `name`: 6-character isotope name
- `mat`: ENDF material number
- `sigma_zero_index`: which sigma-zero was selected (1-based)
- `temperature`: temperature in Kelvin
- `layout`: table layout parameters
- `neutron_tables`: vector of DTFNeutronTable (one per Legendre order)
- `photon_tables`: vector of DTFPhotonTable
- `capture_ssf`: capture self-shielding factors per group (empty if none)
- `fission_ssf`: fission self-shielding factors per group (empty if none)
- `chi`: fission spectrum per group (empty if non-fissile)
- `nu_sigf`: nu*sigma_f per group (empty if non-fissile)
"""
struct DTFMaterial
    name::String
    mat::Int
    sigma_zero_index::Int
    temperature::Float64
    layout::DTFLayout
    neutron_tables::Vector{DTFNeutronTable}
    photon_tables::Vector{DTFPhotonTable}
    capture_ssf::Vector{Float64}
    fission_ssf::Vector{Float64}
    chi::Vector{Float64}
    nu_sigf::Vector{Float64}
end

# =========================================================================
# Validation
# =========================================================================
"""
    validate(mat::DTFMaterial)

Validate internal consistency. Throws on error.
"""
function validate(mat::DTFMaterial)
    validate(mat.layout)
    lay = mat.layout
    mat.mat > 0 || error("DTFMaterial: mat must be > 0")
    mat.sigma_zero_index > 0 ||
        error("DTFMaterial: sigma_zero_index must be > 0")
    mat.temperature >= 0.0 ||
        error("DTFMaterial: temperature must be >= 0")
    length(mat.neutron_tables) <= lay.nlmax ||
        error("DTFMaterial: too many neutron tables")
    for nt in mat.neutron_tables
        validate(nt, lay)
    end
    for pt in mat.photon_tables
        size(pt.data, 1) == lay.ngroups ||
            error("DTFPhotonTable: row count must equal ngroups")
        lay.ngamma > 0 && size(pt.data, 2) == lay.ngamma ||
            error("DTFPhotonTable: column count must equal ngamma")
    end
    if !isempty(mat.chi)
        length(mat.chi) == lay.ngroups ||
            error("DTFMaterial: chi length must equal ngroups")
    end
    if !isempty(mat.nu_sigf)
        length(mat.nu_sigf) == lay.ngroups ||
            error("DTFMaterial: nu_sigf length must equal ngroups")
    end
    return true
end

# =========================================================================
# Writer -- DTF format output
# =========================================================================
"""
    write_dtf(io::IO, mat::DTFMaterial; format=:internal)

Write a DTFMaterial in DTF format.
  - `format=:internal` (iedit=0): tables with edits in-line
  - `format=:claw`     (iedit=1): separate CLAW/TD6 format edits
"""
function write_dtf(io::IO, mat::DTFMaterial; format::Symbol=:internal)
    validate(mat)
    lay = mat.layout

    if format == :internal
        _write_dtf_internal(io, mat)
    elseif format == :claw
        _write_dtf_claw(io, mat)
    else
        error("write_dtf: unknown format $format, use :internal or :claw")
    end
end

function _write_dtf_internal(io::IO, mat::DTFMaterial)
    lay = mat.layout
    for nt in mat.neutron_tables
        il = nt.legendre_order + 1
        @printf(io, "\n il= %d table%3d gp%3d pos, mat=%5d iz=%2d temp=%12.5E\n",
                il, lay.ngroups, lay.itabl, mat.mat, mat.sigma_zero_index,
                mat.temperature)
        nw = lay.itabl * lay.ngroups
        # Write column-major flat data: groups vary first (row), then columns
        k = 0
        for ig in 1:lay.ngroups
            for jpos in 1:lay.itabl
                k += 1
                @printf(io, "%12.4E", nt.data[ig, jpos])
                if k % 6 == 0
                    println(io)
                end
            end
        end
        if k % 6 != 0
            println(io)
        end
    end

    # Photon tables
    for pt in mat.photon_tables
        ip = pt.legendre_order + 1
        @printf(io, "\n ip= %d table%3d gp%3d pos, mat=%5d iz=%2d temp=%12.5E\n",
                ip, lay.ngroups, lay.ngamma, mat.mat, mat.sigma_zero_index,
                mat.temperature)
        k = 0
        for ig in 1:lay.ngroups
            for igp in 1:lay.ngamma
                k += 1
                @printf(io, "%12.4E", pt.data[ig, igp])
                if k % 6 == 0
                    println(io)
                end
            end
        end
        if k % 6 != 0
            println(io)
        end
    end
end

function _write_dtf_claw(io::IO, mat::DTFMaterial)
    lay = mat.layout
    hname = rpad(mat.name, 6)[1:min(6, length(rpad(mat.name, 6)))]
    ned = length(lay.edits)

    # Only process if we have at least P0
    if isempty(mat.neutron_tables)
        return
    end

    nt0 = mat.neutron_tables[1]  # P0 table

    # Edit cross sections header
    @printf(io, "%s     edit xsec (%3dx%3d) proc by njoy.jl\n",
            hname, lay.ngroups, ned)

    # Write each edit column
    for edit in lay.edits
        iseq = 0
        k = 0
        dat = zeros(6)
        for ig in 1:lay.ngroups
            k += 1
            jpos = edit.position
            dat[k] = (jpos <= size(nt0.data, 2)) ? nt0.data[ig, jpos] : 0.0
            if k >= 6 || ig >= lay.ngroups
                iseq += 1
                for ii in (k+1):6
                    dat[ii] = 0.0
                end
                @printf(io, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E%2s%5s%1d\n",
                        dat[1], dat[2], dat[3], dat[4], dat[5], dat[6],
                        hname[1:min(2,end)], edit.name[1:min(5,end)], iseq)
                k = 0
                fill!(dat, 0.0)
            end
        end
    end

    # Write reduced neutron tables for higher Legendre orders
    for idx in 2:length(mat.neutron_tables)
        nt = mat.neutron_tables[idx]
        l = nt.legendre_order
        med = lay.iptotl - 3
        ltabn = lay.itabl - lay.iptotl + 3
        @printf(io, "%s l=%d n-n table (%3dx%3d)\n",
                hname, l, ltabn, lay.ngroups)
        iseq = 0
        k = 0
        dat = zeros(6)
        for ig in 1:lay.ngroups
            for jpos in 1:lay.itabl
                if jpos > med
                    k += 1
                    dat[k] = nt.data[ig, jpos]
                    if k >= 6 || (ig == lay.ngroups && jpos == lay.itabl)
                        iseq += 1
                        for ii in (k+1):6
                            dat[ii] = 0.0
                        end
                        @printf(io, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E%2s%2d%4d\n",
                                dat[1], dat[2], dat[3], dat[4], dat[5], dat[6],
                                hname[1:min(2,end)], l, iseq)
                        k = 0
                        fill!(dat, 0.0)
                    end
                end
            end
        end
    end
end

"""
    write_dtf(filename::AbstractString, mat::DTFMaterial; kwargs...)

Write a DTFMaterial to a file.
"""
function write_dtf(filename::AbstractString, mat::DTFMaterial; kwargs...)
    open(filename, "w") do io
        write_dtf(io, mat; kwargs...)
    end
end

# =========================================================================
# Accessor helpers
# =========================================================================
"""
    total_xs(mat::DTFMaterial, group::Int) -> Float64

Return the total cross section for the given group from the P0 table.
"""
function total_xs(mat::DTFMaterial, group::Int)
    isempty(mat.neutron_tables) && return 0.0
    mat.neutron_tables[1].data[group, mat.layout.iptotl]
end

"""
    absorption_xs(mat::DTFMaterial, group::Int) -> Float64

Return the absorption cross section for the given group from the P0 table.
Absorption = total - total_scattering.
"""
function absorption_xs(mat::DTFMaterial, group::Int)
    isempty(mat.neutron_tables) && return 0.0
    mat.neutron_tables[1].data[group, mat.layout.iptotl - 2]
end

"""
    ingroup_xs(mat::DTFMaterial, group::Int) -> Float64

Return the in-group (self-scatter) cross section.
"""
function ingroup_xs(mat::DTFMaterial, group::Int)
    isempty(mat.neutron_tables) && return 0.0
    mat.neutron_tables[1].data[group, mat.layout.ipingp]
end

# =========================================================================
# Convenience: write_dtf from MultiGroupXS
# =========================================================================

"""
    write_dtf(io::IO, multigroup_data::MultiGroupXS; format=:anisn, options...)

Write multigroup cross sections in DTF-IV or ANISN format for discrete
ordinates transport codes.  Constructs a DTFMaterial automatically from
the MultiGroupXS data.

# Keyword Arguments
- `format`: output style (:internal or :claw; default :internal)
- `mat_id`: ENDF MAT number (default 1)
- `name`: 6-character isotope name (default "mat   ")
- `temperature`: temperature in K (default 300.0)
- `sigma_zero_index`: sigma-zero selection index (default 1)
- `nlmax`: number of Legendre orders (default 1, P0 only)
- `iptotl`: position of total cross section (default 4)
- `ipingp`: position of in-group scattering (default 5)
"""
function write_dtf(io::IO, multigroup_data::MultiGroupXS;
                   format::Symbol = :internal,
                   mat_id::Int = 1,
                   name::AbstractString = "mat   ",
                   temperature::Float64 = 300.0,
                   sigma_zero_index::Int = 1,
                   nlmax::Int = 1,
                   iptotl::Int = 4,
                   ipingp::Int = 5)
    ng = length(multigroup_data.group_bounds) - 1
    itabl = max(ipingp, ng)

    # Build MT lookup
    mt_list = multigroup_data.mt_list
    xs = multigroup_data.xs
    mt_to_col = Dict{Int,Int}()
    for (c, mt) in enumerate(mt_list)
        mt_to_col[mt] = c
    end

    # Build P0 table
    tdata = zeros(Float64, ng, itabl)
    for ig in 1:ng
        # Total at iptotl
        sig_tot = haskey(mt_to_col, 1) ? xs[ig, mt_to_col[1]] : 0.0
        tdata[ig, iptotl] = sig_tot

        # Elastic / in-group at ipingp
        sig_el = haskey(mt_to_col, 2) ? xs[ig, mt_to_col[2]] : 0.0
        tdata[ig, ipingp] = sig_el

        # Absorption at iptotl - 2 (= total - scattering)
        if iptotl > 2
            tdata[ig, iptotl - 2] = sig_tot - sig_el
        end

        # Nu*sigf at iptotl - 1
        if iptotl > 1 && haskey(mt_to_col, 18)
            tdata[ig, iptotl - 1] = xs[ig, mt_to_col[18]]
        end
    end

    p0 = DTFNeutronTable(0, tdata)

    lay = DTFLayout(nlmax, ng, iptotl, ipingp, itabl, 0, 0, 0, 0,
                    DTFEdit[], 0, 0)

    mat = DTFMaterial(
        rpad(name, 6)[1:6],
        mat_id,
        sigma_zero_index,
        temperature,
        lay,
        [p0],
        DTFPhotonTable[],
        Float64[], Float64[],     # SSF
        Float64[], Float64[]      # chi, nu_sigf
    )

    write_dtf(io, mat; format=format)
end
