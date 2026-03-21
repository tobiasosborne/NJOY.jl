# CCCCR -- CCCC standard interface files (ISOTXS, BRKOXS, DLAYXS)
# Correspondence to NJOY2016 ccccr.f90: cisotx->write_isotxs, cbrkxs->write_brkoxs, cdlyxs->write_dlayxs
# CCCC-IV binary format: records bracketed by 4-byte length markers.

"""Write a Fortran unformatted record: 4-byte length prefix, data, 4-byte length suffix."""
function _write_record(io::IO, data::Vector{UInt8})
    n = Int32(length(data))
    write(io, htol(n)); write(io, data); write(io, htol(n))
end

"""Serialize items (Int32, Float32, String, Vector) into a byte buffer for a Fortran record."""
function _record_buf(items...)
    buf = IOBuffer()
    for item in items
        if item isa Int32;           write(buf, htol(item))
        elseif item isa Float32;     write(buf, htol(item))
        elseif item isa Float64;     write(buf, htol(item))
        elseif item isa AbstractString; write(buf, codeunits(rpad(item, 8)[1:8]))
        elseif item isa AbstractVector
            for v in item
                if v isa Int32;           write(buf, htol(v))
                elseif v isa Float32;     write(buf, htol(v))
                elseif v isa Float64;     write(buf, htol(v))
                elseif v isa AbstractString; write(buf, codeunits(rpad(v, 8)[1:8]))
                end
            end
        end
    end
    take!(buf)
end

_pad8(s::AbstractString) = rpad(s, 8)[1:8]

# =========================================================================
# ISOTXS writer
# =========================================================================
"""
    write_isotxs(io::IO, multigroup_data::MultiGroupXS; title, ivers, isotope_name, atom_mass, classification)

Write an ISOTXS (CCCC-IV) binary file. Records: file identification, file control,
file data (group structure), isotope control, principal cross sections, scattering sub-block.
"""
function write_isotxs(io::IO, multigroup_data::MultiGroupXS;
                      title::AbstractString="NJOY.jl", ivers::Integer=0,
                      isotope_name::AbstractString="iso1",
                      atom_mass::Real=1.0, classification::Integer=0)
    gb = multigroup_data.group_bounds
    ngroup = length(gb) - 1
    nscmax = 1; nsblok = 1; maxord = 0
    # Record 1: File identification -- hname(a8), huse(2*a8), ivers(i4)
    _write_record(io, _record_buf("isotxs ", "NJOY.jl ", "        ", Int32(ivers)))
    # Record 2: File control
    _write_record(io, _record_buf(Int32[ngroup, 1, 0, ngroup-1, 0, 0, nscmax, nsblok]))
    # Record 3: File data -- hsetid(12*a8), hisonm, vel, emax, emin, loca
    buf3 = IOBuffer()
    pt = rpad(title, 72)[1:72]
    for i in 1:9; write(buf3, codeunits(pt[(i-1)*8+1:i*8])); end
    for _ in 10:12; write(buf3, codeunits("        ")); end
    write(buf3, codeunits(_pad8(isotope_name)))
    for _ in 1:ngroup; write(buf3, htol(Float32(0.0))); end
    for g in 1:ngroup; write(buf3, htol(Float32(gb[ngroup-g+2]))); end
    write(buf3, htol(Float32(gb[1])))
    write(buf3, htol(Int32(1)))
    _write_record(io, take!(buf3))
    # Record 4: Isotope control
    buf4 = IOBuffer()
    for s in [isotope_name, "endf", "mat"]; write(buf4, codeunits(_pad8(s))); end
    write(buf4, htol(Float32(atom_mass)))
    for _ in 1:5; write(buf4, htol(Float32(0.0))); end
    has_fis = any(==(18), multigroup_data.mt_list)
    has_cap = any(==(102), multigroup_data.mt_list)
    ifis = has_fis ? 1 : 0; ialf = has_cap ? 1 : 0
    for f in Int32[classification, 0, ifis, ialf, 0, 0, 0, 0, 1, 1, 0]
        write(buf4, htol(f))
    end
    write(buf4, htol(Int32(2))); write(buf4, htol(Int32(1)))
    for _ in 1:nscmax*ngroup; write(buf4, htol(Int32(1))); end
    for g in 1:ngroup; write(buf4, htol(Int32(g))); end
    _write_record(io, take!(buf4))
    # Record 5: Principal cross sections
    buf5 = IOBuffer()
    mt2c = Dict{Int,Int}(mt => c for (c, mt) in enumerate(multigroup_data.mt_list))
    xs = multigroup_data.xs
    for g in 1:ngroup
        st = haskey(mt2c, 1) ? Float32(xs[g, mt2c[1]]) : Float32(0.0)
        write(buf5, htol(st)); write(buf5, htol(st))
        if ialf == 1
            sc = haskey(mt2c, 102) ? Float32(xs[g, mt2c[102]]) : Float32(0.0)
            write(buf5, htol(sc))
        end
        if ifis == 1
            sf = haskey(mt2c, 18) ? Float32(xs[g, mt2c[18]]) : Float32(0.0)
            write(buf5, htol(sf))
            write(buf5, htol(Float32(g == 1 ? 1.0 : 0.0)))
        end
    end
    _write_record(io, take!(buf5))
    # Record 6: Scattering sub-block (elastic, P0, diagonal)
    buf6 = IOBuffer()
    ecol = get(mt2c, 2, 0)
    for g in 1:ngroup
        write(buf6, htol(ecol > 0 ? Float32(xs[g, ecol]) : Float32(0.0)))
    end
    _write_record(io, take!(buf6))
    return nothing
end

# =========================================================================
# BRKOXS writer
# =========================================================================
"""
    write_brkoxs(io::IO, shielded_data; sigma0_values, temperatures, ...)

Write BRKOXS (CCCC-IV) binary for Bondarenko self-shielding factors.
`shielded_data` is a 4D array (ngroups, nreact, nsigma0, ntemp).
"""
function write_brkoxs(io::IO, shielded_data::AbstractArray{<:Real,4};
                      sigma0_values::AbstractVector{<:Real},
                      temperatures::AbstractVector{<:Real},
                      title::AbstractString="NJOY.jl", ivers::Integer=0,
                      isotope_name::AbstractString="iso1",
                      group_bounds::AbstractVector{<:Real}=Float64[])
    ng, nreact, nsig, ntemp = size(shielded_data)
    _write_record(io, _record_buf("brkoxs ", "NJOY.jl ", "        ", Int32(ivers)))
    _write_record(io, _record_buf(Int32[ng, 1, nsig, ntemp, nreact, 0]))
    # File data
    buf3 = IOBuffer()
    write(buf3, codeunits(_pad8(isotope_name)))
    for s in sigma0_values; write(buf3, htol(Float32(s))); end
    for t in temperatures;  write(buf3, htol(Float32(t))); end
    if !isempty(group_bounds)
        for g in 1:ng; write(buf3, htol(Float32(group_bounds[ng-g+2]))); end
        write(buf3, htol(Float32(group_bounds[1])))
    else
        for _ in 1:ng+1; write(buf3, htol(Float32(0.0))); end
    end
    for v in Int32[1, ng, nsig, ntemp]; write(buf3, htol(v)); end
    _write_record(io, take!(buf3))
    # Self-shielding factors
    buf4 = IOBuffer()
    for r in 1:nreact, g in 1:ng, t in 1:ntemp, s in 1:nsig
        write(buf4, htol(Float32(shielded_data[g, r, s, t])))
    end
    _write_record(io, take!(buf4))
    # Infinite-dilution cross sections (6*ng)
    buf5 = IOBuffer()
    for r in 1:min(6, nreact), g in 1:ng
        write(buf5, htol(Float32(shielded_data[g, r, 1, 1])))
    end
    for _ in nreact+1:6, _ in 1:ng; write(buf5, htol(Float32(0.0))); end
    _write_record(io, take!(buf5))
    return nothing
end

# =========================================================================
# DLAYXS writer
# =========================================================================
"""
    write_dlayxs(io::IO, delayed_data; title, ivers)

Write DLAYXS (CCCC-III) binary for delayed neutron data.
`delayed_data` has fields: ngroup, nfamilies, decay_constants, spectra, group_bounds, yields.
"""
function write_dlayxs(io::IO, delayed_data; title::AbstractString="NJOY.jl", ivers::Integer=0)
    ngroup = delayed_data.ngroup; nfam = delayed_data.nfamilies
    decay = delayed_data.decay_constants; spectra = delayed_data.spectra
    gb = delayed_data.group_bounds; yields = delayed_data.yields
    iso_name = hasproperty(delayed_data, :isotope_name) ? delayed_data.isotope_name : "iso1"
    _write_record(io, _record_buf("dlayxs ", "NJOY.jl ", "        ", Int32(ivers)))
    _write_record(io, _record_buf(Int32[ngroup, 1, nfam, nfam]))
    # File data
    buf3 = IOBuffer()
    write(buf3, codeunits(_pad8(iso_name)))
    for f in 1:nfam; write(buf3, htol(Float32(decay[f]))); end
    for f in 1:nfam, g in 1:ngroup; write(buf3, htol(Float32(spectra[g, f]))); end
    for g in 1:ngroup+1; write(buf3, htol(Float32(gb[g]))); end
    write(buf3, htol(Int32(nfam))); write(buf3, htol(Int32(1)))
    _write_record(io, take!(buf3))
    # Precursor yield data
    buf4 = IOBuffer()
    for f in 1:nfam, g in 1:ngroup; write(buf4, htol(Float32(yields[g, f]))); end
    for f in 1:nfam; write(buf4, htol(Int32(f))); end
    _write_record(io, take!(buf4))
    return nothing
end
