# CCCCR -- CCCC (Committee on Computer Code Coordination) standard interface files
#
# Produces ISOTXS (isotope cross sections), BRKOXS (Bondarenko self-shielding
# factors), and DLAYXS (delayed neutron data) binary output files.
#
# Correspondence to NJOY2016 ccccr.f90:
#   cisotx  -> write_isotxs   (ISOTXS file writer)
#   cbrkxs  -> write_brkoxs   (BRKOXS file writer)
#   cdlyxs  -> write_dlayxs   (DLAYXS file writer)
#
# CCCC-IV binary format: fixed-length records, Fortran unformatted I/O style.
# Each record is bracketed by 4-byte record-length markers.

# =========================================================================
# Binary record helpers (Fortran unformatted I/O emulation)
# =========================================================================

"""
    _write_record(io, data::Vector{UInt8})

Write a single Fortran unformatted record: 4-byte length prefix, data, 4-byte
length suffix, all in little-endian byte order.
"""
function _write_record(io::IO, data::Vector{UInt8})
    n = Int32(length(data))
    write(io, htol(n))
    write(io, data)
    write(io, htol(n))
end

"""
    _record_buf(items...) -> Vector{UInt8}

Serialize items into a byte buffer for a Fortran record.  Accepts Int32,
Float32, Float64, and String (padded/truncated to 8 bytes, Hollerith-style).
"""
function _record_buf(items...)
    buf = IOBuffer()
    for item in items
        if item isa Int32
            write(buf, htol(item))
        elseif item isa Float32
            write(buf, htol(item))
        elseif item isa Float64
            write(buf, htol(item))
        elseif item isa AbstractString
            s = rpad(item, 8)[1:8]
            write(buf, codeunits(s))
        elseif item isa AbstractVector
            for v in item
                if v isa Int32
                    write(buf, htol(v))
                elseif v isa Float32
                    write(buf, htol(v))
                elseif v isa Float64
                    write(buf, htol(v))
                elseif v isa AbstractString
                    s = rpad(v, 8)[1:8]
                    write(buf, codeunits(s))
                end
            end
        end
    end
    take!(buf)
end

"""
    _pad8(s::AbstractString) -> String

Pad or truncate a string to exactly 8 characters (Hollerith word).
"""
_pad8(s::AbstractString) = rpad(s, 8)[1:8]

# =========================================================================
# ISOTXS writer
# =========================================================================

"""
    write_isotxs(io::IO, multigroup_data::MultiGroupXS;
                 title="NJOY.jl", ivers=0, isotope_name="iso1",
                 atom_mass=1.0, classification=0)

Write an ISOTXS (CCCC-IV) binary file from multigroup cross section data.

The ISOTXS file contains: file identification, file control, file data
(group structure, set ID, isotope names), and per-isotope records
(isotope control, principal cross sections).

# Arguments
- `io`: output IO stream (must be opened in write mode)
- `multigroup_data`: `MultiGroupXS` with group_bounds, mt_list, xs, flux
- `title`: set identification string (up to 72 characters)
- `ivers`: file version number
- `isotope_name`: 6-character isotope label
- `atom_mass`: gram atomic weight
- `classification`: isotope classification (kbr)
"""
function write_isotxs(io::IO, multigroup_data::MultiGroupXS;
                      title::AbstractString="NJOY.jl",
                      ivers::Integer=0,
                      isotope_name::AbstractString="iso1",
                      atom_mass::Real=1.0,
                      classification::Integer=0)
    gb = multigroup_data.group_bounds
    ngroup = length(gb) - 1
    niso = 1
    maxord = 0       # P0 only
    maxup = 0
    maxdn = ngroup - 1
    ichist = 0       # no chi histogram
    nscmax = 1       # one scattering type (elastic)
    nsblok = 1       # single sub-block

    # Record 1: File identification
    # hname(a8), huse(2*a8), ivers(i4)
    rec1 = _record_buf("isotxs ", "NJOY.jl ", "        ", Int32(ivers))
    _write_record(io, rec1)

    # Record 2: File control
    # ngroup, niso, maxup, maxdn, maxord, ichist, nscmax, nsblok
    rec2 = _record_buf(Int32[ngroup, niso, maxup, maxdn, maxord,
                             ichist, nscmax, nsblok])
    _write_record(io, rec2)

    # Record 3: File data
    # hsetid(12*a8), hisonm(niso*a8), vel(ngroup), emax(ngroup), emin, loca(niso)
    buf3 = IOBuffer()
    # Set identification (12 Hollerith words = 72 chars from title)
    padded_title = rpad(title, 72)[1:72]
    for i in 1:9
        s = padded_title[(i-1)*8+1 : i*8]
        write(buf3, codeunits(s))
    end
    # Pad to 12 words
    for _ in 10:12
        write(buf3, codeunits("        "))
    end
    # Isotope names (niso Hollerith words)
    write(buf3, codeunits(_pad8(isotope_name)))
    # Group velocities (ngroup float32, set to zero = not provided)
    for _ in 1:ngroup
        write(buf3, htol(Float32(0.0)))
    end
    # Emax: upper boundary of each group (descending energy order)
    for g in 1:ngroup
        write(buf3, htol(Float32(gb[ngroup - g + 2])))
    end
    # Emin: minimum energy boundary
    write(buf3, htol(Float32(gb[1])))
    # loca: location pointer for each isotope (first isotope at offset 1)
    write(buf3, htol(Int32(1)))
    _write_record(io, take!(buf3))

    # --- Per-isotope records ---

    # Record 4: Isotope control and group-independent data
    buf4 = IOBuffer()
    # habsid, hident, hmat (3 Hollerith words)
    write(buf4, codeunits(_pad8(isotope_name)))
    write(buf4, codeunits(_pad8("endf")))
    write(buf4, codeunits(_pad8("mat")))
    # amass, efiss, ecapt, temp, sigpot, adens (6 reals)
    write(buf4, htol(Float32(atom_mass)))
    write(buf4, htol(Float32(0.0)))  # efiss
    write(buf4, htol(Float32(0.0)))  # ecapt
    write(buf4, htol(Float32(0.0)))  # temp
    write(buf4, htol(Float32(0.0)))  # sigpot
    write(buf4, htol(Float32(0.0)))  # adens
    # Integer flags: kbr, ichi, ifis, ialf, inp, in2n, ind, int_, ltot, ltrn, istrpd
    n_mt = length(multigroup_data.mt_list)
    has_fission = any(mt -> mt == 18, multigroup_data.mt_list)
    has_capture = any(mt -> mt == 102, multigroup_data.mt_list)
    ichi = 0
    ifis = has_fission ? 1 : 0
    ialf = has_capture ? 1 : 0
    ltot_flag = 1
    ltrn_flag = 1
    int_flags = Int32[classification, ichi, ifis, ialf, 0, 0, 0, 0,
                      ltot_flag, ltrn_flag, 0]
    for f in int_flags
        write(buf4, htol(f))
    end
    # idsct, lord arrays (2*nscmax) -- scattering type IDs and orders
    write(buf4, htol(Int32(2)))   # elastic scattering idsct
    write(buf4, htol(Int32(1)))   # lord = 1 (P0 only)
    # jband(nscmax*ngroup): bandwidth for scattering (set to 1 for diagonal)
    for _ in 1:nscmax*ngroup
        write(buf4, htol(Int32(1)))
    end
    # ijj(nscmax*ngroup): lowest sink group in band
    for g in 1:ngroup
        write(buf4, htol(Int32(g)))
    end
    _write_record(io, take!(buf4))

    # Record 5: Principal cross sections
    # Layout: (transport, total, then reaction-specific) * ngroup
    # Each column is ngroup values. Columns present depend on flags.
    buf5 = IOBuffer()
    # Find relevant MT columns
    mt_to_col = Dict{Int,Int}()
    for (c, mt) in enumerate(multigroup_data.mt_list)
        mt_to_col[mt] = c
    end
    xs = multigroup_data.xs
    for g in 1:ngroup
        # Transport cross section (ltrn=1): use total if available, else zero
        sig_tot = haskey(mt_to_col, 1) ? Float32(xs[g, mt_to_col[1]]) : Float32(0.0)
        write(buf5, htol(sig_tot))  # transport
        # Total cross section (ltot=1)
        write(buf5, htol(sig_tot))  # total
        # Capture (ialf=1)
        if ialf == 1
            sig_c = haskey(mt_to_col, 102) ? Float32(xs[g, mt_to_col[102]]) : Float32(0.0)
            write(buf5, htol(sig_c))
        end
        # Fission (ifis=1): nusigf and chi
        if ifis == 1
            sig_f = haskey(mt_to_col, 18) ? Float32(xs[g, mt_to_col[18]]) : Float32(0.0)
            write(buf5, htol(sig_f))    # nu*sigma_f
            write(buf5, htol(Float32(g == 1 ? 1.0 : 0.0)))  # chi (all in group 1)
        end
    end
    _write_record(io, take!(buf5))

    # Record 6: Scattering sub-block (elastic, P0 only)
    # One value per group (diagonal scattering matrix)
    buf6 = IOBuffer()
    elastic_col = get(mt_to_col, 2, 0)
    for g in 1:ngroup
        sig_el = elastic_col > 0 ? Float32(xs[g, elastic_col]) : Float32(0.0)
        write(buf6, htol(sig_el))
    end
    _write_record(io, take!(buf6))

    return nothing
end

# =========================================================================
# BRKOXS writer
# =========================================================================

"""
    write_brkoxs(io::IO, shielded_data;
                 sigma0_values, temperatures,
                 title="NJOY.jl", ivers=0,
                 isotope_name="iso1", group_bounds=Float64[])

Write a BRKOXS (CCCC-IV) binary file for Bondarenko self-shielding factors.

# Arguments
- `io`: output IO stream
- `shielded_data`: 4D array (ngroups, nreact, nsigma0, ntemp) of f-factors
- `sigma0_values`: background cross section values [barns]
- `temperatures`: temperature values [K]
- `title`: set identification string
- `group_bounds`: energy group boundaries (ascending)
- `isotope_name`: isotope label
"""
function write_brkoxs(io::IO, shielded_data::AbstractArray{<:Real,4};
                      sigma0_values::AbstractVector{<:Real},
                      temperatures::AbstractVector{<:Real},
                      title::AbstractString="NJOY.jl",
                      ivers::Integer=0,
                      isotope_name::AbstractString="iso1",
                      group_bounds::AbstractVector{<:Real}=Float64[])
    ng, nreact, nsig, ntemp = size(shielded_data)

    # Record 1: File identification
    rec1 = _record_buf("brkoxs ", "NJOY.jl ", "        ", Int32(ivers))
    _write_record(io, rec1)

    # Record 2: File control
    # ngroup, niso, nsigpt, ntempt, nreact, iblk
    iblk = 0  # all reactions in one block
    rec2 = _record_buf(Int32[ng, 1, nsig, ntemp, nreact, iblk])
    _write_record(io, rec2)

    # Record 3: File data
    # hisonm(niso*a8), sigpo(nsig), temp(ntemp), emax(ng), emin, jbl/jbh/ntap/ntat
    buf3 = IOBuffer()
    write(buf3, codeunits(_pad8(isotope_name)))
    for s in sigma0_values
        write(buf3, htol(Float32(s)))
    end
    for t in temperatures
        write(buf3, htol(Float32(t)))
    end
    # Group structure: emax(ng) in descending order, then emin
    if !isempty(group_bounds)
        for g in 1:ng
            write(buf3, htol(Float32(group_bounds[ng - g + 2])))
        end
        write(buf3, htol(Float32(group_bounds[1])))
    else
        for g in 1:ng+1
            write(buf3, htol(Float32(0.0)))
        end
    end
    # Integer arrays: jbl, jbh, ntap, ntat (each niso=1 values)
    write(buf3, htol(Int32(1)))      # jbl: first group
    write(buf3, htol(Int32(ng)))     # jbh: last group
    write(buf3, htol(Int32(nsig)))   # ntap: number of sigma0 values
    write(buf3, htol(Int32(ntemp)))  # ntat: number of temperatures
    _write_record(io, take!(buf3))

    # Record 4: Self-shielding factors
    # Dimensions: ntap*ntat*(jbh-jbl+1)*nreact
    buf4 = IOBuffer()
    for r in 1:nreact
        for g in 1:ng
            for t in 1:ntemp
                for s in 1:nsig
                    write(buf4, htol(Float32(shielded_data[g, r, s, t])))
                end
            end
        end
    end
    _write_record(io, take!(buf4))

    # Record 5: Infinite-dilution cross sections (6*ngroup values)
    # Use first sigma0 / first temperature slice as infinite dilution
    buf5 = IOBuffer()
    for r in 1:min(6, nreact)
        for g in 1:ng
            write(buf5, htol(Float32(shielded_data[g, r, 1, 1])))
        end
    end
    # Pad remaining reaction slots to 6
    for _ in nreact+1:6
        for _ in 1:ng
            write(buf5, htol(Float32(0.0)))
        end
    end
    _write_record(io, take!(buf5))

    return nothing
end

# =========================================================================
# DLAYXS writer
# =========================================================================

"""
    write_dlayxs(io::IO, delayed_data;
                 title="NJOY.jl", ivers=0)

Write a DLAYXS (CCCC-III) binary file for delayed neutron data.

# Arguments
- `io`: output IO stream
- `delayed_data`: NamedTuple or Dict with fields:
    - `ngroup`: number of energy groups
    - `nfamilies`: number of delayed neutron precursor families
    - `decay_constants`: Vector of decay constants [1/s] (length nfamilies)
    - `spectra`: Matrix (ngroup, nfamilies) of emission spectra
    - `group_bounds`: energy group boundaries
    - `yields`: Matrix (ngroup, nfamilies) of precursor yields
    - `isotope_name`: isotope label (optional, default "iso1")
"""
function write_dlayxs(io::IO, delayed_data;
                      title::AbstractString="NJOY.jl",
                      ivers::Integer=0)
    ngroup = delayed_data.ngroup
    nfam = delayed_data.nfamilies
    nisod = 1  # single isotope
    decay = delayed_data.decay_constants
    spectra = delayed_data.spectra
    gb = delayed_data.group_bounds
    yields = delayed_data.yields
    iso_name = hasproperty(delayed_data, :isotope_name) ?
               delayed_data.isotope_name : "iso1"

    # Record 1: File identification
    rec1 = _record_buf("dlayxs ", "NJOY.jl ", "        ", Int32(ivers))
    _write_record(io, rec1)

    # Record 2: File control
    # ngroup, nisod, nfam, nfam (nfam appears twice per CCCC-III spec)
    rec2 = _record_buf(Int32[ngroup, nisod, nfam, nfam])
    _write_record(io, rec2)

    # Record 3: File data
    # hisonm(nisod*a8), decay(nfam), spectra(ngroup*nfam),
    # group_bounds(ngroup+1), nkfam(nisod), spec_type(nisod)
    buf3 = IOBuffer()
    write(buf3, codeunits(_pad8(iso_name)))
    # Decay constants
    for f in 1:nfam
        write(buf3, htol(Float32(decay[f])))
    end
    # Emission spectra (ngroup values for each family)
    for f in 1:nfam
        for g in 1:ngroup
            write(buf3, htol(Float32(spectra[g, f])))
        end
    end
    # Group boundaries
    for g in 1:ngroup+1
        write(buf3, htol(Float32(gb[g])))
    end
    # nkfam and spec_type per isotope
    write(buf3, htol(Int32(nfam)))  # nkfam
    write(buf3, htol(Int32(1)))     # spec_type
    _write_record(io, take!(buf3))

    # Record 4: Precursor yield data per isotope
    # yields(ngroup, nkfam), family_indices(nkfam)
    buf4 = IOBuffer()
    for f in 1:nfam
        for g in 1:ngroup
            write(buf4, htol(Float32(yields[g, f])))
        end
    end
    for f in 1:nfam
        write(buf4, htol(Int32(f)))
    end
    _write_record(io, take!(buf4))

    return nothing
end
