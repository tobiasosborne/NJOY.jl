# groupr module runner -- Group-averaged cross sections from PENDF
#
# Matches Fortran groupr.f90: read pointwise data from PENDF,
# compute flux-weighted group averages, write GENDF tape.

"""
    groupr_module(tapes::TapeManager, params::GrouprParams)

Run GROUPR: compute multigroup cross sections from PENDF data.
Write output GENDF tape. Matches Fortran groupr.f90 interface.

Input tapes: nendf (ENDF), npend (PENDF from reconr)
Output tapes: nout (GENDF multigroup output)
"""
function groupr_module(tapes::TapeManager, params::GrouprParams)
    @info "groupr: MAT=$(params.mat) ign=$(params.ign) → tape $(params.nout)"

    endf_path = resolve(tapes, params.nendf)
    pendf_path = resolve(tapes, params.npend)
    nout_path = resolve(tapes, params.nout)

    # Get group structure
    egn = _groupr_group_structure(params)
    ngn = length(egn) - 1
    @info "groupr: $ngn groups, $(length(params.mt_list)) MT requests"

    # Read ZA and AWR
    za, awr = _read_za_awr(endf_path, params.mat)

    # Temperature and sigz
    temp = isempty(params.temperatures) ? 0.0 : params.temperatures[1]
    sigz = isempty(params.sigz) ? 1.0e10 : params.sigz[1]

    # Weight function
    wfn = _groupr_weight_function(params.iwt)

    # Read PENDF data for fission XS (needed for nubar weighting)
    tape = read_pendf(pendf_path)
    mf3 = extract_mf3_all(tape, params.mat)

    # Expand Fortran auto-reaction sentinel. Ref: groupr.f90:622,628-632 +
    # nextr at groupr.f90:1104-1123. A deck card like `3 /` (mfd=3, mtd
    # absent) leaves mtdp=-1000; Fortran then walks MF=3 on the PENDF and
    # yields every MT passing: mt<=200 OR 203<=mt<=207 OR mt>300 (thermal
    # 201-202 and derived 208-300 excluded; those must be named explicitly).
    mt_list = _groupr_expand_auto(params.mt_list, mf3)

    # Process each MT request
    mt_results = Vector{NamedTuple{(:mfd,:mt,:name,:records),
                                   Tuple{Int,Int,String,Vector{NTuple{3,Float64}}}}}()

    for (mfd, mtd, name) in mt_list
        if mfd == 3
            if mtd == 452 || mtd == 455 || mtd == 456
                # Nubar: read from ENDF MF1, weight by fission XS from PENDF
                nubar_data = _read_nubar(endf_path, params.mat, mtd)
                if nubar_data !== nothing
                    # Get fission XS from PENDF (MT=18 or MT=19)
                    fis_mt = haskey(mf3, 18) ? 18 : (haskey(mf3, 19) ? 19 : 0)
                    records = _groupr_nubar_records(nubar_data, mf3, fis_mt,
                                                    egn, wfn)
                    push!(mt_results, (mfd=mfd, mt=mtd, name=name, records=records))
                end
            else
                # Standard XS from PENDF
                if haskey(mf3, mtd)
                    records = _groupr_xs_records(mf3[mtd], egn, wfn)
                    push!(mt_results, (mfd=mfd, mt=mtd, name=name, records=records))
                end
            end
        end
    end

    # Write GENDF output tape
    open(nout_path, "w") do io
        _write_groupr_tape(io, params.mat, za, awr, egn, mt_results,
                           params.title, temp, sigz)
    end
    register!(tapes, params.nout, nout_path)

    lines = countlines(nout_path)
    @info "groupr: wrote $nout_path ($lines lines, $ngn groups)"
    nothing
end

# =========================================================================
# Auto-reaction expansion (Fortran nextr sentinel)
# =========================================================================

"""
    _nextr_filter(mt) -> Bool

Inclusion predicate for Fortran `nextr`'s MF<=3 branch
(groupr.f90:1115-1117). Returns true iff `mt` would be yielded by the
auto-walk of MF=3 on the PENDF tape. Excludes thermal (201-202) and
engineering/derived (208-300) reactions — those must be named
explicitly in the deck.
"""
_nextr_filter(mt::Integer) =
    mt > 0 && (mt <= 200 || (203 <= mt <= 207) || mt > 300)

"""
    _groupr_expand_auto(mt_list, mf3) -> Vector{Tuple{Int,Int,String}}

Expand any `(mfd, -1000, name)` sentinel entries into the full set of
auto-yielded MTs from the PENDF MF=3 dict. Scoped to mfd=3, matching
Fortran nextr's MF<=3 branch; other mfd sentinels are currently passed
through unchanged (the mfd=6/8/10/16-36 auto-paths depend on `conver`
lists not yet ported).

Ref: groupr.f90:622 (detect), 628-635 (iauto=1), 1087-1123 (nextr body).
"""
function _groupr_expand_auto(mt_list::Vector{Tuple{Int,Int,String}},
                             mf3::AbstractDict)
    expanded = Tuple{Int,Int,String}[]
    auto_mts = sort!(filter(_nextr_filter, collect(keys(mf3))))
    for (mfd, mtd, name) in mt_list
        if mtd == -1000 && mfd == 3
            for mt in auto_mts
                push!(expanded, (mfd, mt, ""))
            end
        else
            push!(expanded, (mfd, mtd, name))
        end
    end
    return expanded
end

# =========================================================================
# Group structure
# =========================================================================

function _groupr_group_structure(params::GrouprParams)
    ign = params.ign
    if ign == 3
        return collect(Float64, LANL_30)
    elseif ign == 5
        return collect(Float64, RRD_50)
    elseif ign == 9
        return collect(Float64, WIMS_69)
    else
        # Try the generic get_group_structure
        try
            return collect(Float64, get_group_structure(ign))
        catch
            @warn "groupr: unsupported ign=$ign, using LANL-30"
            return collect(Float64, LANL_30)
        end
    end
end

function _groupr_weight_function(iwt::Int)
    iwt == 2 && return constant_weight
    return inv_e_weight  # iwt=3 and default
end

# =========================================================================
# Nubar reading from ENDF MF1
# =========================================================================

"""Read nubar (MF1/MT452, 455, or 456) from ENDF file as (energies, values)."""
function _read_nubar(endf_path::String, mat::Int, mt::Int)
    open(endf_path) do io
        seekstart(io)
        find_section(io, 1, mt; target_mat=mat) || return nothing
        head = read_cont(io)
        lnu = Int(head.L2)  # LNU: 1=polynomial, 2=tabulated
        if lnu == 1
            # Polynomial representation: ν(E) = c₁ + c₂E + c₃E² + ...
            lst = read_list(io)
            nc = Int(lst.N1)
            coeffs = lst.data[1:nc]
            # Generate dense tabulation for proper integration
            # Use logarithmic spacing to cover full energy range
            energies = Float64[]
            for dec in -5:8
                for sub in [1.0, 2.0, 5.0]
                    e = sub * 10.0^dec
                    push!(energies, e)
                end
            end
            push!(energies, 2e7)
            sort!(unique!(energies))
            values = [sum(coeffs[k] * e^(k-1) for k in 1:nc) for e in energies]
            return (energies=energies, values=values)
        elseif lnu == 2
            # Tabulated. For MT=455 (delayed), a LIST of NNF decay
            # constants precedes the TAB1 (ENDF-6 §1.5; Fortran groupr.f90
            # getyld label 110 at line 6472).
            mt == 455 && read_list(io)
            tab = read_tab1(io)
            return (energies=collect(Float64, tab.x), values=collect(Float64, tab.y))
        end
        nothing
    end
end

# =========================================================================
# Group averaging using existing group_integrate infrastructure
# =========================================================================

"""Compute GENDF records for nubar: (flux, nubar_g, sigf_avg) per group.
Nubar weighted by fission XS: ν_g = ∫ν·σf·W dE / ∫σf·W dE.
Uses trapezoidal with 1% stepping matching Fortran groupr panel+getwtf."""
function _groupr_nubar_records(nubar_data, mf3::Dict, fis_mt::Int,
                                egn::Vector{Float64}, wfn)
    ngn = length(egn) - 1
    nu_e = nubar_data.energies; nu_v = nubar_data.values

    if fis_mt > 0 && haskey(mf3, fis_mt)
        fis_e, fis_xs = mf3[fis_mt]
    else
        fis_e = nu_e; fis_xs = ones(Float64, length(nu_e))
    end

    # Merge grids for panel boundaries
    merged_e = sort!(unique!(vcat(nu_e, fis_e, egn)))
    nu_merged = _interp_linlin(merged_e, nu_e, nu_v)
    fis_merged = _interp_linlin(merged_e, fis_e, fis_xs)

    s101 = 1.01  # Fortran stepping factor
    records = NTuple{3,Float64}[]
    ip = 1

    for g in 1:ngn
        elo, ehi = egn[g], egn[g+1]
        flux_g = 0.0; sigf_int_g = 0.0; nu_sigf_int_g = 0.0

        while ip < length(merged_e) - 1 && merged_e[ip+1] <= elo; ip += 1; end

        jp = ip
        while jp < length(merged_e) && merged_e[jp] < ehi
            jn = min(jp + 1, length(merged_e))
            e1, e2 = merged_e[jp], merged_e[jn]
            sf1, sf2 = fis_merged[jp], fis_merged[jn]
            nu1, nu2 = nu_merged[jp], nu_merged[jn]
            ea = max(e1, elo); eb = min(e2, ehi)
            ea >= eb && (jp += 1; continue)

            # 1% sub-stepping matching Fortran getwtf s101=1.01
            e_lo = ea
            while e_lo < eb
                e_hi = min(s101 * e_lo, eb)
                de = e2 > e1 ? e2 - e1 : 1.0
                frac_lo = (e_lo - e1) / de; frac_hi = (e_hi - e1) / de
                sf_lo = sf1 + frac_lo * (sf2 - sf1)
                sf_hi = sf1 + frac_hi * (sf2 - sf1)
                nu_lo = nu1 + frac_lo * (nu2 - nu1)
                nu_hi = nu1 + frac_hi * (nu2 - nu1)

                bq = (e_hi - e_lo) / 2.0
                fl_lo = 1.0 / e_lo; fl_hi = 1.0 / e_hi
                flux_g += (fl_lo + fl_hi) * bq
                sigf_int_g += (sf_lo * fl_lo + sf_hi * fl_hi) * bq
                nu_sigf_int_g += (nu_lo * sf_lo * fl_lo + nu_hi * sf_hi * fl_hi) * bq
                e_lo = e_hi
            end
            jp += 1
        end

        nubar_g = sigf_int_g > 0 ? nu_sigf_int_g / sigf_int_g : 0.0
        sigf_avg = flux_g > 0 ? sigf_int_g / flux_g : 0.0
        push!(records, (round_sigfig(flux_g, 7),
                       round_sigfig(nubar_g, 7),
                       round_sigfig(sigf_avg, 7)))
    end
    records
end

"""Compute GENDF records for standard XS: (flux, sigma_g, sigma_flux) per group."""
function _groupr_xs_records(xs_data::Tuple, egn::Vector{Float64}, wfn)
    energies, xs = xs_data
    ngn = length(egn) - 1

    w_vals = [wfn(energies[i]) for i in 1:length(energies)]
    sig_w = [xs[i] * w_vals[i] for i in 1:length(energies)]

    flux = group_integrate(energies, w_vals, egn)
    sig_int = group_integrate(energies, sig_w, egn)

    records = NTuple{3,Float64}[]
    for g in 1:ngn
        sigma_g = flux[g] > 0 ? sig_int[g] / flux[g] : 0.0
        push!(records, (flux[g], sigma_g, sig_int[g]))
    end
    records
end

"""Linear interpolation of (x_data, y_data) at points x_out."""
function _interp_linlin(x_out::Vector{Float64}, x_data::Vector{Float64},
                        y_data::Vector{Float64})
    n = length(x_out)
    y_out = Vector{Float64}(undef, n)
    j = 1
    for i in 1:n
        x = x_out[i]
        if x <= x_data[1]
            y_out[i] = y_data[1]
        elseif x >= x_data[end]
            y_out[i] = y_data[end]
        else
            while j < length(x_data) - 1 && x_data[j+1] < x
                j += 1
            end
            frac = (x - x_data[j]) / (x_data[j+1] - x_data[j])
            y_out[i] = y_data[j] + frac * (y_data[j+1] - y_data[j])
        end
    end
    y_out
end

# =========================================================================
# GENDF output tape writer
# =========================================================================

"""Write GROUPR output tape in GENDF format matching Fortran groupr."""
function _write_groupr_tape(io::IO, mat::Int, za::Float64, awr::Float64,
                            egn::Vector{Float64}, mt_results, title::String,
                            temp::Float64, sigz::Float64)
    ngn = length(egn) - 1

    # TPID record
    @printf(io, "%-66s%4d%2d%3d%5d\n", title, 0, 0, 0, 0)

    seq = 1
    # MF1/MT451 HEAD: ZA, AWR, 0, ntemp=1, NFC=-1, nmod=1
    _write_cont_line(io, za, awr, 0, 1, -1, 1, mat, 1, 451, seq); seq += 1
    # CONT: 0, 0, NGN, 0, NW, 0
    nw_dir = ngn + 4  # 2 (elow,ehigh) + ngn+1 (bounds) + 1 (trailing zero)
    _write_cont_line(io, 0.0, 0.0, ngn, 0, nw_dir, 0, mat, 1, 451, seq); seq += 1

    # Write elow=0, ehigh=max, then group boundaries, then trailing zero
    all_vals = Float64[0.0, 1e10]
    append!(all_vals, egn)
    push!(all_vals, 0.0)

    idx = 1
    while idx <= length(all_vals)
        buf = ""
        for col in 1:6
            idx > length(all_vals) && break
            buf *= format_endf_float(all_vals[idx]); idx += 1
        end
        _write_data_line(io, buf, mat, 1, 451, seq); seq += 1
    end

    # Fortran groupr (see referenceTape24): MF1/MT451 goes straight to FEND,
    # no SEND. MF3 sections use SEND=99999 but skip FEND; MEND follows directly.
    _write_fend_line(io, mat)

    for res in mt_results
        mt = res.mt
        records = res.records
        seq = 1

        _write_cont_line(io, za, 0.0, 1, 1, 0, ngn, mat, 3, mt, seq); seq += 1

        for (g, rec) in enumerate(records)
            _write_cont_line(io, 0.0, 0.0, 3, 1, 3, g, mat, 3, mt, seq); seq += 1
            buf = format_endf_float(rec[1]) *
                  format_endf_float(rec[2]) *
                  format_endf_float(rec[3])
            _write_data_line(io, buf, mat, 3, mt, seq); seq += 1
        end

        _write_send_line(io, mat, 3)
    end

    # MEND + TEND
    _write_fend_line(io, 0)
    @printf(io, "%66s%4d%2d%3d%5d\n", "", -1, 0, 0, 0)
end
