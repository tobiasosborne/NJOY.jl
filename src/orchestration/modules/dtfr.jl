# dtfr module runner -- DTF-format cross section tables from GENDF
#
# Matches Fortran dtfr.f90: read GENDF tape (from gaminr/groupr),
# build DTF tables, write DTF output and viewr plot tape.

"""
    dtfr_module(tapes::TapeManager, params::DtfrParams)

Run DTFR: produce DTF-format cross section tables from GENDF data.
Reads GENDF tape (from gaminr), writes DTF output tape and viewr plot tape.
"""
function dtfr_module(tapes::TapeManager, params::DtfrParams)
    @info "dtfr: nin=$(params.nin) nout=$(params.nout) npend=$(params.npend) nplot=$(params.nplot)"
    @info "dtfr: $(length(params.materials)) materials, nlmax=$(params.nlmax), ng=$(params.ng)"

    gendf_path = resolve(tapes, params.nin)
    nout_path = resolve(tapes, params.nout)
    nplot_path = params.nplot > 0 ? resolve(tapes, params.nplot) : ""

    # Read GENDF tape
    gendf = _read_gendf_tape(gendf_path)

    # Write DTF output tape (tape34)
    if params.nout > 0
        open(nout_path, "w") do io
            _write_dtf_output(io, gendf, params)
        end
        register!(tapes, params.nout, nout_path)
        @info "dtfr: wrote DTF $nout_path ($(countlines(nout_path)) lines)"
    end

    # Write viewr plot tape (tape36)
    if params.nplot > 0
        pendf_path = params.npend > 0 ? resolve(tapes, params.npend) : ""
        open(nplot_path, "w") do io
            _write_dtfr_plot_tape(io, gendf, params, pendf_path)
        end
        register!(tapes, params.nplot, nplot_path)
        @info "dtfr: wrote plot tape $nplot_path ($(countlines(nplot_path)) lines)"
    end

    nothing
end

# =========================================================================
# GENDF tape reader
# =========================================================================

struct GendfSection
    mt::Int
    data::Vector{Tuple{Float64,Float64}}  # (flux, sigma) per group
end

struct GendfMaterial
    mat::Int; za::Float64; awr::Float64
    ngg::Int; egg::Vector{Float64}
    sections::Vector{GendfSection}
end

"""Read a GENDF tape produced by gaminr. Returns vector of GendfMaterial."""
function _read_gendf_tape(path::String)
    materials = GendfMaterial[]
    lines = readlines(path)
    isempty(lines) && return materials

    # Parse ENDF-format records
    i = 2  # skip TPID
    while i <= length(lines)
        line = rpad(get(lines, i, ""), 80)
        mat_n = _parse_int(line[67:70])
        mf = _parse_int(line[71:72])
        mt = _parse_int(line[73:75])

        if mat_n > 0 && mf == 1 && mt == 451
            # MF1/MT451 header — parse material
            za = parse_endf_float(strip(line[1:11]))
            awr = parse_endf_float(strip(line[12:22]))
            i += 1
            # Next line: CONT with ngg
            line2 = rpad(get(lines, i, ""), 80)
            ngg = _parse_int(line2[23:33])
            nw = _parse_int(line2[45:55])
            i += 1

            # Read nw data values (group boundaries etc.)
            vals = Float64[]
            while length(vals) < nw && i <= length(lines)
                dline = rpad(lines[i], 80)
                for col in 0:5
                    length(vals) >= nw && break
                    s = strip(dline[col*11+1:col*11+11])
                    isempty(s) && continue
                    push!(vals, parse_endf_float(s))
                end
                i += 1
            end

            # Extract group boundaries from vals
            # Format: title_real, emax, egg(1..ngg+1), 0.0
            egg = Float64[]
            if length(vals) >= 3
                # Skip the first value (title packed as real) and emax
                start_idx = 3  # vals[1]=title, vals[2]=emax, vals[3..3+ngg]=egg
                for j in start_idx:min(start_idx+ngg, length(vals))
                    push!(egg, vals[j])
                end
            end

            # Now read MF3/MF23 sections until MEND
            sections = GendfSection[]
            while i <= length(lines)
                fline = rpad(get(lines, i, ""), 80)
                fm = _parse_int(fline[67:70])
                fmf = _parse_int(fline[71:72])
                fmt = _parse_int(fline[73:75])

                if fm == 0 && fmf == 0 && fmt == 0
                    i += 1; break  # MEND
                end

                if fm == mat_n && (fmf == 3 || fmf == 23 || fmf == 6 || fmf == 26) && fmt > 0
                    # Section HEAD
                    i += 1
                    # Read group data records
                    data = Tuple{Float64,Float64}[]
                    while i <= length(lines)
                        gline = rpad(get(lines, i, ""), 80)
                        gmat = _parse_int(gline[67:70])
                        gmf = _parse_int(gline[71:72])
                        gmt = _parse_int(gline[73:75])

                        if gmt == 0  # SEND
                            i += 1; break
                        end

                        # CONT line: ng2, ig2lo, nw, ig
                        ng2 = _parse_int(gline[23:33])
                        ig = _parse_int(gline[45:55])
                        i += 1

                        # Data line(s)
                        dvals = Float64[]
                        target_nw = ng2  # number of data words
                        while length(dvals) < max(target_nw, 2) && i <= length(lines)
                            dl = rpad(get(lines, i, ""), 80)
                            dlmat = _parse_int(dl[67:70])
                            dlmf = _parse_int(dl[71:72])
                            dlmt = _parse_int(dl[73:75])
                            (dlmat != gmat || dlmf != gmf || dlmt != fmt) && break
                            for col in 0:5
                                s = strip(dl[col*11+1:col*11+11])
                                isempty(s) && continue
                                push!(dvals, parse_endf_float(s))
                            end
                            i += 1
                        end

                        # Extract flux and sigma from data
                        flux_val = length(dvals) >= 1 ? dvals[1] : 0.0
                        sig_val  = length(dvals) >= 2 ? dvals[2] : 0.0
                        push!(data, (flux_val, sig_val))
                    end

                    push!(sections, GendfSection(fmt, data))
                else
                    i += 1
                end
            end

            push!(materials, GendfMaterial(mat_n, za, awr, ngg, egg, sections))
        else
            i += 1
        end
    end

    materials
end

# =========================================================================
# DTF output writer — matches Fortran dtfr dtout
# =========================================================================

function _write_dtf_output(io::IO, gendf::Vector{GendfMaterial}, params::DtfrParams)
    for dmat in params.materials
        # Find this material in GENDF
        gmat = nothing
        for gm in gendf
            gm.mat == dmat.mat && (gmat = gm; break)
        end
        gmat === nothing && continue

        ng = params.ng
        itabl = params.itabl

        # Build the cross section table for this material
        # Layout: itabl columns × ng rows per Legendre order
        for il in 1:params.nlmax
            @printf(io, "\n il= %1d table%3d gp%3d pos, mat=%5d iz=%2d temp=%12.5E\n",
                    il, ng, itabl, dmat.mat, dmat.jsigz, dmat.dtemp)

            # Build data matrix
            data = zeros(Float64, ng, itabl)

            if il == 1  # P0 only for now
                # Fill from GENDF sections
                for sec in gmat.sections
                    col = _mt_to_dtf_column(sec.mt, params)
                    col == 0 && continue
                    for g in 1:min(ng, length(sec.data))
                        _, sigma = sec.data[g]
                        data[g, col] += sigma
                    end
                end

                # Total XS at iptotl
                # Absorption at iptotl-2
                # In-group at ipingp
                # Scattering matrix fills columns ipingp+1 .. itabl
            end

            # Write data in column-major order, 6 per line
            # Fortran: write(nout,'(1p,6e12.4)') (sig(j),j=1,nw)
            nw = ng * itabl
            idx = 1
            for g in 1:ng
                for j in 1:itabl
                    if idx == 1
                        # Start new line
                    end
                    @printf(io, "%12.4E", data[g, j])
                    if idx % 6 == 0
                        println(io)
                    end
                    idx += 1
                end
            end
            idx % 6 != 1 && println(io)
        end
    end
end

function _mt_to_dtf_column(mt::Int, params::DtfrParams)
    # Map GENDF MT to DTF table column position
    mt == 501 && return params.iptotl      # total
    mt == 502 && return params.ipingp      # coherent → in-group
    mt == 504 && return params.ipingp      # incoherent → in-group
    # Edit specs
    for spec in params.edit_specs
        spec.mt == mt && return spec.jpos
    end
    return 0
end

# =========================================================================
# Plot tape writer — matches Fortran dtfr dtfplt
# =========================================================================

function _write_dtfr_plot_tape(io::IO, gendf::Vector{GendfMaterial},
                               params::DtfrParams, pendf_path::String)
    # Write NJOY plot tape format matching Fortran dtfplt exactly
    # This format is read by viewr

    # Card 1: page setup (landscape, Swiss font, size 0.300)
    @printf(io, "1 2 0.300 /\n")

    plot_idx = 0
    for dmat in params.materials
        gmat = nothing
        for gm in gendf
            gm.mat == dmat.mat && (gmat = gm; break)
        end
        gmat === nothing && continue

        # For each MT that has data, write a plot
        for sec in gmat.sections
            plot_idx += 1
            ngg = gmat.ngg
            egg = gmat.egg

            # Build step-function data for multigroup XS
            xs_data = Float64[]
            for g in 1:min(ngg, length(sec.data))
                _, sigma = sec.data[g]
                push!(xs_data, sigma)
            end
            isempty(xs_data) && continue

            # Write plot frame
            # Card 2: new plot
            @printf(io, " 1/\n")

            # Card 3: title
            mt_name = _photon_mt_label(sec.mt)
            @printf(io, " '<%8s         %s '/\n", rpad(dmat.name, 8), mt_name)

            # Card 3a: subtitle (blank)
            @printf(io, " /\n")

            # Card 4: axis types (log-log, outside tics, legend)
            # Determine axis ranges from data
            xmin_log = floor(log10(egg[1]))
            xmax_log = ceil(log10(egg[end]))
            ymin_log = 2.0  # 100 barns default min
            ymax_log = 6.0  # 1e6 barns default max

            # Find actual Y range
            valid_xs = filter(x -> x > 0, xs_data)
            if !isempty(valid_xs)
                ymin_log = floor(log10(minimum(valid_xs)))
                ymax_log = ceil(log10(maximum(valid_xs)))
            end

            @printf(io, " 4 0 3 1%10.2E%10.2E/\n", 10.0^xmin_log * 1.2, 10.0^ymax_log * 1.1)
            @printf(io, "%10.2E%10.2E%10.2E/\n", egg[1], egg[end], 1.0)
            @printf(io, " '<e>nergy (e<v>)'/\n")
            @printf(io, "%10.2E%10.2E%10.2E/\n", 10.0^ymin_log, 10.0^ymax_log, 1.0)
            @printf(io, " '<c>ross <s>ection (barns)'/\n")

            # Card 8: dummy
            @printf(io, "/\n")

            # Card 9: curve style (solid, no symbols)
            @printf(io, " 0 0 0 0 1/\n")

            # Card 10: legend
            @printf(io, " '<m.g.>'/\n")

            # Card 12: format code (0 = free 2D)
            @printf(io, " 0\n")

            # Card 13: data points — step function
            for g in 1:min(ngg, length(xs_data))
                g > length(egg) - 1 && break
                e_lo = egg[g]
                e_hi = egg[g+1]
                sigma = xs_data[g]
                xs_plot = sigma > 0 ? sigma : 100.0  # floor for log plot
                @printf(io, "%14.4E%14.4E/\n", e_lo, xs_plot)
                @printf(io, "%14.4E%14.4E/\n", e_hi, xs_plot)
            end
            @printf(io, " /\n")

            # PENDF overlay if available
            if !isempty(pendf_path) && isfile(pendf_path)
                _write_pendf_overlay(io, pendf_path, dmat.mat, sec.mt, plot_idx)
            end
        end
    end

    # Terminator
    @printf(io, " 99/\n")
end

function _write_pendf_overlay(io::IO, pendf_path::String, mat::Int, mt::Int, plot_idx::Int)
    # Read PENDF MF23 section for this MT and write as overlay curve
    mf23 = _read_pendf_mf23(pendf_path, mat)
    haskey(mf23, mt) || return

    energies, xs_vals = mf23[mt]
    isempty(energies) && return

    # Card 9: overlay curve style (line, no symbols, dashed)
    @printf(io, "%3d 0 0 0 1/\n", plot_idx + 1)
    @printf(io, " /\n")

    # Card 12: format code
    @printf(io, " 0\n")

    # Data points
    for i in 1:length(energies)
        xs_vals[i] > 0 || continue
        @printf(io, "%14.4E%14.4E/\n", energies[i], xs_vals[i])
    end
    @printf(io, " /\n")
end

function _photon_mt_label(mt::Int)
    mt == 501 && return "total"
    mt == 502 && return "coherent"
    mt == 504 && return "incoherent"
    mt == 516 && return "pair production"
    mt == 522 && return "photoelectric"
    mt == 602 && return "photoelectric"
    return "mt$mt"
end
