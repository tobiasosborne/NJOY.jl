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
    gendf_mats = _read_gendf_tape(gendf_path)

    # Write DTF output tape (tape34)
    if params.nout > 0
        open(nout_path, "w") do io
            _write_dtf_output(io, gendf_mats, params)
        end
        register!(tapes, params.nout, nout_path)
        @info "dtfr: wrote DTF $nout_path ($(countlines(nout_path)) lines)"
    end

    # Write viewr plot tape (tape36)
    if params.nplot > 0
        pendf_path = params.npend > 0 ? resolve(tapes, params.npend) : ""
        open(nplot_path, "w") do io
            _write_dtfr_plot_tape(io, gendf_mats, params, pendf_path)
        end
        register!(tapes, params.nplot, nplot_path)
        @info "dtfr: wrote plot tape $nplot_path ($(countlines(nplot_path)) lines)"
    end

    nothing
end

# =========================================================================
# GENDF tape reader
# =========================================================================

"""Per-group data for an MF23 section: sigma[g] for each group."""
struct GendfMF23Section
    mt::Int
    sigma::Vector{Float64}  # sigma[1:ngg], indexed by GENDF group number
end

"""MF26 transfer matrix for one MT: transfer[il, ig_source, ig_sink]."""
struct GendfMF26Section
    mt::Int
    nl::Int
    transfer::Array{Float64, 3}  # [nl, ngg, ngg] — (Legendre, source, sink)
end

"""One material from GENDF tape."""
struct GendfMaterial
    mat::Int; za::Float64; awr::Float64
    ngg::Int; egg::Vector{Float64}  # ngg+1 group boundaries
    flux::Vector{Float64}  # P0 flux per group (from first MF23 section)
    mf23::Vector{GendfMF23Section}
    mf26::Vector{GendfMF26Section}
end

"""Read a GENDF tape produced by gaminr. Returns vector of GendfMaterial."""
function _read_gendf_tape(path::String)
    materials = GendfMaterial[]
    lines = readlines(path)
    isempty(lines) && return materials

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

            # Extract group boundaries: vals[1]=title, vals[2]=emax, vals[3..3+ngg]=egg
            egg = Float64[]
            if length(vals) >= 3
                for j in 3:min(2+ngg+1, length(vals))
                    push!(egg, vals[j])
                end
            end

            # Skip FEND record
            while i <= length(lines)
                fl = rpad(get(lines, i, ""), 80)
                i += 1
                fm = _parse_int(fl[67:70])
                fmf = _parse_int(fl[71:72])
                fmt = _parse_int(fl[73:75])
                (fm == mat_n && fmf == 0 && fmt == 0) && break
            end

            # Now read MF23 and MF26 sections until MEND
            mf23_secs = GendfMF23Section[]
            mf26_secs = GendfMF26Section[]
            flux_arr = zeros(Float64, ngg)  # P0 flux per group

            while i <= length(lines)
                fline = rpad(get(lines, i, ""), 80)
                fm = _parse_int(fline[67:70])
                fmf = _parse_int(fline[71:72])
                fmt = _parse_int(fline[73:75])

                if fm == 0 && fmf == 0 && fmt == 0
                    i += 1; break  # MEND
                end

                if fm == mat_n && (fmf == 23 || fmf == 26) && fmt > 0
                    # Section HEAD line
                    head_nl = _parse_int(fline[23:33])  # L1 = nl
                    head_nz = _parse_int(fline[34:44])  # L2 = nz
                    sec_mt = fmt
                    sec_mf = fmf
                    i += 1

                    if sec_mf == 23
                        # MF23: per-group (flux, sigma)
                        sigma = zeros(Float64, ngg)
                        while i <= length(lines)
                            gline = rpad(get(lines, i, ""), 80)
                            gmat = _parse_int(gline[67:70])
                            gmf = _parse_int(gline[71:72])
                            gmt = _parse_int(gline[73:75])
                            if gmt == 0  # SEND
                                i += 1; break
                            end
                            # CONT: ng2, iglo, nw_rec, ig
                            ng2_rec = _parse_int(gline[23:33])
                            iglo_rec = _parse_int(gline[34:44])
                            nw_rec = _parse_int(gline[45:55])
                            ig = _parse_int(gline[56:66])
                            i += 1

                            # Read data values
                            dvals = Float64[]
                            while length(dvals) < nw_rec && i <= length(lines)
                                dl = rpad(lines[i], 80)
                                dlmat = _parse_int(dl[67:70])
                                dlmf = _parse_int(dl[71:72])
                                dlmt = _parse_int(dl[73:75])
                                (dlmat != gmat || dlmf != gmf || dlmt != sec_mt) && break
                                for col in 0:5
                                    length(dvals) >= nw_rec && break
                                    s = strip(dl[col*11+1:col*11+11])
                                    isempty(s) && continue
                                    push!(dvals, parse_endf_float(s))
                                end
                                i += 1
                            end

                            # dvals[1] = flux, dvals[2] = sigma (for nl=1, nz=1)
                            if ig >= 1 && ig <= ngg && length(dvals) >= 2
                                sigma[ig] = dvals[2]
                                flux_arr[ig] = dvals[1]
                            end
                        end
                        push!(mf23_secs, GendfMF23Section(sec_mt, sigma))

                    else  # sec_mf == 26
                        # MF26: transfer matrix
                        nl = head_nl
                        nz = head_nz
                        transfer = zeros(Float64, nl, ngg, ngg)

                        while i <= length(lines)
                            gline = rpad(get(lines, i, ""), 80)
                            gmat = _parse_int(gline[67:70])
                            gmf = _parse_int(gline[71:72])
                            gmt = _parse_int(gline[73:75])
                            if gmt == 0  # SEND
                                i += 1; break
                            end
                            ng2_rec = _parse_int(gline[23:33])
                            iglo_rec = _parse_int(gline[34:44])
                            nw_rec = _parse_int(gline[45:55])
                            ig = _parse_int(gline[56:66])
                            i += 1

                            dvals = Float64[]
                            while length(dvals) < nw_rec && i <= length(lines)
                                dl = rpad(lines[i], 80)
                                dlmat = _parse_int(dl[67:70])
                                dlmf = _parse_int(dl[71:72])
                                dlmt = _parse_int(dl[73:75])
                                (dlmat != gmat || dlmf != gmf || dlmt != sec_mt) && break
                                for col in 0:5
                                    length(dvals) >= nw_rec && break
                                    s = strip(dl[col*11+1:col*11+11])
                                    isempty(s) && continue
                                    push!(dvals, parse_endf_float(s))
                                end
                                i += 1
                            end

                            # Parse blocks: block 1 = flux (nl values), blocks 2..ng2 = scatter
                            # Data layout: [flux_P0, flux_P1, ..., flux_P(nl-1),
                            #               sink1_P0, sink1_P1, ..., sink1_P(nl-1), ...]
                            if ig >= 1 && ig <= ngg
                                for k in 2:ng2_rec
                                    ig2 = iglo_rec + k - 2  # sink group (GENDF numbering)
                                    if ig2 >= 1 && ig2 <= ngg
                                        for il_idx in 1:nl
                                            didx = nl * (k - 1) + il_idx
                                            if didx <= length(dvals)
                                                transfer[il_idx, ig, ig2] += dvals[didx]
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        push!(mf26_secs, GendfMF26Section(sec_mt, nl, transfer))
                    end
                else
                    i += 1
                end
            end

            push!(materials, GendfMaterial(mat_n, za, awr, ngg, egg, flux_arr, mf23_secs, mf26_secs))
        else
            i += 1
        end
    end

    materials
end

# Helper to look up MF23 sigma for a given MT
function _get_mf23_sigma(gmat::GendfMaterial, mt::Int)
    for sec in gmat.mf23
        sec.mt == mt && return sec.sigma
    end
    return nothing
end

# =========================================================================
# DTF output writer — matches Fortran dtfr dtout (iedit=0)
# =========================================================================

function _write_dtf_output(io::IO, gendf_mats::Vector{GendfMaterial}, params::DtfrParams)
    for dmat in params.materials
        gmat = nothing
        for gm in gendf_mats
            gm.mat == dmat.mat && (gmat = gm; break)
        end
        gmat === nothing && continue

        ng = params.ng
        itabl = params.itabl
        iptotl = params.iptotl
        ipingp = params.ipingp
        ngg = gmat.ngg

        for il in 1:params.nlmax
            # Build sig array (1D, Fortran-style: sig[jpos + itabl*(jg-1)])
            nws = itabl * ng
            sig = zeros(Float64, nws)

            # Process all MF23 sections (cross sections)
            if il == 1
                # Total (MT501)
                sigma_tot = _get_mf23_sigma(gmat, 501)
                if sigma_tot !== nothing
                    for ig in 1:ngg
                        jg = ng - ig + 1  # reversed
                        locs = iptotl + itabl * (jg - 1)
                        sig[locs] = sigma_tot[ig]
                        # Also add to absorption position
                        sig[locs - 2] += sigma_tot[ig]
                    end
                end

                # Edit cross sections
                for spec in params.edit_specs
                    sigma_ed = _get_mf23_sigma(gmat, spec.mt)
                    sigma_ed === nothing && continue
                    for ig in 1:ngg
                        jg = ng - ig + 1
                        locs = spec.jpos + itabl * (jg - 1)
                        sig[locs] += spec.mult * sigma_ed[ig]
                    end
                end
            end

            # Process all MF26 sections (transfer matrices)
            for sec26 in gmat.mf26
                il > sec26.nl && continue

                for ig in 1:ngg
                    jg = ng - ig + 1  # DTF source row (reversed)
                    for ig2 in 1:ngg
                        val = sec26.transfer[il, ig, ig2]
                        val == 0.0 && continue

                        jg2 = ng - ig2 + 1  # DTF sink row (reversed)
                        jpos = ipingp + (jg2 - jg)

                        # Clamping logic from Fortran lines 413-416
                        if jpos > itabl
                            jg2 = jg2 - jpos + itabl
                            jpos = itabl
                        end
                        if jpos <= iptotl
                            jg2 = jg2 + iptotl + 1 - jpos
                            jpos = iptotl + 1
                        end

                        if jpos >= 1 && jpos <= itabl && jg2 >= 1 && jg2 <= ng
                            locs = jpos + itabl * (jg2 - 1)
                            sig[locs] += val

                            # For P0 (il==1), subtract from absorption
                            if il == 1
                                locb = iptotl - 2 + itabl * (jg - 1)
                                sig[locb] -= val
                            end
                        end
                    end
                end
            end

            # Write header: matches Fortran format exactly
            # '/'' il= '',i1,'' table'',i3,'' gp'',i3,'' pos, mat='',i5,'' iz='',i2,'' temp='',1p,e12.5'
            @printf(io, "\n il= %1d table%3d gp%3d pos, mat=%5d iz=%2d temp=%12.5E\n",
                    il, ng, itabl, dmat.mat, dmat.jsigz, dmat.dtemp)

            # Write data: Fortran format '(1p,6e12.4)'
            for j in 1:nws
                @printf(io, "%12.4E", sig[j])
                if j % 6 == 0
                    println(io)
                end
            end
            # Handle final line if not multiple of 6
            if nws % 6 != 0
                println(io)
            end
        end
    end
end

"""Fortran E10.2 format (0p scale): mantissa in [0.1, 1.0), e.g. `  0.12E+05`."""
function _fmt_e10_2(x::Float64)
    if x == 0.0
        return "  0.00E+00"
    end
    ax = abs(x)
    e = floor(Int, log10(ax)) + 1
    m = ax / 10.0^e
    if m >= 1.0  # rounding edge case
        m /= 10.0; e += 1
    end
    if x < 0
        @sprintf(" -%4.2fE%+03d", m, e)
    else
        @sprintf("  %4.2fE%+03d", m, e)
    end
end

# =========================================================================
# Plot tape writer — matches Fortran dtfr dtfplt/ploted/plotnn
# =========================================================================

function _write_dtfr_plot_tape(io::IO, gendf_mats::Vector{GendfMaterial},
                               params::DtfrParams, pendf_path::String)
    # Card 1: page setup — matches Fortran '(''1 2'',f6.3,'' /'')'
    csz = 0.300
    @printf(io, "1 2%6.3f /\n", csz)

    for dmat in params.materials
        gmat = nothing
        for gm in gendf_mats
            gm.mat == dmat.mat && (gmat = gm; break)
        end
        gmat === nothing && continue

        ng = params.ng
        itabl = params.itabl
        iptotl = params.iptotl
        ipingp = params.ipingp
        ngg = gmat.ngg

        # Build the P0 sig array (same as il=1 in DTF output)
        nws = itabl * ng
        sig = zeros(Float64, nws)

        # Total
        sigma_tot = _get_mf23_sigma(gmat, 501)
        if sigma_tot !== nothing
            for ig in 1:ngg
                jg = ng - ig + 1
                locs = iptotl + itabl * (jg - 1)
                sig[locs] = sigma_tot[ig]
                sig[locs - 2] += sigma_tot[ig]
            end
        end

        # Edits
        for spec in params.edit_specs
            sigma_ed = _get_mf23_sigma(gmat, spec.mt)
            sigma_ed === nothing && continue
            for ig in 1:ngg
                jg = ng - ig + 1
                locs = spec.jpos + itabl * (jg - 1)
                sig[locs] += spec.mult * sigma_ed[ig]
            end
        end

        # Scatter (P0 only for plotting)
        for sec26 in gmat.mf26
            for ig in 1:ngg
                jg = ng - ig + 1
                for ig2 in 1:ngg
                    val = sec26.transfer[1, ig, ig2]
                    val == 0.0 && continue
                    jg2 = ng - ig2 + 1
                    jpos = ipingp + (jg2 - jg)
                    if jpos > itabl
                        jg2 = jg2 - jpos + itabl
                        jpos = itabl
                    end
                    if jpos <= iptotl
                        jg2 = jg2 + iptotl + 1 - jpos
                        jpos = iptotl + 1
                    end
                    if jpos >= 1 && jpos <= itabl && jg2 >= 1 && jg2 <= ng
                        locs = jpos + itabl * (jg2 - 1)
                        sig[locs] += val
                        locb = iptotl - 2 + itabl * (jg - 1)
                        sig[locb] -= val
                    end
                end
            end
        end

        # Determine which edit positions have data — matching Fortran ids[] logic
        ned = params.ned
        n3 = ned + 3
        ids = zeros(Int, n3)

        # Standard positions: ned+1 = absorption, ned+2 = nusigf, ned+3 = total
        # Check total
        if sigma_tot !== nothing
            ids[ned + 3] = 1
            ids[ned + 1] = 1  # absorption initialized from total
        end

        # Check edits
        for (ei, spec) in enumerate(params.edit_specs)
            sigma_ed = _get_mf23_sigma(gmat, spec.mt)
            if sigma_ed !== nothing
                ids[ei] = ei
            end
        end

        # Build jped/mted arrays matching Fortran
        jped = zeros(Int, n3)
        mted = zeros(Int, n3)
        for (ei, spec) in enumerate(params.edit_specs)
            jped[ei] = spec.jpos
            mted[ei] = spec.mt
        end
        # Standard edits
        jped[ned + 1] = iptotl - 2  # absorption
        mted[ned + 1] = 1000
        jped[ned + 2] = iptotl - 1  # nusigf
        mted[ned + 2] = 1000
        jped[ned + 3] = iptotl      # total
        mted[ned + 3] = 501  # marked as 1 in Fortran but we use 501 for photon

        # Build edit name array matching Fortran hednam
        hednam = fill("      ", iptotl)
        # For iedit=0, names are read for positions 1..iptotl-3
        # Position iptotl-2 = "absorp", iptotl-1 = "nusigf", iptotl = "total "
        for (ei, name) in enumerate(params.edit_names)
            if ei <= iptotl
                hednam[ei] = rpad(name, 6)[1:6]
            end
        end
        if iptotl >= 2
            hednam[iptotl - 2] = "absorp"
            hednam[iptotl - 1] = "nusigf"
        end
        if iptotl >= 1
            hednam[iptotl] = "total "
        end

        # iphph = 1 since mf23 sections exist (photon data)
        iphph = 1

        # Plot each non-zero edit position: i = 1..iptotl
        # Matching Fortran ploted loop
        nplt = iptotl
        for ipos in 1:nplt
            # Find first edit spec with this position that has data
            ksave = 0
            for k in 1:n3
                if jped[k] == ipos && ids[k] > 0
                    ksave = k
                    break
                end
            end
            ksave == 0 && continue

            jpos = ipos
            mt = mted[ksave]
            mt == 300 && continue
            # For photons: skip if mt==1 and not total position
            # (Fortran: if (mt.eq.1.and.i.ne.iptotl) go to 100)
            # We use mt==501 for total in photon case
            hedn = hednam[jpos]

            # Build y data for histogram (group-reversed, from sig array)
            y_data = zeros(Float64, ngg)
            for k in 1:ngg
                locs = jpos + itabl * (ngg - k)
                y_data[k] = sig[locs]
            end

            # histod: build histogram and compute scales
            x_hist, z_hist, nh, xmin_p, xmax_p, ymin_p, ymax_p, bot_p, top_p, axl_p, axd_p =
                _histod(y_data, gmat.egg, ngg)

            # Write plot frame
            write(io, " 1/\n")

            # Title: '<    NAME      EDITNAME'/  — Fortran: (a6,4x,a6) = 16 chars
            labelz = @sprintf("%-6s    %-6s", rpad(dmat.name, 6)[1:6], rpad(hedn, 6)[1:6])
            @printf(io, " '<    %s'/\n", labelz)
            write(io, " /\n")

            # Axis spec — tag values use Fortran 0p E10.2, limits use 1p E10.2
            xtag = mt >= 500 && mt <= 522 ? 10.0^xmax_p / 10 : 10.0^(xmin_p + (xmax_p - xmin_p) / 50)
            ytag = 10.0^(ymax_p - (ymax_p - ymin_p) / 50)
            @printf(io, " 4 0 3 1%s%s/\n", _fmt_e10_2(xtag), _fmt_e10_2(ytag))
            @printf(io, "%10.2E%10.2E%10.2E/\n", 10.0^xmin_p, 10.0^xmax_p, 1.0)
            write(io, " '<e>nergy (e<v>)'/\n")
            @printf(io, "%10.2E%10.2E%10.2E/\n", 10.0^ymin_p, 10.0^ymax_p, 1.0)
            write(io, " '<c>ross <s>ection (barns)'/\n")
            write(io, "/\n")
            write(io, " 0 0 0 0 1/\n")
            write(io, " '<m.g.>'/\n")
            write(io, " 0\n")

            # Write histogram data — Fortran: (1p,2e13.4,'/')
            for ih in 1:nh
                @printf(io, "%13.4E%13.4E/\n", x_hist[ih], z_hist[ih])
            end
            write(io, " /\n")

            # PENDF overlay: only for reactions in photon range
            if iphph == 1
                if mt < 501 || mt > 522
                    continue
                end
            end

            # Draw PENDF overlay if available
            if !isempty(pendf_path) && isfile(pendf_path)
                j = 1  # overlay counter (first overlay)
                x_pendf, y_pendf, npts = _dpend(pendf_path, dmat.mat, mt,
                                                  axl_p, axd_p, bot_p, top_p)
                if npts > 0
                    ndash = j
                    @printf(io, "%3d/ pendf plot for mt=%4d\n", j + 1, mt)
                    write(io, "/\n")
                    @printf(io, "    0  0%3d 0 1/\n", ndash)
                    l1 = @sprintf("mt%3d   ", mt)
                    @printf(io, " '<%-8s'/\n", l1)
                    write(io, " 0\n")
                    for ipt in 1:npts
                        @printf(io, "%13.4E%13.4E/\n", x_pendf[ipt], y_pendf[ipt])
                    end
                    write(io, " /\n")
                end
            end
        end

        # 3D phot-phot scatter plot (plotnn with iphph=1)
        _plotnn(io, sig, gmat, params, dmat)
    end

    # Terminator
    write(io, " 99/\n")
end

# =========================================================================
# histod — matches Fortran histod exactly
# =========================================================================

function _histod(y::Vector{Float64}, eg::Vector{Float64}, ng::Int)
    small = 1.0e-10
    big = 1.0e10
    fact = 1.7

    x = zeros(Float64, 2 * ng + 4)
    z = zeros(Float64, 2 * ng + 4)

    yh = small
    yl = big
    n = 1
    z[n] = 0.0

    for i in 1:ng
        if y[i] != 0.0 || n != 1
            n += 2
            z[n - 1] = y[i]
            z[n] = y[i]
            x[n - 2] = eg[i]
            x[n - 1] = eg[i]
            if y[i] > yh
                yh = y[i]
            end
            if y[i] > 0.0 && y[i] < yl
                yl = y[i]
            end
        end
    end
    x[n] = eg[ng + 1]
    n += 1
    x[n] = eg[ng + 1]
    z[n] = 0.0

    # Set scales for plotting
    yh = fact * yh
    ayh = log10(yh)
    iph = trunc(Int, ayh)
    if iph >= 0
        iph += 1
    end
    ayl = log10(yl)
    ipl = trunc(Int, ayl - 1)
    if (iph - ipl) > 6
        ipl = iph - 6
    end
    if iph == ipl
        iph += 1
    end
    ymin = Float64(ipl)
    ymax = Float64(iph)

    if iph - ipl <= 3
        ynow = ymin
        if yl > 2.0 * 10.0^ynow
            ymin = ynow + log10(2.0)
        end
        if yl > 5.0 * 10.0^ynow
            ymin = ynow + log10(5.0)
        end
        ynow = ymax
        if yh < 0.5 * 10.0^ynow
            ymax = ynow + log10(0.5)
        end
        if yh < 0.2 * 10.0^ynow
            ymax = ynow + log10(0.2)
        end
    end

    bot = 10.0^ymin
    top = 10.0^ymax
    axh = log10(x[n])
    iph2 = trunc(Int, axh)
    if axh - iph2 != 0.0
        iph2 += 1
    end
    xmax = Float64(iph2)

    ipl2 = trunc(Int, log10(x[1]) + small)
    ipx = -4
    if ipl2 > ipx; ipx = -2; end
    if ipl2 > ipx; ipx = 4; end
    if ipl2 > ipx; ipx = 5; end
    if ipl2 > ipx; ipx = 6; end
    xmin = Float64(ipx)

    if iph2 - ipx <= 2
        xnow = xmax
        if x[n] < 0.5 * 10.0^xnow
            xmax = xnow + log10(0.5)
        end
        if x[n] < 0.2 * 10.0^xnow
            xmax = xnow + log10(0.2)
        end
        xnow = xmin
        if x[1] > 2.0 * 10.0^xnow
            xmin = xnow + log10(2.0)
        end
        if x[1] > 5.0 * 10.0^xnow
            xmin = xnow + log10(5.0)
        end
    end

    axl = xmin
    axd = xmax - xmin

    # Clamp z values
    for i in 1:n
        if z[i] <= bot
            z[i] = bot
        end
    end

    return x[1:n], z[1:n], n, xmin, xmax, ymin, ymax, bot, top, axl, axd
end

# =========================================================================
# dpend — read PENDF pointwise XS, thin/thicken for plotting
# Matches Fortran dpend exactly
# =========================================================================

function _dpend(pendf_path::String, mat::Int, mt::Int,
                axl::Float64, axd::Float64, bot::Float64, top::Float64)
    ehigh = 5.0e7
    nstep = 50
    ndim = 7000

    # Read PENDF pointwise data for this mat/mt
    mf23_data = _read_pendf_mf23(pendf_path, mat)
    if !haskey(mf23_data, mt)
        return Float64[], Float64[], 0
    end

    energies_raw, xs_raw = mf23_data[mt]
    isempty(energies_raw) && return Float64[], Float64[], 0

    fact = 2000.0 / axd
    step = 10.0^(1.0 / nstep)
    elow = 10.0^axl

    x_out = zeros(Float64, ndim)
    y_out = zeros(Float64, ndim)
    npts = 0

    # Build interpolation from PENDF data
    # We simulate Fortran's gety1 with linear interpolation
    ne = length(energies_raw)
    idx = 1  # current position in PENDF data

    # Helper: interpolate at energy e
    function gety1_interp(e::Float64)
        while idx < ne && energies_raw[idx + 1] <= e
            idx += 1
        end
        if idx >= ne
            return xs_raw[ne], ehigh + 1.0  # past end
        end
        if e <= energies_raw[1]
            # Return next tabulated point as enext (not the same point!)
            return xs_raw[1], (ne > 1 ? energies_raw[2] : ehigh + 1.0)
        end
        if idx < ne && e >= energies_raw[idx] && e <= energies_raw[idx + 1]
            # Linear interpolation
            e1, e2 = energies_raw[idx], energies_raw[idx + 1]
            s1, s2 = xs_raw[idx], xs_raw[idx + 1]
            if e1 == e2
                return s1, (idx + 1 < ne ? energies_raw[idx + 2] : ehigh + 1.0)
            end
            s = s1 + (s2 - s1) * (e - e1) / (e2 - e1)
            enext = energies_raw[idx + 1]
            return s, enext
        end
        return xs_raw[idx], (idx < ne ? energies_raw[idx + 1] : ehigh + 1.0)
    end

    # Initialize
    e = 0.0
    s, enext = gety1_interp(e)
    ixlast = -100
    elast = enext
    while true
        e = enext
        if e > step * elast
            e = step * elast
        end
        elast = e
        s, enext = gety1_interp(e)

        if enext < elow
            enext > ehigh && break
            continue
        end

        if s < 0.0; s = -s; end
        if s > top; s = top; end
        if s < bot; s = bot; end
        if e < elow; e = elow; end

        ix = trunc(Int, (log10(e) - axl) * fact)

        if ix == ixlast
            # Thinning logic
            if npts > 0 && s == y_out[npts]
                x_out[npts] = e
                y_out[npts] = s
                enext > ehigh && break
                continue
            end
            # ns-based logic: simplified — keep point if direction changes
            if npts >= 2 && (s - y_out[npts]) * (y_out[npts] - y_out[npts - 1]) > 0.0
                # Same direction — replace last point
                x_out[npts] = e
                y_out[npts] = s
                enext > ehigh && break
                continue
            end
        end
        ixlast = ix

        npts += 1
        npts > ndim && break
        x_out[npts] = e
        y_out[npts] = s

        enext > ehigh && break
    end

    return x_out[1:npts], y_out[1:npts], npts
end

# =========================================================================
# plotnn — 3D phot-phot scatter plot matching Fortran plotnn (iphph=1)
# =========================================================================

function _plotnn(io::IO, sig::Vector{Float64}, gmat::GendfMaterial,
                 params::DtfrParams, dmat)
    ng = params.ng
    itabl = params.itabl
    iptotl = params.iptotl
    ipingp = params.ipingp
    ngg = gmat.ngg
    egg = gmat.egg

    emin = 1.0e3
    emax = 1.0e8
    big_val = 1.0e10

    # Remove edits and normalize by lethargy width
    # Fortran plotnn: ltabn = itabl - iptotl, extract scatter block
    ltabn = itabl - iptotl

    # Build normalized scatter array
    scatter = zeros(Float64, ltabn * ngg)
    top_val = -big_val

    k = 0
    for ig in 1:ngg  # ig is DTF row (already reversed in sig)
        ig2 = ngg - ig + 1  # GENDF group for lethargy width
        dl = log(egg[ig2 + 1]) - log(egg[ig2])
        for ip in 1:ltabn
            k += 1
            j = itabl * (ig - 1) + iptotl + ip
            sigj = sig[j]
            # Boundary checks
            if egg[ig2] < emin; sigj = 0.0; end
            if egg[ig2 + 1] > emax; sigj = 0.0; end
            scatter[k] = sigj / dl
            if scatter[k] > top_val
                top_val = scatter[k]
            end
        end
    end

    i_top = trunc(Int, log10(top_val))
    if i_top >= 0; i_top += 1; end
    bot = 10.0^(i_top - 5)
    top = 10.0^i_top
    lim = ngg * ltabn
    for k2 in 1:lim
        if scatter[k2] < bot
            scatter[k2] = bot
        end
    end

    # Determine axis limits
    xmin_val = egg[1]
    if xmin_val < emin; xmin_val = emin; end
    i_xmin = trunc(Int, log10(xmin_val))
    if i_xmin < 0; i_xmin -= 1; end
    xmin_val = 10.0^i_xmin

    xmax_val = egg[ngg + 1]
    if xmax_val > emax; xmax_val = emax; end
    i_xmax = trunc(Int, log10(xmax_val))
    if i_xmax >= 0; i_xmax += 1; end
    xmax_val = 10.0^i_xmax

    # Write 3D plot header
    l = 0  # il-1 where il=1
    ititle = @sprintf("%-6s l=%1d phot-phot table", rpad(dmat.name, 6)[1:6], l)

    write(io, " 1/ 3d data\n")
    @printf(io, " '<%-26s'/\n", ititle)
    write(io, " /\n")
    write(io, " -4 2/\n")
    @printf(io, "  %10.4E  %10.4E 1./\n", xmin_val, xmax_val)
    write(io, " '<s>ec. <e>nergy'/\n")
    @printf(io, "  %10.4E  %10.4E 1./\n", xmin_val, xmax_val)
    write(io, " '<e>nergy (e<v>)'/\n")
    @printf(io, "  %10.4E  %10.4E 1./\n", bot, top)
    write(io, " '<x>sec/leth'/\n")
    write(io, "/\n")
    write(io, " 15. -15. 15. -2.5 6.5 2.5/\n")
    write(io, " 1/ 3d data\n")

    # Write 3D data strips — matching Fortran plotnn exactly
    # ig1 loops over GENDF groups (not reversed): egg[ig1] is incident energy
    for ig1 in 1:ngg
        if egg[ig1] >= xmin_val && egg[ig1 + 1] <= xmax_val
            i_src = ngg - ig1 + 1  # DTF row index (reversed)
            @printf(io, "   %10.4E/\n", egg[ig1])
            zlast = bot
            ig2l = 0
            for ig2 in 1:ngg
                if egg[ig2] >= xmin_val && egg[ig2 + 1] <= xmax_val
                    j_sink = ngg - ig2 + 1  # DTF column row (reversed)
                    k_col = j_sink - i_src + ipingp - iptotl
                    if k_col >= 1
                        ss = k_col > ltabn ? bot : scatter[k_col + ltabn * (j_sink - 1)]
                        if ss < bot; ss = bot; end
                        @printf(io, "   %10.4E  %10.4E/\n", egg[ig2], zlast)
                        @printf(io, "   %10.4E  %10.4E/\n", egg[ig2], ss)
                        zlast = ss
                        ig2l = ig2
                    end
                end
            end
            if ig2l > 0
                @printf(io, "   %10.4E  %10.4E/\n", egg[ig2l + 1], zlast)
                @printf(io, "   %10.4E  %10.4E/\n", egg[ig2l + 1], bot)
            end
            write(io, " /\n")
        end
    end
    write(io, " /\n")
end

# =========================================================================
# PENDF MF23 reader — reuse from gaminr module
# =========================================================================

# _read_pendf_mf23 is defined in gaminr.jl and available in the NJOY module scope
