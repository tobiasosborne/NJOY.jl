using NJOY

function run_t01()
    r = reconr("/home/tobiasosborne/Projects/NJOY.jl/test/validation/oracle_cache/test01/run_reconr/tape20"; mat=1306, err=0.005)
    A = Float64(r.mf2.AWR)
    alpha = A / (NJOY.PhysicsConstants.bk * 296.0)
    thnmax = 4.81207e6
    # CRITICAL: broadn_grid must use PARTIALS only (elastic, capture), NOT total.
    # Fortran broadr's broadn uses nreac=2 for C-nat (elastic + capture).
    # Including the total makes convergence stricter → different grid.
    xs_partials = hcat(r.elastic, r.capture)
    b_e, b_xs = NJOY.broadn_grid(r.energies, xs_partials, alpha, 0.005, 0.05, 0.005/20000, thnmax)
    println("broadr: $(length(b_e)) pts")

    # Broaden total separately on the same grid (Fortran broadens all MTs independently)
    # Below thnmax: sigma1_at. Above thnmax: copy reconr values.
    b_total = Vector{Float64}(undef, length(b_e))
    for i in eachindex(b_e)
        if b_e[i] <= thnmax
            b_total[i] = NJOY.sigma1_at(b_e[i], r.energies, r.total, alpha)
        else
            # Above thnmax: interpolate from reconr grid (no broadening)
            idx = searchsortedfirst(r.energies, b_e[i])
            if idx <= 1
                b_total[i] = r.total[1]
            elseif idx > length(r.energies)
                b_total[i] = r.total[end]
            else
                f = (b_e[i] - r.energies[idx-1]) / (r.energies[idx] - r.energies[idx-1])
                b_total[i] = r.total[idx-1] + f * (r.total[idx] - r.total[idx-1])
            end
        end
    end

    override_mf3 = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}(
        1 => (b_e, b_total), 2 => (b_e, b_xs[:,1]), 102 => (b_e, b_xs[:,2]))

    # Read MF12/MT=102 gamma data from ENDF file (C-nat: 3 gammas)
    # Fortran heatr reads these for gheat photon recoil computation
    endf_file = "/home/tobiasosborne/Projects/NJOY.jl/test/validation/oracle_cache/test01/run_reconr/tape20"
    gamma_102 = Tuple{Float64,Float64}[]
    begin
        io = open(endf_file)
        in_mf12_102 = false; subsec = 0; nk = 0
        for line in eachline(io)
            length(line) < 75 && continue
            p = rpad(line, 80)
            mat = NJOY._parse_int(p[67:70])
            mf = NJOY._parse_int(p[71:72])
            mt = NJOY._parse_int(p[73:75])
            mat != 1306 && continue
            if mf == 12 && mt == 102 && !in_mf12_102
                in_mf12_102 = true
                # HEAD: ZA, AWR, LO, 0, NK, 0
                nk = NJOY._parse_int(p[45:55])
                subsec = 0
                continue
            end
            if in_mf12_102
                if mf != 12 || mt != 102
                    break
                end
                # Parse subsections: first is total yield, rest are per-gamma
                # Each subsection starts with a CONT/TAB1 header
                # Look for gamma energy in first field of subsection CONT
                eg = NJOY.parse_endf_float(p[1:11])
                es = NJOY.parse_endf_float(p[12:22])
                lp = NJOY._parse_int(p[23:33])
                lf = NJOY._parse_int(p[34:44])
                nr = NJOY._parse_int(p[45:55])
                np_val = NJOY._parse_int(p[56:66])
                if subsec == 0 && np_val > 0
                    # First subsection = total yield TAB1, skip it
                    subsec += 1
                    # Skip interpolation + data lines
                    n_skip = 1 + cld(2*np_val, 6)
                    for _ in 1:n_skip
                        readline(io)
                    end
                    continue
                elseif eg > 0 && np_val > 0
                    # Gamma subsection: read yield from TAB1 data
                    subsec += 1
                    readline(io) # interpolation line
                    yline = readline(io) # first data line
                    yp = rpad(yline, 80)
                    yield_val = NJOY.parse_endf_float(yp[12:22]) # yield at first energy
                    push!(gamma_102, (eg, yield_val))
                    # Skip remaining data lines
                    n_data = cld(2*np_val, 6) - 1
                    for _ in 1:n_data
                        readline(io)
                    end
                end
            end
        end
        close(io)
    end
    println("MF12/MT=102 gammas: $(length(gamma_102))")
    for (eg, yld) in gamma_102
        println("  E_gamma=$(eg) eV, yield=$yld")
    end

    # Build PointwiseMaterial from partials (elastic=col1, capture=col2)
    pm_broad = PointwiseMaterial(Int32(1306), b_e,
        hcat(b_xs[:,1], zeros(length(b_e)), b_xs[:,2]), [2, 18, 102])
    gd = Dict{Int, Vector{Tuple{Float64,Float64}}}(102 => gamma_102)
    kr = compute_kerma(pm_broad; awr=A, Z=6, gamma_data=gd)
    extra_mf3 = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}(
        301 => (b_e, kr.total_kerma), 444 => (b_e, kr.damage_energy))

    # === THERMR 1: Free gas MT=221 ===
    # KEY INSIGHT: For iinc=1 (free gas), Fortran thermr calcem line 2455: ex(3)=ex(2)
    # MT=221 MF3 = the BROADENED ELASTIC XS on the elastic grid. No calcem integration.
    emax_thermr = 1.2
    # Fortran thermr line 1913-1914: for free gas, smz=1, sb=((az+1)/az)^2
    sb = ((A + 1) / A)^2

    mt221_e = Float64[]; mt221_xs = Float64[]
    for i in eachindex(b_e)
        b_e[i] > emax_thermr && continue
        push!(mt221_e, b_e[i])
        push!(mt221_xs, b_xs[i, 1])  # elastic column (col 1 = elastic)
    end
    # If emax not in grid, add it with interpolated elastic
    if isempty(mt221_e) || mt221_e[end] < emax_thermr
        idx = searchsortedfirst(b_e, emax_thermr)
        if idx <= 1
            xs_emax = b_xs[1, 1]
        elseif idx > length(b_e)
            xs_emax = b_xs[end, 1]
        else
            f = (emax_thermr - b_e[idx-1]) / (b_e[idx] - b_e[idx-1])
            xs_emax = b_xs[idx-1, 1] + f * (b_xs[idx, 1] - b_xs[idx-1, 1])
        end
        push!(mt221_e, emax_thermr); push!(mt221_xs, xs_emax)
    end
    # Cutoff: sigfig(emax, 7, +1) with XS=0, then etop=2e7
    push!(mt221_e, NJOY.round_sigfig(emax_thermr, 7, 1)); push!(mt221_xs, 0.0)
    push!(mt221_e, 2e7); push!(mt221_xs, 0.0)
    extra_mf3[221] = (mt221_e, mt221_xs)
    println("MT=221: $(length(mt221_e)) pts ($(length(mt221_e)-2) below cutoff)")

    # === THERMR 2: S(α,β) from tape26 ===
    sab = NJOY.read_mf7_mt4("/home/tobiasosborne/Projects/NJOY.jl/njoy-reference/tests/resources/t322", 1065, 296.0)

    # Build Bragg data first (needed for thermal grid)
    bragg = NJOY.build_bragg_data(
        a=2.4573e-8, c=6.7e-8, sigma_coh=5.50, A_mass=12.011,
        natom=1, debye_waller=2.1997, emax=emax_thermr, lat=1)

    # Build thermal grid using BROADENED elastic grid (Fortran coh reads from iold = broadened)
    # NOT the reconr grid — broadr thins some points.
    # Fortran iold also contains emax and sigfig(emax,7,+1) as boundary points
    # (from thermr's grid preparation), needed for the 5-point sliding window.
    broadened_grid = sort(unique(vcat(
        Float64[e for e in b_e if e <= emax_thermr],
        [emax_thermr, NJOY.round_sigfig(emax_thermr, 7, 1)])))
    thermal_e = NJOY.build_thermal_grid(bragg, broadened_grid, emax_thermr; tol=0.05)
    # Fortran tpend adds emax and sigfig(emax,7,+1) as final grid points (NP = ne+1)
    push!(thermal_e, emax_thermr)
    push!(thermal_e, NJOY.round_sigfig(emax_thermr, 7, 1))
    push!(thermal_e, 2e7)  # Fortran also includes 2e7 as sentinel
    sort!(unique!(thermal_e))
    println("Thermal grid: $(length(thermal_e)) pts (Fortran coh: 569)")

    # MT=229: Use calcem's total XS (sigl integration), interpolated to thermal grid
    println("Computing calcem for SAB...")
    esi_sab, xsi_sab, mf6_229 = NJOY.calcem(sab, 296.0, emax_thermr, 8; tol=0.05)
    # Order-5 Lagrangian interpolation matching Fortran terp(esi,xsi,nne,enow,nlt=5)
    function terp_lagrange(xs, ys, arg, il=5)
        n = length(xs)
        n == 0 && return 0.0
        il = min(il, n)
        n == il && return _lagrange_eval(xs, ys, 1, il, arg)
        il2 = il ÷ 2
        iadd = il % 2
        # Increasing sequence (Fortran lines 1459-1468)
        ilow = il2 + 1
        ihi = n - il2 - iadd
        iadd = 0  # Fortran line 1468: overwrite for increasing
        iuseh = n - il + 1
        ibeg = ilow + 1
        iend = ihi - 1
        last = iend - il2 + 1
        # Boundary checks (Fortran lines 1481-1498)
        abs(arg - xs[ilow]) < 1e-10 * abs(arg) && return ys[ilow]
        arg <= xs[ilow] && return _lagrange_eval(xs, ys, 1, il, arg)
        abs(xs[ihi] - arg) < 1e-10 * abs(arg) && return ys[ihi]
        arg >= xs[ihi] && return _lagrange_eval(xs, ys, iuseh, il, arg)
        # Search (Fortran lines 1500-1520)
        l = last
        for m in ibeg:iend
            abs(xs[m] - arg) < 1e-10 * abs(arg) && return ys[m]
            if xs[m] > arg
                l = m - il2 + iadd
                break
            end
        end
        l = max(1, min(l, n - il + 1))
        return _lagrange_eval(xs, ys, l, il, arg)
    end
    function _lagrange_eval(xs, ys, l, il, arg)
        s = 0.0
        for i in 1:il
            p = 1.0; pk = 1.0
            idx_n = l + i - 1
            for ip in 1:il
                ip == i && continue
                idx_p = l + ip - 1
                p *= (arg - xs[idx_p])
                pk *= (xs[idx_n] - xs[idx_p])
            end
            s += p * ys[idx_n] / pk
        end
        return s
    end
    mt229_xs = Float64[]
    for e in thermal_e
        if e > emax_thermr; push!(mt229_xs, 0.0)
        else
            push!(mt229_xs, terp_lagrange(esi_sab, xsi_sab, e, 5))
        end
    end
    extra_mf3[229] = (thermal_e, mt229_xs)

    # MT=230: Bragg edges on thermal grid — sigma_coh=5.50 (Fortran gr4), DW=2.1997 at 296K
    # Round to 7 sigfigs: Fortran sigcoh outputs 7-sigfig precision (a11 format)
    mt230_xs = [e > emax_thermr ? 0.0 : NJOY.round_sigfig(NJOY.bragg_edges(e, bragg), 7, 0) for e in thermal_e]
    extra_mf3[230] = (thermal_e, mt230_xs)

    # === MF6 ===
    esi_fg, xsi_fg, mf6_221 = NJOY.calcem_free_gas(A, 296.0, 1.2, 8; sigma_b=sb, tol=0.05)
    mf6 = Dict{Int, Any}(221 => mf6_221)
    mf6[229] = mf6_229
    mf6_xsi = Dict{Int, Vector{Float64}}(221 => xsi_fg, 229 => xsi_sab)
    mf6_emax = Dict{Int, Float64}(221 => 1.2, 229 => 1.2)
    mf6_stubs = Dict{Int, NamedTuple}(230 => (nbragg=bragg.n_edges, emin=1e-5, emax=emax_thermr))

    # === MF12/MF13 ===
    mf12 = String[]; mf13 = String[]
    for line in readlines("/home/tobiasosborne/Projects/NJOY.jl/test/validation/oracle_cache/test01/after_reconr.pendf")
        length(line) < 75 && continue
        p = rpad(line, 80); mf = NJOY._parse_int(p[71:72]); mt = NJOY._parse_int(p[73:75])
        mf == 12 && mt > 0 && push!(mf12, line)
        mf == 13 && mt > 0 && push!(mf13, line)
    end

    # === WRITE ===
    open("/tmp/t01_tape25.pendf", "w") do io
        NJOY.write_full_pendf(io, r; mat=1306, label="pendf tape for c-nat from endf/b tape 511",
            err=0.005, tempr=296.0, override_mf3=override_mf3, extra_mf3=extra_mf3,
            mf6_records=mf6, mf6_stubs=mf6_stubs, mf12_lines=mf12, mf13_lines=mf13,
            mf6_xsi=mf6_xsi, mf6_emax=mf6_emax,
            thermr_mts=Set([221, 229, 230]),
            descriptions=["6-c-nat from tape 511",
                          "processed by the njoy nuclear data processing system",
                          "see original endf/b-v tape for details of evaluation"])
    end

    # === COMPARE ===
    ref = "/home/tobiasosborne/Projects/NJOY.jl/njoy-reference/tests/01/referenceTape25"
    println("\ntape25: $(countlines("/tmp/t01_tape25.pendf")) (ref: $(countlines(ref)))")

    function list_sections(fn)
        s = Dict{Tuple{Int,Int}, Int}()
        for line in readlines(fn)
            length(line) < 75 && continue
            p = rpad(line, 80); mf = NJOY._parse_int(p[71:72]); mt = NJOY._parse_int(p[73:75])
            mat = NJOY._parse_int(p[67:70])
            mf > 0 && mt > 0 && mat > 0 && (s[(mf,mt)] = get(s, (mf,mt), 0) + 1)
        end; s
    end
    j = list_sections("/tmp/t01_tape25.pendf"); f = list_sections(ref)
    n_match = 0; n_total = 0
    println("\nSection comparison:")
    for k in sort(collect(union(keys(j), keys(f))))
        jc = get(j, k, 0); fc = get(f, k, 0); n_total += 1
        jc == fc ? (n_match += 1) :
            println("  MF=$(k[1]) MT=$(lpad(k[2],3)): Julia=$jc Fortran=$fc")
    end
    println("MATCH: $n_match / $n_total sections")

    # Data comparison — skip headers (first 3 lines), compare data only
    println("\n=== MF3 data comparison (data lines only, skipping headers) ===")
    jall = readlines("/tmp/t01_tape25.pendf"); fall = readlines(ref)
    total_match = 0; total_lines = 0
    for mt in [221, 229, 230, 1, 2, 102, 301, 444, 4, 51, 91, 103]
        jmt = String[]; fmt = String[]
        for line in jall
            length(line) < 75 && continue
            p = rpad(line, 80)
            NJOY._parse_int(p[71:72]) == 3 && NJOY._parse_int(p[73:75]) == mt &&
                NJOY._parse_int(p[67:70]) > 0 && push!(jmt, p[1:66])
        end
        for line in fall
            length(line) < 75 && continue
            p = rpad(line, 80)
            NJOY._parse_int(p[71:72]) == 3 && NJOY._parse_int(p[73:75]) == mt &&
                NJOY._parse_int(p[67:70]) > 0 && push!(fmt, p[1:66])
        end
        isempty(jmt) && isempty(fmt) && continue
        if length(jmt) == length(fmt)
            # Skip first 3 lines (headers)
            nh = min(3, length(jmt))
            data_j = jmt[nh+1:end]; data_f = fmt[nh+1:end]
            exact = count(i -> data_j[i] == data_f[i], 1:length(data_j))
            hdr_match = count(i -> jmt[i] == fmt[i], 1:nh)
            total_match += exact + hdr_match; total_lines += length(jmt)
            pct = length(data_j) > 0 ? round(100*exact/length(data_j), digits=1) : 100.0
            if exact == length(data_j)
                println("  MT=$(lpad(mt,3)): DATA $(exact)/$(length(data_j)) PERFECT (hdr $hdr_match/$nh)")
            else
                println("  MT=$(lpad(mt,3)): DATA $exact/$(length(data_j)) ($pct%) hdr=$hdr_match/$nh")
                for i in 1:length(data_j)
                    if data_j[i] != data_f[i]
                        println("    first diff @$i: J=$(data_j[i])")
                        println("                    F=$(data_f[i])")
                        break
                    end
                end
            end
        else
            println("  MT=$(lpad(mt,3)): SIZE MISMATCH $(length(jmt)) vs $(length(fmt))")
        end
    end
    println("\nTotal data: $total_match / $total_lines ($(round(100*total_match/total_lines, digits=1))%)")
end

run_t01()
