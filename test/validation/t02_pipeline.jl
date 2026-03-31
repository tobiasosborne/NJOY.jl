using NJOY

function run_t02()
    base = @__DIR__
    project_root = dirname(dirname(base))
    oracle_dir = joinpath(base, "oracle_cache", "test02")
    tape20 = joinpath(oracle_dir, "run_reconr", "tape20")

    # ═══════════════════════════════════════════════════════════
    # 1. RECONR — 0K pointwise reconstruction
    # ═══════════════════════════════════════════════════════════
    # T02: Pu-238, MAT=1050, err=0.005, SLBW + URR mode=11 (LSSF=0)
    # Expected: 17/17 MTs BIT-IDENTICAL, 3567 energy points
    println("=== RECONR ===")
    r = reconr(tape20; mat=1050, err=0.005)
    A = Float64(r.mf2.AWR)
    println("  AWR = $A")
    println("  energies: $(length(r.energies)) pts (Fortran: 3567)")
    println("  MTs: $(length(r.mf3_sections))")

    # Verify reconr against oracle
    reconr_oracle = joinpath(oracle_dir, "after_reconr.pendf")
    write_pendf_file("/tmp/t02_reconr.pendf", r; mat=1050, err=0.005)
    _compare_mf3("RECONR", "/tmp/t02_reconr.pendf", reconr_oracle)

    # ═══════════════════════════════════════════════════════════
    # 2. BROADR — Doppler broadening at 3 temperatures
    # ═══════════════════════════════════════════════════════════
    # From Fortran output:
    #   thnmax = 200.0 eV (resolved resonance range upper limit)
    #   Broadened MTs: 2, 18, 102 (elastic, fission, capture)
    #   Pass 1: 0K → 300K:   3567 → 2925 pts
    #   Pass 2: 300K → 900K:  2925 → 2592 pts
    #   Pass 3: 900K → 2100K: 2592 → 2418 pts
    # CRITICAL: broadn_grid uses PARTIALS only for convergence (Trap 53).
    # Total is broadened separately via sigma1_at on the same grid.
    println("\n=== BROADR ===")
    thnmax = 200.0
    temps = [300.0, 900.0, 2100.0]
    bk = NJOY.PhysicsConstants.bk

    # Reconr data — partials for broadn convergence: elastic, fission, capture
    prev_e = r.energies
    prev_partials = hcat(r.elastic, r.fission, r.capture)
    prev_total = r.total
    T_old = 0.0

    broadened_results = []  # (temp, b_e, b_total, b_elastic, b_fission, b_capture)

    for T_new in temps
        T_eff = T_new - T_old
        alpha = A / (bk * T_eff)
        println("\n  Broadening $(T_old)K → $(T_new)K  (T_eff=$(T_eff), alpha=$(round(alpha, sigdigits=6)))")

        # Broaden partials via broadn_grid (Trap 53: partials only)
        b_e, b_xs = NJOY.broadn_grid(prev_e, prev_partials, alpha,
                                       0.005, 0.05, 0.005/20000, thnmax)
        println("    Grid: $(length(prev_e)) → $(length(b_e)) pts")

        # Broaden total separately via sigma1_at on the same grid
        b_total = Vector{Float64}(undef, length(b_e))
        for i in eachindex(b_e)
            if b_e[i] <= thnmax
                b_total[i] = NJOY.sigma1_at(b_e[i], prev_e, prev_total, alpha)
            else
                # Above thnmax: interpolate from previous grid (no broadening)
                idx = searchsortedfirst(prev_e, b_e[i])
                if idx <= 1
                    b_total[i] = prev_total[1]
                elseif idx > length(prev_e)
                    b_total[i] = prev_total[end]
                else
                    f = (b_e[i] - prev_e[idx-1]) / (prev_e[idx] - prev_e[idx-1])
                    b_total[i] = prev_total[idx-1] + f * (prev_total[idx] - prev_total[idx-1])
                end
            end
        end

        # Fortran broadr writes XS through standard a11 which rounds to 7 sigfigs.
        # The reconr output has 9-sigfig precision, but broadened values from sigma1_at
        # have full Float64 precision. Rounding to 7 sigfigs matches Fortran's output.
        b_total_r = NJOY.round_sigfig.(b_total, 7, 0)
        b_el_r = NJOY.round_sigfig.(b_xs[:, 1], 7, 0)
        b_fi_r = NJOY.round_sigfig.(b_xs[:, 2], 7, 0)
        b_ca_r = NJOY.round_sigfig.(b_xs[:, 3], 7, 0)

        push!(broadened_results, (T_new, b_e, b_total_r, b_el_r, b_fi_r, b_ca_r))

        # Prepare for next pass: input = this pass's output
        prev_e = b_e
        prev_partials = b_xs
        prev_total = b_total
        T_old = T_new
    end

    # ═══════════════════════════════════════════════════════════
    # 3. COMPARE against broadr oracle — per-temperature, per-MT
    # ═══════════════════════════════════════════════════════════
    println("\n=== COMPARISON vs after_broadr.pendf ===")
    broadr_oracle = joinpath(oracle_dir, "after_broadr.pendf")

    # Parse oracle: extract each temperature block's MF3 data
    oracle_sections = _parse_broadr_oracle(broadr_oracle, 1050)

    # For each temperature, compare broadened MTs (1, 2, 18, 102)
    # and non-broadened MTs (4, 16, 17, 51-59, 91)
    broadened_mts = Set([1, 2, 18, 102])
    total_match = 0; total_lines = 0

    for (tidx, (T, b_e, b_total, b_el, b_fi, b_ca)) in enumerate(broadened_results)
        println("\n  --- Temperature $tidx: $(T)K ---")
        println("    Grid: $(length(b_e)) pts")

        for mt in sort(collect(keys(oracle_sections[tidx])))
            oracle_data = oracle_sections[tidx][mt]
            if mt in broadened_mts
                # Get Julia data for this MT
                julia_xs = if mt == 1
                    b_total
                elseif mt == 2
                    b_el
                elseif mt == 18
                    b_fi
                elseif mt == 102
                    b_ca
                end
                julia_e = b_e
                julia_lines = _format_tab1_data(julia_e, julia_xs)
            else
                # Non-broadened MT: use reconr data (same at all temps)
                legacy = NJOY._get_legacy_section(r, mt)
                if legacy !== nothing
                    julia_lines = _format_tab1_data(legacy[1], legacy[2])
                else
                    println("    MT=$(lpad(mt,3)): MISSING in Julia")
                    continue
                end
            end

            # Compare data lines (skip headers)
            n_oracle = length(oracle_data)
            n_julia = length(julia_lines)
            if n_oracle == n_julia
                exact = count(i -> julia_lines[i] == oracle_data[i], 1:n_oracle)
                total_match += exact; total_lines += n_oracle
                if exact == n_oracle
                    println("    MT=$(lpad(mt,3)): DATA $exact/$n_oracle PERFECT")
                else
                    pct = round(100*exact/n_oracle, digits=1)
                    println("    MT=$(lpad(mt,3)): DATA $exact/$n_oracle ($pct%)")
                    # Show first diff
                    for i in 1:n_oracle
                        if julia_lines[i] != oracle_data[i]
                            println("      @$i: J=$(julia_lines[i])")
                            println("          F=$(oracle_data[i])")
                            break
                        end
                    end
                end
            else
                total_lines += max(n_oracle, n_julia)
                println("    MT=$(lpad(mt,3)): SIZE MISMATCH Julia=$n_julia Fortran=$n_oracle")
            end
        end
    end
    println("\n  Total data: $total_match / $total_lines ($(round(100*total_match/max(total_lines,1), digits=1))%)")

    # ═══════════════════════════════════════════════════════════
    # 4. Tolerance test — simulate execute.py methodology
    # ═══════════════════════════════════════════════════════════
    println("\n=== TOLERANCE TEST ===")
    fp = r"([-+]?\d+(\.\d+)?(([eE])?[-+]\d+)?)"
    function makeFloat(m)
        s = m.match
        try return parse(Float64, s)
        catch
            # Handle ENDF compact format: "1.234567+4" → "1.234567E+4"
            for i in 2:length(s)-1
                if s[i] in ('+', '-') && s[i-1] != 'E' && s[i-1] != 'e' && isdigit(s[i-1])
                    return parse(Float64, s[1:i-1] * "E" * s[i:end])
                end
            end
            return nothing
        end
    end
    # Collect all broadened data lines for tolerance comparison
    # Build parallel arrays of (julia_val, fortran_val) for each broadened line
    for rtol in [1e-9, 1e-7, 1e-5, 1e-4, 1e-3]
        fails = 0
        for tidx in 1:3
            T, b_e, b_total, b_el, b_fi, b_ca = broadened_results[tidx]
            for mt in sort(collect(keys(oracle_sections[tidx])))
                oracle_data = oracle_sections[tidx][mt]
                if mt in broadened_mts
                    julia_xs = mt == 1 ? b_total : mt == 2 ? b_el : mt == 18 ? b_fi : b_ca
                    julia_lines = _format_tab1_data(b_e, julia_xs)
                else
                    legacy = NJOY._get_legacy_section(r, mt)
                    legacy === nothing && continue
                    julia_lines = _format_tab1_data(legacy[1], legacy[2])
                end
                n = min(length(julia_lines), length(oracle_data))
                for i in 1:n
                    julia_lines[i] == oracle_data[i] && continue
                    jm = collect(eachmatch(fp, julia_lines[i]))
                    fm = collect(eachmatch(fp, oracle_data[i]))
                    length(jm) != length(fm) && (fails += 1; continue)
                    ok = true
                    for (jf, ff) in zip(jm, fm)
                        jv = makeFloat(jf); fv = makeFloat(ff)
                        (jv === nothing || fv === nothing) && continue
                        if !isapprox(jv, fv; rtol=rtol, atol=1e-10)
                            ok = false; break
                        end
                    end
                    ok || (fails += 1)
                end
            end
        end
        println("  rel_tol=$(rpad(string(rtol), 5)): $fails lines fail")
    end
end

# ═══════════════════════════════════════════════════════════
# Helper: format energy/XS arrays as ENDF TAB1 data lines (cols 1-66)
# ═══════════════════════════════════════════════════════════
function _format_tab1_data(energies, xs)
    lines = String[]
    vals = Float64[]
    for i in eachindex(energies)
        push!(vals, energies[i])
        push!(vals, xs[i])
    end
    # Format 6 values per line
    for start in 1:6:length(vals)
        line = ""
        for j in start:min(start+5, length(vals))
            line *= NJOY.format_endf_float(vals[j])
        end
        push!(lines, rpad(line, 66))
    end
    return lines
end

# ═══════════════════════════════════════════════════════════
# Helper: parse broadr oracle into per-temperature, per-MT data lines
# ═══════════════════════════════════════════════════════════
function _parse_broadr_oracle(filename, target_mat)
    lines = readlines(filename)
    # Structure: 3 temperature blocks, each has MF1 + MF2 + MF3 sections
    # We only care about MF3 data lines (cols 1-66)
    # Each MF3/MT section: HEAD(1) + CONT(1) + interp(1) + data lines
    # We extract data lines only (after first 3 lines of each section)

    # First pass: find temperature block boundaries
    # MF1/MT=451 marks the start of each temperature block
    temp_starts = Int[]
    for (i, line) in enumerate(lines)
        length(line) < 75 && continue
        p = rpad(line, 80)
        mat = NJOY._parse_int(p[67:70])
        mf = NJOY._parse_int(p[71:72])
        mt = NJOY._parse_int(p[73:75])
        if mat == target_mat && mf == 1 && mt == 451
            # Check if this is the first line of MF1/MT=451 (has NXC in field 6)
            # Each MF1/MT=451 block is ~25 lines; use a large gap to avoid double-counting
            if isempty(temp_starts) || i > temp_starts[end] + 100
                push!(temp_starts, i)
            end
        end
    end
    push!(temp_starts, length(lines) + 1)  # sentinel
    println("  Oracle: $(length(temp_starts)-1) temperature blocks found")

    # Second pass: extract MF3 data lines per temperature block
    result = [Dict{Int, Vector{String}}() for _ in 1:length(temp_starts)-1]

    for tidx in 1:length(temp_starts)-1
        block_start = temp_starts[tidx]
        block_end = temp_starts[tidx+1] - 1

        # Parse MF3 sections within this block
        i = block_start
        while i <= block_end
            length(lines[i]) < 75 && (i += 1; continue)
            p = rpad(lines[i], 80)
            mat = NJOY._parse_int(p[67:70])
            mf = NJOY._parse_int(p[71:72])
            mt = NJOY._parse_int(p[73:75])

            if mat == target_mat && mf == 3 && mt > 0
                # Found MF3/MT section start — skip 3 header lines
                i += 3
                data = String[]
                while i <= block_end
                    length(lines[i]) < 75 && (i += 1; continue)
                    pd = rpad(lines[i], 80)
                    mfc = NJOY._parse_int(pd[71:72])
                    mtc = NJOY._parse_int(pd[73:75])
                    matc = NJOY._parse_int(pd[67:70])
                    (matc != target_mat || mfc != 3 || mtc != mt) && break
                    push!(data, pd[1:66])
                    i += 1
                end
                result[tidx][mt] = data
            else
                i += 1
            end
        end
    end
    return result
end

# ═══════════════════════════════════════════════════════════
# Helper: compare MF3 data (reconr verification)
# ═══════════════════════════════════════════════════════════
function _compare_mf3(label, julia_file, fortran_file)
    function parse_all_mf3(fn)
        result = Dict{Int, Vector{String}}()
        lines = readlines(fn); idx = 1
        while idx <= length(lines)
            length(lines[idx]) < 75 && (idx += 1; continue)
            p = rpad(lines[idx], 80)
            mf = NJOY._parse_int(p[71:72]); mt = NJOY._parse_int(p[73:75])
            mat = NJOY._parse_int(p[67:70])
            if mf == 3 && mt > 0 && mat > 0
                idx += 3; data = String[]
                while idx <= length(lines)
                    length(lines[idx]) < 75 && (idx += 1; continue)
                    pd = rpad(lines[idx], 80)
                    mfc = NJOY._parse_int(pd[71:72]); mtc = NJOY._parse_int(pd[73:75])
                    (mfc != 3 || mtc != mt) && break
                    push!(data, pd[1:66]); idx += 1
                end; result[mt] = data
            else; idx += 1; end
        end; return result
    end
    j = parse_all_mf3(julia_file); f = parse_all_mf3(fortran_file)
    mts = sort(collect(union(keys(j), keys(f)))); np = 0
    for mt in mts
        jd = get(j, mt, String[]); fd = get(f, mt, String[])
        if jd == fd; np += 1
        elseif isempty(jd); println("  $label MT=$(lpad(mt,3)): MISSING")
        else
            d = findfirst(i -> i > length(jd) || i > length(fd) || jd[i] != fd[i],
                          1:max(length(jd),length(fd)))
            println("  $label MT=$(lpad(mt,3)): DIFF $(length(jd))v$(length(fd)) @$d")
            d !== nothing && d <= min(length(jd),length(fd)) &&
                println("    J: $(jd[d])\n    F: $(fd[d])")
        end
    end
    println("  $label: $np/$(length(mts)) PERFECT")
end

run_t02()
