# wimsr resint port — effective resonance integrals for WIMS.
#
# Ref: njoy-reference/src/wimsr.f90:425-676 (subroutine resint).
#
# Walks the GENDF tape per (temperature × sigma_zero × resonance_group)
# and accumulates capture (MT=102) + fission (MT=18) into raw absorption.
# Applies the Goldstein-Cohen normalization:
#
#     σ_b = σ_0 + λ·σ_p
#     σ_abs,RI = σ_b · σ_a / (σ_b + σ_a)
#     σ_nuf,RI = σ_b · ν·σ_f / (σ_b + σ_a)
#
# (Fortran lines 656-668; siglam = spot[nfg+jg]·glam[jg].)
#
# Output: per resonance group, a flat layout of (ires × nsigz) integrals
# plus the temperatures + sigma_zero list, suitable for WIMSResonanceTable.

"""
    wimsr_extract_resint(in_path, p::WimsrParams, xs_data::WimsrXSecs)
        -> (abs_tables, fis_tables)

Compute resonance integral tables. Returns two `Vector{WIMSResonanceTable}`
of length `nrg`: absorption tables (ifis>=1) and fission tables (ifis==3).
For `ifis<3`, the second vector is empty.

Pre-conditions:
- `xs_data.tempr` has at least `p.ires` entries (one per RI temperature)
- `xs_data.spot[nfg+1..nfg+nrg]` is populated (sigp = potential xs)
- `xs_data.snu[nfg+1..nfg+nrg]` is populated (nu-bar from MF=3/MT=452)
- `xs_data.abs2[nfg+1..nfg+nrg]` carries the non-MT=102 absorption seed
"""
function wimsr_extract_resint(in_path::AbstractString, p::WimsrParams,
                              xs_data::WimsrXSecs)
    ngnd = p.ngnd
    nfg  = p.nfg
    nrg  = p.nrg
    ires = p.ires

    # Validate temperature axis BEFORE slicing — converts the otherwise opaque
    # `BoundsError [1:ires]` on `xs_data.tempr[1:ires]` (line ~275 below) into
    # a self-explaining error that names the upstream gap. Currently
    # src/orchestration/modules/groupr.jl:31 uses only
    # `params.temperatures[1]`, so multi-temperature decks (T11 Pu-238 with
    # 3 broadr temperatures) hit this path. Fixing requires porting groupr's
    # multi-temperature loop (Fortran groupr.f90 outer-temp do-loop).
    length(xs_data.tempr) >= ires ||
        error("wimsr_extract_resint: GENDF tape has $(length(xs_data.tempr)) " *
              "temperature(s) but params.ires=$ires. " *
              "Julia groupr (src/orchestration/modules/groupr.jl:31) currently " *
              "writes only the first temperature; the multi-temperature " *
              "groupr port is required before wimsr can run on multi-T decks.")

    # Resonance group range in GENDF top-down order:
    #   jg = ngnd - ig + 1, jg ∈ [nfg+1, nfg+nrg] (resonance window)
    # Equivalent ig: ig ∈ [ngnd-nfg-nrg+1, ngnd-nfg]
    nghi = ngnd - nfg              # max ig in resonance range
    nglo = nghi - nrg + 1          # min ig in resonance range

    # Determine nsigz from first MF=1/MT=451 LIST (use xs_data egb/tempr already populated).
    nsigz_eff = p.nsigz == 0 ? 7 : p.nsigz   # T11 uses 7

    # Per-(temp, sig0, group) buffers (1-based, jg in 1..nrg)
    sabs = zeros(ires * nsigz_eff * nrg)   # absorption raw (capture + fission)
    snsf = zeros(ires * nsigz_eff * nrg)   # fission raw
    elas = zeros(ires * nsigz_eff * nrg)   # elastic raw (for diagnostic)

    # Seed sabs with non-MT=102 background absorption (Fortran line 476)
    for jg in 1:nrg, it in 1:ires, jz in 1:nsigz_eff
        idx = jz + nsigz_eff * ((it - 1) + ires * (jg - 1))
        sabs[idx] = xs_data.abs2[nfg + jg]
    end

    # Sigma_zero list (constant across temps) — descending in groupr GENDF
    # (sig0[1] = inf-dil = 1e10, sig0[nz] = most-shielded). Fortran reverses
    # so sig0[1] becomes most-shielded after the `sigz(nsigz-i+1)=scr(...)`
    # at line 497.
    sigz = zeros(nsigz_eff)

    has_fission = false

    lines = readlines(in_path)
    cur_temp = 0
    section_nl = 1
    section_nz = 1
    section_mt = 0
    section_mf = 0
    in_target = false

    idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        row = rpad(lines[idx], 80)
        mv = tryparse(Int, strip(row[67:70]))
        mf = tryparse(Int, strip(row[71:72]))
        mt = tryparse(Int, strip(row[73:75]))
        seq = tryparse(Int, strip(row[76:80]))
        (mv === nothing || mf === nothing || mt === nothing || seq === nothing) && (idx += 1; continue)

        if mv == p.mat && mf == 1 && mt == 451 && seq == 1
            cur_temp += 1
            in_target = (cur_temp <= ires)
            cur_temp > ires && break
        end

        # MF=1/MT=451 row 2 contains nsigz at L2; subsequent LIST has the sig0 list.
        # Fortran reads `sigz(nsigz-i+1)=scr(6+ntw+i)` — reverses order so
        # sigz[1] = most-shielded, sigz[nsigz] = inf-dil.
        if mv == p.mat && mf == 1 && mt == 451 && seq == 2 && in_target && cur_temp == 1
            # The actual sig0 values live in the LIST body below; read them.
            # GENDF body: pos 1 = sentinel 0, pos 2 = inf-dil sigref, pos 3..nz+1 = remaining sig0.
            body_start_idx = idx + 1
            body = Float64[]
            bidx = body_start_idx
            target_count = nsigz_eff + 1   # nz sig0 + sentinel
            while bidx <= length(lines) && length(body) < target_count + 5
                brow = rpad(lines[bidx], 80)
                mt_check = tryparse(Int, strip(brow[73:75]))
                mt_check == 451 || break
                seq_check = tryparse(Int, strip(brow[76:80]))
                seq_check == 99999 && break
                for col in 0:5
                    s = brow[1+11*col : 11+11*col]
                    all(isspace, s) && break
                    push!(body, parse_endf_float(s))
                end
                bidx += 1
            end
            # Extract sig0 list: positions 2..nsigz+1 are sig0 in descending order.
            # Fortran reverses → sigz[i] = sig0[nz-i+1] (Fortran 1-based).
            if length(body) >= nsigz_eff + 1
                for i in 1:nsigz_eff
                    sigz[nsigz_eff - i + 1] = body[1 + i]   # body[2..nz+1] reversed
                end
            end
            idx += 1
            continue
        end

        if mv == p.mat && mf in (3, 6) && mt > 0 && seq == 1 && in_target
            section_nl = something(tryparse(Int, strip(row[23:33])), 1)
            section_nz = something(tryparse(Int, strip(row[34:44])), 1)
            section_mt = mt
            section_mf = mf
            if mt == 18; has_fission = true; end
            idx += 1
            continue
        end

        # Per-group LIST records in MF=3 (resint only consumes MF=3, line 519)
        if mv == p.mat && mf == 3 && mt > 0 && seq >= 2 && in_target &&
           mt == section_mt && mf == section_mf
            ig    = tryparse(Int, strip(row[56:66]))
            nw    = tryparse(Int, strip(row[45:55]))
            (ig === nothing || nw === nothing) && (idx += 1; continue)
            (!(1 <= ig <= ngnd) || nw <= 0) && (idx += 1; continue)

            # Restrict to resonance-group range (Fortran line 533)
            if !(nglo <= ig <= nghi)
                idx += 1; continue
            end
            # Skip irrelevant MTs (Fortran line 540)
            if !(mt == 1 || mt == 2 || mt == 18 || mt == 102)
                idx += 1; continue
            end

            body, idx_after = _wimsr_read_list_body(lines, idx + 1, mt, nw)
            if length(body) < nw
                idx = idx_after; continue
            end

            kg = ngnd - ig + 1
            jg = kg - nfg                  # 1..nrg
            (jg < 1 || jg > nrg) && (idx = idx_after; continue)

            nl = section_nl
            nz_rec = section_nz
            lim = min(nsigz_eff, nz_rec)

            # Body indexing in Fortran: scr(nl*jz + loca) where loca = l+lz+nl*(nz-1)
            # for absorption/fission/elastic. In our body coords (1-based, no header):
            # body[nl*jz + (lz+nl*(nz-1)) - lz] = body[nl*jz + nl*(nz-1)]
            # Actually Fortran scr(loca) → our body[loca-l-lz+1] = body[loca-l-5].
            # For loca = l+lz+nl*(nz-1) and the Fortran read scr(nl*jz+loca):
            #   our body index = nl*jz + nl*(nz-1) + 1 = nl*(jz + nz - 1) + 1
            # That maps to (i_l=1, iz=jz, ig2=2) in our _bidx convention:
            #   _bidx(1, jz, 2, nl, nz_rec) = 1 + nl*(jz-1) + nl*nz*1
            #                               = 1 + nl*(jz-1) + nl*nz
            # Different from the loca read! Let me re-derive.
            # Actually Fortran's MF=3 layout: 1st nl*nz block = flux (ig2=1),
            # 2nd nl*nz block = xs (ig2=2). So scr(nl*jz+loca) with
            # loca = l + lz + nl*(nz-1) means:
            #   scr(l+lz+nl*(nz-1)+nl*jz) — NOTE this is offset by nl*(nz-1)
            #   into the flux block! Plus nl*jz pushes to the xs block.
            #   Fortran code reads NL elements per (jz, ig2=2) starting from
            #   the LAST sig0's flux entry, then jumping by nl per jz.
            # Our 1-based body index: pos = nl*jz + nl*(nz-1) + 1
            # For nl=2 nz=7 jz=1: pos = 2 + 12 + 1 = 15 — first xs entry.
            # For jz=2: pos = 17. So body[15], body[17], ... — every nl
            # position. That's body[(i_l=1, iz=jz, ig2=2)] for jz=1..lim.

            # Initial offset before any jz iter: pos = nl*(nz-1) + 1 = body[nl*(nz-1)+1]
            # Then for each jz (1..lim), pos += nl. body[(jz=1)] = body[nl*nz+1].
            # Wait: pos for jz=1 = nl*1 + nl*(nz-1) + 1 = nl*nz + 1. body[nl*nz+1] = first
            # entry of ig2=2 (xs block). ✓.

            offset_base = nl * (nz_rec - 1)
            iadd = nsigz_eff + nsigz_eff * ((cur_temp - 1) + ires * (jg - 1))

            if mt == 102
                for jz in 1:lim
                    pos = nl * jz + offset_base + 1   # 1-based body index
                    if pos <= length(body)
                        sabs[iadd - jz + 1] += body[pos]
                    end
                end
            elseif mt == 18
                has_fission = true
                for jz in 1:lim
                    pos = nl * jz + offset_base + 1
                    if pos <= length(body)
                        v = body[pos]
                        snsf[iadd - jz + 1] += v
                        sabs[iadd - jz + 1] += v   # fission counts as abs
                    end
                end
            elseif mt == 2
                for jz in 1:lim
                    pos = nl * jz + offset_base + 1
                    if pos <= length(body)
                        elas[iadd - jz + 1] += body[pos]
                    end
                end
            end
            # MT=1 flux extraction (Fortran lines 552-567) — TODO Phase 58c
            # (needed for true cross-temperature flux ratio in resint).

            idx = idx_after
            continue
        end

        idx += 1
    end

    # Compute snux = snsf * snu (Fortran lines 607-617) — only if fission seen
    snux = zeros(ires * nsigz_eff * nrg)
    if has_fission
        for jg in 1:nrg
            locn = nfg + jg
            for jz in 1:nsigz_eff, jtem in 1:ires
                ioff = jz + nsigz_eff * ((jtem - 1) + ires * (jg - 1))
                snux[ioff] = snsf[ioff] * xs_data.snu[locn]
            end
        end
    end

    # Goldstein-Cohen normalization (Fortran lines 656-668)
    for jg in 1:nrg
        siglam = xs_data.spot[nfg + jg] * (jg <= length(xs_data.l1) ? 1.0 : 1.0)
        # spot * glam — but glam is per-resonance-group; xs_data has lambdas
        # only inside the WIMSMaterial; resint takes glam from p.glam directly
        glam_jg = jg <= length(p.glam) ? p.glam[jg] : 1.0
        siglam = xs_data.spot[nfg + jg] * glam_jg

        for it in 1:ires
            for iz in 1:nsigz_eff
                idx = iz + nsigz_eff * ((it - 1) + ires * (jg - 1))
                sigb = sigz[iz] + siglam
                siga = sabs[idx]
                sig  = snux[idx]
                if sigb + siga > 0
                    sabs[idx] = sigb * siga / (sigb + siga)
                    snux[idx] = sigb * sig  / (sigb + siga)
                end
            end
        end
    end

    # Build WIMSResonanceTable vectors
    abs_tables = WIMSResonanceTable[]
    fis_tables = WIMSResonanceTable[]
    for jg in 1:nrg
        # Per-group integrals: (ires × nsigz) flat in (it, iz) order
        absorb = Vector{Float64}(undef, ires * nsigz_eff)
        for it in 1:ires, iz in 1:nsigz_eff
            src_idx = iz + nsigz_eff * ((it - 1) + ires * (jg - 1))
            dst_idx = iz + nsigz_eff * (it - 1)
            absorb[dst_idx] = sabs[src_idx]
        end
        push!(abs_tables, WIMSResonanceTable(p.rdfid, ires, nsigz_eff,
                                              xs_data.tempr[1:ires],
                                              copy(sigz),
                                              absorb))
        if has_fission && p.ires > 0
            fis = Vector{Float64}(undef, ires * nsigz_eff)
            for it in 1:ires, iz in 1:nsigz_eff
                src_idx = iz + nsigz_eff * ((it - 1) + ires * (jg - 1))
                dst_idx = iz + nsigz_eff * (it - 1)
                fis[dst_idx] = snux[src_idx]
            end
            push!(fis_tables, WIMSResonanceTable(p.rdfid, ires, nsigz_eff,
                                                  xs_data.tempr[1:ires],
                                                  copy(sigz),
                                                  fis))
        end
    end

    return (abs_tables, fis_tables)
end
