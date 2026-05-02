# wimsr xsecs port — extract per-group cross sections from a GENDF tape.
#
# Ref: njoy-reference/src/wimsr.f90:863-1422 (subroutine xsecs).
#
# Walks the GENDF temperature/MAT/MF/MT structure produced by groupr and
# populates the per-group arrays the WIMS writer consumes:
#
#   spot[g]   potential scattering xs (init = sigp; if sigp==0, spot[g]=scat[g])
#   abs1[g]   MT=102 capture xs (Fortran wimsr.f90:1051-1053)
#   abs2[g]   other absorption MT=103..150 (Fortran line 1055)
#   sf0[g]    fission xs, inf-dilution (Fortran lines 1066-1071)
#   sfi[g]    fission xs, sigma_zero-iz indexed
#   sn2n[g]   (n,2n)+(n,3n)*2+(n,4n) (lines 1089-1097)
#   scat[g]   accumulated MF=6 transfer-out (lines 1177, 1233)
#   xi[g]     log slowing-down power (MT=252; default = 1+log(α)·α/(1-α))
#   snu[g]    nu-bar (MF=3/MT=452; line 1124)
#   xtr[g]    ← scat[g] + ab0[g] - csp1[g]               (line 1342)
#   ab0[g]    ← sf0[g] + abs1[g] + abs2[g] - sn2n[g]     (line 1329)
#   sdp[g]    ← scat[g] * xi[g] / log(egb[g]/egb[g+1])   (line 1359)
#   xs[g, g'] ngnd×ngnd nonthermal P0 transfer matrix accumulated from MF=6
#             (lines 1166-1191, 1218-1245). Sparse-emitted via l1[g]/l2[g] trim.
#   l1[g], l2[g]  per-row first/last nonzero column indices in xs.
#   nth       first thermal group (set when MF=3/MT=mti shows nonzero data)
#   tempr[t]  temperature list (one per temperature on tape)
#   egb[1..ngnd+1]  group boundaries from MF=1/MT=451
#
# Phase 58b coverage:
#   - MF=3 dispatch: MT=1 (skipped — p1flx default), MT=18-21,38 (fission),
#     MT=102-150 (absorption), MT=16/17/24/25/875-891 (n2n), MT=mti (nth),
#     MT=252 (xi), MT=452 (nu-bar)
#   - MF=6 dispatch: MT=2 (elastic), MT=mti, MT=18-21,38 (fission spectrum),
#     other MTs (inelastic-like) — temperature-independent path only
#   - First-temperature only (multi-temp thermal section is Phase 58c)
#   - Inf-dilution sigma_zero only (sgref >= 1e10 ⇒ isg=0)
#
# Parts marked TODO will diff against referenceTape27 until ported.

"""
    WimsrXSecs

Per-group cross sections extracted from a GENDF tape, in the layout
required by the WIMS writer. All vectors are length `ngnd` unless noted.
"""
struct WimsrXSecs
    ngnd::Int
    nfg::Int
    nrg::Int
    nth::Int                       # first thermal group (1-based)
    awr::Float64
    za::Float64
    # Temperature-independent (Fortran nscr2)
    spot::Vector{Float64}
    abs1::Vector{Float64}
    abs2::Vector{Float64}
    sfi::Vector{Float64}
    sf0::Vector{Float64}
    sn2n::Vector{Float64}
    scat::Vector{Float64}
    xi::Vector{Float64}
    snu::Vector{Float64}
    xtr::Vector{Float64}
    ab0::Vector{Float64}
    sdp::Vector{Float64}
    xs::Matrix{Float64}            # ngnd × ngnd nonthermal scatter matrix
    l1::Vector{Int}                # per-row first nonzero col (init ngnd)
    l2::Vector{Int}                # per-row last nonzero col (init 1)
    # Temperature-dependent (Fortran nscr3) — outer index = temperature
    xtr_t::Vector{Vector{Float64}}    # length ngnd per temp (only thermal slice used)
    ab0_t::Vector{Vector{Float64}}
    snus_t::Vector{Vector{Float64}}
    sf0_t::Vector{Vector{Float64}}
    xs_t::Vector{Matrix{Float64}}     # ngnd×ngnd thermal scatter matrix per temp
    l1_t::Vector{Vector{Int}}
    l2_t::Vector{Vector{Int}}
    egb::Vector{Float64}
    tempr::Vector{Float64}
    has_fission::Bool
end

"""
    wimsr_extract_xsecs(in_path, p::WimsrParams) -> WimsrXSecs

Walk the GENDF tape and populate the WIMS-output buffers per the routing
dispatch in `xsecs` (wimsr.f90:863-1422).
"""
function wimsr_extract_xsecs(in_path::AbstractString, p::WimsrParams)
    ngnd = p.ngnd
    nfg  = p.nfg
    nrg  = p.nrg
    nnt  = nfg + nrg

    awr, za, fission_present = _wimsr_scan_header(in_path, p.mat)
    egb, tempr = _wimsr_read_gendf_metadata(in_path, p.mat)
    n_temps = max(length(tempr), 1)

    # --- Temp-independent persistent state (only updated on temp 1)
    spot = fill(p.sigp, ngnd)
    abs1 = zeros(ngnd); abs2 = zeros(ngnd)
    sfi  = zeros(ngnd); sf0  = zeros(ngnd)
    sn2n = zeros(ngnd); snu  = zeros(ngnd)

    α = ((awr - 1) / (awr + 1))^2
    xxi = 1 + log(α) * α / (1 - α)
    xi = fill(xxi, ngnd)

    # Temp-1 nonthermal scatter matrix (frozen at end of temp 1)
    xs1 = zeros(ngnd, ngnd)
    l1_save = fill(ngnd, ngnd)
    l2_save = ones(Int, ngnd)
    scat1 = zeros(ngnd)
    csp1_save = zeros(ngnd)

    # Per-temperature thermal accumulators
    xtr_t   = [zeros(ngnd)        for _ in 1:n_temps]
    ab0_t   = [zeros(ngnd)        for _ in 1:n_temps]
    snus_t  = [zeros(ngnd)        for _ in 1:n_temps]
    sf0_t_v = [zeros(ngnd)        for _ in 1:n_temps]
    xs_t    = [zeros(ngnd, ngnd)  for _ in 1:n_temps]
    l1_t    = [fill(ngnd, ngnd)   for _ in 1:n_temps]
    l2_t    = [ones(Int, ngnd)    for _ in 1:n_temps]

    isg = (p.sgref < 1e10) ? 1 : 0
    nth1_state = Ref(0)
    has_fission = Ref(fission_present)

    _wimsr_dispatch_multitemp!(in_path, p, isg,
        abs1, abs2, sfi, sf0, sn2n, snu, xi,
        spot, scat1, xs1, l1_save, l2_save, csp1_save,
        xtr_t, ab0_t, snus_t, sf0_t_v, xs_t, l1_t, l2_t,
        nth1_state, has_fission, n_temps)

    nth = nth1_state[] > 0 ? (ngnd - nth1_state[] + 1) : (ngnd - nnt)

    # Finalise temp-indep (using temp-1 csp1 + scat)
    ab0 = [sf0[i] + abs1[i] + abs2[i] - sn2n[i] for i in 1:ngnd]
    xtr = [scat1[i] + ab0[i] - csp1_save[i] for i in 1:ngnd]
    sdp = zeros(ngnd)
    for i in 1:ngnd
        if i + 1 <= length(egb) && egb[i] > 0 && egb[i+1] > 0 && egb[i] != egb[i+1]
            sdp[i] = scat1[i] * xi[i] / log(egb[i] / egb[i+1])
        end
    end
    if p.sigp == 0
        for i in 1:ngnd
            spot[i] = scat1[i]
        end
    end

    return WimsrXSecs(ngnd, nfg, nrg, nth, awr, za,
                     spot, abs1, abs2, sfi, sf0, sn2n, scat1, xi, snu,
                     xtr, ab0, sdp, xs1, l1_save, l2_save,
                     xtr_t, ab0_t, snus_t, sf0_t_v, xs_t, l1_t, l2_t,
                     egb, tempr, has_fission[])
end

# ============================================================================
# MF=1/MT=451 minimal scanner — AWR, ZA, fission-present flag
# ============================================================================

function _wimsr_scan_header(in_path::AbstractString, target_mat::Int)
    awr = 0.0; za = 0.0; fission_present = false
    open(in_path, "r") do io
        for line in eachline(io)
            length(line) < 75 && continue
            row = rpad(line, 80)
            mv = tryparse(Int, strip(row[67:70])); mv === nothing && continue
            mv != target_mat && continue
            mf = tryparse(Int, strip(row[71:72]))
            mt = tryparse(Int, strip(row[73:75]))
            seq = tryparse(Int, strip(row[76:80]))
            (mf === nothing || mt === nothing || seq === nothing) && continue
            if mf == 1 && mt == 451 && seq == 1 && awr == 0.0
                za = parse_endf_float(row[1:11])
                awr = parse_endf_float(row[12:22])
            end
            if mf == 3 && mt == 18 && seq == 1
                fission_present = true
            end
        end
    end
    return (awr, za, fission_present)
end

# ============================================================================
# MF=1/MT=451 group-bounds + temperature reader
# ============================================================================

function _wimsr_read_gendf_metadata(in_path::AbstractString, target_mat::Int)
    egb = Float64[]
    tempr = Float64[]
    seen_temps = Set{Int}()

    lines = readlines(in_path)
    idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        row = rpad(lines[idx], 80)
        mv  = tryparse(Int, strip(row[67:70])); mv === nothing && (idx += 1; continue)
        mf  = tryparse(Int, strip(row[71:72]))
        mt  = tryparse(Int, strip(row[73:75]))
        seq = tryparse(Int, strip(row[76:80]))
        if mv == target_mat && mf == 1 && mt == 451 && seq == 1
            # GENDF MF=1/MT=451 HEAD has L2=NZ at cols 34-44, N2=NTW at 56-66.
            # NGN lives on seq=2's L1 field (cols 23-33). Walk to seq=2 to get it.
            idx += 1
            ngn = 0
            if idx <= length(lines)
                row2 = rpad(lines[idx], 80)
                t = parse_endf_float(row2[1:11])
                rounded = round(Int, t)
                if !(rounded in seen_temps)
                    push!(seen_temps, rounded); push!(tempr, t)
                end
                ngn = something(tryparse(Int, strip(row2[23:33])), 0)
                idx += 1
            end
            body = Float64[]
            while idx <= length(lines)
                rrow = rpad(lines[idx], 80)
                mt_check = tryparse(Int, strip(rrow[73:75]))
                mt_check == 451 || break
                seq_check = tryparse(Int, strip(rrow[76:80]))
                seq_check == 99999 && break
                seq_check == 0 && break
                for col in 0:5
                    s = rrow[1+11*col : 11+11*col]
                    all(isspace, s) && break
                    push!(body, parse_endf_float(s))
                end
                idx += 1
            end
            if length(body) >= ngn + 1 && isempty(egb)
                # GENDF stores bounds ascending; wimsr stores descending
                # (Fortran wimsr.f90:346 `egb(i)=scr(ngnd1-i+i1)`).
                # Find the first monotonic increasing run of length ngn+1
                # and reverse it.
                for start in 1:length(body)-ngn
                    run = body[start:start+ngn]
                    if all(run[i] <= run[i+1] for i in 1:ngn) && run[end] > run[1]
                        egb = reverse(run); break
                    end
                end
            end
        else
            idx += 1
        end
    end
    return (egb, tempr)
end

# ============================================================================
# Tape dispatch — walk records and route per MT
# ============================================================================

# Position of body element (i, iz, ig2) in the LIST body (1-based body[]).
# Fortran convention: scr[lz + i + nl*(iz-1) + nl*nz*(ig2-1)] where lz=6 is
# the LIST CONT-record stride before body. In our split (header parsed
# separately, body[] is body-only), the index is just:
#   body[i + nl*(iz-1) + nl*nz*(ig2-1)]
@inline _bidx(i, iz, ig2, nl, nz) = i + nl*(iz - 1) + nl*nz*(ig2 - 1)

"""
Multi-temperature dispatcher: walks the GENDF tape once per temperature,
routing records to the temp-1 buffers (for the temp-indep nscr2 layout)
on the first pass, and to per-temperature thermal accumulators on every
pass. Mirrors Fortran's outer temp-loop at wimsr.f90:947 with per-temp
reset of `scat`, `xs`, `l1`, `l2`, `csp1`.
"""
function _wimsr_dispatch_multitemp!(in_path, p::WimsrParams, isg,
        abs1_out, abs2_out, sfi_out, sf0_out, sn2n_out, snu_out, xi_out,
        spot, scat1, xs1, l1_save, l2_save, csp1_save,
        xtr_t, ab0_t, snus_t, sf0_t_v, xs_t, l1_t, l2_t,
        nth1_state, has_fission, n_temps::Int)
    ngnd = p.ngnd

    lines = readlines(in_path)

    α = ((spot[1] != p.sigp ? 1.0 : 1.0) - 1)  # noop, kept for clarity
    # Fortran resets all per-group accumulators each temp (line 948).
    # So we allocate fresh arrays per pass; only outputs at temp==1
    # populate the temp-indep buffers.
    for ti in 1:n_temps
        abs1 = zeros(ngnd); abs2 = zeros(ngnd)
        sfi  = zeros(ngnd); sf0  = zeros(ngnd)
        sn2n = zeros(ngnd); snu  = zeros(ngnd); snus = zeros(ngnd)
        scat = zeros(ngnd); csp1 = zeros(ngnd)
        # xi is temp-indep but per-loop default; MT=252 (if present) overwrites
        xi = copy(xi_out)
        xs   = zeros(ngnd, ngnd)
        l1   = fill(ngnd, ngnd)
        l2   = ones(Int, ngnd)

        _wimsr_dispatch_one_temp!(lines, p, ti, isg,
            abs1, abs2, sfi, sf0, sn2n, snu, snus, xi,
            scat, xs, l1, l2, csp1,
            nth1_state, has_fission)

        # Per-temp finalisation (Fortran wimsr.f90:1339-1343):
        #   if (ip1opt > 0): xs[i, i] -= csp1[i]   ← diagonal P1 correction
        #   xtr[i] = scat[i] + ab0[i] - csp1[i]
        # The diagonal correction is what makes the WIMS scatter matrix
        # transport-corrected (P0 + P1-self-removal).
        if p.ip1opt > 0
            for i in 1:ngnd
                xs[i, i] -= csp1[i]
            end
        end
        ab0_full = [sf0[i] + abs1[i] + abs2[i] - sn2n[i] for i in 1:ngnd]
        xtr_full = [scat[i] + ab0_full[i] - csp1[i] for i in 1:ngnd]

        for i in 1:ngnd
            xtr_t[ti][i]   = xtr_full[i]
            ab0_t[ti][i]   = ab0_full[i]
            snus_t[ti][i]  = snu[i] * sf0[i]
            sf0_t_v[ti][i] = sf0[i]
        end
        xs_t[ti] .= xs
        l1_t[ti] .= l1
        l2_t[ti] .= l2

        # Save temp-1 state for the temp-indep nscr2 block
        if ti == 1
            scat1     .= scat
            xs1       .= xs
            l1_save   .= l1
            l2_save   .= l2
            csp1_save .= csp1
            abs1_out  .= abs1
            abs2_out  .= abs2
            sfi_out   .= sfi
            sf0_out   .= sf0
            sn2n_out  .= sn2n
            snu_out   .= snu
            xi_out    .= xi
        end
    end
end

function _wimsr_dispatch_one_temp!(lines, p::WimsrParams, target_temp_idx::Int, isg::Int,
                                     abs1, abs2, sfi, sf0, sn2n, snu, snus, xi,
                                     scat, xs, l1, l2, csp1,
                                     nth1_state, has_fission)
    ngnd = p.ngnd
    nfg  = p.nfg
    nrg  = p.nrg
    nnt  = nfg + nrg
    mti  = p.mti
    mtc  = p.mtc

    idx = 1
    cur_temp_idx = 0
    section_nl = 1
    section_nz = 1
    section_mt = 0
    section_mf = 0
    in_target_temp = false

    # Locate iz_target if specific sigma_zero requested (sgref < 1e10).
    # For T11 with sgref=1e10, isg=0 ⇒ first sig0 column suffices.
    iz_target = 1

    i318 = false  # MT=18 seen flag (Fortran line 1061)
    jn2n = 0      # n2n MT seen indicator (line 1036)

    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        row = rpad(lines[idx], 80)
        mv = tryparse(Int, strip(row[67:70]))
        mf = tryparse(Int, strip(row[71:72]))
        mt = tryparse(Int, strip(row[73:75]))
        seq = tryparse(Int, strip(row[76:80]))
        if mv === nothing || mf === nothing || mt === nothing || seq === nothing
            idx += 1; continue
        end

        # Track temperature blocks via MF=1/MT=451 HEAD records
        if mv == p.mat && mf == 1 && mt == 451 && seq == 1
            cur_temp_idx += 1
            in_target_temp = (cur_temp_idx == target_temp_idx)
            # Stop early when we move past the target temp
            if cur_temp_idx > target_temp_idx
                break
            end
        end

        # Section HEAD (seq==1) for MF=3/6 with mt>0: capture NL, NZ
        if mv == p.mat && mf in (3, 6) && mt > 0 && seq == 1 && in_target_temp
            section_nl = something(tryparse(Int, strip(row[23:33])), 1)
            section_nz = something(tryparse(Int, strip(row[34:44])), 1)
            section_mt = mt
            section_mf = mf
            if mt == 18; i318 = true; has_fission[] = true; end
            if mt == 16; jn2n = 16; end
            idx += 1
            continue
        end

        # Per-group LIST record header (mid-section)
        if mv == p.mat && mf in (3, 6) && mt > 0 && seq >= 2 && in_target_temp &&
           mt == section_mt && mf == section_mf
            ng2_in   = tryparse(Int, strip(row[23:33]))
            ig2lo    = tryparse(Int, strip(row[34:44]))
            nw       = tryparse(Int, strip(row[45:55]))
            ig       = tryparse(Int, strip(row[56:66]))
            if ng2_in === nothing || ig2lo === nothing || nw === nothing || ig === nothing
                idx += 1; continue
            end
            if !(1 <= ig <= ngnd) || nw <= 0
                idx += 1; continue
            end
            body, idx_after = _wimsr_read_list_body(lines, idx + 1, mt, nw)
            if length(body) < nw
                idx = idx_after; continue
            end

            jg = ngnd - ig + 1   # Fortran flips to top-down
            nl = section_nl; nz = section_nz

            if mf == 3
                _dispatch_mf3!(mt, mti, mtc, jg, ig, nl, nz, ng2_in, ig2lo,
                                body, isg, iz_target, abs1, abs2, sfi, sf0,
                                sn2n, snu, xi, jn2n, i318, nth1_state)
            elseif mf == 6
                _dispatch_mf6!(mt, mti, mtc, jg, ig, nl, nz, ng2_in, ig2lo,
                                body, isg, iz_target, scat, xs, l1, l2,
                                csp1, ngnd, nfg, nrg, nth1_state)
            end
            idx = idx_after
            continue
        end

        idx += 1
    end
    return nothing
end

# Read NW floats from subsequent body lines until the section MT changes.
function _wimsr_read_list_body(lines, start_idx::Int, expect_mt::Int, nw::Int)
    body = Float64[]
    idx = start_idx
    while idx <= length(lines) && length(body) < nw
        brow = rpad(lines[idx], 80)
        mt_check = tryparse(Int, strip(brow[73:75]))
        mt_check == expect_mt || break
        for col in 0:5
            length(body) >= nw && break
            s = brow[1+11*col : 11+11*col]
            all(isspace, s) && break
            push!(body, parse_endf_float(s))
        end
        idx += 1
    end
    return (body, idx)
end

# MF=3 dispatch (Fortran wimsr.f90:1029-1114)
function _dispatch_mf3!(mt, mti, mtc, jg, ig, nl, nz, ng2_in, ig2lo,
                          body, isg, iz_target, abs1, abs2, sfi, sf0,
                          sn2n, snu, xi, jn2n, i318, nth1_state)
    # For MF=3 with ig2_in=2, position = body[1 + nl*nz*(2-1)] = body[1 + nl*nz]
    # Specific sig0: body[1 + nl*(iz-1) + nl*nz*1]
    pos_infdil = _bidx(1, 1, 2, nl, nz)            # (i=1, iz=1, ig2=2)
    pos_sigsel = _bidx(1, iz_target, 2, nl, nz)    # (i=1, iz=iz_target, ig2=2)
    xs_infdil  = body[pos_infdil]
    xs_iz      = (isg > 0 && nz >= iz_target) ? body[pos_sigsel] : xs_infdil

    if 102 <= mt <= 150
        if mt == 102
            abs1[jg] = xs_iz
        else
            abs2[jg] += xs_infdil
        end
    elseif (18 <= mt <= 21) || mt == 38
        if mt == 18 || i318
            sfi[jg] += xs_infdil
            sf0[jg] += xs_iz
        end
    elseif mt == 16 || mt == 24 || (875 <= mt <= 891 && jn2n != 16)
        sn2n[jg] += xs_infdil
    elseif mt == 17 || mt == 25
        sn2n[jg] += 2 * xs_infdil
    elseif mt == 252
        xi[jg] = xs_infdil
    elseif mt == mti
        # Thermal inelastic boundary update (Fortran line 1102)
        if ig > nth1_state[] && xs_infdil != 0.0
            nth1_state[] = ig
        end
    elseif mt == 452
        # MF=3/MT=452: nu-bar (Fortran line 1124, loca = l+lz+1)
        # That's body position 2 (1-based) in our split
        snu[jg] = body[2]
    end
    return nothing
end

# MF=6 dispatch (Fortran wimsr.f90:1160-1250) — temperature-independent path
# (Fortran `300 continue` branch at line 1205, gated by jtemp==1 at line 1206).
function _dispatch_mf6!(mt, mti, mtc, jg, ig, nl, nz, ng2_in, ig2lo,
                          body, isg, iz_target, scat, xs, l1, l2,
                          csp1, ngnd, nfg, nrg, nth_state)
    # Fortran routing (wimsr.f90:1043-1046, 1164-1166, 1214-1217):
    #   MT in {2, mti, mtc} → 270 block (temp-dep): we're called from per-temp
    #     dispatch already, so accumulate into the passed scat/xs.
    #     Inside 270: MT=2 with jg<nth falls back to 301 (temp-indep style);
    #     MT=2 with jg>=nth is skipped (thermal elastic is handled by mti).
    #   MT in {18-21, 38} → 315 (fission spectrum, NOT scat). Skip here.
    #   MT in {221-250} \ {mti} → skip entirely (line 1217).
    if (18 <= mt <= 21) || mt == 38
        return nothing  # TODO: fission-spectrum accumulation (Phase 58c)
    end
    if 221 <= mt <= 250 && mt != mti
        return nothing
    end
    # Fortran line 1164: MT=2 thermal-source skip (jg >= nth → goto 295).
    # Thermal elastic is handled by MT=mti accumulation; double-counting here
    # makes scat[thermal] ~2× too big.
    # nth_state holds nth1 (highest GENDF ig where MT=mti has data); nth in
    # WIMS top-down = ngnd - nth1 + 1. Fall back to nominal first-thermal
    # boundary (ngnd - nfg - nrg + 1 in jg = nfg+nrg+1) if nth1 unset.
    nth1 = nth_state[]
    nth = nth1 > 0 ? (ngnd - nth1 + 1) : (nfg + nrg + 1)
    if mt == 2 && jg >= nth
        return nothing
    end

    il = 1   # Legendre order index for P0 scattering (Fortran's `il` is 1)
    max_g = 0; min_g = ngnd
    for i in 2:ng2_in  # secondary group loop (Fortran line 1220)
        ig2 = ig2lo + i - 2
        if 1 <= ig2 <= ngnd
            jg2 = ngnd - ig2 + 1
            # Position: (il, nz_sel, i)  (Fortran lines 1223-1229)
            # For MT=2 (elastic) inf-dilution branch: use last sig0 (iz=nz)
            # otherwise iz=1 default per (il-1) + nl*nz*(i-1).
            if mt == 2
                if isg > 0 && nz >= iz_target
                    pos = _bidx(il, iz_target, i, nl, nz)
                else
                    pos = _bidx(il, nz, i, nl, nz)  # (i=il, iz=nz, ig2=i)
                end
            else
                pos = _bidx(il, 1, i, nl, nz)
            end
            if pos < 1 || pos > length(body)
                continue
            end
            v = body[pos]
            xs[jg, jg2] += v
            scat[jg]    += v
            if nl != 1
                # csp1 = P1-current transport correction (Fortran line 1241):
                #   csp1[jg2] += body[pos+1] * p1flx[jg]/p1flx[jg2]
                # Without populated p1flx (Phase 58c — MT=1 dispatch), use
                # ratio=1. This produces xtr correct to ~0.1% rather than bit;
                # the residual is the missing p1flx weighting.
                if jg2 < ngnd - nfg - nrg + 1   # jg2 < nth
                    if pos + 1 <= length(body)
                        csp1[jg2] += body[pos + 1]
                    end
                end
            end
            # Fortran 270 block (MT in {2, mti, mtc}) gates l1/l2 on value!=0
            # (line 1184-1186). 300 block (other MTs) extends unconditionally
            # (line 1234-1235). Effect: MT=mti zero transfers don't extend
            # l2 into thermal sinks where they'd produce trailing zeros.
            in_270_block = (mt == 2 || mt == mti || mt == mtc)
            if !in_270_block || v != 0.0
                if jg2 > max_g; max_g = jg2; end
                if jg2 < min_g; min_g = jg2; end
            end
        end
    end
    if min_g < l1[jg]; l1[jg] = min_g; end
    if max_g > l2[jg]; l2[jg] = max_g; end
    return nothing
end
