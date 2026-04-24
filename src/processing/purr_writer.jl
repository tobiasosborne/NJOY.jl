# PURR PENDF writer -- emits MT152 (Bondarenko self-shielded XS) and
# MT153 (probability tables) into the URR region of a PENDF tape.
#
# Refs: njoy-reference/src/purr.f90
#   MT152 emission : lines 415-470
#   MT153 emission : lines 472-576
#   DICT rewrite   : lines 316-400
#   PENDF splicing : around the material/temperature loop
#
# Architecture: line-by-line state-machine splice of the input PENDF,
# inserting two new MF2 sections after MT151's SEND. DICT (MF1/MT451) is
# rewritten to include both MT152 and MT153 entries with Fortran's exact
# NC counts. Pattern mirrors src/orchestration/modules/unresr.jl's
# _write_unresr_pendf but inserts two MTs instead of one.

using Printf

"""
    PurrBlock

Per-temperature container for purr output that will be serialised into
one (MT152, MT153) pair in the PENDF tape.

- `sigu[nunr, 5*nsigz]` — Bondarenko-weighted self-shielded cross
  sections, packed per Fortran purr.f90:448-452 as
  `sigu[ie, (ix-1)*nsigz + isigz]` with `ix in 1..5` = (total, elastic,
  fission, capture, transport/MT4).
- `tabl[nbin, 5, nunr]` — probability-table cross sections per bin and
  reaction, already `sigfig(_,7,0)`-rounded.
- `heating[nbin, nunr]` — conditional heating probabilities (zero when
  no heatr precedes purr; matches `ihave=0` branch in Fortran).
"""
struct PurrBlock
    temp::Float64
    sigu::Matrix{Float64}
    tabl::Array{Float64,3}
    heating::Matrix{Float64}
end

# Fortran purr.f90 DICT NC formulas (lines 365-368).
_purr_mt152_nc(nsigz::Int, nunr::Int) = 2 + div(nsigz + nunr * (1 + 5 * nsigz) - 1, 6)
_purr_mt153_nc(nbin::Int, nunr::Int)  = 2 + div((1 + 6 * nbin) * nunr - 1, 6)

# Pack a Vector{Float64} 6-per-line with MAT/MF/MT/SEQ trailers and an
# appended SEND record. Mirrors Fortran listio/moreio/asend.
function _purr_emit_list(io::IO, data::Vector{Float64},
                         mat::Int, mf::Int, mt::Int, seq::Int)
    idx = 1
    while idx <= length(data)
        buf = ""
        for _ in 1:6
            if idx <= length(data)
                buf *= format_endf_float(data[idx]; extended=false)
                idx += 1
            end
        end
        seq += 1
        println(io, rpad(buf, 66) * @sprintf("%4d%2d%3d%5d", mat, mf, mt, seq))
    end
    println(io, rpad("", 66) * @sprintf("%4d%2d%3d%5d", mat, mf, 0, 99999))
end

"""
MT152 block — Bondarenko self-shielded XS.
Ref: njoy-reference/src/purr.f90:415-470
"""
function _write_mt152_purr(io::IO, mat::Int, urr, sigz::Vector{Float64},
                           block::PurrBlock)
    nunr = length(urr.energies)
    nsigz = length(sigz)
    trailer(seq) = @sprintf("%4d%2d%3d%5d", mat, 2, 152, seq)

    # HEAD CONT — ZA, AWR, LSSF, 0, 0, INTUNR
    seq = 1
    cont = format_endf_float(urr.za; extended=false) *
           format_endf_float(urr.awr; extended=false) *
           lpad(string(urr.lssf), 11) * lpad("0", 11) *
           lpad("0", 11) * lpad(string(urr.intunr), 11)
    println(io, cont * trailer(seq))

    # LIST header — TEMZ, 0.0, NREACT=5, NSIGZ, NPL, NE
    npl = nsigz + nunr * (1 + 5 * nsigz)
    seq += 1
    list_head = format_endf_float(block.temp; extended=false) *
                format_endf_float(0.0; extended=false) *
                lpad("5", 11) * lpad(string(nsigz), 11) *
                lpad(string(npl), 11) * lpad(string(nunr), 11)
    println(io, list_head * trailer(seq))

    # LIST data: [sigz values] then per-energy [E, sigu(1..5, 1..nsigz)]
    # Per Fortran 434-453: sigz first (as raw), then for each energy
    # a raw E (line 447 — no sigfig) and 5*nsigz sigma values.
    data = Float64[]
    append!(data, sigz)
    for ie in 1:nunr
        push!(data, urr.energies[ie])
        for k in 1:5*nsigz
            push!(data, block.sigu[ie, k])
        end
    end
    _purr_emit_list(io, data, mat, 2, 152, seq)
end

"""
MT153 block — probability tables + heating conditionals.
Ref: njoy-reference/src/purr.f90:472-576
"""
function _write_mt153_purr(io::IO, mat::Int, urr, nbin::Int,
                           iinel::Int, iabso::Int, block::PurrBlock)
    nunr = length(urr.energies)
    trailer(seq) = @sprintf("%4d%2d%3d%5d", mat, 2, 153, seq)

    # HEAD CONT — ZA, AWR, IINEL, IABSO, INTUNR, NBIN
    seq = 1
    cont = format_endf_float(urr.za; extended=false) *
           format_endf_float(urr.awr; extended=false) *
           lpad(string(iinel), 11) * lpad(string(iabso), 11) *
           lpad(string(urr.intunr), 11) * lpad(string(nbin), 11)
    println(io, cont * trailer(seq))

    # LIST header — TEMZ, 0.0, LSSF, 0, NPL, NE
    npl = (1 + 6 * nbin) * nunr
    seq += 1
    list_head = format_endf_float(block.temp; extended=false) *
                format_endf_float(0.0; extended=false) *
                lpad(string(urr.lssf), 11) * lpad("0", 11) *
                lpad(string(npl), 11) * lpad(string(nunr), 11)
    println(io, list_head * trailer(seq))

    # LIST data: per energy -> [E, tabl(1..nbin, 1..5), heating(1..nbin)].
    # Layout per Fortran 499-524: E, then 5 reaction blocks of nbin
    # values each, then nbin heating conditionals. Total 1+6*nbin per E.
    data = Float64[]
    for ie in 1:nunr
        push!(data, urr.energies[ie])
        for ix in 1:5
            for j in 1:nbin
                push!(data, block.tabl[j, ix, ie])
            end
        end
        for j in 1:nbin
            push!(data, block.heating[j, ie])
        end
    end
    _purr_emit_list(io, data, mat, 2, 153, seq)
end

"""
    _write_purr_pendf(pendf_in, pendf_out, params, urr, blocks; iinel, iabso)

Splice the purr-generated MT152 and MT153 sections into the input
PENDF tape, rewriting the MF1/MT451 directory and MT451 self-reference
NC count per Fortran purr.f90:316-413.

State machine:
  (a) At each MF1/MT451 HEAD, rewrite HEAD+CONT NXC, then copy NWD
      hollerith lines, then rewrite the directory inserting MT152 and
      MT153 entries just after MT151. Renumber the MF1 SEND.
  (b) At each MF2/MT151 SEND (detected by the next line being the MF2
      FEND), append our MT152 and MT153 sections before continuing.
  (c) Everything else is copied verbatim.
"""
function _write_purr_pendf(pendf_in::String, pendf_out::String,
                           params::PurrParams, urr,
                           blocks::Vector{PurrBlock};
                           iinel::Int=-1, iabso::Int=-1)
    in_lines = readlines(pendf_in)
    mat = params.mat

    function _lp(line)
        nl = length(line)
        nl < 14 && return (0, 0, 0)
        (_parse_int(line[nl-13:nl-10]), _parse_int(line[nl-9:nl-8]), _parse_int(line[nl-7:nl-5]))
    end

    function _dir_vals(line)
        p = rpad(line, 80)[1:66]
        vals = Int[]
        for f in 1:6
            v = Int(_parse_int(p[(f-1)*11+1:f*11]))
            v != 0 && push!(vals, v)
        end
        vals
    end

    mt152_nc = _purr_mt152_nc(length(params.sigz), length(urr.energies))
    mt153_nc = _purr_mt153_nc(params.nbin, length(urr.energies))

    block_idx = 0
    expect_new_block = true

    open(pendf_out, "w") do out
        i = 1
        while i <= length(in_lines)
            line = in_lines[i]
            pm, pf, pt = _lp(line)

            # MF1/MT451 HEAD -> rewrite directory to include MT152/MT153
            if pm == mat && pf == 1 && pt == 451 && expect_new_block
                expect_new_block = false
                block_idx += 1

                cont_line = in_lines[i + 1]
                cont_p = rpad(cont_line, 80)[1:80]
                nwd_blk = Int(_parse_int(cont_p[45:55]))
                nxc_blk = Int(_parse_int(cont_p[56:66]))

                dir_start = i + 2 + nwd_blk
                dir_end = dir_start + nxc_blk - 1
                blk_has_152 = false
                blk_has_153 = false
                mt151_dir_idx = 0
                for j in dir_start:min(dir_end, length(in_lines))
                    vals = _dir_vals(in_lines[j])
                    length(vals) >= 2 || continue
                    vals[1] == 2 && vals[2] == 152 && (blk_has_152 = true)
                    vals[1] == 2 && vals[2] == 153 && (blk_has_153 = true)
                    vals[1] == 2 && vals[2] == 151 && (mt151_dir_idx = j)
                end
                # Canonical rule: always emit MT152+MT153 fresh after MT151.
                # Upstream entries (reconr sometimes seeds a stub MT152) get
                # skipped here and their LIST sections replaced in phase (b).
                new_count = 2 - (blk_has_152 ? 1 : 0) - (blk_has_153 ? 1 : 0)

                # Rewrite HEAD line: NXC field (cols 56:66) += new_count
                if new_count > 0
                    hp = rpad(line, 80)[1:80]
                    old_nxc_head = Int(_parse_int(hp[56:66]))
                    println(out, hp[1:55] * lpad(string(old_nxc_head + new_count), 11) * hp[67:80])
                else
                    println(out, line)
                end

                # Rewrite CONT line: NXC += new_count (cols 56:66)
                if new_count > 0
                    seq_str = rpad(cont_line, 80)[67:80]
                    println(out, cont_p[1:44] * lpad(string(nwd_blk), 11) *
                                  lpad(string(nxc_blk + new_count), 11) * seq_str)
                else
                    println(out, cont_line)
                end

                # Copy description (hollerith) lines
                for j in (i + 2):(i + 1 + nwd_blk)
                    println(out, in_lines[j])
                end

                # Rewrite directory: copy non-152/153 entries, then insert fresh
                # 152+153 pair after MT151. 1/451 self-ref NC is bumped so the
                # downstream count stays consistent.
                seq_offset = 0
                for j in dir_start:dir_end
                    dln = in_lines[j]
                    vals = _dir_vals(dln)

                    # Skip upstream MT152/MT153 DICT entries (we insert fresh)
                    if length(vals) >= 2 && vals[1] == 2 && (vals[2] == 152 || vals[2] == 153)
                        seq_offset -= 1  # this line disappears -> bump subsequent seqs down
                        continue
                    end

                    if new_count > 0 && length(vals) >= 3 && vals[1] == 1 && vals[2] == 451
                        old_nc = vals[3]
                        trailer = rpad(dln, 80)[67:80]
                        println(out, @sprintf("%22s%11d%11d%11d%11d", "", 1, 451, old_nc + new_count, 0) * trailer)
                    elseif seq_offset != 0
                        old_seq = _parse_int(rpad(dln, 80)[76:80])
                        println(out, rpad(dln, 80)[1:75] * @sprintf("%5d", old_seq + seq_offset))
                    else
                        println(out, dln)
                    end

                    if j == mt151_dir_idx
                        base_seq = _parse_int(rpad(dln, 80)[76:80]) + seq_offset
                        seq_offset += 1
                        @printf(out, "%22s%11d%11d%11d%11d%4d%2d%3d%5d\n",
                                "", 2, 152, mt152_nc, 0, mat, 1, 451, base_seq + 1)
                        seq_offset += 1
                        @printf(out, "%22s%11d%11d%11d%11d%4d%2d%3d%5d\n",
                                "", 2, 153, mt153_nc, 0, mat, 1, 451, base_seq + 2)
                    end
                end

                # MF1 SEND: fixed sentinel seq=99999 (not renumbered)
                send_idx = dir_end + 1
                if send_idx <= length(in_lines)
                    println(out, in_lines[send_idx])
                end
                fend_idx = send_idx + 1
                fend_idx <= length(in_lines) && println(out, in_lines[fend_idx])
                i = fend_idx + 1
                continue
            end

            # MF2/MT151 SEND detected when next line is MF2 FEND (or material FEND)
            if pm == mat && pf == 2 && pt == 0 && block_idx > 0 && block_idx <= length(blocks)
                next_pm, next_pf, next_pt = i + 1 <= length(in_lines) ? _lp(in_lines[i+1]) : (0, 0, 0)
                is_fend_next = (next_pm == mat && next_pf == 0 && next_pt == 0) || (next_pm == 0 && next_pf == 0)
                is_mt152_next = (next_pm == mat && next_pf == 2 && next_pt == 152)
                is_mt153_next = (next_pm == mat && next_pf == 2 && next_pt == 153)

                if is_fend_next
                    println(out, line)
                    _write_mt152_purr(out, mat, urr, params.sigz, blocks[block_idx])
                    _write_mt153_purr(out, mat, urr, params.nbin, iinel, iabso, blocks[block_idx])
                    i += 1
                    continue
                elseif is_mt152_next || is_mt153_next
                    # Replace existing MT152/153 sections in-place
                    println(out, line)
                    _write_mt152_purr(out, mat, urr, params.sigz, blocks[block_idx])
                    _write_mt153_purr(out, mat, urr, params.nbin, iinel, iabso, blocks[block_idx])
                    i += 1
                    # Skip old 152/153 content + their SEND records
                    while i <= length(in_lines)
                        pm2, pf2, pt2 = _lp(in_lines[i])
                        if pm2 == mat && pf2 == 2 && (pt2 == 152 || pt2 == 153)
                            i += 1
                        elseif pm2 == mat && pf2 == 2 && pt2 == 0
                            next2 = i + 1 <= length(in_lines) ? _lp(in_lines[i+1]) : (0, 0, 0)
                            if next2[1] == mat && next2[2] == 2 && (next2[3] == 152 || next2[3] == 153)
                                i += 1
                            else
                                i += 1; break
                            end
                        else
                            break
                        end
                    end
                    continue
                end
            end

            # MEND -> next material block expected
            if pm == 0 && pf == 0 && pt == 0
                expect_new_block = true
            end

            println(out, line)
            i += 1
        end
    end
end
