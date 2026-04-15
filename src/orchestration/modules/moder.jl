# moder module runner -- Tape format conversion / copy
#
# In Fortran, MODER converts between binary and ASCII tapes.
# In Julia, all tapes are ASCII. MODER reduces to file copy/registration.
#
# Two roles in a typical chain:
#   1. Initial: `moder 20 -21` — register ENDF tape path for downstream modules
#   2. Final: `moder -23 25` — copy final PENDF to output tape

"""
    moder_module(tapes::TapeManager, mc::ModuleCall)

Run MODER: copy/convert tape. In Julia this is essentially a file copy
since we work exclusively in ASCII format.
"""
function moder_module(tapes::TapeManager, mc::ModuleCall)
    isempty(mc.raw_cards) && return nothing

    card1 = mc.raw_cards[1]
    length(card1) < 2 && return nothing

    nin  = _parse_int_token(card1[1])
    nout = _parse_int_token(card1[2])

    # Fortran moder extract/merge mode: if 1 <= |nin| <= 19, card 1's nin
    # is a flag (1=endf/pendf, 2=gendf, 3=errorr) — not a real tape unit.
    # Card 2 is a tape-ID label; card 3+ are (real_nin, matd) pairs
    # terminated by real_nin=0. Single-material stub: cp the first
    # referenced tape to the output. Material filtering / multi-material
    # merging is not yet implemented.
    if 1 <= abs(nin) <= 19
        return _moder_extract_stub!(tapes, mc.raw_cards, nin, nout)
    end

    nin_path = resolve(tapes, nin)
    nout_path = resolve(tapes, nout)

    @info "moder: tape $(abs(nin)) → tape $(abs(nout))"

    if isfile(nin_path)
        if nin_path != nout_path
            cp(nin_path, nout_path; force=true)
        end
        register!(tapes, nout, nout_path)
    else
        @warn "moder: input tape $nin_path not found"
    end

    nothing
end

function _moder_extract_stub!(tapes::TapeManager, cards, nin_flag::Int, nout::Int)
    if length(cards) < 3 || isempty(cards[3])
        @warn "moder: extract-mode missing card 3 (nin_real, matd) — skipping"
        return nothing
    end
    real_nin = _parse_int_token(cards[3][1])
    if real_nin == 0
        @warn "moder: extract-mode card 3 starts with 0 sentinel — nothing to extract"
        return nothing
    end
    nin_path  = resolve(tapes, real_nin)
    nout_path = resolve(tapes, nout)
    @info "moder: extract-mode (flag=$nin_flag) tape $(abs(real_nin)) → tape $(abs(nout)) " *
          "(material filter not implemented)"
    if isfile(nin_path)
        nin_path != nout_path && cp(nin_path, nout_path; force=true)
        register!(tapes, nout, nout_path)
    else
        @warn "moder: extract-mode input tape $nin_path not found"
    end
    nothing
end

"""
    final_assembly!(tapes::TapeManager, output_unit::Int,
                    reconr_result, override_mf3, extra_mf3,
                    mf6_records, mf6_xsi, mf6_emax, mf6_stubs,
                    mf12_lines, mf13_lines, descriptions, thermr_mts,
                    thermr_coh_ne; mat, err, tempr, label)

Assemble the final PENDF tape using write_full_pendf with all accumulated
module outputs. This is called instead of a simple file copy when the
output requires MF6, MF12/MF13, and proper MF1 directory formatting.
"""
function final_assembly!(tapes::TapeManager, output_unit::Int,
                         reconr_result, override_mf3, extra_mf3,
                         mf6_records, mf6_xsi, mf6_emax, mf6_stubs,
                         mf12_lines, mf13_lines, descriptions, thermr_mts,
                         thermr_coh_ne; mat, err, tempr, label="")
    output_path = resolve(tapes, output_unit)
    @info "final_assembly: writing $output_path"

    open(output_path, "w") do io
        write_full_pendf(io, reconr_result;
            mat=mat, label=label, err=err, tempr=tempr,
            override_mf3=override_mf3, extra_mf3=extra_mf3,
            mf6_records=mf6_records, mf6_stubs=mf6_stubs,
            mf12_lines=mf12_lines, mf13_lines=mf13_lines,
            mf6_xsi=mf6_xsi, mf6_emax=mf6_emax,
            thermr_mts=thermr_mts,
            thermr_coh_ne=thermr_coh_ne,
            descriptions=descriptions)
    end

    register!(tapes, output_unit, output_path)
    lines = countlines(output_path)
    @info "final_assembly: $lines lines written"
    nothing
end
