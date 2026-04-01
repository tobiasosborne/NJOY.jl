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

    nin_path = resolve(tapes, nin)
    nout_path = resolve(tapes, nout)

    @info "moder: tape $(abs(nin)) → tape $(abs(nout))"

    # If input tape exists, copy it to output location
    if isfile(nin_path)
        if nin_path != nout_path
            cp(nin_path, nout_path; force=true)
        end
        # Register the output tape so downstream modules can find it
        register!(tapes, nout, nout_path)
    else
        @warn "moder: input tape $nin_path not found"
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
