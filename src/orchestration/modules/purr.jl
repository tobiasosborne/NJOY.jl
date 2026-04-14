# purr module runner — Probability tables for the unresolved resonance range.
#
# Matches Fortran purr.f90 pipeline entry point: reads ENDF (MF2 unresolved)
# and PENDF (MF3 backgrounds), generates MC probability tables (MT152),
# writes augmented PENDF.
#
# Current coverage: STUB. The low-level `generate_ptable` algorithm is ported
# in src/processing/purr.jl, but turning a `ProbabilityTable` into the
# Fortran's MT152 ENDF section format is not yet wired up. For now: copy
# input PENDF to output PENDF unchanged so downstream modules (acer, viewr)
# find a tape to read.

"""
    purr_module(tapes::TapeManager, params::PurrParams)

STUB: Copy input PENDF (`params.npendf_in`) to output PENDF
(`params.npendf_out`) so downstream modules don't crash on missing tapes.
A future implementation will run `generate_ptable` and emit MT152 sections.
"""
function purr_module(tapes::TapeManager, params::PurrParams)
    @info "purr: STUB nendf=$(params.nendf) npendf_in=$(params.npendf_in) " *
          "npendf_out=$(params.npendf_out) mat=$(params.mat) " *
          "(probability tables not yet emitted)"

    if params.npendf_out <= 0
        @warn "purr: no output unit — nothing to do"
        return nothing
    end

    out_path = resolve(tapes, params.npendf_out)

    if params.npendf_in > 0
        in_path = resolve(tapes, params.npendf_in)
        if isfile(in_path)
            cp(in_path, out_path; force=true)
            register!(tapes, params.npendf_out, out_path)
            @info "purr: copied PENDF $(in_path) → $(out_path)"
            return nothing
        else
            @warn "purr: input PENDF (unit $(params.npendf_in)) not found at $in_path"
        end
    end

    # Fallback: empty file so downstream open() doesn't SystemError
    touch(out_path)
    register!(tapes, params.npendf_out, out_path)
    nothing
end
