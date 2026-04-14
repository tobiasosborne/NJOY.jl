# leapr module runner — S(α,β) generation from phonon DOS.
#
# Matches Fortran leapr.f90 pipeline entry point: reads phonon density of
# states from input deck, runs `generate_sab` + post-processing, writes
# ENDF MF7 file.
#
# Current coverage: STUB. The phonon-expansion algorithm is ported in
# src/processing/leapr.jl (`generate_sab`), but parsing the full leapr
# input deck (oscillators, T-effective, Bragg cards, descriptions) and
# emitting MF7/MT2+MT4 ENDF format is not yet wired up. For now: write
# an empty output tape so downstream modules find a file at the expected
# unit number.

"""
    leapr_module(tapes::TapeManager, params::LeaprParams)

STUB: Write an empty output tape at `params.nout` so downstream modules
(thermr, moder) don't crash on missing files. A future implementation
will run `generate_sab` and emit MF7 ENDF sections.
"""
function leapr_module(tapes::TapeManager, params::LeaprParams)
    @info "leapr: STUB nout=$(params.nout) (S(α,β) generation not yet emitted)"

    if params.nout <= 0
        @warn "leapr: no output unit — nothing to do"
        return nothing
    end

    out_path = resolve(tapes, params.nout)
    touch(out_path)
    register!(tapes, params.nout, out_path)
    nothing
end
