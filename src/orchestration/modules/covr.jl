# covr module runner — Covariance plot-tape generation.
#
# Matches Fortran covr.f90 pipeline entry point: reads errorr's covariance
# tape and writes a plotr-format tape consumed by viewr to produce the
# PostScript covariance-matrix plots.
#
# Current coverage: STUB. Low-level covariance → correlation conversion
# and boxer/csv formatting live in src/processing/covr.jl, but the
# plotr-format emitter (the big block of plot commands in Fortran covr.f90
# that viewr reads) is not yet wired up. For now: write an empty output
# tape at params.nout so downstream modules (viewr, moder) find a file to
# open rather than crashing with SystemError. Classification moves from
# CRASH to DIFFS in the reference sweep.

"""
    covr_module(tapes::TapeManager, params::CovrParams)

STUB: Write an empty output tape at `params.nout` so downstream modules
(viewr, moder) don't crash on a missing file. A future implementation
will consume the errorr covariance tape at `params.nin` and emit plotr-
format commands matching Fortran covr.f90's plot-tape output.
"""
function covr_module(tapes::TapeManager, params::CovrParams)
    @info "covr: STUB nin=$(params.nin) npend=$(params.npend) nout=$(params.nout) " *
          "(plot-tape generation not yet emitted)"

    if params.nout <= 0
        @warn "covr: no output unit — nothing to do"
        return nothing
    end

    out_path = resolve(tapes, params.nout)
    touch(out_path)
    register!(tapes, params.nout, out_path)
    nothing
end
