# plotr module runner — Plot-tape generation for later viewr rendering.
#
# Matches Fortran plotr.f90 pipeline entry point: reads ENDF (or covariance
# tape) and emits a plot-command tape that viewr turns into PostScript.
#
# Current coverage: STUB. The plot-command format (BRK/FIT/CURVE cards that
# viewr reads) is not yet emitted. For now: touch the output unit so
# downstream viewr finds a file rather than SystemError-ing.

"""
    plotr_module(tapes::TapeManager, params::PlotrParams)

STUB: Write an empty output tape at `params.nplt` so downstream `viewr`
can open it (it will produce an empty PostScript which still trips line-
count diff but no longer crashes the pipeline). A future implementation
will emit plot-command cards matching Fortran `plotr.f90`.
"""
function plotr_module(tapes::TapeManager, params::PlotrParams)
    @info "plotr: STUB nplt=$(params.nplt) (plot-command tape not yet emitted)"

    if params.nplt <= 0
        @warn "plotr: no output unit — nothing to do"
        return nothing
    end

    out_path = resolve(tapes, params.nplt)
    touch(out_path)
    register!(tapes, params.nplt, out_path)
    nothing
end
