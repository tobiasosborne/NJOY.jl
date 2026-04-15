# gaspr module runner — Gas production cross sections.
#
# Matches Fortran gaspr.f90 pipeline entry point: reads ENDF + PENDF,
# computes MT=203/204/205/206/207 gas-production cross sections from
# MF6 distributions and MT thresholds, writes an augmented PENDF.
#
# Current coverage: STUB. The gas-production algorithm is ported in
# src/processing/gaspr.jl (`gas_production`, `accumulate_gas`), but
# the "splice MT20x sections into a PENDF tape" step is not wired up.
# For now: copy npendf_in → npendf_out so downstream modules (acer,
# moder) find a readable file.

"""
    gaspr_module(tapes::TapeManager, params::GasprParams)

STUB: Copy input PENDF (`params.npendf_in`) to output PENDF
(`params.npendf_out`) so downstream modules don't crash on missing
tapes. A future implementation will compute gas-production MT20x
sections and splice them into the output.
"""
function gaspr_module(tapes::TapeManager, params::GasprParams)
    @info "gaspr: STUB nendf=$(params.nendf) npendf_in=$(params.npendf_in) " *
          "npendf_out=$(params.npendf_out) (MT20x gas production not yet emitted)"

    if params.npendf_out <= 0
        @warn "gaspr: no output unit — nothing to do"
        return nothing
    end

    out_path = resolve(tapes, params.npendf_out)

    if params.npendf_in > 0
        in_path = resolve(tapes, params.npendf_in)
        if isfile(in_path)
            cp(in_path, out_path; force=true)
            register!(tapes, params.npendf_out, out_path)
            @info "gaspr: copied PENDF $(in_path) → $(out_path)"
            return nothing
        else
            @warn "gaspr: input PENDF (unit $(params.npendf_in)) not found at $in_path"
        end
    end

    touch(out_path)
    register!(tapes, params.npendf_out, out_path)
    nothing
end
