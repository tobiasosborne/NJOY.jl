# viewr module runner -- PostScript rendering from NJOY plot tapes
#
# Matches Fortran viewr.f90 + graph.f90: reads plot tape produced by
# dtfr/plotr/covr, renders to PostScript output.
#
# The actual rendering engine lives in src/viewr/*.jl and is assembled
# by including those files in order. This orchestration module just
# wires the tape manager interface to viewr_render!.

"""
    viewr_module(tapes::TapeManager, params::ViewrParams)

Run VIEWR: render NJOY plot tape to PostScript.
Reads plot commands from infile, writes PostScript to nps.
"""
function viewr_module(tapes::TapeManager, params::ViewrParams)
    @info "viewr: infile=$(params.infile) nps=$(params.nps)"

    infile_path = resolve(tapes, params.infile)
    nps_path = resolve(tapes, params.nps)

    # Read the plot tape and render to PostScript
    open(nps_path, "w") do ps_io
        viewr_render!(ps_io, infile_path)
    end
    register!(tapes, params.nps, nps_path)

    lines = countlines(nps_path)
    @info "viewr: wrote $nps_path ($lines lines)"
    nothing
end
