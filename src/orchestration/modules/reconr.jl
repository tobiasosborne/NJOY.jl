# reconr module runner -- Reconstruct pointwise cross sections from ENDF
#
# Thin wrapper: calls existing reconr(), writes PENDF via write_pendf_file().
# All 19 bit-identical reconr tests validate this path.

"""
    reconr_module(tapes::TapeManager, params::ReconrParams)

Run RECONR: read ENDF tape, reconstruct pointwise cross sections,
write PENDF tape. Matches Fortran reconr.f90 interface.

Input tapes: nendf (ENDF evaluation)
Output tapes: npend (pointwise ENDF = PENDF)
"""
function reconr_module(tapes::TapeManager, params::ReconrParams)
    endf_path = resolve(tapes, params.nendf)
    pendf_path = resolve(tapes, params.npend)

    @info "reconr: MAT=$(params.mat) err=$(params.err) → $pendf_path"

    r = reconr(endf_path; mat=params.mat, err=params.err)

    write_pendf_file(pendf_path, r; mat=params.mat, err=params.err)

    @info "reconr: $(length(r.energies)) points written"
    nothing
end
