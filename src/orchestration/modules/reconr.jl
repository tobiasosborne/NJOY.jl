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

    materials = params.materials
    isempty(materials) && (materials = [ReconrMatSpec(params.mat, params.err, "")])

    if length(materials) == 1
        mspec = materials[1]
        @info "reconr: MAT=$(mspec.mat) err=$(mspec.err) → $pendf_path"
        r = reconr(endf_path; mat=mspec.mat, err=mspec.err)
        write_pendf_file(pendf_path, r; mat=mspec.mat, err=mspec.err)
        @info "reconr: $(length(r.energies)) points written"
    else
        open(pendf_path, "w") do io
            _write_tpid_line(io, "PENDF tape produced by NJOY.jl", materials[1].mat)
            for mspec in materials
                @info "reconr: MAT=$(mspec.mat) err=$(mspec.err)"
                r = reconr(endf_path; mat=mspec.mat, err=mspec.err)
                _write_reconr_mat_block(io, r, mspec.mat, mspec.err)
                @info "reconr: MAT=$(mspec.mat) $(length(r.energies)) points"
            end
            blanks = repeat(" ", 66)
            @printf(io, "%s%4d%2d%3d%5d\n", blanks, -1, 0, 0, 0)
        end
        @info "reconr: wrote $pendf_path"
    end
    nothing
end

"""Write a single MAT block (MF1+MF2+MF3+MEND) for reconr output, no TPID/TEND."""
function _write_reconr_mat_block(io::IO, result::NamedTuple, mat::Int, err::Float64)
    mf2 = result.mf2
    actual_mat = mat > 0 ? Int32(mat) : Int32(max(1, round(Int, mf2.ZA / 10)))
    ns = Ref(1)

    reactions = _collect_reactions(result)
    mt_mf = Dict{Int,Int}()
    for sec in result.mf3_sections
        mt_mf[Int(sec.mt)] = Int(sec.mf)
    end

    _write_legacy_mf1(io, mf2, actual_mat, err, 0.0, length(result.energies), reactions, ns;
                      mt_to_mf=mt_mf)
    _write_legacy_mf2(io, mf2, actual_mat, ns)
    _write_legacy_mf3(io, result, actual_mat, reactions, ns)
    _write_fend_zero(io, 0)  # MEND
end
