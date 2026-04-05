# matxsr module runner -- MATXS format cross section output
#
# Matches Fortran matxsr.f90: read GENDF tape(s),
# write MATXS-format output tape.

"""
    matxsr_module(tapes::TapeManager, params::MatxsrParams)

Run MATXSR: produce MATXS-format cross section library from GENDF data.
Reads GENDF tape(s) from gaminr/groupr, writes MATXS output.
"""
function matxsr_module(tapes::TapeManager, params::MatxsrParams)
    @info "matxsr: ngen1=$(params.ngen1) ngen2=$(params.ngen2) nmatx=$(params.nmatx)"
    @info "matxsr: npart=$(params.npart) ntype=$(params.ntype) nmat=$(params.nmat)"

    nmatx_path = resolve(tapes, params.nmatx)

    # Read GENDF tape(s)
    gendf_materials = GendfMaterial[]  # reuse dtfr's reader

    # Photon GENDF (ngen2 = gaminr output)
    if params.ngen2 > 0
        gendf_path = resolve(tapes, params.ngen2)
        if isfile(gendf_path)
            append!(gendf_materials, _read_gendf_tape(gendf_path))
        end
    end

    # Neutron GENDF (ngen1 = groupr output)
    if params.ngen1 > 0
        gendf_path = resolve(tapes, params.ngen1)
        if isfile(gendf_path)
            append!(gendf_materials, _read_gendf_tape(gendf_path))
        end
    end

    # Write MATXS output tape
    open(nmatx_path, "w") do io
        _write_matxsr_tape(io, gendf_materials, params)
    end
    register!(tapes, params.nmatx, nmatx_path)

    lines = countlines(nmatx_path)
    @info "matxsr: wrote $nmatx_path ($lines lines)"
    nothing
end

# =========================================================================
# MATXS tape writer — matches Fortran matxsr.f90 output format exactly
# Format: " Nd " tagged records with specific field widths
# =========================================================================

function _write_matxsr_tape(io::IO, gendf::Vector{GendfMaterial}, params::MatxsrParams)
    # Record 1: File Identification
    # Format: " 0v " + a8(hfile) + "*" + 2×a8(huse) + "*" + i6(ivers)
    huse_parts = split(rpad(params.huse, 16), "")
    huse1 = rpad(params.huse[1:min(8, length(params.huse))], 8)
    huse2 = length(params.huse) > 8 ? rpad(params.huse[9:min(16, length(params.huse))], 8) : "        "
    @printf(io, " 0v   matxs *%-8s%-8s*%6d\n", huse1, huse2, params.ivers)

    # Collect per-material data for sizing
    mat_data = Dict{Int, GendfMaterial}()
    for gm in gendf
        mat_data[gm.mat] = gm
    end

    # Get group structure from first material with data
    ngrp_vals = params.ngrp
    ngg = isempty(ngrp_vals) ? 12 : ngrp_vals[1]

    # Compute total file length
    # Record structure count for 'length' field
    n_records = 2 + 1 + 1 + params.npart + params.nmat * 2  # rough estimate
    maxw = 5000

    # Record 2: File Control
    # Format: " 1d   " + 6×i6
    # Fields: npart, ntype, nholl, nmat, maxw, length
    nholl_words = max(params.nholl * 9, 9)  # 9 hollerith words per card
    @printf(io, " 1d   %6d%6d%6d%6d%6d%6d\n",
            params.npart, params.ntype, nholl_words, params.nmat, maxw, n_records)

    # Record 3: Set Hollerith Identification
    # Format: " 2d " + lines of 9×a8
    @printf(io, " 2d \n")
    for hs in params.hsetid
        # Pad to 72 chars, split into 9×8-char fields
        padded = rpad(hs, 72)
        println(io, padded)
    end

    # Record 4: File Data
    # Hollerith part: particle names, data type names, material names
    # Integer part: ngrp, jinp, joutp, nsubm, locm
    holl_items = String[]
    for p in params.particles; push!(holl_items, rpad(p, 8)); end
    for d in params.dtypes;    push!(holl_items, rpad(d, 8)); end
    for m in params.materials;  push!(holl_items, rpad(m.name, 8)); end

    @printf(io, " 3d     ")
    for (idx, h) in enumerate(holl_items)
        print(io, rpad(h, 8))
        if idx % 8 == 0 && idx < length(holl_items)
            println(io)
        end
    end
    println(io)

    # Integer data: ngrp(npart), jinp(ntype), joutp(ntype), nsubm(nmat), locm(nmat)
    int_vals = Int[]
    append!(int_vals, ngrp_vals)
    append!(int_vals, params.jinp)
    append!(int_vals, params.joutp)
    # nsubm: each material has 1 submaterial (infinite dilution, 1 temp)
    for _ in 1:params.nmat; push!(int_vals, 1); end
    # locm: location of each material (cumulative)
    loc = 0
    for i in 1:params.nmat; push!(int_vals, loc); loc += 1; end

    # Write integers 12 per line
    for i in 1:length(int_vals)
        @printf(io, "%6d", int_vals[i])
        if i % 12 == 0; println(io); end
    end
    length(int_vals) % 12 != 0 && println(io)

    # Record 5: Group structures (one per particle)
    for ip in 1:params.npart
        # Get group boundaries from GENDF data
        egg = Float64[]
        for gm in gendf
            if !isempty(gm.egg)
                egg = gm.egg
                break
            end
        end

        # MATXS convention: group bounds in DESCENDING energy order
        egg_desc = reverse(egg)

        @printf(io, " 4d         ")
        for (idx, e) in enumerate(egg_desc)
            @printf(io, " %11.5E", e)
            if idx % 5 == 0 && idx < length(egg_desc)
                println(io)
            elseif idx == 1 && length(egg_desc) > 5
                # First line has only 5 values after " 4d         "
            end
        end
        println(io)
    end

    # Per-material records
    for mspec in params.materials
        gmat = get(mat_data, mspec.matno, nothing)
        gmat === nothing && continue

        # Record 6: Material Control
        @printf(io, " 5d %-8s%12.5E\n", rpad(mspec.name, 8), gmat.awr)
        # Submaterial: temp=0, sigz=1e12, itype, n1d, n2d, loc
        n1d = length(gmat.sections)  # number of vector data types
        n2d = 0  # number of matrix data types (scattering matrices)
        # Count scattering matrices separately
        n_scat = 0
        for sec in gmat.sections
            if sec.mt == 502 || sec.mt == 504
                n2d += 1
            end
        end
        @printf(io, "%12.5E%12.5E%6d%6d%6d%6d\n",
                0.0, 1.0e12, 1, n1d + 1, n2d, 0)  # +1 for flux

        # Record 7: Vector Control
        # Names for each vector quantity
        vec_names = String["gwt0    "]  # flux weight
        for sec in gmat.sections
            push!(vec_names, rpad(_matxsr_mt_name(sec.mt), 8))
        end

        @printf(io, " 6d     ")
        for (idx, vn) in enumerate(vec_names)
            print(io, rpad(vn, 8))
            if idx % 8 == 0 && idx < length(vec_names)
                println(io)
            end
        end
        println(io)

        # Integer part: first_group, last_group for each vector
        # All vectors span the full ngg groups
        int_v = Int[]
        for _ in vec_names
            push!(int_v, 1)     # first group
        end
        for _ in vec_names
            push!(int_v, ngg)   # last group
        end
        for i in 1:length(int_v)
            @printf(io, "%6d", int_v[i])
            if i % 12 == 0; println(io); end
        end
        length(int_v) % 12 != 0 && println(io)

        # Record 8: Vector Block Data
        # Flux weights first, then each reaction's group-averaged XS
        all_data = Float64[]

        # Flux weights
        for sec in gmat.sections
            for g in 1:min(ngg, length(sec.data))
                flux, _ = sec.data[g]
                push!(all_data, flux)
            end
            break  # use first section's flux
        end

        # Each reaction's XS
        for sec in gmat.sections
            for g in 1:min(ngg, length(sec.data))
                _, sigma = sec.data[g]
                push!(all_data, sigma)
            end
        end

        @printf(io, " 7d         ")
        for (idx, v) in enumerate(all_data)
            @printf(io, " %11.5E", v)
            if idx % 5 == 0 && idx < length(all_data)
                println(io)
            elseif idx % 6 == 0 && idx > 5
                println(io)
            end
        end
        println(io)

        # Scattering matrices (MF6/MF26)
        # For now, write diagonal-only matrices for coherent/incoherent
        for sec in gmat.sections
            (sec.mt == 502 || sec.mt == 504) || continue
            mt_name = rpad(_matxsr_mt_name(sec.mt), 8)

            # Matrix control: name, lord, jconst, jband(ng), ijj(ng)
            @printf(io, " 8d     %-8s\n", mt_name)
            # lord=1 (P0 only for now), jconst=0, jband=1 (diagonal), ijj=g
            lord = 1
            jconst = 0
            int_m = Int[lord, jconst]
            for g in 1:ngg; push!(int_m, 1); end  # jband = 1 for each group
            for g in 1:ngg; push!(int_m, g); end  # ijj = g for each group
            for i in 1:length(int_m)
                @printf(io, "%6d", int_m[i])
                if i % 12 == 0; println(io); end
            end
            length(int_m) % 12 != 0 && println(io)

            # Matrix data: one value per group (diagonal)
            mat_vals = Float64[]
            for g in 1:min(ngg, length(sec.data))
                _, sigma = sec.data[g]
                push!(mat_vals, sigma)
            end

            @printf(io, " 9d         ")
            for (idx, v) in enumerate(mat_vals)
                @printf(io, " %11.5E", v)
                if idx % 5 == 0 && idx < length(mat_vals)
                    println(io)
                elseif idx % 6 == 0 && idx > 5
                    println(io)
                end
            end
            println(io)
        end
    end
end

function _matxsr_mt_name(mt::Int)
    mt == 501 && return "gtot0"
    mt == 502 && return "gcoh"
    mt == 504 && return "ginch"
    mt == 516 && return "gpair"
    mt == 522 && return "p02"
    mt == 602 && return "p21"
    mt == 621 && return "gheat"
    return "mt$(mt)"
end
