# RESXSR -- Resonance cross section file (RESXS format)
#
# Extracts epithermal resonance-range cross sections (elastic, fission,
# capture) from PENDF data and writes them in a tabular format suitable
# for hyper-fine flux calculations.
#
# Correspondence to NJOY2016 resxsr.f90:
#   energy selection + MT extraction (lines 296-347) -> extract_resxs
#   thinning loop (lines 260-397)                   -> thin_resxs
#   RESXS output (lines 436-501)                    -> write_resxs
#
# Design: extraction and thinning are pure functions; only write_resxs
# performs I/O. Thinning uses only arithmetic (AD-compatible).

# ==========================================================================
# Result type
# ==========================================================================

"""
    RESXSRecord

Extracted resonance cross sections for one material.
"""
struct RESXSRecord
    name::String
    mat::Int
    awr::Float64
    nreac::Int                   # 2 (non-fissile) or 3 (fissile)
    energies::Vector{Float64}
    elastic::Vector{Float64}     # MT2
    fission::Vector{Float64}     # MT18 (zeros if not fissile)
    capture::Vector{Float64}     # MT102
end

# ==========================================================================
# Thinning  (resxsr.f90 lines 260-397)
# ==========================================================================

"""
    thin_resxs(energies, xs_cols::NTuple{N,Vector{Float64}}, eps) ->
        (thinned_e, thinned_cols)

Thin pointwise cross sections using linear-interpolation tolerance `eps`.
All `N` reaction columns must satisfy the tolerance for a point to be removed.
Pure function, AD-compatible.
"""
function thin_resxs(energies::AbstractVector{<:Real},
                    xs_cols::NTuple{N, AbstractVector{<:Real}},
                    eps::Real) where N
    ne = length(energies)
    ne <= 2 && return (Float64.(energies), Tuple(Float64.(c) for c in xs_cols))

    keep = Int[1]  # always keep first point
    base = 1

    for i in 2:(ne - 1)
        can_skip = true
        e_base = Float64(energies[base])
        e_end  = Float64(energies[i + 1])
        e_i    = Float64(energies[i])

        if e_end == e_base
            can_skip = false
        else
            frac = (e_i - e_base) / (e_end - e_base)
            for k in 1:N
                y_base  = Float64(xs_cols[k][base])
                y_end   = Float64(xs_cols[k][i + 1])
                y_interp = y_base + frac * (y_end - y_base)
                y_actual = Float64(xs_cols[k][i])
                if abs(y_interp - y_actual) > eps * abs(y_actual)
                    can_skip = false
                    break
                end
            end
        end

        if !can_skip
            push!(keep, i)
            base = i
        end
    end

    push!(keep, ne)  # always keep last

    out_e = Float64[energies[i] for i in keep]
    out_cols = Tuple(Float64[c[i] for i in keep] for c in xs_cols)
    (out_e, out_cols)
end

# ==========================================================================
# Extraction (pure function)
# ==========================================================================

"""
    extract_resxs(pendf::PointwiseMaterial;
                  elow=1.0, ehigh=1.0e5, eps=0.005,
                  name="", awr=1.0) -> RESXSRecord

Extract resonance-range cross sections from a PointwiseMaterial.
Selects MT2, MT18, MT102 within [elow, ehigh] and thins with tolerance `eps`.
Pure function (no I/O).
"""
function extract_resxs(pendf::PointwiseMaterial;
                       elow::Real=1.0, ehigh::Real=1.0e5,
                       eps::Real=0.005,
                       name::String="",
                       awr::Real=1.0)
    col_el  = findfirst(==(2),   pendf.mt_list)
    col_fis = findfirst(==(18),  pendf.mt_list)
    col_cap = findfirst(==(102), pendf.mt_list)

    # Energy range selection
    idx = findall(e -> elow <= e <= ehigh, pendf.energies)

    if isempty(idx)
        nreac = col_fis !== nothing ? 3 : 2
        return RESXSRecord(name, Int(pendf.mat), Float64(awr), nreac,
                           Float64[elow], Float64[0.0], Float64[0.0], Float64[0.0])
    end

    e_sel   = pendf.energies[idx]
    el_sel  = col_el  !== nothing ? pendf.cross_sections[idx, col_el]  : zeros(length(idx))
    fis_sel = col_fis !== nothing ? pendf.cross_sections[idx, col_fis] : zeros(length(idx))
    cap_sel = col_cap !== nothing ? pendf.cross_sections[idx, col_cap] : zeros(length(idx))

    nreac = col_fis !== nothing ? 3 : 2

    # Thin using all reaction columns
    cols = nreac == 3 ? (el_sel, fis_sel, cap_sel) : (el_sel, cap_sel)
    e_thin, cols_thin = thin_resxs(e_sel, cols, Float64(eps))

    if nreac == 3
        return RESXSRecord(name, Int(pendf.mat), Float64(awr), nreac,
                           e_thin, cols_thin[1], cols_thin[2], cols_thin[3])
    else
        return RESXSRecord(name, Int(pendf.mat), Float64(awr), nreac,
                           e_thin, cols_thin[1], zeros(length(e_thin)), cols_thin[2])
    end
end

# ==========================================================================
# Dict-based extraction (composable, AD-compatible)
# ==========================================================================

"""
    extract_resxs_dict(energies, reactions;
                       elow=1.0, ehigh=1.0e5, eps=0.005)
        -> (thinned_e, thinned_reactions::Dict{Int,Vector{Float64}})

Extract and thin resonance cross sections from a Dict-based representation.
Only considers MT2, MT18, MT102. Pure function.
"""
function extract_resxs_dict(energies::AbstractVector{<:Real},
                            reactions::Dict{Int,<:AbstractVector{<:Real}};
                            elow::Real=1.0, ehigh::Real=1.0e5,
                            eps::Real=0.005)
    idx = findall(e -> elow <= e <= ehigh, energies)
    isempty(idx) && return (Float64[elow], Dict{Int,Vector{Float64}}())

    e_sel = Float64.(energies[idx])
    resxs_mts = [2, 18, 102]
    present_mts = Int[]
    cols = Vector{Float64}[]

    for mt in resxs_mts
        xs = get(reactions, mt, nothing)
        xs === nothing && continue
        push!(present_mts, mt)
        push!(cols, Float64.(xs[idx]))
    end

    isempty(cols) && return (e_sel, Dict{Int,Vector{Float64}}())

    e_thin, cols_thin = thin_resxs(e_sel, Tuple(cols), Float64(eps))

    result = Dict{Int,Vector{Float64}}()
    for (k, mt) in enumerate(present_mts)
        result[mt] = cols_thin[k]
    end

    (e_thin, result)
end

# ==========================================================================
# Writer (I/O)
# ==========================================================================

"""
    write_resxs(io::IO, records::Vector{RESXSRecord};
                user_id="NJOY.jl", ivers=1,
                description="Resonance cross sections")

Write one or more RESXS records in text format.
"""
function write_resxs(io::IO, records::Vector{RESXSRecord};
                     user_id::String="NJOY.jl",
                     ivers::Int=1,
                     description::String="Resonance cross sections")
    isempty(records) && return nothing
    elow  = minimum(r.energies[1] for r in records)
    ehigh = maximum(r.energies[end] for r in records)
    nmat  = length(records)

    # File identification
    println(io, "resxs  ", rpad(user_id, 16), " ", ivers)

    # File control
    @Printf.printf(io, "%14.6e %14.6e %6d %6d\n", elow, ehigh, 1, nmat)

    # Description
    println(io, rpad(description, 72))

    # Material directory
    for r in records
        println(io, rpad(r.name, 8), "  ntemp=1  nreac=", r.nreac,
                "  nener=", length(r.energies))
    end

    # Material data blocks
    for r in records
        ne = length(r.energies)
        @Printf.printf(io, "# %s  awr=%12.5e  nreac=%d  nener=%d\n",
                       rpad(r.name, 8), r.awr, r.nreac, ne)
        for i in 1:ne
            if r.nreac == 3
                @Printf.printf(io, "%14.6e %14.6e %14.6e %14.6e\n",
                               r.energies[i], r.elastic[i],
                               r.fission[i], r.capture[i])
            else
                @Printf.printf(io, "%14.6e %14.6e %14.6e\n",
                               r.energies[i], r.elastic[i], r.capture[i])
            end
        end
    end

    nothing
end

"""
    write_resxs(io::IO, pendf::PointwiseMaterial; kwargs...)

Convenience: extract and write a single material.
"""
function write_resxs(io::IO, pendf::PointwiseMaterial;
                     elow::Real=1.0, ehigh::Real=1.0e5,
                     eps::Real=0.005,
                     name::String="material",
                     awr::Real=1.0,
                     user_id::String="NJOY.jl",
                     ivers::Int=1,
                     description::String="Resonance cross sections")
    rec = extract_resxs(pendf; elow=Float64(elow), ehigh=Float64(ehigh),
                        eps=Float64(eps), name=name, awr=Float64(awr))
    write_resxs(io, [rec]; user_id=user_id, ivers=ivers, description=description)
end

"""
    write_resxs(filename::AbstractString, args...; kwargs...)

Write RESXS data to a file.
"""
function write_resxs(filename::AbstractString, args...; kwargs...)
    open(filename, "w") do io
        write_resxs(io, args...; kwargs...)
    end
end
