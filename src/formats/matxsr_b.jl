# MATXSR-B: Type-safe MATXS binary writer (Proposer-B)
# Records: FileId, FileControl, HollID, FileData, GroupStruct, MatCtrl, Vec, Matrix

using Printf

# Record types
struct MATXSParticleB
    name::String; ngrp::Int32
    function MATXSParticleB(name::AbstractString, ngrp::Integer)
        ngrp > 0 || throw(ArgumentError("ngrp must be > 0, got $ngrp"))
        new(rpad(name, 8)[1:8], Int32(ngrp))
    end
end

struct MATXSDataTypeB
    name::String; jinp::Int32; joutp::Int32
    function MATXSDataTypeB(name::AbstractString, jinp::Integer, joutp::Integer)
        jinp >= 0 || throw(ArgumentError("jinp must be >= 0"))
        joutp > 0 || throw(ArgumentError("joutp must be > 0"))
        new(rpad(name, 8)[1:8], Int32(jinp), Int32(joutp))
    end
end

struct MATXSVectorB
    name::String; nfg::Int32; nlg::Int32; data::Vector{Float32}
    function MATXSVectorB(name::AbstractString, nfg::Integer, nlg::Integer,
                          data::AbstractVector{<:Real})
        nfg > 0 || throw(ArgumentError("nfg must be > 0"))
        nlg >= nfg || throw(ArgumentError("nlg=$nlg must be >= nfg=$nfg"))
        expected = nlg - nfg + 1
        length(data) == expected || throw(ArgumentError(
            "data length $(length(data)) != expected $expected"))
        new(rpad(name, 8)[1:8], Int32(nfg), Int32(nlg), Float32.(data))
    end
end

struct MATXSMatrixBlockB
    name::String; lord::Int32; jconst::Int32
    jband::Vector{Int32}; ijj::Vector{Int32}
    data::Vector{Float32}; const_spec::Vector{Float32}; const_prod::Vector{Float32}
    function MATXSMatrixBlockB(name, lord, jconst, jband, ijj, data, spec, prod)
        lord > 0 || throw(ArgumentError("lord must be > 0"))
        jconst >= 0 || throw(ArgumentError("jconst must be >= 0"))
        length(jband) == length(ijj) || throw(ArgumentError("jband/ijj length mismatch"))
        new(rpad(name, 8)[1:8], Int32(lord), Int32(jconst),
            Int32.(jband), Int32.(ijj), Float32.(data), Float32.(spec), Float32.(prod))
    end
end

struct MATXSSubmaterialB
    temp::Float32; sigz::Float32; itype::Int32
    vectors::Vector{MATXSVectorB}; matrices::Vector{MATXSMatrixBlockB}
end

struct MATXSMaterialB
    name::String; amass::Float32; submaterials::Vector{MATXSSubmaterialB}
    function MATXSMaterialB(name::AbstractString, amass::Real,
                            subs::AbstractVector{MATXSSubmaterialB})
        amass > 0 || throw(ArgumentError("amass must be > 0"))
        new(rpad(name, 8)[1:8], Float32(amass), subs)
    end
end

"""Complete MATXS file with validated structure."""
struct MATXSFileB
    huse1::String; huse2::String; ivers::Int32; hsetid::Vector{String}
    particles::Vector{MATXSParticleB}; data_types::Vector{MATXSDataTypeB}
    materials::Vector{MATXSMaterialB}; group_bounds::Vector{Vector{Float32}}
    maxw::Int32
    function MATXSFileB(h1, h2, iv, hs, parts, dtypes, mats, gbs, mw)
        np, nt, nm = length(parts), length(dtypes), length(mats)
        np > 0 || throw(ArgumentError("need >= 1 particle"))
        nt > 0 || throw(ArgumentError("need >= 1 data type"))
        nm > 0 || throw(ArgumentError("need >= 1 material"))
        length(gbs) == np || throw(ArgumentError("group_bounds count != npart"))
        for (j, gb) in enumerate(gbs)
            length(gb) == Int(parts[j].ngrp)+1 || throw(ArgumentError(
                "particle $j: bounds length != ngrp+1"))
        end
        for (k, dt) in enumerate(dtypes)
            dt.jinp > 0 && dt.jinp > np && throw(ArgumentError("dtype $k: jinp > npart"))
            dt.joutp > np && throw(ArgumentError("dtype $k: joutp > npart"))
        end
        for (i, mat) in enumerate(mats)
            for (j, sub) in enumerate(mat.submaterials)
                sub.itype > nt && throw(ArgumentError("mat $i sub $j: itype > ntype"))
            end
        end
        mw > 0 || throw(ArgumentError("maxw must be > 0"))
        new(rpad(h1,8)[1:8], rpad(h2,8)[1:8], Int32(iv),
            [rpad(s,8)[1:8] for s in hs], parts, dtypes, mats,
            [Float32.(gb) for gb in gbs], Int32(mw))
    end
end

# File length in records
function _matxs_length_b(f::MATXSFileB)
    nrec = 4 + length(f.particles)
    for mat in f.materials
        nrec += 1
        for sub in mat.submaterials
            if length(sub.vectors) > 0; nrec += 2; end
            for blk in sub.matrices
                nrec += 2
                if blk.jconst > 0; nrec += 1; end
            end
        end
    end
    nrec
end

# Binary helpers
function _write_matxs_rec_b(io::IO, buf::Vector{UInt8})
    n = Int32(length(buf)); write(io, n); write(io, buf); write(io, n)
end
_mx_emit!(buf, x::Int32)   = append!(buf, reinterpret(UInt8, [x]))
_mx_emit!(buf, x::Float32) = append!(buf, reinterpret(UInt8, [x]))
_mx_emit!(buf, s::String)  = append!(buf, UInt8.(collect(rpad(s, 8)[1:8])))
_mx_emit_ints!(buf, vs)  = for v in vs; _mx_emit!(buf, Int32(v)); end
_mx_emit_reals!(buf, vs) = for v in vs; _mx_emit!(buf, Float32(v)); end

"""
    write_matxs_b(io::IO, f::MATXSFileB)

Write MATXS in CCCC binary format with Fortran record markers.
"""
function write_matxs_b(io::IO, f::MATXSFileB)
    np, nt, nm = length(f.particles), length(f.data_types), length(f.materials)
    nholl = length(f.hsetid)
    # Record 1: File identification
    buf = UInt8[]
    _mx_emit!(buf, rpad("matxs",8)[1:8]); _mx_emit!(buf, f.huse1)
    _mx_emit!(buf, f.huse2); _mx_emit!(buf, f.ivers)
    _write_matxs_rec_b(io, buf)
    # Record 2: File control
    buf = UInt8[]
    for v in (Int32(np), Int32(nt), Int32(nholl), Int32(nm), f.maxw, Int32(_matxs_length_b(f)))
        _mx_emit!(buf, v)
    end
    _write_matxs_rec_b(io, buf)
    # Record 3: Hollerith ID
    buf = UInt8[]
    for s in f.hsetid; _mx_emit!(buf, s); end
    _write_matxs_rec_b(io, buf)
    # Record 4: File data
    buf = UInt8[]
    for p in f.particles; _mx_emit!(buf, p.name); end
    for dt in f.data_types; _mx_emit!(buf, dt.name); end
    for m in f.materials; _mx_emit!(buf, m.name); end
    for p in f.particles; _mx_emit!(buf, p.ngrp); end
    for dt in f.data_types; _mx_emit!(buf, dt.jinp); end
    for dt in f.data_types; _mx_emit!(buf, dt.joutp); end
    for m in f.materials; _mx_emit!(buf, Int32(length(m.submaterials))); end
    loc = 1
    for m in f.materials
        _mx_emit!(buf, Int32(loc)); loc += 1
        for sub in m.submaterials
            if length(sub.vectors) > 0; loc += 2; end
            loc += 2*length(sub.matrices)
            loc += count(b->b.jconst>0, sub.matrices)
        end
    end
    _write_matxs_rec_b(io, buf)
    # Record 5: Group structures
    for (j, p) in enumerate(f.particles)
        buf = UInt8[]; _mx_emit_reals!(buf, f.group_bounds[j]); _write_matxs_rec_b(io, buf)
    end
    # Material records
    for mat in f.materials; _write_material_b(io, mat, f); end
    return nothing
end

function _write_material_b(io::IO, mat::MATXSMaterialB, f::MATXSFileB)
    buf = UInt8[]; _mx_emit!(buf, mat.name); _mx_emit!(buf, mat.amass)
    for sub in mat.submaterials
        _mx_emit!(buf, sub.temp); _mx_emit!(buf, sub.sigz); _mx_emit!(buf, sub.itype)
        _mx_emit!(buf, Int32(length(sub.vectors)))
        _mx_emit!(buf, Int32(length(sub.matrices))); _mx_emit!(buf, Int32(0))
    end
    _write_matxs_rec_b(io, buf)
    for sub in mat.submaterials; _write_sub_b(io, sub); end
end

function _write_sub_b(io::IO, sub::MATXSSubmaterialB)
    if length(sub.vectors) > 0
        buf = UInt8[]
        for v in sub.vectors; _mx_emit!(buf, v.name); end
        for v in sub.vectors; _mx_emit!(buf, v.nfg); end
        for v in sub.vectors; _mx_emit!(buf, v.nlg); end
        _write_matxs_rec_b(io, buf)
        buf = UInt8[]
        for v in sub.vectors; _mx_emit_reals!(buf, v.data); end
        _write_matxs_rec_b(io, buf)
    end
    for blk in sub.matrices
        buf = UInt8[]
        _mx_emit!(buf, blk.name); _mx_emit!(buf, blk.lord); _mx_emit!(buf, blk.jconst)
        _mx_emit_ints!(buf, blk.jband); _mx_emit_ints!(buf, blk.ijj)
        _write_matxs_rec_b(io, buf)
        if !isempty(blk.data)
            buf = UInt8[]; _mx_emit_reals!(buf, blk.data); _write_matxs_rec_b(io, buf)
        end
        if blk.jconst > 0 && (!isempty(blk.const_spec) || !isempty(blk.const_prod))
            buf = UInt8[]
            _mx_emit_reals!(buf, blk.const_spec); _mx_emit_reals!(buf, blk.const_prod)
            _write_matxs_rec_b(io, buf)
        end
    end
end

# MT name mapping
const _MATXS_MT_NAMES_B = Dict{Int,String}(
    1=>"ntot", 2=>"nelas", 4=>"ninel", 16=>"n2n", 17=>"n3n", 18=>"nfiss",
    102=>"ngamm", 103=>"np", 104=>"nd", 105=>"nt", 106=>"nhe3", 107=>"nalpha",
    251=>"mubar", 252=>"xi", 253=>"gamma", 452=>"nubar")
_mt_name_b(mt::Integer) = get(_MATXS_MT_NAMES_B, Int(mt), @sprintf("mt%d", mt))

"""
    build_matxs_b(mg::MultiGroupXS; kwargs...) -> MATXSFileB

Build validated MATXS from multigroup cross section data.
"""
function build_matxs_b(mg::MultiGroupXS;
        huse::AbstractString="", hsetid::AbstractVector{<:AbstractString}=["NJOY.jl"],
        ivers::Integer=0, mat_name::AbstractString="mat1",
        amass::Real=1.0, maxw::Integer=5000)
    gb = mg.group_bounds; ng = length(gb) - 1
    gbounds = Float32[Float32(gb[ng-g+2]) for g in 1:ng]
    push!(gbounds, Float32(gb[1]))
    mt_idx = Dict(mt=>c for (c,mt) in enumerate(mg.mt_list))
    vectors = MATXSVectorB[]
    for (mt, col) in mt_idx
        data = Float32[Float32(mg.xs[g,col]) for g in 1:ng]
        fnz, lnz = findfirst(!=(0), data), findlast(!=(0), data)
        fnz !== nothing && lnz !== nothing &&
            push!(vectors, MATXSVectorB(_mt_name_b(mt), fnz, lnz, data[fnz:lnz]))
    end
    matrices = MATXSMatrixBlockB[]
    if haskey(mt_idx, 2)
        sd = Float32[Float32(mg.xs[g,mt_idx[2]]) for g in 1:ng]
        push!(matrices, MATXSMatrixBlockB("nelas", Int32(1), Int32(0),
            ones(Int32,ng), Int32.(1:ng), sd, Float32[], Float32[]))
    end
    sub = MATXSSubmaterialB(Float32(300), Float32(1e10), Int32(1), vectors, matrices)
    mat = MATXSMaterialB(mat_name, amass, [sub])
    h1 = length(huse) >= 8 ? huse[1:8] : huse
    h2 = length(huse) > 8 ? huse[9:min(end,16)] : ""
    MATXSFileB(h1, h2, ivers, hsetid, [MATXSParticleB("n",ng)],
               [MATXSDataTypeB("nscat",1,1)], [mat], [gbounds], maxw)
end
