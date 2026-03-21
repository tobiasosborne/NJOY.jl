# CCCCR-B: Type-safe CCCC-IV ISOTXS binary writer (Proposer-B)
# Records: FileId, FileControl, FileData, per-isotope control/XS/scatter

using Printf

# Hollerith8: 8-byte padded string for CCCC Hollerith words
struct Hollerith8
    data::NTuple{8,UInt8}
end
Hollerith8(s::AbstractString) = Hollerith8(NTuple{8,UInt8}(UInt8.(collect(rpad(s, 8)[1:8]))))
Hollerith8() = Hollerith8("        ")
Base.String(h::Hollerith8) = String(UInt8[h.data...])
Base.show(io::IO, h::Hollerith8) = print(io, "H\"", rstrip(String(h)), "\"")

# Record 1: File identification
struct ISOTXSFileIdB
    hname::Hollerith8; huse1::Hollerith8; huse2::Hollerith8; ivers::Int32
    function ISOTXSFileIdB(huse::AbstractString="", ivers::Integer=0)
        h1 = length(huse) >= 8 ? huse[1:8] : huse
        h2 = length(huse) > 8 ? huse[9:min(end,16)] : ""
        new(Hollerith8("ISOTXS"), Hollerith8(h1), Hollerith8(h2), Int32(ivers))
    end
end

# Record 2: File control (validated)
struct ISOTXSFileControlB
    ngroup::Int32; niso::Int32; maxup::Int32; maxdn::Int32
    maxord::Int32; ichist::Int32; nscmax::Int32; nsblok::Int32
    function ISOTXSFileControlB(ng, ni, mu, md, mo, ich, nsc, nsb)
        ng > 0   || throw(ArgumentError("ngroup must be > 0, got $ng"))
        ni > 0   || throw(ArgumentError("niso must be > 0, got $ni"))
        mu >= 0  || throw(ArgumentError("maxup must be >= 0"))
        md >= 0  || throw(ArgumentError("maxdn must be >= 0"))
        mo >= 0  || throw(ArgumentError("maxord must be >= 0"))
        nsc >= 0 || throw(ArgumentError("nscmax must be >= 0"))
        nsb > 0  || throw(ArgumentError("nsblok must be > 0, got $nsb"))
        new(Int32(ng), Int32(ni), Int32(mu), Int32(md),
            Int32(mo), Int32(ich), Int32(nsc), Int32(nsb))
    end
end

# Per-isotope control and group-independent data
struct IsotopeControlB
    hisonm::Hollerith8; habsid::Hollerith8; hident::Hollerith8
    amass::Float32; efiss::Float32; ecapt::Float32
    temp::Float32; sigpot::Float32; adens::Float32
    kbr::Int32; ichi::Int32; ifis::Int32; ialf::Int32; inp::Int32
    in2n::Int32; ind::Int32; int_::Int32
    ltot::Int32; ltrn::Int32; istrpd::Int32
end

# Principal (1-D) cross sections (validated: non-empty vectors same length)
struct PrincipalXSB
    transport::Vector{Float32}; total::Vector{Float32}; ngamma::Vector{Float32}
    fission::Vector{Float32}; nubar::Vector{Float32}; chi::Vector{Float32}
    function PrincipalXSB(tr, tot, ng, fis, nu, ch)
        lens = [length(v) for v in (tr,tot,ng,fis,nu,ch) if !isempty(v)]
        if !isempty(lens) && !all(==(lens[1]), lens)
            throw(ArgumentError("non-empty XS vectors must have equal length; got $lens"))
        end
        new(Float32.(tr), Float32.(tot), Float32.(ng),
            Float32.(fis), Float32.(nu), Float32.(ch))
    end
end

# Scattering sub-block
struct ScatterSubBlockB
    idsct::Int32; lord::Int32
    jband::Vector{Int32}; ijj::Vector{Int32}; data::Vector{Float32}
end

# Complete isotope data
struct IsotopeDataB
    control::IsotopeControlB; principal::PrincipalXSB
    scatter::Vector{ScatterSubBlockB}
end

# Complete ISOTXS file (validated on construction)
struct ISOTXSFileB
    file_id::ISOTXSFileIdB; file_ctrl::ISOTXSFileControlB
    hsetid::Vector{Hollerith8}; emax::Vector{Float32}
    emin::Float32; vel::Vector{Float32}; isotopes::Vector{IsotopeDataB}
    function ISOTXSFileB(fid, fc, hsetid, emax, emin, vel, isotopes)
        ng, ni = Int(fc.ngroup), Int(fc.niso)
        length(emax) == ng || throw(ArgumentError("emax length != ngroup"))
        length(vel) == ng || throw(ArgumentError("vel length != ngroup"))
        length(isotopes) == ni || throw(ArgumentError("isotope count != niso"))
        length(hsetid) <= 12 || throw(ArgumentError("hsetid has >12 entries"))
        for i in 1:ng-1
            emax[i] >= emax[i+1] || throw(ArgumentError(
                "emax must descend: emax[$i]=$(emax[i]) < emax[$(i+1)]"))
        end
        ng < 1 || emax[ng] > emin || throw(ArgumentError("emax[end] must exceed emin"))
        all(>(0), vel) || throw(ArgumentError("all velocities must be positive"))
        new(fid, fc, hsetid, Float32.(emax), Float32(emin), Float32.(vel), isotopes)
    end
end

# Binary record helpers
function _write_cccc_rec_b(io::IO, buf::Vector{UInt8})
    n = Int32(length(buf)); write(io, n); write(io, buf); write(io, n)
end
_emit_b!(buf, x::Int32)       = append!(buf, reinterpret(UInt8, [x]))
_emit_b!(buf, x::Float32)     = append!(buf, reinterpret(UInt8, [x]))
_emit_b!(buf, h::Hollerith8)  = append!(buf, UInt8[h.data...])
_emit_b!(buf, v::AbstractVector{Int32})   = for x in v; _emit_b!(buf, x); end
_emit_b!(buf, v::AbstractVector{Float32}) = for x in v; _emit_b!(buf, x); end

"""
    write_isotxs_b(io::IO, f::ISOTXSFileB)

Write ISOTXS in CCCC-IV binary format with Fortran record markers.
"""
function write_isotxs_b(io::IO, f::ISOTXSFileB)
    fc = f.file_ctrl; ng = Int(fc.ngroup); ni = Int(fc.niso)
    # Record 1: File identification
    buf = UInt8[]
    _emit_b!(buf, f.file_id.hname); _emit_b!(buf, f.file_id.huse1)
    _emit_b!(buf, f.file_id.huse2); _emit_b!(buf, f.file_id.ivers)
    _write_cccc_rec_b(io, buf)
    # Record 2: File control
    buf = UInt8[]
    for v in (fc.ngroup,fc.niso,fc.maxup,fc.maxdn,fc.maxord,fc.ichist,fc.nscmax,fc.nsblok)
        _emit_b!(buf, v)
    end
    _write_cccc_rec_b(io, buf)
    # Record 3: File data
    buf = UInt8[]
    for i in 1:12; _emit_b!(buf, i<=length(f.hsetid) ? f.hsetid[i] : Hollerith8()); end
    for iso in f.isotopes; _emit_b!(buf, iso.control.hisonm); end
    _emit_b!(buf, f.vel); _emit_b!(buf, f.emax); _emit_b!(buf, f.emin)
    for i in 1:ni; _emit_b!(buf, Int32(i)); end
    _write_cccc_rec_b(io, buf)
    # Per-isotope records
    for iso in f.isotopes; _write_isotope_b(io, iso, fc); end
    return nothing
end

function _write_isotope_b(io::IO, iso::IsotopeDataB, fc::ISOTXSFileControlB)
    ctrl = iso.control; ng = Int(fc.ngroup); nsc = Int(fc.nscmax)
    buf = UInt8[]
    _emit_b!(buf, ctrl.hisonm); _emit_b!(buf, ctrl.habsid); _emit_b!(buf, ctrl.hident)
    for v in (ctrl.amass,ctrl.efiss,ctrl.ecapt,ctrl.temp,ctrl.sigpot,ctrl.adens)
        _emit_b!(buf, v)
    end
    for v in (ctrl.kbr,ctrl.ichi,ctrl.ifis,ctrl.ialf,ctrl.inp,
              ctrl.in2n,ctrl.ind,ctrl.int_,ctrl.ltot,ctrl.ltrn,ctrl.istrpd)
        _emit_b!(buf, v)
    end
    for blk in iso.scatter; _emit_b!(buf, blk.idsct); end
    for _ in length(iso.scatter)+1:nsc; _emit_b!(buf, Int32(0)); end
    for blk in iso.scatter; _emit_b!(buf, blk.lord); end
    for _ in length(iso.scatter)+1:nsc; _emit_b!(buf, Int32(0)); end
    for blk in iso.scatter; _emit_b!(buf, blk.jband); end
    for blk in iso.scatter; _emit_b!(buf, blk.ijj); end
    _write_cccc_rec_b(io, buf)
    # Principal XS
    p = iso.principal; buf = UInt8[]
    for g in 1:ng
        _emit_b!(buf, p.ngamma[g])
        if !isempty(p.transport); _emit_b!(buf, p.transport[g]); end
        if !isempty(p.total); _emit_b!(buf, p.total[g]); end
        if !isempty(p.fission); _emit_b!(buf, p.fission[g]); end
        if !isempty(p.nubar); _emit_b!(buf, p.nubar[g]); end
        if !isempty(p.chi); _emit_b!(buf, p.chi[g]); end
    end
    isempty(buf) || _write_cccc_rec_b(io, buf)
    # Scatter sub-blocks
    for blk in iso.scatter
        isempty(blk.data) && continue
        buf = UInt8[]; _emit_b!(buf, blk.data); _write_cccc_rec_b(io, buf)
    end
end

"""
    build_isotxs_b(mg::MultiGroupXS; kwargs...) -> ISOTXSFileB

Build validated ISOTXS from multigroup cross section data.
"""
function build_isotxs_b(mg::MultiGroupXS;
        huse::AbstractString="", hsetid::AbstractVector{<:AbstractString}=String[],
        ivers::Integer=0, name::AbstractString="iso1",
        amass::Real=1.0, vel::Union{Nothing,AbstractVector{<:Real}}=nothing)
    gb = mg.group_bounds; ng = length(gb) - 1
    emax_desc = Float32[Float32(gb[ng-g+2]) for g in 1:ng]
    emin_val = Float32(gb[1])
    if vel === nothing
        vel_vec = Float32[Float32(sqrt(2*0.5*(gb[g]+gb[g+1])*1.602e-19/1.675e-27))
                          for g in ng:-1:1]
    else vel_vec = Float32.(vel) end
    mt_idx = Dict(mt=>c for (c,mt) in enumerate(mg.mt_list))
    has_fis, has_cap, has_el, has_tot = haskey(mt_idx,18), haskey(mt_idx,102),
                                        haskey(mt_idx,2), haskey(mt_idx,1)
    _gv(mt) = Float32[Float32(mg.xs[g,mt_idx[mt]]) for g in 1:ng]
    transport = has_tot ? _gv(1) : Float32[]
    total     = has_tot ? _gv(1) : Float32[]
    ngamma    = has_cap ? _gv(102) : zeros(Float32, ng)
    fission   = has_fis ? _gv(18) : Float32[]
    nubar     = has_fis ? fill(Float32(2.5), ng) : Float32[]
    nscmax = has_el ? 1 : 0
    scatter = ScatterSubBlockB[]
    if has_el
        push!(scatter, ScatterSubBlockB(Int32(2), Int32(1),
            ones(Int32, ng), Int32.(1:ng), _gv(2)))
    end
    ctrl = IsotopeControlB(Hollerith8(name), Hollerith8(), Hollerith8("endf"),
        Float32(amass), Float32(0), Float32(0), Float32(300), Float32(0), Float32(0),
        Int32(0), Int32(0), Int32(has_fis), Int32(has_cap), Int32(0),
        Int32(0), Int32(0), Int32(0), Int32(!isempty(total)), Int32(!isempty(transport)), Int32(0))
    pxs = PrincipalXSB(transport, total, ngamma, fission, nubar, Float32[])
    iso = IsotopeDataB(ctrl, pxs, scatter)
    ISOTXSFileB(ISOTXSFileIdB(huse, ivers), ISOTXSFileControlB(ng,1,0,ng-1,0,0,nscmax,1),
                [Hollerith8(s) for s in hsetid], emax_desc, emin_val, vel_vec, [iso])
end
