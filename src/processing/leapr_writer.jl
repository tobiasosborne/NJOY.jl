# =================================================================
# endout — MF1/MT451 + MF7/MT2 + MF7/MT4 ENDF tape writer.
# Ref: leapr.f90:2972-3623 (endout).
# =================================================================

"Fortran a11-style 11-char float format, via the existing ENDF formatter."
_leapr_a11(x::Real) = format_endf_float(Float64(x))

"Trailer bytes for an ENDF line: cols 67-80 = MAT/MF/MT/NS."
@inline _leapr_trailer(mat, mf, mt, ns) = @sprintf("%4d%2d%3d%5d", mat, mf, mt, ns)

"Write a CONT-shaped line: 6 floats + trailer. Returns ns+1."
function _leapr_write_cont(io::IO, c1, c2, l1::Integer, l2::Integer,
                           n1::Integer, n2::Integer,
                           mat::Integer, mf::Integer, mt::Integer, ns::Integer)
    line = _leapr_a11(c1) * _leapr_a11(c2) *
           lpad(string(Int(l1)), 11) * lpad(string(Int(l2)), 11) *
           lpad(string(Int(n1)), 11) * lpad(string(Int(n2)), 11) *
           _leapr_trailer(mat, mf, mt, ns)
    println(io, line)
    return ns + 1
end

"Write a line of raw text for MF1/MT451 hollerith data (66-char). Returns ns+1."
function _leapr_write_text_line(io::IO, text::AbstractString,
                                mat::Integer, mf::Integer, mt::Integer, ns::Integer)
    t66 = rpad(text, 66)[1:66]
    println(io, t66 * _leapr_trailer(mat, mf, mt, ns))
    return ns + 1
end

"Write a DICT entry: 22 spaces + 4 integers + trailer (Fortran dictio convention)."
function _leapr_write_dict(io::IO, mf::Integer, mt::Integer, nc::Integer, mod::Integer,
                           mat::Integer, ns::Integer)
    line = repeat(" ", 22) *
           lpad(string(Int(mf)), 11)  * lpad(string(Int(mt)), 11) *
           lpad(string(Int(nc)), 11)  * lpad(string(Int(mod)), 11) *
           _leapr_trailer(mat, 1, 451, ns)
    println(io, line)
    return ns + 1
end

"Write LIST-style data (6 values/line) returning updated ns."
function _leapr_write_data_6perline(io::IO, vals::AbstractVector{<:Real},
                                    mat::Integer, mf::Integer, mt::Integer, ns::Integer)
    n = length(vals); k = 0
    while k < n
        parts = String[]
        for j in 1:6
            if k + j <= n
                push!(parts, _leapr_a11(vals[k+j]))
            else
                push!(parts, "           ")   # blank 11 chars
            end
        end
        println(io, join(parts) * _leapr_trailer(mat, mf, mt, ns))
        ns += 1; k += 6
    end
    return ns
end

"Write (x,y) pairs 3/line returning updated ns (TAB1-style)."
function _leapr_write_pairs_3perline(io::IO, x::AbstractVector{<:Real},
                                     y::AbstractVector{<:Real},
                                     mat::Integer, mf::Integer, mt::Integer, ns::Integer)
    n = length(x); k = 0
    while k < n
        parts = String[]
        for j in 1:3
            if k + j <= n
                push!(parts, _leapr_a11(x[k+j]), _leapr_a11(y[k+j]))
            else
                push!(parts, "           ", "           ")
            end
        end
        println(io, join(parts) * _leapr_trailer(mat, mf, mt, ns))
        ns += 1; k += 3
    end
    return ns
end

"Write integer interp table (NBT, INT) pairs 3-per-line, returning ns."
function _leapr_write_interp(io::IO, nbt::Vector{Int}, ilaw::Vector{Int},
                             mat::Integer, mf::Integer, mt::Integer, ns::Integer)
    n = length(nbt); k = 0
    while k < n
        parts = String[]
        for j in 1:3
            if k + j <= n
                push!(parts, lpad(string(nbt[k+j]), 11),
                             lpad(string(ilaw[k+j]), 11))
            else
                push!(parts, "           ", "           ")
            end
        end
        println(io, join(parts) * _leapr_trailer(mat, mf, mt, ns))
        ns += 1; k += 3
    end
    return ns
end

"Blank separator record (SEND/FEND/MEND/TEND)."
function _leapr_write_sep(io::IO, mat::Integer, mf::Integer, mt::Integer, ns::Integer)
    println(io, repeat(" ", 66) * _leapr_trailer(mat, mf, mt, ns))
end

"""
Count the MF1/MT451 section NC lines (for its own DICT entry).

    4 header CONT records + NWD hollerith lines + NXC dict entries + 1 SEND = 4 + NWD + NXC + 1
"""
_leapr_mf1_nc(nwd::Int, nxc::Int) = 4 + nwd + nxc + 1

"""
Count the MF7/MT4 section NC (for DICT).

Structure: HEAD + LIST (B-coeffs, 1 header + ceil(NPL/6) data = 2 lines for NPL=6) +
TAB2 (1 header + ceil(NR*2/6)=1 interp line) + per-β blocks + SEND.

For T22 (ntempr=1, isym=1, nbeta=105, nalpha=59, NPL=6):
  Per-β block = TAB1 (1 header + ceil(NR*2/6)=1 interp + ceil(NP/3)=20 pair lines) = 22 lines
  Total β-rows = 2·nbeta - 1 = 209 for isym ∈ {1,3}, else nbeta.
"""
function _leapr_mf7_mt4_nc(p::LeaprParams, isym::Int)
    nbeta  = p.nbeta
    nalpha = p.nalpha
    ntempr = p.ntempr
    nbt_rows = (isym == 1 || isym == 3) ? 2*nbeta - 1 : nbeta
    # Per-β TAB1 (first temp) = 1 header + 1 interp + ceil(nalpha/3) pair lines
    # Per-β LIST (subsequent temps) = 1 header + ceil(nalpha/6) data lines
    first_t = 1 + 1 + cld(nalpha, 3)
    per_extra_t = 1 + cld(nalpha, 6)
    per_beta = first_t + (ntempr - 1) * per_extra_t
    # HEAD + LIST(B-coefs: 1 header + ceil(NPL/6) data)
    npl = 6 * (p.nss + 1)
    b_block = 1 + 1 + cld(npl, 6)
    tab2_block = 1 + 1                   # header + 1 interp line (NR=1)
    beta_header = cld(nbt_rows, 6)       # β-grid list in TAB2 body
    # But TAB2 itself doesn't carry the β values as a list — the β values are the
    # C2 field of each nested TAB1. So no separate β list block.
    # Tail: effective-temperature TAB1s (1 each for tempf, plus tempf1 if nss>0 && b7<=0)
    tail = 1 + 1 + cld(ntempr, 3)        # tempf TAB1 (conservative)
    return b_block + tab2_block + per_beta * nbt_rows + tail + 1  # +1 SEND
end

"""
    write_leapr_tape(path, p::LeaprParams, ssm, ssp;
                     bragg=Float64[], nedge=0,
                     tempf=[p.temperatures[1]], tempf1=Float64[],
                     dwpix=[0.0], dwp1=Float64[],
                     isym=determined-from-deck, iel=p.iel)

Emit the full ENDF leapr output tape: TPID + MF1/MT451 descriptive header
+ (optional MF7/MT2 elastic if iel≠0) + MF7/MT4 inelastic + SEND/FEND/MEND/TEND.

Array-layout contract: `ssm[β, α, T]` and `ssp[β, α, T]` (β as the fastest
axis, matching the Fortran ssm(nbeta,nalpha,ntempr)).

Ref: leapr.f90:2972-3623 (endout).
"""
function write_leapr_tape(path::AbstractString, p::LeaprParams,
                          ssm::AbstractArray{Float64,3},
                          ssp::AbstractArray{Float64,3};
                          tempf::Vector{Float64}=copy(p.temperatures),
                          tempf1::Vector{Float64}=Float64[],
                          dwpix::Vector{Float64}=zeros(p.ntempr),
                          dwp1::Vector{Float64}=Float64[],
                          isym::Int=_leapr_infer_isym(p),
                          iel::Int=p.iel,
                          bragg::Vector{Float64}=Float64[],
                          nedge::Int=0)
    mat = p.mat
    open(path, "w") do io
        # ------------- TPID -------------
        # Per NJOY convention: blank text, MAT=1, MF=0, MT=0, NS=0
        _leapr_write_text_line(io, "", 1, 0, 0, 0)

        # ------------- MF1/MT451 -------------
        ns = 1
        # HEAD1: ZA, AWR, LRP=-1, 0, 0, 0
        ns = _leapr_write_cont(io, p.za, p.awr, -1, 0, 0, 0, mat, 1, 451, ns)
        # HEAD2: 0, 0, 0, 0, 0, 6
        ns = _leapr_write_cont(io, 0.0, 0.0, 0, 0, 0, 6, mat, 1, 451, ns)
        # HEAD3: 1.0, 0.0, 0, 0, 12, 6  (AWI=1, EMAX=0, LREL=0, 0, NSUB=12, NVER=6)
        ns = _leapr_write_cont(io, 1.0, 0.0, 0, 0, 12, 6, mat, 1, 451, ns)

        # HEAD4: (0, 0, 0, 0, NWD, NXC) — emit after NWD+NXC known.
        nwd = length(p.comments)
        nxc = iel == 0 ? 2 : 3
        ns = _leapr_write_cont(io, 0.0, 0.0, 0, 0, nwd, nxc, mat, 1, 451, ns)

        # Hollerith descriptive lines (66-char text)
        for txt in p.comments
            ns = _leapr_write_text_line(io, txt, mat, 1, 451, ns)
        end

        # DICT — section NC counts must match final tape.
        mf1_nc = _leapr_mf1_nc(nwd, nxc)
        mf7_mt4_nc = _leapr_mf7_mt4_nc(p, isym)
        ns = _leapr_write_dict(io, 1, 451, mf1_nc, 0, mat, ns)
        if iel != 0
            # Placeholder MT=2 NC when coher is wired
            ns = _leapr_write_dict(io, 7, 2, 1, 0, mat, ns)
        end
        ns = _leapr_write_dict(io, 7, 4, mf7_mt4_nc, 0, mat, ns)

        # SEND for MF1
        _leapr_write_sep(io, mat, 1, 0, 99999)
        # FEND for MF1 (zero-line sep with MF=0, MT=0, NS=0)
        _leapr_write_sep(io, mat, 0, 0, 0)

        # ------------- MF7/MT4 inelastic -------------
        ns = _write_mf7_mt4(io, p, ssm, ssp, isym, tempf, tempf1, dwpix, dwp1)

        # SEND for MF7
        _leapr_write_sep(io, mat, 7, 0, 99999)
        _leapr_write_sep(io, mat, 0, 0, 0)

        # MEND, TEND
        _leapr_write_sep(io, 0, 0, 0, 0)
        _leapr_write_sep(io, -1, 0, 0, 0)
    end
    return path
end

"Infer isym from the deck: 1 → cold-H (ncold>0), 0 → symmetric."
_leapr_infer_isym(p::LeaprParams) = p.ncold > 0 ? 1 : 0

"""
    _write_mf7_mt4(io, p, ssm, ssp, isym, tempf, tempf1, dwpix, dwp1) -> ns

Emit the MF=7 MT=4 inelastic section exactly per leapr.f90:3291-3577.
Returns the next ns. T22-aware (ntempr=1, isym=1, nss=0 short-cuts apply).
"""
function _write_mf7_mt4(io::IO, p::LeaprParams,
                        ssm::AbstractArray{Float64,3},
                        ssp::AbstractArray{Float64,3},
                        isym::Int,
                        tempf::Vector{Float64},
                        tempf1::Vector{Float64},
                        dwpix::Vector{Float64},
                        dwp1::Vector{Float64})
    mat = p.mat
    nbeta  = p.nbeta
    nalpha = p.nalpha
    ntempr = p.ntempr
    ns = 1

    # HEAD: ZA, AWR, 0, LAT, LASYM, 0
    ns = _leapr_write_cont(io, p.za, p.awr, 0, p.lat, isym, 0, mat, 7, 4, ns)

    # LIST: (0, 0, ilog, 0, NPL, nss) + data = [σ_b, β_max·therm, AWR, β_max·therm, 0, npr, (+nss secondary entries)]
    therm = 0.0253
    sc = p.lat == 1 ? therm / (PhysicsConstants.bk * p.temperatures[1]) : 1.0
    β_max_energy = p.beta[end] * (p.lat == 1 ? therm : therm)  # β_max·therm (eV)
    sb = p.spr * ((1 + p.awr) / p.awr)^2                  # bound scattering XS (unused but computed)
    npl = 6 * (p.nss + 1)
    list_data = Float64[p.spr * p.npr, β_max_energy, p.awr, β_max_energy, 0.0, Float64(p.npr)]
    if p.nss > 0
        sbs = p.sps * ((1 + p.aws) / p.aws)^2
        append!(list_data, [p.b7, p.sps * p.mss, p.aws, 0.0, 0.0, Float64(p.mss)])
    end
    ns = _leapr_write_cont(io, 0.0, 0.0, p.ilog, 0, npl, p.nss, mat, 7, 4, ns)
    ns = _leapr_write_data_6perline(io, list_data, mat, 7, 4, ns)

    # TAB2: (0, 0, 0, 0, NR=1, NZ=nbt_rows) + interp (nbt_rows, INT=4)
    nbt_rows = (isym == 1 || isym == 3) ? (2*nbeta - 1) : nbeta
    ns = _leapr_write_cont(io, 0.0, 0.0, 0, 0, 1, nbt_rows, mat, 7, 4, ns)
    ns = _leapr_write_interp(io, [nbt_rows], [4], mat, 7, 4, ns)

    # Per-β nested TAB1 (first temp) + LIST (subsequent temps)
    for jj in 1:nbt_rows
        # β-index + direction from isym (Fortran :3346-3369)
        β_val, k_side, i_idx = _leapr_beta_for_row(jj, nbeta, isym, p.beta)
        for nt in 1:ntempr
            tev = PhysicsConstants.bk * p.temperatures[nt]
            be_scaled = β_val * (p.lat == 1 ? therm/tev : 1.0)
            yvals = [_leapr_s_value(ssm, ssp, k_side, i_idx, j, nt, be_scaled, isym,
                                    p.ilog, p.smin)
                     for j in 1:nalpha]
            if nt == 1
                # TAB1: (T, β, LT=ntempr-1, 0, NR=1, NP=nalpha)
                ns = _leapr_write_cont(io, p.temperatures[nt], β_val,
                                       ntempr - 1, 0, 1, nalpha, mat, 7, 4, ns)
                ns = _leapr_write_interp(io, [nalpha], [4], mat, 7, 4, ns)
                ns = _leapr_write_pairs_3perline(io, p.alpha, yvals, mat, 7, 4, ns)
            else
                # LIST: (T, β, LI=1, 0, NP=nalpha, 0)
                ns = _leapr_write_cont(io, p.temperatures[nt], β_val,
                                       1, 0, nalpha, 0, mat, 7, 4, ns)
                ns = _leapr_write_data_6perline(io, yvals, mat, 7, 4, ns)
            end
        end
    end

    # Effective-temperature TAB1s (tempf always; tempf1 when nss>0 and b7<=0)
    if p.nss > 0 && p.b7 <= 0 && !isempty(tempf1)
        ns = _leapr_write_cont(io, 0.0, 0.0, 0, 0, 1, ntempr, mat, 7, 4, ns)
        ns = _leapr_write_interp(io, [ntempr], [2], mat, 7, 4, ns)
        ns = _leapr_write_pairs_3perline(io, p.temperatures, tempf1, mat, 7, 4, ns)
    end
    ns = _leapr_write_cont(io, 0.0, 0.0, 0, 0, 1, ntempr, mat, 7, 4, ns)
    ns = _leapr_write_interp(io, [ntempr], [2], mat, 7, 4, ns)
    ns = _leapr_write_pairs_3perline(io, p.temperatures, tempf, mat, 7, 4, ns)

    return ns
end

"""Given the extended-β row index `jj` ∈ [1, nbt_rows] and symmetry flag `isym`,
return `(β_val, k_side, i_idx)` where `k_side ∈ {:m, :p}` selects ssm or ssp
and `i_idx` is the β-axis index into the selected array.

Fortran ref: leapr.f90:3341-3369.
"""
function _leapr_beta_for_row(jj::Int, nbeta::Int, isym::Int, beta::Vector{Float64})
    if iseven(isym)
        return (beta[jj], :m, jj)
    elseif jj < nbeta
        i = nbeta - jj + 1
        return (-beta[i], :m, i)
    else
        i = jj - nbeta + 1
        return (beta[i], :p, i)
    end
end

"""Fetch S(α,β,T) with exp(±β/2) detailed-balance transform per isym.

Fortran ref: leapr.f90:3354-3453. isym semantics:
- 0 → symmetric S(α,|β|); emit `S · exp(-β/2)`
- 1 → asymmetric, cold-H (ssm,ssp); emit `S · exp(+β/2)` on both branches
- 2 → asymmetric, negative-β only; emit raw S
- 3 → asymmetric both branches; emit raw S

`ilog` controls log₁₀ emission (unused for T22 which has ilog=0).
"""
function _leapr_s_value(ssm::AbstractArray{Float64,3}, ssp::AbstractArray{Float64,3},
                        k_side::Symbol, i::Int, j::Int, nt::Int, be::Float64,
                        isym::Int, ilog::Int, smin::Float64)
    raw = k_side === :m ? ssm[i, j, nt] : ssp[i, j, nt]
    val = if isym == 0
        raw * exp(-be/2)
    elseif isym == 1
        raw * exp(+be/2)
    else           # isym ∈ {2, 3}
        raw
    end
    val < smin && return 0.0
    ilog == 1 && return val > 0 ? log10(val) : -999.0
    return val
end
