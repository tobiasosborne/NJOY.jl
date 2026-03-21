# ACE Type 1 (ASCII) writer for MCNP continuous-energy neutron tables
# Matches NJOY2016 aceout()/typen()/change() format exactly.
# aceout -> write_ace_table; typen -> _ace_fmt_int/_ace_fmt_real;
# acelod -> build_ace_from_pendf; change -> build_xss (int/real tracking)

using Printf

# =========================================================================
# Low-level Type 1 format helpers (matching NJOY2016 typen subroutine)
# =========================================================================

"""Format an integer for XSS data block: right-justified in 20 chars (i20)."""
_ace_fmt_int(val::Integer) = lpad(string(Int(val)), 20)

"""Format a real for XSS data block: scientific in 20 chars (1pe20.11)."""
_ace_fmt_real(val::Real) = @sprintf("%20.11E", Float64(val))

# Also keep Proposer-A naming for compatibility
_format_xss_value(x::Float64) = _ace_fmt_real(x)
_format_xss_integer(x::Int) = _ace_fmt_int(x)

# =========================================================================
# Build NXS/JXS/XSS from ACENeutronTable (Proposer-B: structured -> flat)
# =========================================================================

"""
    build_xss(table::ACENeutronTable) -> (nxs, jxs, xss, is_int)

Serialize an `ACENeutronTable` into flat NXS(16), JXS(32), XSS arrays.
Returns a boolean mask `is_int[i]` indicating whether XSS[i] is an integer
field (for correct Type 1 formatting).
"""
function build_xss(table::ACENeutronTable)
    n_es = length(table.energy_grid)
    n_react = length(table.reactions)
    nr_count = count(r -> abs(r.ty) > 0, table.reactions)

    nxs = zeros(Int32, 16)
    jxs = zeros(Int32, 32)
    xss = Float64[]
    is_int = Bool[]

    push_real!(v) = (push!(xss, Float64(v)); push!(is_int, false))
    push_reals!(vs) = for v in vs; push_real!(v); end
    push_int!(v) = (push!(xss, Float64(v)); push!(is_int, true))

    # --- ESZ block ---
    jxs[JXS_ESZ] = 1
    push_reals!(table.energy_grid)
    push_reals!(table.total_xs)
    push_reals!(table.absorption_xs)
    push_reals!(table.elastic_xs)
    push_reals!(table.heating_numbers)

    # --- NU block (not yet populated) ---
    jxs[JXS_NU] = 0

    # --- Reaction data blocks ---
    if n_react > 0
        jxs[JXS_MTR] = length(xss) + 1
        for r in table.reactions; push_int!(r.mt); end

        jxs[JXS_LQR] = length(xss) + 1
        for r in table.reactions; push_real!(r.q_value); end

        jxs[JXS_TYR] = length(xss) + 1
        for r in table.reactions; push_int!(r.ty); end

        # LSIG: locators relative to SIG start
        jxs[JXS_LSIG] = length(xss) + 1
        sig_offset = 1
        sig_offsets = Int32[]
        for r in table.reactions
            push!(sig_offsets, sig_offset)
            sig_offset += 2 + length(r.xs)
        end
        for off in sig_offsets; push_int!(off); end

        # SIG block
        jxs[JXS_SIG] = length(xss) + 1
        for r in table.reactions
            push_int!(r.ie_start)
            push_int!(length(r.xs))
            push_reals!(r.xs)
        end
    end

    # --- LAND block ---
    jxs[JXS_LAND] = length(xss) + 1
    if table.angular_elastic !== nothing
        push_int!(1)   # data present at offset 1
    else
        push_int!(0)   # isotropic
    end
    for r in table.reactions
        abs(r.ty) > 0 || continue
        push_int!(haskey(table.angular, Int(r.mt)) ? 0 : 0)
    end

    # --- AND block ---
    jxs[JXS_AND] = length(xss) + 1
    if table.angular_elastic !== nothing
        ab = table.angular_elastic
        ne_a = length(ab.energies)
        push_int!(ne_a)
        push_reals!(ab.energies)
        loc_start = length(xss) + 1
        and_base = jxs[JXS_AND]
        for _ in 1:ne_a; push_int!(0); end
        for i in 1:ne_a
            d = ab.distributions[i]
            if d isa EquiprobableBins
                rel = length(xss) + 1 - and_base + 1
                xss[loc_start + i - 1] = Float64(rel)
                push_reals!(d.cosines)
            elseif d isa TabulatedAngular
                rel = -(length(xss) + 1 - and_base + 1)
                xss[loc_start + i - 1] = Float64(rel)
                push_int!(d.iflag)
                push_int!(length(d.cosine))
                push_reals!(d.cosine)
                push_reals!(d.pdf)
                push_reals!(d.cdf)
            end
        end
    end

    # --- Remaining blocks (0 = not present) ---
    jxs[JXS_LDLW]  = 0
    jxs[JXS_DLW]   = 0
    jxs[JXS_GPD]   = 0
    for j in JXS_MTRP:JXS_YP; jxs[j] = 0; end
    jxs[JXS_FIS]   = 0
    jxs[JXS_END]   = length(xss) + 1
    jxs[JXS_IURPT] = 0
    jxs[JXS_NUD]   = 0
    jxs[JXS_DNDAT] = 0
    jxs[JXS_LDND]  = 0
    jxs[JXS_DND]   = 0

    # --- Fill NXS ---
    za, _ = parse_zaid(strip(table.header.hz))
    nxs[NXS_LEN2]  = length(xss)
    nxs[NXS_IZAID] = za
    nxs[NXS_NES]   = n_es
    nxs[NXS_NTR]   = n_react
    nxs[NXS_NR]    = nr_count
    nxs[NXS_IZ]    = za ÷ 1000
    nxs[NXS_IA]    = za % 1000

    (nxs, jxs, xss, is_int)
end

# =========================================================================
# Type 1 ACE file writer (Proposer-B: from ACENeutronTable)
# =========================================================================

"""
    write_ace_table(io::IO, table::ACENeutronTable; format=:ascii)

Write an ACENeutronTable as a Type 1 (ASCII) ACE file.
Output format matches NJOY2016 aceout() exactly for MCNP compatibility.
"""
function write_ace_table(io::IO, table::ACENeutronTable; format::Symbol=:ascii)
    format === :ascii || throw(ArgumentError("only :ascii (Type 1) supported"))
    nxs, jxs, xss, is_int = build_xss(table)
    h = table.header

    # Line 1: hz(a10) aw0(f12.6) ' ' tz(1pe11.4) ' ' hd(a10)
    @printf(io, "%s%12.6f %11.4E %s\n",
            rpad(h.hz, 10)[1:10], h.aw0, h.tz, rpad(h.hd, 10)[1:10])

    # Line 2: hk(a70) hm(a10)
    print(io, rpad(h.hk, 70)[1:70], rpad(h.hm, 10)[1:10], "\n")

    # Lines 3-6: IZ/AW pairs, 4(i7,f11.0) per line
    for row in 1:4
        for col in 1:4
            idx = (row - 1) * 4 + col
            @printf(io, "%7d%11.0f", table.iz[idx], table.aw[idx])
        end
        print(io, "\n")
    end

    # Lines 7-12: NXS(16)+JXS(32) = 48 integers, 8i9 per line
    all_ij = vcat(Int.(nxs), Int.(jxs))
    for row in 1:6
        for col in 1:8
            @printf(io, "%9d", all_ij[(row-1)*8 + col])
        end
        print(io, "\n")
    end

    # Data: XSS values, 4 per line, 20 chars each
    for i in 1:length(xss)
        if is_int[i]
            print(io, _ace_fmt_int(round(Int, xss[i])))
        else
            print(io, _ace_fmt_real(xss[i]))
        end
        if i % 4 == 0 || i == length(xss)
            print(io, "\n")
        end
    end
end

"""
    write_ace_directory(io::IO, table::ACENeutronTable, nxs;
                        itype=1, filepath="filename", route="route")

Write the MCNP xsdir directory line for an ACE table.
"""
function write_ace_directory(io::IO, table::ACENeutronTable,
                              nxs::AbstractVector{<:Integer};
                              itype::Int=1, filepath::String="filename",
                              route::String="route")
    h = table.header
    @printf(io, "%s%12.6f %s %s%2d%4d %8d%6d%6d%10.3E\n",
            rpad(strip(h.hz), 10)[1:10], h.aw0,
            filepath, route, itype, 1, nxs[NXS_LEN2], 0, 0, h.tz)
end

# write_ace dispatches to write_ace_table for ACENeutronTable
write_ace(io::IO, table::ACENeutronTable) = write_ace_table(io, table)

# =========================================================================
# Build ACENeutronTable from PointwiseMaterial (PENDF output)
# =========================================================================
"""
    build_ace_from_pendf(pendf::PointwiseMaterial;
                         suffix="80c", temp_kelvin=300.0,
                         comment="", mat_id=0) -> ACENeutronTable

Construct an ACENeutronTable from a PointwiseMaterial.
Energies are converted from eV to MeV. Cross sections remain in barns.
"""
function build_ace_from_pendf(pendf::PointwiseMaterial;
                               suffix::AbstractString = "80c",
                               temp_kelvin::Float64 = 300.0,
                               comment::AbstractString = "",
                               mat_id::Integer = 0)
    za = Int(pendf.mat)
    actual_mat = mat_id > 0 ? mat_id : za
    zaid_str = format_zaid(za, suffix)
    temp_mev = temp_to_mev(temp_kelvin)
    mat_str = @sprintf("   mat%4d", actual_mat)

    header = ACEHeader(zaid=zaid_str, awr=Float64(za % 1000),
                       temp_mev=temp_mev, comment=comment,
                       mat_string=mat_str)

    n = length(pendf.energies)
    energy_mev = pendf.energies .* 1e-6

    mt_to_col = Dict(mt => i for (i, mt) in enumerate(pendf.mt_list))
    col_total   = get(mt_to_col, 1, 0)
    col_elastic = get(mt_to_col, 2, 0)

    total_xs = col_total > 0 ? pendf.cross_sections[:, col_total] : zeros(n)
    elastic_xs = col_elastic > 0 ? pendf.cross_sections[:, col_elastic] : zeros(n)
    absorption_xs = max.(total_xs .- elastic_xs, 0.0)
    heating = zeros(n)

    reactions = ReactionXS[]
    for (col, mt) in enumerate(pendf.mt_list)
        mt in (1, 2) && continue
        xs_col = pendf.cross_sections[:, col]
        ie_start = findfirst(x -> x > 0.0, xs_col)
        ie_start === nothing && continue
        push!(reactions, ReactionXS(
            Int32(mt), 0.0, Int32(0), Int32(ie_start),
            xs_col[ie_start:end]))
    end

    ACENeutronTable(header=header,
                    energy_grid=energy_mev, total_xs=total_xs,
                    absorption_xs=absorption_xs, elastic_xs=elastic_xs,
                    heating_numbers=heating, reactions=reactions)
end

# =========================================================================
# Write flat ACETable (NXS/JXS/XSS representation)
# =========================================================================
"""
    write_ace(io::IO, table::ACETable)

Write a flat ACETable (NXS/JXS/XSS) as a Type 1 (ASCII) ACE file.
Output format matches NJOY2016 aceout() for itype=1.
"""
function write_ace(io::IO, table::ACETable)
    # Line 1: ZAID(a10), AWR(f12.6), ' ', TZ(e11.4), ' ', date(a10)
    @printf(io, "%s%12.6f %11.4E %s\n",
            rpad(table.zaid, 10)[1:10], table.awr, table.temp,
            rpad(table.date, 10)[1:10])

    # Line 2: comment(a70), mat_id(a10)
    print(io, rpad(table.comment, 70)[1:70], rpad(table.mat_id, 10)[1:10], "\n")

    # Lines 3-6: IZ/AW pairs, 4(i7,f11.0) per line
    for row in 1:4
        for col in 1:4
            idx = (row - 1) * 4 + col
            @printf(io, "%7d%11.0f", table.pairs_iz[idx], table.pairs_aw[idx])
        end
        print(io, "\n")
    end

    # Lines 7-12: NXS(16)+JXS(32) = 48 integers, 8i9 per line
    all_ints = Int32[table.nxs..., table.jxs...]
    for row in 1:6
        for col in 1:8
            @printf(io, "%9d", all_ints[(row - 1) * 8 + col])
        end
        print(io, "\n")
    end

    # Data: XSS values, 4 per line, 20 chars each (1pe20.11)
    n = length(table.xss)
    for i in 1:n
        print(io, _ace_fmt_real(table.xss[i]))
        if i % 4 == 0 || i == n
            print(io, "\n")
        end
    end
    return nothing
end

# =========================================================================
# build_ace: flexible entry point for PointwiseMaterial -> ACETable (flat)
# =========================================================================
"""
    build_ace(pendf::PointwiseMaterial;
              suffix="80c", awr=0.0, temperature=300.0,
              comment="", date="")

Construct a flat `ACETable` from pointwise cross section data.
This is the Proposer-A interface that produces an ACETable with
NXS/JXS/XSS arrays ready for `write_ace`.
"""
function build_ace(pendf::PointwiseMaterial;
                   suffix::String = "80c",
                   awr::Float64 = 0.0,
                   temperature::Float64 = 300.0,
                   comment::String = "",
                   date::String = "")
    mat = Int(pendf.mat)
    n_es = length(pendf.energies)
    izaid = mat
    iz_val = div(izaid, 1000)
    ia_val = mod(izaid, 1000)

    if awr <= 0.0
        awr = Float64(ia_val)
    end

    bk_ev = 8.617333262e-5
    tz = temperature * bk_ev / 1.0e6

    zaid = @sprintf("%d", izaid) * "." * suffix
    if isempty(date)
        date = Dates.format(Dates.today(), "mm/dd/yyyy")
    end
    if isempty(comment)
        comment = "NJOY.jl ACE output"
    end
    mat_id_str = @sprintf("   mat%4d", mat)

    mt_list = pendf.mt_list
    col_total   = findfirst(==(1), mt_list)
    col_elastic = findfirst(==(2), mt_list)
    extra_cols = Int[]
    extra_mts  = Int[]
    for (j, mt) in enumerate(mt_list)
        if mt != 1 && mt != 2
            push!(extra_cols, j)
            push!(extra_mts, mt)
        end
    end
    ntr = length(extra_mts)

    xss = Float64[]

    # ESZ block
    for i in 1:n_es; push!(xss, pendf.energies[i] * 1.0e-6); end
    if col_total !== nothing
        for i in 1:n_es; push!(xss, pendf.cross_sections[i, col_total]); end
    else
        for i in 1:n_es
            s = sum(pendf.cross_sections[i, j] for j in 1:length(mt_list))
            push!(xss, s)
        end
    end
    for i in 1:n_es
        total_v = xss[n_es + i]
        elast_v = col_elastic !== nothing ? pendf.cross_sections[i, col_elastic] : 0.0
        push!(xss, max(0.0, total_v - elast_v))
    end
    if col_elastic !== nothing
        for i in 1:n_es; push!(xss, pendf.cross_sections[i, col_elastic]); end
    else
        append!(xss, zeros(n_es))
    end
    append!(xss, zeros(n_es))  # heating

    # Reaction blocks
    mtr_loc = Int32(0); lqr_loc = Int32(0); tyr_loc = Int32(0)
    lsig_loc = Int32(0); sig_loc = Int32(0)
    if ntr > 0
        mtr_loc = Int32(length(xss) + 1)
        for mt in extra_mts; push!(xss, Float64(mt)); end

        lqr_loc = Int32(length(xss) + 1)
        append!(xss, zeros(ntr))

        tyr_loc = Int32(length(xss) + 1)
        for mt in extra_mts
            push!(xss, mt == 102 ? 0.0 : mt == 18 ? 19.0 : 1.0)
        end

        lsig_loc = Int32(length(xss) + 1)
        sig_loc  = Int32(length(xss) + ntr + 1)
        sig_data = Float64[]
        for (_, col) in enumerate(extra_cols)
            push!(xss, Float64(length(sig_data) + 1))
            ie_s = 1
            for i in 1:n_es
                if pendf.cross_sections[i, col] > 0.0; ie_s = i; break; end
            end
            ne_r = n_es - ie_s + 1
            push!(sig_data, Float64(ie_s))
            push!(sig_data, Float64(ne_r))
            for i in ie_s:n_es; push!(sig_data, pendf.cross_sections[i, col]); end
        end
        append!(xss, sig_data)
    end

    land_loc = Int32(length(xss) + 1)
    push!(xss, -1.0)
    and_loc = Int32(length(xss) + 1)

    len2 = Int32(length(xss))

    nxs = ntuple(i -> Int32(0), 16)
    nxs = Base.setindex(nxs, len2, NXS_LEN2)
    nxs = Base.setindex(nxs, Int32(izaid), NXS_IZAID)
    nxs = Base.setindex(nxs, Int32(n_es), NXS_NES)
    nxs = Base.setindex(nxs, Int32(ntr), NXS_NTR)
    nxs = Base.setindex(nxs, Int32(0), NXS_NR)
    nxs = Base.setindex(nxs, Int32(iz_val), NXS_IZ)
    nxs = Base.setindex(nxs, Int32(ia_val), NXS_IA)

    jxs = ntuple(i -> Int32(0), 32)
    jxs = Base.setindex(jxs, Int32(1), JXS_ESZ)
    jxs = Base.setindex(jxs, mtr_loc, JXS_MTR)
    jxs = Base.setindex(jxs, lqr_loc, JXS_LQR)
    jxs = Base.setindex(jxs, tyr_loc, JXS_TYR)
    jxs = Base.setindex(jxs, lsig_loc, JXS_LSIG)
    jxs = Base.setindex(jxs, sig_loc, JXS_SIG)
    jxs = Base.setindex(jxs, land_loc, JXS_LAND)
    jxs = Base.setindex(jxs, and_loc, JXS_AND)

    pairs_iz = ntuple(i -> Int32(0), 16)
    pairs_aw = ntuple(i -> 0.0, 16)

    ACETable(zaid, awr, tz, rpad(date, 10), comment, mat_id_str,
             pairs_iz, pairs_aw, nxs, jxs, xss)
end
