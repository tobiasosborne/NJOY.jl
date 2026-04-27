# ACE table builders: construct ACE data from PENDF/pointwise data

using Printf

"""
    build_xss(table::ACENeutronTable) -> (nxs, jxs, xss, is_int)

Serialize an `ACENeutronTable` into flat NXS(16), JXS(32), XSS arrays.
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

    # ESZ block
    jxs[JXS_ESZ] = 1
    push_reals!(table.energy_grid)
    push_reals!(table.total_xs)
    push_reals!(table.absorption_xs)
    push_reals!(table.elastic_xs)
    push_reals!(table.heating_numbers)

    # NU block (not yet populated)
    jxs[JXS_NU] = 0

    # Reaction data blocks. Fortran NJOY aceout always sets JXS[3..7] to a
    # valid position — when NTR=0 the tables are zero-length but the pointers
    # still mark the position (end of ESZ + 1). Only the JXS[2] (NU/fission)
    # pointer is 0 for truly-absent block types. Ref: acefc.f90:5367-5371.
    mtr_pos = length(xss) + 1
    jxs[JXS_MTR]  = mtr_pos
    jxs[JXS_LQR]  = mtr_pos
    jxs[JXS_TYR]  = mtr_pos
    jxs[JXS_LSIG] = mtr_pos
    jxs[JXS_SIG]  = mtr_pos
    if n_react > 0
        jxs[JXS_MTR] = length(xss) + 1
        for r in table.reactions; push_int!(r.mt); end

        jxs[JXS_LQR] = length(xss) + 1
        for r in table.reactions; push_real!(r.q_value); end

        jxs[JXS_TYR] = length(xss) + 1
        for r in table.reactions; push_int!(r.ty); end

        # LSIG: locators
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

    # LAND block
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

    # AND block
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

    # LDLW/DLW — energy-distribution blocks. These are always "structurally
    # present" in ACE fast tables; when NR=0 (no reactions with secondary
    # neutrons) they're zero-length and Fortran points them at the
    # end-of-XSS+1 sentinel (= length+1). Ref: acefc.f90:5891-5893.
    ldlw_sentinel = length(xss) + 1
    jxs[JXS_LDLW] = ldlw_sentinel
    jxs[JXS_DLW]  = ldlw_sentinel

    # Conditional blocks — 0 means "block type not present in this file"
    jxs[JXS_GPD]   = 0
    for j in JXS_MTRP:JXS_YP; jxs[j] = 0; end
    jxs[JXS_FIS]   = 0
    jxs[JXS_IURPT] = 0
    jxs[JXS_NUD]   = 0
    jxs[JXS_DNDAT] = 0
    jxs[JXS_LDND]  = 0
    jxs[JXS_DND]   = 0

    # JXS[22] (END) = last word of basic data, i.e. length(xss). Fortran
    # aceout assigns `end=next-1` after filling XSS (acefc.f90:5893).
    jxs[JXS_END] = length(xss)

    # Fill NXS
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

"""
    build_ace_from_pendf(pendf; suffix="80c", temp_kelvin=300.0, comment="", mat_id=0)

Construct an ACENeutronTable from a PointwiseMaterial (eV -> MeV conversion).
"""
function build_ace_from_pendf(pendf::PointwiseMaterial;
                               suffix::AbstractString = "80c",
                               temp_kelvin::Float64 = 300.0,
                               comment::AbstractString = "",
                               mat_id::Integer = 0,
                               za::Integer = 0,
                               awr::Float64 = NaN,
                               date::AbstractString = "",
                               charged_elastic::Union{Nothing, NamedTuple} = nothing)
    # Prefer explicit ZA (from MF1/MT451) over pendf.mat (which is MAT, not ZA).
    # Same for AWR — pendf often has mass number A in place of real AWR.
    true_za = za > 0 ? Int(za) : Int(pendf.mat)
    actual_mat = mat_id > 0 ? mat_id : Int(pendf.mat)
    zaid_str = format_zaid(true_za, suffix)
    temp_mev = temp_to_mev(temp_kelvin)
    mat_str = @sprintf("   mat%4d", actual_mat)
    header_awr = isnan(awr) ? Float64(true_za % 1000) : awr

    header = ACEHeader(zaid=zaid_str, awr=header_awr,
                       temp_mev=temp_mev, comment=comment,
                       mat_string=mat_str, date=date)

    n = length(pendf.energies)
    energy_mev = pendf.energies .* 1e-6

    mt_to_col = Dict(mt => i for (i, mt) in enumerate(pendf.mt_list))
    col_total   = get(mt_to_col, 1, 0)
    col_elastic = get(mt_to_col, 2, 0)

    elastic_xs = col_elastic > 0 ? pendf.cross_sections[:, col_elastic] : zeros(n)
    # Total XS: prefer MT=1 if present; otherwise reconstruct as sum of all
    # non-redundant reaction MTs (matches Fortran acer when the PENDF has no
    # explicit MT=1 — e.g. charged-particle evaluations with only MT=2).
    # Ref: njoy-reference/src/acefc.f90 unionx() lines 1538-1702.
    total_xs = if col_total > 0
        pendf.cross_sections[:, col_total]
    else
        s = zeros(n)
        for (col, mt) in enumerate(pendf.mt_list)
            # skip redundant sums if present (3, 4, 103-107)
            mt in (3, 4, 103, 104, 105, 106, 107) && continue
            # skip MT=251-300 (derived quantities like mu_bar, xi)
            251 <= mt <= 300 && continue
            s .+= pendf.cross_sections[:, col]
        end
        s
    end
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

    # Charged-particle elastic override: replace total/elastic/heating with
    # Coulomb-corrected values from acer_charged_elastic, attach the LAW=14
    # angular block. Note: Fortran acecpe modifies xss(esz+nes..2*nes-1)
    # (total) and xss(esz+3*nes..4*nes-1) (elastic) but leaves the
    # absorption column xss(esz+2*nes..3*nes-1) untouched — see
    # acefc.f90:6646-6667. We therefore preserve `absorption_xs` from the
    # pre-Coulomb computation rather than recomputing it as total-elastic
    # (which would introduce ~1e-14 FP noise from the asymmetric 9/7-sigfig
    # rounding of total vs elastic).
    angular_elastic = nothing
    if charged_elastic !== nothing
        angular_elastic = charged_elastic.angular
        total_xs    = charged_elastic.total
        elastic_xs  = charged_elastic.elastic
        heating     = charged_elastic.heating
    end

    ACENeutronTable(header=header,
                    energy_grid=energy_mev, total_xs=total_xs,
                    absorption_xs=absorption_xs, elastic_xs=elastic_xs,
                    heating_numbers=heating, reactions=reactions,
                    angular_elastic=angular_elastic)
end

"""
    build_ace(pendf; suffix="80c", awr=0.0, temperature=300.0, comment="", date="")

Construct a flat `ACETable` from pointwise cross section data.
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
    if awr <= 0.0; awr = Float64(ia_val); end
    bk_ev = 8.617333262e-5
    tz = temperature * bk_ev / 1.0e6
    zaid = @sprintf("%d", izaid) * "." * suffix
    if isempty(date); date = Dates.format(Dates.today(), "mm/dd/yyyy"); end
    if isempty(comment); comment = "NJOY.jl ACE output"; end
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
        for mt in extra_mts; push!(xss, mt == 102 ? 0.0 : mt == 18 ? 19.0 : 1.0); end
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
