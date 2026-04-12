# ERRORR -- Covariance processing for multigroup cross sections from ENDF MF33.
# Correspondence: covcal->expand_covariance_block, covout->multigroup_covariance,
# gridd->_build_union_grid, errorr->process_covariance

"One NI-type covariance sub-subsection from ENDF MF33."
struct CovarianceBlock
    mt1::Int; mt2::Int; lb::Int; lt::Int
    energies::Vector{Float64}
    data::Vector{Float64}
end

"Symmetric multigroup covariance matrix for a reaction pair."
struct CovarianceMatrix
    mt1::Int; mt2::Int
    energy_groups::Vector{Float64}
    matrix::Matrix{Float64}    # (ng x ng), relative covariance
    is_relative::Bool
end

Base.show(io::IO, cm::CovarianceMatrix) = print(io, "CovarianceMatrix(MT",
    cm.mt1, "/MT", cm.mt2, ", ", size(cm.matrix,1), " groups, ",
    cm.is_relative ? "relative" : "absolute", ")")

"""Parsed MF33 covariance section for one reaction pair."""
struct CovarianceData
    mat::Int; mt::Int; mt1::Int
    nc_subsections::Vector{NamedTuple{(:coeffs,:mts),Tuple{Vector{Float64},Vector{Int}}}}
    ni_subsections::Vector{CovarianceBlock}
end

# Block expansion: LB=0/1/2 (diagonal), LB=5 (matrix), LB=6 (asymmetric)

"Expand LB=0/1: diagonal covariance on energy intervals."
function expand_lb1(block::CovarianceBlock, egrid::AbstractVector{<:Real})
    ne = length(egrid) - 1; C = zeros(Float64, ne, ne)
    ek, fk = block.energies, block.data
    nk = min(length(fk), length(ek)) - 1
    for k in 1:nk
        elo, ehi = ek[k], ek[k+1]
        for i in 1:ne
            emid = 0.5 * (egrid[i] + egrid[i+1])
            emid >= elo && emid < ehi && (C[i,i] = fk[k])
        end
    end; C
end

"""Expand LB=2: fully correlated fractional variance, C[i,j] = F_k * F_l where
intervals are looked up by the group's LOWER bound (Fortran `un(jg)`/`un(jh)`
are lower bounds of union intervals; see errorr.f90:2087-2090)."""
function expand_lb2(block::CovarianceBlock, egrid::AbstractVector{<:Real})
    ne = length(egrid) - 1; C = zeros(Float64, ne, ne)
    ek, fk = block.energies, block.data
    nk = min(length(fk), length(ek)) - 1
    f_lo = zeros(Float64, ne)
    for i in 1:ne
        e_lo = egrid[i]
        for k in 1:nk
            if e_lo >= ek[k] && e_lo < ek[k+1]
                f_lo[i] = fk[k]; break
            end
        end
    end
    for i in 1:ne, j in 1:ne
        C[i,j] = f_lo[i] * f_lo[j]
    end
    C
end

"Expand LB=5, LT=1: symmetric upper-triangular covariance."
function expand_lb5_symmetric(block::CovarianceBlock, egrid::AbstractVector{<:Real})
    ne = length(egrid) - 1; ek = block.energies; nk = length(ek) - 1
    C = zeros(Float64, ne, ne)
    for i in 1:ne
        ki = _find_bin(0.5*(egrid[i]+egrid[i+1]), ek); ki == 0 && continue
        for j in i:ne
            kj = _find_bin(0.5*(egrid[j]+egrid[j+1]), ek); kj == 0 && continue
            k_lo, k_hi = minmax(ki, kj)
            idx = _upper_tri_index(k_lo, k_hi, nk)
            1 <= idx <= length(block.data) && (C[i,j] = C[j,i] = block.data[idx])
        end
    end; C
end

"Expand LB=5/6 full rectangular covariance."
function expand_lb5_full(block::CovarianceBlock, egrid::AbstractVector{<:Real})
    ne = length(egrid) - 1
    C = zeros(Float64, ne, ne)
    if block.lb == 6 && block.lt > 0
        # LB=6 asymmetric: energies stores [row_grid..., col_grid...]
        # row grid has NE energies (block header ne), col grid has LT energies
        nek = length(block.energies) - block.lt  # NE (row grid count)
        nel = block.lt                            # NE' (col grid count)
        ek_row = block.energies[1:nek]
        ek_col = block.energies[nek+1:nek+nel]
        nk_row = nek - 1
        nk_col = nel - 1
        for i in 1:ne
            ki = _find_bin(0.5*(egrid[i]+egrid[i+1]), ek_row); ki == 0 && continue
            for j in 1:ne
                kj = _find_bin(0.5*(egrid[j]+egrid[j+1]), ek_col); kj == 0 && continue
                idx = (ki-1)*nk_col + kj
                1 <= idx <= length(block.data) && (C[i,j] = block.data[idx])
            end
        end
    else
        # LB=5 full or LB=6 with identical grids: single shared energy grid
        ek = block.energies; nk = length(ek) - 1
        for i in 1:ne
            ki = _find_bin(0.5*(egrid[i]+egrid[i+1]), ek); ki == 0 && continue
            for j in 1:ne
                kj = _find_bin(0.5*(egrid[j]+egrid[j+1]), ek); kj == 0 && continue
                idx = (ki-1)*nk + kj
                1 <= idx <= length(block.data) && (C[i,j] = block.data[idx])
            end
        end
    end
    C
end

"""Dispatch on LB flag to expand a CovarianceBlock into a matrix on `egrid`."""
function expand_covariance_block(block::CovarianceBlock, egrid::AbstractVector{<:Real})
    block.lb in (0,1) && return expand_lb1(block, egrid)
    block.lb == 2 && return expand_lb2(block, egrid)
    block.lb == 5 && return block.lt == 1 ? expand_lb5_symmetric(block, egrid) :
                                            expand_lb5_full(block, egrid)
    block.lb == 6 && return expand_lb5_full(block, egrid)
    error("Unsupported LB=$(block.lb). Supported: 0,1,2,5,6.")
end

"""Collapse NI-type CovarianceBlocks into a multigroup covariance matrix.

When `is_relative` is not explicitly provided, it is inferred from the LB flags:
LB=0 indicates absolute covariance, so if all blocks have LB=0 the result is
marked absolute. Otherwise the result is marked relative (the common case for
LB=1..6). An optional `flux` vector (one value per fine-grid bin) enables
flux-weighted collapse matching the Fortran ERRORR module."""
function multigroup_covariance(blocks::AbstractVector{CovarianceBlock},
                                group_bounds::AbstractVector{<:Real};
                                is_relative::Union{Bool,Nothing}=nothing,
                                flux::Union{Nothing,AbstractVector{<:Real}}=nothing)
    isempty(blocks) && error("no covariance blocks provided")
    ng = length(group_bounds) - 1
    mt1, mt2 = blocks[1].mt1, blocks[1].mt2
    gb = collect(Float64, group_bounds)

    # Determine is_relative from LB flags when not explicitly provided.
    # LB=0 is absolute covariance per the ENDF manual; all other LB values
    # are relative.
    if is_relative === nothing
        is_relative = !all(b -> b.lb == 0, blocks)
    end

    # Check total fine grid size to avoid OOM
    total_e = sum(length(b.energies) for b in blocks) + ng + 1
    if total_e > 10000
        # Direct expansion onto group grid (no fine-grid intermediate)
        C_group = zeros(Float64, ng, ng)
        for block in blocks; C_group .+= expand_covariance_block(block, gb); end
    else
        ugrid = _build_union_grid(blocks, gb)
        nu = length(ugrid) - 1
        C_fine = zeros(Float64, nu, nu)
        for block in blocks; C_fine .+= expand_covariance_block(block, ugrid); end
        T = _collapse_matrix(ugrid, gb; flux=flux)
        C_group = T * C_fine * T'
    end
    C_group = 0.5 * (C_group + C_group')
    CovarianceMatrix(mt1, mt2, gb, C_group, is_relative)
end

"""Read one MF33 covariance section. Seeks to MF=33/MT=mt from file start."""
function read_mf33(io::IO, mat::Integer, mt::Integer)
    seekstart(io)
    find_section(io, 33, mt) || error("read_mf33: MF33/MT$mt not found")
    head = read_cont(io); nl = Int(head.N2)
    nc_subs = NamedTuple{(:coeffs,:mts),Tuple{Vector{Float64},Vector{Int}}}[]
    ni_subs = CovarianceBlock[]
    for _ in 1:nl
        sh = read_cont(io); mt1 = Int(sh.L2); mt1 == 0 && (mt1 = mt)
        nc = Int(sh.N1); ni = Int(sh.N2)
        for _ in 1:nc
            lst = read_list(io); nci = Int(lst.N1)
            coeffs = Float64[]; mts_nc = Int[]
            k = 1
            while k + 1 <= length(lst.data) && k <= nci
                push!(mts_nc, round(Int, lst.data[k]))
                push!(coeffs, lst.data[k+1])
                k += 2
            end
            push!(nc_subs, (coeffs=coeffs, mts=mts_nc))
        end
        for _ in 1:ni
            lst = read_list(io)
            lt = Int(lst.L1); lb = Int(lst.L2)
            np = Int(lst.N1)  # total data items in LIST
            ne = Int(lst.N2)  # NE (number of energies for LB=5,6)
            if lb in (0,1,2,3,4)
                # (Ek, Fk) pairs: NP items total = 2*NK where NK = num intervals
                nk = div(np, 2)
                push!(ni_subs, CovarianceBlock(mt, mt1, lb, lt,
                    lst.data[1:2:2nk], lst.data[2:2:2nk]))
            elseif lb == 5
                # NE energies followed by covariance matrix values
                ek = lst.data[1:ne]
                fvals = lst.data[ne+1:np]
                push!(ni_subs, CovarianceBlock(mt, mt1, lb, lt, ek, fvals))
            elseif lb == 6
                # LB=6: NE energies for rows, NE' energies for cols, then matrix
                # ne = NE (rows), lt = NE' (cols) per ENDF manual
                nek = ne; nel = lt
                ek = vcat(lst.data[1:nek], lst.data[nek+1:nek+nel])
                fvals = lst.data[nek+nel+1:np]
                push!(ni_subs, CovarianceBlock(mt, mt1, lb, lt, ek, fvals))
            end
        end
    end
    CovarianceData(Int(mat), Int(mt), Int(mt), nc_subs, ni_subs)
end

"""Compute relative covariance from an NI-type sub-subsection."""
ni_covariance(ni_data::CovarianceBlock, energy_grid::AbstractVector{<:Real}) =
    expand_covariance_block(ni_data, energy_grid)

"""Compute covariance from NC-type (derived reaction linear combination).
Returns zero matrix; full NC requires external cross section data."""
nc_covariance(nc_data::NamedTuple, energy_grid::AbstractVector{<:Real};
              xs_funcs::Union{Nothing,Dict}=nothing) =
    zeros(Float64, length(energy_grid)-1, length(energy_grid)-1)

"""Process MF33 covariance data into multigroup covariance matrices."""
function process_covariance(endf_file::AbstractString,
                             group_structure::AbstractVector{<:Real}; mts=:all)
    results = CovarianceMatrix[]; gb = collect(Float64, group_structure)
    open(endf_file, "r") do io
        seekstart(io); avail = Int[]
        while !eof(io)
            line = readline(io); p = rpad(line, 80)
            mf = Int(_parse_int(p[71:72])); mt = Int(_parse_int(p[73:75]))
            mf == 33 && mt > 0 && !(mt in avail) && push!(avail, mt)
        end
        target = mts === :all ? avail : intersect(collect(Int, mts), avail)
        matnum = _detect_mat_from_io(io)
        for mt in target
            try
                cd = read_mf33(io, matnum, mt)
                !isempty(cd.ni_subsections) && push!(results,
                    multigroup_covariance(cd.ni_subsections, gb))
            catch e
                @warn "process_covariance: skipping MT=$mt" exception=(e, catch_backtrace())
            end
        end
    end; results
end

"""Cov_g = J * Cov_params * J' (sandwich rule). Result is symmetric."""
function sandwich_covariance(J::AbstractMatrix{<:Real}, cov_params::AbstractMatrix{<:Real})
    C = J * cov_params * J'; 0.5*(C+C')
end

"""Finite-difference Jacobian d(sigma_g)/d(params)."""
function sensitivity_jacobian(sigma_func, params::AbstractVector{<:Real},
                               group_bounds::AbstractVector{<:Real}; dp::Real=1e-6)
    ng = length(group_bounds) - 1; np_p = length(params)
    J = zeros(Float64, ng, np_p)
    sigma_base = _group_avg_simple(sigma_func, params, group_bounds)
    for j in 1:np_p
        p_pert = copy(params); h = max(abs(params[j])*dp, dp); p_pert[j] += h
        J[:,j] = (_group_avg_simple(sigma_func, p_pert, group_bounds) .- sigma_base) ./ h
    end; J
end

"Check symmetry within tolerance."
is_symmetric(cm::CovarianceMatrix; atol::Real=1e-12) =
    maximum(abs.(cm.matrix - cm.matrix')) < atol

"Check positive semi-definiteness (all eigenvalues >= -atol)."
is_psd(cm::CovarianceMatrix; atol::Real=1e-10) =
    minimum(eigvals(Symmetric(cm.matrix))) >= -atol

_find_bin(x::Real, grid::AbstractVector{<:Real}) = begin
    for k in 1:length(grid)-1; x >= grid[k] && x < grid[k+1] && return k; end; 0
end

_upper_tri_index(i::Int, j::Int, n::Int) = n*(i-1) - div((i-1)*(i-2), 2) + (j-i+1)

function _build_union_grid(blocks::AbstractVector{CovarianceBlock}, gb::AbstractVector{<:Real})
    pts = Set{Float64}(Float64.(gb))
    for b in blocks, e in b.energies; push!(pts, Float64(e)); end
    sort!(collect(pts))
end

function _collapse_matrix(ugrid::AbstractVector{<:Real}, gb::AbstractVector{<:Real};
                          flux::Union{Nothing,AbstractVector{<:Real}}=nothing)
    ng = length(gb)-1; nu = length(ugrid)-1; T = zeros(Float64, ng, nu)
    if flux !== nothing && length(flux) == nu
        # Flux-weighted collapse (matching Fortran ERRORR)
        # T[g,k] = flux[k] * overlap / (sum of flux*overlap in group g)
        group_flux = zeros(Float64, ng)
        for k in 1:nu
            u_lo, u_hi = ugrid[k], ugrid[k+1]; u_w = u_hi - u_lo; u_w <= 0 && continue
            for g in 1:ng
                ov = min(u_hi, gb[g+1]) - max(u_lo, gb[g])
                if ov > 0
                    T[g,k] = flux[k] * ov / u_w
                    group_flux[g] += flux[k] * ov / u_w
                end
            end
        end
        for g in 1:ng
            group_flux[g] > 0 && (T[g,:] ./= group_flux[g])
        end
    else
        # Area-weighted collapse (flat flux assumption)
        for k in 1:nu
            u_lo, u_hi = ugrid[k], ugrid[k+1]; u_w = u_hi - u_lo; u_w <= 0 && continue
            for g in 1:ng
                ov = min(u_hi, gb[g+1]) - max(u_lo, gb[g]); ov > 0 && (T[g,k] = ov/u_w)
            end
        end
    end
    T
end

function _group_avg_simple(sigma_func, params, gb)
    ng = length(gb)-1; result = zeros(Float64, ng); npts = 32
    for g in 1:ng
        elo, ehi = gb[g], gb[g+1]; de = (ehi-elo)/npts
        s = 0.5*(sigma_func(elo, params) + sigma_func(ehi, params))
        for k in 1:npts-1; s += sigma_func(elo+k*de, params); end
        result[g] = s*de/(ehi-elo)
    end; result
end

function _detect_mat_from_io(io::IO)
    seekstart(io)
    while !eof(io)
        line = readline(io); p = rpad(line, 80)
        mat = Int(_parse_int(p[67:70])); mf = Int(_parse_int(p[71:72]))
        mat > 0 && mf > 0 && return mat
    end
    error("could not detect MAT number")
end
