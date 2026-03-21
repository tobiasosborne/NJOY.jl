# Doppler broadening pipeline (BROADR) -- PROPOSAL B
# Kernel math lives in sigma1.jl; this file does orchestration.
# Reuses the generic adaptive_reconstruct from adaptive_grid.jl.
# Correspondence: broadr->doppler_broaden, broadn->adaptive_reconstruct,
#                 bsigma->sigma1_at, thinb->thin_xs

"""
    doppler_broaden(energies, xs, T, awr; tol=0.001, ...) -> (new_e, new_xs)
Doppler-broaden a piecewise-linear cross section from 0K to temperature T.
Uses SIGMA1 exact kernel + adaptive grid refinement.
"""
function doppler_broaden(energies::AbstractVector{<:Real},
                         xs::AbstractVector{<:Real},
                         T::Real, awr::Real;
                         tol::Real = 0.001, do_thin::Bool = true,
                         errmax::Real = 10*tol, errint::Real = tol/20000,
                         max_depth::Int = 30)
    n = length(energies)
    @assert n == length(xs) && n >= 2 && T > 0 && awr > 0
    @assert issorted(energies) "energies must be sorted"
    alpha = awr / (PhysicsConstants.bk * T)
    seg_e = Float64.(energies); seg_xs = Float64.(xs)
    broadened_eval = E -> (sigma1_at(E, seg_e, seg_xs, alpha),)
    initial_grid = _prepare_grid(seg_e, seg_xs)
    config = AdaptiveConfig(Float64(tol); errmax=Float64(errmax),
                            errint=Float64(errint), max_depth=max_depth)
    out_e, out_v = adaptive_reconstruct(broadened_eval, initial_grid, config)
    new_xs = vec(out_v)
    do_thin && ((out_e, new_xs) = thin_xs(out_e, new_xs; tol=Float64(tol)))
    return (out_e, new_xs)
end

"""
    doppler_broaden_multi(energies, xs_matrix, T, awr; ...) -> (new_e, new_xs)
Broaden multiple cross sections simultaneously on a shared adaptive grid.
`xs_matrix` is (N, n_reactions). Returns matrix output.
"""
function doppler_broaden_multi(energies::AbstractVector{<:Real},
                               xs_matrix::AbstractMatrix{<:Real},
                               T::Real, awr::Real;
                               tol::Real = 0.001, do_thin::Bool = true,
                               errmax::Real = 10*tol, errint::Real = tol/20000,
                               max_depth::Int = 30)
    n_pts, n_reac = size(xs_matrix)
    @assert n_pts == length(energies) && T > 0
    alpha = awr / (PhysicsConstants.bk * T)
    seg_e = Float64.(energies)
    broadened_eval = E -> ntuple(i -> sigma1_at(E, seg_e, @view(xs_matrix[:,i]), alpha), n_reac)
    initial_grid = filter(e -> e > 0.0, seg_e); sort!(initial_grid); unique!(initial_grid)
    config = AdaptiveConfig(Float64(tol); errmax=Float64(errmax),
                            errint=Float64(errint), max_depth=max_depth)
    out_e, out_v = adaptive_reconstruct(broadened_eval, initial_grid, config)
    do_thin && ((out_e, out_v) = thin_xs(out_e, out_v; tol=Float64(tol)))
    return (out_e, out_v)
end

"""
    doppler_broaden(pendf::PointwiseMaterial, T_new; T_old=0.0, awr=1.0,
                    tol=0.001, thnmax=Inf) -> PointwiseMaterial
Broaden a PointwiseMaterial from T_old to T_new.
"""
function doppler_broaden(pendf::PointwiseMaterial, T_new::Float64;
                         T_old::Float64=0.0, awr::Float64=1.0,
                         tol::Float64=0.001, thnmax::Float64=Inf)
    T_eff = T_new - T_old
    T_eff == 0.0 && return pendf
    T_eff < 0.0 && error("doppler_broaden: T_new must be >= T_old")
    emax = isfinite(thnmax) && thnmax > 0.0 ? thnmax : pendf.energies[end]
    idx_broad = findall(e -> e <= emax, pendf.energies)
    idx_copy  = findall(e -> e >  emax, pendf.energies)
    (isempty(idx_broad) || length(idx_broad) < 2) && return pendf
    out_e, out_v = doppler_broaden_multi(
        pendf.energies[idx_broad], pendf.cross_sections[idx_broad, :],
        T_eff, awr; tol=tol, do_thin=true)
    if !isempty(idx_copy)
        out_e = vcat(out_e, pendf.energies[idx_copy])
        out_v = vcat(out_v, pendf.cross_sections[idx_copy, :])
    end
    return PointwiseMaterial(pendf.mat, out_e, out_v, copy(pendf.mt_list))
end

# ==========================================================================
# Thinning (matches NJOY2016 thinb)
# ==========================================================================

"""
    thin_xs(energies, xs; tol=0.001, step_max=1.24) -> (thinned_e, thinned_xs)
Remove redundant points where linear interpolation suffices within `tol`.
"""
function thin_xs(energies::AbstractVector{<:Real}, xs::AbstractVector{<:Real};
                 tol::Real=0.001, step_max::Real=1.24)
    n = length(energies); @assert n == length(xs)
    n <= 2 && return (collect(Float64, energies), collect(Float64, xs))
    keep = _compute_thin_mask(energies, reshape(Float64.(xs), n, 1),
                              Float64(tol), Float64(step_max))
    idx = findall(keep)
    return (Float64.(energies[idx]), Float64.(xs[idx]))
end

function thin_xs(energies::AbstractVector{<:Real}, xs::AbstractMatrix{<:Real};
                 tol::Real=0.001, step_max::Real=1.24)
    n = size(xs, 1); @assert n == length(energies)
    n <= 2 && return (collect(Float64, energies), Float64.(xs))
    keep = _compute_thin_mask(energies, Float64.(xs), Float64(tol), Float64(step_max))
    idx = findall(keep)
    return (Float64.(energies[idx]), Float64.(xs[idx, :]))
end

# Core thinning: test ALL reactions simultaneously for union grid preservation.
function _compute_thin_mask(energies, xs_matrix::Matrix{Float64},
                            tol::Float64, step_max::Float64)
    n, nr = size(xs_matrix)
    keep = falses(n); keep[1] = true; keep[n] = true
    last_kept = 1
    for i in 3:n
        can_thin = true
        e_lo = energies[last_kept]; denom = energies[i] - e_lo
        denom <= 0.0 && continue
        for j in (last_kept+1):(i-1)
            energies[j] >= step_max*e_lo && (can_thin = false; break)
            for r in 1:nr
                slope = (xs_matrix[i,r] - xs_matrix[last_kept,r]) / denom
                sp = xs_matrix[last_kept,r] + slope*(energies[j] - e_lo)
                abs(sp - xs_matrix[j,r]) > tol*abs(xs_matrix[j,r]) && (can_thin = false; break)
            end
            can_thin || break
        end
        if !can_thin; keep[i-1] = true; last_kept = i-1; end
    end
    return keep
end

# ==========================================================================
# Grid enrichment: add midpoints near sharp slope changes
# ==========================================================================

function _prepare_grid(seg_e, seg_xs)
    grid = filter(e -> e > 0.0, seg_e); sort!(grid); unique!(grid)
    _enrich_broadr_grid(grid, seg_e, seg_xs)
end

function _enrich_broadr_grid(grid::Vector{Float64}, seg_e, seg_xs)
    enriched = copy(grid); n = length(seg_e)
    n <= 2 && return enriched
    for i in 2:(n-1)
        de_lo = seg_e[i] - seg_e[i-1]; de_hi = seg_e[i+1] - seg_e[i]
        (de_lo <= 0.0 || de_hi <= 0.0) && continue
        sl = (seg_xs[i] - seg_xs[i-1])/de_lo; sh = (seg_xs[i+1] - seg_xs[i])/de_hi
        ms = max(abs(sl), abs(sh))
        if ms > 0 && abs(sh - sl) > 0.5*ms
            push!(enriched, 0.5*(seg_e[i-1] + seg_e[i]))
            push!(enriched, 0.5*(seg_e[i] + seg_e[i+1]))
        end
    end
    sort!(enriched); unique!(enriched)
    return enriched
end
