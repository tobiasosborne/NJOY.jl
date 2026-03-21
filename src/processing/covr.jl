# COVR -- Covariance post-processing (correlation matrices, relative std dev,
# and boxer-format output from ERRORR covariance data).
#
# Correspondence to NJOY2016 covr.f90:
#   corr        -> covariance_to_correlation (normalize C to [-1,1])
#   covard      -> (reads from CovarianceMatrix directly)
#   press       -> format_covariance_output (boxer BCD format)
#   truncg      -> (optional truncation in format_covariance_output)
#   plotit      -> (plot-ready data via CorrelationMatrix)
#   relative_std_dev -> sqrt(diag), matching rsdx/rsdy in covr.f90

# ============================================================================
# Types
# ============================================================================

"""
    CorrelationMatrix

Normalized correlation matrix with values in [-1, 1], plus associated
relative standard deviations and group structure. Produced from a
CovarianceMatrix via `covariance_to_correlation`.

Fields:
- `mt1`, `mt2`: reaction identifiers
- `groups`: energy group boundaries (length ng+1)
- `correlation`: (ng x ng) matrix with entries in [-1, 1]
- `std_dev_1`, `std_dev_2`: relative standard deviation vectors for mt1, mt2
"""
struct CorrelationMatrix
    mt1::Int
    mt2::Int
    groups::Vector{Float64}
    correlation::Matrix{Float64}
    std_dev_1::Vector{Float64}
    std_dev_2::Vector{Float64}
end

function Base.show(io::IO, cm::CorrelationMatrix)
    ng = size(cm.correlation, 1)
    print(io, "CorrelationMatrix(MT", cm.mt1, "/MT", cm.mt2, ", ", ng, " groups)")
end

# ============================================================================
# Core functions
# ============================================================================

"""
    relative_std_dev(cov::CovarianceMatrix) -> Vector{Float64}

Extract the relative standard deviation from the diagonal of a covariance
matrix: rsd[i] = sqrt(C[i,i]). Negative diagonal entries are clamped to zero.
This matches the Fortran: `rsdy(i) = sqrt(cf(ixmax*(i-1)+i))`.
"""
function relative_std_dev(cov::CovarianceMatrix)
    ng = size(cov.matrix, 1)
    rsd = Vector{Float64}(undef, ng)
    for i in 1:ng
        d = cov.matrix[i, i]
        rsd[i] = d > 0.0 ? sqrt(d) : 0.0
    end
    rsd
end

"""
    covariance_to_correlation(cov::CovarianceMatrix) -> CorrelationMatrix
    covariance_to_correlation(cov::CovarianceMatrix, cov2::CovarianceMatrix) -> CorrelationMatrix

Convert a relative covariance matrix to a correlation matrix. The diagonal
provides relative standard deviations and the off-diagonal elements are
normalized to the range [-1, 1].

For a cross-reaction covariance (mt1 != mt2), pass the two self-covariance
matrices to get the correct normalization for each reaction.

Matches the Fortran `corr` subroutine in covr.f90 (lines 578-718):
    cf(i,j) = cf(i,j) / (rsdx(i) * rsdy(j))
"""
function covariance_to_correlation(cov::CovarianceMatrix)
    rsd = relative_std_dev(cov)
    _build_correlation(cov, rsd, rsd)
end

function covariance_to_correlation(cov_cross::CovarianceMatrix,
                                    cov_mt1::CovarianceMatrix,
                                    cov_mt2::CovarianceMatrix)
    rsd1 = relative_std_dev(cov_mt1)
    rsd2 = relative_std_dev(cov_mt2)
    _build_correlation(cov_cross, rsd1, rsd2)
end

function _build_correlation(cov::CovarianceMatrix,
                            rsd_row::Vector{Float64},
                            rsd_col::Vector{Float64})
    ng = size(cov.matrix, 1)
    @assert length(rsd_row) == ng && length(rsd_col) == ng
    corr = zeros(Float64, ng, ng)
    for i in 1:ng, j in 1:ng
        denom = rsd_row[i] * rsd_col[j]
        if denom > 0.0 && cov.matrix[i, j] != 0.0
            corr[i, j] = cov.matrix[i, j] / denom
        end
    end
    # Clamp to [-1, 1] for numerical safety
    clamp!(corr, -1.0, 1.0)
    CorrelationMatrix(cov.mt1, cov.mt2, copy(cov.energy_groups),
                      corr, copy(rsd_row), copy(rsd_col))
end

"""
    max_abs_correlation(cm::CorrelationMatrix) -> Float64

Return the maximum absolute off-diagonal correlation coefficient.
Useful for identifying strongly correlated groups.
"""
function max_abs_correlation(cm::CorrelationMatrix)
    ng = size(cm.correlation, 1)
    mx = 0.0
    for i in 1:ng, j in 1:ng
        i != j && (mx = max(mx, abs(cm.correlation[i, j])))
    end
    mx
end

"""
    is_valid_correlation(cm::CorrelationMatrix; atol::Real=1e-10) -> Bool

Verify that a correlation matrix has diagonal entries of 1.0 (for self-
correlations) and all entries in [-1, 1].
"""
function is_valid_correlation(cm::CorrelationMatrix; atol::Real=1e-10)
    ng = size(cm.correlation, 1)
    all(x -> -1.0 - atol <= x <= 1.0 + atol, cm.correlation) || return false
    # For self-correlations (mt1==mt2), diagonal should be 1.0
    if cm.mt1 == cm.mt2
        for i in 1:ng
            abs(cm.correlation[i, i] - 1.0) > atol && cm.std_dev_1[i] > 0 && return false
        end
    end
    true
end

# ============================================================================
# Boxer format output (BCD covariance library)
# ============================================================================

"""
    format_covariance_output(covs::Vector{CovarianceMatrix};
                             format=:boxer, matype=3, libid="njoy",
                             description="covariance library") -> String

Format covariance data for output. Supported formats:
- `:boxer` -- compressed BCD format matching NJOY COVR library output
- `:csv`   -- comma-separated values (one matrix per block)

`matype` controls whether covariances (3) or correlations (4) are written
in boxer format, matching the Fortran `press` subroutine.
"""
function format_covariance_output(covs::Vector{CovarianceMatrix};
                                   format::Symbol=:boxer,
                                   matype::Int=3,
                                   libid::String="njoy",
                                   description::String="covariance library")
    format === :boxer && return _format_boxer(covs; matype, libid, description)
    format === :csv   && return _format_csv(covs; matype)
    error("format_covariance_output: unsupported format=$format (use :boxer or :csv)")
end

function _format_boxer(covs::Vector{CovarianceMatrix};
                        matype::Int=3, libid::String="njoy",
                        description::String="covariance library")
    isempty(covs) && return ""
    ng = size(covs[1].matrix, 1)
    hlibid = rpad(libid, 12)[1:12]
    if matype == 3
        hlibid = rpad(libid[1:min(6,end)] * "-a-" * lpad(string(ng), 3), 12)[1:12]
    elseif matype == 4
        hlibid = rpad(libid[1:min(6,end)] * "-b-" * lpad(string(ng), 3), 12)[1:12]
    end
    lines = String[]
    push!(lines, " " * hlibid * "  " * rpad(description, 21)[1:21])
    # Group boundaries
    gb = covs[1].energy_groups
    push!(lines, @sprintf(" %4d group bounds (eV)", ng))
    for i in 1:length(gb)
        push!(lines, @sprintf(" %13.6e", gb[i]))
    end
    # Each covariance matrix
    for cov in covs
        M = matype == 4 ? _to_correlation_matrix(cov) : cov.matrix
        rsd = relative_std_dev(cov)
        push!(lines, @sprintf(" mat/mt %4d/%3d  mat1/mt1 %4d/%3d  type=%d  ng=%d",
                              0, cov.mt1, 0, cov.mt2, matype, ng))
        # Cross sections / std devs
        push!(lines, " rsd")
        for i in 1:ng
            push!(lines, @sprintf(" %13.6e", rsd[i]))
        end
        # Matrix (upper triangle for symmetric, full otherwise)
        is_sym = cov.mt1 == cov.mt2
        push!(lines, is_sym ? " matrix (upper triangle)" : " matrix (full)")
        for i in 1:ng
            jstart = is_sym ? i : 1
            for j in jstart:ng
                push!(lines, @sprintf(" %4d %4d %13.6e", i, j, M[i, j]))
            end
        end
    end
    join(lines, "\n") * "\n"
end

function _to_correlation_matrix(cov::CovarianceMatrix)
    rsd = relative_std_dev(cov)
    ng = size(cov.matrix, 1)
    corr = zeros(Float64, ng, ng)
    for i in 1:ng, j in 1:ng
        d = rsd[i] * rsd[j]
        corr[i, j] = d > 0.0 ? cov.matrix[i, j] / d : 0.0
    end
    clamp!(corr, -1.0, 1.0)
    corr
end

function _format_csv(covs::Vector{CovarianceMatrix}; matype::Int=3)
    lines = String[]
    for cov in covs
        M = matype == 4 ? _to_correlation_matrix(cov) : cov.matrix
        ng = size(M, 1)
        push!(lines, "# MT$(cov.mt1)/MT$(cov.mt2) $(ng)x$(ng) " *
                     (matype == 4 ? "correlation" : "covariance"))
        for i in 1:ng
            push!(lines, join([@sprintf("%.6e", M[i,j]) for j in 1:ng], ","))
        end
        push!(lines, "")
    end
    join(lines, "\n")
end

# ============================================================================
# Batch processing
# ============================================================================

"""
    process_covr(endf_file::AbstractString,
                 group_bounds::AbstractVector{<:Real};
                 mts=:all, output_format=:boxer, matype=3) ->
        (matrices=Vector{CovarianceMatrix},
         correlations=Vector{CorrelationMatrix},
         formatted=String)

End-to-end COVR pipeline: reads ENDF MF33 covariance data via ERRORR processing,
then computes correlation matrices and formatted output.
"""
function process_covr(endf_file::AbstractString,
                      group_bounds::AbstractVector{<:Real};
                      mts=:all, output_format::Symbol=:boxer, matype::Int=3)
    covs = process_covariance(endf_file, group_bounds; mts=mts)
    corrs = CorrelationMatrix[]
    for cov in covs
        push!(corrs, covariance_to_correlation(cov))
    end
    formatted = format_covariance_output(covs; format=output_format, matype=matype)
    (matrices=covs, correlations=corrs, formatted=formatted)
end

"""
    covr_summary(corrs::Vector{CorrelationMatrix}) -> String

Print a summary of correlation matrices: reaction pairs, group counts,
max absolute off-diagonal correlation, and mean relative std dev.
"""
function covr_summary(corrs::Vector{CorrelationMatrix})
    lines = String[]
    push!(lines, "COVR Summary: $(length(corrs)) reaction pair(s)")
    push!(lines, "-"^60)
    for cm in corrs
        ng = size(cm.correlation, 1)
        max_corr = max_abs_correlation(cm)
        mean_rsd = sum(cm.std_dev_1) / max(ng, 1)
        push!(lines, @sprintf("  MT%d/MT%d: %d groups, max|rho|=%.4f, mean_rsd=%.4e",
                              cm.mt1, cm.mt2, ng, max_corr, mean_rsd))
    end
    join(lines, "\n") * "\n"
end
