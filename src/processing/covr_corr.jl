# covr_corr.jl -- Correlation and grid-truncation steps for covr.
#
# Ports two Fortran covr.f90 subroutines:
#   corr   (covr.f90:578-718)  -- convert covariance → std dev + correlation
#   truncg (covr.f90:939-1011) -- truncate energy grids on zero-xs ends

"""
    CorrResult

State carried out of `corr` and consumed by the plot writer. Mirrors the
Fortran-side variables `x, y, xx, xy, rsdx, rsdy, cf, ixmin, ixmax,
izero, ismall, izap, einc` (covr.f90:578-718).

After `corr` returns,
- `cf` is the *correlation* matrix (cov[i,j] / (rsdx[i]·rsdy[j])), or zero
  where either standard deviation vanished.
- `rsdx` is sqrt(diag) for the row reaction (mt) — relative std dev unless
  the source matrix was absolute and the user passed irelco=0.
- `rsdy` is sqrt(diag) for the column reaction (mt1).
- `xx`/`xy` retain the per-group cross sections (used by plotit when
  mat==mat1 ∧ mt==mt1 and mtflg=1; covr.f90:1048-1059).
- `izero=0` ⇒ all-null matrix (skip plot, register 'null' on nscr1).
- `ismall=0` ⇒ no |corr| ≥ tlev[1] (skip matrix, register 'small' on nscr1).
"""
struct CorrResult
    x::Vector{Float64}    # row group bounds (length ixmax+1) — used as plot x-axis
    y::Vector{Float64}    # column group bounds — equal to x for self/aligned tapes
    xx::Vector{Float64}   # MT cross section per group
    xy::Vector{Float64}   # MT1 cross section per group
    cf::Vector{Float64}   # ixmax*ixmax row-major correlation matrix
    rsdx::Vector{Float64}
    rsdy::Vector{Float64}
    ixmax::Int
    izero::Int
    ismall::Int
    izap::Int
    einc::Float64
    epmin::Float64        # possibly auto-bumped value (covr.f90:628)
end

"""
    corr(tape::ErrorrTape, mat, mt, mat1, mt1; irelco=1, epmin0=0.0,
         tlev_first::Real=0.001, plot_mode::Bool=true) -> CorrResult

Convert covariance values from `tape` for the (mat,mt,mat1,mt1) reaction
pair into standard deviations and correlation coefficients. Mirrors
covr.f90:578-718 (`corr`).

Arguments:
- `irelco`     : 0/1 absolute/relative on input (covr.f90:929 ratio path)
- `epmin0`     : baseline plot lower limit (post-`*rdn` from parse_covr)
- `tlev_first` : `xlev(1)` for the ismall trigger (covr.f90:683)
- `plot_mode`  : selects covr.f90's `nout.le.0` branch — matters for the
                 epmin auto-reset (covr.f90:620-634), which is plot-only.

Notes
- For self-covariance (`mat==mat1 && mt==mt1`) the Fortran skips the
  cross-pre-check (covr.f90:599-604) and computes rsdx==rsdy in one pass.
- For cross blocks the Fortran calls `covard` four times; here we mirror
  exactly that flow because each call's `cf` ends up holding a different
  matrix (final = cross), and `rsdx`/`rsdy` come from the two self diags.
"""
function corr(tape::ErrorrTape,
              mat::Integer, mt::Integer,
              mat1::Integer, mt1::Integer;
              irelco::Integer=1, epmin0::Real=0.0,
              tlev_first::Real=0.001,
              plot_mode::Bool=true)::CorrResult
    mat1 = mat1 == 0 ? mat : mat1
    mt1  = mt1  == 0 ? mt  : mt1
    self = (mat == mat1 && mt == mt1)
    ixmax = tape.ixmax

    # ---- pre-check for cross-reaction null matrix (covr.f90:599-604) ----
    if !self
        cross0 = covard(tape, mat, mt, mat1, mt1; irelco=irelco)
        if cross0.izero == 0
            return _null_corr_result(cross0, epmin0)
        end
    end

    # ---- rsdy from self(mat1, mt1) (covr.f90:606-650) ----
    self_y = covard(tape, mat1, mt1, mat1, mt1; irelco=irelco)
    rsdy = _diag_sqrt(self_y.cf, ixmax)

    if self
        rsdx = copy(rsdy)
        cf   = self_y.cf
        xx   = self_y.xx
        xy   = self_y.xy
        izero = self_y.izero
        izap  = self_y.izap
        einc  = self_y.einc
        x = copy(tape.groups)
        y = copy(tape.groups)
    else
        # ---- rsdx from self(mat, mt) (covr.f90:651-670) ----
        self_x = covard(tape, mat, mt, mat, mt; irelco=irelco)
        rsdx = _diag_sqrt(self_x.cf, ixmax)
        # ---- final cf = cross matrix (covr.f90:669) ----
        cross1 = covard(tape, mat, mt, mat1, mt1; irelco=irelco)
        cf   = cross1.cf
        xx   = cross1.xx
        xy   = cross1.xy
        izero = cross1.izero
        izap  = cross1.izap
        einc  = cross1.einc
        x = copy(tape.groups)
        y = copy(tape.groups)
    end

    # ---- epmin auto-reset for plot mode (covr.f90:620-634) ----
    epmin = Float64(epmin0)
    if plot_mode
        epmn = 0.499999e-4   # covr.f90:592
        zp4  = 0.4
        xsize = 4.25         # covr.f90:591
        emax  = x[ixmax + 1]
        xcycle = xsize / log10(emax / epmn)
        while xcycle < zp4
            epmn *= 2
            xcycle = xsize / log10(emax / epmn)
        end
        if epmn > epmin
            epmin = epmn
            @info @sprintf("corr: epmin reset to %12.4e, xcycle=%12.4e", epmin, xcycle)
        end
    end

    # ---- normalize cov → correlation, set ismall (covr.f90:672-688) ----
    ismall = 0
    if izero != 0
        @inbounds for i in 1:ixmax
            ind = ixmax * (i - 1)
            for j in 1:ixmax
                v = cf[ind + j]
                d = rsdx[i] * rsdy[j]
                if v != 0.0 && d != 0.0
                    c = v / d
                    cf[ind + j] = c
                    abs(c) >= tlev_first && (ismall = 1)
                else
                    cf[ind + j] = 0.0
                end
            end
        end
    end

    CorrResult(x, y, xx, xy, cf, rsdx, rsdy, ixmax,
               izero, ismall, izap, einc, epmin)
end

# Pull the standard-deviation vector from the diagonal of a flat row-major
# matrix. Negative diagonal entries clamp to zero (covr.f90:636-642).
function _diag_sqrt(cf::Vector{Float64}, ixmax::Int)::Vector{Float64}
    rsd = zeros(Float64, ixmax)
    @inbounds for i in 1:ixmax
        d = cf[ixmax * (i - 1) + i]
        rsd[i] = d > 0.0 ? sqrt(d) : 0.0
    end
    rsd
end

function _null_corr_result(c::CovardResult, epmin0::Real)::CorrResult
    ixmax = c.ixmax
    x = copy(c.groups)
    y = copy(c.groups)
    CorrResult(x, y, copy(c.xx), copy(c.xy), c.cf,
               zeros(Float64, ixmax), zeros(Float64, ixmax),
               ixmax, 0, 0, c.izap, c.einc, Float64(epmin0))
end

# ---------------------------------------------------------------------------
# truncg
# ---------------------------------------------------------------------------

"""
    truncg!(x, y, xx, xy, ixmax, epmin) -> ixmin::Int

Truncate the lower end of the MT and MT1 energy group structures to
eliminate zero-xs regions from plots, mutating `x[1]` and `y[1]` to the
new lower bound (Fortran in-place update). Returns the index of the new
lower group on the (1-based) row grid.

Direct port of covr.f90:939-1011.

Behaviour mirror:
- Defines `rlim_x` / `rlim_y` as `xslim*∫xs dE / (E_max - E_threshold)`.
- Walks groups from low to high; first index whose group lower bound
  exceeds `epmin` and whose xs is significant becomes `ixmin`.
- For thermal groups (group 1 below 1e-6 eV): clamps the lower bound to a
  decade if the next group is more than 10× wider.
- Final clamp `x[ixmin] >= epmin`, `y[ixmin] >= epmin`.

Constants identical to Fortran:
- `emin = 0.9999e6` (1 MeV roof for the threshold/`go to 120` short-cut)
- `xslim = 1e-4`    (relative xs cutoff)
"""
function truncg!(x::Vector{Float64}, y::Vector{Float64},
                 xx::Vector{Float64}, xy::Vector{Float64},
                 ixmax::Int, epmin::Float64)::Int
    emin  = 0.9999e6
    xslim = 1.0e-4

    # Build the rlim_x / rlim_y normalisation from cumulative ∫xs dE
    rlimx = 0.0
    rlimy = 0.0
    ethrx = x[1]
    ethry = y[1]
    for i in 2:ixmax
        rlimx += xx[i-1] * (x[i] - x[i-1])
        rlimy += xy[i-1] * (y[i] - y[i-1])
        xx[i-1] <= 0.0 && (ethrx = x[i])
        xy[i-1] <= 0.0 && (ethry = y[i])
    end
    if rlimx != 0.0 && rlimy != 0.0
        rlimx = xslim * rlimx / (y[ixmax + 1] - ethry)
        rlimy = xslim * rlimy / (x[ixmax + 1] - ethrx)
    end

    # Find ixmin from xx-side
    ixmin = 1
    for i in 1:ixmax
        if x[1 + i] > epmin
            (xx[i] >= xslim || xx[i] >= rlimx) && @goto found_x
            x[1 + i] > emin && @goto found_x
        end
        ixmin = i + 1
    end
    error("truncg!: bad data — could not locate ixmin (covr.f90:982)")
    @label found_x

    # Same logic for yy-side
    iymin = 1
    for i in 1:ixmax
        if y[1 + i] > epmin
            (xy[i] >= xslim || xy[i] >= rlimy) && @goto found_y
            y[1 + i] > emin && @goto found_y
        end
        iymin = i + 1
    end
    error("truncg!: bad data — could not locate iymin (covr.f90:992)")
    @label found_y

    iymin < ixmin && (ixmin = iymin)

    # Thermal-group decade clamp (covr.f90:996-1006)
    if ixmin == 1
        if 10 * x[1] < x[2] && 10 * y[1] < y[2]
            elo = log10(x[2] / 10)
            ielo = round(Int, elo)
            x[1] < 10.0^ielo && (x[1] = 10.0^ielo)
            elo = log10(y[2] / 10)
            ielo = round(Int, elo)
            y[1] < 10.0^ielo && (y[1] = 10.0^ielo)
        end
    end
    x[ixmin] < epmin && (x[ixmin] = epmin)
    y[ixmin] < epmin && (y[ixmin] = epmin)

    ixmin
end
