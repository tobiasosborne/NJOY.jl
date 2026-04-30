# covr_module -- Top-level orchestration of the covr Fortran subroutine.
#
# Mirrors covr.f90:48-505 (`covr` main). Reads the user's input cards
# (already parsed into CovrParams), opens the errorr-output covariance
# tape, loops over (mat,mt,mat1,mt1) cases, and emits either:
#   * a viewr plot tape at unit `nplot` (plot mode, nout ≤ 0); or
#   * a boxer-format covariance/correlation library at unit `nout`
#     (library mode, nout > 0).
#
# Both modes coexist on the same input deck. The Fortran handles them in
# two disjoint code paths (covr.f90:287-306 vs. 425-459); this Julia
# orchestrator follows the same dispatch.
#
# Inputs travelled per case:
#   1. read_errorr_tape(path)         -- one ErrorrTape per `mat = imat(n)`
#   2. expndo (covr.f90:508-576)      -- expand `imt(n) ≤ 0` into the full
#                                        list of MT pairs in ENDF order
#   3. corr   (covr.f90:578-718)      -- (cov → corr, std-dev) per pair
#   4. truncg (covr.f90:939-1011)     -- drop zero-xs leading groups
#   5. plotit + matshd (plot mode)    -- 3 frames per pair to `nplot`
#      OR
#      press (library mode)           -- compressed boxer to `nout`
#
# All Fortran constants and FP semantics are matched exactly; deviations
# are flagged in inline comments where they appear.

using Printf

"""
    covr_module(tapes::TapeManager, params::CovrParams)

Run one covr invocation.  Reads the errorr-output tape mapped to
`params.nin`, and writes either a plot-tape (`params.nplot`) or a
boxer-format library (`params.nout`), per the input-card configuration.

Errors propagate as Julia `ErrorException`s with Fortran citation in the
message. No silent fallbacks (per CLAUDE.md Rule 6 / fail-fast-fail-loud):
malformed sections raise rather than emit empty output.
"""
function covr_module(tapes::TapeManager, params::CovrParams)
    plot_mode = params.nout <= 0
    in_path = resolve(tapes, params.nin)
    isfile(in_path) ||
        error("covr: input tape unit $(params.nin) → $in_path does not exist " *
              "(covr.f90:184-185 `openz(nin,0)`)")

    plot_path = plot_mode ? resolve(tapes, params.nplot) : ""
    out_path  = plot_mode ? "" : resolve(tapes, params.nout)

    # Always ensure the registered output paths exist on the work_dir,
    # even if no frames are written (matches Fortran's habit of touching
    # the file via openz).
    if plot_mode
        params.nplot != 0 && touch(plot_path)
    else
        touch(out_path)
    end

    cases = params.cases
    isempty(cases) && return _covr_finish(tapes, params, plot_path, out_path)

    if plot_mode
        if params.nplot == 0
            @info "covr: plot mode but nplot=0 — no plot tape produced"
            return nothing
        end
        open(plot_path, "w") do plotio
            write_covr_header!(plotio, params.icolor)
            _covr_run_cases!(plotio, params, cases, in_path; plot_mode=true)
            write_covr_footer!(plotio)
        end
        register!(tapes, params.nplot, plot_path)
    else
        open(out_path, "w") do outio
            _covr_run_cases!(outio, params, cases, in_path; plot_mode=false)
        end
        register!(tapes, params.nout, out_path)
    end
    nothing
end

# Final no-op shutdown when there were zero cases. We still register the
# empty tapes so downstream modules don't trip on missing files.
function _covr_finish(tapes::TapeManager, params::CovrParams,
                      plot_path::String, out_path::String)
    if params.nout <= 0
        params.nplot != 0 && register!(tapes, params.nplot, plot_path)
    else
        register!(tapes, params.nout, out_path)
    end
    nothing
end

# Run the full per-case loop. Reads the errorr tape once per unique MAT
# encountered (mirrors Fortran's `repoz(nin); finds(mat,...)` pattern,
# which re-seeks rather than re-opens, but avoids re-parsing).
function _covr_run_cases!(io::IO, params::CovrParams,
                          cases::Vector{CovrCase}, in_path::String;
                          plot_mode::Bool)
    state = nothing
    boxer = nothing
    tapes_by_mat = Dict{Int, ErrorrTape}()
    iomit_log = String[]

    for (n, case) in pairs(cases)
        mat = case.mat
        mat == 0 && continue

        tape = get!(tapes_by_mat, mat) do
            read_errorr_tape(in_path; mat=mat)
        end
        mf35 = tape.mfflg == -12 ? 5 : 3

        # Expand the MT list when imt(n) ≤ 0 (covr.f90:369).
        mt_pairs = if case.mt <= 0
            mts = expand_mt_list(tape;
                                 strip_mt   = case.mt,
                                 strip_mat1 = case.mat1,
                                 strip_mt1  = case.mt1)
            # all (i,j) with i ≤ j (covr.f90:558-568)
            ps = Tuple{Int,Int,Int,Int}[]
            for im in 1:length(mts), jm in im:length(mts)
                push!(ps, (mat, mts[im], mat, mts[jm]))
            end
            ps
        else
            mat1 = case.mat1 == 0 ? mat : case.mat1
            mt1  = case.mt1  == 0 ? case.mt : case.mt1
            [(mat, case.mt, mat1, mt1)]
        end

        if plot_mode && state === nothing
            iza = tape.iza
            state = PlotState(params; iverf=tape.iverf)
        end
        if !plot_mode && boxer === nothing
            boxer = BoxerWriter(io;
                                matype=params.matype,
                                hlibid=isempty(params.hlibid) ? "njoy  " : params.hlibid,
                                hdescr=params.hdescr)
            # Per covr.f90:443: write the group structure as itype=0 first.
            x = collect(Float64, tape.groups)
            boxer_press!(boxer, mat, 0, 0, 0, 0,
                         reshape(x, length(x), 1);
                         nrow=length(x), ncol=1, nvf=10, ncf=3)
        end

        for (ne, (mat_, mt, mat1, mt1)) in pairs(mt_pairs)
            cr = corr(tape, mat_, mt, mat1, mt1;
                      irelco=params.irelco, epmin0=params.epmin,
                      tlev_first=state === nothing ? 0.001 : state.xlev[1],
                      plot_mode=plot_mode)

            if plot_mode
                _covr_plot_one!(io, state::PlotState, cr, tape, mat_, mt, mat1, mt1, ne, iomit_log)
            else
                _covr_press_one!(boxer::BoxerWriter, params, cr, mat_, mt, mat1, mt1, ne, n)
            end
        end
    end

    plot_mode && !isempty(iomit_log) &&
        @info "covr: $(length(iomit_log)) plots omitted (null/small)"
end

# One plot-mode subcase. Returns nothing; updates `state.nfig` and may
# append to `iomit_log` (used by Fortran's nscr1 statistics tape).
function _covr_plot_one!(io::IO, state::PlotState, cr::CorrResult,
                          tape::ErrorrTape,
                          mat::Integer, mt::Integer,
                          mat1::Integer, mt1::Integer, ne::Integer,
                          iomit_log::Vector{String})
    # Fortran flags (covr.f90:386-391): if cf is null OR ismall=0, log + skip.
    if cr.izero == 0 || cr.ismall == 0
        reason = cr.izero == 0 ? "null" : "small"
        push!(iomit_log,
              @sprintf("%4d  %3d  %4d  %3d  %s", mat, mt, mat1, mt1, reason))
        return
    end

    # Truncate the leading zero-xs groups (covr.f90:404).
    ixmin = truncg!(cr.x, cr.y, cr.xx, cr.xy, cr.ixmax, cr.epmin)

    # Self-cov sanity short-circuit (covr.f90:411): if the truncate left
    # the matrix empty, mark izero=2 and skip.
    cr.izero == 2 && return

    mfflg = tape.mfflg
    mf35  = mfflg == -12 ? 5 : 3
    iza   = tape.iza
    iza1  = tape.iza   # single-MAT tape ⇒ iza1 == iza
    ishade = write_covr_frames!(io, state, cr,
                                 mat, mt, mat1, mt1,
                                 mfflg, mf35, cr.izap, ne,
                                 iza, iza1, ixmin)
    if ishade == 0
        push!(iomit_log,
              @sprintf("%4d  %3d  %4d  %3d  empty", mat, mt, mat1, mt1))
        return
    end
    state.nfig += 1
    nothing
end

# One library-mode subcase. Mirrors covr.f90:425-459.
function _covr_press_one!(boxer::BoxerWriter, params::CovrParams,
                           cr::CorrResult,
                           mat::Integer, mt::Integer,
                           mat1::Integer, mt1::Integer,
                           ne::Integer, ncase_idx::Integer)
    cr.izero == 0 && return
    ixmax = cr.ixmax

    if ne == 1
        suffix = params.matype == 4 ? "-b-" : "-a-"
        boxer.hlibid = rpad(rstrip(params.hlibid) * suffix * lpad(string(ixmax), 3), 12)[1:12]
    end

    # Cross sections + std dev for self-cov pairs (covr.f90:444-450).
    if mat == mat1 && mt == mt1
        boxer_press!(boxer, mat, mt, mat1, mt1, 1,
                     reshape(cr.xx, length(cr.xx), 1);
                     nrow=ixmax, ncol=1, nvf=10, ncf=3)
        boxer_press!(boxer, mat, mt, mat1, mt1, 2,
                     reshape(cr.rsdx, length(cr.rsdx), 1);
                     nrow=ixmax, ncol=1, nvf=10, ncf=3)
    end

    # The matrix itself.  ncol=0 marks symmetric upper-triangle storage.
    M = reshape(cr.cf, ixmax, ixmax)
    nvf = params.matype == 3 ? 10 : 7
    ncf = ixmax <= 30 ? 4 : (ixmax <= 100 ? 5 : 6)
    ncol = (mat == mat1 && mt == mt1) ? 0 : ixmax
    boxer_press!(boxer, mat, mt, mat1, mt1, params.matype, M;
                 nrow=ixmax, ncol=ncol, nvf=nvf, ncf=ncf)
    nothing
end
