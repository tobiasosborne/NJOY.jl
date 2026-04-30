# covr_io.jl -- Read covariance data from an errorr-output ENDF tape.
#
# Ports covr.f90:720-937 (`covard`) — the data ingest for covr. The Fortran
# subroutine streams an errorr-format covariance tape, extracting:
#   * MF1/MT451 group boundaries  → y(1:ixp)
#   * MF3 per-group cross sections for the active reaction → xx(1:ixmax)
#                                            and the cross reaction → xy(1:ixmax)
#   * MF33/MF34/MF35/MF40 covariance subsections for the (mat,mt,mat1,mt1)
#     pair, which errorr writes one matrix-row per LIST record (L2=lgp1
#     starting column, N2=kl row, NPL covariance values along the row),
#     terminated when kl ≥ kgp (kgp = NI from the sub-section header).
#
# Post-processing inside covard (covr.f90:892-934) eliminates spurious
# off-diagonal covariances in zero-xs regions, optionally converts
# absolute → relative (irelco != 1), and sets izero=0 if every covariance
# value ended up zero (null-matrix flag, used to skip empty plots).
#
# This file produces a CovardResult for one (mat,mt,mat1,mt1) request. The
# orchestration layer calls into here per (ne) iteration of the outer loop
# in covr.f90 (lines 372-470).

using Printf

"""
    ErrorrTape

Lightweight in-memory view of an errorr-output ENDF tape. Holds the parsed
group structure, MF1/451 metadata, MF3 cross-section vectors, and raw MF33/
MF34/MF35/MF40 line streams for later subsection extraction.

Built once per case (covr.f90 outer `do 130 n=1,ncase` loop) by
`read_errorr_tape(path)`. Callers query the per-MT covariance matrix via
`extract_covariance_matrix(tape, mat, mf3x, mt, mat1, mt1)`.

Fields:
- `tpid`     : 66-char tape ID line
- `mat`      : MAT number (single-material errorr tapes are the norm)
- `mfflg`    : `scr(5)` from MF1/MT451 head — marks tape kind
                  -11/-14 → MF35 maps to MF3 (normal covariance)
                  -12     → MF35 maps to MF5 (spectrum / nubar)
                  any other → covr aborts (covr.f90:325-327)
- `iverf`    : ENDF format version from MF1/MT451 (`l1h` of second contio)
- `iza`      : ZA from the MF1/MT451 HEAD record (1000*Z + A)
- `awr`      : AWR from MF1/MT451 HEAD
- `ixmax`    : number of groups (= L1 of MF1/MT451 LIST head)
- `ixp`      : number of group boundaries (= ixmax + 1, = N1 of LIST head)
- `groups`   : group boundaries (ixp values, monotonically increasing)
- `mf3_xs`   : MT → length-`ixmax` Vector{Float64} (per-group cross section)
- `mf3_einc` : MT → C2 of the LIST head (the "incident energy" stored in
               MF5-style records; zero for normal MF3)
- `mf33_raw` : (mt, mat1, mt1) → CovarianceSubsection block, lazily filled
               on first request via `extract_covariance_subsection`
- `sections` : Map (mf, mt) → Vector{String} of raw lines for late access
"""
struct ErrorrTape
    tpid::String
    mat::Int
    mfflg::Int
    iverf::Int
    iza::Int
    awr::Float64
    ixmax::Int
    ixp::Int
    groups::Vector{Float64}
    mf3_xs::Dict{Int, Vector{Float64}}
    mf3_einc::Dict{Int, Float64}
    sections::Dict{Tuple{Int,Int}, Vector{String}}
end

"""
    read_errorr_tape(path::AbstractString; mat::Integer=0) -> ErrorrTape

Read an errorr-output tape and parse all MF1/MF3 records into the
ErrorrTape struct. MF33/MF34/MF35/MF40 raw lines are kept for later
subsection extraction.

If `mat == 0` (the default) the tape's first material is used; this
matches covr.f90 line 339 (`mat = imat(n)`) — covr requires the user to
specify the MAT explicitly via card 4, so multi-material tapes are
handled by the caller.
"""
function read_errorr_tape(path::AbstractString; mat::Integer=0)::ErrorrTape
    tape = read_pendf(path)
    isempty(tape.materials) && error("read_errorr_tape: empty tape at $path")
    target_mat = mat == 0 ? tape.materials[1].mat : Int(mat)
    mat_idx = findfirst(m -> m.mat == target_mat, tape.materials)
    mat_idx === nothing &&
        error("read_errorr_tape: MAT=$target_mat not on tape $path")
    material = tape.materials[mat_idx]

    isempty(material.mf1_lines) &&
        error("read_errorr_tape: MAT=$target_mat has no MF1/MT451 (covr.f90:319)")

    iza, awr, mfflg, iverf, ixmax, ixp, groups =
        _parse_errorr_mf1_451(material.mf1_lines)
    ixmax = Int(ixmax); ixp = Int(ixp)
    iverf = Int(iverf); mfflg = Int(mfflg); iza = Int(iza)

    # Sanity check (covr.f90:319-327)
    mfflg == -11 || mfflg == -12 || mfflg == -14 ||
        error("read_errorr_tape: illegal mfflg=$mfflg (expected -11/-12/-14, covr.f90:325-327)")

    mf3_xs   = Dict{Int, Vector{Float64}}()
    mf3_einc = Dict{Int, Float64}()
    sections = Dict{Tuple{Int,Int}, Vector{String}}()

    for sec in material.sections
        if sec.mf == 3 || sec.mf == 5
            xs, einc = _parse_errorr_mf3_section(sec.lines, ixmax)
            mf3_xs[sec.mt] = xs
            mf3_einc[sec.mt] = einc
        end
        sections[(sec.mf, sec.mt)] = sec.lines
    end

    ErrorrTape(tape.tpid, target_mat, mfflg, iverf, iza, awr,
               ixmax, ixp, groups, mf3_xs, mf3_einc, sections)
end

# Parse the MF1/MT451 of an errorr tape. Two records:
#   1. HEAD CONT  : C1=ZA, C2=AWR, L1=iverf (≥5 since errorr stamps ENDF version),
#                   L2=lrp, N1=mfflg, N2=0
#       (covr.f90:317-320 `mfflg=nint(scr(5))` reads the N1 of the FIRST contio)
#   2. LIST       : C1=0, C2=0, L1=ixmax (#groups), L2=0, N1=ixp (#boundaries),
#                   N2=0; payload = ixp group boundaries.
function _parse_errorr_mf1_451(lines::Vector{String})
    length(lines) < 2 &&
        error("MF1/MT451 too short (need head + list head minimum)")
    # ---- HEAD ----
    f1, _, _, _, _ = parse_endf_line(lines[1])
    iza   = round(Int, parse_endf_float(f1[1]))
    awr   = parse_endf_float(f1[2])
    iverf = _parse_int(f1[3])
    mfflg = _parse_int(f1[5])
    # ---- LIST head ----
    f2, _, _, _, _ = parse_endf_line(lines[2])
    ixmax = _parse_int(f2[3])
    ixp   = _parse_int(f2[5])
    # Payload starts at lines[3]
    groups = Vector{Float64}(undef, ixp)
    k = 0
    li = 3
    while k < ixp
        li > length(lines) &&
            error("MF1/MT451: ran out of lines reading $ixp group bounds (got $k)")
        flds, _, _, _, _ = parse_endf_line(lines[li])
        for j in 1:6
            k >= ixp && break
            k += 1
            groups[k] = parse_endf_float(flds[j])
        end
        li += 1
    end
    (iza, awr, mfflg, iverf, ixmax, ixp, groups)
end

# Parse an MF3 (or MF5) section of an errorr tape into an `ixmax`-element
# cross-section vector. errorr writes one LIST record per MT:
#   CONT/LIST head : C1=ZA, C2=einc (=0 for MF3, incident-E for MF5),
#                    L1=0, L2=0, N1=ixmax, N2=0
#   payload        : ixmax data values
# (see covr.f90:786-812)
function _parse_errorr_mf3_section(lines::Vector{String}, ixmax::Integer)
    isempty(lines) && error("empty MF3 section")
    f1, _, _, _, _ = parse_endf_line(lines[1])
    einc = parse_endf_float(f1[2])
    npl  = _parse_int(f1[5])
    # Some flavours (MF=5) keep ixmax in N1; if smaller than expected, we still
    # try to read what's there — match covr.f90's permissive read-into-scr.
    n = min(npl, ixmax)
    xs = zeros(Float64, ixmax)
    k = 0
    li = 2
    while k < n
        li > length(lines) && break
        flds, _, _, _, _ = parse_endf_line(lines[li])
        for j in 1:6
            k >= n && break
            k += 1
            xs[k] = parse_endf_float(flds[j])
        end
        li += 1
    end
    (xs, einc)
end

# ----------------------------------------------------------------------
# MF33/MF34/MF35/MF40 sub-section extraction
# ----------------------------------------------------------------------

"""
    CovardResult

Result of reading + post-processing one covariance sub-section, mirroring
the `cf, xx, xy, ...` outputs of covr.f90:720-937 (`covard`).

Fields:
- `cf`     : flat row-major covariance matrix of size ixmax*ixmax
- `xx`     : per-group cross section of MT
- `xy`     : per-group cross section of MT1
- `groups` : group boundaries (= tape.groups; copied for self-containedness)
- `ixmax`  : matrix dimension
- `izero`  : 0 if cf is all-zero (null matrix), 1 otherwise (covr.f90:897)
- `izap`   : ZAP for MF40 sub-section (covr.f90:838-843); 0 for MF33
- `einc`   : MF5 incident energy (covr.f90:791); 0 for MF3
"""
struct CovardResult
    cf::Vector{Float64}
    xx::Vector{Float64}
    xy::Vector{Float64}
    groups::Vector{Float64}
    ixmax::Int
    izero::Int
    izap::Int
    einc::Float64
end

"""
    covard(tape::ErrorrTape, mat, mt, mat1, mt1; irelco=1) -> CovardResult

Port of covr.f90:720-937. Returns the (mat,mt) × (mat1,mt1) covariance
matrix on the errorr group grid, flagged with `izero` for null matrices.

For self-covariance (mat==mat1 and mt==mt1), the matrix is the same
that errorr emits in its MF33 self-section. For cross-reaction blocks,
the asymmetric LIST stream is reshaped into a full ixmax×ixmax matrix.

`irelco != 1` triggers the absolute → relative conversion at line 929
(`cf(ind+n)=cf(ind+n)/(xx(k)*xy(n))`).
"""
function covard(tape::ErrorrTape,
                mat::Integer, mt::Integer,
                mat1::Integer, mt1::Integer;
                irelco::Integer=1)::CovardResult
    mat == tape.mat ||
        error("covard: tape MAT=$(tape.mat), requested mat=$mat (only single-MAT tapes for now)")
    mat1 = mat1 == 0 ? mat : mat1
    mt1  = mt1  == 0 ? mt  : mt1

    ixmax = tape.ixmax
    groups = copy(tape.groups)

    # ---- choose covariance MF (covr.f90:820-823) ----
    mf3x = if tape.mfflg == -14
        40
    elseif tape.mfflg == -12       # MF5 = nubar / spectrum
        35
    elseif mt == 251               # mubar covariance lives in MF34
        34
    else
        33
    end

    # ---- locate cross sections (covr.f90:786-812) ----
    haskey(tape.mf3_xs, mt) ||
        error("covard: MF3/MT=$mt missing for MAT=$mat (covr.f90:788 `finds`)")
    xx = copy(tape.mf3_xs[mt])
    einc = get(tape.mf3_einc, mt, 0.0)
    xy = if mat1 == mat && mt1 == mt
        copy(xx)
    else
        haskey(tape.mf3_xs, mt1) ||
            error("covard: MF3/MT=$mt1 missing for MAT=$mat (covr.f90:811 `finds`)")
        copy(tape.mf3_xs[mt1])
    end

    # ---- locate covariance section ----
    raw = get(tape.sections, (mf3x, mt), nothing)
    raw === nothing &&
        error("covard: MF$mf3x/MT=$mt not present (covr.f90:824 `finds`)")

    cf, izap = _read_mf33_subsection(raw, mat, mt, mat1, mt1, ixmax;
                                     mf3x=mf3x, mfflg=tape.mfflg)

    # ---- post-process: zero-xs cleanup, abs→rel, null-flag (covr.f90:892-934) ----
    izero = 0
    @inbounds for k in 1:ixmax
        ind = ixmax * (k - 1)
        for n in 1:ixmax
            if xx[k] == 0.0 || xy[n] == 0.0
                cf[ind + n] = 0.0
                continue
            end
            irelco != 1 && (cf[ind + n] = cf[ind + n] / (xx[k] * xy[n]))
            cf[ind + n] != 0.0 && (izero = 1)
        end
    end

    CovardResult(cf, xx, xy, groups, ixmax, izero, izap, einc)
end

# Read one MF33 (or MF34/MF35/MF40) section's payload, find the
# sub-section matching (mat1,mt1), and unpack its per-row LIST stream into
# a flat ixmax×ixmax matrix.
#
# Section layout (matches errorr.f90 `cwrite` and ENDF-6 manual §33.2):
#   line 1 : HEAD  -- C1=ZA, C2=AWR, L1=0, L2=MTL, N1=0, N2=NL (#subsections)
#   per sub-section :
#     line k : CONT -- C1=XMF1, C2=XLFS1, L1=MAT1 (=0 → use math), L2=MT1,
#                      N1=NC, N2=NI    (NC=0 in covr-compatible output)
#                      For MF40: C2 holds IZAP (covr.f90:841 `c2h`).
#     lines  : NI  LIST records, each -- C1=0, C2=0, L1=LT, L2=lgp1 (start col),
#                      N1=NPL (#data values), N2=kl (row index 1..NI=ixmax).
#                      Loop until kl ≥ NI (== ixmax for errorr output).
#
# Returns (cf::Vector{Float64} length ixmax^2, izap::Int).
function _read_mf33_subsection(lines::Vector{String},
                               mat::Integer, mt::Integer,
                               mat1::Integer, mt1::Integer,
                               ixmax::Integer;
                               mf3x::Integer=33, mfflg::Integer=-11)
    isempty(lines) && error("MF$mf3x/MT=$mt: empty section")
    # ---- HEAD ----
    f1, _, _, _, _ = parse_endf_line(lines[1])
    nl = _parse_int(f1[6])
    nl == 0 && return (zeros(Float64, ixmax * ixmax), 0)
    li = 2

    izap = 0
    # ---- locate matching sub-section ----
    for sub in 1:nl
        li > length(lines) &&
            error("MF$mf3x/MT=$mt: ran off end while searching sub-section $sub of $nl")
        f, math, _, _, _ = parse_endf_line(lines[li])
        # Sub-section CONT: MAT1=L1; for MF34, MAT1 lives in math (cols 67-70)
        # per covr.f90:846-847 (`if mf3x.eq.34: mat1x=math; if mat1x.eq.0: mat1x=math`).
        l1   = _parse_int(f[3])
        l2   = _parse_int(f[4])
        nc   = _parse_int(f[5])
        ni   = _parse_int(f[6])
        mat1x = mf3x == 34 ? Int(math) : (l1 == 0 ? Int(math) : l1)
        mtx   = mf3x == 34 ? l1 : l2
        # MF40 carries IZAP in C2 of the matched (self) sub-section.
        sub_izap = (mfflg == -14 && mat == mat1 && mt == mt1) ?
                   round(Int, parse_endf_float(f[2])) : 0
        li += 1
        # NB: covr.f90 covard (lines 826-887) never inspects N1 of the
        # sub-section CONT — only N2 (= NI = `kgp`). For MF33/MF40, errorr
        # writes N1=0; for MF34/MT=251, errorr stores NSS=1 there. We must
        # NOT consume any extra LIST records on the strength of N1 — doing
        # so eats the first row of MF34 covariance data.
        if mat1x == mat1 && mtx == mt1
            izap = sub_izap
            cf = zeros(Float64, ixmax * ixmax)
            kgp = ni
            kl = 0
            while kl < kgp
                li > length(lines) &&
                    error("MF$mf3x/MT=$mt: ran off end inside matched sub-section")
                fl, _, _, _, _ = parse_endf_line(lines[li])
                lgp1 = _parse_int(fl[4])
                npl  = _parse_int(fl[5])
                kl   = _parse_int(fl[6])
                li += 1
                # data: NPL values, 6 per line
                k = 0
                while k < npl
                    li > length(lines) &&
                        error("MF$mf3x/MT=$mt: ran off end reading $npl LIST data values")
                    fd, _, _, _, _ = parse_endf_line(lines[li])
                    for j in 1:6
                        k >= npl && break
                        k += 1
                        v = parse_endf_float(fd[j])
                        col = lgp1 + k - 1
                        if 1 <= kl <= ixmax && 1 <= col <= ixmax
                            cf[ixmax * (kl - 1) + col] = v
                        end
                    end
                    li += 1
                end
            end
            return (cf, izap)
        else
            # not a match — skip this sub-section's NI LIST records
            kgp = ni
            kl = 0
            while kl < kgp
                li > length(lines) && error("MF$mf3x/MT=$mt: skip overran end")
                fl, _, _, _, _ = parse_endf_line(lines[li])
                npl = _parse_int(fl[5])
                kl  = _parse_int(fl[6])
                li += 1 + cld(npl, 6)
            end
        end
    end
    # No matching sub-section — return null. covr.f90:836 errors out, but the
    # plot-mode logic later flags izero=0 and skips. Mirror Fortran by raising.
    error("MF$mf3x/MT=$mt: no sub-section matching mat1=$mat1, mt1=$mt1 (covr.f90:836)")
end

# Skip exactly one LIST record (used for NC sub-subsections we don't read).
function _skip_one_list(lines::Vector{String}, li::Integer)::Int
    li > length(lines) && error("_skip_one_list: index past end")
    f, _, _, _, _ = parse_endf_line(lines[li])
    npl = _parse_int(f[5])
    li + 1 + cld(npl, 6)
end

"""
    expand_mt_list(tape::ErrorrTape; mfflg::Integer=tape.mfflg,
                   strip_mt::Integer=0, strip_mat1::Integer=0,
                   strip_mt1::Integer=0) -> Vector{Int}

Build the list of MT-numbers carried by the covariance MF on the tape,
applying the negative-`imt` strip rules (covr.f90:546-555). Used by the
expndo path when card 4 has imt(n) ≤ 0 (covr.f90:369).

The strip rules — enforced verbatim:
  * if mtn == 1: strip
  * if mtn == mstrip[i]: strip
  * if mtn == 3 and mstrip[i] == 4: strip
  * if mstrip[i] in 51..90 and mtn > mstrip[i] and mtn ≤ 90: strip
"""
function expand_mt_list(tape::ErrorrTape;
                        mfflg::Integer=tape.mfflg,
                        strip_mt::Integer=0,
                        strip_mat1::Integer=0,
                        strip_mt1::Integer=0)::Vector{Int}
    mf35 = mfflg == -12 ? 5 : 3
    mts  = sort!([mt for ((mf, mt), _) in tape.sections if mf == mf35])
    mstrip = (-strip_mt, -strip_mat1, -strip_mt1)
    keep = Int[]
    for mtn in mts
        skip = false
        for s in mstrip
            s == 0 && continue
            if mtn == 1 || mtn == s ||
               (mtn == 3 && s == 4) ||
               (51 <= s <= 90 && s < mtn <= 90)
                skip = true; break
            end
        end
        skip || push!(keep, mtn)
    end
    keep
end
