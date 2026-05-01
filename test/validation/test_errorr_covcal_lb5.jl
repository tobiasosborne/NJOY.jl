# RED→GREEN test for errorr covcal LB=5 XS·flux-weighted collapse
# (NJOY.jl-f8k Bug B — worklog/T15_covcal_mt77_diagnosis.md).
#
# T15 U-238 JENDL MF33 MT=77 (n,n7' inelastic) canary: LANL-30 group 20
# self-cov diagonal element C[20,20].
#
# Pre-fix: Julia value 4.000000e-2 — direct midpoint sample of a single
# input LB=5 bin by `expand_lb5_symmetric`, bypassing the union-grid +
# XS·flux-weighted collapse path the Fortran uses.
#
# Post-fix: 2.987998e-2 — XS·flux-weighted average across 9 union
# sub-groups covering [1.353e6, 1.738e6] eV, dragged down by a
# sub-threshold bin at 1e-5..1.4e6 (σ ≈ 7e-5 vs 5e-3 in the on-threshold
# bins).
#
# Reference algorithm — Fortran covcal (errorr.f90:1770-2417):
#   LB=5 branch (2208-2235):  cov(jh) += fvals[k(jg),k(jh)] · sig·sig1
#   per-row write (2336-2353): b[jg,jh] = cov(jh) · flx(jg)·flx(jh)
# Then covout (errorr.f90:7431-7438) collapses union→output by sum, and
# writes relative cov = Σ b / (sig_out·flx_out)(ig) · (sig_out·flx_out)(igp).
# For piecewise-constant σ within an output group this collapses to a
# σ·flx-weighted mean of the LB=5 fvals over the (ig, igp) rectangle.
#
# Run:  julia --project=. test/validation/test_errorr_covcal_lb5.jl

using Test
using NJOY

const REPO      = normpath(joinpath(@__DIR__, "..", ".."))
const T15_CACHE = joinpath(REPO, "test", "validation", "oracle_cache", "test15")
const T15_INPUT = joinpath(REPO, "njoy-reference", "tests", "15", "input")
const T15_REF   = joinpath(REPO, "njoy-reference", "tests", "15",
                           "referenceTape26")

"Run T15 groupr + second errorr (MF33) into a temp dir. Returns tape26 path."
function _run_t15_errorr_mf33(; work_dir::AbstractString = mktempdir())
    txt = read(T15_INPUT, String)
    calls = NJOY.tokenise_njoy_input(txt)
    gcall = calls[findfirst(c -> c.name === :groupr, calls)]
    gparams = NJOY.parse_groupr(gcall)
    tapes = NJOY.TapeManager(Dict{Int,String}(), work_dir)
    NJOY.register!(tapes, 21, joinpath(T15_CACHE, "run_broadr", "tape20"))
    NJOY.register!(tapes, 23, joinpath(T15_CACHE, "after_broadr.pendf"))
    NJOY.groupr_module(tapes, gparams)
    errorrs = findall(c -> c.name === :errorr, calls)
    ecall = calls[errorrs[2]]
    eparams = NJOY.parse_errorr(ecall)
    NJOY.errorr_module(tapes, eparams)
    return NJOY.resolve(tapes, abs(eparams.nout))
end

"""Parse MF=33 MT=mt self-cov sub-section into a symmetric matrix.

Reads the ENDF-format sparse row output from Julia's _write_mfcov_rows
and Fortran covout (errorr.f90:7530-7605). Self-cov is the sub-section
whose CONT header has L2 == mt."""
function _parse_mf33_self_cov(path::AbstractString, mat::Int, mt::Int, ngn::Int)
    C = zeros(Float64, ngn, ngn)
    open(path, "r") do io
        NJOY.find_section(io, 33, mt; target_mat=mat) ||
            error("_parse_mf33_self_cov: MF=33/MT=$mt not found at MAT=$mat")
        head = NJOY.read_cont(io)
        nk = Int(head.N2)
        for _ in 1:nk
            sub_hdr = NJOY.read_cont(io)
            mt2 = Int(sub_hdr.L2); mt2 == 0 && (mt2 = mt)
            sub_ngn = Int(sub_hdr.N2)
            is_self = (mt2 == mt)
            while true
                row = NJOY.read_list(io)
                ig = Int(row.N2); ig2lo = Int(row.L2); cnt = Int(row.N1)
                if is_self
                    for k in 1:cnt
                        igp = ig2lo + k - 1
                        if 1 <= ig <= ngn && 1 <= igp <= ngn
                            C[ig, igp] = row.data[k]
                        end
                    end
                end
                ig >= sub_ngn && break
            end
        end
    end
    return C
end

"Return NK (the sub-section count from the MF=33/MT=mt HEAD record's N2)."
function _read_mf33_nk(path::AbstractString, mat::Int, mt::Int)
    open(path, "r") do io
        NJOY.find_section(io, 33, mt; target_mat=mat) ||
            error("_read_mf33_nk: MF=33/MT=$mt not found at MAT=$mat")
        return Int(NJOY.read_cont(io).N2)
    end
end

"Return the sorted list of mt2 values in MF=33/MT=mt sub-section CONT.L2 fields."
function _read_mf33_mt2_list(path::AbstractString, mat::Int, mt::Int)
    mt2s = Int[]
    open(path, "r") do io
        NJOY.find_section(io, 33, mt; target_mat=mat) ||
            error("_read_mf33_mt2_list: MF=33/MT=$mt not found at MAT=$mat")
        head = NJOY.read_cont(io)
        nk = Int(head.N2)
        for _ in 1:nk
            sub_hdr = NJOY.read_cont(io)
            mt2 = Int(sub_hdr.L2); mt2 == 0 && (mt2 = mt)
            push!(mt2s, mt2)
            sub_ngn = Int(sub_hdr.N2)
            while true
                row = NJOY.read_list(io)
                ig = Int(row.N2)
                ig >= sub_ngn && break
            end
        end
    end
    return sort(mt2s)
end

@testset "errorr covcal LB=5 XS·flux-weighted collapse (T15 MT=77)" begin
    tape26 = _run_t15_errorr_mf33()
    @test isfile(tape26)

    ref_mt77 = _parse_mf33_self_cov(T15_REF, 9237, 77, 30)
    jul_mt77 = _parse_mf33_self_cov(tape26,  9237, 77, 30)

    # Sanity: the reference itself contains the Phase-50 trace canary
    # at [20, 20]. Guards against parser bugs.
    @test abs(ref_mt77[20, 20] - 2.987998e-2) < 1e-7

    # Canary: Julia must now match the reference at [20, 20].
    @info "MT=77 self-cov C[20,20]: jul=$(jul_mt77[20,20]) ref=$(ref_mt77[20,20])"
    @test abs(jul_mt77[20, 20] - 2.987998e-2) < 1e-6

    # 5×5 nonzero block rows/cols 20..24 must match ref to ~1e-5.
    max_diff = 0.0
    for ig in 20:24, igp in 20:24
        max_diff = max(max_diff, abs(jul_mt77[ig, igp] - ref_mt77[ig, igp]))
    end
    @info "MT=77 self-cov max |jul - ref| on [20:24, 20:24] = $max_diff"
    @test max_diff < 1e-5
end

# Bug A — writer NK count for cross-pairs.
# Ref: njoy-reference/src/errorr.f90:7244 (`scr(6)=nmts-ix+1`) +
# inner loop `do 180 ixp=ix,nmts` at 7254. Fortran emits one sub-section
# per (mt, mt2) with mt2 >= mt in the active reactions list, even when
# the matrix is empty (writes a 2-line zero stub via the iabort=1 path
# at lines 7350-7356, terminated at label 390).
# Pre-fix: Julia synthesises cross-pairs only for nc_derived_mts, so
# every non-NC mt collapses to NK=1 (self only).
@testset "errorr covcal Bug A — NK cross-pair stubs (T15 MT=77)" begin
    tape26 = _run_t15_errorr_mf33()
    @test isfile(tape26)

    # Reference NK and mt2 list for MT=77.
    ref_nk  = _read_mf33_nk(T15_REF, 9237, 77)
    ref_mt2 = _read_mf33_mt2_list(T15_REF, 9237, 77)
    @info "Reference T15 MT=77: NK=$ref_nk  mt2=$ref_mt2"
    @test ref_nk == 3
    @test ref_mt2 == [77, 91, 102]

    # Julia must match.
    jul_nk  = _read_mf33_nk(tape26, 9237, 77)
    jul_mt2 = _read_mf33_mt2_list(tape26, 9237, 77)
    @info "Julia     T15 MT=77: NK=$jul_nk  mt2=$jul_mt2"
    @test jul_nk == 3
    @test jul_mt2 == [77, 91, 102]
end
