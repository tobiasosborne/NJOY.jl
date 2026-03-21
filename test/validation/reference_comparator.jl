# reference_comparator.jl -- Compare NJOY.jl output against NJOY2016 reference tapes
#
# Error metrics:
#   - Filtered relative error: |ours - ref| / |ref| only where |ref| >= ABS_FLOOR
#   - Classification: PASS (<1%), MARGINAL (1-10%), FAIL (>10%)

using NJOY
using Printf

const ABS_FLOOR = 1.0e-2
const THRESH_PASS = 0.01
const THRESH_MARGINAL = 0.10

const MT_FIELD_MAP = Dict{Int, Tuple{Symbol, String}}(
    1 => (:total, "total"), 2 => (:elastic, "elastic"),
    18 => (:fission, "fission"), 102 => (:capture, "capture"),
)

# =========================================================================
# PENDF tape parser
# =========================================================================

function _parse_endf_data_pairs(line::AbstractString)
    p = rpad(line, 80)
    pairs = Tuple{Float64,Float64}[]
    for i in 0:2
        f1, f2 = p[i*22+1:i*22+11], p[i*22+12:i*22+22]
        (strip(f1) == "" && strip(f2) == "") && break
        push!(pairs, (parse_endf_float(f1), parse_endf_float(f2)))
    end
    pairs
end

function read_ref_pendf(filename::AbstractString)
    result = Dict{Int, @NamedTuple{energies::Vector{Float64},
                                    xs::Vector{Float64}, n_points::Int}}()
    isfile(filename) || return result
    local lines
    try; lines = readlines(filename); catch; return result; end
    isempty(lines) && return result
    i = 1
    while i <= length(lines)
        length(lines[i]) < 75 && (i += 1; continue)
        p = rpad(lines[i], 80)
        local mf, mt, mat_val
        try
            mf = NJOY._parse_int(p[71:72])
            mt = NJOY._parse_int(p[73:75])
            mat_val = NJOY._parse_int(p[67:70])
        catch; i += 1; continue; end
        if mf == 3 && mt > 0 && mat_val > 0
            i += 1; i > length(lines) && break
            length(lines[i]) < 66 && (i += 1; continue)
            p2 = rpad(lines[i], 80)
            local n_points
            try; n_points = NJOY._parse_int(p2[56:66]); catch; i += 1; continue; end
            i += 2; i > length(lines) && break
            energies, xs_vals = Float64[], Float64[]
            while i <= length(lines)
                length(lines[i]) < 75 && (i += 1; continue)
                pd = rpad(lines[i], 80)
                local mfc, mtc
                try; mfc = NJOY._parse_int(pd[71:72]); mtc = NJOY._parse_int(pd[73:75])
                catch; i += 1; continue; end
                (mfc != mf || mtc != mt) && break
                try
                    for (e, x) in _parse_endf_data_pairs(pd)
                        push!(energies, e); push!(xs_vals, x)
                    end
                catch; end
                i += 1
            end
            !isempty(energies) && (result[mt] = (energies=energies, xs=xs_vals, n_points=n_points))
        else
            i += 1
        end
    end
    result
end

# =========================================================================
# Interpolation and comparison
# =========================================================================

function _interp_at(energies::Vector{Float64}, xs::Vector{Float64}, e::Float64)
    e <= energies[1] && return xs[1]
    e >= energies[end] && return xs[end]
    idx = searchsortedlast(energies, e)
    idx >= length(energies) && return xs[end]
    frac = (e - energies[idx]) / (energies[idx+1] - energies[idx])
    xs[idx] + frac * (xs[idx+1] - xs[idx])
end

struct ReactionComparison
    mt::Int; label::String; n_ref::Int; n_ours::Int
    max_rel_err::Float64; max_abs_err::Float64; worst_energy::Float64
    filt_rel_err::Float64; filt_abs_err::Float64; filt_worst_energy::Float64
    n_significant::Int; mean_rel_err::Float64; classification::Symbol
end

function compare_reaction(ref_data, our_e::Vector{Float64},
                          our_xs::Vector{Float64}, mt::Int, label::String)
    re, rxs, nr = ref_data.energies, ref_data.xs, length(ref_data.energies)
    max_r = 0.0; max_a = 0.0; we = 0.0
    fr = 0.0; fa = 0.0; fe = 0.0; ns = 0; sr = 0.0
    for k in 1:nr
        (re[k] <= 0 || !isfinite(re[k])) && continue
        xo = _interp_at(our_e, our_xs, re[k])
        ae = abs(xo - rxs[k]); den = max(abs(rxs[k]), 1e-30); rr = ae / den
        rr > max_r && (max_r = rr; max_a = ae; we = re[k])
        if abs(rxs[k]) >= ABS_FLOOR
            ns += 1; sr += rr
            rr > fr && (fr = rr; fa = ae; fe = re[k])
        end
    end
    mr = ns > 0 ? sr / ns : 0.0
    cl = fr <= THRESH_PASS ? :pass : (fr <= THRESH_MARGINAL ? :marginal : :fail)
    ReactionComparison(mt, label, nr, length(our_e), max_r, max_a, we,
                       fr, fa, fe, ns, mr, cl)
end

# =========================================================================
# Full comparison
# =========================================================================

struct ComparisonReport
    test_number::Int; ref_tape::String
    reactions::Vector{ReactionComparison}
    overall::Symbol; n_ref_total_points::Int; n_ours_points::Int; grid_ratio::Float64
end

function compare_all(result, ref_tape_path::AbstractString, test_number::Int)
    ref = read_ref_pendf(ref_tape_path)
    isempty(ref) && return ComparisonReport(test_number, basename(ref_tape_path),
        ReactionComparison[], :no_data, 0, 0, 0.0)
    reactions = ReactionComparison[]
    for (mt, (field, label)) in sort(collect(MT_FIELD_MAP))
        haskey(ref, mt) || continue
        hasfield(typeof(result), field) || continue
        push!(reactions, compare_reaction(ref[mt], result.energies,
              getfield(result, field), mt, label))
    end
    nrt = haskey(ref, 1) ? ref[1].n_points : 0
    no = hasfield(typeof(result), :energies) ? length(result.energies) : 0
    ratio = nrt > 0 ? no / nrt : 0.0
    ov = isempty(reactions) ? :no_data :
         all(r -> r.classification == :pass, reactions) ? :pass :
         any(r -> r.classification == :fail, reactions) ? :fail : :marginal
    ComparisonReport(test_number, basename(ref_tape_path), reactions, ov, nrt, no, ratio)
end

function compare_best_ref(result, ref_tapes::Dict{String,String}, test_num::Int)
    best = nothing; best_n = -1
    for (_, path) in ref_tapes
        rpt = compare_all(result, path, test_num)
        length(rpt.reactions) > best_n && (best = rpt; best_n = length(rpt.reactions))
    end
    best === nothing ? ComparisonReport(test_num, "", ReactionComparison[],
                                        :no_data, 0, 0, 0.0) : best
end

# =========================================================================
# Formatting
# =========================================================================

# Format a comparison report as a multi-line string.
function format_comparison(report::ComparisonReport)
    buf = IOBuffer()
    st = uppercase(string(report.overall))
    @printf(buf, "Test %02d [%s]:  (ref=%s, grid_ratio=%.2fx)\n",
            report.test_number, st, report.ref_tape, report.grid_ratio)
    @printf(buf, "  Points: ref=%d  ours=%d\n",
            report.n_ref_total_points, report.n_ours_points)
    for rc in report.reactions
        fl = uppercase(string(rc.classification))
        @printf(buf, "  MT=%-3d %-8s  ref=%5d  sig=%5d  filt=%.4f%%  mean=%.4f%%  [%s]\n",
                rc.mt, rc.label, rc.n_ref, rc.n_significant,
                rc.filt_rel_err * 100, rc.mean_rel_err * 100, fl)
        rc.filt_rel_err > 0 && @printf(buf, "         worst at E=%.6e eV\n", rc.filt_worst_energy)
    end
    String(take!(buf))
end

# One-line summary for a comparison report.
function summary_line(report::ComparisonReport)
    st = uppercase(string(report.overall))
    isempty(report.reactions) && return @sprintf("%-8s Test %02d  (no MF3 data)", st, report.test_number)
    w = maximum(r -> r.filt_rel_err, report.reactions)
    @sprintf("%-8s Test %02d  worst=%.3f%%  pts=%d/%d  (%.1fx)",
             st, report.test_number, w*100, report.n_ours_points,
             report.n_ref_total_points, report.grid_ratio)
end
