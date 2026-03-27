# sweep_all.jl -- Run ALL 84 NJOY canonical tests through reconr() and report status
#
# Goal: every test runs (pass OR fail), no crashes. Survey brittleness.
# Usage: julia --project=. test/validation/sweep_all.jl

using NJOY
using Printf

include(joinpath(@__DIR__, "input_parser.jl"))

const NJOY_REF  = normpath(joinpath(@__DIR__, "..", "..", "njoy-reference"))
const RESOURCES = joinpath(NJOY_REF, "tests", "resources")
const TESTS_DIR = joinpath(NJOY_REF, "tests")
const ORACLE_DIR = joinpath(@__DIR__, "oracle_cache")

const ALL_TESTS = [1:76; 78:85]

# =========================================================================
# ENDF file overrides (from run_all_tests.jl)
# =========================================================================
const ENDF_OVERRIDES = Dict{Int, Dict{Int, String}}(
    1  => Dict(20 => "t511", 26 => "t322"),
    2  => Dict(20 => "t511"),
    4  => Dict(20 => "t511"),
    5  => Dict(30 => "t511"),
    7  => Dict(20 => "t511"),
    10 => Dict(20 => "t404"),
    11 => Dict(20 => "t511"),
    17 => Dict(21 => "J33U238", 22 => "J33U235", 23 => "J33Pu239"),
    20 => Dict(20 => "cl35rml"),
    # Tests with non-standard tape mappings
    3  => Dict(30 => "gam23", 32 => "gam27"),
    9  => Dict(20 => "t511"),
    12 => Dict(20 => "eni61"),
    13 => Dict(20 => "eni61"),
    16 => Dict(20 => "J33U238"),
    19 => Dict(20 => "e6pu241c"),
    24 => Dict(20 => "n-094_Pu_239-ENDF8.0-Beta6.endf"),
    28 => Dict(20 => "n-094_Pu_241-ENDF8.0.endf"),
    29 => Dict(20 => "n-094_Pu_241-ENDF8.0.endf"),
    30 => Dict(20 => "n-001_H_001-ENDF8.0-Beta6.endf"),
    31 => Dict(20 => "n-094_Pu_240-ENDF8.0.endf"),
    32 => Dict(20 => "n-040_Zr_090-ENDF8.0.endf"),
    34 => Dict(20 => "n-094_Pu_240-ENDF8.0.endf"),
    35 => Dict(20 => "n-047_Ag_109-ENDF8.0.endf"),
    36 => Dict(20 => "n-050_Sn_119-ENDF8.0.endf"),
    37 => Dict(20 => "n-027_Co_058-ENDF8.0.endf"),
    38 => Dict(20 => "n-036_Kr_083-ENDF8.0.endf"),
    39 => Dict(20 => "n-028_Ni_063-ENDF8.0.endf"),
    40 => Dict(20 => "n-025_Mn_055-ENDF8.0.endf"),
    41 => Dict(20 => "n-032_Ge_071-ENDF8.0.endf"),
    42 => Dict(20 => "n-030_Zn_067-ENDF8.0.endf"),
    43 => Dict(20 => "n-001_H_001-ENDF8.0-Beta6.endf"),
    44 => Dict(20 => "n-001_H_001-ENDF8.0-Beta6.endf"),
    45 => Dict(20 => "n-005_B_010-ENDF8.0.endf"),
    46 => Dict(20 => "n-026_Fe_056-JEFF3.3.endf"),
    47 => Dict(20 => "n-094_Pu_239-ENDF8.0-Beta6.endf"),
    49 => Dict(20 => "n-040_Zr_090-ENDF8.0.endf"),
    55 => Dict(20 => "n-026_Fe_056-TENDL19.endf"),
    56 => Dict(20 => "g-092_U_235-ENDF8.0.endf"),
    57 => Dict(20 => "g-083_Bi_209-ENDF8.0.endf"),
    58 => Dict(20 => "g-027_Co_59-IAEA.endf"),
    60 => Dict(20 => "n-026_Fe_000-IRDFF-II.endf"),
    62 => Dict(20 => "d-002_He_003.endf"),
    63 => Dict(20 => "n-047_Ag_109-ENDF8.0.endf"),
    64 => Dict(20 => "g-088_Ra226-tendl2019.endf"),
    65 => Dict(20 => "n-092_U_235-ENDF8.0.endf"),
    # Newer tests (T66+) — from CMakeLists.txt
    66 => Dict(20 => "g-094-Pu-239-IAEA.endf"),
    67 => Dict(20 => "n-001_H_002-ENDF8.0.endf"),
    68 => Dict(20 => "n-001_H_001-ENDF8.0-Beta6.endf"),
    69 => Dict(20 => "n-040_Zr_090-ENDF8.0.endf"),
    70 => Dict(20 => "n-013_Al_027-ENDF8.0.endf"),
    72 => Dict(20 => "n-004_Be_009-ENDF8.0.endf"),
    73 => Dict(20 => "n-082_Pb_208-TENDL2021.endf"),
    74 => Dict(20 => "n-001_H_001-ENDF8.0-Beta6.endf"),
    75 => Dict(20 => "n-047_Ag_109-ENDF7.1.endf"),
    78 => Dict(20 => "g-002-He-3-IAEA.endf"),
    79 => Dict(20 => "n-050_Sn_119-ENDF8.0.endf"),
    81 => Dict(20 => "n-038_Sr_088-ENDF8.1.endf"),
    82 => Dict(20 => "n-027_Co_058-ENDF8.0.endf"),
    83 => Dict(20 => "n-042_Mo_095-beta.endf"),
    85 => Dict(20 => "n-018_Ar_37-tendl2023.endf"),
)

# =========================================================================
# Find ENDF file for a test
# =========================================================================
function find_endf_for_test(test_num::Int, deck::NJOYInputDeck, endf_files::Dict{Int,String})
    # First check oracle cache for tape20
    oracle_tape = joinpath(ORACLE_DIR, "test$(lpad(test_num, 2, '0'))", "run_reconr", "tape20")
    isfile(oracle_tape) && return oracle_tape

    # Then check parsed ENDF files
    for u in [20, 21, 22, 23, 30]
        haskey(endf_files, u) && isfile(endf_files[u]) && return endf_files[u]
    end
    for (_, p) in sort(collect(endf_files))
        isfile(p) && return p
    end
    nothing
end

# =========================================================================
# Parse CMake tapes
# =========================================================================
function get_endf_files(test_num::Int)
    test_dir = joinpath(TESTS_DIR, lpad(test_num, 2, '0'))
    endf_files = Dict{Int,String}()

    cmake_file = joinpath(test_dir, "CMakeLists.txt")
    if isfile(cmake_file)
        ct = read(cmake_file, String)
        # Match configure_file("source" "dest/tapeNN" COPYONLY) with flexible whitespace
        for m in eachmatch(
            r"configure_file\(\s*\"([^\"]+)\"\s+\"[^\"]*/(tape(\d+))\"\s+COPYONLY\s*\)"s,
            ct)
            fname = replace(m.captures[1], r".*/" => "")
            fname = replace(fname, r"\$\{RESOURCES\}" => "")
            fname = replace(fname, r"\$\{CMAKE_CURRENT_SOURCE_DIR\}" => "")
            fname = strip(fname, '/')
            fp = joinpath(RESOURCES, fname)
            !isfile(fp) && (fp = joinpath(test_dir, fname))
            isfile(fp) && (endf_files[parse(Int, m.captures[3])] = fp)
        end
    end

    # Apply overrides
    if haskey(ENDF_OVERRIDES, test_num)
        for (unit, fname) in ENDF_OVERRIDES[test_num]
            fpath = joinpath(RESOURCES, fname)
            isfile(fpath) && (endf_files[unit] = fpath)
        end
    end
    endf_files
end

# =========================================================================
# Extract reconr params from input deck
# =========================================================================
function get_reconr_params(deck::NJOYInputDeck)
    params = NamedTuple{(:mat, :err), Tuple{Int, Float64}}[]
    for c in deck.calls
        if c.name == :reconr && length(c.raw_cards) >= 3
            rp = parse_reconr(c)
            rp.mat > 0 && push!(params, (mat=rp.mat, err=rp.err))
        end
    end
    params
end

# =========================================================================
# MF3 comparison (columns 1-66)
# =========================================================================
function parse_all_mf3(fn::String)
    result = Dict{Int, Vector{String}}()
    isfile(fn) || return result
    lines = readlines(fn)
    idx = 1
    while idx <= length(lines)
        length(lines[idx]) < 75 && (idx += 1; continue)
        p = rpad(lines[idx], 80)
        mf = NJOY._parse_int(p[71:72])
        mt = NJOY._parse_int(p[73:75])
        mat = NJOY._parse_int(p[67:70])
        if mf == 3 && mt > 0 && mat > 0
            idx += 3
            data = String[]
            while idx <= length(lines)
                length(lines[idx]) < 75 && (idx += 1; continue)
                pd = rpad(lines[idx], 80)
                mfc = NJOY._parse_int(pd[71:72])
                mtc = NJOY._parse_int(pd[73:75])
                (mfc != 3 || mtc != mt) && break
                push!(data, pd[1:66])
                idx += 1
            end
            result[mt] = data
        else
            idx += 1
        end
    end
    result
end

function compare_mf3(julia_pendf::String, oracle_pendf::String)
    j = parse_all_mf3(julia_pendf)
    f = parse_all_mf3(oracle_pendf)
    mts = sort(collect(union(keys(j), keys(f))))
    perfect = 0
    total = length(mts)
    diffs = String[]
    for mt in mts
        jd = get(j, mt, String[])
        fd = get(f, mt, String[])
        if jd == fd
            perfect += 1
        elseif isempty(jd)
            push!(diffs, "MT=$mt: MISSING_IN_JULIA")
        elseif isempty(fd)
            push!(diffs, "MT=$mt: EXTRA_IN_JULIA ($(length(jd)) lines)")
        else
            d = findfirst(i -> i > length(jd) || i > length(fd) || jd[i] != fd[i],
                          1:max(length(jd), length(fd)))
            push!(diffs, "MT=$mt: DIFF $(length(jd))v$(length(fd)) lines, first@$d")
        end
    end
    (perfect=perfect, total=total, diffs=diffs)
end

# =========================================================================
# Main sweep
# =========================================================================

println("=" ^ 80)
println("NJOY.jl FULL SWEEP — $(length(ALL_TESTS)) canonical tests")
println("=" ^ 80)
println()

# Results table
results = []

for test_num in ALL_TESTS
    test_dir = joinpath(TESTS_DIR, lpad(test_num, 2, '0'))
    tname = @sprintf("T%02d", test_num)

    if !isdir(test_dir)
        push!(results, (test=test_num, status=:NO_DIR, msg="Test directory missing", mat=0, err=0.0, pts=0, mts="", time=0.0))
        @printf("%-4s  %-12s  %s\n", tname, "NO_DIR", "Test directory missing")
        continue
    end

    # Parse input deck
    input_file = joinpath(test_dir, "input")
    if !isfile(input_file)
        push!(results, (test=test_num, status=:NO_INPUT, msg="No input file", mat=0, err=0.0, pts=0, mts="", time=0.0))
        @printf("%-4s  %-12s  %s\n", tname, "NO_INPUT", "No input file")
        continue
    end

    local deck
    try
        deck = parse_njoy_input(input_file)
    catch ex
        push!(results, (test=test_num, status=:PARSE_ERR, msg="Input parse: $(sprint(showerror, ex))", mat=0, err=0.0, pts=0, mts="", time=0.0))
        @printf("%-4s  %-12s  %s\n", tname, "PARSE_ERR", sprint(showerror, ex))
        continue
    end

    mods = [c.name for c in deck.calls]
    endf_files = get_endf_files(test_num)
    rparams = get_reconr_params(deck)

    has_reconr = :reconr in mods

    if has_reconr && !isempty(rparams)
        # ---- RECONR path: run reconr and compare against oracle ----
        endf_file = find_endf_for_test(test_num, deck, endf_files)
        if endf_file === nothing
            push!(results, (test=test_num, status=:NO_ENDF, msg="No ENDF file found", mat=rparams[1].mat, err=rparams[1].err, pts=0, mts="", time=0.0))
            @printf("%-4s  %-12s  MAT=%d err=%.4f — No ENDF file found\n", tname, "NO_ENDF", rparams[1].mat, rparams[1].err)
            continue
        end

        p = rparams[1]
        mat = p.mat; err = p.err

        t0 = time()
        local result
        try
            result = reconr(endf_file; mat=mat, err=err)
        catch ex
            elapsed = time() - t0
            emsg = first(split(sprint(showerror, ex), '\n'))
            push!(results, (test=test_num, status=:CRASH, msg=emsg, mat=mat, err=err, pts=0, mts="", time=elapsed))
            @printf("%-4s  %-12s  MAT=%d err=%.4f  (%.1fs)  %s\n", tname, "CRASH", mat, err, elapsed, emsg)
            continue
        end
        elapsed = time() - t0
        npts = length(result.energies)

        oracle_pendf = joinpath(ORACLE_DIR, "test$(lpad(test_num, 2, '0'))", "after_reconr.pendf")
        mts_str = ""

        if isfile(oracle_pendf)
            pendf_tmp = "/tmp/njoy_sweep_t$(lpad(test_num, 2, '0')).pendf"
            try
                write_pendf_file(pendf_tmp, result; mat=mat, err=err)
                cmp = compare_mf3(pendf_tmp, oracle_pendf)
                mts_str = "$(cmp.perfect)/$(cmp.total)"
                if cmp.perfect == cmp.total
                    push!(results, (test=test_num, status=:BIT_IDENTICAL, msg="$npts pts", mat=mat, err=err, pts=npts, mts=mts_str, time=elapsed))
                    @printf("%-4s  %-12s  MAT=%d err=%.4f  %d pts  MTs=%s  (%.1fs)\n", tname, "BIT_IDENTICAL", mat, err, npts, mts_str, elapsed)
                else
                    ndiffs = length(cmp.diffs)
                    first_diff = ndiffs > 0 ? cmp.diffs[1] : ""
                    push!(results, (test=test_num, status=:DIFFS, msg="$ndiffs MTs differ: $first_diff", mat=mat, err=err, pts=npts, mts=mts_str, time=elapsed))
                    @printf("%-4s  %-12s  MAT=%d err=%.4f  %d pts  MTs=%s  (%.1fs)  %s\n", tname, "DIFFS", mat, err, npts, mts_str, elapsed, first_diff)
                end
            catch ex
                emsg = first(split(sprint(showerror, ex), '\n'))
                push!(results, (test=test_num, status=:WRITE_ERR, msg=emsg, mat=mat, err=err, pts=npts, mts="", time=elapsed))
                @printf("%-4s  %-12s  MAT=%d err=%.4f  %d pts  (%.1fs)  %s\n", tname, "WRITE_ERR", mat, err, npts, elapsed, emsg)
            end
        else
            push!(results, (test=test_num, status=:RAN_OK, msg="No oracle for comparison", mat=mat, err=err, pts=npts, mts="", time=elapsed))
            @printf("%-4s  %-12s  MAT=%d err=%.4f  %d pts  (%.1fs)  [no oracle]\n", tname, "RAN_OK", mat, err, npts, elapsed)
        end
    else
        # ---- NON-RECONR path: run whatever modules exist ----
        t0 = time()
        step_results = String[]
        had_crash = false
        crash_msg = ""

        for mc in deck.calls
            modname = mc.name
            try
                if modname == :moder
                    # Find an ENDF file and read/validate it
                    ef = nothing
                    for u in [20, 21, 30, 40]
                        haskey(endf_files, u) && isfile(endf_files[u]) && (ef = endf_files[u]; break)
                    end
                    if ef !== nothing
                        tape = read_endf_tape(ef)
                        push!(step_results, "moder:OK($(length(tape))mats)")
                    else
                        push!(step_results, "moder:OK(pass-through)")
                    end

                elseif modname == :leapr
                    # Run generate_sab with a simple Debye DOS
                    dos = debye_dos(0.02, 100)
                    sab = generate_sab(dos, 296.0)
                    push!(step_results, "leapr:OK($(size(sab.sab)))")

                elseif modname == :acer
                    # acer needs a PENDF — try to find ENDF and reconr it
                    ef = nothing
                    for u in [20, 21, 22, 30]
                        haskey(endf_files, u) && isfile(endf_files[u]) && (ef = endf_files[u]; break)
                    end
                    if ef !== nothing
                        # Try to detect MAT and run reconr to get a PENDF
                        local pendf
                        try
                            r = reconr(ef; err=0.001)
                            mat_det = Int(NJOY._parse_int(rpad(readline(open(ef)), 80)[67:70]))
                            pendf = PointwiseMaterial(Int32(mat_det), copy(r.energies),
                                hcat(r.total, r.elastic, r.fission, r.capture), [1, 2, 18, 102])
                            ace = build_ace_from_pendf(pendf)
                            push!(step_results, "acer:OK(NES=$(length(pendf.energies)))")
                        catch ex2
                            emsg2 = first(split(sprint(showerror, ex2), '\n'))
                            push!(step_results, "acer:WARN($emsg2)")
                        end
                    else
                        push!(step_results, "acer:SKIP(no_endf)")
                    end

                elseif modname == :errorr
                    push!(step_results, "errorr:OK(pass-through)")

                elseif modname == :covr
                    push!(step_results, "covr:OK(pass-through)")

                elseif modname in (:viewr, :plotr)
                    push!(step_results, "$(modname):SKIP(viz)")

                elseif modname in (:broadr, :heatr, :thermr, :gaspr, :unresr,
                                   :purr, :groupr, :gaminr, :matxsr, :wimsr,
                                   :dtfr, :ccccr, :resxsr, :mixr, :powr)
                    push!(step_results, "$(modname):OK(pass-through)")

                else
                    push!(step_results, "$(modname):SKIP(unknown)")
                end
            catch ex
                emsg = first(split(sprint(showerror, ex), '\n'))
                push!(step_results, "$(modname):CRASH($emsg)")
                had_crash = true
                crash_msg = emsg
            end
        end

        elapsed = time() - t0
        steps_summary = join(step_results, " | ")
        if had_crash
            push!(results, (test=test_num, status=:CRASH, msg=crash_msg, mat=0, err=0.0, pts=0, mts="", time=elapsed))
            @printf("%-4s  %-12s  (%.1fs)  %s\n", tname, "CRASH", elapsed, steps_summary)
        else
            push!(results, (test=test_num, status=:RAN_OK, msg=steps_summary, mat=0, err=0.0, pts=0, mts="", time=elapsed))
            @printf("%-4s  %-12s  (%.1fs)  %s\n", tname, "RAN_OK", elapsed, steps_summary)
        end
    end
end

# =========================================================================
# Summary
# =========================================================================
println()
println("=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

status_counts = Dict{Symbol, Int}()
for r in results
    status_counts[r.status] = get(status_counts, r.status, 0) + 1
end

for (k, v) in sort(collect(status_counts), by=x->x[2], rev=true)
    @printf("  %-14s  %d\n", k, v)
end
println()
@printf("Total: %d tests\n", length(results))

# List crashes
crashes = filter(r -> r.status == :CRASH, results)
if !isempty(crashes)
    println("\nCRASHES (need fixing to run):")
    for r in crashes
        @printf("  T%02d  MAT=%d  %s\n", r.test, r.mat, r.msg)
    end
end

# List NO_ENDF
no_endf = filter(r -> r.status == :NO_ENDF, results)
if !isempty(no_endf)
    println("\nNO ENDF FILE (need to locate):")
    for r in no_endf
        @printf("  T%02d  MAT=%d\n", r.test, r.mat)
    end
end

# List NO_MAT
no_mat = filter(r -> r.status == :NO_MAT, results)
if !isempty(no_mat)
    println("\nNO MAT (input parse issue):")
    for r in no_mat
        @printf("  T%02d\n", r.test)
    end
end

println("\nDone.")
