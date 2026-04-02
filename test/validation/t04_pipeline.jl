# T04 Pipeline Validation: moder → reconr → errorr → groupr → errorr
#
# Test 04: U-235 (MAT=1395) covariance processing chain
#   - reconr with err=0.10 (10% tolerance for errorr test problem)
#   - errorr #1: MF33 covariance → tape23
#   - groupr: 30-group nubar (MT=452) → tape24
#   - errorr #2: MF31 covariance with groupr input → tape25
#
# Comparison targets: referenceTape23, referenceTape24, referenceTape25

using NJOY

const TEST_DIR = joinpath(@__DIR__, "../../njoy-reference/tests/04")
const INPUT_PATH = joinpath(TEST_DIR, "input")

function run_t04(; work_dir = mktempdir())
    println("=" ^ 72)
    println("T04 Pipeline: moder → reconr → errorr → groupr → errorr")
    println("=" ^ 72)

    # Run the full pipeline via run_njoy
    tapes = run_njoy(INPUT_PATH; work_dir=work_dir)

    println("\nOutput tapes:")
    for unit in sort(collect(keys(tapes.unit_to_path)))
        path = tapes.unit_to_path[unit]
        if isfile(path)
            lines = countlines(path)
            println("  tape $unit: $path ($lines lines)")
        end
    end

    # Compare output tapes against reference
    println("\n" * "=" ^ 72)
    println("Comparing output tapes against Fortran reference")
    println("=" ^ 72)

    ref_tapes = Dict(
        23 => "referenceTape23",
        24 => "referenceTape24",
        25 => "referenceTape25",
    )

    total_pass = true
    for (unit, ref_name) in sort(collect(ref_tapes))
        ref_path = joinpath(TEST_DIR, ref_name)
        if !isfile(ref_path)
            println("\n  tape$unit: SKIP (no reference)")
            continue
        end
        julia_path = haskey(tapes.unit_to_path, unit) ? tapes.unit_to_path[unit] : ""
        if !isfile(julia_path)
            println("\n  tape$unit: FAIL (not produced)")
            total_pass = false
            continue
        end

        pass = compare_tapes(julia_path, ref_path, "tape$unit")
        total_pass = total_pass && pass
    end

    println("\n" * "=" ^ 72)
    println(total_pass ? "T04: ALL TAPES MATCH" : "T04: DIFFS FOUND (see above)")
    println("=" ^ 72)
    total_pass
end

"""Compare two ENDF-format tapes line by line with numeric tolerance."""
function compare_tapes(julia_path::String, ref_path::String, label::String;
                       rel_tol=1e-5, abs_tol=1e-10)
    jlines = readlines(julia_path)
    flines = readlines(ref_path)

    println("\n  $label: Julia=$(length(jlines)) lines, Fortran=$(length(flines)) lines")

    if length(jlines) != length(flines)
        println("  $label: LINE COUNT DIFF ($(length(jlines)) vs $(length(flines)))")
    end

    nlines = min(length(jlines), length(flines))
    ndiff = 0; npass = 0; first_diff = 0

    for i in 1:nlines
        jp = rpad(jlines[i], 80)
        fp = rpad(flines[i], 80)

        # Compare columns 1-66 (data portion)
        jdata = jp[1:min(66, length(jp))]
        fdata = fp[1:min(66, length(fp))]

        if jdata == fdata
            npass += 1
            continue
        end

        # Try numeric comparison
        jfloats = extract_floats(jdata)
        ffloats = extract_floats(fdata)

        if length(jfloats) == length(ffloats) &&
           all(isclose(jfloats[k], ffloats[k]; rel_tol=rel_tol, abs_tol=abs_tol)
               for k in 1:length(jfloats))
            npass += 1
        else
            ndiff += 1
            if first_diff == 0; first_diff = i; end
            if ndiff <= 5
                println("  $label line $i:")
                println("    J: $(strip(jdata))")
                println("    F: $(strip(fdata))")
            end
        end
    end

    # Count extra lines
    extra = abs(length(jlines) - length(flines))

    if ndiff == 0 && extra == 0
        println("  $label: PERFECT ($npass/$nlines lines match)")
        return true
    else
        println("  $label: $ndiff diffs in $nlines common lines" *
                (extra > 0 ? ", $extra extra lines" : "") *
                (ndiff > 5 ? " (showing first 5)" : ""))
        return false
    end
end

function extract_floats(s::AbstractString)
    floats = Float64[]
    for m in eachmatch(r"[+-]?\d+\.?\d*[eE+\-]\d+|[+-]?\d+\.\d+", s)
        v = tryparse(Float64, replace(m.match, r"([0-9])([+-])(\d)" => s"\1e\2\3"))
        v !== nothing && push!(floats, v)
    end
    floats
end

function isclose(a::Float64, b::Float64; rel_tol=1e-5, abs_tol=1e-10)
    abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
end

# Run if invoked directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_t04()
end
