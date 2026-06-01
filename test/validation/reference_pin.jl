# reference_pin.jl -- Preflight check that the Fortran reference oracle
# (njoy-reference/) is checked out at the SHA pinned in REFERENCE_PIN.
#
# njoy-reference/ tracks upstream NJOY2016 `develop`, which rebases — so the
# oracle is a moving target unless pinned. REFERENCE_PIN (repo root) freezes it
# to one commit; this check warns LOUDLY when the local oracle has drifted off
# the pin, with the remedy (scripts/setup_reference.sh).
#
# THIS IS THE ONE PLACE WE DO NOT FAIL-LOUD: a missing REFERENCE_PIN, a missing
# njoy-reference/, or a git error emits ONE @warn and returns. It must never
# throw or abort the validation run — it is advisory only.

"""
    _repo_root_from_here() -> String

Repo root derived from this file's location (test/validation/ → ../..).
"""
_repo_root_from_here() = normpath(joinpath(@__DIR__, "..", ".."))

"""
    _read_pinned_sha(pin_path) -> Union{String,Nothing}

Parse the `SHA=...` line from REFERENCE_PIN. Returns the SHA string, or
`nothing` if the file is missing or has no parseable SHA line.
"""
function _read_pinned_sha(pin_path::AbstractString)
    isfile(pin_path) || return nothing
    for line in eachline(pin_path)
        m = match(r"^SHA=([0-9a-fA-F]+)", strip(line))
        m === nothing || return m.captures[1]
    end
    nothing
end

"""
    check_reference_pin(; repo_root=_repo_root_from_here(), quiet=false)

Verify that njoy-reference/ HEAD matches the SHA pinned in REFERENCE_PIN.

  - On match: silent (or one `@info` when `!quiet`).
  - On mismatch: one LOUD `@warn` naming both SHAs and pointing at
    `scripts/setup_reference.sh`.
  - On any setup problem (missing REFERENCE_PIN, missing njoy-reference/, git
    error): exactly ONE advisory `@warn`, then return. NEVER throws — this is
    the deliberate exception to the repo's fail-loud rule, since a validation
    run must not be aborted merely because the pin check could not run.

Returns `true` when the oracle is confirmed at the pin, `false` otherwise.
"""
function check_reference_pin(; repo_root::AbstractString = _repo_root_from_here(),
                              quiet::Bool = false)
    pin_path = joinpath(repo_root, "REFERENCE_PIN")
    ref_dir  = joinpath(repo_root, "njoy-reference")

    pinned = _read_pinned_sha(pin_path)
    if pinned === nothing
        @warn "reference pin check skipped: REFERENCE_PIN missing or unparseable; \
               run scripts/setup_reference.sh to align the oracle" pin_path
        return false
    end

    if !isdir(joinpath(ref_dir, ".git")) && !isfile(joinpath(ref_dir, ".git"))
        @warn "reference pin check skipped: njoy-reference/ missing (not a git \
               clone); run scripts/setup_reference.sh to align the oracle" ref_dir pinned
        return false
    end

    head = nothing
    try
        head = strip(read(`git -C $ref_dir rev-parse HEAD`, String))
    catch ex
        @warn "reference pin check skipped: could not read njoy-reference HEAD \
               ($(sprint(showerror, ex))); run scripts/setup_reference.sh to align the oracle" ref_dir
        return false
    end

    if head == pinned
        quiet || @info "reference oracle at pin" sha=pinned
        return true
    else
        @warn "REFERENCE ORACLE DRIFT: njoy-reference is NOT at the pinned SHA — \
               validation may compare against a moving target. Run \
               scripts/setup_reference.sh to align the oracle." pinned head
        return false
    end
end

# Allow standalone invocation: `julia test/validation/reference_pin.jl`
if abspath(PROGRAM_FILE) == @__FILE__
    check_reference_pin()
end
