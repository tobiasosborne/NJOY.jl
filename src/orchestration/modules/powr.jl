# powr module runner — Phase A scaffold (parser + dispatch only).
#
# powr produces input for the EPRI-CELL codes (GAMTAP, LIBRAR) and the
# EPRI-CPM code (CLIB). The Fortran (njoy-reference/src/powr.f90, 4902 LOC,
# 26 subroutines) splits into three independent driver paths selected by
# `lib`:
#   lib=1  →  fast()   — GAMTAP fast-neutron library    (~290 LOC + helpers)
#   lib=2  →  therm()  — LIBRAR thermal-neutron library (~370 LOC + helpers)
#   lib=3  →  cpm()    — CLIB EPRI-CPM library          (~270 LOC + sf*out helpers)
#
# Each path emits a different binary library format. None of the 84 reference
# tests in njoy-reference/tests/ exercises powr, so a real implementation
# will need standalone oracles per mode (built by running Fortran NJOY's
# powr against a custom GENDF and capturing the binary output).
#
# Phase A landing (this commit): parser recognises and validates cards 1-2,
# dispatch is wired in pipeline.jl, and `powr_module` fails LOUDLY (Rule 6)
# when called — naming the missing port phase. This makes powr "wired" in
# the project sense (23/23 modules dispatched) without producing a silent
# stub that downstream code might mistake for real output.
#
# Phase B (next): port `fast` (lib=1) end-to-end against a self-built oracle.
# Phase C: port `therm` (lib=2). Phase D: port `cpm` (lib=3).
#
# Ref: njoy-reference/src/powr.f90:63-296 (driver), :362-651 (fast),
#      :1544-1912 (therm), :1992-2263 (cpm).

"""
    powr_module(tapes::TapeManager, params::PowrParams)

Phase A scaffold: validates inputs and reports the missing implementation
phase. Will be replaced by real per-`lib` ports in Phases B-D.
"""
function powr_module(tapes::TapeManager, params::PowrParams)
    @info "powr: ngendf=$(params.ngendf) nout=$(params.nout) " *
          "lib=$(params.lib) iprint=$(params.iprint) iclaps=$(params.iclaps)"

    mode = params.lib == 1 ? "fast (GAMTAP)" :
           params.lib == 2 ? "thermal (LIBRAR)" :
                             "cpm (CLIB)"
    error(
        "powr lib=$(params.lib) [$mode] not yet ported. Phase A scaffold " *
        "only — see src/orchestration/modules/powr.jl header for the multi-" *
        "phase plan. To unblock a deck that exercises powr, implement the " *
        "matching Fortran driver in njoy-reference/src/powr.f90.")
end
