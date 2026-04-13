# Phase 9 — Dispatch `acer` in `run_njoy`

## Date: 2026-04-13 (continuing from Phase 8 sweep baseline)

## Summary

Wired `:acer` into the pipeline dispatcher. The low-level ACE format work
(types, builders, writer) was already ported in `src/formats/ace_*.jl`
(~846 lines); the missing piece was the orchestration wrapper that glues
the input-deck parser to those builders and writes output tapes.

## Files

**Added:** `src/orchestration/modules/acer.jl` (~150 lines).

**Modified:**
- `src/NJOY.jl` — `include("orchestration/modules/acer.jl")`.
- `src/orchestration/pipeline.jl` — `:acer` branch in the dispatch loop.

## Behavior

### iopt=1 (fast continuous-energy neutron)

Full path:
1. Resolve PENDF input tape (`params.npendf`) via `TapeManager`.
2. `read_pendf(path)` → `PENDFTape`.
3. `extract_mf3_all(tape, mat)` → `Dict{Int,(energies, xs)}`.
4. Pick MT=1's energy grid as the master (reconr/broadr output has all
   partials on the lunion grid — verified by grid-length check with
   linear-interp fallback if they disagree).
5. Assemble an `xs_matrix` column-per-MT in canonical order
   (1, 2, 18, 102, then others sorted).
6. Construct `PointwiseMaterial(mat, master_e, xs_matrix, mt_list)`.
7. Extract ZA from MF1/MT451 HEAD record for the ZAID header string
   (differs from `mat` for most isotopes).
8. `build_ace_from_pendf(pm; suffix, temp_kelvin, mat_id)` →
   `ACENeutronTable`.
9. `write_ace_table(io, table)` to `tapes[params.nace]`.
10. `write_ace_directory(io, table, nxs)` to `tapes[params.ndir]`.

### iopt≠1 (thermal, photoatomic, dosimetry, continuous-energy photon, etc.)

`@warn` + write **empty** stub tapes at `params.nace` / `params.ndir` so
downstream modules (viewr, moder) don't `SystemError` when opening them.
The sweep then classifies these as `STRUCTURAL_FAIL` (line count
mismatch) rather than `CRASH`, which is useful signal: we know the deck
executed end-to-end, we just haven't produced real data for that iopt.

## Verification

T14 and T50 transitioned from `CRASH: SystemError: opening .../tape33`
to `DIFFS` — pipeline runs through both acer calls + viewr. Output
doesn't byte-match Fortran yet because:
- T14 second acer is iopt=7 (photoatomic) → stub.
- T14 first acer writes xsdir in neutron format (`.80c`) but the
  reference expects photoatomic `.00h` format.
- Those would only match once iopt=7 lands.

Full sweep not yet re-run at time of commit.

## Missing for bit-identical acer output

Tracked here so a future session can pick up directly.

### iopt modes still stubbed

| iopt | Description | Count of reference tests affected |
|------|-------------|-----------------------------------|
| 2    | thermal S(α,β) ACE | T50–T54 (via leapr chain), T62 |
| 3    | dosimetry | T48 |
| 4    | incident photon continuous-energy | ??? |
| 5    | incident proton continuous-energy | T14 (second call) |
| 6    | ??? |  |
| 7    | photoatomic | T14 (second call), T50+ |
| 8    | dosimetry (alternate form) |  |
| 9    | charged-particle ACE | T52, T54 |

### iopt=1 tuning needed

Even iopt=1 won't immediately byte-match the Fortran reference because:
- `format_zaid` needs verification against `acer.f90:acelod` conventions.
- AWR in the ACE header is computed as `za % 1000` in
  `build_ace_from_pendf` — should be read from MF1/MT451 AWR field.
- Angular and energy distributions (MF4, MF5, MF6) aren't carried
  through: `build_ace_from_pendf` only uses MF3 cross sections. Fortran
  acer's full output includes MU, PDF, CDF tables for every reaction.
  These are substantial (MF6 tabulated distributions → ACE LAW-3/LAW-4
  tables).

## Sweep impact — estimate pending re-sweep

Conservative: **~15 tests CRASH → DIFFS/STRUCTURAL_FAIL**. These tests
previously failed to open tape33/34 because acer never wrote them. Now
they get written (empty for iopt≠1, real but format-incomplete for
iopt=1). The CRASH bucket drops; the DIFFS/STRUCTURAL_FAIL bucket grows
by the same amount. No new BIT_IDENTICAL or NUMERIC_PASS expected until
the iopt=1 header/AWR fixes land.

## Recommendations (priority-ordered)

### Immediate (before next sweep, ideally today)

1. **Dispatch `purr`** in `pipeline.jl`. PURR (probability-table
   generation for MCNP) is structurally similar to acer — reads PENDF,
   writes a tape with URR probability tables. `purr_module` doesn't
   exist yet; source code for the PURR algorithm may or may not be
   ported. If the core code isn't ported, write a stub module that
   writes an empty output tape so T35–T42, T28, T34, T63, T65, T71, T72
   move from CRASH → STRUCTURAL_FAIL. Conservative: **~10 tests
   unblock**.

2. **Dispatch `leapr`** in `pipeline.jl`. LEAPR generates S(α,β)
   thermal scattering data. `src/processing/leapr_*` — **check if
   ported before writing a stub**. If ported, wire up; if not, stub.
   Unblocks T09, T22, T23, T33, T80 plus several downstream thermr
   chains. Conservative: **~5 tests**.

3. **Re-sweep** (2h) after all three dispatches land. Compare new
   REFERENCE_SWEEP.md to `REFERENCE_SWEEP_opaque_baseline.md` to measure
   exact impact.

### Short-term (next session)

4. **Fix `UndefVarError: h` in heatr.jl** — 6 tests (T08, T13, T21,
   T26, T49, T79) crash on this same variable-scope bug. Should be a
   one-liner.

5. **Add Bragg lattice entries** for MAT=1, 7, 15, 53, 58 in
   `src/orchestration/auto_params.jl::BRAGG_LATTICE_PARAMS`. Materials:
   H-1 (MAT=1 is nonsense — probably Be-9 at MAT=1 in some library),
   Li-7 (MAT=7), N-15 (MAT=15), Cr-53 (MAT=53), Fe-58 (MAT=58). 7 tests
   unblock (T25, T32, T67, T68, T69, T70, T74).

6. **Accept INT=0 in MF TAB1 reader** for U-238 JENDL-3.3 (T15, T17).
   This is a real bug: the ENDF standard says INT=0 is invalid, but
   JENDL-3.3 emits it for some sections, and Fortran NJOY silently
   accepts it. Check `src/endf/io.jl` for `InterpolationLaw` enum and
   add INT=0 handling (treat as INT=2 linear?).

### Medium-term

7. **Tune iopt=1 ACE output** to byte-match references for the simpler
   tests (T43/T44 H-1, T48 H-3): fix ZAID formatting, AWR source,
   header date convention. Then the ~15 acer tests start earning
   `NUMERIC_PASS`/`BIT_IDENTICAL` statuses.

8. **Implement iopt=2/7** for T14/T50–T54 thermal and photoatomic ACE
   — requires more substantial work: S(α,β) ACE tables (iopt=2) and
   photoatomic tables (iopt=7) each have their own XSS layouts.

9. **Grind the 12 existing DIFFS cases** to reduce their diff counts.
   Most are a few lines from passing at 1e-5 or 1e-3.

## Framework health

The Phase 8 framework is doing its job: every module added to the
dispatcher immediately becomes testable under Fortran-identical
conditions, with classified output. Adding acer was a one-commit,
~150-line change with immediate sweep-visible results. The same
pattern applies to purr and leapr.

Re-sweep cost is 2h each time. For this reason, batching multiple
dispatches per sweep is more efficient than iterating one-at-a-time
unless the changes interact (they don't — acer, purr, leapr are
independent pipeline entries).
