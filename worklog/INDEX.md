# worklog/ — Chronological Navigation Index

This directory contains per-session debug and progress journals for NJOY.jl,
a Julia port of NJOY2016. Each file records root causes found, fixes landed,
and reference-test outcomes for one development session or closely-related
phase cluster.

**Naming is inconsistent by design** (files were created organically across
sessions). Patterns present: `T##_phase##_handoff.md` (early sessions),
`phase##_topic.md` (mid/late), `T##_topic.md` (test-focused). A future
cleanup pass could normalise to `phase##_YYYY-MM-DD_topic.md`, but this
index intentionally avoids renaming any file (HANDOFF.md references them
by their current names).

---

## Chronological Index (oldest → newest)

| Date | Phase / Test | File | Outcome (1 line) |
|------|-------------|------|-----------------|
| 2026-04-05 | Phase 1 / T03 | [T03_phase1_handoff.md](T03_phase1_handoff.md) | T03 orchestration layer scaffolded (input parsers for gaminr/dtfr/matxsr/viewr added; no tapes passing yet). |
| 2026-04-07 | Phase 2 / T03 | [T03_phase2_handoff.md](T03_phase2_handoff.md) | Fixed `_dpend` infinite loop + 3 plot-tape formatting bugs; recovered from crashed subagent session. |
| 2026-04-08 | Phase 3 / T03 | [T03_phase3_handoff.md](T03_phase3_handoff.md) | MT=621 total-heating accumulation fixed (critical); Gauss-Lobatto incoherent-scatter matrix ported. |
| 2026-04-09 | Phase 4 / T03 | [T03_phase4_handoff.md](T03_phase4_handoff.md) | 5 bugs fixed (KN factor, MT=516 NL=1, GL p-space quadrature, igzero skip, pair-production threshold). |
| 2026-04-10 | Phase 5 / T03 | [T03_phase5_handoff.md](T03_phase5_handoff.md) | Coherent angular integration rewritten in x-space with GL-6/GL-10 quadrature; ~200 tape33 lines corrected. |
| 2026-04-11 | Phase 6 / T03 | [T03_phase6_handoff.md](T03_phase6_handoff.md) | Fortran `gpanel` `rndoff`/`delta` boundary nudges applied to all 6 integration paths. |
| 2026-04-12 | Phase 7 / T03+T04 | [T03_phase7_T04_handoff.md](T03_phase7_T04_handoff.md) | **T03 BIT_IDENTICAL** (tape34/36/37); T04 tape23 NUMERIC_PASS @1e-7, tape24 @1e-5. |
| 2026-04-13 | Phase 8 | [T08_reference_test_framework.md](T08_reference_test_framework.md) | Built Fortran-faithful 84-test reference framework (`reference_test.jl`, `sweep_reference_tests.jl`, `module_coverage.jl`). |
| 2026-04-13 | Phase 9 | [T09_acer_dispatch.md](T09_acer_dispatch.md) | Wired `:acer` iopt=1 into pipeline; `acer_module` orchestration wrapper added (~150 LOC). |
| 2026-04-14 | Phase 10 | [T10_phase10_batch_dispatch.md](T10_phase10_batch_dispatch.md) | 5-fix batch: purr/leapr stubs, heatr uninit, Bragg gate, INT=0 coerce; CRASH 16→2. |
| 2026-04-15 | Phase 11 | [T11_covr_dispatch.md](T11_covr_dispatch.md) | Wired `:covr` stub into pipeline; `CovrParams` + `parse_covr` added. |
| 2026-04-15 | Phase 12 | [T12_small_batch.md](T12_small_batch.md) | plotr stub + T20 errorr fix + T43 broadr MT=1 fallback; CRASH 2→1. |
| 2026-04-15 | Phase 13 | [T13_moder_gaspr.md](T13_moder_gaspr.md) | moder extract-mode stub (iopt=1 cp) + gaspr dispatch stub; T24 CRASH→DIFFS. |
| 2026-04-15 | Phase 14 | [T14_drain_crash.md](T14_drain_crash.md) | Drained structural CRASHes: heatr plot-tape stub (T24), broadr MT=1 fallback (T60), thermr SAB fallback (T09). |
| 2026-04-17 | Phase ~15 / T49 | [T49_urr_sentinel.md](T49_urr_sentinel.md) | URR upper-boundary sentinel drop matches Fortran `lunion` ig≥ngo path; T49 41→44/46 MTs BIT_IDENTICAL. |
| 2026-04-18 | Phase ~16 | [T15_T17_mt455_crash.md](T15_T17_mt455_crash.md) | MT=455 LIST-skip crash fixed (BoundsError on empty vector); T15/T17 CRASH→DIFFS; sweep CRASH 2→0. |
| 2026-04-18 | Phase 44 | [T15_T17_errorr_size.md](T15_T17_errorr_size.md) | errorr output-group structure fixed (`ign` dispatch); 2.5 GB / 30M-line tape26 → 7953 lines. |
| 2026-04-19 | Phase 45 | [T15_groupr_auto_expand.md](T15_groupr_auto_expand.md) | groupr `3 /` auto-expand sentinel (Fortran `mtdp=-1000` path); T15 GENDF MF=3 7→39 MTs. |
| 2026-04-19 | Phase 46 | [T15_errorr_gendf_readback.md](T15_errorr_gendf_readback.md) | errorr GENDF MF3 readback fixed (`npend=0` path); T15 tape26 MF=3 0→36 MTs. |
| 2026-04-20 | Phase 47 | [T15_errorr_mf33_sparse.md](T15_errorr_mf33_sparse.md) | errorr MF33 sparse per-row emission (Fortran `ng2==0` row skip); tape26 8205→1859 lines. |
| 2026-04-20 | Phase 48 | [T15_T17_errorr_nc_expansion.md](T15_T17_errorr_nc_expansion.md) | LTY=0 NC-block expansion (LTY=0 single-NC cross-pairs); T15 tape26 1859→4178 lines. |
| 2026-04-21 | Phase 49 | [T15_T17_errorr_nc_expansion_v2.md](T15_T17_errorr_nc_expansion_v2.md) | Double-NC cross-pairs (Cov(2,4) 3→65 lines); T15 tape26 4178→4240. |
| 2026-04-21 | Phase 50 | [T15_covcal_mt77_diagnosis.md](T15_covcal_mt77_diagnosis.md) | T15 MT=77 diagnosis: Bug A (NK count) + Bug B (XS·flux-weighted union-grid) pinned; no code change. |
| 2026-04-22 | Phase 51 | [T15_covcal_lb5_weighted.md](T15_covcal_lb5_weighted.md) | LB=5 σ·flux-weighted collapse landed; MT=77 C[20,20] 0.04→0.02987998 (exact ref). |
| 2026-04-23 | Phase ~52 / T22+T80 | [T22_T80_leapr_wiring.md](T22_T80_leapr_wiring.md) | **T22 BIT_IDENTICAL** (4636/4636); T80 MISSING_TAPE→DIFFS (structural match, 76.8% S(α,β) BI). |
| 2026-04-24 | Phase 52-53 / T50 | [T50_acer_phase1_scaffolding.md](T50_acer_phase1_scaffolding.md) | T50 ACER Phases 1-6: parser fix, ACE header BI, ESZ grid, acecpe Coulomb; NUMERIC_PASS @1e-5 (143/163 BI). |
| 2026-04-24 | Phase ~38 / T38 | [T38_purr_wiring.md](T38_purr_wiring.md) | purr MT152/MT153 emission wired; T38 2849→3642/3480 lines; T21/T31 URR graceful fallback. |
| 2026-04-28 | Phase 54 | [T05_T16_covr_full_port.md](T05_T16_covr_full_port.md) | covr.f90 fully ported (2249 lines, 18 subroutines); **3/3 isolation tapes BIT_IDENTICAL** (T05/T16). |
| 2026-05-01 | Phase 55 | [T15_covcal_bug_a_nk_writer.md](T15_covcal_bug_a_nk_writer.md) | covcal Bug A fixed: errorr writer NK = all (mt,mt2) pairs; T15 tape26 4275→5964 lines (ref 5958). |
| 2026-05-02 | Phase 57 / T80 | [T80_leapr_contin_phase_a_lat_sc.md](T80_leapr_contin_phase_a_lat_sc.md) | leapr `contin` lat/sc rescaling ported; T80 75.7%→NUMERIC_PASS 99.95% @1e-5 (47 lines remain). |
| 2026-05-02 | Phase 58a / T11 | [T11_phase58a_wimsr_scaffold.md](T11_phase58a_wimsr_scaffold.md) | wimsr scaffold wired; T11 MISSING_TAPE→partial (2/3 assertions), `wimsr_module` dispatch live. |
| 2026-05-02 | Phase 58b-d / T11 | [T11_phase58_wimsr_complete.md](T11_phase58_wimsr_complete.md) | wimsr xsecs/resint/p1flx ported; T11 tape27 0→130/1169 BI (11.1%); 1039 residual diffs are ULP-class. |
| 2026-05-02 | Phase 59 | [phase59_mixr_resxsr_powr_wiring.md](phase59_mixr_resxsr_powr_wiring.md) | mixr BIT_IDENTICAL (2 oracles), resxsr BIT_IDENTICAL (CCCC binary 7968/7968 B), powr Phase-A scaffold. |
| 2026-05-02 | Phase 60 | [phase60_powr_lib1_carbon.md](phase60_powr_lib1_carbon.md) | **powr lib=1 carbon BIT_IDENTICAL** (tape50 16/16 lines; 68-group GAM-I). |
| 2026-05-02 | Phase 61-62 | [phase61_63_powr_lib1_matrices_fission.md](phase61_63_powr_lib1_matrices_fission.md) | powr lib=1 MF=6 elastic + inelastic matrices BIT_IDENTICAL (carbon 84 + 113 lines). |
| 2026-05-02 | Phase 63 | [phase61_63_powr_lib1_matrices_fission.md](phase61_63_powr_lib1_matrices_fission.md) | powr lib=1 fission path (MF=3/MT=18 + MF=6/MT=18 chi/nu) **BIT_IDENTICAL** on U-235. |
| 2026-05-03 | Phase 64-65 | [phase64_65_powr_lib1_selfshielding.md](phase64_65_powr_lib1_selfshielding.md) | powr lib=1 multi-temp + gamff f-factors **BIT_IDENTICAL** (Fe-56 ntemp=3 nsigz=4, 177×80-char tape50). |
| 2026-05-03 | Phase 66 | [phase66_powr_lib1_delayed.md](phase66_powr_lib1_delayed.md) | powr lib=1 delayed-neutron spectra (MF=5/MT=455 U-235 nfs=6) **BIT_IDENTICAL** (131-line tape50). |
| 2026-05-03 | Phase 67 | [phase67_full_sweep_post_phase66.md](phase67_full_sweep_post_phase66.md) | Fresh 84-test sweep: 2 BI (T03,T22), 2 NUMERIC_PASS (T01,T80), 78 DIFFS, 1 CRASH (T20); errorr/covr `mfflg` fix. |
| 2026-05-04 | Phase 68 | [phase68_errorr_body_mf_dispatch.md](phase68_errorr_body_mf_dispatch.md) | errorr mfcov∈{31,34} remap + MT=1 PENDF unfilter; T15+T16+T65 CRASH→runs; 47-assertion unit test. |
| 2026-05-05 | Phase 69 | [phase69_groupr_multi_T.md](phase69_groupr_multi_T.md) | Multi-temperature groupr port (T11 CRASH→DIFFS); all 4 Phase-67 CRASHes resolved. |
| 2026-05-06 | Phase 70 | [phase70_musigc_mt251_derivation.md](phase70_musigc_mt251_derivation.md) | musigc MT=251 derivation + per-mfcov MF=3 echo restriction; T15+T16+T65 CRASH→runs, sweep CRASH 4→0. |
| 2026-05-07 | Phase 71 | [phase71_rescon_diagnosis.md](phase71_rescon_diagnosis.md) | T15 MT=102 rescon diagnosis: MF=32→MF=33 sandwich-rule identified as missing piece; RED canary + warn-scaffold. |
| 2026-05-08 | Phase 72 | [phase72_mf32_reader.md](phase72_mf32_reader.md) | MF=32 reader (213 LOC, 10 ranges/317 resonances/951 params) + apply_rescon! skeleton; 57 reader tests PASS. |
| 2026-05-08 | Phase 72b | [phase72b_rescon_sandwich.md](phase72b_rescon_sandwich.md) | Full rescon sensitivity + sandwich builder (~400 LOC); MT=102 row-1 10/10 nonzero, signs correct. |
| 2026-05-16 | Phase 72c | [phase72c_rescon_iwt_fix.md](phase72c_rescon_iwt_fix.md) | iwt=6 rpxgrp/egtwtf port; C[1,1] 3.580e-4→2.658408e-4 (ref 2.658914e-4); `@test_broken` flipped GREEN. |
| 2026-05-16 | Phase 73 | [phase73_nc_xmat_and_sigma_ratio.md](phase73_nc_xmat_and_sigma_ratio.md) | Cross-material sub-section pollution fix + NC σ-ratio weighting; T15 tape26 gap +150→−68. |
| 2026-05-18 | Phase 74 | [phase74_mt1_nc_lower_ref.md](phase74_mt1_nc_lower_ref.md) | NC cross-pair to lower-MT handled (ref_j < mt); T15 tape26 gap −68→−6; MT=1/mt2=2 geometry BI. |
| 2026-05-19 | Phase 75 | [phase75_urr_rescon.md](phase75_urr_rescon.md) | URR rescon (ggunr1/rpxunr) ported; T15 tape26 gap −6→**0** (5958=ref); all MF=33 row layouts match. |
| 2026-05-28 | Phase 77 | [phase77_broadr_acer_sweep.md](phase77_broadr_acer_sweep.md) | broadr O(n²)→O(n) (~26× speedup, BI); T50 BI 143→163/163; T20 CRASH→DIFFS; fresh sweep 3 BI, 0 CRASH. |
| 2026-05-29 | Phase 78 | [phase78_acer_charged_arc.md](phase78_acer_charged_arc.md) | **4 NEW BIT_IDENTICAL: T52, T62, T61, T53**; ptlegc/coul, charged NTR>0, iopt=7 xsdir, acelcp particle-type blocks. |
