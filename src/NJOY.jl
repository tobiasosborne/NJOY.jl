"""
    NJOY

Julia reimplementation of the NJOY Nuclear Data Processing System.

Provides ENDF-6 format I/O, resonance cross section reconstruction,
Doppler broadening, and related nuclear data processing capabilities.
"""
module NJOY

using Printf
using Dates
using Random
using SpecialFunctions
using StaticArrays
using LinearAlgebra

# Physics constants (CODATA 2014, matching NJOY2016 phys.f90)
include("constants.jl")

# ENDF format types, I/O, and interpolation
include("endf/types.jl")
include("endf/io.jl")
include("endf/interpolation.jl")
include("endf/readers.jl")

# Resonance formalisms
include("resonances/types.jl")
include("resonances/penetrability.jl")
include("resonances/faddeeva_exact.jl")
include("resonances/faddeeva_table.jl")
include("resonances/reader.jl")
include("resonances/slbw.jl")
include("resonances/mlbw.jl")
include("resonances/reich_moore.jl")
include("resonances/sammy.jl")

# Processing modules (RECONR pipeline)
include("processing/adaptive_grid.jl")
include("processing/reconr_types.jl")
include("resonances/unresolved.jl")  # after reconr_types.jl (needs MF3Section)
include("processing/reconr_evaluator.jl")
include("processing/reconr_grid.jl")
include("processing/reconr.jl")
include("processing/pendf_writer.jl")

# Processing modules (BROADR pipeline -- Proposal B: sigma1 kernel + pipeline)
include("processing/sigma1.jl")
include("processing/broadr.jl")

# Processing modules (HEATR -- KERMA coefficients and damage energy)
include("processing/heatr.jl")

# Processing modules (THERMR -- thermal scattering cross sections)
include("processing/thermr.jl")

# Processing modules (LEAPR -- S(alpha,beta) generation from phonon DOS)
include("processing/leapr.jl")

# Processing modules (UNRESR/PURR -- unresolved resonance self-shielding)
include("processing/unresr.jl")
include("processing/purr.jl")

# Processing modules (GROUPR -- group-averaged cross sections)
include("processing/group_structures.jl")
include("processing/weight_functions.jl")
include("processing/groupr.jl")

# Processing modules (GASPR -- gas production cross sections)
include("processing/gaspr.jl")

# Processing modules (MIXR -- cross section mixing)
include("processing/mixr.jl")

# Processing modules (RESXSR -- resonance cross section file)
include("processing/resxsr.jl")

# Processing modules (MODER -- tape management and material extraction)
include("processing/moder.jl")

# Processing modules (ERRORR -- covariance processing)
include("processing/errorr.jl")

# Processing modules (COVR -- covariance visualization and library output)
include("processing/covr.jl")

# Processing modules (GAMINR -- photon interaction cross sections)
include("processing/gaminr.jl")

# Output formats (ACER -- ACE format for MCNP)
include("formats/ace_types.jl")
include("formats/ace_neutron.jl")
include("formats/ace_builder.jl")
include("formats/ace_writer.jl")

# Output formats (CCCCR -- CCCC standard interface files)
include("formats/ccccr.jl")
include("formats/ccccr_b.jl")

# Output formats (MATXSR -- MATXS interface format)
include("formats/matxsr.jl")
include("formats/matxsr_b.jl")

# Output formats (WIMSR -- WIMS-D/WIMS-E library format)
include("formats/wimsr.jl")

# Output formats (DTFR -- DTF-IV/ANISN format)
include("formats/dtfr.jl")

# Output formats (POWR -- EPRI-CELL/EPRI-CPM library format)
include("formats/powr.jl")

# Visualization (replaces plotr/viewr -- zero-dependency plot specs + renderers)
include("visualization/plotting.jl")
include("visualization/backends.jl")

# Orchestration layer -- module-level dispatch matching Fortran main.f90
include("orchestration/types.jl")
include("orchestration/input_parser.jl")
include("orchestration/auto_params.jl")
include("orchestration/pendf_io.jl")
include("orchestration/modules/moder.jl")
include("orchestration/modules/reconr.jl")
include("orchestration/modules/broadr.jl")
include("orchestration/modules/heatr.jl")
include("orchestration/modules/thermr.jl")
include("orchestration/modules/errorr.jl")
include("orchestration/modules/groupr.jl")
include("orchestration/modules/unresr.jl")
include("orchestration/modules/ccccr.jl")
include("orchestration/modules/gaminr.jl")
include("orchestration/modules/dtfr.jl")
include("orchestration/modules/matxsr.jl")
include("orchestration/modules/viewr.jl")
include("orchestration/pipeline.jl")

# Public API -- constants
export PhysicsConstants, CODATA2014

# Public API -- ENDF types
export ContRecord, ListRecord, Tab1Record, Tab2Record
export InterpolationLaw, Histogram, LinLin, LinLog, LogLin, LogLog, CoulombPen
export InterpolationTable, TabulatedFunction, MaterialId

# Public API -- ENDF I/O
export parse_endf_float, format_endf_float
export read_cont, read_list, read_tab1, read_tab2, read_head, read_tpid
export write_cont, write_list, write_tab1, write_tab2
export find_section

# Public API -- Interpolation
export terp1, interpolate, integrate

# Public API -- Resonance types
export AbstractResonanceFormalism
export SLBWParameters, MLBWParameters, ReichMooreParameters
export SAMMYParameters, SAMMYParticlePair, SAMMYSpinGroup
export AdlerAdlerParameters, UnresolvedParameters
export ResonanceRange, CrossSections

# Public API -- Penetrability
export penetrability, shift_factor, phase_shift

# Public API -- Faddeeva / Doppler broadening
export faddeeva_w, FaddeevaTable, build_faddeeva_table, quickw
export psi_chi, faddeeva_w_julia

# Public API -- MF2 reader
export MF2Data, IsotopeData, read_mf2

# Public API -- Cross section evaluation
export cross_section, cross_section_slbw, cross_section_mlbw, cross_section_rm, cross_section_sammy
export cwaven_constant, channel_radius

# Public API -- Adaptive grid reconstruction
export AdaptiveConfig, adaptive_reconstruct, round_sigfig

# Public API -- RECONR processing
export MF3Section, ENDFMaterial, PointwiseMaterial
export read_mf3_sections, build_grid, build_evaluator, linearize_one_over_v!
export merge_background!, reconstruct, reconr, sigma_mf2

# Public API -- PENDF writer
export write_pendf, write_pendf_file

# Public API -- BROADR (Doppler broadening -- Proposal B)
export f_func, f_all, h_func, h_all, h_taylor
export sigma1_at, doppler_broaden, doppler_broaden_multi, thin_xs

# Public API -- HEATR (KERMA coefficients and damage energy -- Proposal B)
export KERMAResult, LindharParams, FissionQComponents
export displacement_energy, lindhard_params, lindhard_damage
export elastic_heating, elastic_heating_aniso, elastic_damage
export capture_heating, capture_recoil, capture_damage
export fission_heating, inelastic_heating, nxn_heating
export compute_kerma, verify_kerma_sum_rule

# Public API -- THERMR (thermal scattering cross sections)
export SABData, BraggData, ThermalResult
export free_gas_xs, free_gas_kernel
export read_thermal_data, sab_kernel, sab_xs
export bragg_edges, bragg_edge_energies, build_bragg_data, structure_factor
export incoh_elastic_xs
export compute_thermal_xs, compute_thermal
export THERMR_EGRID

# Public API -- LEAPR (S(alpha,beta) generation from phonon DOS)
export PhononDOS, DiscreteOscillator, SABTable
export debye_waller_factor, phonon_expansion
export generate_sab, add_discrete_oscillators!
export debye_dos, default_alpha_grid, default_beta_grid
export sab_table_to_thermr

# Public API -- UNRESR (Bondarenko self-shielding)
export URRSpinSequence, URRStatModel
export urr_penetrability, ajku, bondarenko_xs, infinite_dilution_xs
export HWANG_QW, HWANG_QP

# Public API -- PURR (probability tables)
export ProbabilityTable
export chi2_sample, wigner_spacing
export generate_ladder, generate_ptable, bondarenko_from_ptable
export CHI2_QUANTILES

# Public API -- GROUPR (group-averaged cross sections)
export LANL_30, RRD_50, WIMS_69, VITAMINJ_175, SANDII_620, XMAS_172, ECCO_33
export GroupStructureId, IGN_LANL30, IGN_RRD50, IGN_WIMS69, IGN_SANDII620
export IGN_VITAMINJ, IGN_XMAS172, IGN_ECCO33
export get_group_structure, num_groups, validate_group_bounds, find_group
export MultiGroupXS
export constant_weight, inv_e_weight, maxwell_inv_e_fission
export vitamin_e_weight, thermal_fission_fusion, tabulated_weight
export get_weight_function
export weight_flat, weight_inv_e, weight_maxwell_fission
export group_integrate, group_average, group_average_shielded

# Public API -- GASPR (gas production cross sections)
export GasProductionResult
export gas_multiplicity, gas_yield, accumulate_gas, gas_production
export gas_production_dict, compute_gas_production

# Public API -- MIXR (cross section mixing)
export MixComponent, MixInput, union_energy_grid
export interpolate_column, mix_reactions, mix_materials

# Public API -- RESXSR (resonance cross section file)
export RESXSRecord, extract_resxs, extract_resxs_dict
export thin_resxs, write_resxs

# Public API -- MODER (tape management)
export TapeEntry, TapeDirectory, ENDFTapeSection, ENDFTapeMaterial
export read_tape_directory, extract_material, merge_tapes
export write_tpid, write_tend, moder_copy, materials, sections
export validate_tape, read_endf_tape, write_endf_tape

# Public API -- ERRORR (covariance processing)
export CovarianceBlock, CovarianceMatrix, CovarianceData
export expand_covariance_block, multigroup_covariance
export sandwich_covariance, sensitivity_jacobian
export is_symmetric, is_psd
export read_mf33, ni_covariance, nc_covariance, process_covariance

# Public API -- COVR (covariance visualization and library output)
export CorrelationMatrix
export relative_std_dev, covariance_to_correlation
export max_abs_correlation, is_valid_correlation
export format_covariance_output, process_covr, covr_summary

# Public API -- GAMINR (photon interaction cross sections)
export PhotonMultiGroupXS, PhotonGroupId
export IGG_CSEWG94, IGG_LANL12, IGG_STEINER21, IGG_STRAKER22
export IGG_LANL48, IGG_LANL24, IGG_VITAMINC36, IGG_VITAMINE38, IGG_VITAMINJ42
export CSEWG_94, LANL_12_GAMMA, STEINER_21, STRAKER_22
export LANL_48_GAMMA, LANL_24_GAMMA, VITAMINC_36, VITAMINE_38, VITAMINJ_42_GAMMA
export get_photon_group_structure, get_photon_weight_function
export photon_weight_constant, photon_weight_inv_e
export photon_group_average, photon_heating_kerma
export gaminr, group_photon_xs, PHOTON_REACTIONS

# Public API -- ACER (ACE format for MCNP)
export ACEHeader, ACENeutronTable, ACETable
export ReactionXS, EquiprobableBins, TabulatedAngular, AngularBlock
export NXS_LEN2, NXS_IZAID, NXS_NES, NXS_NTR, NXS_NR
export NXS_NTRP, NXS_NTYPE, NXS_NDNF, NXS_IS, NXS_IZ, NXS_IA
export JXS_ESZ, JXS_NU, JXS_MTR, JXS_LQR, JXS_TYR
export JXS_LSIG, JXS_SIG, JXS_LAND, JXS_AND
export JXS_LDLW, JXS_DLW, JXS_GPD, JXS_FIS, JXS_END
export ESZ_ENERGY, ESZ_TOTAL, ESZ_DISAP, ESZ_ELASTIC, ESZ_HEATING
export nxs_length, nxs_nes, nxs_ntr
export esz_energies, esz_total, esz_elastic
export format_zaid, parse_zaid, temp_to_mev, mev_to_temp
export write_ace, write_ace_table, build_ace, build_ace_from_pendf, build_xss
export write_ace_directory, ace_nes, ace_ntr

# Public API -- CCCCR (CCCC standard interface files)
export write_isotxs, write_brkoxs, write_dlayxs

# Public API -- CCCCR-B (type-safe ISOTXS writer)
export Hollerith8, ISOTXSFileIdB, ISOTXSFileControlB
export IsotopeControlB, PrincipalXSB, ScatterSubBlockB
export IsotopeDataB, ISOTXSFileB
export write_isotxs_b, build_isotxs_b

# Public API -- MATXSR (MATXS interface format)
export write_matxs, write_matxs_binary

# Public API -- MATXSR-B (type-safe MATXS writer)
export MATXSParticleB, MATXSDataTypeB, MATXSVectorB
export MATXSMatrixBlockB, MATXSSubmaterialB, MATXSMaterialB, MATXSFileB
export write_matxs_b, build_matxs_b

# Public API -- WIMSR (WIMS-D/WIMS-E library format)
export WIMSMaterial, WIMSBurnup, BurnupProduct, WIMSResonanceTable, WIMSP1Block
export WIMS_VERSION_D, WIMS_VERSION_E, WIMS_DEFAULT_NGROUPS
export write_wims

# Public API -- DTFR (DTF-IV/ANISN format)
export DTFMaterial, DTFLayout, DTFNeutronTable, DTFPhotonTable, DTFEdit
export claw_layout, write_dtf
export total_xs, ingroup_xs

# Public API -- POWR (EPRI-CELL/EPRI-CPM library format)
export EPRIFastIsotope, EPRIFastScatterMatrix
export EPRIThermalIsotope
export CPMLibrary, CPMLibraryHeader, CPMNuclideSpec, CPMResonanceData, CPMBurnupIsotope
export POWROutput, POWR_LIB_FAST, POWR_LIB_THERMAL, POWR_LIB_CPM
export POWR_NGNF, POWR_NGND_FAST, POWR_NGMIN, POWR_NGMAX
export write_powr, write_powr_fast, write_powr_cpm
export write_epri_cell, write_epri_cpm

# Public API -- Visualization (plotr/viewr replacement)
export AxisScale, LINLIN, LINLOG, LOGLIN, LOGLOG
export GridStyle, GRID_NONE, GRID_LINES, GRID_TICKS_OUT, GRID_TICKS_IN
export LineStyle, LINE_SOLID, LINE_DASHED, LINE_CHAIN_DASH, LINE_CHAIN_DOT, LINE_DOT, LINE_INVISIBLE
export MarkerShape, MARKER_NONE, MARKER_SQUARE, MARKER_OCTAGON, MARKER_TRIANGLE
export MARKER_CROSS, MARKER_EX, MARKER_DIAMOND, MARKER_CIRCLE
export CurveColor, COLOR_BLACK, COLOR_RED, COLOR_GREEN, COLOR_BLUE
export COLOR_MAGENTA, COLOR_CYAN, COLOR_BROWN, COLOR_PURPLE, COLOR_ORANGE
export CurveData, AxisSpec, PlotSpec, HeatmapSpec
export needs_autoscale, add_curve
export auto_scale_linear, auto_scale_log, resolve_axes
export extract_curve, plot_material, plot_multigroup, plot_covariance, plot_probability_table
export render_ascii, render_postscript, to_plot_recipe

# Public API -- ENDF readers (MF12 gammas, MF5 evaporation)
export read_mf12_gammas, read_mf5_evaporation

# Public API -- Orchestration (run_njoy, TapeManager)
export TapeManager, PENDFTape, PENDFMaterial, PENDFSection
export resolve, register!
export run_njoy, build_tape_manager, RunContext, final_assembly!
export read_pendf, write_pendf_tape, extract_mf3, extract_mf3_all, copy_with_modifications
export reconr_module, broadr_module, heatr_module, thermr_module, moder_module
export compute_thnmax, resolve_thnmax, select_broadr_partials
export lookup_bragg_params, extract_Z
export ModuleCall, NJOYInputDeck, parse_njoy_input
export ReconrParams, BroadrParams, HeatrParams, ThermrParams
export parse_reconr, parse_broadr, parse_heatr, parse_thermr

end # module NJOY
