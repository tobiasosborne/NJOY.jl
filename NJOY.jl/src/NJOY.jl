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

# Resonance formalisms
include("resonances/types.jl")
include("resonances/penetrability.jl")
include("resonances/faddeeva_exact.jl")
include("resonances/faddeeva_table.jl")
include("resonances/reader.jl")
include("resonances/slbw.jl")
include("resonances/mlbw.jl")
include("resonances/reich_moore.jl")

# Processing modules (RECONR pipeline)
include("processing/adaptive_grid.jl")
include("processing/reconr_types.jl")
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

# Processing modules (UNRESR/PURR -- unresolved resonance self-shielding)
include("processing/unresr.jl")
include("processing/purr.jl")

# Processing modules (GROUPR -- group-averaged cross sections)
include("processing/group_structures.jl")
include("processing/weight_functions.jl")
include("processing/groupr.jl")

# Processing modules (MODER -- tape management and material extraction)
include("processing/moder.jl")

# Processing modules (ERRORR -- covariance processing)
include("processing/errorr.jl")

# Output formats (ACER -- ACE format for MCNP)
include("formats/ace_types.jl")
include("formats/ace_neutron.jl")
include("formats/ace_builder.jl")
include("formats/ace_writer.jl")

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
export cross_section, cross_section_slbw, cross_section_mlbw, cross_section_rm
export cwaven_constant, channel_radius

# Public API -- Adaptive grid reconstruction
export AdaptiveConfig, adaptive_reconstruct, round_sigfig

# Public API -- RECONR processing
export MF3Section, ENDFMaterial, PointwiseMaterial
export read_mf3_sections, build_grid, build_evaluator
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
export bragg_edges, bragg_edge_energies, build_bragg_data
export incoh_elastic_xs
export compute_thermal_xs, compute_thermal
export THERMR_EGRID

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
export LANL_30, WIMS_69, VITAMINJ_175, SANDII_620, XMAS_172, ECCO_33
export GroupStructureId, IGN_LANL30, IGN_WIMS69, IGN_SANDII620
export IGN_VITAMINJ, IGN_XMAS172, IGN_ECCO33
export get_group_structure, num_groups, validate_group_bounds, find_group
export MultiGroupXS
export constant_weight, inv_e_weight, maxwell_inv_e_fission
export vitamin_e_weight, thermal_fission_fusion, tabulated_weight
export get_weight_function
export weight_flat, weight_inv_e, weight_maxwell_fission
export group_integrate, group_average, group_average_shielded

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

end # module NJOY
