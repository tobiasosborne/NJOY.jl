"""
    NJOY

Julia reimplementation of the NJOY Nuclear Data Processing System.

Provides ENDF-6 format I/O, resonance cross section reconstruction,
Doppler broadening, and related nuclear data processing capabilities.
"""
module NJOY

using Printf
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
include("resonances/faddeeva.jl")
include("resonances/reader.jl")
include("resonances/breit_wigner.jl")
include("resonances/reich_moore.jl")

# Processing modules (RECONR pipeline)
include("processing/adaptive_grid.jl")
include("processing/reconr.jl")
include("processing/pendf_writer.jl")

# Processing modules (BROADR pipeline -- Proposal B: sigma1 kernel + pipeline)
include("processing/sigma1.jl")
include("processing/broadr.jl")

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

end # module NJOY
