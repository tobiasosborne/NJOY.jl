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

end # module NJOY
