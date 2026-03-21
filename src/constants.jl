# Physics constants matching NJOY2016 phys.f90 — SI units (J, kg, m)
# Source: CODATA 2014 as given in ENDF-102 Appendix H (Feb 2018 edition)
#
# These are the *exact* values used by NJOY2016 v2016.78.
# Variable names use NJOY Fortran conventions for traceability.
#
# IMPORTANT: SI units throughout. All formulas ported from NJOY2016 Fortran
# assume SI. Previous versions of this file used CGS (ergs, grams, cm/s),
# causing cwaven_constant() to be 100x too small and pifac 10,000x too large.

"""
    PhysicsConstants

Collection of fundamental physics constants used throughout NJOY.
All values match NJOY2016's `phys.f90` (CODATA 2014 / ENDF-102 Appendix H).
Units: SI (Joules, kilograms, metres, seconds).
"""
module PhysicsConstants

# Mathematical constants
const pi      = 3.141592653589793238    # pi
const euler   = 0.57721566490153286     # Euler-Mascheroni constant

# Fundamental constants (CODATA 2014 via ENDF-102 Appendix H) — SI
const bk      = 8.617333262e-5          # Boltzmann constant [eV/K]
const ev      = 1.602176634e-19         # electron-volt [J/eV]
const clight  = 2.99792458e8            # speed of light [m/s]
const amu     = 1.66053906660e-27       # atomic mass unit [kg]
const hbar    = 1.054571817e-34         # reduced Planck constant [J·s]

# Inverse fine-structure constant α⁻¹ ≈ 137.036
const finstri = 137.035999084            # CODATA 2018 exact value

# ------------------------------------------------------------------
# Light particle masses in atomic mass units (amu)
# These are *particle* masses, not atomic masses.
# Source: ENDF-102 Appendix H
# ------------------------------------------------------------------
const amassn  = 1.00866491595           # neutron mass [amu]
const amassp  = 1.007276466621          # proton mass [amu]
const amassd  = 2.013553212745          # deuteron mass [amu]
const amasst  = 3.01550071621           # triton mass [amu]
const amassh  = 3.014932247175          # helion (He-3 nucleus) mass [amu]
const amassa  = 4.001506179127          # alpha particle mass [amu]
const amasse  = 5.48579909065e-4        # electron mass [amu]

# Mass ratios to neutron mass (derived)
const pnratio = amassp / amassn         # proton/neutron
const dnratio = amassd / amassn         # deuteron/neutron
const tnratio = amasst / amassn         # triton/neutron
const hnratio = amassh / amassn         # helion/neutron
const anratio = amassa / amassn         # alpha/neutron

# Pair production threshold energy [eV]
# epair = m_e [amu] * amu [kg] * c² [m²/s²] / ev [J/eV]
const epair   = amasse * amu * clight * clight / ev

end # module PhysicsConstants

# Convenience alias
const CODATA2014 = PhysicsConstants
