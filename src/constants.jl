# Physics constants matching NJOY2016 phys.f90
# Source: CODATA 2014 as given in ENDF-102 Appendix H (Feb 2018 edition)
#
# These are the *exact* values used by NJOY2016 v2016.78.
# Variable names use NJOY Fortran conventions for traceability.

"""
    PhysicsConstants

Collection of fundamental physics constants used throughout NJOY.
All values match NJOY2016's `phys.f90` (CODATA 2014 / ENDF-102 Appendix H).
"""
module PhysicsConstants

# Mathematical constants
const pi      = 3.141592653589793238    # pi
const euler   = 0.57721566490153286     # Euler-Mascheroni constant

# Fundamental constants (CODATA 2014 via ENDF-102 Appendix H)
const bk      = 8.617333262e-5          # Boltzmann constant [eV/K]
const ev      = 1.602176634e-12         # electron-volt [erg/eV]
const clight  = 2.99792458e10           # speed of light [cm/s]
const amu     = 931.49410242e6 * ev /
                (clight * clight)        # atomic mass unit [g]
const hbar    = 6.582119569e-16 * ev    # reduced Planck constant [erg*s]

# Inverse fine-structure constant: alpha^{-1} = hbar*c / (e^2)
# NJOY formula: finstri = 1e16 * hbar / (ev^2 * clight)
const finstri = 1.0e16 * hbar /
                (ev * ev * clight)       # inverse fine structure constant

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
# epair = m_e * amu_in_grams * c^2 / eV_in_ergs
const epair   = amasse * amu * clight * clight / ev

end # module PhysicsConstants

# Convenience alias
const CODATA2014 = PhysicsConstants
