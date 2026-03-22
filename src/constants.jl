# Physics constants matching NJOY2016 phys.f90 — CGS units (erg, g, cm)
# Source: CODATA 2014 as given in ENDF-102 Appendix H (Feb 2018 edition)
#
# These are the *exact* values used by NJOY2016 v2016.78 (phys.f90).
# Variable names use NJOY Fortran conventions for traceability.
#
# IMPORTANT: CGS units throughout, matching the Fortran exactly.
# All formulas in the codebase (cwaven, channel_radius, penetrability,
# phase shift, etc.) were ported from Fortran and assume CGS.
# An earlier session incorrectly changed these to SI, breaking cwaven
# by 100x (0.002197 → 0.2197) and all resonance cross sections.

"""
    PhysicsConstants

Collection of fundamental physics constants used throughout NJOY.
All values match NJOY2016's `phys.f90` (CODATA 2014 / ENDF-102 Appendix H).
Units: CGS (ergs, grams, centimetres, seconds) — matching Fortran exactly.
"""
module PhysicsConstants

# Mathematical constants
const pi      = 3.141592653589793238    # pi
const euler   = 0.57721566490153286     # Euler-Mascheroni constant

# Fundamental constants (CODATA 2014 via ENDF-102 Appendix H) — CGS
# These values are copied verbatim from njoy-reference/src/phys.f90
const bk      = 8.617333262e-5          # Boltzmann constant [eV/K]
const ev      = 1.602176634e-12         # electron-volt [erg/eV]
const clight  = 2.99792458e10           # speed of light [cm/s]

# Derived constants matching phys.f90 exactly:
#   amu  = 931.49410242e6 * ev / clight^2  [g]
#   hbar = 6.582119569e-16 * ev            [erg·s]
const amu     = 931.49410242e6 * ev / (clight * clight)  # atomic mass unit [g]
const hbar    = 6.582119569e-16 * ev                      # reduced Planck [erg·s]

# Inverse fine-structure constant (derived in phys.f90 as 1e16*hbar/(ev*ev*clight))
const finstri = 1.0e16 * hbar / (ev * ev * clight)       # α⁻¹ ≈ 137.036

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
