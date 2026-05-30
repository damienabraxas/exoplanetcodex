"""
abundances_derive.py
====================
Derive stellar elemental abundances from measured EWs.

Method: 1D LTE stellar atmosphere modeling
  - Model: ATLAS9 (Kurucz) plane-parallel, LTE
  - Radiative transfer: MOOG (Sneden 1973) or equivalent
  - Iterative convergence on Teff, log g, vturb using:
      * Excitation equilibrium (Fe I lines vs. excitation potential)
      * Ionization equilibrium (Fe I abundance = Fe II abundance)
      * EW trend removal (abundance vs. reduced EW)

Linear issue: RYA-TBD
"""

# TODO: Implement abundance analysis
# Inputs:  EW table from lines_fit.py, stellar params from params_stellar.py
# Outputs: A(X) for each element, [X/H], [X/Fe]
# Key decisions:
#   - MOOG interface vs. pyMOOGi vs. iSpec
#   - Model atmosphere interpolation scheme
#   - Convergence criteria for stellar parameter self-consistency

from config.constants import STAR_55CNC, SOLAR_ASPLUND2021, MODEL, PATHS


def run(star_id: str) -> None:
    """Called by run_pipeline.py. TODO: implement abundance derivation."""
    raise NotImplementedError("abundances_derive.run() is not yet implemented.")
