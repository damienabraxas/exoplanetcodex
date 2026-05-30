"""
params_stellar.py
=================
Determine self-consistent stellar parameters via spectroscopic equilibria.

Method: iterative LTE convergence on three criteria:
  1. Excitation equilibrium  — Fe I abundance must NOT correlate with
     excitation potential → constrains Teff
  2. Ionization equilibrium  — Fe I abundance must equal Fe II abundance
     → constrains log g
  3. Microturbulence         — Fe I abundance must NOT correlate with
     reduced EW (log EW/λ) → constrains vturb

Starting values: STAR_55CNC from config/constants.py (CHARA interferometry)
Output: refined (Teff, log g, [Fe/H], vturb) with formal uncertainties

Linear issue: RYA-TBD
"""

# TODO: Implement stellar parameter determination
# Inputs:  EW table from lines_fit.py, initial params from STAR_55CNC
# Outputs: converged (Teff, logg, feh, vturb), delta values, convergence plot
# Key decisions:
#   - Convergence criterion: slope < 0.005 dex/eV for excitation equilibrium
#   - Step size for grid search vs gradient descent
#   - How many Fe I/II lines are required for reliable convergence

from config.constants import STAR_55CNC, PIPELINE, PATHS


def run(star_id: str) -> dict:
    """Called by run_pipeline.py. TODO: implement stellar parameter determination."""
    raise NotImplementedError("params_stellar.run() is not yet implemented.")
