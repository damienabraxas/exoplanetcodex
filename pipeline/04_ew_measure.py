"""
04_ew_measure.py
================
Measure equivalent widths (EW) of absorption lines.

Equivalent width definition:
  EW = ∫ (1 - F_norm(λ)) dλ   [measured in milliAngstroms]

For each line in the filtered line list:
  1. Extract a local window around the line center (±1.5 Å)
  2. Fit a Gaussian or Voigt profile to the absorption
  3. Integrate the fit to get EW
  4. Apply quality cuts: PIPELINE['ew_min_mA'] < EW < PIPELINE['ew_max_mA']

Linear issue: RYA-TBD
"""

# TODO: Implement EW measurement
# Key decisions:
#   - Gaussian vs Voigt profile (Voigt better for strong lines)
#   - Local continuum re-normalization around each line
#   - Blend rejection strategy
#   - Uncertainty estimation (fitting error + continuum uncertainty)

from config.constants import PIPELINE, PATHS
