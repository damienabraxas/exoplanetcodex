"""
spectra_normalize.py
====================
Continuum normalization of the HARPS spectrum.

Goal: divide out the stellar continuum so that absorption line depths
are measured relative to a local flux of ~1.0.

Approaches:
  - Iterative polynomial fitting (mask absorption lines, fit continuum)
  - Spline fitting with sigma-clipping
  - Order-by-order normalization using the S2D format

Linear issue: RYA-TBD
"""

# TODO: Implement continuum normalization
# Inputs:  wavelength, flux arrays from spectra_acquire.py
# Outputs: wavelength, flux_normalized (flux / continuum)
# Key decisions:
#   - Polynomial degree per echelle order
#   - Sigma-clipping iterations and threshold
#   - How to handle emission features (chromospheric Ca H&K)

from config.constants import PIPELINE, PATHS


def run(star_id: str) -> None:
    """Called by run_pipeline.py. TODO: implement continuum normalization."""
    raise NotImplementedError("spectra_normalize.run() is not yet implemented.")
