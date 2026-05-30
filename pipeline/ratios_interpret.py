"""
ratios_interpret.py
===================
Interpret derived abundances in the context of exoplanet science.

Outputs:
  - Solar-normalized abundance ratios [X/Fe] for 55 Cancri A
  - Comparison with solar abundances (Asplund 2021)
  - Comparison with literature values for 55 Cnc
  - Rock-forming element ratios: Mg/Si, C/O (informs planet bulk composition)
  - Summary table and publication-quality figures

Linear issue: RYA-TBD
"""

# TODO: Implement interpretation and plotting
# Key science outputs:
#   - [Fe/H], [Mg/Fe], [Si/Fe], [Ca/Fe], [O/Fe], [Ni/Fe]
#   - Mg/Si ratio (mantle-to-core mass fraction proxy)
#   - C/O ratio (if C lines available; affects mineralogy)
#   - Comparison with Valenti & Fischer (2005) for validation

from config.constants import SOLAR_ASPLUND2021, STAR_55CNC, PATHS


def run(star_id: str) -> None:
    """Called by run_pipeline.py. TODO: implement abundance interpretation."""
    raise NotImplementedError("ratios_interpret.run() is not yet implemented.")
