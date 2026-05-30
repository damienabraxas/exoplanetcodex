"""
lines_load.py
=============
Load and validate the master line list for EW analysis.

Sources:
  - VALD3 (Ryabchikova et al. 2015): atomic line data, log gf values
  - NIST ASD: transition probability grades (A+/A/B/C/D/E)

The merged master CSV lives at data/linelists/linelist_master.csv.
See data/linelists/README.md for the column schema.

Linear issue: RYA-62 (complete)
"""

# TODO: Implement pipeline wrapper around data/linelists/loader.py
# The core loading logic is in data/linelists/loader.py
# This script wraps it for the pipeline context:
#   - Load master list filtered to PIPELINE['min_nist_grade']
#   - Cross-match against the normalized spectrum wavelength range
#   - Flag potential blends
#   - Export filtered list for lines_fit.py

from config.constants import PIPELINE, PATHS
from data.linelists.loader import load_linelist, summarize_linelist


def run(star_id: str) -> None:
    """Called by run_pipeline.py. TODO: implement pipeline line list wrapper."""
    raise NotImplementedError("lines_load.run() is not yet implemented.")
