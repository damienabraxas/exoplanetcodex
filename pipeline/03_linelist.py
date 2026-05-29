"""
03_linelist.py
==============
Load and validate the master line list for EW analysis.

Sources:
  - VALD3 (Ryabchikova et al. 2015): atomic line data, log gf values
  - NIST ASD: transition probability grades (A+/A/B/C/D/E)

The merged master CSV lives at data/linelists/master_linelist.csv.
See data/linelists/README.md for the column schema.

Linear issue: RYA-TBD
"""

# TODO: Implement line list loading and validation
# The core loading logic is already in data/linelists/loader.py
# This script wraps it for the pipeline context:
#   - Load master list filtered to PIPELINE['min_nist_grade']
#   - Cross-match against the normalized spectrum wavelength range
#   - Flag potential blends
#   - Export filtered list for EW measurement

from config.constants import PIPELINE, PATHS
from data.linelists.loader import load_linelist, summarize_linelist
