"""
Forward model layer:
θ → synthetic equivalent widths via MOOG
"""

import subprocess
import numpy as np
import pandas as pd
from pathlib import Path

from config.constants import PATHS


def run_moog(params, line_list):
    """
    Run MOOG with given stellar parameters.

    Parameters
    ----------
    params : dict
        Teff, logg, vturb, abundances
    line_list : DataFrame
        Selected spectral lines

    Returns
    -------
    ew_model : DataFrame
        Predicted equivalent widths (mÅ)
    """

    # 1. write MOOG input file
    input_file = _write_moog_input(params, line_list)

    # 2. call MOOG
    subprocess.run(["moog", input_file], check=True)

    # 3. parse output
    ew_model = _parse_moog_output()

    return ew_model


def _write_moog_input(params, line_list):
    """
    Convert θ → MOOG input format
    """
    # placeholder
    return Path("moog_input.in")


def _parse_moog_output():
    """
    Parse MOOG output EW table
    """
    # placeholder parser
    return pd.DataFrame()
