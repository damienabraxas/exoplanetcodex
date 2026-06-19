"""
Astrophysical priors for stellar parameters
"""

import numpy as np

def log_prior(theta):
    """
    θ = [Teff, logg, vturb, AFe, AC, AO, ...]
    """

    Teff, logg, vturb = theta[:3]

    # hard physical bounds
    if not (4000 < Teff < 7000):
        return -np.inf
    if not (3.0 < logg < 5.0):
        return -np.inf
    if not (0.5 < vturb < 3.0):
        return -np.inf

    # weak Gaussian priors (optional refinement later)
    lp = 0.0
    lp += -0.5 * ((Teff - 5772) / 200)**2  # RYA-298: canonical solar Teff (IAU/GBS)

    return lp
