"""
Nested sampling inference engine
"""

import numpy as np
import dynesty

from .forward_model import run_moog
from .likelihood import log_likelihood
from .priors import log_prior


class StellarInference:
    def __init__(self, ew_obs, ew_err, line_list):
        self.ew_obs = ew_obs
        self.ew_err = ew_err
        self.line_list = line_list

    def prior_transform(self, u):
        """
        Transform unit cube → physical parameters
        """

        Teff  = 4000 + u[0] * 3000
        logg  = 3.0   + u[1] * 2.0
        vturb = 0.5   + u[2] * 2.5

        # abundances (simple starting point)
        AFe = -1.0 + u[3] * 2.0
        AC  = -1.0 + u[4] * 2.0
        AO  = -1.0 + u[5] * 2.0

        return np.array([Teff, logg, vturb, AFe, AC, AO])

    def log_likelihood(self, theta):
        ew_model = run_moog(theta, self.line_list)
        return log_likelihood(self.ew_obs, ew_model, self.ew_err)

    def run(self, nlive=500):
        """
        Run nested sampling
        """

        sampler = dynesty.NestedSampler(
            self.log_likelihood,
            self.prior_transform,
            ndim=6,
            nlive=nlive
        )

        sampler.run_nested()

        return sampler.results
