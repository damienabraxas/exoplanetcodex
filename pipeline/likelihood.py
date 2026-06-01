import numpy as np

def log_likelihood(ew_obs, ew_model, ew_err):
    """
    Gaussian likelihood in EW space
    """

    residual = ew_obs - ew_model
    sigma2 = ew_err ** 2 + 1e-6  # prevent division issues

    return -0.5 * np.sum(residual**2 / sigma2 + np.log(sigma2))
