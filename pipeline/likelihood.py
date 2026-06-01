import numpy as np

def log_likelihood(ew_obs, ew_model, ew_err):
    return -0.5 * np.sum(
        ((ew_obs - ew_model) / ew_err) ** 2
    )
