"""Screen numerical extrapolation using freshly evaluated physical rate balance."""
import numpy as np


def select_fresh_population(plain,proposal,evaluate):
    """Leave the evaluator at the chosen state, including after rejection.

    `evaluate(populations)` must solve radiation at those fixed populations and
    return a finite nonnegative `residual` from the physical rate equations.
    Its state may be mutable. Re-evaluation restores the direct-map state after
    a rejected proposal; no radiative snapshot from another state is reused.
    """
    direct=evaluate(plain)
    if not np.isfinite(direct['residual'])or direct['residual']<0:raise ValueError('Invalid direct physical residual')
    info=dict(plain_fresh_residual=direct['residual'],proposal_fresh_residual=None,accepted=False,restored_plain=False)
    if np.array_equal(plain,proposal):return plain,direct,info
    trial=evaluate(proposal)
    if not np.isfinite(trial['residual'])or trial['residual']<0:raise ValueError('Invalid proposal physical residual')
    info['proposal_fresh_residual']=trial['residual']
    if trial['residual']<direct['residual']:
        info['accepted']=True;return proposal,trial,info
    restored=evaluate(plain)
    if not np.isfinite(restored['residual'])or restored['residual']<0:raise ValueError('Invalid restored physical residual')
    info['restored_plain']=True
    return plain,restored,info
