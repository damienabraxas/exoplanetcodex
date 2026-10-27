"""Physical flow accounting, separate from the accelerated SE operator."""
import numpy as np


def process_flow(matrix, populations):
    """Off-diagonal destination/source rates -> signed and gross level flows.

    Rates are s^-1; populations m^-3; outputs m^-3 s^-1. Diagonals are ignored.
    """
    if matrix.shape != (len(populations), len(populations), populations.shape[1]):
        raise ValueError('Rate/population shapes differ')
    rates = matrix.copy()
    for i in range(len(populations)):
        rates[i, i] = 0.
    if not np.isfinite(rates).all() or np.any(rates < 0) or not np.isfinite(populations).all() or np.any(populations <= 0):
        raise ValueError('Positive populations and nonnegative finite physical rates required')
    incoming = np.einsum('ijd,jd->id', rates, populations)
    outgoing = rates.sum(axis=0) * populations
    return incoming-outgoing, incoming+outgoing


def collision_classes(atom, atmosphere, eqpops):
    """Evaluate each already-initialized native collision class separately."""
    groups = {}
    shape = (len(atom.levels), len(atom.levels), atmosphere.Nspace)
    for collision in atom.collisions:
        key = type(collision).__name__
        matrix = groups.setdefault(key, np.zeros(shape))
        collision.compute_rates(atmosphere, eqpops, matrix)
    return groups
