"""Declared input-difference diagnostics, never an atmosphere accuracy estimate."""
import numpy as np


def diagnostic_depth_masks(native_logtau):
    """Select rows on the source grid, never a rounded interpolated coordinate."""
    tau = np.asarray(native_logtau)
    if tau.ndim != 1 or len(tau) < 2 or not np.isfinite(tau).all() or np.any(np.diff(tau) <= 0):
        raise ValueError('Require ordered finite native depth coordinates')
    subset = (tau >= -3) & (tau <= 0)
    if not np.any(subset):
        raise ValueError('Declared depth subset is empty')
    return [('all_depths', np.ones(len(tau), dtype=bool)), ('logtau_-3_0', subset)]


def field_difference(candidate, baseline, mask):
    """Max absolute T and fractional Pe/rho/ne differences on explicit native rows.

    ne is proportional to Pe/T; the shared Boltzmann constant cancels. Baseline
    denominators are intentional. Neither input is asserted to be physical truth.
    """
    fields = ('temperature_K', 'electron_pressure_Pa', 'rho_kg_m3')
    n = len(baseline['temperature_K'])
    mask = np.asarray(mask)
    if mask.dtype != bool or mask.shape != (n,) or not np.any(mask):
        raise ValueError('Require explicit nonempty boolean depth mask')
    for structure in (candidate, baseline):
        for name in fields:
            values = np.asarray(structure[name])
            if values.shape != (n,) or not np.isfinite(values).all() or np.any(values <= 0):
                raise ValueError('Invalid input field for sensitivity diagnostic')
    result = dict(max_abs_temperature_difference_K=float(np.max(abs(candidate['temperature_K'][mask]-baseline['temperature_K'][mask]))))
    for name, label in (('electron_pressure_Pa', 'Pe'), ('rho_kg_m3', 'rho')):
        result[f'max_abs_fractional_{label}_difference'] = float(np.max(abs(candidate[name][mask]/baseline[name][mask]-1)))
    ne_ratio = (candidate['electron_pressure_Pa']/candidate['temperature_K']) / (baseline['electron_pressure_Pa']/baseline['temperature_K'])
    result['max_abs_fractional_ne_difference'] = float(np.max(abs(ne_ratio[mask]-1)))
    return result
