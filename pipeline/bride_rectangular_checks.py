"""Diagnostic rectangular atmosphere interpolation; no solver input admission."""
import copy

import numpy as np


def rectangular_weights(points, target):
    """Return bilinear weights in input order; refuse incomplete cells/extrapolation."""
    p = np.asarray(points, dtype=float)
    t = np.asarray(target, dtype=float)
    if p.shape != (4, 2) or t.shape != (2,) or not np.isfinite(p).all() or not np.isfinite(t).all():
        raise ValueError('Require four finite 2D corners and a finite 2D target')
    x, y = np.unique(p[:, 0]), np.unique(p[:, 1])
    if len(x) != 2 or len(y) != 2 or len(set(map(tuple, p))) != 4:
        raise ValueError('Require all four corners of a nondegenerate rectangle')
    if np.any(t < [x[0], y[0]]) or np.any(t > [x[1], y[1]]):
        raise ValueError('Target outside rectangular cell; no extrapolation')
    u, v = (t - [x[0], y[0]]) / [x[1] - x[0], y[1] - y[0]]
    return np.array([(u if a == x[1] else 1-u) * (v if b == y[1] else 1-v) for a, b in p])


def rectangular_structure(nodes, target, *, thermodynamic_scale):
    """Bilinear common-tau diagnostic with declared linear/log field conventions.

    T remains linear. Positive Pe/rho and geometric layer thickness may use
    log interpolation, matching the existing experimental simplex conventions.
    No EOS/equilibrium solve, fitted offset or column metadata is supplied.
    """
    weights = rectangular_weights([p for p, s in nodes], target)
    reference = nodes[0][1]
    names = ('height_m', 'temperature_K', 'electron_pressure_Pa', 'rho_kg_m3', 'logtau_ross', 'logtau_500')
    for p, s in nodes:
        if s['ndep'] != reference['ndep'] or not np.array_equal(s['logtau_ross'], reference['logtau_ross']):
            raise ValueError('Atmospheres do not share native tauRoss rows')
        if s['abundances_dex'] != reference['abundances_dex'] or s['hydrogen_mass_fraction'] != reference['hydrogen_mass_fraction']:
            raise ValueError('Unmodeled composition change')
    if thermodynamic_scale not in ('linear', 'log_positive', 'log_positive_layers'):
        raise ValueError('Unknown thermodynamic interpolation convention')
    result = copy.deepcopy(reference)
    result.pop('column_mass_kg_m2', None)
    for name in names:
        values = np.stack([s[name] for p, s in nodes])
        if not np.isfinite(values).all():
            raise ValueError('Nonfinite source atmosphere')
        logarithmic = thermodynamic_scale.startswith('log_positive') and name in ('electron_pressure_Pa', 'rho_kg_m3')
        if logarithmic and np.any(values <= 0):
            raise ValueError('Nonpositive pressure/density')
        result[name] = np.exp(weights @ np.log(values)) if logarithmic else weights @ values
    if thermodynamic_scale == 'log_positive_layers':
        thickness = np.stack([-np.diff(s['height_m']) for p, s in nodes])
        if np.any(thickness <= 0):
            raise ValueError('Nonpositive layer thickness')
        dz = np.exp(weights @ np.log(thickness))
        anchor = int(np.argmin(abs(reference['logtau_ross'])))
        height = result['height_m']
        for i in range(anchor, len(height)-1):
            height[i+1] = height[i]-dz[i]
        for i in range(anchor-1, -1, -1):
            height[i] = height[i+1]+dz[i]
    if any(not np.isfinite(result[k]).all() for k in names) or any(np.any(result[k] <= 0) for k in ('temperature_K', 'electron_pressure_Pa', 'rho_kg_m3')) or np.any(np.diff(result['height_m']) >= 0):
        raise ValueError('Invalid candidate atmosphere')
    result['logg_cgs'] = float(target[1])
    result['model_id'] = f'DIAGNOSTIC_RECTANGULAR_T{target[0]}_g{target[1]}'
    return result, dict(method='bilinear diagnostic, common native tauRoss; no EOS re-solve',
                        weights=weights.tolist(), source_indices=list(range(4)), target=list(target))
