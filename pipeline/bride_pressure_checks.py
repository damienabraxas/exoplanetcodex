"""Plane-parallel diagnostic pressure/column checks; no equilibrium or cap admission."""
from pathlib import Path
import numpy as np


def read_pressure_table(path, structure):
    lines = Path(path).read_text().splitlines()
    if not any('No mass for plane-parallel models' in line for line in lines):
        raise ValueError('Pressure audit requires an explicit plane-parallel model')
    index = next(i for i, line in enumerate(lines) if line.split() == ['k', 'lgTauR', 'lgTau5', 'Depth', 'T', 'Pe', 'Pg', 'Prad', 'Pturb'])
    n = structure['ndep']
    a = np.asarray([[float(x) for x in line.split()] for line in lines[index+1:index+1+n]])
    if a.shape != (n, 9) or not np.isfinite(a).all() or not np.array_equal(a[:, 0], np.arange(1, n+1)):
        raise ValueError('Malformed native pressure rows')
    if not np.array_equal(a[:, 1], structure['logtau_ross']) or not np.array_equal(a[:, 4], structure['temperature_K']) or not np.array_equal(a[:, 5]*.1, structure['electron_pressure_Pa']):
        raise ValueError('Pressure/native field row identities differ')
    if np.any(a[:, 6] <= 0) or np.any(a[:, 7:] < 0):
        raise ValueError('Invalid native pressures')
    return dict(gas_pressure_Pa=a[:, 6]*.1, radiation_pressure_Pa=a[:, 7]*.1, turbulent_pressure_Pa=a[:, 8]*.1)


def pressure_column_metrics(pressure_Pa, column_kg_m2, gravity_m_s2):
    """Differences eliminate an unknown pressure intercept; not top-column recovery.

    Plane-parallel d(Pg+Prad+Pturb)/dm=g, from Gustafsson2008 Eqs1-3.
    Finite adjacent intervals and top-to-row differences are separate diagnostics.
    """
    p = np.asarray(pressure_Pa, dtype=float)
    m = np.asarray(column_kg_m2, dtype=float)
    if p.ndim != 1 or len(p) < 2 or m.shape != p.shape or not np.isfinite(p).all() or not np.isfinite(m).all() or np.any(p <= 0) or np.any(m <= 0) or np.any(np.diff(m) <= 0) or not np.isfinite(gravity_m_s2) or gravity_m_s2 <= 0:
        raise ValueError('Invalid pressure/column/gravity inputs')
    adjacent = np.diff(p)/(gravity_m_s2*np.diff(m))-1
    anchored = (p[1:]-p[0])/(gravity_m_s2*(m[1:]-m[0]))-1
    return dict(max_abs_adjacent_hydrostatic_difference=float(abs(adjacent).max()),
                max_abs_top_anchored_hydrostatic_difference=float(abs(anchored).max()))


def pressure_column_candidate(pressure_Pa, gravity_m_s2, top_column_kg_m2):
    """Construct a hydrostatic column with an external top anchor; not validation."""
    p = np.asarray(pressure_Pa, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or np.any(p <= 0) or np.any(np.diff(p) <= 0) or not np.isfinite(gravity_m_s2) or gravity_m_s2 <= 0 or not np.isfinite(top_column_kg_m2) or top_column_kg_m2 <= 0:
        raise ValueError('Invalid pressure-column construction inputs')
    return top_column_kg_m2+(p-p[0])/gravity_m_s2
