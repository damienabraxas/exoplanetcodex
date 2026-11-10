"""Read MARCS MRXF continuum tables using the official readopa.f contract.

Wavelengths retain their source Angstrom coordinate: the inspected reader does
not state its air/vacuum medium. This reader does not admit a spectral comparison.
"""
import gzip
from pathlib import Path

import numpy as np


def read_marcs_opacity(path):
    path = Path(path)
    raw = path.read_bytes()
    text = (gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw).decode('ascii')
    fields = text.split()
    if len(fields) < 4 or fields[0] != 'MRXF':
        raise ValueError('Expected MARCS MRXF opacity header')
    n, standard, nwave = int(fields[1]), float(fields[2]), int(fields[3])
    if not (1 <= n <= 400 and 2 <= nwave <= 10000 and standard > 0):
        raise ValueError('Invalid MARCS opacity dimensions')
    values = np.asarray([float(x.replace('D', 'E')) for x in fields[4:]])
    if len(values) != nwave + n * (8 + 2*nwave) + 92 or not np.isfinite(values).all():
        raise ValueError('Truncated/extra/nonfinite MARCS opacity data')
    wave = values[:nwave]
    if np.any(wave <= 0) or np.any(np.diff(wave) <= 0):
        raise ValueError('Non-increasing opacity wavelengths')
    depth = values[nwave:-92].reshape(n, 8+2*nwave)
    structure = depth[:, :8]
    relative = depth[:, 8:].reshape(n, nwave, 2)
    if (np.any(structure[:, :6] <= 0) or np.any(structure[:, 6] < 0)
            or np.any(structure[:, 7] <= 0) or np.any(np.diff(structure[:, 1]) <= 0)
            or np.any(relative < 0) or values[-92] != 12.):
        raise ValueError('Invalid MARCS opacity physical fields')
    # File coefficients are ratios to OPS, not absolute cm2/g. Multiplication
    # by OPS is required by readopa.f lines 125-131; 1 cm2/g = 0.1 m2/kg.
    absolute = relative * structure[:, 7, None, None] * .1
    return dict(ndep=n, nwave=nwave, standard_wavelength_source_A=standard,
        wavelength_source_A=wave, wavelength_medium='UNVERIFIED',
        radius_outer_normalized=structure[:, 0], tau_standard=structure[:, 1],
        temperature_K=structure[:, 2], electron_pressure_Pa=structure[:, 3]*.1,
        gas_pressure_Pa=structure[:, 4]*.1, rho_kg_m3=structure[:, 5]*1000,
        vturb_m_s=structure[:, 6]*1000, standard_opacity_m2_kg=structure[:, 7]*.1,
        absorption_m2_kg=absolute[:, :, 0].T, scattering_m2_kg=absolute[:, :, 1].T,
        abundances_dex=values[-92:])
