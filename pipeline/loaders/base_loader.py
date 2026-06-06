"""
pipeline/loaders/base_loader.py
================================
Abstract base class and standard output container for all instrument loaders.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np


@dataclass
class SpectrumData:
    """
    Standard internal representation of a reduced spectrum.

    All loaders return this object. Downstream pipeline (normalize,
    lines_fit, abundances_derive) only ever sees SpectrumData — never
    raw FITS structures.

    Attributes
    ----------
    wave_A : np.ndarray
        Wavelength array in Angstroms, barycentric frame.
        Shape: (n_pixels,) for merged 1D; (n_orders, n_pixels) for echelle.
    flux : np.ndarray
        Flux array (same shape as wave_A). Telluric-corrected where applicable.
        NaN marks regions where correction failed (deep telluric bands, gaps).
    err : np.ndarray or None
        Flux uncertainty (same shape). None if not available from instrument.
    meta : dict
        Provenance and observation metadata. Required keys documented below.
    """
    wave_A : np.ndarray
    flux   : np.ndarray
    err    : np.ndarray | None
    meta   : dict = field(default_factory=dict)

    _REQUIRED_META = [
        'instrument',
        'object',
        'date_obs',
        'exptime_s',
        'snr_summary',
        'apero_version',
        'berv_kms',
        'wave_units_raw',
        'telluric_corrected',
        'wapiti_applied',
        'filepath',
    ]

    def __post_init__(self):
        assert self.wave_A.shape == self.flux.shape, \
            f"wave/flux shape mismatch: {self.wave_A.shape} vs {self.flux.shape}"
        if self.err is not None:
            assert self.err.shape == self.flux.shape, \
                f"err shape mismatch: {self.err.shape} vs {self.flux.shape}"
        missing = [k for k in self._REQUIRED_META if k not in self.meta]
        if missing:
            raise ValueError(f"SpectrumData missing meta keys: {missing}")

    @property
    def n_pixels(self) -> int:
        return self.flux.size

    @property
    def wave_range_A(self) -> tuple[float, float]:
        valid = self.wave_A[np.isfinite(self.wave_A)]
        return float(valid.min()), float(valid.max())

    @property
    def nan_fraction(self) -> float:
        return float(np.isnan(self.flux).sum() / self.flux.size)
