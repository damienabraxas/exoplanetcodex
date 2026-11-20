"""Hydrogen density an LTE equation of state infers from (T, Pe) for a given abundance mix.

RYA-1282, 2026-11-20. Diagnostic for the Gerber+2023 MULTI 2.3 reference grid. MULTI reads
an atmosphere as T and n_e on a tau5000 scale and gets n(H) from its own EOS with the
abundances in its ABUND file. If ABUND is solar for every model (pymulti, the MULTI wrapper
our MULTI source comes from, always writes solar Asplund+2009), the inferred n(H) differs
from the model's by  r = e(model mix) / e(solar mix),  e = free electrons per H nucleus at
the model's (T, n_e). On a tau scale the opacities and the oxygen count all scale with
n(H), so their ratios cancel; what remains is the H-collision rate, scaled by r.

Measured (docs/audit/rya1282/lightweaver_differences.md): scaling only the Bride's CH rates
by r, with nothing fitted, reproduces the grid's 777-triplet departures over [Fe/H] +1..-5
for 4500-5750 K dwarfs (median error 0.023 -> 0.005 dex over 120 nodes).

Saha with TS's Irwin partition functions and ionisation energies; H- included. Reproduces
MARCS's own n(H) from (T, Pe) to 1-1.5% (tests/test_bride_reference_eos_rya1282.py).
"""
from __future__ import annotations

import numpy as np

_K = 1.380649e-16
_H = 6.62607015e-27
_ME = 9.1093837e-28
_EV = 1.602176634e-12


def electrons_per_H(T, ne_cm3, abundances_dex, irwin_path=None):
    """Free electrons per H nucleus at (T, n_e) for the mix A(X) (Z = 1..len)."""
    from pipeline.bride_line_haze import _stage_fractions
    if irwin_path is None:
        from pipeline.bride_ts_background import ts_data_dir, TsMetalBackground
        irwin_path = ts_data_dir() / TsMetalBackground.IRWIN
    T = np.asarray(T, dtype=float); ne = np.asarray(ne_cm3, dtype=float)
    e = np.zeros_like(T)
    for Z in range(1, len(abundances_dex) + 1):
        A = abundances_dex[Z - 1]
        if not np.isfinite(A) or A < -5:
            continue
        try:
            f, _ = _stage_fractions(irwin_path, Z, T, ne)
        except ValueError:                      # no Irwin data for this element
            continue
        e += 10. ** (A - 12.) * sum(s * f[s] for s in range(len(f)))
    phi = (2 * np.pi * _ME * _K * T / _H ** 2) ** 1.5
    fH, UH = _stage_fractions(irwin_path, 1, T, ne)
    h_minus = fH[0] * ne / (2. * phi / UH[0]) * np.exp(0.754 * _EV / (_K * T))
    return e - h_minus


def nh_from_pe(T, ne_cm3, abundances_dex, irwin_path=None):
    """n(H nuclei) (cm^-3) implied by (T, n_e) and the mix."""
    return np.asarray(ne_cm3, dtype=float) / electrons_per_H(T, ne_cm3, abundances_dex, irwin_path)


def solar_mix_nh_ratio(T, ne_cm3, model_abundances_dex, solar_abundances_dex, irwin_path=None):
    """n(H) a solar-mix EOS would infer / n(H) of the model's own mix, at the same (T, n_e)."""
    return (electrons_per_H(T, ne_cm3, model_abundances_dex, irwin_path)
            / electrons_per_H(T, ne_cm3, solar_abundances_dex, irwin_path))
