"""
tests/test_procyon_ew_integrity_rya273.py
=========================================
RYA-273 — generalize the RYA-458 EW-integrity QA pass beyond solar.

Pins: the measured-EW pool resolves per star (solar -> the promoted canonical pool;
Procyon -> its own lines_fit output), the EW-verify gate is no longer solar-only, and
the layer flags an F-star's saturated profiles via COG_FLAG while never mutating an EW.
The solar charter behaviour (RYA-458) is unaffected by the generalization.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import abundances_derive as AD       # noqa: E402
import pipeline.ew_integrity as ei                  # noqa: E402


def test_measured_ew_path_solar_is_canonical():
    assert AD._measured_ew_path('solar') == AD.PATHS['solar_ew_canonical']
    assert AD._measured_ew_path('sun')   == AD.PATHS['solar_ew_canonical']


def test_measured_ew_path_procyon_prefers_star_pool(tmp_path):
    """Procyon resolves to its own EW staging output (procyon_ew.csv) when present —
    the file where the per-line chi2/profile/blend_flag live — not the solar pool."""
    p = AD.PATHS['procyon_ew']
    existed = p.exists()
    if not existed:
        p.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({'element': ['Fe'], 'ion': ['I'],
                      'wavelength_air_A': [5000.0], 'ew_mA': [50.0]}).to_csv(p, index=False)
    try:
        assert AD._measured_ew_path('procyon') == p
        assert p.name == 'procyon_ew.csv'
    finally:
        if not existed:
            p.unlink()


def test_unknown_star_falls_back_to_solar_canonical():
    # No silent crash for an unmapped star — documented fallback to the solar pool.
    assert AD._measured_ew_path('betelgeuse') == AD.PATHS['solar_ew_canonical']


def test_ew_integrity_flags_fstar_saturation_without_mutation():
    """F-star strong lines (REW above the linear COG knee, -4.90) -> COG_FLAG, and the
    measured EWs are returned byte-identical (flag-only; RYA-451/458 cardinal rule)."""
    df = pd.DataFrame({
        'element': ['Fe', 'Fe', 'Fe', 'Fe'],
        'ion':     ['I', 'I', 'I', 'I'],
        'wavelength_air_A': [5000.0, 5500.0, 6000.0, 6200.0],
        'ew_mA':   [120.0, 150.0, 200.0, 8.0],      # first three saturated, last weak
        'a_lte':   [7.5, 7.5, 7.5, 7.5],
    })
    before = df['ew_mA'].tolist()
    out = ei.flag_ew_integrity(df)
    assert out['ew_integrity'].astype(str).str.contains('COG_FLAG').any()
    # the weak line (8 mA) is below the knee -> not COG-flagged
    weak = out[out['wavelength_air_A'] == 6200.0].iloc[0]
    assert 'COG_FLAG' not in str(weak['ew_integrity'])
    # cardinal guard: no EW mutated
    assert out['ew_mA'].tolist() == before
    ei.assert_no_ew_mutation(df, out)


def test_ew_verify_gate_is_no_longer_solar_only():
    """The run() EW-verify call must not be solar-gated — the source line is
    'if ew_verify:' (RYA-273), so Procyon/others run through the layer too."""
    import inspect
    src = inspect.getsource(AD.run)
    assert 'if ew_verify and' not in src     # the old solar-only gate is gone
    assert 'if ew_verify:' in src
