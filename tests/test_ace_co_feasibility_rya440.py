"""
tests/test_ace_co_feasibility_rya440.py
=======================================
RYA-440 — guard the ACE-FTS solar CO feasibility module. These are the FAST,
deterministic checks (line-list inventory, the air<->vac LOUD-FAIL boundary, the
ACE pointing geometry, constants wiring). The full Turbospectrum synthesis +
A(C) fit (~45 s, needs the iSpec install) is exercised by the CLI smoke test
`python -m pipeline.ace_co_feasibility --validate`, not in CI.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline import ace_co_feasibility as ace  # noqa: E402


def test_constants_single_source():
    # A(C)/A(O) and the C/O tolerance come from constants.py / SOLAR_VIS_GATES,
    # not round-defaults.
    assert ace.A_C_REF == 8.46
    assert ace.A_O_REF == 8.69
    assert ace.GEOM_TOL_DEX == 0.05
    p = ace._solar_params()
    assert p['teff_K'] == 5778.0 and p['logg'] == 4.438
    assert p['feh'] == 0.0 and p['vturb_kms'] == 1.0


def test_step0_linelists():
    s0 = ace.step0_linelists()
    assert s0['co12_present'] is True
    assert s0['co12_species_token'].startswith('0608.012016')   # 12C16O
    assert s0['co12_lines_in_overtone_band'] > 100              # CO band populated
    # 13C16O is absent on disk -> the fast-follow flag must be False/empty
    assert s0['co13_present'] is False
    assert 'band-scoped symlink' in s0['wiring_note']


def test_step5_geometry_collapses_to_mu1():
    s5 = ace.step5_geometry()
    lo, hi = s5['mu_range_sampled']
    assert hi == 1.0 and lo > 0.98          # disk-center sampling
    assert 'Hase' in s5['reference']
    assert 'disk-CENTER' in s5['finding'] or 'disk-center' in s5['finding'].lower()


def _fake_ace(tmp_path, monkeypatch, vac_head, air_head):
    """Write a tiny ACE-like CSV with a single absorption dip at a chosen vac/air
    bandhead, point the module at it."""
    w_vac = np.linspace(22925.0, 22965.0, 401)
    w_air = w_vac - (vac_head - air_head)        # constant air/vac offset
    flux = np.ones_like(w_vac)
    flux[np.argmin(np.abs(w_vac - vac_head))] = 0.3   # the dip = bandhead
    import pandas as pd
    df = pd.DataFrame({'wavelength_vac_A': w_vac, 'wavelength_air_A': w_air,
                       'intensity': flux})
    csv = tmp_path / 'ace.csv'
    df.to_csv(csv, index=False)
    monkeypatch.setattr(ace, 'ACE_CSV', csv)


def test_airvac_matched_vacuum(tmp_path, monkeypatch):
    # Synthesis (vacuum) bandhead at 22935.3; vac column aligned, air offset ~6.3 A.
    _fake_ace(tmp_path, monkeypatch, vac_head=22935.3, air_head=22929.0)
    out = ace.step1_load_ace(synth_bandhead_A=22935.3)
    assert out['matched_frame'] == 'wavelength_vac_A'
    assert abs(out['airvac_slip_kms']) > 50          # the ~83 km/s slip on the air column


def test_airvac_loud_fail_on_total_mismatch(tmp_path, monkeypatch):
    # Bandhead nowhere near the synthesis fiducial -> neither column aligns -> STOP.
    _fake_ace(tmp_path, monkeypatch, vac_head=22950.0, air_head=22944.0)
    with pytest.raises(ValueError, match='AIR/VAC LOUD-FAIL'):
        ace.step1_load_ace(synth_bandhead_A=22935.3)
