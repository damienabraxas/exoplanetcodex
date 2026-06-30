"""
tests/test_solar_n_bank_rya369.py
=================================
RYA-369 follow-on — solar N I measurement. The verdict is that solar N I red is
continuum-limited (a ~0.3-dex lever on 2-4% deep lines), so the bank is provisional
and the prior phase_c 8.202 is the raw/un-localized upper bracket, not a bug. These
tests pin the FINDING's premises (CI-safe — no synthesis): the KPNO pseudo-continuum
is blanketed below 1.0, the lines are weak, and the solar N I NLTE is small.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KPNO = ROOT / 'data' / 'solar_reference' / 'kpno_flux_atlas'
N_LINES = [(7468.31, 'kpno_NI_7442_7468'), (8216.34, 'kpno_NI_8216_8223'),
           (8683.40, 'kpno_NI_8680_8718')]


def test_kpno_red_pseudocontinuum_is_blanketed_below_unity():
    # No true continuum point in-window -> the atlas-1.0 (raw) treatment over-reads A(N);
    # this is why the measurement is continuum-limited.
    for wl, f in N_LINES:
        d = pd.read_csv(KPNO / f'{f}.csv')
        m = (d['wavelength_air_A'] >= wl - 1.0) & (d['wavelength_air_A'] <= wl + 1.0)
        assert float(d['residual_flux'][m].max()) < 0.999      # never reaches 1.0


def test_n_lines_are_weak():
    # 2-4% deep -> a 0.5-1% continuum offset is a large fractional depth error.
    for wl, f in N_LINES:
        d = pd.read_csv(KPNO / f'{f}.csv')
        core = (d['wavelength_air_A'] > wl - 0.15) & (d['wavelength_air_A'] < wl + 0.15)
        depth = 1.0 - float(d['residual_flux'][core].min())
        assert 0.0 < depth < 0.10                              # shallow line


def test_solar_n_nlte_is_small_near_lte():
    import json
    j = json.loads((ROOT / 'data' / 'results' / 'n_grid_loadtest_rya369.json').read_text())
    deltas = [abs(float(v)) for v in j['resolve']['per_line'].values()]
    assert max(deltas) < 0.05                                  # solar N I red is near-LTE


def test_provisional_bank_is_continuum_limited():
    # The run records the lever as a systematic, not a clean millidex bank.
    import json
    p = ROOT / 'data' / 'results' / 'solar_n_bank_rya369.json'
    if not p.exists():
        return                                                 # artifact optional in CI
    o = json.loads(p.read_text())
    assert o['continuum_limited'] is True
    assert o['sigma_continuum'] > 0.1                          # lever folded into sigma
    assert o['a_n_local'] < o['a_n_raw']                       # local re-norm shallower -> lower A
