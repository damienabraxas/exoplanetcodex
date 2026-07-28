"""
tests/test_fe_leg_placeholder_guard_rya417.py
=============================================
RYA-417 — extend the RYA-413 placeholder-zero protection to the SEPARATE primary-Fe
Bergemann NLTE leg (_load_mpia_fe_grid / apply_fe_nlte_corrections), the path the
registry guard did not cover (where the 31 Fe placeholder-zeros + Si 6253 survived).
The guard must be symmetric: a line identically 0 across all nodes is refused on the
Fe leg exactly as on the registry leg — and a GENUINE near-zero (real, varying) is kept.
"""
import pandas as pd
import pytest

from pipeline import nlte_corrections as N
from pipeline.nlte_corrections import detect_placeholder_zero_lines

_COLS = ['element', 'ion', 'wave_A', 'teff_K', 'logg', 'feh', 'delta_nlte']
_NODES = [(5400, 4.3, 0.0), (5800, 4.4, 0.0), (6200, 4.5, 0.2), (6600, 4.6, 0.35)]
# non-degenerate cube (8 nodes) for cases that must build a LinearND interpolator
_CUBE = [(t, g, f) for t in (5400, 6200) for g in (4.0, 4.6) for f in (-0.2, 0.35)]


def _fe_rows(ion, wave, vals):
    return [['Fe', ion, wave, t, g, f, v] for (t, g, f), v in zip(_NODES, vals)]


def _fe_cube(ion, wave, val):
    return [['Fe', ion, wave, t, g, f, val] for (t, g, f) in _CUBE]


def test_detect_is_ion_aware_on_fe_grid():
    # same wave in Fe I (real) and Fe II (placeholder) — only the Fe II one is flagged
    df = pd.DataFrame(_fe_rows('I', 5225.5, [0.01, 0.02, 0.03, 0.02]) +
                      _fe_rows('II', 5225.5, [0.0, 0.0, 0.0, 0.0]), columns=_COLS)
    assert detect_placeholder_zero_lines(df[df.ion == 'II']) == [5225.5]
    assert detect_placeholder_zero_lines(df[df.ion == 'I']) == []


def test_fe_leg_loader_refuses_placeholder_grid(tmp_path, monkeypatch):
    bad = pd.DataFrame(_fe_rows('I', 5269.5, [0.01, 0.02, 0.03, 0.02]) +
                       _fe_rows('I', 6430.846, [0.0, 0.0, 0.0, 0.0]), columns=_COLS)
    p = tmp_path / 'Fe_bad.csv'
    p.write_text(bad.to_csv(index=False))
    monkeypatch.setattr(N, '_MPIA_FE_GRID', p)
    N._mpia_cache.clear()
    with pytest.raises(ValueError, match='PLACEHOLDER_ZERO'):
        N._load_mpia_fe_grid()
    N._mpia_cache.clear()


def test_fe_leg_loader_keeps_genuine_near_zero(tmp_path, monkeypatch):
    # a small-but-real correction with at least one nonzero node must NOT be refused
    ok = pd.DataFrame(_fe_cube('I', 5269.5, 0.001)[:7] + [['Fe', 'I', 5269.5, 6200, 4.6, 0.35, 0.003]] +
                      _fe_cube('II', 6432.0, 0.002)[:7] + [['Fe', 'II', 6432.0, 6200, 4.6, 0.35, 0.004]],
                      columns=_COLS)
    p = tmp_path / 'Fe_ok.csv'
    p.write_text(ok.to_csv(index=False))
    monkeypatch.setattr(N, '_MPIA_FE_GRID', p)
    N._mpia_cache.clear()
    c = N._load_mpia_fe_grid()           # must not raise
    assert 1 in c['waves'] and 2 in c['waves']
    N._mpia_cache.clear()


def test_committed_fe_grid_loads_clean():
    # RYA-417 Part C dropped the 31 verified placeholders — the live grid loads clean.
    N._mpia_cache.clear()
    c = N._load_mpia_fe_grid()
    assert len(c['waves'].get(1, [])) > 100 and len(c['waves'].get(2, [])) > 30
    # none of the dropped placeholder waves remain
    for w in (5225.533, 5855.077, 6430.846, 6369.375):
        assert round(w, 0) not in [round(float(x), 0) for x in c['waves'].get(1, [])
                                   if abs(float(x) - w) < 0.05]
    N._mpia_cache.clear()


def test_curate_si_leg_refuses_placeholder(tmp_path, monkeypatch):
    import pipeline.curate_nonfe_pools as CU
    bad = pd.DataFrame([['Si', None, 6253.6, t, g, f, 0.0] for (t, g, f) in _NODES],
                       columns=_COLS).drop(columns=['ion'])
    (tmp_path / 'Si_bad.csv').write_text(bad.to_csv(index=False))
    monkeypatch.setattr(CU, '_NLTE_GRID_DIR', tmp_path)
    monkeypatch.setitem(CU._NLTE_GRID_FILES, 'Si', 'Si_bad.csv')
    with pytest.raises(ValueError, match='PLACEHOLDER_ZERO'):
        CU.solar_nlte_delta('Si')
