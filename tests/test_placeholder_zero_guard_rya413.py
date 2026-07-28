"""
tests/test_placeholder_zero_guard_rya413.py
===========================================
RYA-413 — the placeholder-zero guard. MPIA OFFERS some lines in its dropdown but RETURNS a
literal 0 for them (it serves but does not NLTE-model the line; e.g. Ca 6166). A line that
is identically 0 across ALL (teff,logg,feh) nodes is NOT a physical NLTE correction — it must
be refused, not silently registered as a correction of exactly zero everywhere.
"""
import numpy as np
import pandas as pd
import pytest

from pipeline import nlte_corrections as N
from pipeline.nlte_corrections import detect_placeholder_zero_lines


def _grid(rows):
    return pd.DataFrame(rows, columns=['element', 'wave_A', 'teff_K', 'logg', 'feh', 'delta_nlte'])


def _nodes(wave, vals):
    tlgf = [(5500, 4.0, 0.0), (5772, 4.4, 0.0), (6000, 4.5, 0.0), (5772, 4.4, 0.3)]
    return [['X', wave, t, g, f, v] for (t, g, f), v in zip(tlgf, vals)]


def test_detects_all_zero_line():
    df = _grid(_nodes(6166.44, [0.0, 0.0, 0.0, 0.0]))
    assert detect_placeholder_zero_lines(df) == [6166.44]


def test_ignores_real_varying_line():
    df = _grid(_nodes(5867.57, [0.027, -0.011, 0.022, 0.018]))
    assert detect_placeholder_zero_lines(df) == []


def test_mixed_grid_flags_only_the_placeholder():
    df = _grid(_nodes(5867.57, [0.027, -0.011, 0.022, 0.018]) + _nodes(6166.44, [0.0, 0.0, 0.0, 0.0]))
    assert detect_placeholder_zero_lines(df) == [6166.44]


def test_single_node_not_flagged():
    # a genuinely tiny single-node value is not a placeholder claim (need >=2 nodes)
    df = _grid([['X', 6000.0, 5772, 4.4, 0.0, 0.0]])
    assert detect_placeholder_zero_lines(df) == []


def test_fixed_ca_grid_has_no_placeholder():
    # RYA-413 Part A dropped Ca 6166 -> the registered Ca grid must now load clean.
    N._mpia_element_cache.pop('Ca', None)
    c = N._load_mpia_element_grid('Ca')
    assert 6166.44 not in [round(float(w), 2) for w in c['waves']]
    assert len(c['waves']) >= 5


def test_loader_refuses_a_placeholder_grid(tmp_path, monkeypatch):
    import config.constants as C
    bad = _grid(_nodes(5000.0, [0.05, 0.06, 0.04, 0.05]) + _nodes(6166.44, [0.0, 0.0, 0.0, 0.0]))
    (tmp_path / 'Bad_placeholder.csv').write_text(bad.to_csv(index=False))
    # the loader does `from config.constants import NLTE_CORRECTION_ELEMENTS` per call,
    # so patch the canonical dict there.
    monkeypatch.setitem(C.NLTE_CORRECTION_ELEMENTS, 'Xx',
                        {'ion': 1, 'grid': 'Bad_placeholder.csv', 'ref': 'test'})
    monkeypatch.setattr(N, '_MPIA_ELEMENT_DIR', tmp_path)
    N._mpia_element_cache.pop('Xx', None)
    with pytest.raises(ValueError, match='PLACEHOLDER_ZERO'):
        N._load_mpia_element_grid('Xx')
