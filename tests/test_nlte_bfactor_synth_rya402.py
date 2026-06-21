"""
tests/test_nlte_bfactor_synth_rya402.py
=======================================
RYA-402 — the testable core of the b-factor NLTE-synthesis capability: the Amarsi
grid reader, the b-factor->delta extraction math, and the fail-loud/no-fake
contract. The live Turbospectrum synthesis + Na/Al/K results are gated on Ryan's
Amarsi-2020 fetch and are validated by the --validate-against Na STOP-gate, not here.
"""
from pathlib import Path

import numpy as np
import pytest

from pipeline import nlte_bfactor_synth as bf


# ── b-factor -> delta extraction (pure math) ─────────────────────────────────

def test_delta_recovers_known_shift_identity_cog():
    # ew_lte_at(A) = A  (monotonic stub). NLTE EW above the LTE EW at a_lte0 by +0.107
    # must yield delta = +0.107; below by 0.107 -> -0.107 (the INSPECT sign).
    assert bf.delta_from_ews(lambda A: A, ew_nlte=5.107, a_lte0=5.0) == pytest.approx(0.107, abs=1e-3)
    assert bf.delta_from_ews(lambda A: A, ew_nlte=4.893, a_lte0=5.0) == pytest.approx(-0.107, abs=1e-3)


def test_delta_with_curved_cog():
    # a realistic monotone EW(A) = 10**(A-7): recover a known delta
    ew = lambda A: 10.0 ** (A - 7.0)
    a0 = 6.20
    target = 6.40
    delta = bf.delta_from_ews(ew, ew_nlte=ew(target), a_lte0=a0, a_lo=5.0, a_hi=7.5)
    assert delta == pytest.approx(target - a0, abs=2e-3)


def test_delta_unbracketed_raises():
    with pytest.raises(ValueError, match="not bracketed"):
        bf.delta_from_ews(lambda A: A, ew_nlte=99.0, a_lte0=5.0, a_lo=4.0, a_hi=6.0)


# ── Amarsi grid reader (synthetic fixture) ───────────────────────────────────

def _write_fixture(d: Path, element: str, nodes, b_per_node):
    (d / f'label_{element}.txt').write_text(
        "# teff logg feh a_x\n" + "\n".join(" ".join(f"{x:g}" for x in row) for row in nodes) + "\n")
    out = []
    for arr in b_per_node:           # arr shape (nlevel, ndepth)
        nl, nd = arr.shape
        out.append(f"{nd} {nl}")
        out.append(" ".join(f"{v:g}" for v in arr.flatten()))
    (d / f'nlte_{element}_test.txt').write_text("\n".join(out) + "\n")


def test_reader_parses_synthetic_grid(tmp_path):
    nodes = [[5772, 4.44, 0.0, 6.43], [5000, 4.5, -0.5, 6.0]]
    b0 = np.array([[1.0, 1.1, 0.9], [1.0, 0.8, 1.2]])     # (nlevel=2, ndepth=3)
    b1 = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])     # unity node
    _write_fixture(tmp_path, 'Al', nodes, [b0, b1])
    g = bf.read_amarsi_grid('Al', base_dir=tmp_path)
    assert g.element == 'Al' and g.n_nodes == 2
    assert g.nlevel == 2 and g.ndepth == 3
    # nearest node snap + departures
    np.testing.assert_allclose(g.departures(5772, 4.44, 0.0), b0)
    np.testing.assert_allclose(g.departures(5000, 4.5, -0.5), b1)


def test_reader_node_outside_tol_raises(tmp_path):
    _write_fixture(tmp_path, 'K', [[5772, 4.44, 0.0, 5.07]],
                   [np.ones((2, 3))])
    g = bf.read_amarsi_grid('K', base_dir=tmp_path)
    with pytest.raises(ValueError, match="outside grid tol"):
        g.departures(3000, 1.0, 0.0)            # far from the only node


def test_reader_missing_files_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="not found"):
        bf.read_amarsi_grid('Al', base_dir=tmp_path)


# ── fail-loud / no-fake contract ─────────────────────────────────────────────

def test_live_synth_fails_loud_without_depgrid():
    # No dep-grid installed -> RuntimeError, never a silent LTE number.
    with pytest.raises((RuntimeError, NotImplementedError)):
        bf.synth_ew_nlte_vs_lte('Al', 3961.52, {'teff': 5772, 'logg': 4.44, 'feh': 0.0})


def test_validate_against_stops_without_data():
    # Na anchor exists, but the grid isn't fetched -> loud stop (no fabricated pass).
    with pytest.raises((RuntimeError, FileNotFoundError)):
        bf.validate_against('Na')


def test_known_delta_anchor_is_the_inspect_na():
    assert bf._KNOWN_DELTA['Na']['delta'] == pytest.approx(-0.107)
    assert bf.TARGET_ELEMENTS == ['Al', 'K']
