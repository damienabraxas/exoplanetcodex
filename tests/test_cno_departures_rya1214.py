"""RYA-1214 — C I / O I per-line departure legs for the band-product route.

Amarsi, Nissen & Skuladottir 2019 tabulate two legs per grid line: 1D-NLTE - 1D-LTE
(-> ENGINE-A) and 3D-NLTE - 1D-LTE (-> ENGINE-A-3DNLTE). These pin the three things that
could go wrong silently: a neighbour's correction handed to an off-grid line, a sign flip,
and an element the grid does not tabulate being "served".
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import nlte_cno as nc  # noqa: E402

SUN = dict(teff=5772.0, logg=4.438, feh=0.0, vturb=1.0)


def _deltas(element, waves, A, leg):
    used = [SimpleNamespace(wavelength_air_A=w, abundance=A) for w in waves]
    return nc.departure_deltas(element, "I", used, SUN, leg)


def test_grid_lines_are_served_on_both_legs_with_the_physical_sign():
    for leg in ("1D", "3D"):
        c = _deltas("C", [5052.167, 8335.15], 8.45, leg)
        o = _deltas("O", [7771.94, 7774.17, 7775.39], 8.80, leg)
        assert set(c) == {5052.167, 8335.15}
        assert all(v < 0 for v in c.values())
        # the 777 triplet routes to ONE multiplet-averaged label: identical corrections
        assert len(set(o.values())) == 1 and next(iter(o.values())) < -0.1


def test_the_solar_o_i_777_3d_correction_reproduces_amarsi_2019_table7():
    """Table 7 Sun, 3N - 1L on O I 777: -0.182. The interpolated grid gives -0.171."""
    d = _deltas("O", [7774.17], 8.80, "3D")[7774.17]
    assert d == pytest.approx(-0.182, abs=0.05)


def test_an_off_grid_line_is_NOT_SERVED_not_given_its_neighbours_correction():
    """0.10 A from C I 965.8nm and 0.07 A from 906.1nm: both refused. The resolver's own
    1.5 A default would have served both with a neighbour's value."""
    assert _deltas("C", [9658.53, 9061.5], 8.45, "1D") == {}
    assert nc.resolve_line("CI", 9658.53) is not None      # the default WOULD match
    assert nc.resolve_line("CI", 9658.53, tol_A=nc.CNO_GRID_MATCH_TOL_A) is None


def test_untabulated_species_are_never_served():
    assert _deltas("N", [8216.34], 7.9, "1D") == {}
    used = [SimpleNamespace(wavelength_air_A=5052.167, abundance=8.45)]
    assert nc.departure_deltas("C", "II", used, SUN, "1D") == {}


def test_the_source_string_names_the_leg():
    assert "1D non-LTE" in nc.departure_source("1D")
    assert "3D non-LTE" in nc.departure_source("3D")
