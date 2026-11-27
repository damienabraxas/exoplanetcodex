"""RYA-1395 v72: hydrogenic G92RWF writer/reader (scripts/rya1395_hydrogenic_rwfn.py) and the v72 gate inputs."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'scripts'))
import rya1395_hydrogenic_rwfn as H  # noqa: E402

RNT, HH, N = math.exp(-65 / 16) / 8, 0.03125, 439   # the v59+ half grid


@pytest.fixture(scope='module')
def orbs(tmp_path_factory):
    p = tmp_path_factory.mktemp('rwfn') / 'hyd.w'
    H.write_rwfn(p, H.build(['3d', '4f', '5f', '5s'], RNT, HH, N))
    return {H.label(o): o for o in H.read_rwfn(p)}


def test_grid_is_grasp_log_grid():
    r = H.grasp_grid(RNT, HH, N)
    assert r[0] == 0.0
    assert np.allclose(np.log1p(r[1:] / RNT) / HH, np.arange(1, N))


@pytest.mark.parametrize('nl,eps,mean_r', [('3d', 1 / 18, 10.5), ('4f', 1 / 32, 18.0), ('5f', 1 / 50, 31.5), ('5s', 1 / 50, 37.5)])
def test_zeff1_energy_and_mean_radius_are_hydrogenic(orbs, nl, eps, mean_r):
    for lab in (nl, nl + '-') if not nl.endswith('s') else (nl,):
        o = orbs[lab]
        assert o['eps'] == pytest.approx(eps, rel=1e-12)
        # <r> = (3n^2 - l(l+1)) / (2 Z): the relationship, not a stored constant
        assert H.mean_r(o) == pytest.approx(mean_r, rel=2e-3)


def test_round_trip_normalised_and_orthogonal(orbs):
    for o in orbs.values():
        assert H.overlap(o, o) == pytest.approx(1.0, abs=2e-4)
    assert abs(H.overlap(orbs['4f-'], orbs['5f-'])) < 1e-4
    assert abs(H.overlap(orbs['4f'], orbs['5f'])) < 1e-4


def test_kappa_convention_matches_grasp_nak(orbs):
    assert orbs['4f-']['nak'] == 3 and orbs['4f']['nak'] == -4
    assert orbs['3d-']['nak'] == 2 and orbs['3d']['nak'] == -3
    assert orbs['5s']['nak'] == -1


def test_small_component_is_order_alpha(orbs):
    o = orbs['4f']
    assert np.abs(o['Q']).max() < 1e-2 * np.abs(o['P']).max()


def test_bare_z8_reproduces_the_v71_failure():
    """v71: GRASP screened hydrogenic with no core subshells = bare Z=8 (4f eps 2.000, <r> 2.25)."""
    r = H.grasp_grid(RNT, HH, N)
    k, eps, my, p0, P, Q, R = H.orbital(4, 3, False, r, zeff=8.0)
    o = dict(n=4, nak=k, eps=eps, my=my, p0=p0, P=P, Q=Q, R=R)
    assert eps == pytest.approx(2.0)
    assert H.mean_r(o) == pytest.approx(2.25, rel=2e-3)


def test_nist_binding_and_start_zeff():
    import rya1395_mcdhf_s0_v72 as V
    b = V.nist_binding()
    V.BINDING = b
    # Rydberg relationship: start Z_eff gives the NIST binding for the hydrogenic n
    for nl in ('5s', '6s', '3d', '4d', '5d'):
        n = int(nl[:-1])
        assert V.start_zeff(nl) ** 2 / (2 * n * n) == pytest.approx(b[nl], rel=1e-12)
    assert V.start_zeff('4f') == 1.0 and V.start_zeff('5f') == 1.0
    assert b['4f'] == pytest.approx(1 / 32, rel=0.01)       # nf is hydrogenic in O I
    assert b['3s'] > b['4s'] > b['5s'] > b['6s']
