"""RYA-1282: the Bride applies MULTI collision semantics (mul23_opacu.f GENCL), not RH/Lightweaver's.

CH  : coefficient x n(H I) is the DOWNWARD rate; upward by detailed balance.
CH0 : coefficient x n(H I) runs FROM the first written index TO the second, with its
      detailed-balance reverse, interpolated in log(rate).
Lightweaver's Cmat is [destination, source].
"""
import numpy as np
import pytest
from types import SimpleNamespace as S

lw = pytest.importorskip('lightweaver')
from pipeline.multi_atom import MultiCHDownward, MultiCH0Directed  # noqa: E402


class _Pops:
    def __init__(self, nstar, nh0):
        self.atomicPops = {'O': S(nStar=nstar)}
        self._h = np.array([[nh0]])

    def __getitem__(self, key):
        assert key == 'H'
        return self._h


def _atom(n):
    return S(element='O', levels=[S(g=1.0, E=0.0, stage=0)] * n)


def test_CH_coefficient_is_the_downward_rate_and_balances_lte():
    r = MultiCHDownward(j=1, i=0, temperature=[4000., 6000.], rates=[1e-16, 1e-16])
    r.setup(_atom(2))
    m = np.zeros((2, 2, 1)); lte = np.array([[10.], [2.]])
    r.compute_rates(S(temperature=np.array([5000.])), _Pops(lte, 1e20), m)
    np.testing.assert_allclose(m[0, 1], 1e4)          # upper(1) -> lower(0): q * nH
    np.testing.assert_allclose(m[1, 0], 1e4 * 2 / 10)  # upward from detailed balance
    np.testing.assert_allclose(m[0, 1] * lte[1], m[1, 0] * lte[0])


@pytest.mark.parametrize('src,dst', [(0, 2), (2, 0)])
def test_CH0_runs_in_written_direction_with_detailed_balance_reverse(src, dst):
    r = MultiCH0Directed(j=src, i=dst, temperature=[4000., 6000.], rates=[1e-20, 1e-16])
    r.setup(_atom(3))
    m = np.zeros((3, 3, 1)); lte = np.array([[5.], [1.], [3.]])
    r.compute_rates(S(temperature=np.array([5000.])), _Pops(lte, 1e20), m)
    fwd = np.sqrt(1e-20 * 1e-16) * 1e20                # log interpolation midpoint
    np.testing.assert_allclose(m[dst, src, 0], fwd, rtol=1e-12)
    np.testing.assert_allclose(m[src, dst, 0] * lte[dst], m[dst, src, 0] * lte[src], rtol=1e-12)
    assert m[1].sum() == 0 and m[:, 1].sum() == 0     # touches only its own pair


def test_CH0_refuses_non_positive_coefficients_for_log_interpolation():
    r = MultiCH0Directed(j=0, i=1, temperature=[4000., 6000.], rates=[0., 1e-16])
    with pytest.raises(ValueError):
        r.setup(_atom(2))


def test_held_oxygen_atom_converts_CH0_as_ionising_and_keeps_its_order():
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom, to_lightweaver
    try:
        path = require_codex_path('experiments.bride_oxygen_atom')
    except Exception:
        pytest.skip('held O atom not on this host')
    src = read_multi_atom(path)
    model = to_lightweaver(src)
    ch0 = [c for c in model.collisions if isinstance(c, MultiCH0Directed)]
    raw = [c for c in src.collisions if c.kind == 'CH0']
    assert len(ch0) == len(raw) == 24
    for c, r in zip(ch0, raw):
        assert (c.src, c.dst) == (r.j - 1, r.i - 1) and c.upward
        assert model.levels[c.dst].stage == model.levels[c.src].stage + 1
