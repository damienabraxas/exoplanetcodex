"""Bride QA (RYA-1282): the physics grades pass on the solar benchmark AND catch the bugs they exist for.

Positive: Gerber+2023 node p5750/4.5/A8.70 through the Bride default path.
Negative controls re-inject each fixed defect; a grade that cannot fail is not a test:
  * RH/Lightweaver CH (coefficient read as UPWARD)  -> Q6 triplet departures vs the benchmark
    must fail. NOT Q2: a CH misread is still detailed-balanced, it only over-thermalises,
    so the LTE limit cannot see it (measured: Q2 9.7e-7 with the bug). The direction itself
    is fixed by the MULTI format contract (mul23_opacu.f GENCL) and pinned in
    test_bride_collision_semantics_rya1282.py.
  * CH0 one-way, no detailed-balance reverse        -> Q5 detailed balance AND Q2 must fail

Bounds, and why (not tuned to pass):
  Q2 <= 1e-5 : the solve stops at 1e-6 relative change; 10x that is the declared ceiling
               for "deep layers are in LTE". Measured 9.7e-7.
  Q3 <= 1e-12: particle number is an exact constraint of the linear SE solve (float64).
  Q4 <= 1e-5 : physical rate residual after a FRESH formal solution, same 10x-rtol basis.
               Measured 4.9e-7.
  Q5 <= 1e-12: detailed balance is built into every collision class algebraically.
Requires the held O atom and MARCS grid (Mac/Sirius); skipped elsewhere.
"""
from pathlib import Path
import pytest

pytest.importorskip('lightweaver')
from config.constants import codex_path  # noqa: E402

BENCH = codex_path('grids.marcs_standard')  # presence probe only
RECORD = Path.home() / 'codex/bride_inputs/rya1282/bench_gerber2023_O_MARCS/p5750_g+4.5_+0.00_a+0.00__A8.7.json'
pytestmark = pytest.mark.skipif(not RECORD.exists(), reason='benchmark record not on this host')


def _grade():
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
    from bride_qa_benchmark import grade
    return grade(str(RECORD))


def test_solar_benchmark_passes_every_physics_grade():
    g = _grade()
    assert g['Q1_converged'], g.get('error')
    assert g['Q2_lte_limit_max_abs_b_minus_1'] <= 1e-5
    assert g['Q3_conservation_max'] <= 1e-12
    assert g['Q4_rate_residual_fresh'] <= 1e-5
    assert g['Q5_detailed_balance_max_rel'] <= 1e-12
    assert all(r > 1.0 for r in g['nlte_over_lte'])        # 777 strengthens in NLTE


def test_QA_catches_the_RH_CH_direction_bug(monkeypatch):
    import pipeline.multi_atom as MA
    from lightweaver.collisional_rates import CH
    monkeypatch.setattr(MA, 'MultiCHDownward', CH)
    g = _grade()
    assert not g['Q1_converged'] or g['Q6_triplet_phot_max_dlogb'] > 0.05   # fixed: 0.0047


def test_QA_catches_one_way_CH0(monkeypatch):
    import numpy as np
    import pipeline.multi_atom as MA

    class OneWay(MA.MultiCH0Directed):
        def compute_rates(self, atmos, eqPops, Cmat):
            q = np.exp(np.interp(atmos.temperature, self.temperature, self.log_rates))
            Cmat[self.dst, self.src, :] += q * eqPops['H'][0, :]
    OneWay.__module__ = MA.__name__
    monkeypatch.setattr(MA, 'MultiCH0Directed', OneWay)
    g = _grade()
    assert not g['Q1_converged'] or (g['Q5_detailed_balance_max_rel'] > 0.5
                                     and g['Q2_lte_limit_max_abs_b_minus_1'] > 0.05)
