"""
tests/test_grade_weighted_aggregation_rya329.py
===============================================
RYA-329 — line_score-weighted median as primary, GATED on line_score discrimination.

Guards: (1) the weighted-median estimator is correct; (2) the discrimination gate fails
on a flat/circular signal (so the primary is NOT weighted on an unvalidated score) and
passes on a genuinely discriminating one; (3) inert sub-scores are detected.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from config.constants import ISPEC_DIR  # noqa: E402  RYA-1140: not ROOT.parent
sys.path.insert(0, str(ISPEC_DIR))

import pipeline.abundances_derive as ad  # noqa: E402


def test_weighted_median_matches_plain_when_weights_equal():
    v = [7.0, 7.2, 7.4, 7.6, 7.9]
    assert ad._weighted_median(v, [1, 1, 1, 1, 1]) == np.median(v)
    # a heavy weight pulls the crossing point onto that value
    assert ad._weighted_median([7.0, 7.5, 8.0], [1, 1, 50]) == 8.0


def _pool(a, line_score, subscores):
    df = pd.DataFrame({'a_1dlte': a, 'line_score': line_score,
                       'mq_grade': ['MQ-B'] * len(a)})
    for k, vals in subscores.items():
        df[k] = vals
    return df


def test_gate_fails_on_flat_and_circular_signal():
    # line_score varies ONLY via the circular abundance_outlier_score; the independent
    # sub-scores are inert (constant). This is the solar-Fe situation → must NOT weight.
    rng = np.random.default_rng(0)
    a = 7.5 + rng.normal(0, 0.15, 40)
    dev = np.abs(a - np.median(a))
    outlier = 1.0 - dev / dev.max()                 # circular: derived from |a-median|
    ls = 0.5 + 0.3 * outlier
    subs = {'ew_snr_score': np.full(40, 0.6), 'fit_chi2_score': np.full(40, 0.5),
            'saturation_score': np.zeros(40), 'abundance_outlier_score': outlier,
            'nlte_correction_score': np.zeros(40)}
    rep = ad._line_score_discriminates(_pool(a, ls, subs))
    assert rep['discriminates'] is False
    assert set(rep['inert_subscores']) >= {'fit_chi2_score', 'saturation_score',
                                           'nlte_correction_score'}
    # circular ρ looks strong, non-circular ρ is ~0 → the whole point of the gate
    assert rep['rho_noncircular'] is None or abs(rep['rho_noncircular']) < 0.3


def test_gate_passes_when_independent_subscore_discriminates():
    rng = np.random.default_rng(1)
    a = 7.5 + rng.normal(0, 0.15, 40)
    dev = np.abs(a - np.median(a))
    snr = 1.0 - dev / dev.max() + rng.normal(0, 0.02, 40)   # INDEPENDENT signal tracks quality
    ls = np.clip(0.2 + 0.7 * snr, 0, 1)
    subs = {'ew_snr_score': snr, 'fit_chi2_score': rng.uniform(0.3, 0.9, 40),
            'saturation_score': rng.uniform(0.3, 0.9, 40),
            'abundance_outlier_score': rng.uniform(0, 1, 40),
            'nlte_correction_score': rng.uniform(0.3, 0.9, 40)}
    rep = ad._line_score_discriminates(_pool(a, ls, subs))
    assert rep['discriminates'] is True
    assert rep['inert_subscores'] == []


def test_aggregation_diagnostics_emits_all_estimators():
    a = [7.4, 7.5, 7.5, 7.6, 8.2]
    df = _pool(a, [0.7, 0.7, 0.7, 0.7, 0.3],
               {s: np.full(5, 0.5) for s in ad._LS_SUBSCORES})
    d = ad._aggregation_diagnostics(df, delta_nlte=0.01)
    for k in ('unweighted_mean', 'plain_median', 'weighted_mean', 'weighted_median',
              'ab_cut', 'wmed_minus_median', 'plain_median_nlte', 'weighted_median_nlte'):
        assert k in d
    assert d['plain_median'] == 7.5
