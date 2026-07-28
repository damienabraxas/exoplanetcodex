"""
tests/test_solar_s_costasilva_rya492.py
=======================================
RYA-492 — adopt Costa Silva+2020 atlas-tuned log gf for S I Multiplet 8, and pin the
finding that the gf is NOT the fix for solar S's +0.40 offset. CI-safe (no synthesis):
the abundance re-measure lives in scripts/measure_solar_s_costasilva_rya492.py.
"""
from pathlib import Path

import numpy as np

from pipeline import gf_resolver as gr

ROOT = Path(__file__).resolve().parents[1]

# Costa Silva+2020 Table 1 (arXiv:1912.08659, verified) — 6743 triplet per-component gf.
CS_6743 = [-1.27, -0.95, -0.93]
CS_6743_TOTAL = float(np.log10(np.sum(10.0 ** np.array(CS_6743))))   # -0.5476


def test_canonical_6743_gf_is_costa_silva_total():
    gr._index.cache_clear()
    v = gr.resolve((16, 1), 6743.5466, 7.866)            # S I Mult-8 6743 aggregate
    assert abs(v - CS_6743_TOTAL) < 1e-3                  # banked to the CS triplet total
    assert abs(v - (-0.5476)) < 1e-3


def test_canonical_row_cites_costa_silva():
    row = next(l for l in (ROOT / 'data' / 'linelists' / 'canonical_gf.csv')
               .read_text().splitlines() if l.startswith('gf_113247,'))
    assert 'Costa Silva' in row
    assert 'RYA-492' in row


def test_gf_lever_is_too_small_to_close_the_offset():
    # The finding: adopting the published gf moves the 6743 total by only +0.063 dex
    # (old GES total -0.6103 -> CS -0.5476). A <0.1 dex lever cannot close a 0.40 offset,
    # so S's residual is an astrophysical zero-point (RYA-161), not a gf error.
    old_total = -0.6103
    lever = CS_6743_TOTAL - old_total
    assert 0.0 < lever < 0.10                             # small, strengthening
    assert lever < 0.40                                   # cannot explain the +0.40 offset


def test_solar_S_not_banked_as_pass_in_frozen_reference():
    # S stays curation-owed / LOW_CONFIDENCE in the frozen gold ref — this fix does not
    # bank a solar S value; the residual escalates to RYA-161.
    s_row = next(l for l in (ROOT / 'data' / 'reference' / 'solar' / 'solar_abundances_v1.csv')
                 .read_text().splitlines() if l.startswith('S,I,'))
    assert 'LOW_CONFIDENCE' in s_row
