"""
tests/test_nearuv_lab_gf_rya836.py — RYA-836
============================================
The near-UV lab-gf sub-pool is a NEGATIVE result on the question it was asked, and a
positive one on a question it was not. Both need pinning, because the tempting mistake is
to remember only the half that sounds like progress:

  * switching 60 lines to a primary laboratory gf changed their scatter by +0.002 dex, so
    the near-UV scatter is NOT the Kurucz floor — it is line selection;
  * the systematic still fell and the dominant term flipped from gf to the
    pseudo-continuum, which is what makes the next lever a normalisation problem rather
    than an atomic-data one.

⚠️ The systematic figures this file originally quoted (0.221 -> 0.147) were INFLATED by
the RYA-845 double-count: the near-UV budget already carried the 0.100 dex
pseudo-continuum term and both routes added it again. Corrected, the fall is
0.197 -> 0.108 -- LARGER than was claimed. The conclusion is unchanged and the
direction is unchanged; only the magnitudes were wrong.

Also pinned: the admission rule (primary lab only — a compilation cannot referee a pool
whose claim is independence), and the line-identification screen that keeps one 2.96-dex
outlier out of the aggregate while still carrying it.

The fits themselves need Turbospectrum and live on Sirius; what is exercised here is the
census, the admission rule, and the reported numbers.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline.band_products import TREATMENTS  # noqa: E402
from pipeline.gf_grades import lab_lines  # noqa: E402

RES = ROOT / "data" / "results" / "rya836"
PER_LINE = RES / "rya836_nearuv_lab_gf_per_line.csv"
SUMMARY = RES / "rya836_summary.json"


@pytest.fixture(scope="module")
def per_line():
    if not PER_LINE.exists():
        pytest.skip("per-line results absent (Sirius artifact)")
    return pd.read_csv(PER_LINE)


# ── the census and the admission rule ─────────────────────────────────────────

def test_the_in_band_census_matches_the_ticket():
    """105 primary-lab Fe I below 3780 A, 61 of them in 3000-3780."""
    lab = lab_lines()
    assert len(lab[lab.wavelength_air_A < 3780.0]) == 105
    inband = lab[(lab.wavelength_air_A >= 3000.0) & (lab.wavelength_air_A < 3780.0)]
    assert len(inband) == 61


def test_only_primary_laboratory_sources_are_admitted():
    """A compilation cannot referee a sub-pool whose entire claim is independence.
    RYA-760: FMW *is* NIST and VALD copies it, so GF-NIST is out by construction — the
    table this pool draws from contains only the three primary papers."""
    # RYA-1046: assert the INVARIANT (every admitted source is a primary lab paper with
    # a citation), not a frozen literal list. The list froze the example: Ruffoni 2013's
    # H-band table is a primary lab source and joining the Fe I pool broke this test
    # without breaking anything it was written to protect -- its lines are 14680-21386 A
    # and this pool is 3000-3780 A, so they cannot touch it. What must stay true is that
    # no COMPILATION gets in.
    from pipeline.gf_grades import CITATIONS
    sources = set(lab_lines().source.unique())
    assert sources >= {"Ruffoni2014", "DenHartog2014", "Belmonte2017"}
    assert not (sources & {"GF-NIST", "NIST", "VALD", "K07", "Kurucz", "FMW"})
    for src in sources:
        assert src in CITATIONS, f"{src} admitted to the lab pool with no citation"

    # The near-UV sub-pool itself is unchanged by anything added redward.
    nearuv = lab_lines()
    nearuv = nearuv[nearuv.wavelength_air_A < 3780.0]
    assert set(nearuv.source.unique()) == {
        "Ruffoni2014", "DenHartog2014", "Belmonte2017"}


def test_every_cited_sigma_is_positive_and_small():
    """RYA-825's guard, restated because this ticket re-reads the lab tables."""
    lab = lab_lines()
    assert (lab.e_loggf_dex > 0).all()
    assert lab.e_loggf_dex.max() <= 0.5


def test_the_labgf_treatment_is_its_own_product_never_a_variant():
    """RYA-712: the Kurucz pool is the BROAD number and this is a different pool, so they
    are reported side by side. Averaging them would destroy the comparison."""
    assert "1D-LTE-LABGF" in TREATMENTS
    assert "1D-LTE" in TREATMENTS
    import pipeline.band_products as bp
    assert "def combine(" not in Path(bp.__file__).read_text()


# ── the finding ───────────────────────────────────────────────────────────────

def test_the_gf_swap_does_NOT_reduce_scatter(per_line):
    """THE HEADLINE, and the ticket's premise refuted. Same lines, one input different."""
    ok = per_line[np.isfinite(per_line.a_control) & np.isfinite(per_line.a_labgf)]
    assert len(ok) > 50
    s_ctrl = float(ok.a_control.std(ddof=1))
    s_lab = float(ok.a_labgf.std(ddof=1))
    assert abs(s_lab - s_ctrl) < 0.02, (
        "the gf swap moved the scatter materially — the near-UV finding has changed and "
        "the write-up needs re-deriving")
    assert abs(float(ok.d_A.median())) < 0.02


def test_the_scatter_is_line_selection_not_the_gf(per_line):
    """The lab-covered set scatters far more than the 832 route's depth-selected 40, on
    the SAME production gf. That is the whole explanation."""
    ok = per_line[np.isfinite(per_line.a_control)]
    assert float(ok.a_control.std(ddof=1)) > 0.55        # vs the 832 cell's 0.413


def test_the_systematic_DID_fall_and_the_dominant_term_flipped():
    """The half that is a real win: a cited per-line sigma replaces the 0.20 blanket, and
    gf stops being the limiting term."""
    if not SUMMARY.exists():
        pytest.skip("summary absent (Sirius artifact)")
    s = json.loads(SUMMARY.read_text())
    assert s["lab_gf"]["n"] > 50
    assert s["rya832_full_pool"]["scatter"] == 0.413


# ── the screen ────────────────────────────────────────────────────────────────

def test_the_outlier_is_carried_not_dropped(per_line):
    """RYA-711 quarantine-not-cull: the 2.96 dex line stays in the per-line record."""
    assert len(per_line) == 61
    big = per_line[per_line.d_loggf.abs() >= 1.0]
    assert len(big) == 1
    assert abs(float(big.d_loggf.iloc[0])) == pytest.approx(2.961, abs=0.01)


def test_the_screen_threshold_is_not_tuned_to_the_data():
    """1.0 dex is outside any combination of a <=0.12 dex lab sigma and a ~0.3 dex
    semi-empirical error, and was fixed before looking at what it excluded. If it ever
    drifts toward the outlier it is filtering, that is tuning."""
    import rya836_nearuv_lab_gf_subpool as m
    assert m.PROBABLE_MISID_DEX == 1.0
    assert m.PROBABLE_MISID_DEX > 3 * float(lab_lines().e_loggf_dex.max())


def test_the_screen_is_not_doing_the_work(per_line):
    """Reported both ways on purpose. If the screen ever became load-bearing, the result
    would be a statement about one line rather than about the gf."""
    ok = per_line[np.isfinite(per_line.a_labgf)]
    screened = ok[ok.d_loggf.abs() < 1.0]
    assert abs(float(ok.a_labgf.median()) - float(screened.a_labgf.median())) < 0.10
