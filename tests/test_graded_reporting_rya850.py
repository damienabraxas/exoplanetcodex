"""
tests/test_graded_reporting_rya850.py — RYA-850
===============================================
Promoting the lab-gf pool to the reported value is a WIRING change, so what needs pinning
is the wiring and the claims it rests on — not a new measurement.

The load-bearing claim is that a graded pool reports a SMALLER total even though RYA-836
measured it scattering WORSE (0.652 vs 0.413) and RYA-842 ratified that gf does not buy a
tighter spread. Both still hold; the resolution is that scatter enters only through
stat = scatter/sqrt(N), which averages down, while the gf term enters sys, which does not.
That is arithmetic, and it is checked per cell rather than assumed — a graded cell that
does not beat its twin must be flagged, not promoted.
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

from pipeline.error_budget import (GRADED_GF_SYSTEMATIC_DEX,  # noqa: E402
                                   UNGRADED_GF_SYSTEMATIC_DEX, build)
from pipeline.graded_reporting import (base_treatment, element_table,  # noqa: E402
                                       format_value, is_graded, pair_products,
                                       stat_from_scatter, total_sigma)

RES = ROOT / "data" / "results" / "rya850"
SUMMARY = RES / "rya850_summary.json"
MATRIX = RES / "rya850_two_product_matrix.csv"


@pytest.fixture(scope="module")
def summary():
    if not SUMMARY.exists():
        pytest.skip("RYA-850 summary absent")
    return json.loads(SUMMARY.read_text())


# ── the uncertainty construction: standard, and nothing more ──────────────────

def test_the_only_combination_is_stat_and_sys_in_quadrature():
    """No new statistics (the ticket is explicit). total = sqrt(stat^2 + sys^2)."""
    assert total_sigma(0.03, 0.04) == pytest.approx(0.05)
    assert stat_from_scatter(0.4, 4) == pytest.approx(0.2)
    assert np.isnan(stat_from_scatter(0.4, 0))


def test_the_scatter_enters_only_through_stat_which_averages_down():
    """THE RESOLUTION OF THE APPARENT CONFLICT WITH RYA-836/842. A pool may scatter worse
    and still report a smaller total, because scatter/sqrt(N) shrinks with N while the gf
    systematic does not shrink at all."""
    worse_scatter_graded = stat_from_scatter(0.652, 59)      # RYA-836 lab pool
    better_scatter_ungraded = stat_from_scatter(0.413, 40)   # RYA-832 Kurucz pool
    assert worse_scatter_graded > better_scatter_ungraded
    g = total_sigma(worse_scatter_graded, GRADED_GF_SYSTEMATIC_DEX)
    u = total_sigma(better_scatter_ungraded, UNGRADED_GF_SYSTEMATIC_DEX)
    assert g < u, "the graded pool no longer reports a smaller total — re-derive RYA-850"


def test_the_two_gf_terms_are_the_ones_already_in_the_budget():
    """This ticket wires terms that existed; it does not invent them."""
    assert GRADED_GF_SYSTEMATIC_DEX == 0.041
    assert UNGRADED_GF_SYSTEMATIC_DEX == 0.17
    b_g = build("Fe", 5000.0, 10, scatter_dex=0.1, gf_graded=True,
                harness_residual_dex=0.0129, handler="ProfileFitHandler")
    b_u = build("Fe", 5000.0, 10, scatter_dex=0.1, gf_graded=False,
                harness_residual_dex=0.0129, handler="ProfileFitHandler")
    assert b_g.total()[1] < b_u.total()[1]


# ── the pairing rule ──────────────────────────────────────────────────────────

def test_a_graded_product_is_recognised_by_its_pool_not_its_engine():
    """RYA-712 keys a product on what produced it; `1D-LTE` and `1D-LTE-LABGF` are the
    SAME engine on DIFFERENT pools, which is exactly what makes them comparable."""
    assert is_graded("1D-LTE-LABGF") and not is_graded("1D-LTE")
    assert base_treatment("1D-LTE-LABGF") == base_treatment("1D-LTE") == "1D-LTE"
    assert base_treatment("ENGINE-B-NLTE") == "ENGINE-B-NLTE"


def test_a_cell_with_no_lab_gf_pool_falls_back_and_says_so():
    """15 of 18 cells have no lab-gf pool measured. They must NOT be relabelled graded on
    the strength of a term they are not entitled to."""
    df = pd.DataFrame([
        dict(element="Fe", ion="II", band="VIS", treatment="ENGINE-A", A=7.5,
             n_lines=3, stat_dex=0.08, syst_dex=0.17),
    ])
    p = pair_products(df)[0]
    assert p.graded is None and p.primary_is_graded is False
    assert p.primary is not None, "the ungraded product is KEPT, never dropped"
    assert p.graded_beats_ungraded is None


def test_the_ungraded_product_is_never_discarded(summary):
    """An ungraded line is usually good data with a worse-known gf (the ticket is
    explicit). Every cell keeps both pools it has."""
    if not MATRIX.exists():
        pytest.skip("matrix absent")
    d = pd.read_csv(MATRIX)
    assert (~d.treatment.map(is_graded)).sum() >= 15


# ── the reported result ───────────────────────────────────────────────────────

def test_every_graded_cell_actually_tightens_its_bar(summary):
    """Measured per cell, not assumed. A graded cell that WIDENED the bar would mean its
    scatter penalty outran the gf gain, and promoting it would publish a worse number."""
    cells = summary["graded_cells"]
    assert len(cells) == 3
    for c in cells:
        assert c["graded_beats_ungraded"] is True, (
            f"{c['band']}: promoting the graded pool widens the bar "
            f"({c['ungraded_total_dex']} -> {c['total_dex']}) — do not promote it")
        assert c["total_dex"] < c["ungraded_total_dex"]


def test_the_graded_bands_agree_with_each_other(summary):
    """If the graded bands disagreed by more than their bars, 'the graded value' would not
    be well defined and picking one would be a choice dressed as a measurement."""
    for c in summary["graded_consistency"]:
        assert c["n_sigma"] < 2.0, (
            f"{c['a']} vs {c['b']} now differ by {c['n_sigma']:.1f} sigma — the graded "
            f"cells no longer agree and the headline choice needs revisiting")


def test_no_combined_cross_band_value_is_computed(summary):
    """RYA-712: separate products, never combined. `band_products` has no combine() and
    this reporting layer must not smuggle one in."""
    bp = (ROOT / "pipeline" / "band_products.py").read_text()
    assert "def combine(" not in bp
    gr = (ROOT / "pipeline" / "graded_reporting.py").read_text()
    assert "def combine(" not in gr
    assert "graded_band_spread_dex" in summary, (
        "the band-to-band spread must be REPORTED rather than collapsed")


def test_the_generic_graded_term_understates_the_measured_lab_sigma(summary):
    """The cited laboratory sigma of EVERY graded pool exceeds the generic NIST grade-B
    bound, so wiring the generic term publishes a bar tighter than the pool's own
    oscillator-strength uncertainty supports.

    RYA-853 refereed this table against the source papers — Ruffoni 142/142 and Den Hartog
    203/203 perfect on value AND cited sigma — so the cited figure is audited data rather
    than a plausible-looking alternative. It is derived from each pool's own lines here,
    not hardcoded, and reproduces RYA-824's published 0.0524 / 0.0600."""
    cited = summary["cited_gf_sigma_by_band"]
    assert set(cited) >= {"VIS", "red-optical", "near-UV"}, (
        "a graded band lost its cited-sigma derivation")
    for band, c in cited.items():
        assert c["cited_sigma"] > GRADED_GF_SYSTEMATIC_DEX, (
            f"{band}: the cited lab sigma no longer exceeds the generic term")
    assert cited["red-optical"]["cited_sigma"] == pytest.approx(0.0524, abs=5e-4)
    assert cited["VIS"]["cited_sigma"] == pytest.approx(0.0600, abs=5e-4)


def test_the_cited_sigma_variant_is_reported_for_every_graded_cell(summary):
    """All three graded cells get the comparison, not just the two wired from RYA-824.
    The near-UV moves only +2.6% because its 0.100 dex continuum term dominates — which is
    itself worth seeing, since it says where the choice actually matters."""
    cells = {c["band"]: c for c in summary["cells_with_cited_sigma"]}
    assert set(cells) == {"VIS", "red-optical", "near-UV"}
    for c in cells.values():
        assert c["total_cited"] > c["total_generic"], (
            f"{c['band']}: the cited sigma no longer widens the bar — re-derive")
    assert cells["near-UV"]["pct"] < 5.0 < cells["VIS"]["pct"]


def test_the_pre_847_caveat_travels_with_the_numbers(summary):
    """The graded pool becomes lab-gf INTERSECT constrained once RYA-847 lands, so these
    values are provisional. A number that outlives its caveat is how a provisional result
    becomes a published one."""
    joined = " ".join(summary["caveats"]).lower()
    assert "847" in joined
    assert "provisional" in joined or "regenerated" in joined


def test_the_plus_minus_is_visible_in_the_reported_string():
    """Spec item 4: the main table shows the value WITH its uncertainty."""
    assert format_value(7.577, 0.1375) == "7.577 +/- 0.138"
    assert format_value(7.577, float("nan")).endswith("n/a")
