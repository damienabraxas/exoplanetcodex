"""
tests/test_validate_element_rya813.py — RYA-813
===============================================
The per-element validate-and-report stage, and the one property the whole RYA-812
design rests on: THIS PATH NEVER READS GOLD.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.litscan import (  # noqa: E402
    LitscanError, has_litscan, literature_range)
from pipeline.validate_element import (  # noqa: E402
    FAIL, PASS, PASS_WITH_EXCEPTION, REPORT, UN_ANCHORABLE, BandProduct,
    GoldReadInValidationError, appendix, assert_no_gold_read, validate_element)


# ── the load-bearing guarantee ───────────────────────────────────────────────
def test_validation_path_never_reads_gold_statically():
    """
    STATIC proof: neither module in the validation path may even MENTION a gold
    accessor. This is the structural half of RYA-812 principle 5 — the loop
    gold -> verdict -> gold cannot form if the validator cannot name gold.
    """
    forbidden = ("read_solar_reference", "differential_denominator",
                 "data_namespace", "solar_abundances_v")
    for rel in ("pipeline/validate_element.py", "pipeline/litscan.py"):
        src = (ROOT / rel).read_text()
        tree = ast.parse(src)
        # imports are the thing that matters: a docstring may DISCUSS gold (and
        # validate_element's docstring does, at length) without touching it.
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported |= {a.name for a in node.names}
            elif isinstance(node, ast.ImportFrom):
                imported.add(node.module or "")
                imported |= {a.name for a in node.names}
        for bad in forbidden:
            assert not any(bad in i for i in imported), (
                f"{rel} imports {bad!r} — gold is OUTPUT ONLY (RYA-812 principle 5)")


def test_assert_no_gold_read_passes_in_the_real_module():
    assert_no_gold_read()          # must not raise


def test_assert_no_gold_read_actually_fires(monkeypatch):
    """The guard must be able to FAIL, or it proves nothing."""
    import pipeline.validate_element as ve
    monkeypatch.setitem(ve.__dict__, "read_solar_reference", lambda *a, **k: None)
    with pytest.raises(GoldReadInValidationError, match="OUTPUT ONLY"):
        ve.assert_no_gold_read()


# ── Fe: the acceptance case ──────────────────────────────────────────────────
def _fe_products():
    return [
        BandProduct(band="VIS", engine="EW-1D-NLTE+3D", value=7.466,
                    sigma_statistical=0.062, sigma_systematic=0.040,
                    n_lines=62, scale="3D-NLTE",
                    note="reported anchor, Magic-2013 1D->3D applied once (RYA-553)."),
        BandProduct(band="NIR", engine="Engine-A", value=7.58,
                    sigma_statistical=0.11, sigma_systematic=0.09,
                    n_lines=31, scale="3D-NLTE",
                    note="EW-inversion note (RYA-809): 12 flagged problem children "
                         "under per-line RCA; bars are honest, not tuned."),
    ]


def test_fe_vis_validates_and_passes():
    ev = validate_element("Fe", _fe_products(), asplund=7.46)
    vis = [b for b in ev.bands if b.band == "VIS"][0]
    assert vis.tier == "validate"
    assert vis.verdict == PASS, vis.reason
    assert vis.literature_range == (7.42, 7.50)   # RYA-714 §4: 7.46 ± σ_ext 0.04
    assert abs(vis.offset_vs_literature - 0.006) < 1e-9
    assert not ev.un_anchorable


def test_fe_ir_reports_and_is_never_failed():
    """A frontier band outside the literature range is a NOTED GAP, not a failure."""
    ev = validate_element("Fe", _fe_products(), asplund=7.46)
    nir = [b for b in ev.bands if b.band == "NIR"][0]
    assert nir.tier == "report"
    assert nir.verdict == REPORT
    assert nir.excess_beyond_band > 0          # 7.58 IS outside [7.41, 7.51]
    assert "NOT a failure" in nir.reason
    assert nir.verdict != FAIL
    assert ev.failures == []
    # 7.58 is +0.12 from 7.46, i.e. BEYOND 2*sigma_ext -- a DEVIATE-level
    # disagreement. In a frontier band that is still REPORTED, never failed. This
    # is principle 7 made concrete: no gate that punishes being at the frontier.
    assert nir.deviate is True
    assert nir.is_loud is True
    assert "RYA-809" in nir.reason or "RYA-809" in appendix(ev)


def test_ir_carries_honest_two_part_bars():
    ev = validate_element("Fe", _fe_products(), asplund=7.46)
    nir = [b for b in ev.bands if b.band == "NIR"][0]
    assert "stat" in nir.bars and "sys" in nir.bars   # never collapsed into one


# ── the loud cases ───────────────────────────────────────────────────────────
def test_vis_miss_without_a_reason_fails_loudly():
    prods = [BandProduct(band="VIS", engine="E", value=7.90, scale="3D-NLTE")]
    ev = validate_element("Fe", prods, asplund=7.46)
    vis = ev.bands[0]
    assert vis.verdict == FAIL
    assert vis.is_loud
    assert "NO exception was documented" in vis.reason
    assert ev.failures


def test_vis_miss_with_a_documented_reason_is_pass_with_exception():
    prods = [BandProduct(band="VIS", engine="E", value=7.90, scale="3D-NLTE")]
    ev = validate_element("Fe", prods, asplund=7.46,
                          exceptions={"VIS": "blend-dominated pool, RYA-XXX"})
    vis = ev.bands[0]
    assert vis.verdict == PASS_WITH_EXCEPTION
    assert vis.exception_reason
    assert vis.is_loud                      # an exception is never quiet
    assert ev.exceptions and "appendix" in vis.reason.lower()
    assert "DOCUMENTED EXCEPTION" in appendix(ev)


def test_un_anchorable_element_is_report_only_everywhere_and_loud():
    prods = [BandProduct(band="VIS", engine="E", value=1.23),
             BandProduct(band="NIR", engine="E", value=1.30)]
    ev = validate_element("Nh", prods, asplund=None)     # no litscan exists
    assert ev.un_anchorable
    assert {b.verdict for b in ev.bands} == {UN_ANCHORABLE}
    assert all(b.is_loud for b in ev.bands)
    assert ev.failures == []                              # a gap is not a failure
    assert "UN-ANCHORABLE" in appendix(ev)


def test_scale_mismatch_is_surfaced_not_silently_compared():
    prods = [BandProduct(band="VIS", engine="E", value=7.516, scale="1D-NLTE")]
    ev = validate_element("Fe", prods, asplund=7.46)
    vis = ev.bands[0]
    assert vis.scale_mismatch
    assert vis.is_loud
    assert "SCALE MISMATCH" in appendix(ev)


def test_unknown_band_is_reported_not_dropped():
    prods = [BandProduct(band="mid-IR", engine="E", value=7.5)]
    ev = validate_element("Fe", prods, asplund=7.46)
    assert len(ev.bands) == 1
    assert ev.bands[0].verdict == REPORT
    assert "not in the declared vocabulary" in ev.bands[0].reason


# ── litscan contract ─────────────────────────────────────────────────────────
def test_fe_litscan_exists_and_declares_its_derivation_honestly():
    assert has_litscan("Fe")
    lit = literature_range("Fe")
    assert lit is not None
    # RYA-714 §4 ratifies PASS as |delta| <= sigma_ext against a NAMED best
    # external, so the band is Amarsi-2021 7.46 +/- 0.04 -- NOT the FE_GATE policy
    # window [7.41, 7.51]. The two answer different questions and the litscan must
    # not silently inherit the gate.
    assert (lit.min, lit.central, lit.max) == (7.42, 7.46, 7.50)
    assert lit.sigma_external == 0.04
    assert lit.deviate_beyond == 0.08                 # 2*sigma_ext
    assert lit.derivation == "best-external-plus-sigma"
    assert "Amarsi" in lit.best_external
    assert lit.verification_owed
    assert any("Amarsi" in c for c in lit.citations)
    assert any("Caffau" in c for c in lit.citations)  # the full scan, not a stub


def test_band_is_not_the_fe_gate_policy_window():
    """Regression: the first build of this file used FE_GATE (+/-0.05) as the band.

    That validated Fe against a POLICY window while claiming to validate it against
    the literature. The litscan band is agreement-with-literature; FE_GATE is a
    physical acceptance window on the anchor. Different questions."""
    lit = literature_range("Fe")
    assert (lit.min, lit.max) != (7.41, 7.51)


def test_deviate_is_a_distinct_code_from_a_near_miss():
    """RYA-714 §4 codes DEVIATE at >2 sigma; a 1.2 sigma miss is not the same claim."""
    lit = literature_range("Fe")
    assert lit.is_deviate(7.46 + 0.09) is True
    assert lit.is_deviate(7.46 + 0.05) is False


def test_missing_litscan_is_none_not_an_exception():
    """Absent is DATA (un-anchorable), never an error."""
    assert literature_range("Nh") is None


def test_incomplete_range_refuses_rather_than_inventing_a_band(tmp_path, monkeypatch):
    import pipeline.litscan as ls
    monkeypatch.setattr(ls, "LITSCAN_DIR", tmp_path)
    (tmp_path / "Xx.yaml").write_text("element: Xx\nrange:\n  central: 1.0\n  min: 0.9\n")
    with pytest.raises(LitscanError, match="incomplete band"):
        ls.literature_range("Xx")
