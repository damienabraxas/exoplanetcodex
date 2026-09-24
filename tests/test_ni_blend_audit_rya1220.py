"""RYA-1220 -- the N I diagnostic set must be selected on measured contamination.

RYA-1220 expanded the native N I set to five lines and reported a five-line median.
The two it added are the two most contaminated of the five, both by CN -- the competing
N diagnostic. These tests pin the measurement and the circularity rule, so the set
cannot be widened again without someone re-measuring what sits next to the lines.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "rya1220_ni_blend_audit.py"
ARTIFACT = ROOT / "data" / "output" / "rya1220" / "ni_blend_contamination.json"


def _mod():
    spec = importlib.util.spec_from_file_location("rya1220_blend", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_five_line_set_is_the_amarsi_selection():
    assert _mod().NI_SET == [7442.298, 7468.312, 8216.336, 8629.235, 8683.403]


def test_contamination_is_measured_from_the_preserved_extract_not_the_built_list():
    """The built linelist_*.csv has no neighbours and no term/J; the question is neighbours."""
    src = SCRIPT.read_text()
    assert "vald_raw_archive" in src, "audit no longer reads the preserved VALD extract"
    assert "linelist_solar.csv" not in src, (
        "the built line list cannot answer a neighbour question -- it carries no neighbours")


def test_ratio_is_summed_depth_not_a_nearest_neighbour_count():
    """RYA-1189: rank isolation by summed catalogued depth; a count hides a deep blend."""
    m = _mod()
    rows = [{"species": "N 1", "wl": 8629.235, "depth": 0.004},
            {"species": "CN 1", "wl": 8629.2426, "depth": 0.003},
            {"species": "CN 1", "wl": 8629.2100, "depth": 0.001}]
    got = m.audit(rows)[3]
    assert got["blend_depth_ratio"] == pytest.approx(1.0, abs=1e-6)  # (0.003+0.001)/0.004
    assert got["n_contaminants"] == 2


def test_a_cn_blended_line_cannot_adjudicate_against_cn():
    m = _mod()
    rows = [{"species": "N 1", "wl": 8629.235, "depth": 0.004},
            {"species": "CN 1", "wl": 8629.2426, "depth": 0.003}]
    got = m.audit(rows)[3]
    assert got["cn_circular"] is True
    assert got["admit"] is False


def test_an_iron_blend_is_not_circular():
    """Only CN contamination is circular for an N I vs CN adjudication."""
    m = _mod()
    rows = [{"species": "N 1", "wl": 7468.312, "depth": 0.005},
            {"species": "Fe 1", "wl": 7468.2857, "depth": 0.002}]
    got = m.audit(rows)[1]
    assert got["cn_circular"] is False
    assert got["admit"] is True


@pytest.mark.skipif(not ARTIFACT.exists(), reason="audit artifact not generated")
def test_admitted_set_equals_the_production_registry():
    """The measured-clean set must be the set production actually runs."""
    doc = json.loads(ARTIFACT.read_text())
    assert doc["admitted"] == [7468.312, 8216.336, 8683.403]
    assert doc["rejected"] == [7442.298, 8629.235]

    import sys
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from pipeline import pysme_nlte
    registered = sorted(row[0] for row in pysme_nlte.NLTE_LINES["N"])
    assert registered == [7468.31, 8216.34, 8683.40], (
        "NLTE_LINES['N'] drifted from the measured-clean set")
    for a, b in zip(sorted(doc["admitted"]), registered):
        assert abs(a - b) < 0.01, f"admitted {a} does not match registered {b}"
