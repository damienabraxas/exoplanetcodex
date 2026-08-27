from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/audit/rya1061_uv_completeness"


def test_audit_reproduces_committed_outputs(tmp_path):
    subprocess.run([sys.executable, str(ROOT / "scripts/rya1061_uv_completeness_audit.py"), "--out", str(tmp_path)], check=True)
    for path in OUT.iterdir():
        assert (tmp_path / path.name).read_bytes() == path.read_bytes()


def test_canonical_uv_gap_is_explicit():
    summary = json.loads((OUT / "summary.json").read_text())
    assert summary["canonical_by_regime"] == {"FUV": 0, "NUV": 0, "NEAR_UV": 21279}
    assert summary["canonical_by_tier"]["LAB"] == 71
    assert summary["mutation_policy"] == "READ_ONLY_AUDIT_NO_CANONICAL_OR_HOLDING_MUTATIONS"


def test_science_lanes_and_authority_are_not_collapsed():
    papers = pd.read_csv(OUT / "stellar_uv_paper_census.csv")
    assert "PHOTOSPHERIC_ABUNDANCE" in set(papers.science_lane)
    assert "CHROMOSPHERE_TRANSITION_REGION_CORONA" in set(papers.science_lane)
    non_abundance = papers[papers.science_lane != "PHOTOSPHERIC_ABUNDANCE"]
    assert not (non_abundance.authority == "CANDIDATE_ABUNDANCE").any()


def test_source_counts_are_not_misrepresented_as_in_window_counts():
    sources = pd.read_csv(OUT / "laboratory_source_census.csv")
    ward = sources[sources.source == "Ward et al. 2023"].iloc[0]
    assert ward.published_transition_count == 268
    assert "paper" not in str(ward.audit_note).lower() or "window clips" in ward.audit_note
    assert ward.codex_status == "ACQUIRE_TABLE"


def test_holdings_are_derived_from_registry():
    rows = pd.read_csv(OUT / "target_uv_holdings.csv")
    assert {"alpha_cen_a_stis", "alpha_cen_b_stis", "procyon_stis", "procyon_cos"} <= set(rows.holding_id)
    assert set(rows.instrument_id) <= {"hst_stis", "hst_cos"}


def test_procyon_spectrum_without_abundance_paper_is_visible():
    rows = pd.read_csv(OUT / "target_literature_overlap.csv")
    procyon = rows[rows.system_id == "procyon"].iloc[0]
    assert procyon.abundance_literature_status == "NO_PHOTOSPHERIC_ABUNDANCE_PAPER_IN_CENSUS"
