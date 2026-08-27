from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/audit/rya1062_uv_ion_roadmap"


def test_outputs_reproduce(tmp_path):
    subprocess.run([sys.executable, str(ROOT / "scripts/rya1062_uv_ion_roadmap_audit.py"), "--out", str(tmp_path)], check=True)
    for p in OUT.iterdir():
        assert (tmp_path / p.name).read_bytes() == p.read_bytes()


def test_read_only_policy_and_required_gap_classes():
    s = json.loads((OUT / "summary.json").read_text())
    assert s["mutation_policy"] == "READ_ONLY_NO_STATE_CHANGES"
    rec = pd.read_csv(OUT / "uv_source_reconciliation.csv")
    assert {"SOURCE_KNOWN_NOT_INGESTED", "NEW_LAB_SOURCE", "NEW_EVALUATED_SOURCE", "ROADMAP_ION_MISSING"} <= set(rec.gap_class)


def test_primary_and_evaluated_sources_stay_distinct():
    rec = pd.read_csv(OUT / "uv_source_reconciliation.csv")
    v = rec[rec.ion == "V II"]
    assert "PRIMARY_LAB" in set(v.evidence_type)
    assert "EVALUATED" in set(v.evidence_type)
    assert len(set(v.evidence_type)) == 2


def test_known_label_mismatches_are_carried_not_fixed():
    m = pd.read_csv(OUT / "ion_label_mismatches.csv")
    assert {"Sr", "Y", "Zr", "Ba", "Eu"} == set(m.element)
    assert "Fe" not in set(m.element)  # Fe II is an arbiter beside the Fe I headline.
    assert set(m.gap_class) == {"ION_LABEL_MISMATCH"}


def test_missing_ions_require_feasibility_before_promotion():
    m = pd.read_csv(OUT / "roadmap_missing_ions.csv")
    assert {"P", "S", "Ni", "Cu", "Zn"} == set(m.element)
    assert set(m.gap_class) == {"ROADMAP_ION_MISSING"}


def test_crii_opportunity_is_not_claimed_as_ingested():
    rec = pd.read_csv(OUT / "uv_source_reconciliation.csv")
    ward = rec[(rec.ion == "Cr II") & (rec.source == "Ward et al.")].iloc[0]
    assert ward.published_transition_count == 268
    assert ward.already_represented_in_canonical_gf == "no"
    assert ward.gap_class == "SOURCE_KNOWN_NOT_INGESTED"
