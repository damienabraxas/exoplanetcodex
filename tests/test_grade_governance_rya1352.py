"""Characterize the unresolved policy conflict; never silently correct production."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.audit_grade_governance_rya1352 import classifier


def test_current_classifier_preserves_theory_and_evaluated_admission():
    f = classifier((ROOT / "pipeline/gf_grades.py").read_text())
    rows = pd.DataFrame([
        {"gf_tier": "LAB", "gf_sigma_dex": None, "nist_grade": None},
        {"gf_tier": "PR2024-CALC", "gf_sigma_dex": .05, "nist_grade": None},
        {"gf_tier": "NIST-C", "gf_sigma_dex": None, "nist_grade": "C"},
        {"gf_tier": "VALD3", "gf_sigma_dex": None, "nist_grade": None},
        {"gf_tier": "OTHER", "gf_sigma_dex": None, "nist_grade": "D"},
    ])
    assert f(rows).tolist() == [True, True, True, False, False]


def test_snapshot_proves_no_protected_changes():
    import json
    import hashlib
    r = json.loads((ROOT / "docs/audit/rya1352/baseline.json").read_text())
    assert r["pass"] and r["classification_unchanged"]
    assert r["changed_protected_paths"] == []
    assert r["accepted_not_exact_LAB"] > 0  # observed conflict, not an origin inference
    for name, expected in r["protected_file_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    rows = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    mask = classifier((ROOT / "pipeline/gf_grades.py").read_text())(rows)
    ids = "\n".join(rows.loc[mask, "line_id"]).encode()
    assert hashlib.sha256(ids).hexdigest() == r["accepted_line_ids_sha256"]


def test_current_predicate_does_not_validate_sigma():
    # A characterization, not acceptance of invalid scientific uncertainties.
    f = classifier((ROOT / "pipeline/gf_grades.py").read_text())
    rows = pd.DataFrame({"gf_tier": ["OTHER"] * 3,
                         "gf_sigma_dex": [0., -1., float("inf")],
                         "nist_grade": [None] * 3})
    assert f(rows).tolist() == [True, True, True]
