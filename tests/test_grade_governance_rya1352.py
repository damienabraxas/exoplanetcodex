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


def test_snapshot_records_no_protected_changes():
    import json
    r = json.loads((ROOT / "docs/audit/rya1352/baseline.json").read_text())
    assert r["pass"] and r["classification_unchanged"]
    assert r["changed_protected_paths"] == []
    assert r["accepted_not_exact_LAB"] > 0  # observed conflict, not an origin inference


def test_audit_detects_changed_product_and_classifier(tmp_path, monkeypatch):
    import json
    import subprocess
    from scripts import audit_grade_governance_rya1352 as audit_module
    for directory in ("pipeline", "data/linelists", "data/products"):
        (tmp_path / directory).mkdir(parents=True)
    code = ROOT / "pipeline/gf_grades.py"
    (tmp_path / "pipeline/gf_grades.py").write_text(code.read_text())
    pd.DataFrame([{"line_id": "test", "gf_tier": "LAB", "gf_sigma_dex": .1,
                   "nist_grade": None}]).to_csv(tmp_path / "data/linelists/canonical_gf.csv", index=False)
    product = tmp_path / "data/products/test.json"
    product.write_text(json.dumps({"products": [{"grade": "Codex Grade", "A": 1.}]}))
    def git(*args):
        subprocess.run(["git", *args], cwd=tmp_path, check=True, capture_output=True)
    git("init")
    git("add", ".")
    git("-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-m", "fixture")
    monkeypatch.setattr(audit_module, "ROOT", tmp_path)
    assert audit_module.audit("HEAD")["pass"]
    product.write_text(json.dumps({"products": [{"grade": "Deep Grade", "A": 2.}]}))
    result = audit_module.audit("HEAD")
    assert not result["pass"]
    assert result["classification_unchanged"]
    assert result["changed_protected_paths"] == ["data/products/test.json"]
    source = tmp_path / "pipeline/gf_grades.py"
    source.write_text(source.read_text().replace('tier.str.contains("LAB", na=False)', 'tier.str.contains("NEVER", na=False)').replace('sig.notna()', 'sig.isna()'))
    result = audit_module.audit("HEAD")
    assert not result["classification_unchanged"]


def test_current_predicate_does_not_validate_sigma():
    # A characterization, not acceptance of invalid scientific uncertainties.
    f = classifier((ROOT / "pipeline/gf_grades.py").read_text())
    rows = pd.DataFrame({"gf_tier": ["OTHER"] * 3,
                         "gf_sigma_dex": [0., -1., float("inf")],
                         "nist_grade": [None] * 3})
    assert f(rows).tolist() == [True, True, True]
