from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import audit_si_uncertainty_rya1218 as audit


def test_si_uncertainty_audit_keeps_all_components_on_hold():
    doc = audit.audit(ROOT / "data/results/rya1218")
    assert len(doc["products"]) >= 20
    assert doc["admitted_products"] == 0
    assert len(doc["components"]) == 16
    assert {c["state"] for c in doc["components"]} == {"HOLD"}


def test_si_uncertainty_audit_reports_real_line_scatter():
    doc = audit.audit(ROOT / "data/results/rya1218")
    rows = [p for p in doc["products"] if "crires_abundance_run" in p["artifact"] and p["accepted_line_count"] > 1]
    assert rows
    assert all(p["line_scatter_dex"] is not None for p in rows)
    assert all(p["scatter_sem_dex"] is not None for p in rows)
