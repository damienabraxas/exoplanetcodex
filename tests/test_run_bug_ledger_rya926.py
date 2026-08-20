"""RYA-926 — the committed defect ledger has one stable, honest schema."""
from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "audit" / "run_bug_ledger.csv"

REQUIRED = {
    "ticket_id", "run_id", "discovered_at_utc", "element", "system",
    "instrument", "holding_id", "product_identity", "stage", "severity",
    "category", "symptom", "expected_behavior", "actual_behavior",
    "root_cause", "detected_by", "affected_files", "affected_artifacts",
    "science_impact", "status", "fix_ticket_id", "regression_test",
    "base_sha", "fix_sha", "notes_provenance",
}


def _rows():
    with LEDGER.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return reader.fieldnames, list(reader)


def test_schema_and_rows_are_well_formed():
    fields, rows = _rows()
    assert fields is not None
    assert set(fields) == REQUIRED
    assert rows, "historical Kitt Peak→HARPS regression seed is required"
    for row in rows:
        assert None not in row, "row has more cells than the canonical header"
        assert re.fullmatch(r"RYA-\d+", row["ticket_id"])
        assert row["severity"] in {"CRITICAL", "WARNING", "SUGGESTION"}
        assert row["status"] in {
            "open", "characterized", "fixed", "wont_fix",
            "expected_behavior", "duplicate",
        }
        assert row["symptom"] and row["science_impact"]


def test_cross_instrument_regression_is_seeded_without_fabricated_fields():
    _, rows = _rows()
    seeds = [r for r in rows if r["category"] == "wrong_holding"
             and r["instrument"] == "harps"]
    assert len(seeds) == 1
    seed = seeds[0]
    assert seed["fix_ticket_id"] == "RYA-913"
    assert seed["regression_test"] == "tests/test_instrument_agnostic_rya913.py"
    assert seed["run_id"] == ""
    assert seed["base_sha"] == ""
    assert "UNKNOWN" in seed["affected_artifacts"]
