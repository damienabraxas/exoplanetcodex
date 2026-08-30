import csv
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HOLDING = ROOT / "data/reference/amarsi2021_cno"
AUDIT = ROOT / "data/audit/rya1136_cno_intake"


def test_ingest_is_reproducible_and_complete():
    subprocess.run(
        [sys.executable, "scripts/ingest_amarsi2021_cno_rya1136.py"],
        cwd=ROOT,
        check=True,
    )
    rows = list(csv.DictReader(
        (HOLDING / "derived/amarsi2021_cno_molecular_lines.csv").open()
    ))
    assert len(rows) == 408
    assert Counter(row["species"] for row in rows) == {
        "12C16O": 80, "C2": 39, "CH": 54, "CN": 59, "NH": 31, "OH": 145,
    }
    assert Counter(row["source_band"] for row in rows) == {
        "VIS": 45, "NIR": 122, "IR": 241,
    }
    assert all(row["canonical_join_status"] == "CROSSMATCH_REVIEW" for row in rows)


def test_coverage_matrix_is_explicit_across_all_reporting_bands():
    rows = list(csv.DictReader((AUDIT / "molecular_coverage_matrix.csv").open()))
    assert len(rows) == 6 * 6
    assert {row["band"] for row in rows} == {
        "FUV", "NUV", "VIS", "RED_OPTICAL", "NIR", "IR",
    }
    assert sum(int(row["used"]) for row in rows) == 408
    assert sum(int(row["ambiguous"]) for row in rows) == 408


def test_summary_refuses_premature_freeze():
    summary = json.loads((AUDIT / "summary.json").read_text())
    assert summary["molecular_used_rows"] == 408
    assert summary["canonical_matched"] == 0
    assert summary["crossmatch_review"] == 408
    assert summary["frozen_ready"] is False
    assert summary["verdict"] == "CROSSMATCH_REVIEW"
