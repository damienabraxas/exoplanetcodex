import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/audit/rya1169_si_intake"


def rows(name):
    with (OUT / name).open(newline="") as f:
        return list(csv.DictReader(f))


def setup_module():
    subprocess.run([sys.executable, str(ROOT / "scripts/build_si_intake_rya1169.py")], check=True)


def test_reference_membership_and_negative_selections():
    data = rows("si_agss21_reference_lines.csv")
    assert len(data) == 12
    assert sum(r["reference_status"] == "used" for r in data) == 10
    assert sum(r["reference_status"] == "explicitly_rejected" for r in data) == 2


def test_primary_lab_rows_are_not_promoted_without_depth():
    data = rows("si_primary_lab_gf_census.csv")
    assert len(data) == 22
    assert {r["grade"] for r in data} == {"ungraded"}
    assert {r["depth_status"] for r in data} == {"NOT_MEASURED"}


def test_freeze_gate_is_honest():
    summary = json.loads((OUT / "summary.json").read_text())
    assert summary["freeze_status"] == "BLOCKED"
    assert summary["active_graded_rows"] == 0


def test_si2_keeps_its_own_experimental_provenance():
    """Scott 2015 section 5.4 and Table 2 note 8, distinct from note 7."""
    row = next(r for r in rows("si_agss21_reference_lines.csv") if r["species"] == "Si II")
    assert float(row["published_loggf"]) == -0.044
    assert "Garz" not in row["published_gf_source"]
    for source in ("Schulz-Gulde1969", "Blanco1995", "Matheron2001"):
        assert source in row["published_gf_source"]
    assert float(row["published_gf_sigma_dex"]) == 0.02
    assert all(not r["published_gf_sigma_dex"] for r in rows("si_agss21_reference_lines.csv")
               if r["species"] == "Si I")
