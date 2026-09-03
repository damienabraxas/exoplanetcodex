import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/results/rya1169"


def rows(name):
    with (OUT / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def setup_module():
    subprocess.run([sys.executable, str(ROOT / "scripts/build_si_intake_rya1169.py")], check=True)
    subprocess.run([sys.executable, str(ROOT / "scripts/grade_si_asplund_rya1169.py")], check=True)


def test_source_membership_is_not_conflated_with_depth():
    data = rows("si_asplund_line_test.csv")
    assert sum(x["asplund_grade"] == "asplund" for x in data) == 10
    assert sum(x["asplund_grade"] == "excluded" for x in data) == 2
    assert {x["depth_route"] for x in data if x["asplund_grade"] == "asplund"} == {"ew_valid"}
    assert {x["depth_route"] for x in data if x["asplund_grade"] == "excluded"} == {"deep_synthesis"}


def test_all_bands_are_reported_and_empty_cells_not_backfilled():
    data = rows("si_asplund_all_band_matrix.csv")
    assert [x["band"] for x in data] == ["FUV", "NUV", "VIS", "red-optical", "NIR", "H", "K"]
    active = {x["band"]: int(x["asplund_active"]) for x in data}
    assert active == {"FUV": 0, "NUV": 0, "VIS": 8, "red-optical": 2, "NIR": 0, "H": 0, "K": 0}
    assert all(x["source_scope_note"] for x in data if x["band"] in {"NIR", "H", "K"})


def test_all_ten_source_used_lines_pass_depth_test():
    summary = json.loads((OUT / "si_asplund_test_summary.json").read_text())
    assert summary["source_used"] == 10
    assert summary["active_test_passed"] == 10
    assert summary["ew_valid"] == 10
    assert summary["deep_synthesis"] == 0
    assert summary["all_band_verdict"] == "PASS_WITH_SOURCE_DEFINED_GAPS"
