import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LINESET = ROOT / "data/linelists/reference_sets/agss21_solar_al_rya1134.csv"
PROVENANCE = ROOT / "data/reference/asplund2021_al/provenance_rya1134.json"


def rows():
    with LINESET.open(newline="") as stream:
        return list(csv.DictReader(stream))


def test_lineage_has_six_adopted_and_one_explicit_telluric_rejection():
    data = rows()
    assert len(data) == 7
    assert sum(row["nordlander2017_solar_selected"] == "Y" for row in data) == 6
    rejected = [row for row in data if row["nordlander2017_solar_selected"] == "N"]
    assert len(rejected) == 1
    assert rejected[0]["canonical_line_id"] == "alphys_I_10891.7360_0405"
    assert rejected[0]["agss21_adoption_status"] == "EXCLUDED_TELLURIC"


def test_every_line_is_physically_crossmatched_to_frozen_intake():
    data = rows()
    assert all(row["canonical_line_id"].startswith("alphys_I_") for row in data)
    assert all(float(row["match_distance_mA"]) <= 10.0 for row in data)
    assert all(row["manifest_intake_status"] in {"FROZEN", "CROSSMATCH_REVIEW"} for row in data)


def test_agss21_source_count_discrepancy_is_not_erased():
    provenance = json.loads(PROVENANCE.read_text())
    assert provenance["counts"] == {
        "scott2015_retained": 7,
        "nordlander2017_after_telluric_exclusion": 6,
    }
    discrepancy = provenance["published_count_discrepancy"]
    assert discrepancy["agss21_prose_count"] == 5
    assert discrepancy["traceable_selection_count"] == 6
    assert discrepancy["status"] == "OPEN_DOCUMENTED_SOURCE_INCONSISTENCY"
