import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/literature/igrins_nandakumar_2024"


def rows(name):
    with (DATA / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def test_complete_cds_package_is_preserved_with_hashes():
    provenance = json.loads((DATA / "provenance.json").read_text())
    assert provenance["raw_files_immutable"] is True
    assert len(provenance["files"]) == 23
    assert all(len(meta["sha256"]) == 64 for meta in provenance["files"].values())


def test_all_published_lines_and_observations_are_normalized():
    lines = rows("igrins_nandakumar_2024_lines.csv")
    observations = rows("igrins_nandakumar_2024_observations.csv")
    assert len(lines) == 76
    assert len(observations) == 76 * 50
    assert {row["element"] for row in lines} == {"F", "Na", "Mg", "Al", "Si", "S", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Co", "Ni", "Cu", "Zn", "Y", "Ce", "Nd", "Yb"}
    assert len([row for row in lines if row["element"] == "Al"]) == 7


def test_wavelength_only_matches_are_never_promoted():
    lines = rows("igrins_nandakumar_2024_lines.csv")
    assert all(row["promotion_allowed"] == "False" for row in lines)
    assert all(row["crossmatch_status"] in {"HOLD_TRANSITION_IDENTITY_INCOMPLETE", "NO_CANONICAL_WAVELENGTH_CANDIDATE"} for row in lines)
    for row in lines:
        assert row["excitation_potential_eV"] == ""
        assert row["log_gf_igrins"] == ""


def test_delta_matrix_refuses_unearned_gradeability():
    delta = rows("igrins_element_delta.csv")
    assert len(delta) == 21
    assert sum(int(row["igrins_lines"]) for row in delta) == 76
    assert all(row["transition_confirmed_matches"] == "0" for row in delta)
    assert all(row["primary_lab_gradeable"] == "0" for row in delta)


def test_al_audit_is_complete_and_non_promotional():
    audit = rows("igrins_al_completeness_audit.csv")
    assert len(audit) == 7
    assert all(row["promotion_allowed"] == "False" for row in audit)
