import json
from pathlib import Path

import pytest

from pipeline import molecular_canonical as mc

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "data/schemas/canonical_molecular_transition.schema.json"


def _fixture(tmp_path: Path, name="CO_IR_Li2015.dat") -> Path:
    path = tmp_path / name
    path.write_text("' code' 1 1\n'ExoMol Li2015'\n"
                    " 4515.025 0.01716 -16.314 0.000 19.0 8.35E-10 'X' 'X' 9.0 8.0 'X' 'X' 0.0 0.0 'v11-0_J9-8_Li2015'\n")
    return path


def _complete_metadata():
    return {
        "database": "ExoMol", "release": "Li2015", "retrieved_at": "2026-08-30",
        "license": "CC BY 4.0", "partition_function_source": "source:pf",
        "dissociation_energy_source": "source:de", "component_convention": "one row per component",
        "isotopologue_convention": "explicit isotopologue", "wavelength_medium": "vacuum",
        "wavelength_frame": "laboratory", "transition_probability_source": "Li2015",
        "probability_kind": "ab_initio", "conversion_provenance": "test conversion",
    }


def test_schema_is_separate_from_atomic_gf_and_has_required_domains():
    schema = json.loads(SCHEMA.read_text())
    text = SCHEMA.read_text()
    assert schema["properties"]["schema"]["const"] == mc.SCHEMA
    assert "canonical_gf" in schema["description"]
    assert "Codex Grade" not in text and "Deep Grade" not in text
    for field in ("isotopologue", "transition_identity", "source", "thermochemistry",
                  "normalization", "precedence", "band_reach", "intake_status"):
        assert field in schema["properties"]
        assert field in schema["required"]


def test_complete_provenance_yields_source_verified_transition(tmp_path):
    row = next(mc.canonical_rows(_fixture(tmp_path), "CO", _complete_metadata()))
    assert row["intake_status"] == "SOURCE_VERIFIED"
    assert row["isotopologue"] == "12C16O"
    assert row["transition_identity"]["v_upper"] == "11"
    assert row["transition_identity"]["v_lower"] == "0"
    assert row["transition_identity"]["j_upper"] == "9"
    assert row["transition_identity"]["j_lower"] == "8"
    assert row["source"]["source_sha256"] == mc.sha256(_fixture(tmp_path))


def test_missing_typed_provenance_is_loudly_blocked(tmp_path):
    row = next(mc.canonical_rows(_fixture(tmp_path, "12CH_400-450.bsyn"), "CH", {}))
    assert row["intake_status"] == "BLOCKED_PROVENANCE"
    assert row["transition_probability"]["source_id"] == "UNRESOLVED"
    assert row["source"]["license"] is None


def test_raw_transition_identity_is_never_discarded(tmp_path):
    path = _fixture(tmp_path, "12CH_400-450.bsyn")
    row = next(mc.canonical_rows(path, "CH", _complete_metadata()))
    assert row["transition_identity"]["raw_label"] == "v11-0_J9-8_Li2015"


def test_malformed_rows_fail(tmp_path):
    path = tmp_path / "12CH_400-450.bsyn"
    path.write_text("' h'\n's'\n4200.0 0.1 -1.0\n")
    with pytest.raises(ValueError, match="malformed"):
        list(mc.parse_bsyn_rows(path))


def test_source_omitted_identity_is_preserved_as_review_not_invented(tmp_path):
    path = tmp_path / "12CH_400-450.bsyn"
    path.write_text("' h'\n's'\n"
                    "4200.0 0.1 -1.0 0 1 1 'X' 'X' 0 0 'X' 'X' 0 0 ''\n")
    row = next(mc.canonical_rows(path, "CH", _complete_metadata()))
    assert row["transition_identity"]["raw_label"] == "SOURCE_LABEL_NOT_SUPPLIED"
    assert row["intake_status"] == "CROSSMATCH_REVIEW"


def test_band_reporting_reaches_uv_through_ir():
    assert mc.bands_for_wavelength(1500) == ["FUV"]
    assert mc.bands_for_wavelength(3000) == ["NUV"]
    assert mc.bands_for_wavelength(5000) == ["VIS"]
    assert mc.bands_for_wavelength(8000) == ["RED_OPTICAL"]
    assert mc.bands_for_wavelength(16000) == ["NIR"]
    assert mc.bands_for_wavelength(46000) == ["IR"]


@pytest.mark.parametrize(("filename", "expected"), [
    ("12C-1H__MoLLIST_rovib.bsyn", "12C1H"),
    ("14N-1H__kNigHt_rovib.bsyn", "14N1H"),
    ("16O-1H__MYTHOS_rovib.bsyn", "16O1H"),
])
def test_explicit_rovibrational_isotopologue_filenames(filename, expected):
    assert mc.isotopologue_from_filename(Path(filename)) == expected
