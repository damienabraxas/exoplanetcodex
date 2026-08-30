import json

import pandas as pd

from scripts.build_atomic_intake_rya1129 import build


def test_builds_cno_atomic_manifests_without_abundances(tmp_path):
    build(tmp_path)
    metadata = json.loads((tmp_path / "manifest_metadata.json").read_text())
    assert metadata["abundances_generated"] is False
    assert metadata["registry_n_rows_observed"] == metadata["registry_n_targets_declared"]
    assert metadata["registry_n_unique_atomic_symbols"] == 26

    for element in ("C", "N", "O"):
        manifest = pd.read_csv(tmp_path / f"{element}_atomic_manifest.csv")
        assert len(manifest) > 0
        assert manifest.species.str.startswith(element + " ").all()
        assert set(manifest.medium) == {"air"}
        assert not manifest.canonical_line_id.duplicated().any()
        assert set(manifest.intake_status) <= {
            "INVENTORIED_NOT_FROZEN", "CROSSMATCH_REVIEW", "SOURCE_REVIEW"
        }


def test_ledger_is_live_registry_driven(tmp_path):
    build(tmp_path)
    ledger = pd.read_csv(tmp_path / "intake_status_ledger.csv")
    assert len(ledger) == 26
    assert "Zn" not in set(ledger.element)
    assert set(ledger[ledger.element.isin(["C", "N", "O"])].intake_status) == {"CROSSMATCH_REVIEW"}
