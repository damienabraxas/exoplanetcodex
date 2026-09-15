"""Regression tests for the RYA-1155 coverage plumbing and band gaps."""
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def test_registered_normalized_crires_products_resolve_in_coverage_registry():
    from pipeline.coverage import load_registry

    rows = load_registry("solar")
    ids = {r.instrument_id for r in rows if r.instrument_id == "crires_plus"}
    assert "crires_plus" in ids
    spans = {(round(r.wave_min_A, 2), round(r.wave_max_A, 2))
             for r in rows if r.instrument_id == "crires_plus"}
    assert (15007.11, 17493.69) in spans
    assert any(r.loader == "csv_normalized" and r.reachable
               for r in rows if r.instrument_id == "crires_plus")


def test_raw_crires_inventory_is_still_not_mistaken_for_a_spectrum():
    from pipeline.coverage import load_registry

    rows = load_registry("solar")
    assert all(r.path != str(ROOT / "data/audit/vesta_crires_plus/vesta_crires_plus_idp_manifest.csv")
               for r in rows)


def test_al_band_classifier_covers_the_census_nir_gaps():
    from scripts.build_al_intake_rya1132 import band

    assert band(13123.416) == "NIR"
    assert band(13150.753) == "NIR"
    assert band(17699.094) == "NIR"
