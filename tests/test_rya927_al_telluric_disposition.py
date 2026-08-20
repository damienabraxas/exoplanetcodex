"""RYA-927: the HARPS product gap is not an Al-line contamination claim."""

from pathlib import Path

import pandas as pd

from pipeline.telluric_policy import TELLURIC_BANDS, exclusion


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "data/audit/rya927_harps_telluric/al_line_telluric_disposition.csv"
WEBSITE = ROOT / "data/measured/sol_ew_results_v1.csv"


def test_every_website_al_ew_is_frozen_in_clean_line_audit():
    website = pd.read_csv(WEBSITE)
    website = website.loc[website.element.eq("Al"),
                          ["wavelength_air_A", "ew_mA", "ew_err_mA"]]
    audit = pd.read_csv(AUDIT)
    assert len(website) == len(audit) == 2
    joined = website.merge(audit, on="wavelength_air_A", validate="one_to_one")
    assert (joined.ew_mA == joined.website_ew_mA).all()
    assert (joined.ew_err_mA == joined.website_ew_err_mA).all()


def test_website_al_lines_are_outside_every_registered_telluric_band():
    audit = pd.read_csv(AUDIT)
    for wave in audit.wavelength_air_A:
        assert exclusion(wave, instrument="harps") == ""
        assert all(not (lo <= wave <= hi) for lo, hi, _ in TELLURIC_BANDS)
    assert set(audit.disposition) == {"CLEAN_OUTSIDE_REGISTERED_BANDS"}


def test_distance_to_nearest_band_is_recomputed_not_asserted_by_label():
    audit = pd.read_csv(AUDIT)
    first_blue_edge = min(lo for lo, _, _ in TELLURIC_BANDS)
    expected = first_blue_edge - audit.wavelength_air_A
    assert (abs(expected - audit.distance_to_band_A) < 1e-9).all()
