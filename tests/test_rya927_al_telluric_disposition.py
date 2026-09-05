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
    """⚠️ THE FORMULA ASSUMED THE BLUEST BAND WAS THE NEAREST ONE, AND RYA-1193 BROKE THAT.
    It computed `min(lo) - wavelength`, which is the distance to the bluest band's blue
    edge -- correct only while every audited line sits blueward of every band. RYA-1193
    added the O2 gamma-band at 6270 A, blueward of both Al lines (6631.218, 6696.185), so
    the old expression went NEGATIVE (-361, -426) against a committed +235.782 / +170.815.

    🔴 The DATA was never wrong: the nearest band to both lines is still O2 B at 6867 A.
    Computing the true nearest EDGE reproduces the committed values exactly, and is
    robust to a band being added on either side."""
    audit = pd.read_csv(AUDIT)

    def nearest_edge(w: float) -> float:
        return min(0.0 if lo <= w <= hi else min(abs(lo - w), abs(hi - w))
                   for lo, hi, _ in TELLURIC_BANDS)

    expected = audit.wavelength_air_A.map(nearest_edge)
    assert (abs(expected - audit.distance_to_band_A) < 1e-9).all(), (
        f"expected {list(expected)} vs committed {list(audit.distance_to_band_A)}")
