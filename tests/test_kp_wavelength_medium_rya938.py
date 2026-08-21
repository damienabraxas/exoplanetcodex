"""RYA-938 — Kurucz 2005 is on a VACUUM grid, and that is what RYA-929 hit.

RYA-929's window sweep was a clean GO and its line sweep was nonsense: Kurucz
showed near-zero depth at Al I 6696 and K I 7665/7699 where the 1984 atlas and
the telluric-free IAG atlas agree on deep lines. Both results are consistent
with one cause. Reading a vacuum grid as air displaces every line by ~1.7 A at
6600 A -- roughly 200 sampled pixels -- which empties a +/-0.35 A line window
while barely moving a 17-300 A window statistic.

The invariant pinned here is the RELATIONSHIP, not the numbers: the measured
registration lag must equal the air-vacuum offset predicted independently by
`pipeline.wavelength_util`, and converting must drive it to zero. That keeps
holding if the atlas is re-staged or the windows change.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pipeline.wavelength_util import air_to_vac, vac_to_air

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "data" / "audit" / "rya938_kp_crosscheck"


def _report() -> dict:
    return json.loads((EVIDENCE / "crosscheck.json").read_text())


def test_air_vacuum_round_trip_is_the_reference_the_test_leans_on():
    """Positive control on the instrument itself before trusting its verdict."""
    air = np.array([5180.0, 6130.0, 6625.0, 8216.0])
    assert np.allclose(vac_to_air(air_to_vac(air)), air, atol=1e-6)
    offset = air_to_vac(air) - air
    assert np.all(offset > 0), "vacuum wavelengths exceed air wavelengths"
    assert np.all(np.diff(offset) > 0), "the offset grows with wavelength"


def test_the_measured_lag_equals_the_predicted_air_vacuum_offset():
    """The relationship, not the constant: lag == offset, window by window."""
    for row in _report()["registration"]:
        measured = row["kurucz_as_read_nm_lag_A"]
        predicted = row["air_vac_offset_at_mid_A"]
        assert measured == pytest.approx(predicted, abs=0.01), row["window_A"]


def test_converting_to_air_drives_the_lag_to_zero():
    for row in _report()["registration"]:
        assert row["kurucz_vac_to_air_lag_A"] == pytest.approx(0.0, abs=0.01)
        assert row["kurucz_vac_to_air_corr"] > 0.99


def test_the_iag_control_proves_the_test_discriminates():
    """An already-air product must return zero, or a zero means nothing."""
    for row in _report()["registration"]:
        assert row["iag_vs_kp1984_lag_A"] == pytest.approx(0.0, abs=0.01)


def test_the_verdict_is_recorded_with_a_threshold_declared_before_measuring():
    verdict = _report()["medium_verdict"]
    assert verdict["kurucz2005"] == "vacuum"
    assert verdict["decision_threshold_A"] == pytest.approx(
        0.5 * verdict["median_air_vac_offset_A"], rel=1e-9)
    assert verdict["median_lag_as_read_A"] > verdict["decision_threshold_A"]


def test_the_correction_recovers_the_lines_rya929_could_not_explain():
    """Every flagged line must land on the independent IAG depth once converted."""
    flagged = {"Al I 6696", "K I 7699", "N I 8216", "Fe I 6881"}
    seen = set()
    for line in _report()["diagnostic_lines"]:
        if line["line"] not in flagged:
            continue
        seen.add(line["line"])
        iag = line["iag"]["depth"]
        before = abs(line["kurucz_as_read_nm"]["depth"] - iag)
        after = abs(line["kurucz_vac_to_air"]["depth"] - iag)
        assert after < before, f'{line["line"]}: conversion did not improve agreement'
        assert after < 0.05, f'{line["line"]}: still disagrees with IAG by {after:.4f}'
    assert seen == flagged


def test_kurucz_does_not_cover_the_reddest_registered_telluric_band():
    """Kurucz 2005 stops at 1000 nm, so it cannot replace 1984 beyond it."""
    bands = {(b["lo_A"], b["hi_A"]): b for b in _report()["telluric_bands"]}
    reddest = bands[(11120.0, 11560.0)]
    assert reddest["kurucz2005"]["min"] is None
    assert reddest["iag"]["min"] is None
    assert reddest["kp1984"]["pct_below_0.5"] > 20.0, (
        "the 1984 atlas is the only arm reaching this band, and it is uncorrected")


def test_1984_is_uncorrected_and_2005_tracks_the_telluric_free_reference():
    """The window-level GO from RYA-929 survives, and is restated as a relationship."""
    for band in _report()["telluric_bands"]:
        if band["band"] == "clean control" or band["kurucz2005"]["min"] is None:
            continue
        kp = band["kp1984"]["pct_below_0.5"]
        kur = band["kurucz2005"]["pct_below_0.5"]
        iag = band["iag"]["pct_below_0.5"]
        assert kp > kur, f'{band["band"]}: 1984 must retain more absorption'
        assert abs(kur - iag) < 1.0, (
            f'{band["band"]}: corrected products must agree with each other')


def test_clean_controls_agree_across_all_three_products():
    """If the controls disagreed, the band comparison would mean nothing."""
    for band in _report()["telluric_bands"]:
        if band["band"] != "clean control":
            continue
        values = [band[p]["pct_below_0.5"] for p in ("kp1984", "kurucz2005", "iag")]
        assert max(values) - min(values) < 1.0, band["lo_A"]
