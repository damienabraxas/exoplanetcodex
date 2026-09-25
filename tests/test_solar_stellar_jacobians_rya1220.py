"""RYA-1220 -- solar stellar-parameter Jacobians for the CNO diagnostics.

Closes stellar.teff and stellar.xi, which the CNO gate matrix carried as HOLD on every
route. The discipline worth pinning is the distinction between a SLOPE and CURVATURE:
a leg pair that straddles nominal gives a slope, one that does not does not.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "rya1220_assemble_stellar_jacobians.py"
ART = ROOT / "data" / "output" / "rya1220" / "solar_stellar_jacobians.json"


def _mod():
    spec = importlib.util.spec_from_file_location("jac", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_straddling_legs_give_a_slope():
    m = _mod()
    fits = {"nominal": {"A_X": 8.0, "constrained": True},
            "teff_K_-1": {"A_X": 7.98, "constrained": True},
            "teff_K_+1": {"A_X": 8.02, "constrained": True}}
    got = m.jacobian(fits, "teff_K_-1", "teff_K_+1")
    assert got["state"] == "MEASURED"
    assert got["sigma_dex"] == pytest.approx(0.02, abs=1e-6)


def test_same_side_legs_are_unresolved_not_a_small_slope():
    """Both legs above nominal: the central difference is curvature, not a slope."""
    m = _mod()
    fits = {"nominal": {"A_X": 8.477, "constrained": True},
            "teff_K_-1": {"A_X": 8.485, "constrained": True},
            "teff_K_+1": {"A_X": 8.483, "constrained": True}}
    got = m.jacobian(fits, "teff_K_-1", "teff_K_+1")
    assert got["state"] == "UNRESOLVED"
    assert got["sigma_dex"] == pytest.approx(0.008, abs=1e-6)   # the larger excursion
    assert "curvature, not a slope" in got["why"]


def test_a_leg_sitting_exactly_on_nominal_does_not_count_as_straddling():
    """N's teff_+1 landed exactly on nominal; there is no resolved upward response."""
    m = _mod()
    fits = {"nominal": {"A_X": 7.942, "constrained": True},
            "teff_K_-1": {"A_X": 7.935, "constrained": True},
            "teff_K_+1": {"A_X": 7.942, "constrained": True}}
    got = m.jacobian(fits, "teff_K_-1", "teff_K_+1")
    assert got["state"] == "UNRESOLVED"


def test_a_nonconverged_leg_has_no_slope():
    """CN_red's xi leg returned NON-MINIMUM; a fit that determined nothing has no slope."""
    m = _mod()
    fits = {"nominal": {"A_X": 7.5, "constrained": True},
            "vturb_kms_-1": {"A_X": 7.4, "constrained": True},
            "vturb_kms_+1": {"A_X": 7.126, "constrained": False}}
    got = m.jacobian(fits, "vturb_kms_-1", "vturb_kms_+1")
    assert got["state"] == "HOLD"
    assert got["sigma_dex"] is None
    assert "non-minimum" in got["why"].lower()


def test_nitrogen_is_measured_on_the_ratified_route():
    m = _mod()
    assert m.RUNS["N"][1] == "nir_cn_iag"
    assert m.RUNS["N"][2] == "CN_AX_IR"


@pytest.mark.skipif(not ART.exists(), reason="jacobians not assembled")
def test_every_cno_element_has_both_parameter_components():
    doc = json.loads(ART.read_text())
    by_el = {r["element"]: r for r in doc["elements"]}
    assert set(by_el) == {"C", "N", "O"}
    for element, row in by_el.items():
        assert row.get("state") != "NOT_RUN", f"{element} was not run"
        for comp in ("stellar.teff", "stellar.xi"):
            g = row["components"][comp]
            assert g["state"] in {"MEASURED", "UNRESOLVED"}, f"{element} {comp}"
            assert g["sigma_dex"] is not None, f"{element} {comp} has no value or bound"


@pytest.mark.skipif(not ART.exists(), reason="jacobians not assembled")
def test_the_parameter_terms_are_small_against_the_dominant_systematics():
    """The point of the result: stellar parameters are a MINOR term across CNO.

    Every value or bound is <= 0.016 dex, against a molecular model-form term of
    0.061-0.140 and a 0.225 dex fit start-dependence. Closing these two components
    barely moves the budget, and that is worth stating rather than implying.
    """
    doc = json.loads(ART.read_text())
    for row in doc["elements"]:
        for comp, g in row["components"].items():
            assert g["sigma_dex"] <= 0.02, f"{row['element']} {comp} = {g['sigma_dex']}"
