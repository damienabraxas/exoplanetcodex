"""RYA-785 — the Fe Gerber deck's gate, pinned where it is easiest to get wrong.

These do not re-run Turbospectrum. They pin the two things the second pass changed and the
one thing it deliberately did NOT change, so a later edit cannot quietly undo the
adjudication or manufacture a pass by loosening the test.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
GATE = ROOT / "scripts" / "ts_gerber_gate.py"
PROV = ROOT / "data" / "nlte_grids" / "gerber_ts" / "Fe_gerber2023.prov.json"


def _cfg():
    """Load the gate module's Fe config without importing its Sirius-only dependencies."""
    src = GATE.read_text()
    spec = importlib.util.spec_from_file_location("_g", GATE)
    m = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(m)
    except Exception:
        pytest.skip("gate module needs the Sirius synthesis stack")
    return m


def test_fe_anchor_is_the_decks_own_published_value_not_mpia():
    """THE ADJUDICATION. RYA-785 says validate against the deck's OWN published anchor and
    forbids MPIA as the target; RYA-525/712 make a Gerber-vs-MPIA difference a DIAGNOSTIC.

    The first pass anchored on MPIA (+0.0124) and returned CHECK by 0.001 dex. The anchor is
    now Gerber+2023 A&A 669 A43 Table 3/4 (solar 1D NLTE - 1D LTE = +0.06). If this ever
    reverts toward the MPIA value, the gate is asking the wrong question again.
    """
    m = _cfg()
    fe = m.ELEMENTS["Fe"]
    assert fe["anchor"] == pytest.approx(0.06), "anchor must be Gerber+2023's published +0.06"
    assert "Gerber+2023" in fe["ref"] and "not MPIA" in fe["ref"]
    assert fe["anchor"] > 0.04, "an anchor near +0.012 is the MPIA value — wrong referee"


def test_tolerance_was_not_loosened_to_obtain_the_pass():
    """Validate-don't-tune (RYA-161). Only the anchor and the line set changed; the
    tolerance is the same 0.05 the first pass failed against."""
    assert _cfg().ELEMENTS["Fe"]["tol"] == pytest.approx(0.05)


def test_fe_gate_lines_are_isolated_and_plural():
    """The first pass used 3 lines, 2 of which had an Fe I neighbour inside the 1.2 A EW
    window (5806.725 at 0.098 A, 6027.070 at 0.239 A). A same-species blend does not cancel
    in the EW->A* inversion. Both are gone; the set is wider."""
    waves = _cfg().ELEMENTS["Fe"]["waves"]
    assert len(waves) >= 8, "3 lines is not a determination"
    assert 5806.725 not in waves and 6027.070 not in waves, "contaminated lines are back"
    assert 5905.671 in waves, "the one clean first-pass line should be retained"


def test_weak_line_guard_applies_to_the_synthetic_ew():
    """The Ti lesson enforced on the quantity the inversion actually consumes.

    Selection filtered on MEASURED EW, but A* is recovered from the SYNTHETIC EW, and Fe I
    5607.664 differs between them by 21x (20.1 mA measured, 0.95 mA synthesised).
    """
    m = _cfg()
    assert getattr(m, "SYN_EW_MIN", 0) >= 5.0
    assert m.SYN_EW_MIN < 20.0, "a floor this high would start discarding real gate lines"


def test_prov_records_pass_and_keeps_the_cross_engine_delta_as_a_diagnostic():
    d = json.loads(PROV.read_text())
    g = d["gate"]
    assert g["verdict"].startswith("PASS")
    assert g["tol"] == 0.05
    # the ~5x Gerber-vs-MPIA difference is a FINDING, never folded into the error bar
    assert "cross_engine_diagnostic" in g
    assert "never folded into the error bar" in g["cross_engine_diagnostic"]
    # and the deck's abundance is the deck's property, not a fitted choice
    assert d["deck_abundance"]["a_sun"] == pytest.approx(7.46)


def test_gerber_nlte_is_wired_and_cannot_silently_fall_back_to_lte():
    """RYA-798 SUPERSEDES the refusal this test used to pin.

    RYA-785 left the driver refusing `gerber-nlte` because the deck was validated but had
    no route to a synthesis. RYA-798 built that route, so "it must refuse" is no longer the
    invariant — but the reason the refusal existed still is: **an LTE spectrum must never
    ship under an NLTE label**. iSpec drops an unlabelled element into `nlte_ignored` and
    synthesises LTE without raising (RYA-764), and it never reports `nlte_available` back,
    so the check has to happen before the call. That guard is what this now pins.
    """
    src = (ROOT / "scripts" / "derive_band_products.py").read_text()
    assert "gerber-nlte" in src
    # the two stale refusal reasons must both be gone
    assert "is not provisioned: the TS-native Gerber" not in src
    assert "VALIDATED but NOT WIRED" not in src
    # the product is SEPARATE, never a correction to Engine-B-LTE (RYA-712)
    assert "ENGINE-B-NLTE" in src
    # and the deck must run on the atmosphere family it was computed on
    assert "MARCS.GES" in src, "the Gerber deck is MARCS; ATLAS9 has no usable tau column"

    gn = (ROOT / "pipeline" / "gerber_nlte.py").read_text()
    assert "def assert_linelist_supports_nlte" in gn, "the RYA-764 silent-LTE guard"
    assert "def assert_depth_match" in gn, "iSpec overwrites the departure tau"
