"""RYA-1193 — the IR science window, pinned where it is read and enforced in code.

The decision (RYA-1094) lived only in a closed ticket, so RYA-1192 re-derived CRIRES+'s
raw 53000 Å reach and audited against it. These tests keep the window on the first-read
surface AND in the config, so the next audit cannot wander out there again.
"""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REGISTER = ROOT / "CODEX_STATE_REGISTER.md"


def test_the_register_pins_the_window_where_it_is_read_first():
    """A decision recorded only in a closed ticket is a decision that gets re-derived."""
    t = REGISTER.read_text()
    head = t[:t.index("## How to read this")]
    assert "IR science window (PINNED, RYA-1094)" in head, (
        "the pin is no longer in the first-read block")
    for fact in ("17277.5", "21386", "24985", "25k", "RYA-1094"):
        assert fact in head, f"the pin no longer states {fact}"
    # and the standing row, in the section that already holds "do NOT re-attempt" rulings
    scope = t[t.index("## Scope / drift corrections"):t.index("## Versioning & edit protocol")]
    assert "IR science window" in scope and "do NOT pursue beyond ~25k" in scope


def test_the_ceiling_is_a_constant_the_code_can_enforce():
    from pipeline.telluric_policy import IR_SCIENCE_CEILING_A, TELLURIC_BANDS
    assert IR_SCIENCE_CEILING_A == 25000.0
    over = [b for b in TELLURIC_BANDS if b[1] > IR_SCIENCE_CEILING_A]
    assert not over, f"bands enumerated above the science ceiling: {over}"


def test_the_h_band_is_declared_so_FALSE_stops_reading_as_clean():
    """🔴 The defect this fixes: the policy stopped at 11560 Å, so `in_telluric_band()`
    answered False for every graded Fe line past 1 µm — and False reads as "clean", not as
    "no band declared here". The graded pool runs to 17277.5 Å."""
    from pipeline.telluric_policy import in_telluric_band, TELLURIC_BANDS
    assert in_telluric_band(15750.0), "CO2 window not declared"
    assert in_telluric_band(16050.0), "CO2 window not declared"
    assert not in_telluric_band(30000.0), "a band exists past the science ceiling"
    assert max(hi for _lo, hi, _n in TELLURIC_BANDS) >= 17400.0


def test_the_o2gamma_band_the_manifests_fitted_is_listed():
    """RYA-940 shipped SEVEN fit manifests; the policy listed six. Edges come from
    o2gamma's own `band_A`, not from a typed-in guess."""
    import json
    from pipeline.telluric_policy import TELLURIC_BANDS
    man = json.loads((ROOT / "data/audit/rya940_kp1984_telluric/o2gamma/"
                      "fit_manifest.json").read_text())
    lo, hi = man["band_A"]
    assert any(abs(a - lo) < 1e-9 and abs(b - hi) < 1e-9 for a, b, _ in TELLURIC_BANDS), (
        f"o2gamma {lo}-{hi} from the manifest is not in TELLURIC_BANDS")
    fitted = len(list((ROOT / "data/audit/rya940_kp1984_telluric").glob("*/fit_manifest.json")))
    listed_below_12k = len([b for b in TELLURIC_BANDS if b[1] <= 12000.0])
    assert listed_below_12k == fitted, (
        f"{fitted} bands were fitted but {listed_below_12k} are listed below 12000 A")


def test_the_gate_ceiling_did_NOT_move_with_the_new_bands():
    """🔴 THE REGRESSION THIS AVOIDS, AND IT WAS MEASURED, NOT FEARED.

    `preflight_readiness.TELLURIC_DOMAIN_MAX_A` used to derive from `max(hi)`. The H-band
    entries are three measured 100 Å windows, not a survey of 11560–17500 Å, so letting the
    ceiling follow them would assert we had looked across six thousand Ångström on the
    strength of three hundred — and, measured, flips an UNCORRECTED 9199–12976 Å band from
    `ir-band-uncorrected` (refused) to satisfied. Adding telluric knowledge would have
    RELAXED the telluric gate."""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    import preflight_readiness as pr
    from pipeline.telluric_policy import ENUMERATION_COMPLETE_TO_A, TELLURIC_BANDS
    assert pr.TELLURIC_DOMAIN_MAX_A == ENUMERATION_COMPLETE_TO_A == 11560.0
    assert pr.TELLURIC_DOMAIN_MAX_A < max(hi for _lo, hi, _n in TELLURIC_BANDS)
    # the gate still refuses an uncorrected IR band
    ok, code, _ = pr.telluric_satisfied("kpno_solar_atlas", "NIR", 9199.0, 12976.0,
                                        "corrected", "not-applied", "audited")
    assert ok is False and code == "ir-band-uncorrected", (
        "an uncorrected 9199-12976 A band is no longer refused — the ceiling drifted")


def test_no_band_is_added_above_the_ceiling():
    """The pin's operative half: enumerating a region we will never measure invites the
    next audit back out to 53000 Å."""
    from pipeline.telluric_policy import TELLURIC_BANDS, IR_SCIENCE_CEILING_A
    assert all(hi <= IR_SCIENCE_CEILING_A for _lo, hi, _n in TELLURIC_BANDS)
    src = (ROOT / "pipeline/telluric_policy.py").read_text()
    assert "DO NOT ADD BANDS ABOVE" in src
    assert "RYA-1094" in src
