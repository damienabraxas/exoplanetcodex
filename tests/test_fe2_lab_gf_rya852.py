"""
tests/test_fe2_lab_gf_rya852.py — RYA-852
=========================================
The Fe II arbiter lines cannot be graded, and the interesting part is why: the grade already
stored for two of them says `B` while NIST — the source the label names — says `D` and `E`.
Grade B is 0.041 dex; D and E are 0.176 and 0.301. RYA-850 keys its graded gf term on that
metadata, so a wrong grade here does not stay here.

Most of what follows pins NEGATIVE results and traps, because both are what this ticket
actually produced and both are what gets quietly re-discovered:

  * a vacuum/air default that hides every line in the band, and
  * a cross-match that needs EP on BOTH sides, which I got wrong on the second side after
    flagging it on the first.

The NIST query needs network + astroquery and runs on Sirius; what runs here is the
arithmetic, the recorded findings, and the guards.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

RES = ROOT / "data" / "results" / "rya852"
SUMMARY = RES / "rya852_summary.json"
AUDIT = RES / "rya852_fe2_gf_audit.csv"


@pytest.fixture(scope="module")
def summary():
    if not SUMMARY.exists():
        pytest.skip("RYA-852 summary absent (Sirius artifact)")
    return json.loads(SUMMARY.read_text())


# ── the accuracy arithmetic ───────────────────────────────────────────────────

def test_nist_accuracy_letters_convert_to_the_dex_they_claim():
    """B = 10% = log10(1.10) = 0.041 — the same number error_budget uses for a graded gf.
    D = 50% and E = 100% are 4-7x larger, which is the whole finding."""
    import rya852_fe2_lab_gf_audit as m
    assert m.accuracy_to_dex("B") == pytest.approx(0.0414, abs=5e-4)
    assert m.accuracy_to_dex("D") == pytest.approx(0.1761, abs=5e-4)
    assert m.accuracy_to_dex("E") == pytest.approx(0.3010, abs=5e-4)
    assert m.accuracy_to_dex("B") == pytest.approx(m.GRADED_GF_TERM_DEX, abs=1e-3)


# ── the defect ────────────────────────────────────────────────────────────────

def test_the_stored_grade_disagrees_with_nist(summary):
    """🔴 THE FINDING. Two lines are stored as NIST grade B; NIST grades them D and E, and
    their stored log gf is not NIST's either."""
    bad = summary["stored_grade_disagrees_with_nist"]
    assert len(bad) == 2
    for b in bad:
        assert b["stored"] == "B"
        assert b["nist_accuracy"] in ("D", "E")
        assert b["nist_dex"] > 4 * b["stored_dex"], (
            "the stored grade no longer understates the cited accuracy — re-derive")
        assert abs(b["delta_loggf"]) > 0.1, (
            "the stored value now agrees with NIST — the 'NIST ASD' label may have been "
            "corrected, so re-check the grade claim too")


def test_no_arbiter_line_is_graded(summary):
    """Do not fabricate coverage. Every arbiter line's cited accuracy is worse than the
    ~0.1 dex floor the ticket anticipated, so none is entitled to a graded term."""
    acc = summary["arbiter_nist_accuracy_dex"]
    assert len(acc) == 3
    assert min(acc.values()) > 0.1, (
        "an arbiter line now has a cited accuracy better than 0.1 dex — grading may be "
        "possible and RYA-852's verdict needs revisiting")
    assert "UNGRADED" in summary["verdict"]


def test_the_arbiter_trio_is_read_from_the_products_not_the_ticket(summary):
    """The ticket guessed 6247.6 / 6432.7 / 6456.4; only the first is in the trio."""
    arb = set(round(w, 3) for w in summary["arbiter_lines_air_A"])
    assert arb == {6147.734, 6238.386, 6247.557}
    guess = set(summary["ticket_guess_was"])
    assert len(arb & {6247.557}) == 1
    assert 6432.7 in guess and 6456.4 in guess, "the guess is recorded for contrast"


# ── the pool-wide scale offset ────────────────────────────────────────────────

def test_the_pool_sits_above_nist_as_a_whole(summary):
    """Not just the mislabelled lines — a coherent offset across the pool, its plain-VALD3
    members included, which is a SCALE difference rather than independent errors."""
    off = summary["pool_scale_offset"]
    assert off["median_dex"] > 0.05
    assert off["n_above_nist"] >= 0.7 * off["n_matched"]
    assert "407" in off["note"], "the ionization-balance consequence must travel with it"


def test_the_offset_is_larger_than_the_balance_it_underwrites():
    """Fe I - Fe II is +0.018 on the current scale. The gf offset is ~0.106, so the balance
    is smaller than the systematic that sets it — which is why this is a finding and not a
    footnote."""
    fe1, fe2, offset = 7.586, 7.568, 0.106
    assert abs(fe1 - fe2) < offset


# ── the traps, pinned because both nearly produced a wrong answer ─────────────

def test_the_vacuum_air_trap_is_recorded(summary):
    """astroquery.nist defaults to VACUUM; at 6150 A that is +1.71 A and every arbiter line
    vanishes from the result. Same class as RYA-846's Wallace atlas."""
    traps = " ".join(summary["traps"]).lower()
    assert "vacuum" in traps and "air" in traps
    src = (ROOT / "scripts" / "rya852_fe2_lab_gf_audit.py").read_text()
    assert 'wavelength_type="vac+air"' in src, (
        "the NIST query lost its explicit wavelength_type — the default is vacuum and it "
        "silently hides every line in this band")


def test_both_sides_of_the_cross_match_use_ep():
    """I flagged the wavelength-only trap for canonical_gf and then made the same mistake on
    the NIST side, which produced a -2.298 dex offset for 6432.676. Both sides need EP."""
    src = (ROOT / "scripts" / "rya852_fe2_lab_gf_audit.py").read_text()
    assert "nist_ep_eV" in src, "the NIST side no longer carries an EP to match on"
    assert src.count("EP_TOL_EV") >= 3, (
        "an EP tolerance disappeared from one side of the cross-match")


def test_the_unverified_claims_are_labelled_unverified(summary):
    """RYA-833: an absence in the SEARCH is not an absence in the SOURCE. Den Hartog 2019
    could not be reached, and that must not harden into 'DH19 does not cover these lines'."""
    ne = " ".join(summary["not_established"]).lower()
    assert "den hartog" in ne and ("unverified" in ne or "not the source" in ne)
    assert "s/l" in ne or "firewall" in ne


def test_the_audit_table_keeps_every_pool_line(summary):
    """RYA-711: the lines that cannot be graded are carried with their evidence."""
    if not AUDIT.exists():
        pytest.skip("audit table absent")
    d = pd.read_csv(AUDIT)
    assert len(d) == summary["n_pool_lines"] == 11
    assert d.is_arbiter.sum() == 3
