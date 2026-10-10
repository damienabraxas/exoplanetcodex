"""RYA-1194 — the telluric "is it corrected?" question is answered PER HOLDING.

RYA-1193 found `basis()` resolving per INSTRUMENT, so `solar_kpno` and
`solar_kpno_molecfit_corrected` received one answer. Harmless only by luck: RYA-1192 had
already proved KP-molecfit is RAW, so `line_selection` happened to be right for both.

It stops being luck the moment RYA-1191 actually corrects that holding — an
instrument-keyed read would still say `line_selection` and the pipeline would not
recognise its own fix. `test_the_rya1191_forward_path_is_recognised` is the test that
matters here: it flips the verified entry and asserts the answer follows.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from pipeline import telluric_policy as TP  # noqa: E402

#: Inside the O2 gamma-band (RYA-1193), so every case below is a line the band set catches.
IN_BAND_A = 6271.28


# ── the spec's own two checks ────────────────────────────────────────────────
def test_two_holdings_on_one_instrument_get_different_bases():
    """Spec item 3. HARPS serves a raw atlas and a molecfit-corrected sibling from one
    instrument row; an instrument-keyed resolver cannot tell them apart and this is the
    case that proves the re-key took."""
    raw, corr = "solar_harps", "solar_harps_molecfit_corrected"
    assert TP.applied_state(raw) == "not-applied"
    assert TP.applied_state(corr) == "applied"
    # one instrument, so the OLD axis gives one answer for both
    assert TP.basis("harps") == TP.basis("harps") == "line_selection"
    # the new axis gives two
    assert TP.holding_basis(raw, "harps") == "line_selection"
    assert TP.holding_basis(corr, "harps") == "corrected"
    assert TP.holding_basis(raw, "harps") != TP.holding_basis(corr, "harps")
    # and it reaches the decision the measurement path actually makes
    assert TP.exclusion(IN_BAND_A, "harps", raw)
    assert not TP.exclusion(IN_BAND_A, "harps", corr)


def test_solar_kpno_still_resolves_line_selection():
    """Spec item 3's control: this is a re-key, not a blanket flip to `corrected`."""
    assert TP.holding_basis("solar_kpno", "kpno_solar_atlas") == "line_selection"
    assert TP.holding_basis("solar_kpno") == "line_selection"       # instrument optional
    assert TP.exclusion(IN_BAND_A, "kpno_solar_atlas", "solar_kpno")


# ── the reason this ticket exists ────────────────────────────────────────────
def test_the_corrected_kitt_peak_holding_is_corrected_and_its_raw_sibling_is_not():
    """RYA-1230 (505bdf0f) corrected solar_kpno_molecfit_corrected over 3000-13000 A, 0 A
    raw, and the map now says so. The raw sibling must not move with it (RYA-1194)."""
    h, i = "solar_kpno_molecfit_corrected", "kpno_solar_atlas"
    assert TP.VERIFIED_HOLDING_STATE[h] == "corrected"
    assert TP.holding_basis(h, i) == "corrected"
    assert not TP.exclusion(IN_BAND_A, i, h)
    assert TP.holding_basis("solar_kpno", i) == "line_selection"
    assert TP.exclusion(IN_BAND_A, i, "solar_kpno")


def test_a_measured_raw_verdict_outranks_the_registry_label(monkeypatch):
    """A MEASUREMENT that a product is raw still beats a registry `applied` claim."""
    h = "solar_kpno_molecfit_corrected"
    assert TP.applied_state(h) == "applied"
    monkeypatch.setitem(TP.VERIFIED_HOLDING_STATE, h, "raw")
    assert TP.holding_basis(h, "kpno_solar_atlas") == "line_selection"


def test_the_verified_map_carries_its_measurements():
    """HARPS molecfit is RYA-1192's VERIFIED-CORRECTED. Kitt Peak molecfit was RYA-1192's
    VERIFIED-RAW and is now corrected by RYA-1230's full-coverage molecfit, which
    superseded that verdict -- so it is checked against the corrected product on disk,
    not against the pre-correction artifact."""
    doc = json.loads((ROOT / "data/results/rya1192/rya1192_verification.json").read_text())
    probe = {k: v["probe"] for k, v in doc["availability_on_this_machine"].items()}
    assert probe["solar_harps_molecfit_corrected"] == "VERIFIED-CORRECTED"
    assert TP.VERIFIED_HOLDING_STATE["solar_harps_molecfit_corrected"] == "corrected"
    kp = ROOT / "data/processed/kp1984_telluric_corrected"
    assert any(kp.glob("kp1984_corrected_3000_*.txt")), "RYA-1230's blue correction is gone"
    assert TP.VERIFIED_HOLDING_STATE["solar_kpno_molecfit_corrected"] == "corrected"


def test_a_registry_applied_correction_is_believed():
    """RYA-1233, Ryan 2026-10-02: a telluric-corrected holding is the data; an incomplete
    correction is fixed, never answered by excluding its lines. This test used to pin the
    opposite (`test_a_claimed_but_unverified_correction_is_not_believed`)."""
    h = "solar_crires_plus_h_rya1094"
    assert TP.applied_state(h) == "applied"
    assert h not in TP.VERIFIED_HOLDING_STATE
    assert TP.holding_basis(h, "crires_plus") == "corrected"
    assert not TP.exclusion(IN_BAND_A, "crires_plus", h)


def test_above_the_atmosphere_stays_a_per_instrument_fact():
    """not_applicable is genuinely per-INSTRUMENT — no holding of a space telescope has
    tellurics — and must survive the re-key."""
    for h in ("procyon_stis", "alpha_cen_a_stis"):
        assert TP.holding_basis(h, "hst_stis") == "not_applicable"
        assert not TP.exclusion(IN_BAND_A, "hst_stis", h)


def test_basis_stays_per_instrument_and_says_so():
    """The other half of the re-key: `basis()` answers a DIFFERENT question and re-keying
    it would collapse the two axes this module forbids collapsing. It must also tell a
    caller who passes a holding what they did wrong."""
    assert TP.basis("harps") == "line_selection"
    assert TP.basis("iag_fts_solar_atlas") == "corrected"
    with pytest.raises(KeyError) as e:
        TP.basis("solar_harps_molecfit_corrected")
    assert "is a HOLDING, not an instrument" in str(e.value)
    assert "holding_basis" in str(e.value)


def test_an_unregistered_holding_is_refused_not_assumed(monkeypatch):
    """RYA-806: a holding whose state nobody recorded must not default to anything."""
    with pytest.raises(KeyError):
        TP.holding_basis("no_such_holding_rya1194", "harps")


# ── no value moved ───────────────────────────────────────────────────────────
def test_omitting_the_holding_reproduces_the_old_instrument_behaviour():
    """Every caller that cannot name a holding must be unaffected by this ticket, and the
    fallback can only over-exclude — never unlock a band."""
    cat = pd.read_csv(TP.CATALOG)
    for inst in cat.instrument_id.astype(str):
        old = TP.basis(inst) in ("corrected", "not_applicable")
        assert bool(TP.exclusion(IN_BAND_A, inst)) == (not old), inst


def test_applied_holdings_keep_their_lines_and_raw_ground_holdings_do_not():
    """Holding-keyed: every registry-`applied` holding keeps an in-band line; a raw ground
    holding (not-applied, atmosphere in the path) excludes it."""
    reg = pd.read_csv(TP.HOLDINGS)
    for _, r in reg.iterrows():
        h, i = str(r.holding_id), str(r.instrument_id)
        try:
            b = TP.basis(i)
        except KeyError:
            continue
        if b in ("not_applicable", "corrected"):
            continue
        if str(r.telluric_applied) == "applied" and TP.VERIFIED_HOLDING_STATE.get(h) != "raw":
            assert not TP.exclusion(IN_BAND_A, i, h), h
        elif str(r.telluric_applied) == "not-applied":
            assert TP.exclusion(IN_BAND_A, i, h), h
