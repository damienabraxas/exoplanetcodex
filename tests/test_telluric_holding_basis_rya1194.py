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
def test_the_rya1191_forward_path_is_recognised(monkeypatch):
    """🔴 THE POINT. Today KP-molecfit is verified RAW, so it excludes — the same answer
    the instrument-keyed code gave, for a completely different reason. Once RYA-1191
    corrects the product and re-verifies, the entry flips and the pipeline must follow it
    WITHOUT any further code change. If this test fails, RYA-1191's fix is invisible."""
    h, i = "solar_kpno_molecfit_corrected", "kpno_solar_atlas"
    assert TP.holding_basis(h, i) == "line_selection"
    assert TP.exclusion(IN_BAND_A, i, h), "today: verified RAW, so the line is quarantined"

    monkeypatch.setitem(TP.VERIFIED_HOLDING_STATE, h, "corrected")
    assert TP.holding_basis(h, i) == "corrected", (
        "RYA-1191 corrected the holding and the policy did not notice — this is exactly "
        "the defect RYA-1194 was filed to prevent")
    assert not TP.exclusion(IN_BAND_A, i, h)
    # and the raw sibling must NOT have moved with it
    assert TP.holding_basis("solar_kpno", i) == "line_selection"


def test_a_measurement_outranks_the_registry_label():
    """RYA-1192's finding, as a rule. The registry says this holding is corrected; the
    bytes say it is raw on 146/146 lines. The bytes win."""
    h = "solar_kpno_molecfit_corrected"
    assert TP.applied_state(h) == "applied", "premise gone: the registry no longer claims it"
    assert TP.VERIFIED_HOLDING_STATE[h] == "raw"
    assert TP.holding_basis(h, "kpno_solar_atlas") == "line_selection", (
        "a registry label was believed over a measurement that contradicts it")


def test_the_verified_map_matches_rya1192s_committed_artifact():
    """The constant is a transcription of a measurement, so it must still agree with the
    measurement. Otherwise it is just an opinion with a ticket number on it."""
    doc = json.loads((ROOT / "data/results/rya1192/rya1192_verification.json").read_text())
    probe = {k: v["probe"] for k, v in doc["availability_on_this_machine"].items()}
    expect = {"VERIFIED-RAW": "raw", "VERIFIED-CORRECTED": "corrected"}
    checked = 0
    for holding, state in TP.VERIFIED_HOLDING_STATE.items():
        assert holding in probe, f"{holding} is not in RYA-1192's output at all"
        assert expect.get(probe[holding]) == state, (
            f"{holding}: map says {state!r}, RYA-1192 measured {probe[holding]!r}")
        checked += 1
    assert checked == 2, "both RYA-1192 verdicts must be carried"


# ── the conservative direction ───────────────────────────────────────────────
def test_a_claimed_but_unverified_correction_is_not_believed():
    """`applied` in the registry is a LABEL, and RYA-1192 caught one that was false. An
    unverified claim excludes — and says so distinguishably, so somebody can close it."""
    h = "solar_crires_plus_h_rya1094"
    assert TP.applied_state(h) == "applied"
    assert h not in TP.VERIFIED_HOLDING_STATE
    assert TP.holding_basis(h, "crires_plus") == TP.APPLIED_UNVERIFIED
    why = TP.exclusion(IN_BAND_A, "crires_plus", h)
    assert why, "an unverified correction claim must not unlock a telluric band"
    assert "no measurement has confirmed" in why, (
        "the reason must distinguish 'we measured tellurics here' from 'nobody checked'")
    assert TP.APPLIED_UNVERIFIED != "corrected"


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


def test_only_one_holding_changed_decision_and_it_is_the_verified_one():
    """The measured blast radius, pinned. Every registered holding is compared instrument-
    keyed (old) against holding-keyed (new); exactly one may differ, and only because a
    measurement says so. A second entry here is a value move that needs its own ticket."""
    reg = pd.read_csv(TP.HOLDINGS)
    moved = []
    for _, r in reg.iterrows():
        h, i = str(r.holding_id), str(r.instrument_id)
        try:
            old = bool(TP.exclusion(IN_BAND_A, i))
        except KeyError:
            continue                      # instrument not in the catalog at all
        if old != bool(TP.exclusion(IN_BAND_A, i, h)):
            moved.append(h)
    assert moved == ["solar_harps_molecfit_corrected"], moved
    assert TP.VERIFIED_HOLDING_STATE["solar_harps_molecfit_corrected"] == "corrected"
