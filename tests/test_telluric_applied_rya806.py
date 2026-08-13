"""RYA-806 — the per-holding `telluric_applied` switch.

The thing under test is a REFUSAL, so most of these assert that something does not
happen. The risk with that shape is a test that passes because the code never ran, so
each refusal test also pins the REASON, and the positive cases assert the gate lets real
work through — a gate that refuses everything would satisfy a careless suite.
"""
from __future__ import annotations

import csv
import io
from pathlib import Path

import pytest

from pipeline import telluric_policy as TP
from pipeline.telluric_intake import APPLIED, NOT_APPLIED, UNKNOWN, VALUES

ROOT = Path(__file__).resolve().parent.parent
HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"


def _rows():
    text = HOLDINGS.read_bytes().decode("utf-8")
    return list(csv.DictReader(io.StringIO(text, newline="")))


@pytest.fixture(autouse=True)
def _clear_cache():
    TP._catalog_cache.clear()
    yield
    TP._catalog_cache.clear()


# ── the column exists and every value is legal ───────────────────────────────
def test_every_holding_carries_a_legal_telluric_applied_value():
    rows = _rows()
    assert rows, "holdings registry is empty"
    for r in rows:
        v = (r.get("telluric_applied") or "").strip()
        assert v in VALUES, (
            f"{r['holding_id']} has telluric_applied={v!r}; legal values are {VALUES}. "
            f"An unrecognised value reads as `unknown` and refuses, but it should never "
            f"get into the registry in the first place.")


def test_no_IR_holding_is_blank_or_unknown():
    """The ticket's smoke test: an instrument that NEEDS correction must have a
    determined per-product state, or the switch it gates is meaningless."""
    unresolved = []
    for r in _rows():
        if TP.requires_correction(r["instrument_id"]):
            v = (r.get("telluric_applied") or "").strip()
            if v in ("", UNKNOWN):
                unresolved.append((r["holding_id"], v or "<blank>"))
    assert not unresolved, (
        f"IR holdings with no determined telluric state: {unresolved}. Determine each "
        f"with pipeline.telluric_intake.from_headers() and record it (RYA-806).")


def test_the_ir_holdings_are_actually_covered():
    """Guards the guard: if `requires_correction` ever returned False for everything,
    the test above would pass vacuously."""
    ir = [r["holding_id"] for r in _rows()
          if TP.requires_correction(r["instrument_id"])]
    assert len(ir) >= 5, f"expected the IR holdings to be found, got {ir}"


# ── the two axes are orthogonal and are NOT collapsed ────────────────────────
def test_same_instrument_axis_different_product_states():
    """The case that proves the axes cannot be merged: alpha Cen CRIRES+ and alpha Cen
    NIRPS are BOTH telluric_required=yes, and their products differ."""
    by_id = {r["holding_id"]: r for r in _rows()}
    crires = by_id["alpha_cen_a_crires_plus"]
    nirps = by_id["alpha_cen_b_nirps"]
    assert TP.requires_correction(crires["instrument_id"])
    assert TP.requires_correction(nirps["instrument_id"])
    assert crires["telluric_applied"] == NOT_APPLIED
    assert nirps["telluric_applied"] == APPLIED, (
        "NIRPS S1D_FINAL_A carries FLUX_TELL_* whose ratio to FLUX_EL tracks "
        "1/ATM_TRANSM at r=1.0000 — it IS corrected. If the holding axis were derived "
        "from the instrument axis this would be wrong, which is the point.")


def test_applied_skips_the_stage_and_not_applied_routes_to_it():
    ok, why = TP.gate_holding("alpha_cen_b_nirps")
    assert ok and "SKIPPED" in why

    ok, why = TP.gate_holding("solar_vesta_crires_plus_idp")
    assert not ok
    assert "RYA-424" in why and "not-applied" in why


def test_line_selection_instrument_runs_uncorrected_on_a_stated_basis():
    """Kitt Peak is not-applied and runs anyway — per-line clean-line selection is a
    ratified method (RYA-460/786), not uncorrected data sneaking through."""
    ok, why = TP.gate_holding("solar_kpno")
    assert ok
    assert "clean-line selection" in why
    assert "No telluric_corrected declaration is made" in why


def test_not_applicable_instrument_is_not_gated_by_the_product_axis():
    """HST is above the atmosphere. Its holdings are literally `not-applied`, and that
    must not refuse them — the instrument axis short-circuits."""
    row = next(r for r in _rows() if r["holding_id"] == "procyon_stis")
    assert row["telluric_applied"] == NOT_APPLIED
    ok, why = TP.gate_holding("procyon_stis")
    assert ok and "not_applicable" in why


# ── `unknown` refuses, loudly, and is never defaulted ────────────────────────
def test_unknown_on_an_IR_holding_raises(tmp_path, monkeypatch):
    """No real IR holding is `unknown` today, so this synthesises one. Without it the
    loud-fail path would be untested precisely because the backfill succeeded."""
    text = HOLDINGS.read_bytes().decode("utf-8")
    rows = list(csv.reader(io.StringIO(text, newline="")))
    hdr = rows[0]
    i_id, i_app = hdr.index("holding_id"), hdr.index("telluric_applied")
    for r in rows[1:]:
        if r and r[i_id] == "solar_vesta_crires_plus_idp":
            r[i_app] = UNKNOWN
    fake = tmp_path / "holdings.csv"
    buf = io.StringIO(newline="")
    csv.writer(buf, lineterminator="\n").writerows(rows)
    fake.write_text(buf.getvalue())

    monkeypatch.setattr(TP, "HOLDINGS", fake)
    TP._catalog_cache.clear()
    with pytest.raises(TP.TelluricStateUnknown) as ei:
        TP.gate_holding("solar_vesta_crires_plus_idp")
    msg = str(ei.value)
    assert "UNVERIFIED" in msg
    assert "never defaulted" in msg
    assert "telluric_intake" in msg, "the refusal must say how to resolve itself"


def test_unknown_does_not_ground_a_non_IR_instrument():
    """HARPS Phase-3 ADPs carry no PRO REC chain, so their telluric state is honestly
    unknown. Refusing there would ground the main EW arm on a switch built for the IR —
    and it is not a silent default, because per-line selection is valid on uncorrected
    data whatever the product turns out to be."""
    row = next(r for r in _rows() if r["holding_id"] == "alpha_cen_a_harps")
    assert row["telluric_applied"] == UNKNOWN
    ok, why = TP.gate_holding("alpha_cen_a_harps")
    assert ok
    assert "UNKNOWN and stays unknown" in why


def test_an_unregistered_holding_is_loud():
    with pytest.raises(KeyError):
        TP.gate_holding("no_such_holding_anywhere")


# ── the loader gap RYA-805 found is closed ───────────────────────────────────
def test_load_window_refuses_the_uncorrected_IR_arm():
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    import measure_band_ew as M

    # The gate fires BEFORE any file is opened, which is why this test needs no staged
    # data: a refusal must not depend on the data being reachable, or it would report
    # "not staged" for a dataset we hold but may not measure.
    with pytest.raises(M.TelluricNotCorrected) as ei:
        M.load_window("crires_plus", 15000.0, 1.0)
    assert "RYA-424" in str(ei.value)
    assert "allow_uncorrected=True" in str(ei.value), (
        "the refusal must name the door the correction leg uses, or RYA-373 is locked "
        "out of the data it exists to correct")
    # That door, and the rest-frame gate still standing behind it, are exercised against
    # staged frames in tests/test_crires_window_loader_rya796.py.


def test_load_window_still_serves_the_optical_arms():
    """A gate that refuses everything would pass the refusal tests above."""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    import measure_band_ew as M

    w, f, prov = M.load_window("kpno_solar_atlas", 6300.0, 1.0)
    assert len(w) and len(w) == len(f) and isinstance(prov, str)
