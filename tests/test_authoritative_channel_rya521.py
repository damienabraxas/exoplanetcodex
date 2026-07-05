"""
tests/test_authoritative_channel_rya521.py
RYA-521 — one authoritative channel (phase_c verdict); raw EW is diagnostic-only;
the divergence guard fires on a raw-vs-verdict mismatch (the C=10.260 class).
"""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import authoritative_channel as ac  # noqa: E402


def _write_verdict(tmp, per_element):
    p = tmp / "solar_verdict.json"
    p.write_text(json.dumps({"summary": {"star": "solar"}, "verdicts": per_element}))
    return p


def _write_raw(tmp, rows):
    p = tmp / "solar_abundances.csv"
    lines = ["# DIAGNOSTIC-ONLY test", "element,ion,A_X,authoritative"]
    lines += [f"{el},I,{a},{auth}" for el, a, auth in rows]
    p.write_text("\n".join(lines) + "\n")
    return p


def test_divergence_guard_fires_on_flagship_mismatch(tmp_path):
    # raw EW says C=10.260 (the saturated-line artifact); verdict says C=8.491 PASS.
    verdict = _write_verdict(tmp_path, [
        {"element": "C", "A_measured": 8.491, "verdict": "PASS", "channel": "synthesis"},
        {"element": "Fe", "A_measured": 7.516, "verdict": "PASS", "channel": "EW"},
    ])
    raw = _write_raw(tmp_path, [("C", 10.260, True), ("Fe", 7.516, True)])
    flags = ac.channel_divergence(raw_path=raw, verdict_path=verdict, tol=0.10)
    assert [f["element"] for f in flags] == ["C"]           # C flagged, Fe agrees
    c = flags[0]
    assert c["pass_element"] is True and c["delta"] == pytest.approx(1.769, abs=1e-3)


def test_suppressed_rows_are_ignored(tmp_path):
    # a synthesis-required row marked authoritative=False must not enter the guard.
    verdict = _write_verdict(tmp_path, [
        {"element": "C", "A_measured": 8.491, "verdict": "PASS", "channel": "synthesis"}])
    raw = _write_raw(tmp_path, [("C", 10.260, False)])      # RYA-520 suppressed
    assert ac.channel_divergence(raw_path=raw, verdict_path=verdict) == []


def test_load_verdict_abundances(tmp_path):
    verdict = _write_verdict(tmp_path, [
        {"element": "O", "A_measured": 8.735, "verdict": "PASS", "channel": "synthesis"},
        {"element": "V", "A_measured": None, "verdict": "CURATION-OWED", "channel": "HFS"}])
    ab = ac.load_verdict_abundances(path=verdict)
    assert ab["O"]["A"] == 8.735 and ab["O"]["verdict"] == "PASS"
    assert ab["V"]["A"] is None                              # owed → no point value


def test_assert_authoritative_rejects_raw_ew():
    with pytest.raises(RuntimeError):
        ac.assert_authoritative_is_verdict("data/outputs/solar/solar_abundances.csv")
    # reading the frozen reference (verdict-sourced) is allowed
    ac.assert_authoritative_is_verdict("data/reference/solar/solar_abundances_v2.csv")
