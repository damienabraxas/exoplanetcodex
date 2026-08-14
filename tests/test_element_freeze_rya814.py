"""
tests/test_element_freeze_rya814.py — RYA-814
=============================================
The round-trip gate and per-element immutability.

The round-trip is the guardrail that keeps Lever A a pure REPRESENTATION refactor:
if exploding gold and re-assembling it is not byte-identical, a value moved, and
the decomposition is wrong. That claim has to be unfalsifiable-if-wrong, so it is
asserted on the live frozen artifact rather than on a fixture.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.element_freeze import (  # noqa: E402
    GOLD_COMPLETE_COUNT, ElementFreezeError, assemble, assembly_status,
    freeze_element, frozen_elements, read_element, read_preamble,
    serialise_fields, verify_consistency)

SOLAR = ROOT / "data" / "reference" / "solar"
SCRIPT = ROOT / "scripts" / "decouple_gold_rya814.py"


def _current() -> str:
    return (SOLAR / "CURRENT").read_text().strip()


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


# ── THE GATE ─────────────────────────────────────────────────────────────────
def test_roundtrip_is_byte_identical_to_the_live_frozen_gold():
    """Explode -> reassemble -> sha256 equality, against the REAL current gold."""
    r = subprocess.run([sys.executable, str(SCRIPT), "--roundtrip"],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "PASS — byte-identical" in r.stdout


def test_assembly_from_disk_records_equals_the_monolith():
    """The committed per-element records must reassemble to the frozen file."""
    v = _current()
    original = (SOLAR / f"solar_abundances_{v}.csv").read_bytes()
    rebuilt = assemble(version_label=v).encode("utf-8")
    assert _sha(rebuilt) == _sha(original)
    assert rebuilt == original          # byte equality, not just digest agreement


def test_every_record_fields_reserialise_to_its_verbatim_line():
    """`fields` may never drift from the bytes it claims to describe."""
    cols = read_preamble()["columns"]
    els = frozen_elements()
    assert els, "no frozen element records"
    for el in els:
        verify_consistency(read_element(el), cols)   # raises on drift


def test_every_record_hash_matches_its_own_line():
    for el in frozen_elements():
        rec = read_element(el)
        assert rec.sha256 == _sha(rec.verbatim.encode("utf-8")), el


def test_row_order_is_preserved_because_it_is_part_of_the_bytes():
    v = _current()
    body = [l for l in (SOLAR / f"solar_abundances_{v}.csv").read_text(
        encoding="utf-8").splitlines() if l and not l.startswith("#")][1:]
    original_order = [l.split(",")[0] for l in body]
    recs = sorted((read_element(e) for e in frozen_elements()),
                  key=lambda r: r.row_index)
    assert [r.element for r in recs] == original_order


# ── per-element immutability (RYA-469, decoupled not weakened) ───────────────
def test_refreezing_same_element_same_version_is_refused():
    el = frozen_elements()[0]
    rec = read_element(el)
    with pytest.raises(ElementFreezeError, match="WRITE-ONCE"):
        freeze_element(el, rec.fields, rec.verbatim, version=rec.version)


def test_freezing_one_element_does_not_touch_another(tmp_path, monkeypatch):
    """The whole point of Lever A: closing Fe must not disturb Ba."""
    import pipeline.element_freeze as ef
    monkeypatch.setattr(ef, "ELEMENTS_DIR", tmp_path)
    cols = ["element", "value"]
    monkeypatch.setattr(ef, "read_preamble", lambda: {
        "columns": cols, "comment_lines": ["# t"], "trailing_newline": True})

    ef.freeze_element("Fe", {"element": "Fe", "value": "7.466", "_row_index": 0},
                      "Fe,7.466", version="v1", columns=cols)
    ef.freeze_element("Ba", {"element": "Ba", "value": "2.237", "_row_index": 1},
                      "Ba,2.237", version="v1", columns=cols)
    ba_before = ef.read_element("Ba")

    # close Fe again at a NEW version — Ba must be untouched, bytes and hash
    ef.freeze_element("Fe", {"element": "Fe", "value": "7.470", "_row_index": 0},
                      "Fe,7.470", version="v2", columns=cols, allow_new_version=True)
    ba_after = ef.read_element("Ba")
    assert (ba_after.sha256, ba_after.verbatim, ba_after.version) == (
        ba_before.sha256, ba_before.verbatim, ba_before.version)
    assert ef.read_element("Fe").version == "v2"


def test_a_tampered_record_is_caught_not_silently_assembled(tmp_path, monkeypatch):
    import pipeline.element_freeze as ef
    monkeypatch.setattr(ef, "ELEMENTS_DIR", tmp_path)
    cols = ["element", "value"]
    monkeypatch.setattr(ef, "read_preamble", lambda: {
        "columns": cols, "comment_lines": ["# t"], "trailing_newline": True})
    ef.freeze_element("Fe", {"element": "Fe", "value": "7.466", "_row_index": 0},
                      "Fe,7.466", version="v1", columns=cols)
    p = tmp_path / "Fe.json"
    d = json.loads(p.read_text())
    d["verbatim"] = "Fe,9.999"                 # edit the bytes, leave the hash
    p.write_text(json.dumps(d))
    with pytest.raises(ElementFreezeError, match="altered on disk|re-serialise"):
        ef.assemble()


# ── the completeness bar (principle 6) ───────────────────────────────────────
def test_partial_assembly_labels_itself_partial():
    st = assembly_status()
    assert st["required_for_gold_complete"] == GOLD_COMPLETE_COUNT == 28
    if st["elements"] < GOLD_COMPLETE_COUNT:
        assert not st["gold_complete"]
        assert "PARTIAL" in st["label"]
        assert "NOT 'gold'" in st["label"]


def test_current_gold_is_below_the_28_element_bar():
    """
    26 frozen today against a ratified bar of 28.

    The repo carries 26 everywhere -- gold v4 has 26 rows and
    nlte_two_engine_coverage.csv has the same 26 -- so reaching 28 needs Zn
    (RYA-757, definitely absent) AND one further element the repo does not name.
    This test pins the gap as a FINDING rather than letting the bar be quietly
    lowered to match the data.
    """
    st = assembly_status()
    assert st["elements"] == 26
    assert st["missing_count"] == 2
    assert "Zn" not in st["frozen"]
    assert not st["gold_complete"]


# ── gold stays output-only (RYA-813's invariant must not regress) ────────────
def test_freeze_path_does_not_import_the_validation_stage():
    """Assembly WRITES a view; it must never pull validation back into the loop."""
    src = (ROOT / "pipeline" / "element_freeze.py").read_text()
    assert "validate_element" not in src.split('"""', 2)[-1], (
        "element_freeze must not call the validator — gold is written FROM verdicts, "
        "never read back INTO validation (RYA-812 principle 5)")
