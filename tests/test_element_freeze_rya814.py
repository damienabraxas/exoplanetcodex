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
    freeze_element, frozen_elements, read_element, read_preamble, read_record,
    record_key, serialise_fields, split_key, verify_consistency)

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
        verify_consistency(read_record(el), cols)   # raises on drift


def test_every_record_hash_matches_its_own_line():
    for key in frozen_elements():
        rec = read_record(key)
        assert rec.sha256 == _sha(rec.verbatim.encode("utf-8")), key


def test_row_order_is_preserved_because_it_is_part_of_the_bytes():
    v = _current()
    body = [l for l in (SOLAR / f"solar_abundances_{v}.csv").read_text(
        encoding="utf-8").splitlines() if l and not l.startswith("#")][1:]
    original_order = [l.split(",")[0] for l in body]
    recs = sorted((read_record(k) for k in frozen_elements()),
                  key=lambda r: r.row_index)
    assert [r.element for r in recs] == original_order


# ── per-element immutability (RYA-469, decoupled not weakened) ───────────────
def test_refreezing_same_species_same_version_is_refused():
    key = frozen_elements()[0]
    el, ion = split_key(key)
    rec = read_record(key)
    with pytest.raises(ElementFreezeError, match="WRITE-ONCE"):
        freeze_element(el, rec.fields, rec.verbatim, version=rec.version, ion=ion)


def test_passing_a_record_key_where_an_element_is_expected_is_refused():
    """
    Caught by the suite: frozen_elements() returns KEYS, these take an ELEMENT.
    Passing 'Al_I' produced a spurious 'Al_I_I' record and inflated the species
    count from 26 to 27. Refuse the shape rather than accept both.
    """
    rec = read_record(frozen_elements()[0])
    with pytest.raises(ElementFreezeError, match="not a record key"):
        freeze_element("Al_I", rec.fields, rec.verbatim, version="v9")


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


def test_the_26_of_28_gap_is_exactly_fe_ii_and_zn():
    """
    RESOLVED (RYA-109 + RYA-757): the canonical 27 counts Fe I and Fe II
    SEPARATELY (Fe II is #2 at priority 1 -- log g via ionisation equilibrium),
    and Zn is the 28th. So 27 + Zn = 28, gold carries 26 = 27 - Fe II, and the
    two missing are Fe II and Zn. Not an off-by-one; a species-vs-element
    distinction.
    """
    st = assembly_status()
    assert st["elements"] == 26
    assert st["missing_count"] == 2
    assert not st["gold_complete"]
    keys = set(st["frozen"])
    assert "Fe_I" in keys and "Fe_II" not in keys   # Fe II adopted (RYA-406), unfrozen
    assert not any(k.startswith("Zn") for k in keys)
    assert "Fe II" in st["label"] and "Zn" in st["label"]


# ── (element, ion) keying — Fe I and Fe II are different species ─────────────
def test_record_key_is_element_and_ion():
    assert record_key("Fe", "I") == "Fe_I"
    assert record_key("Fe", "II") == "Fe_II"
    assert record_key("Fe", "I") != record_key("Fe", "II")


def test_gold_already_carries_an_ion_ii_species():
    """Sc II is frozen today, so the ion-bearing key is not hypothetical."""
    assert "Sc_II" in frozen_elements()


def test_fe_i_and_fe_ii_can_both_be_frozen_without_collision(tmp_path, monkeypatch):
    """
    The defect this keying fixes. Keyed by element alone, freezing Fe II would
    land on Fe.json and silently overwrite Fe I -- a frozen value changing without
    a version bump, exactly what RYA-469 forbids.
    """
    import pipeline.element_freeze as ef
    monkeypatch.setattr(ef, "ELEMENTS_DIR", tmp_path)
    cols = ["element", "ion", "A_X"]
    monkeypatch.setattr(ef, "read_preamble", lambda: {
        "columns": cols, "comment_lines": ["# t"], "trailing_newline": True})

    ef.freeze_element("Fe", {"element": "Fe", "ion": "I", "A_X": "7.466",
                             "_row_index": 0}, "Fe,I,7.466", version="v4", columns=cols)
    ef.freeze_element("Fe", {"element": "Fe", "ion": "II", "A_X": "7.486",
                             "_row_index": 1}, "Fe,II,7.486", version="v4", columns=cols)

    assert sorted(ef.frozen_elements()) == ["Fe_I", "Fe_II"]
    assert ef.read_element("Fe", "I").fields["A_X"] == "7.466"
    assert ef.read_element("Fe", "II").fields["A_X"] == "7.486"   # NOT overwritten
    # and an ambiguous element-only read must refuse rather than guess
    with pytest.raises(ElementFreezeError, match="disambiguate"):
        ef.read_element("Fe")


# ── gold stays output-only (RYA-813's invariant must not regress) ────────────
def test_freeze_path_does_not_import_the_validation_stage():
    """Assembly WRITES a view; it must never pull validation back into the loop."""
    src = (ROOT / "pipeline" / "element_freeze.py").read_text()
    assert "validate_element" not in src.split('"""', 2)[-1], (
        "element_freeze must not call the validator — gold is written FROM verdicts, "
        "never read back INTO validation (RYA-812 principle 5)")
