"""RYA-728: derive's "NOT ONE pooled line carries an NLTE label" refusal is a HOLD, not a FAILED.

8 Cr cells (HARPS Cr II VIS, Cr I red-optical; gerber-nlte / gerber-mean3d, graded pool of
1-2 lines) were refused correctly by derive and then filed as FAILED. The orchestrator
recognises the refusal by a marker string, so the marker must stay in derive's message.
"""
from pathlib import Path

from pipeline import run_matrix as rm

ROOT = Path(__file__).resolve().parents[1]
DERIVE = (ROOT / "scripts" / "derive_band_products.py").read_text(encoding="utf-8")


def test_every_marker_is_in_derives_refusal_message():
    for m in rm.NLTE_UNLABELLED_POOL_MARKERS:
        # derive builds the message across f-string pieces; the marker sits in one piece
        assert m in DERIVE, m


def test_the_marker_does_not_also_classify_as_an_empty_pool():
    msg = ("derive_products: exit 1: NLTE was requested for Cr I but NOT ONE of the 2 pooled "
           "lines carries an NLTE label, so bsyn would apply departure = 1")
    assert any(m in msg for m in rm.NLTE_UNLABELLED_POOL_MARKERS)
    assert not any(m in msg for m in rm.EMPTY_POOL_MARKERS)


# ── the model registry is fingerprinted by its emitting rows only (RYA-728) ──────────────
def _reg_text():
    return rm.MODEL_REGISTRY.read_text(encoding="utf-8")


def test_editing_a_row_with_no_token_does_not_move_the_fingerprint():
    import csv
    import io
    text = _reg_text()
    rows = list(csv.DictReader(io.StringIO(text)))
    idle = next(r for r in rows if not r["stored_token"].strip())
    edited = text.replace(idle["notes"], idle["notes"] + " EDITED", 1) if idle["notes"] else text
    assert edited != text
    assert rm.model_registry_fingerprint(text=edited) == rm.model_registry_fingerprint(text=text)


def test_editing_an_emitting_row_moves_the_fingerprint():
    import csv
    import io
    text = _reg_text()
    live = next(r for r in csv.DictReader(io.StringIO(text)) if r["stored_token"].strip())
    edited = text.replace(live["atmosphere"] + ",", live["atmosphere"] + "-x,", 1)
    assert rm.model_registry_fingerprint(text=edited) != rm.model_registry_fingerprint(text=text)


def test_migration_restamps_only_a_record_whose_emitting_rows_match():
    import hashlib
    import importlib.util
    import json
    spec = importlib.util.spec_from_file_location(
        "rescope", ROOT / "scripts/rya728_rescope_model_registry_fingerprint.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    text = _reg_text()
    old_text = text.replace("RYA-1118/1278", "RYA-1012-old", 1)     # an idle-row edit
    assert old_text != text
    d_old = hashlib.sha256(old_text.encode()).hexdigest()
    versions = {d_old: old_text}
    now = rm.model_registry_fingerprint(text=text)

    def entry(digest):
        rows = [{"kind": "spectrum", "name": "x", "digest": "s"},
                {"kind": "ledger", "name": m.REG, "digest": digest}]
        return {"inputs": rows,
                "inputs_hash": hashlib.sha256(json.dumps(rows, sort_keys=True).encode()).hexdigest()}

    e = entry(d_old)
    new, why = m.rescope(e, versions, now)
    assert why == "" and new == hashlib.sha256(json.dumps(e["inputs"], sort_keys=True).encode()).hexdigest()
    assert {"kind": "ledger_scoped", "name": m.SCOPED_NAME, "digest": now} in e["inputs"]
    # an emitting-row difference stays stale; an unknown digest stays stale
    live_old = text.replace("stagger-mean3d,", "stagger-x,", 1)
    d_live = hashlib.sha256(live_old.encode()).hexdigest()
    assert m.rescope(entry(d_live), {d_live: live_old}, now)[0] is None
    assert m.rescope(entry("nope"), versions, now)[0] is None
