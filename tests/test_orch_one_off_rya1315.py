"""RYA-1315: the orchestrator handles a one-off -- select band / instrument / engine, force the
cells a library change made stale (inputs_hash cannot see it), and publish ONLY those cells."""
import json

import pytest

from pipeline import element_publish as ep
from pipeline import run_matrix as rm


def test_force_needs_a_narrowing_filter():
    with pytest.raises(rm.MatrixError, match="narrow"):
        rm.run("solar", "Fe", force=True, dry_run=True, echo=False)


def test_a_forced_rebuild_gets_fresh_legs():
    assert ep._built_from({"inputs_hash": "h"}) == "h"
    forced = ep._built_from({"inputs_hash": "h", "forced_at": "2026-10-31T00:00:00Z"})
    assert forced != "h" and forced.startswith("h ")
    assert ep._built_from(None) is None


def test_record_inputs_hash_stamps_a_forced_build(tmp_path, monkeypatch):
    monkeypatch.setattr(rm, "LEDGER", tmp_path / "inputs_hashes.json")
    rm.record_inputs_hash("k", "h", product_keys=[], code_commit="c", forced=True)
    rm.record_inputs_hash("j", "h", product_keys=[], code_commit="c")
    led = json.loads(rm.LEDGER.read_text())
    assert led["k"]["forced_at"] == led["k"]["recorded_at"] and "forced_at" not in led["j"]


def test_scope_cells_cover_only_what_the_run_built():
    rep = {"cells": [
        {"ion": "I", "band": "VIS", "holding": "h", "pool": "REFERENCE", "engine": "ts-lte",
         "status": rm.DONE},
        {"ion": "I", "band": "VIS", "holding": "h", "pool": "REFERENCE", "engine": "gerber-nlte",
         "status": rm.FAILED}]}
    assert ep._scope_cells(rep) == {("I", "VIS", "h", "REFERENCE"): set(rm.DECK_EMITS["ts-lte"])}


def _feed(tmp_path, monkeypatch, products):
    import sys
    sys.path.insert(0, str(ep.ROOT / "scripts"))
    import publish_product as pp
    monkeypatch.setattr(pp, "STORE", tmp_path / "store")
    (tmp_path / "store" / "solar").mkdir(parents=True)
    path = tmp_path / "store" / "solar" / "Fe.json"
    path.write_text(json.dumps({"version": "1.0", "products": products}))
    monkeypatch.setattr(pp, "write_feed", lambda p, d: p.write_text(json.dumps(d)))
    return path


def _p(holding, selector, treatment, ion="I"):
    return {"route": "SYNTH", "star": "solar", "ion": ion, "band": "VIS", "instrument": "x",
            "holding": holding, "selector": selector, "treatment": treatment, "A": 7.5,
            "tier": selector if selector in ("REFERENCE", "GRADED") else "ALL"}


def test_a_scoped_publish_withdraws_only_what_it_held_in_scope(tmp_path, monkeypatch):
    held_in = _p("harps", "REFERENCE", "1D-LTE")
    not_reemitted = _p("harps", "REFERENCE", "ENGINE-A-3DNLTE")
    other_holding = _p("kpno", "REFERENCE", "1D-LTE")
    special = _p("harps", "ASPLUND_AGSS21", "ENGINE-A-3DNLTE")
    path = _feed(tmp_path, monkeypatch, [held_in, not_reemitted, other_holding, special])
    (tmp_path / "stage").mkdir()
    gone = ep._withdraw_superseded(
        "solar", "Fe", stage=tmp_path / "stage", keep=[],
        held=[{"product": "VIS harps REFERENCE 1D-LTE", "why": "budget over gate"}],
        why="t", started="2026-10-31", scope={("I", "VIS", "harps", "REFERENCE"):
                                              set(rm.DECK_EMITS["ts-lte"])})
    left = json.loads(path.read_text())["products"]
    assert len(gone) == 1 and held_in not in left
    assert not_reemitted in left and other_holding in left and special in left
