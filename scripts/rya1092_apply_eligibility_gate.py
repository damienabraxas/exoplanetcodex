#!/usr/bin/env python3
"""
scripts/rya1092_apply_eligibility_gate.py
=========================================
RYA-1092 — apply the live-product eligibility gate to an element feed.

    python3 scripts/rya1092_apply_eligibility_gate.py --element Fe            # report
    python3 scripts/rya1092_apply_eligibility_gate.py --element Fe --apply    # write

🔴 THE POOLS EXISTED; THE GATE DID NOT. `data/products/<star>/<El>.json` has carried
`quarantine[]` since RYA-711 and nothing ever checked `products[]` against the discipline
that governs it, so ineligible products stayed LIVE and the public feed rendered them.
This script is the one-time sweep; `tests/test_product_eligibility_rya1092.py` is the
standing enforcement that stops them coming back.

WHAT IT DOES, AND THE THREE THINGS IT REFUSES TO DO
---------------------------------------------------
It MOVES a failing product from `products[]` to `quarantine[]`, adding exactly three
keys: `quarantined_at`, `quarantine_reason` (prose), `quarantine_codes` (machine-readable).

  * It does not DELETE. A withdrawn product still has to be defensible in the appendix,
    so destroying it destroys the evidence for its own rejection (RYA-711).
  * It does not EDIT a published value. Every field of every moved record is carried
    across byte-identical; `test_no_published_value_was_edited_to_reconcile` (RYA-1080)
    checks that across ALL pools, so an edit smuggled in during a move fails CI.
  * It does not RANK. Where the gate empties a cell, that cell is reported as re-run debt
    and left empty. Choosing a replacement is a measurement decision, not a sweep's.

⚠️ THE SWEEP IS DERIVED, NOT LISTED. Nothing here names the Amarsi products, the eight
pre-fix products, or any holding. Everything comes from `pipeline.product_eligibility`
evaluating the record. A hardcoded list would have to be re-edited every time the feed
changes -- and, worse, would freeze today's verdict into the enforcement mechanism, which
is exactly what Ryan's 2026-08-28 correction to this ticket had to undo for Amarsi.
"""
from __future__ import annotations

import argparse
import collections
import datetime as _dt
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import product_eligibility as pe            # noqa: E402

STORE = ROOT / "data" / "products"


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _bump(v: str) -> str:
    major, _, minor = str(v).partition(".")
    return f"{major}.{int(minor or 0) + 1}"


def sweep(doc: dict, *, ticket: str = "RYA-1092") -> tuple[dict, list]:
    """Return (new_doc, moved) without touching the input."""
    doc = json.loads(json.dumps(doc))          # deep copy; the caller's doc is evidence
    verdicts = pe.evaluate_feed(doc)
    keep, moved = [], []
    stamp = _now()
    for p in doc["products"]:
        reasons = verdicts[pe.key_of(p)]
        if not reasons:
            keep.append(p)
            continue
        q = dict(p)
        q["quarantined_at"] = stamp
        q["quarantine_codes"] = [r.code for r in reasons]
        q["quarantine_reason"] = (
            f"{ticket} live-eligibility gate: "
            + " ".join(f"[{r.code}] {r.detail}" for r in reasons)
            + " Moved, not deleted -- every published field is carried across unchanged "
              "and this product returns to products[] the moment it passes the gate.")
        moved.append((p, reasons))
        doc.setdefault("quarantine", []).append(q)
    doc["products"] = keep
    return doc, moved


def report(doc: dict, moved: list) -> None:
    print(f"  live before : {len(moved) + len(doc['products'])}")
    print(f"  live after  : {len(doc['products'])}")
    print(f"  quarantined : {len(moved)}")
    codes = collections.Counter(r.code for _, rs in moved for r in rs)
    print("\n  reason breakdown (a product may carry more than one):")
    for c, n in codes.most_common():
        print(f"    {c:22s} {n}")
    print("\n  moved:")
    for p, rs in moved:
        print(f"    {pe.key_of(p)}")
        print(f"       A={p.get('A')} stat={p.get('sigma_stat')} syst={p.get('sigma_syst')} "
              f"n={p.get('n_lines')} -> {', '.join(r.code for r in rs)}")


def empty_cells(before: list, after: list) -> list:
    """Identities that had a live product and now have none -- re-run debt, reported so it
    cannot be discovered later as a mysterious hole (RYA-833). NOT declared as `gaps[]`:
    a gap means the cell CANNOT exist, and these can."""
    def cell(p):
        return (p["element"], p["ion"], p["band"], p["holding"], p["tier"],
                p["route"], p["treatment"])
    return sorted({cell(p) for p in before} - {cell(p) for p in after})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--star", default="solar")
    ap.add_argument("--restore", action="store_true",
                    help="also move quarantined products back to products[] when they "
                         "now pass the CURRENT gate. Needed after a gate is retired: the "
                         "sweep only ever moved products OUT, so retiring a criterion "
                         "left its victims sitting in quarantine forever.")
    ap.add_argument("--apply", action="store_true",
                    help="write the feed. Without it this only reports -- the default is "
                         "deliberately the harmless one.")
    a = ap.parse_args()

    path = STORE / a.star / f"{a.element}.json"
    if not path.exists():
        raise SystemExit(f"no feed at {path}")
    doc = json.loads(path.read_text())
    before = list(doc["products"])

    restored = []
    if a.restore:
        keep = []
        for q in doc.get("quarantine", []):
            # Re-evaluate against the CURRENT gate, on the record's own fields. The three
            # quarantine bookkeeping keys are not part of the product and are dropped on
            # the way back out; every other field is carried across byte-identical, which
            # test_no_published_value_was_edited_to_reconcile (RYA-1080) enforces.
            probe = {k: v for k, v in q.items()
                     if k not in ("quarantined_at", "quarantine_reason", "quarantine_codes")}
            if pe.evaluate(probe, peers=doc.get("products")):
                keep.append(q)
            else:
                restored.append(probe)
        doc["quarantine"] = keep
        doc["products"] = list(doc["products"]) + restored

    new, moved = sweep(doc)
    print(f"\n=== RYA-1092 eligibility gate — {a.star} {a.element} ===")
    if a.restore:
        print(f"  RESTORED from quarantine (now pass the current gate): {len(restored)}")
        for r in restored[:40]:
            print(f"    {r.get('band')} · {r.get('holding')} · {r.get('selector')} · "
                  f"{r.get('treatment')} · A={r.get('A')}")
        print(f"  still quarantined: {len(doc.get('quarantine', []))}")
    report(new, moved)

    holes = empty_cells(before, new["products"])
    if holes:
        print(f"\n  ⚠️ {len(holes)} cell(s) now have NO live product. This is RE-RUN DEBT, "
              f"not a gap:")
        for c in holes:
            print(f"    {' · '.join(str(x) for x in c)}")

    dups = pe.duplicate_live_cells(new)
    print(f"\n  duplicate live cells after the sweep: {len(dups)}")

    if not a.apply:
        print("\n  (report only -- pass --apply to write)")
        return 0
    new["version"] = _bump(new.get("version") or "1.0")
    new["updated_at"] = _now()
    path.write_text(json.dumps(new, indent=2) + "\n")
    print(f"\n  wrote {path.relative_to(ROOT)}  (version -> {new['version']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
