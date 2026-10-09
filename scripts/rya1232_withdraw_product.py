#!/usr/bin/env python3
"""RYA-1232 -- withdraw one live product to its feed's `archive`, with the reason, through
publish_product.write_feed (the RYA-587 gate still applies).

    python3 scripts/rya1232_withdraw_product.py --element N --holding solar_iag \
        --selector MOL-CN_AX_IR --reason "..." [--apply]

`rya1230_publish_cno` only withdraws held BAND products that have a staged stem; a held
synthesis REGION (IAG CN, beyond the IAG cap) has none, so it would stay live unrevised.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT))

import publish_product as pp  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--element", required=True)
    ap.add_argument("--holding", required=True)
    ap.add_argument("--selector", required=True)
    ap.add_argument("--reason", required=True)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    path = ROOT / "data/products/solar" / f"{a.element}.json"
    doc = json.loads(path.read_text())
    hit = [p for p in doc["products"]
           if p.get("holding") == a.holding and str(p.get("selector")) == a.selector]
    for p in hit:
        print(f"  {pp.key_of(p)}  A={p.get('A')}")
    if not hit:
        raise SystemExit("no live product matches")
    if a.apply:
        now = pp._now()
        for p in hit:
            old = dict(p, superseded_at=now, superseded_reason=a.reason)
            doc.setdefault("archive", []).append(old)
        doc["products"] = [p for p in doc["products"] if p not in hit]
        doc["version"] = pp.bump(doc["version"])
        doc["updated_at"] = now
        pp.write_feed(path, doc)
        print(f"withdrew {len(hit)} -> {path.name} v{doc['version']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
