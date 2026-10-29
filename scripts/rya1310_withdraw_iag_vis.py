#!/usr/bin/env python3
"""RYA-1310 -- withdraw every live IAG product in the VIS band to its feed's `archive`.

Ryan 2026-10-28: "I don't trust the VIS IAG at all ... we were primarily using IAG for only
Near IR, because we have proper VIS representation already." The orchestrator no longer
builds IAG VIS cells (run_matrix.EXCLUDED_HOLDING_BANDS); this removes the products the
earlier runs left live. Archived, never deleted, through publish_product.write_feed.

    python3 scripts/rya1310_withdraw_iag_vis.py            # report only
    python3 scripts/rya1310_withdraw_iag_vis.py --apply
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
from pipeline.run_matrix import EXCLUDED_HOLDING_BANDS  # noqa: E402

REASON = ("RYA-1310: IAG is not used in the VIS band -- VIS is served by HARPS and Kitt Peak, "
          "IAG is the near-IR arm (Ryan 2026-10-28: \"I don't trust the VIS IAG at all\"). "
          "Withdrawn, not re-measured; the value is kept here for the record.")


def excluded(p: dict) -> bool:
    return (str(p.get("instrument")), str(p.get("band"))) in EXCLUDED_HOLDING_BANDS


def withdraw(path: Path, *, apply: bool) -> list[dict]:
    doc = json.loads(path.read_text())
    hit = [p for p in doc.get("products") or [] if excluded(p)]
    if apply and hit:
        now = pp._now()
        doc.setdefault("archive", []).extend(
            dict(p, superseded_at=now, superseded_reason=REASON) for p in hit)
        doc["products"] = [p for p in doc["products"] if not excluded(p)]
        doc["version"] = pp.bump(doc["version"])
        doc["updated_at"] = now
        pp.write_feed(path, doc)
    return hit


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--star", default="solar")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    total = 0
    for path in sorted((pp.STORE / a.star).glob("*.json")):
        hit = withdraw(path, apply=a.apply)
        if hit:
            total += len(hit)
            ver = json.loads(path.read_text())["version"]
            print(f"{path.stem}: {len(hit)} IAG VIS product(s) "
                  f"{'withdrawn -> v' + ver if a.apply else 'live'}")
    print(f"total {total}" + ("" if a.apply else " (report only; --apply to withdraw)"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
