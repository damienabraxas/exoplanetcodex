#!/usr/bin/env python3
"""
RYA-1101 — verify `data/catalog/model_registry.csv` against the CODE.

This is the anti-memory guard. The registry is meant to BE the source of truth, so it must
not be able to launder an unsourced value: every `stored_token` on a live row has to be found
verbatim in `band_products.TREATMENTS`. A token that is in the registry and not in the code
means somebody typed it from memory, and that is a CRITICAL failure, not a warning.

Usage:
    python3 scripts/verify_model_registry.py
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import band_products, model_registry  # noqa: E402

FORBIDDEN_COLUMN_HINTS = ("count", "n_products", "n_lines", "live_count")


def main() -> int:
    rows = model_registry.load()
    fails: list[str] = []

    # (e) row count
    live = [r for r in rows if r["status"] == "live"]
    indev = [r for r in rows if r["status"] == "in-dev"]
    if len(rows) != 8 or len(live) != 7 or len(indev) != 1:
        fails.append(f"(e) expected 8 rows (7 live + 1 in-dev), got {len(rows)} "
                     f"({len(live)} live, {len(indev)} in-dev)")

    # (a) every live token exists verbatim in the code
    code_tokens = set(band_products.TREATMENTS)
    for r in live:
        tok = r["stored_token"].strip()
        if not tok:
            fails.append(f"(a) model {r['model_id']} is live with no stored_token")
        elif tok not in code_tokens:
            fails.append(f"(a) model {r['model_id']} token {tok!r} NOT FOUND in "
                         f"band_products.TREATMENTS — supplied from memory, not source")

    # (b) no token maps to two model_ids, except the documented bare alias
    dupes = [t for t, n in Counter(r["stored_token"].strip() for r in rows
                                   if r["stored_token"].strip()).items() if n > 1]
    for t in dupes:
        fails.append(f"(b) token {t!r} maps to more than one model_id — undocumented "
                     f"collision")
    note = ROOT / "docs" / "catalog" / "model_registry_notes.md"
    for alias in model_registry.AMBIGUOUS_TOKENS:
        if not note.exists() or alias not in note.read_text():
            fails.append(f"(b) ambiguous token {alias!r} is not documented in {note.name}")
        try:
            model_registry.resolve(alias)
            fails.append(f"(b) bare {alias!r} RESOLVED — it must be guarded, not bound")
        except model_registry.ModelRegistryError:
            pass

    # (c) the Frankenstein pair
    frank = [r for r in rows if r["pair_group"] == "frankenstein"]
    if {r["model_id"] for r in frank} != {"5", "6"}:
        fails.append(f"(c) pair_group=frankenstein is "
                     f"{sorted(r['model_id'] for r in frank)}, expected models 5 and 6")

    # (d) no live-count column
    for col in (rows[0].keys() if rows else []):
        if any(h in col.lower() for h in FORBIDDEN_COLUMN_HINTS):
            fails.append(f"(d) column {col!r} looks like a live count — those belong to "
                         f"RYA-1015, never here")

    if fails:
        print("MODEL_REGISTRY FAILED:")
        for f in fails:
            print(f"  🔴 {f}")
        return 1

    print(f"MODEL_REGISTRY OK: {len(rows)} models ({len(live)} live, {len(indev)} in-dev)")
    print("  all live tokens resolved in code")
    print("  no undocumented token collisions (bare ENGINE-B documented + guarded)")
    print("  frankenstein pair intact (rows 5,6)")
    print("  no live-count columns present")

    extra = sorted(code_tokens - {r["stored_token"].strip() for r in rows})
    if extra:
        print(f"\n  ⚠️ live in code, not a model in this roster: {', '.join(extra)}")
        print("     (reported, not silently dropped — see docs/catalog/"
              "model_registry_notes.md)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
