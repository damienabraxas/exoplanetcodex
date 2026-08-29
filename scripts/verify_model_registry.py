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

    # (e) row count. 9 since RYA-1115 added Frankenstein's Dog: 7 live + the Bride
    # (in-dev) + the Dog (not-emitted). The Dog is NOT live — nothing emits it — so it is
    # deliberately outside the (a) token check below rather than exempted from it.
    live = [r for r in rows if r["status"] == "live"]
    indev = [r for r in rows if r["status"] == "in-dev"]
    notemit = [r for r in rows if r["status"] == "not-emitted"]
    if len(rows) != 9 or len(live) != 7 or len(indev) != 1 or len(notemit) != 1:
        fails.append(f"(e) expected 9 rows (7 live + 1 in-dev + 1 not-emitted), got "
                     f"{len(rows)} ({len(live)} live, {len(indev)} in-dev, "
                     f"{len(notemit)} not-emitted)")
    for r in rows:
        if r["status"] not in model_registry.STATUSES:
            fails.append(f"(e) model {r['model_id']} status {r['status']!r} is not in "
                         f"{model_registry.STATUSES}")

    # (f) RYA-1115 — the line_set hook exists and holds only vocabulary values. Opening
    # the column with a guard already on it is what stops RYA-1111 typing a pool name
    # from memory on the day it first populates this.
    if "line_set" not in (rows[0].keys() if rows else ()):
        fails.append("(f) no line_set column — RYA-1111 has no hook to wire into")
    else:
        fails.extend(f"(f) {m}" for m in model_registry.check_line_sets(rows))

    # (g) a non-live row must NOT carry a stored_token. The Bride and the Dog are both in
    # the roster precisely because they do not exist in code yet; a token on either would
    # be an invented one, which is the exact failure (a) guards against on the live side.
    for r in indev + notemit:
        if r["stored_token"].strip():
            fails.append(f"(g) model {r['model_id']} is {r['status']} but carries "
                         f"stored_token {r['stored_token']!r} — a non-live row with a "
                         f"token means it was invented, not read from code")

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

    print(f"MODEL_REGISTRY OK: {len(rows)} models ({len(live)} live, {len(indev)} in-dev, "
          f"{len(notemit)} not-emitted)")
    print("  all live tokens resolved in code")
    print("  no undocumented token collisions (bare ENGINE-B documented + guarded)")
    print("  frankenstein pair intact (rows 5,6)")
    print("  no live-count columns present")
    print("  line_set hook present, all values in vocabulary (RYA-1111 wires them)")
    print("  no non-live row carries an invented token (models 8, 9 blank)")

    extra = sorted(code_tokens - {r["stored_token"].strip() for r in rows})
    if extra:
        print(f"\n  ⚠️ live in code, not a model in this roster: {', '.join(extra)}")
        print("     (reported, not silently dropped — see docs/catalog/"
              "model_registry_notes.md)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
