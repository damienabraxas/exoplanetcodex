#!/usr/bin/env python3
"""RYA-1250: re-stamp the run ledger after the cull registry's fingerprint was scoped.

`run_matrix.input_fingerprints` used to hash data/registry/problem_children.csv WHOLE, so one
culled line (Mg I 24566.30) made every cell of every element stale. It now hashes only the
rows a cell's window can see (`registry_fingerprint`). That changes every cell's inputs hash
once, although no cell's inputs changed.

This script proves, per ledger entry, that the cell was built against the registry exactly as
it stands now (its recorded whole-file digest == today's), then swaps that one input row for
the scoped row and recomputes the hash the way `inputs_hash` does. An entry recorded against
an older registry is left alone: it was already stale and re-runs. Run it BEFORE editing the
registry.

With --element, the element's budget record (`units.prev.txt`, whole-element) is converted to
the per-unit `units.prev.json`, so legs built against the old hashes are kept.

    python scripts/rya1250_rescope_registry_fingerprint.py --star solar --element Mg          # report
    python scripts/rya1250_rescope_registry_fingerprint.py --star solar --element Mg --apply
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import element_publish as ep  # noqa: E402
from pipeline import run_matrix as rm  # noqa: E402

REG = str(rm.PROBLEM_CHILDREN.relative_to(ROOT))
WINDOW = re.compile(r"canonical_gf\.csv\[(-?[0-9]+)-(-?[0-9]+)\]$")


def rescope(entry: dict, whole: str) -> tuple[str | None, str]:
    """(new inputs hash, why) for one ledger entry; None when it is not re-stamped."""
    rows = entry.get("inputs")
    if not rows:
        return None, "no recorded inputs"
    if hashlib.sha256(json.dumps(rows, sort_keys=True).encode()).hexdigest() != entry.get("inputs_hash"):
        return None, "recorded inputs do not reproduce the recorded hash"
    reg = [r for r in rows if r["kind"] == "ledger" and r["name"] == REG]
    if len(reg) != 1:
        return None, "no whole-file registry row"
    if reg[0]["digest"] != whole:
        return None, "built against an older registry (already stale)"
    win = [WINDOW.search(r["name"]) for r in rows if r["kind"] == "ledger_scoped"]
    win = [m for m in win if m]
    if len(win) != 1:
        return None, "no canonical_gf window to scope by"
    lo, hi = float(win[0].group(1)), float(win[0].group(2))
    new = [r for r in rows if r is not reg[0]] + [
        {"kind": "registry_scoped", "name": f"{REG}[{lo:.0f}-{hi:.0f}]",
         "digest": rm.registry_fingerprint(lo, hi)}]
    entry["inputs"] = new
    return hashlib.sha256(json.dumps(new, sort_keys=True).encode()).hexdigest(), ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--star", default="solar")
    ap.add_argument("--element", help="also convert this element's budget leg record")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    ledger = json.loads(rm.LEDGER.read_text())
    whole = rm._file_fingerprint(rm.PROBLEM_CHILDREN)
    moved, why_not = {}, {}
    for k, v in ledger.items():
        old = v.get("inputs_hash")
        new, why = rescope(v, whole)
        if new is None:
            why_not[why] = why_not.get(why, 0) + 1
            continue
        v["inputs_hash"] = new
        moved[k] = (old, new)
    print(f"re-stamped {len(moved)} of {len(ledger)} ledger entries")
    for why, n in sorted(why_not.items()):
        print(f"  left alone {n}: {why}")

    sigs = None
    if a.element:
        work = ep.budget_dir(a.star, a.element)
        txt, report = work / "units.prev.txt", rm.ROOT / "data/results/orchestrator" / f"{a.star}_{a.element}_latest.json"
        if txt.exists() and report.exists():
            body = txt.read_text().splitlines()
            cmds = [ln for ln in body if not ln.startswith("# ")]
            built = dict(ln[2:].rsplit(" ", 1) for ln in body if ln.startswith("# "))
            key_of = {u[1]: u[3] for u in ep.units_for(a.star, a.element,
                                                         report=json.loads(report.read_text()),
                                                         interpreter=sys.executable)}
            sigs, kept = [], 0
            for c in cmds:
                k = key_of.get(c)
                old, new = moved.get(k, (None, None))
                ok = k is not None and new is not None and built.get(k) == old
                sigs.append(ep.unit_signature(c, new if ok else "UNPROVEN"))
                kept += ok
            print(f"budget {work.name}: {kept} of {len(cmds)} units' legs carried to the new hashes")
        else:
            print(f"budget {work.name}: no units.prev.txt or run report -- nothing to convert")
    if not a.apply:
        print("(report only; --apply writes)")
        return 0
    rm.LEDGER.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if sigs is not None:
        (ep.budget_dir(a.star, a.element) / "units.prev.json").write_text(json.dumps(sigs, indent=0) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
