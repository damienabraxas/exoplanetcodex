#!/usr/bin/env python3
"""RYA-728: re-stamp ledger entries after the model registry's fingerprint was scoped.

`run_matrix.input_fingerprints` used to hash data/catalog/model_registry.csv WHOLE, so one
edit to the in-dev Bride row (model 8, RYA-1278) made every built cell of every element stale.
It now hashes only the rows that carry a stored_token (`model_registry_fingerprint`).

Unlike RYA-1250, the cells here were built against the registry BEFORE that edit, so the
recorded whole-file digest is not today's. Per entry, this script finds the registry version
whose bytes give the recorded digest in git history, scopes THAT version, and re-stamps the
entry only when its emitting rows are identical to today's. Proven by construction:
  1. the recorded inputs reproduce the recorded hash;
  2. exactly one whole-file model_registry row is present;
  3. a git revision of the file has exactly that sha256;
  4. that revision's emitting rows fingerprint the same as today's.
Every other input is left as recorded, so a cell whose gf, lines, spectrum or code changed
stays stale and re-runs. Only the named element's entries are touched.

    python scripts/rya728_rescope_model_registry_fingerprint.py --element Cr          # report
    python scripts/rya728_rescope_model_registry_fingerprint.py --element Cr --apply
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import run_matrix as rm  # noqa: E402

REG = str(rm.MODEL_REGISTRY.relative_to(ROOT))
SCOPED_NAME = f"{REG}[emitting]"


def history_by_digest() -> dict[str, str]:
    """{sha256 of the file's bytes: its text} for every committed version of the registry."""
    revs = subprocess.run(["git", "log", "--all", "--format=%H", "--", REG], cwd=ROOT,
                          capture_output=True, text=True, check=True).stdout.split()
    out = {}
    for r in revs:
        b = subprocess.run(["git", "show", f"{r}:{REG}"], cwd=ROOT, capture_output=True).stdout
        if b:
            out.setdefault(hashlib.sha256(b).hexdigest(), b.decode("utf-8"))
    return out


def rescope(entry: dict, versions: dict[str, str], now: str) -> tuple[str | None, str]:
    """(new inputs hash, why) for one ledger entry; None when it is not re-stamped."""
    rows = entry.get("inputs")
    if not rows:
        return None, "no recorded inputs"
    if hashlib.sha256(json.dumps(rows, sort_keys=True).encode()).hexdigest() != entry.get("inputs_hash"):
        return None, "recorded inputs do not reproduce the recorded hash"
    reg = [r for r in rows if r["kind"] == "ledger" and r["name"] == REG]
    if len(reg) != 1:
        return None, "no whole-file model registry row"
    text = versions.get(reg[0]["digest"])
    if text is None:
        return None, "recorded registry digest matches no committed version"
    if rm.model_registry_fingerprint(text=text) != now:
        return None, "emitting rows differ from today's (genuinely stale)"
    new = [r for r in rows if r is not reg[0]] + [
        {"kind": "ledger_scoped", "name": SCOPED_NAME, "digest": now}]
    entry["inputs"] = new
    return hashlib.sha256(json.dumps(new, sort_keys=True).encode()).hexdigest(), ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--star", default="solar")
    ap.add_argument("--element", required=True)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    ledger = json.loads(rm.LEDGER.read_text())
    versions, now = history_by_digest(), rm.model_registry_fingerprint()
    prefix = f"{a.star}|{a.element}|"
    mine = [k for k in ledger if k.startswith(prefix)]
    moved, why_not = 0, {}
    for k in mine:
        new, why = rescope(ledger[k], versions, now)
        if new is None:
            why_not[why] = why_not.get(why, 0) + 1
            continue
        ledger[k]["inputs_hash"] = new
        moved += 1
    print(f"re-stamped {moved} of {len(mine)} {a.element} ledger entries "
          f"({len(ledger) - len(mine)} other entries untouched)")
    for why, n in sorted(why_not.items()):
        print(f"  left alone {n}: {why}")
    if not a.apply:
        print("(report only; --apply writes)")
        return 0
    rm.LEDGER.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
