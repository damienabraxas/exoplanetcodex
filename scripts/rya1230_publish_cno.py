#!/usr/bin/env python3
"""RYA-1230 -- publish every C/N/O product whose RYA-587 budget VALIDATED, through the
normal publisher and its gates (never around them).

    python3 scripts/rya1230_publish_cno.py --report R.json --stage DIR [--apply]

Reads the assembler's report (`scripts/rya1230_cno_budget.py`), and for each staged
products CSV calls `scripts/publish_product.py` with the same identity flags the nominal
run's product carries (holding, tier ALL, route SYNTH, selector). Without --apply every
call is a --dry-run.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REASON = ("RYA-1230 re-run: molecular opacity in every band + the standing local-continuum "
          "rule; RYA-587 budget from paired legs on this pool (xi, core window, envelope "
          "estimator, model grid, telluric, CN blends)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--stage", type=Path, required=True)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    rows = json.loads(a.report.read_text())["products"]
    done, rc_all = set(), 0
    for r in rows:
        if r.get("skip") or r.get("verdict") or not r.get("prod_stem"):
            continue
        src = a.stage / r["prod_stem"]
        if src in done or not src.exists():
            continue
        done.add(src)
        cmd = [sys.executable, "scripts/publish_product.py", "--from", str(src),
               "--holding", r["holding"], "--tier", "ALL", "--route", "SYNTH",
               "--element", r["element"], "--star", "solar", "--reason", REASON]
        if r.get("selector"):
            cmd += ["--selector", r["selector"]]
        if not a.apply:
            cmd.append("--dry-run")
        p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        tail = (p.stdout + p.stderr).strip().splitlines()[-3:]
        print(f"rc={p.returncode} {r['element']} {r['holding'][:26]:26} {str(r.get('selector')):11} "
              f"{r['key_treatment']:16} A={r['A']}\n    " + "\n    ".join(tail))
        rc_all |= p.returncode
    if a.apply and rc_all == 0:
        withdraw_held(rows)
    return rc_all


def withdraw_held(rows) -> None:
    """A live LEGACY row whose key was re-run under RYA-1230 and whose re-run budget HOLDS
    is a pre-correction value (no molecular opacity, no placed continuum). RYA-1097's gated
    precedence would keep it live because the re-run cannot publish -- the gate preserving
    a known defect. Evidence is owed to ASSERT a value, not to WITHDRAW one, so the row moves
    to `archive` with the hold that stopped its successor. Through write_feed (RYA-587 gate)."""
    sys.path.insert(0, str(ROOT / "scripts"))
    sys.path.insert(0, str(ROOT))
    import publish_product as pp
    from publish_product import normalise
    import pandas as pd
    held = [r for r in rows if r.get("verdict") and r.get("prod_stem")]
    by_el: dict = {}
    for r in held:
        by_el.setdefault(r["element"], []).append(r)
    for el, rs in by_el.items():
        path = ROOT / "data/products/solar" / f"{el}.json"
        doc = json.loads(path.read_text())
        keys = {}
        for r in rs:
            stem = r["prod_stem"]
            src = next(iter(Path.home().glob(f"codex/rya1230_runs/legs2/unit{r['unit']:02d}/nominal/{stem}")), None)
            if src is None:
                continue
            row = normalise(pd.read_csv(src), holding=r["holding"], tier="ALL", route="SYNTH",
                            selector=r.get("selector"))[0]
            row["star"] = "solar"
            keys[pp.key_of(row)] = r
        now, moved = pp._now(), []
        keep = []
        for prod in doc["products"]:
            k = pp.key_of(prod)
            if k in keys and "uncertainty" not in prod:
                old = dict(prod)
                old["superseded_at"] = now
                old["superseded_reason"] = (
                    "RYA-1230 withdrawal: pre-correction value (no molecular opacity, no "
                    "placed continuum). Its RYA-1230 re-run HOLDS under RYA-587 -- "
                    f"{keys[k]['verdict']}; {'; '.join(keys[k].get('notes') or [])} -- so no "
                    "corrected value replaces it, and a known-biased value is not kept live "
                    "by the gate that refused its successor.")
                doc.setdefault("archive", []).append(old)
                moved.append(k)
            else:
                keep.append(prod)
        if moved:
            doc["products"] = keep
            doc["version"] = pp.bump(doc["version"])
            doc["updated_at"] = now
            pp.write_feed(path, doc)
            print(f"WITHDRAWN to archive in {el}.json v{doc['version']}:")
            for k in moved:
                print(f"    - {k}")


if __name__ == "__main__":
    sys.exit(main())
