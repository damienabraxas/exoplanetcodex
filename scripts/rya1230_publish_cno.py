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
    return rc_all


if __name__ == "__main__":
    sys.exit(main())
