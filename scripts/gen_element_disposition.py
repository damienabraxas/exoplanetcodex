#!/usr/bin/env python3
"""
Generate the per-element disposition report (RYA-663).

    python scripts/gen_element_disposition.py            # write both artifacts
    python scripts/gen_element_disposition.py --check    # CI mode: fail if stale

``--check`` regenerates in memory and compares against what is committed, the same
contract as the RYA-654 tracker generator: the committed report must equal a fresh
run, so a report that has drifted from the ledger cannot sit in the repo unnoticed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.element_disposition import (  # noqa: E402
    REPO_ROOT,
    build_report,
    render_markdown,
)

JSON_OUT = REPO_ROOT / "data" / "audit" / "element_disposition_rya663.json"
MD_OUT = REPO_ROOT / "docs" / "audit" / "element_disposition_rya663.md"


def _render(report: dict) -> tuple[str, str]:
    # The git provenance of the inputs is part of the report, but it changes on every
    # commit that touches them — including this one. --check compares the SUBSTANCE,
    # so the volatile block is excluded from the comparison, not from the output.
    stable = {k: v for k, v in report.items() if k != "inputs"}
    return json.dumps(report, indent=2, sort_keys=True) + "\n", json.dumps(
        stable, indent=2, sort_keys=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="verify the committed report matches a fresh run; write nothing")
    args = ap.parse_args()

    report = build_report()
    js, stable = _render(report)
    md = render_markdown(report)

    if args.check:
        if not JSON_OUT.exists():
            print(f"MISSING: {JSON_OUT.relative_to(REPO_ROOT)} — run without --check",
                  file=sys.stderr)
            return 1
        _, committed_stable = _render(json.loads(JSON_OUT.read_text(encoding="utf-8")))
        if committed_stable != stable:
            print(f"STALE: {JSON_OUT.relative_to(REPO_ROOT)} does not match a fresh run — "
                  "regenerate it", file=sys.stderr)
            return 1
        print("Element disposition report is up to date.")
        return 0

    JSON_OUT.parent.mkdir(parents=True, exist_ok=True)
    MD_OUT.parent.mkdir(parents=True, exist_ok=True)
    JSON_OUT.write_text(js, encoding="utf-8")
    MD_OUT.write_text(md, encoding="utf-8")
    flip = report["can_flip_now"]
    print(f"Wrote {JSON_OUT.relative_to(REPO_ROOT)} and {MD_OUT.relative_to(REPO_ROOT)}")
    print(f"  can flip now : {', '.join(flip) if flip else 'none'}"
          + ("  (PROVISIONAL — gate 3 read a stale input)" if report["gate3_provisional"]
             and flip else ""))
    print(f"  stale inputs : {len(report['stale_input_evidence'])} contradiction(s)")
    print(f"  value splits : {len(report['value_disagreements'])} element(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
