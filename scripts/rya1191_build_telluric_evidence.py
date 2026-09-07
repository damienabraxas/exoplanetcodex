#!/usr/bin/env python3
"""RYA-1191 (D) — turn the measured telluric verdicts into the table the PIPELINE reads.

    python3 scripts/rya1191_build_telluric_evidence.py

Reads `data/results/rya1191/rya1191_telluric_residual.json` and writes
`data/catalog/telluric_correction_evidence.csv`. Applies no correction and moves no value;
it only makes an existing MEASUREMENT reachable by the code that has to act on it.

🔴 WHY A PER-BAND TABLE, WHEN RYA-1194 ALREADY SHIPPED A PER-HOLDING ONE
------------------------------------------------------------------------
`telluric_policy.VERIFIED_HOLDING_STATE` records one state per holding -- "corrected" or
"raw" -- and its own docstring names this ticket as the thing that would flip the Kitt
Peak entry. It cannot be flipped, because the answer is not one state:

    solar_kpno_kurucz2005_corrected   CLEAN in H2O 7160-7340 and O2 A-band
    solar_kpno_molecfit_corrected     PARTIAL residual in O2 gamma AND in H2O 9280-9600
    solar_harps_molecfit_corrected    essentially UNCORRECTED in O2 gamma 6270-6300

A holding is corrected in some bands and not in others, and collapsing that to one word
loses exactly the distinction a line-selection quarantine needs. So this table is keyed on
(holding, band) and the per-holding map keeps its job as the default.

⚠️ GENERATED, NEVER HAND-TYPED. A hand-maintained copy of a measurement drifts from the
measurement (RYA-686), and the stale side is the one that passes. `--check` re-derives and
diffs instead of writing, so CI can fail on drift rather than discover it later.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/results/rya1191/rya1191_telluric_residual.json"
OUT = ROOT / "data/catalog/telluric_correction_evidence.csv"

FIELDS = ["holding_id", "band_name", "lo_A", "hi_A", "state", "judged_by",
          "statistic", "value", "reference_value", "ticket", "note"]

#: A holding is CLEAN in a band only when the test that judged it had power AND the
#: holding did not trip it. Everything else is recorded with the reason it is not clean --
#: `partial`, `uncorrected`, or `undetermined` -- because the three are different and a
#: quarantine should behave differently for each.
CLEAN, PARTIAL, UNCORRECTED, UNDETERMINED = "clean", "partial", "uncorrected", "undetermined"


def rows_from(doc: dict) -> list[dict]:
    out = []
    for b in doc["bands"]:
        if b.get("state") != "MEASURED":
            continue
        judged = b.get("judged_by")
        gaps = {g.split(":")[0]: g for g in b.get("gaps", [])}
        refs = [c for c in b["holdings"].values()
                if c.get("note") == "UNCORRECTED REFERENCE" and c["state"] == "MEASURED"]
        for hold, c in b["holdings"].items():
            if c.get("note") == "UNCORRECTED REFERENCE" or c.get("state") != "MEASURED":
                continue
            if judged is None:
                state, stat, val, ref = UNDETERMINED, "none", "", ""
            elif judged.startswith("TEMPLATE"):
                stat = "residual_z"
                val = c.get("residual_z")
                ref = max((x.get("residual_z", float("-inf")) for x in refs), default="")
                if hold in gaps:
                    state = (UNCORRECTED if "essentially UNCORRECTED" in gaps[hold]
                             else PARTIAL)
                else:
                    state = CLEAN
            else:
                stat = "depth_over_own_side_windows"
                val = c.get("depth_over_side")
                ref = max((x.get("depth_over_side") or 0 for x in refs), default="")
                state = PARTIAL if hold in gaps else CLEAN
            out.append({
                "holding_id": hold, "band_name": b["band"],
                "lo_A": f"{b['lo_A']:.1f}", "hi_A": f"{b['hi_A']:.1f}",
                "state": state, "judged_by": judged or "NO POWER",
                "statistic": stat, "value": ("" if val is None else val),
                "reference_value": ref, "ticket": "RYA-1191",
                "note": (gaps.get(hold, "") or
                         ("clean: the uncorrected references light this band up and this "
                          "holding does not follow them" if state == CLEAN else
                          "no test had power in this band — NOT a clean verdict (RYA-833)")),
            })
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="re-derive and diff against the committed table; write nothing")
    a = ap.parse_args(argv)
    if not SRC.exists():
        raise SystemExit(f"{SRC.relative_to(ROOT)} is absent — run "
                         f"scripts/rya1191_telluric_residual.py first")
    rows = rows_from(json.loads(SRC.read_text()))
    body = ["# RYA-1191 — MEASURED telluric state per holding per band. GENERATED by",
            "# scripts/rya1191_build_telluric_evidence.py from",
            "# data/results/rya1191/rya1191_telluric_residual.json. Do not hand-edit:",
            "# `--check` diffs this against a fresh derivation and CI fails on drift.",
            ",".join(FIELDS)]
    for r in rows:
        body.append(",".join('"%s"' % str(r[f]).replace('"', "'") for f in FIELDS))
    text = "\n".join(body) + "\n"
    if a.check:
        if not OUT.exists():
            print(f"MISSING {OUT.relative_to(ROOT)}"); return 1
        same = OUT.read_text() == text
        print("IN SYNC" if same else "🔴 DRIFTED from the measurement — re-generate")
        return 0 if same else 1
    OUT.write_text(text)
    print(f"wrote {OUT.relative_to(ROOT)} — {len(rows)} rows")
    for r in rows:
        if r["state"] != "clean":
            print(f"  {r['state'].upper():<13} {r['holding_id']:<34}{r['band_name']:<16}"
                  f"{r['lo_A']}-{r['hi_A']}  {r['statistic']}={r['value']} "
                  f"(refs {r['reference_value']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
