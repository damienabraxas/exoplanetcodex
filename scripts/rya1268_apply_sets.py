#!/usr/bin/env python3
"""RYA-1268 -- adopt every remaining element's published-set gf through the orchestrator's own
prepare -> apply (pipeline.element_prepare), WITHOUT culls (culls stay Ryan's yes/no) and
WITHOUT any adoption that would move a line by more than HOLD_MOVE_DEX ("super bad dex should
not be included", Ryan 2026-10-23) -- those are listed in held_set_adoptions.csv.
Missing published lines are added by apply's own published_line_add plan (clash-free only).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import element_prepare as ep  # noqa: E402

HOLD_MOVE_DEX = 0.5
ELEMENTS = ["Li", "Na", "K", "Ca", "P", "S", "Sc", "Ti", "V", "Cr", "Mn", "Co", "Ni", "Cu", "Zn",
            "Sr", "Y", "Zr", "Ba", "Eu"]


def main() -> None:
    held, summary = [], {}
    for el in (sys.argv[1:] or ELEMENTS):
        try:
            rep = ep.prepare("solar", el)
        except Exception as exc:  # report and move on; one element must not stop the rest
            summary[el] = f"PREPARE FAILED: {type(exc).__name__}: {exc}"
            print(el, summary[el], flush=True)
            continue
        keep = []
        for a in rep["adopt"]:
            d = a.get("d_loggf")
            if d is not None and abs(d) > HOLD_MOVE_DEX:
                held.append({"element": el, **{k: a.get(k) for k in (
                    "set", "species", "wavelength_A", "canonical_loggf", "authors_loggf", "d_loggf")}})
            else:
                keep.append(a)
        rep["adopt"], rep["cull_candidates"] = keep, []
        done = ep.apply("solar", el, rep)
        summary[el] = done
        print(el, done, flush=True)
    out = ROOT / "data" / "audit" / "rya1268_line_coverage"
    pd.DataFrame(held).to_csv(out / "held_set_adoptions.csv", index=False)
    (out / "apply_sets_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
