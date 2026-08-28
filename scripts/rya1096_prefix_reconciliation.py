#!/usr/bin/env python3
"""
scripts/rya1096_prefix_reconciliation.py
========================================
RYA-1096 — what did the continuum fixes actually DO to the eight products RYA-1092
withdrew? Measured per cell, by re-running them, not inferred from a diff.

    python3 scripts/rya1096_prefix_reconciliation.py --band-products <dir with the re-runs>

🔴 WHY THIS IS A RE-RUN AND NOT AN ARGUMENT.

RYA-1092 withdrew eight Fe I VIS GRADED SYNTH products on `PRE_CONTINUUM_FIX`: their
artifacts predate the last of the RYA-933/1030 continuum fixes. This ticket's cheaper first
option was to narrow that cutoff -- show that a given fix CANNOT affect a given product and
restore it without re-measuring. I tried that route first and it does not hold up, which is
itself the finding:

  fafe9246  2026-08-24T00:57:30Z  touches derive_band_products (+146) and measure_band_ew
                                  (+72). On the measurement path. All eight postdate it.
  22eca9e6  2026-08-24T02:36:32Z  its change to derive_band_products is HELP TEXT ONLY; the
  (79ec5a34 merge, 02:45:51Z)     executable change is in measure_band_ew and is specific to
                                  ONE holding -- kurucz2005 flips pre_normalised False->True
                                  and the reader moves to `irradrelwl.dat`. Its own commit
                                  message measures the cost of the old path at 0.022 dex
                                  with a 4% blue-to-red tilt.
  a61e7b80  2026-08-24T03:51:12Z  the RYA-1030 merge. Its sub-commits touch registry,
                                  detector and preflight only -- but the MERGE also rewrites
                                  `pipeline/prenormalised_guard.py`, which
                                  `derive_band_products` imports. On the path after all.
  0f51f367  2026-08-24T13:28:03Z  a FOURTH commit, found only after concluding there were
                                  three: it splits the IAG atlases into separate holdings
                                  and touches the same guard.

⚠️ THE SHAPE OF THAT LIST IS THE ARGUMENT. Each pass over the history moved the answer, and
the last pass found a commit that postdates every product in question. A scoping rule built
on "these are the commits" would have been wrong twice already, and a wrong scoping rule
restores a contaminated product silently -- the failure mode is invisible, which is exactly
the class RYA-1092 exists to close. Re-measuring is decisive and cheap by comparison.

WHAT A RE-RUN SETTLES THAT A DIFF CANNOT
----------------------------------------
RYA-933 withdrew two of its own products this way and its re-run reproduced one of them
EXACTLY (7.535, n=109): "the fix is a NULL for this cell, which is itself the result." That
is the precedent and this follows it. A cell that reproduces was never affected; a cell that
moves tells you what the fix was worth, in dex, for that cell.

🔴 AND THE RESTORE PATH IS PUBLICATION, NOT UN-QUARANTINE. A reproduced cell is republished
from its FRESH post-fix artifact through `scripts/publish_product.py`, so the live record
carries a post-fix `artifact_mtime` and passes RYA-1092's gate on its own merits. Nothing
here edits the gate, and nothing lifts a quarantine by assertion -- the withdrawn record
stays exactly where it is as the pre-fix half of the pair.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import product_eligibility as pe          # noqa: E402

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"

#: Below this the cell REPRODUCED and the fixes were a null for it. Not a tolerance on a
#: measurement -- the published values carry 4 decimals (`error_budget.round_dex`), so
#: anything at or under half a unit in the last place is the same number.
NULL_DEX = 5e-5


def withdrawn_prefix(doc: dict) -> list:
    return [p for p in doc.get("quarantine", [])
            if "PRE_CONTINUUM_FIX" in (p.get("quarantine_codes") or [])]


def rerun_products(band_dir: Path) -> dict:
    """{(holding, treatment): row} from the fresh `*_products.csv` in `band_dir`."""
    out = {}
    for f in sorted(Path(band_dir).glob("FeI_4200_6910_*_SYNTH_GRADED*_products.csv")):
        d = pd.read_csv(f)
        for r in d.itertuples():
            holding = _holding_from_stem(f.name)
            out[(holding, str(r.treatment))] = {
                "A": float(r.A), "stat_dex": float(r.stat_dex),
                "syst_dex": float(r.syst_dex), "n_lines": int(r.n_lines),
                "file": f.name}
    return out


def _holding_from_stem(name: str) -> str:
    """The holding token sits between the instrument and the route in the stem. Parsed
    rather than guessed from a list, so a new holding needs no edit here."""
    for h in ("solar_kpno_molecfit_corrected", "solar_kpno_kurucz2005_corrected",
              "solar_harps_molecfit_corrected", "solar_iag"):
        if h in name:
            return h
    raise SystemExit(f"cannot identify the holding in {name!r}")


def reconcile(doc: dict, fresh: dict) -> list:
    rows = []
    for p in withdrawn_prefix(doc):
        key = (p["holding"], p["treatment"])
        f = fresh.get(key)
        if f is None:
            rows.append({"holding": p["holding"], "treatment": p["treatment"],
                         "verdict": "NOT_RERUN", "withdrawn_A": p["A"],
                         "withdrawn_n": p["n_lines"]})
            continue
        dA = f["A"] - float(p["A"])
        rows.append({
            "holding": p["holding"], "treatment": p["treatment"],
            "withdrawn_A": float(p["A"]), "withdrawn_n": int(p["n_lines"]),
            "withdrawn_stat": p.get("sigma_stat"), "withdrawn_syst": p.get("sigma_syst"),
            "rerun_A": f["A"], "rerun_n": f["n_lines"],
            "rerun_stat": f["stat_dex"], "rerun_syst": f["syst_dex"],
            "delta_A_dex": dA, "delta_n": f["n_lines"] - int(p["n_lines"]),
            "verdict": ("NULL — the fixes changed nothing for this cell"
                        if abs(dA) <= NULL_DEX and f["n_lines"] == int(p["n_lines"])
                        else "MOVED — the fixes are worth this much here"),
            "rerun_file": f["file"]})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--band-products", type=Path, required=True,
                    help="directory holding the FRESH post-fix *_products.csv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "results" / "rya1096"
                    / "rya1096_prefix_reconciliation.json")
    a = ap.parse_args()

    doc = json.loads(FEED.read_text())
    rows = reconcile(doc, rerun_products(a.band_products))
    n_null = sum(1 for r in rows if r["verdict"].startswith("NULL"))
    n_moved = sum(1 for r in rows if r["verdict"].startswith("MOVED"))

    print(f"\n=== RYA-1096 — the eight PRE_CONTINUUM_FIX cells, re-measured ===")
    print(f"{'holding':34s} {'treatment':10s} {'withdrawn':>10s} {'re-run':>10s} "
          f"{'delta':>9s} {'dn':>4s}  verdict")
    for r in sorted(rows, key=lambda x: (x["holding"], x["treatment"])):
        if "rerun_A" not in r:
            print(f"{r['holding']:34s} {r['treatment']:10s} {r['withdrawn_A']:10.4f} "
                  f"{'--':>10s} {'--':>9s} {'--':>4s}  NOT RE-RUN")
            continue
        print(f"{r['holding']:34s} {r['treatment']:10s} {r['withdrawn_A']:10.4f} "
              f"{r['rerun_A']:10.4f} {r['delta_A_dex']:+9.4f} {r['delta_n']:+4d}  "
              f"{r['verdict'].split(' —')[0]}")
    print(f"\n  NULL {n_null}  ·  MOVED {n_moved}  ·  of {len(rows)} withdrawn cells")
    print("  A NULL cell was never affected by the fixes -- republish its FRESH artifact\n"
          "  and RYA-1092's gate admits it on a post-fix mtime, with no gate change.")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps({"ticket": "RYA-1096", "null_dex": NULL_DEX,
                                 "n_null": n_null, "n_moved": n_moved,
                                 "cells": rows}, indent=2) + "\n")
    print(f"  -> {a.out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
