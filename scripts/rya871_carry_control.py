#!/usr/bin/env python3
"""
RYA-871 — is carrying `ep_eV` ADDITIVE? A same-inputs control on the EW artifact.
=================================================================================
    python3 scripts/rya871_carry_control.py --new <dir> --banked <dir> [--stem ...]

The claim this ticket rests on is that carrying the excitation potential changes line
IDENTITY and nothing else — no EW, no abundance, no bar. That claim is cheap to make and
the repo has paid for it before (RYA-848: a banked artifact also carries months of input
drift, so a code fix must be proved with a SAME-INPUTS control).

So the band is re-measured on the same spectrum with the same line selection, and the new
artifact is compared to the banked one COLUMN BY COLUMN:

* the row set must be identical, matched on wavelength — not merely the same LENGTH,
  because two 152-row files can disagree about which 152 lines they are;
* every column the banked file already had must be byte-identical, `ew_method` prose
  included;
* the only difference permitted is the ADDITION of `ep_eV`, which must be non-null on
  every row — a blank identity key is not a key.

A control that only checked "same number of rows and same mean EW" would pass on a file
that measured different lines to the same average, so it is done per row and per column.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "data" / "results" / "rya871"

#: The column this ticket adds. Anything else appearing or disappearing is a finding.
NEW_COLUMNS = {"ep_eV"}


def compare(new_csv: Path, banked_csv: Path) -> dict:
    a = pd.read_csv(banked_csv)
    b = pd.read_csv(new_csv)
    res: dict = {"banked": str(banked_csv), "new": str(new_csv),
                 "n_banked": len(a), "n_new": len(b)}

    added = set(b.columns) - set(a.columns)
    removed = set(a.columns) - set(b.columns)
    res["columns_added"] = sorted(added)
    res["columns_removed"] = sorted(removed)
    res["columns_as_expected"] = bool(added == NEW_COLUMNS and not removed)

    # Row identity by WAVELENGTH, not by position: a reordered file would compare equal
    # row-by-row only by luck, and unequal only by ordering.
    ka = a.wavelength_air_A.round(4)
    kb = b.wavelength_air_A.round(4)
    res["same_line_set"] = bool(sorted(ka) == sorted(kb))
    res["only_in_banked"] = sorted(set(ka) - set(kb))
    res["only_in_new"] = sorted(set(kb) - set(ka))
    if not res["same_line_set"]:
        res["verdict"] = "DIFFERENT LINES — not a same-inputs control; nothing else compared"
        return res

    a = a.set_index(ka).sort_index()
    b = b.set_index(kb).sort_index()
    drift = {}
    for c in a.columns:
        x, y = a[c], b[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            same = ((x.isna() & y.isna())
                    | np.isclose(x.astype(float), y.astype(float),
                                 rtol=0, atol=0, equal_nan=True))
        else:
            same = (x.fillna("").astype(str) == y.fillna("").astype(str))
        n_diff = int((~same).sum())
        if n_diff:
            ex = a.index[~same][:3].tolist()
            drift[c] = {"n_rows_differing": n_diff, "example_wavelengths": ex}
    res["columns_that_drifted"] = drift
    res["carry_is_additive"] = bool(not drift and res["columns_as_expected"])

    ep = b["ep_eV"] if "ep_eV" in b.columns else pd.Series(dtype=float)
    res["ep_present"] = "ep_eV" in b.columns
    res["ep_non_null"] = int(ep.notna().sum())
    res["ep_null"] = int(ep.isna().sum())
    res["ep_range_eV"] = ([round(float(ep.min()), 4), round(float(ep.max()), 4)]
                          if len(ep) and ep.notna().any() else None)
    res["verdict"] = ("ADDITIVE — identity carried, nothing else moved"
                      if res["carry_is_additive"] else "NOT ADDITIVE — see the drift")
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", required=True, type=Path,
                    help="directory holding the freshly measured *_ew.csv")
    ap.add_argument("--banked", required=True, type=Path,
                    help="directory holding the artifact as published")
    ap.add_argument("--stem", default="*_ew.csv")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    results = []
    for new_csv in sorted(a.new.glob(a.stem)):
        banked = a.banked / new_csv.name
        if not banked.exists():
            print(f"  (no banked twin) {new_csv.name}")
            continue
        results.append(compare(new_csv, banked))

    if not results:
        raise SystemExit(f"nothing to compare: {a.new}/{a.stem} vs {a.banked}")

    print("=== same-inputs control: is carrying ep_eV ADDITIVE? ===")
    for r in results:
        name = Path(r["new"]).name
        print(f"\n  {name}")
        print(f"    rows            {r['n_banked']} banked / {r['n_new']} new   "
              f"same line set: {r['same_line_set']}")
        print(f"    columns added   {r['columns_added']}   removed {r['columns_removed']}")
        if r.get("columns_that_drifted"):
            for c, d in r["columns_that_drifted"].items():
                print(f"    🔴 DRIFTED      {c}: {d['n_rows_differing']} rows, "
                      f"e.g. {d['example_wavelengths']}")
        else:
            print(f"    every pre-existing column BYTE-IDENTICAL")
        print(f"    ep_eV           {r['ep_non_null']} non-null / {r['ep_null']} null   "
              f"range {r['ep_range_eV']} eV")
        print(f"    => {r['verdict']}")

    ok = all(r.get("carry_is_additive") and r.get("ep_null") == 0 for r in results)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "rya871_carry_control.json").write_text(json.dumps(
        {"ticket": "RYA-871", "additive_everywhere": ok, "files": results}, indent=2) + "\n")
    print(f"\n  {'ALL ADDITIVE' if ok else '⚠️ NOT ADDITIVE'} — wrote "
          f"{(a.out / 'rya871_carry_control.json').relative_to(ROOT)}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
