#!/usr/bin/env python3
"""
scripts/rya1094_label_nlte_levels.py
====================================
RYA-1094 — write Gerber NLTE level labels onto a VALD-built iSpec linelist.

WHY
---
Departures are per LEVEL. bsyn maps a line to the model atom by its two level
identifications, and falls back to departure = 1 for any line it cannot identify --
silently, so an unlabelled list yields an LTE spectrum under an NLTE label
(RYA-534 Co/Ni, RYA-764). `pipeline/nearuv_linelist.py` writes `nlte_label_* = 'none'`
for every row it builds, so EVERY VALD-built list in this repo is unlabelled: the
near-UV list, RYA-762's IR 9200-13000 list, and the RYA-1094 H list alike. NLTE has
therefore never been reachable in any of those bands -- it works in VIS only because
the iSpec-vendored GESv6 list ships labels (15,706 Fe lines).

WHAT IT DOES
------------
For each line of `--species`, resolves BOTH level endpoints against the deck's model
atom by (energy, J) using `pipeline.model_atom.resolve_level` (RYA-823), and writes
`nlte`, `nlte_level_low/up` and `nlte_label_low/up` where BOTH resolve UNIQUELY.

🔴 BOTH ENDPOINTS OR NEITHER. A line with one identified level is not half-corrected;
bsyn needs the pair, and writing a single label would present an unusable line as
labelled -- inflating the coverage number the product reports.

⚠️ AMBIGUOUS AND ABSENT ARE LEFT UNLABELLED ON PURPOSE. An energy is a weak proxy for a
level (RYA-763), so a near-miss is not promoted by widening the tolerance. Lines that
stay unlabelled run in LTE, which is why the caller must read the POOL coverage the
product reports and not just the aggregate.

Dry-run by default; --apply writes the file in place.
"""
from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import codex_path                      # noqa: E402
from pipeline.model_atom import read_gerber_atom, resolve_level  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--linelist", required=True)
    ap.add_argument("--atom", default=None, help="model atom (default: the deck's atom.fe607a)")
    ap.add_argument("--species", default="Fe 1")
    ap.add_argument("--tol-eV", type=float, default=0.001)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    atom = Path(a.atom) if a.atom else Path(str(codex_path("grids.gerber_ts"))) / "atom.fe607a"
    lv = read_gerber_atom(atom)
    print(f"model atom : {atom.name} — {len(lv)} levels")

    ll = Path(a.linelist)
    d = pd.read_csv(ll, sep="\t", low_memory=False)
    sel = d["element"].astype(str) == a.species
    print(f"linelist   : {ll.parent.name} — {len(d)} rows, {int(sel.sum())} {a.species}")

    term_col = "term" if "term" in lv.columns else None
    verdicts: collections.Counter = collections.Counter()
    n = 0
    for i in d.index[sel]:
        r = d.loc[i]
        vl, il, _ = resolve_level(lv, energy_eV=float(r.lower_state_eV),
                                  j=float(r.lower_j), tol_eV=a.tol_eV)
        vu, iu, _ = resolve_level(lv, energy_eV=float(r.upper_state_eV),
                                  j=float(r.upper_j), tol_eV=a.tol_eV)
        verdicts[f"{vl}/{vu}"] += 1
        if vl != "UNIQUE" or vu != "UNIQUE":
            continue                       # both or neither — see the docstring
        d.at[i, "nlte"] = "T"
        d.at[i, "nlte_level_low"] = int(il)
        d.at[i, "nlte_level_up"] = int(iu)
        if term_col:
            d.at[i, "nlte_label_low"] = str(lv.loc[lv.index[il - 1], term_col])
            d.at[i, "nlte_label_up"] = str(lv.loc[lv.index[iu - 1], term_col])
        else:
            d.at[i, "nlte_label_low"] = f"lev{int(il)}"
            d.at[i, "nlte_label_up"] = f"lev{int(iu)}"
        n += 1

    print(f"\n  labelled (BOTH endpoints UNIQUE): {n} of {int(sel.sum())} "
          f"({100.0 * n / max(1, int(sel.sum())):.1f}%)")
    print("  verdict pairs (top 6):")
    for k, v in verdicts.most_common(6):
        print(f"    {k:22} {v}")

    if a.apply:
        d.to_csv(ll, sep="\t", index=False)
        print(f"\n  ✅ wrote {ll}")
    else:
        print("\n  (dry-run — re-run with --apply)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
