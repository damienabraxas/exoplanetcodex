#!/usr/bin/env python3
"""RYA-1315 (Ryan, 2026-11-03): put Asplund+2021's OWN atmosphere terms on every live product.

The published sigma is Asplund, Amarsi & Grevesse 2021 Sect. 2.1: statistical (+) half the
mean (<3D> - HM74) (+) half the mean (3D - <3D>) (+) NLTE. Live products carried a
MARCS-vs-ATLAS9 leg (some with a 3D-NLTE minus 1D-NLTE term added) in `model_atmosphere`
instead. This replaces that ONE component with the Scott+2015 I/II per-line terms on the
product's own pool (`uncertainty_indicator_ids`, the lines its budget was built on), re-assembles the budget through
the RYA-587 contract, and rewrites the feed through `publish_product.write_feed` (which
re-stamps sigma_asplund and runs the publication gate). No abundance changes; nothing is run.

    python scripts/rya1315_asplund_restamp.py --element Na [--write]

A product whose species Scott does not tabulate (C/N/O), or whose pool ids are not all in canonical_gf,
is left exactly as it is and listed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT)]

from pipeline.uncertainty_contract import assemble  # noqa: E402
from rya1230_cno_budget import scott_atmosphere  # noqa: E402

_GF = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False,
                  usecols=["physical_id", "wavelength_air_A"])
SOURCE = ("Asplund+2021 Sect. 2.1: half the mean (<3D> - Holweger & Mueller 1974) (+) half the "
          "mean (3D - <3D>), per-line LTE abundances of Scott+2015 I/II")


def restamp(product: dict, element: str) -> str | None:
    """Replace the product's model_atmosphere in place; return why it was skipped, or None."""
    u = product.get("uncertainty")
    if not u:
        return "legacy row (no budget)"
    ids = set(product.get("uncertainty_indicator_ids") or [])
    pool = _GF[_GF.physical_id.isin(ids)].drop_duplicates("physical_id")
    if not ids or len(pool) != len(ids):
        return f"pool ids not all in canonical_gf ({len(pool)}/{len(ids)})"
    scott = scott_atmosphere(element, str(product.get("ion") or "").strip(),
                             pool["wavelength_air_A"].astype(float).round(3))
    if scott is None:
        return "species not in Scott+2015"
    comps = [dict(c) for c in u["components"] if c["name"] != "model_atmosphere"]
    comps.append(dict(name="model_atmosphere",
                      sigma_dex=float(np.hypot(scott["temperature_structure_dex"],
                                               scott["inhomogeneity_dex"])),
                      state="MEASURED",
                      source=f"{SOURCE} ({scott['basis']}, n={scott['n_scott_lines']})",
                      evidence={"pool_sha256": u["scope"]["pool_sha256"], **scott}))
    numeric = [c for c in comps if c["state"] in {"MEASURED", "DEFINED"}]
    doc = assemble(u["scope"], comps,
                   covariance=np.diag([c["sigma_dex"] ** 2 for c in numeric]).tolist(),
                   covariance_source=u["covariance_source"], assumptions=u["assumptions"],
                   required_extra=tuple(u.get("required_extra") or ()))
    product["uncertainty"] = doc
    product["sigma_reported"] = doc["sigma_reported"]
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--element", required=True)
    ap.add_argument("--feed", type=Path, default=None)
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    from publish_product import write_feed
    feed = a.feed or ROOT / f"data/products/solar/{a.element}.json"
    before = json.loads(feed.read_text())
    doc = json.loads(feed.read_text())
    skipped, moved = {}, []
    for i, p in enumerate(doc["products"]):
        why = restamp(p, a.element)
        if why:
            skipped[why] = skipped.get(why, 0) + 1
        else:
            moved.append(i)
    if a.write:
        write_feed(feed, doc)
        doc = json.loads(feed.read_text())
    print(f"{a.element}: {len(moved)} re-stamped, skipped {skipped}")
    for i in moved:
        o, p = before["products"][i], doc["products"][i]
        s_old, s_new = o.get("sigma_asplund"), p.get("sigma_asplund")
        print(f"  {p.get('ion'):3} {p.get('band'):12} {str(p.get('holding'))[:32]:32} "
              f"{str(p.get('selector'))[:16]:16} {str(p.get('treatment'))[:20]:20} A={p.get('A')}  "
              f"sigma_asplund {None if s_old is None else round(s_old, 3)} -> "
              f"{None if s_new is None else round(s_new, 3)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
