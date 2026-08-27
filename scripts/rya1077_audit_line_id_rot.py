#!/usr/bin/env python3
"""RYA-1077: which committed artifacts cite a canonical_gf id that has since moved?

🔴 THE CHECK NOBODY HAD. `canonical_gf.line_id` is positional, so a block removal or
replacement shifts every id after it. Only ONE artifact family had a reproducibility test
(`test_chiappino_digest_rya1059`), which is why 1,739 wrong references accumulated in
silence -- the test is the smoke alarm, not the fire.

This resolves every committed row that carries BOTH a `gf_` id AND a wavelength, and asks
whether the id still points where the row says it does. It takes seconds and it is the
measurement the defect was invisible without.

⚠️ IT COMPARES AGAINST THE WAVELENGTH RECORDED BESIDE THE ID, not against a guess. That is
what makes it decidable: the artifact carries its own witness, so a rotted id is provable
rather than suspected.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
WAVE_COLS = ("wavelength_A", "wavelength_air_A", "wave_air_A", "lambda_air_A")
TOL_A = 0.05


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(ROOT / "data"))
    ap.add_argument("--max-bytes", type=int, default=8_000_000)
    ap.add_argument("--out", default=None)
    ap.add_argument("--fail-on-rot", action="store_true",
                    help="exit non-zero if any id has moved -- for CI")
    a = ap.parse_args()

    cg = pd.read_csv(CANON, low_memory=False)
    by = dict(zip(cg.line_id.astype(str),
                  zip(cg.wavelength_air_A, cg.species.astype(str))))
    per, tot, bad, miss = {}, 0, 0, 0
    for f in sorted(glob.glob(os.path.join(a.root, "**", "*.csv"), recursive=True)):
        if os.path.getsize(f) > a.max_bytes:
            continue
        try:
            d = pd.read_csv(f, low_memory=False)
        except Exception:
            continue
        idc = [c for c in d.columns if "line_id" in c]
        wc = [c for c in d.columns if c in WAVE_COLS]
        if not idc or not wc:
            continue
        n = nb = nm = 0
        for _, r in d.iterrows():
            lid = str(r[idc[0]])
            if not lid.startswith("gf_"):
                continue
            if lid not in by:
                nm += 1
                continue
            try:
                w_rec = float(r[wc[0]])
            except (TypeError, ValueError):
                continue
            n += 1
            if abs(by[lid][0] - w_rec) > TOL_A:
                nb += 1
        if n or nm:
            per[os.path.relpath(f, ROOT)] = {"rows": n, "moved": nb, "id_absent": nm}
            tot += n; bad += nb; miss += nm

    doc = {"ticket": "RYA-1077",
           "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "tolerance_A": TOL_A,
           "totals": {"rows_checked": tot, "id_moved": bad, "id_absent": miss},
           "per_file": per}
    for f, s in sorted(per.items(), key=lambda kv: -kv[1]["moved"]):
        if s["moved"] or s["id_absent"]:
            print(f"{s['moved']:5d}/{s['rows']:5d} moved  {s['id_absent']:3d} absent  {f}")
    print(f"\nTOTAL {bad}/{tot} ids moved, {miss} absent")
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        Path(a.out).write_text(json.dumps(doc, indent=2) + "\n")
        print(f"wrote {a.out}")
    if a.fail_on_rot and (bad or miss):
        raise SystemExit(
            f"{bad} committed artifact row(s) cite a canonical_gf id that no longer points "
            f"at the line recorded beside it. A positional id is not an identity "
            f"(RYA-1077); re-resolve by PHYSICAL KEY and re-stamp.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
