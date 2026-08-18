#!/usr/bin/env python3
"""
RYA-877 item 4 — what does dispositioning Fe II 5991.371 cost the published number?

    python3 scripts/rya877_disposition_impact.py

READS ONLY. The delta is computed from the per-line abundances of the last successful
Fe II VIS derivation, which are committed beside this script.

🔴 WHY NOT A FRESH SAME-INPUTS PAIR, WHICH IS WHAT THE TICKET ASKS FOR
The intended control — one tree, one commit, the registry row the only difference — was
built and launched (`rya877_impact.sh`) and the deriver REFUSED the band:

    the EW artifact for Fe II carries no `ep_eV` column, so its lines cannot be identified
    by anything but their wavelength (RYA-871). It predates RYA-871; re-run
    scripts/measure_band_profilefit.py for this band.

That is the same RYA-871 EW-table migration already blocking RYA-874 and 8 of the 9 cells
of the RYA-847 sweep restoration. It is not a defect in this disposition.

WHAT IS AND IS NOT ESTABLISHED HERE
  * The DELTA is EXACT. Removing one line from a fixed set of per-line abundances is
    arithmetic over that set, not a re-fit: the same numbers, one fewer of them.
  * The ABSOLUTE values are that run's, not today's. They come from the pre-RYA-871
    derivation, so they carry whatever input drift has happened since. The question 876
    asks is "does the number move", and that is answered exactly; "what is the number"
    is owed a fresh run once the EW tables migrate.
"""
from __future__ import annotations

import csv
import json
import math
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "results" / "rya877"
LINE = 5991.371
TOL_A = 0.05
#: the RYA-852 ionization arbiter trio. Asserted absent here, never assumed.
ARBITER_A = (6147.734, 6238.386, 6247.557)


def _pool(treatment: str):
    f = OUT / f"FeII_3800_6910_kpno_solar_atlas_PROFILEFIT_{treatment}_lines.csv"
    with open(f, newline="") as fh:
        rows = [r for r in csv.DictReader(fh)
                if r["in_aggregate"] == "True" and r["abundance"] not in ("", "nan")]
    return [(float(r["wavelength_air_A"]), float(r["abundance"])) for r in rows]


def _stats(vals):
    n = len(vals)
    sc = st.stdev(vals) if n > 1 else float("nan")
    return {"n": n, "median": st.median(vals), "scatter": sc,
            "stat_dex": sc / math.sqrt(n) if n > 1 else float("nan")}


def main() -> int:
    report = {"ticket": "RYA-877", "line_air_A": LINE, "products": {}}
    for treat in ("1D-LTE", "ENGINE-A", "ENGINE-B"):
        pool = _pool(treat)
        held = [a for w, a in pool if abs(w - LINE) < TOL_A]
        before = _stats([a for _, a in pool])
        rec = {"carries_the_line": bool(held), "before": before}
        if held:
            after = _stats([a for w, a in pool if abs(w - LINE) >= TOL_A])
            rec.update({
                "line_own_abundance": held[0], "after": after,
                "delta": {"value": after["median"] - before["median"],
                          "scatter": after["scatter"] - before["scatter"],
                          "stat_dex": after["stat_dex"] - before["stat_dex"]}})
        report["products"][treat] = rec
        head = f"{treat:<10} n={before['n']:<3}"
        if held:
            d = rec["delta"]
            print(f"{head} {before['median']:.4f} -> {rec['after']['median']:.4f}  "
                  f"value {d['value']:+.4f}   bar(stat) "
                  f"{before['stat_dex']:.4f} -> {rec['after']['stat_dex']:.4f} "
                  f"({d['stat_dex']:+.4f})   [the line itself sat at "
                  f"{rec['line_own_abundance']:.4f}]")
        else:
            print(f"{head} {before['median']:.4f}   does NOT carry the line — unchanged")

    # ── the arbiter, asserted ─────────────────────────────────────────────────────
    arb = _pool("ENGINE-A")
    arb_waves = sorted(w for w, _ in arb)
    carries = [w for w in arb_waves if abs(w - LINE) < TOL_A]
    report["ionization_arbiter"] = {
        "treatment": "ENGINE-A", "lines_air_A": arb_waves,
        "expected_trio": list(ARBITER_A), "carries_the_dispositioned_line": bool(carries),
        "value_unchanged": not carries}
    print(f"\nionization arbiter (ENGINE-A, n={len(arb_waves)}): "
          f"{', '.join(f'{w:.3f}' for w in arb_waves)}")
    if carries:
        raise SystemExit("🔴 the arbiter carries the dispositioned line — STOP, the "
                         "disposition would move the ionization gate")
    print("  does NOT contain 5991.371 => the arbiter value is UNCHANGED (asserted)")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "rya877_disposition_impact.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\n[out] {OUT}/rya877_disposition_impact.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
