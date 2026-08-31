"""
scripts/rya1163_ir_xi.py — RYA-1163
===================================
dA/dxi for the three IR (NIR) Fe pools, as a PER-LINE PAIRED DIFFERENTIAL.

WHY THIS EXISTS. `data/results/rya1120/xi_sigma_reported.json` covers 46 products and
NONE of them is IR — nor near-UV, nor red-optical. The cause is not a skipped arm:
`~/xi_campaign/run_campaign.py:132` hardcodes `--lo 4200 --hi 6910` on EVERY unit, so all
15 pools x 2 xi measured the VIS window. IAG and Kitt Peak were both in that campaign.
The BAND was pinned, not the arm. This closes the IR third of that gap.

🔴 THE DERIVATIVE IS A PAIRED DIFFERENTIAL OR IT IS NOTHING (RYA-1083). Line acceptance
moves with xi, so `median(A_hi) - median(A_lo)` is a different statistic from
`median(A_hi - A_lo)`. Measured here: they disagree on FOUR of the six pools, by up to
2x. `pipeline.paired_differential` is called, never reimplemented.

Inputs: `~/xi_campaign/out_ir/FeI_<tag>_NIR_GRADED_SYNTH_xi{0.90,1.10}/*_lines.csv`,
produced by `~/xi_campaign/run_ir.py` (snapshot-then-move, so no committed artifact is
touched). The runs are Mac-local — a recorded regenerability gap (RYA-1011), not a
silent one.
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pipeline.paired_differential import paired_differential  # noqa: E402

CAMPAIGN = pathlib.Path.home() / "xi_campaign" / "out_ir"
FEED = REPO / "data" / "products" / "solar" / "Fe.json"
OUT = REPO / "data" / "results" / "rya1163" / "ir_xi_dadxi.json"

#: (xi+0.10) - (xi-0.10). Same convention as `rya1120_sigma_reported.XI_SPAN_KMS`.
XI_SPAN_KMS = 0.20
#: The campaign's delta_xi, from `data/results/rya1120/xi_sigma_reported.json`.
DELTA_XI_KMS = 0.2912
#: Below this many paired lines a slope is not publishable. Two points define a line
#: exactly and carry no dispersion, so the "derivative" would be an artifact of which two
#: lines survived -- report it, never publish it (RYA-907: unmeasured is not a value).
MIN_PAIRED = 3

POOLS = {
    "IAG":     ("FeI_iag_NIR_GRADED_SYNTH",    "iag_fts_solar_atlas"),
    "KP":      ("FeI_kpno_NIR_GRADED_SYNTH",   "kpno_solar_atlas"),
    "CRIRES+": ("FeI_crires_NIR_GRADED_SYNTH", "crires_plus"),
}
TREATMENTS = ("1D-LTE", "ENGINE-A")


def _lines(stem: str, xi: str, treat: str) -> pd.DataFrame | None:
    d = CAMPAIGN / f"{stem}_xi{xi}"
    hits = [p for p in d.glob("*_lines.csv") if p.name.endswith(f"_{treat}_lines.csv")]
    return pd.read_csv(hits[0]) if hits else None


def main() -> int:
    feed = json.loads(FEED.read_text())
    nir = {(x["instrument"], x["treatment"]): x
           for x in feed["products"] if x.get("band") == "NIR"}
    rows = []
    for pool, (stem, inst) in POOLS.items():
        for treat in TREATMENTS:
            hi, lo = _lines(stem, "1.10", treat), _lines(stem, "0.90", treat)
            if hi is None or lo is None:
                print(f"  {pool} {treat}: no per-line pair — SKIP")
                continue
            r = paired_differential(hi, lo)
            p = nir[(inst, treat)]
            d_a_dxi = r.median / XI_SPAN_KMS
            sigma_xi = abs(d_a_dxi) * DELTA_XI_KMS
            publishable = r.n_paired >= MIN_PAIRED
            syst_new = math.sqrt(p["sigma_syst"] ** 2 + sigma_xi ** 2)
            rows.append({
                "pool": pool, "instrument": inst, "holding": p.get("holding"),
                "band": "NIR", "tier": p.get("tier"), "treatment": treat,
                "A": p["A"], "n_lines": p["n_lines"],
                "n_paired": r.n_paired,
                "median_delta_dex": r.median,
                "dA_dxi": round(d_a_dxi, 4),
                "sigma_xi": round(sigma_xi, 4),
                "xi_state": "MEASURED" if publishable else "UNMEASURED",
                "xi_note": (f"|dA/dxi|={abs(d_a_dxi):.4f} x delta_xi={DELTA_XI_KMS} on "
                            f"{r.n_paired} paired lines")
                           if publishable else
                           (f"REFUSED: only {r.n_paired} paired lines (< {MIN_PAIRED}). Two "
                            f"points define a slope exactly and carry no dispersion; the "
                            f"value would be an artifact of which lines survived."),
                "sigma_stat": p["sigma_stat"],
                "sigma_syst_published": p["sigma_syst"],
                "sigma_syst_with_xi": round(syst_new, 4) if publishable else None,
                "sigma_reported_published": round(
                    math.hypot(p["sigma_stat"], p["sigma_syst"]), 4),
                "sigma_reported_with_xi": round(
                    math.hypot(p["sigma_stat"], syst_new), 4) if publishable else None,
                # RYA-1083's trap, carried so it is visible in the artifact.
                "difference_of_aggregates": r.difference_of_aggregates,
                "collision": r.collision,
            })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "ticket": "RYA-1163", "read_only": True, "star": "solar", "species": "Fe I",
        "band": "NIR", "xi_span_kms": XI_SPAN_KMS, "delta_xi_kms": DELTA_XI_KMS,
        "min_paired": MIN_PAIRED, "n_pools": len(rows), "pools": rows,
    }, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(REPO)}  ({len(rows)} pools)")
    for r in rows:
        print(f"  {r['pool']:8s} {r['treatment']:9s} n_paired={r['n_paired']:3d} "
              f"dA/dxi={r['dA_dxi']:+.4f} sigma_xi={r['sigma_xi']:.4f} "
              f"{r['xi_state']}{'  COLLISION' if r['collision'] else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
