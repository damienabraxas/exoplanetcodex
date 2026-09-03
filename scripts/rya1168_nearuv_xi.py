"""
scripts/rya1168_nearuv_xi.py — RYA-1168
=======================================
dA/dxi for the four near-UV (3000-3780 A) DEEPGRADED Fe pools, as a PER-LINE PAIRED
DIFFERENTIAL. Sibling of `scripts/rya1163_ir_xi.py`.

WHY THIS EXISTS. RYA-1113 audited the near-UV bars and recorded, in its F7, that "no
stellar-parameter term reaches the published bar ... Fe II has no row in the stellar
budget at all". It was READ-ONLY, so it named the absence and stopped. This is the
missing run.

🔴 THE ARM WAS NEVER SKIPPED; THE BAND WAS. `~/xi_campaign/run_campaign.py:132` hardcodes
`--lo 4200 --hi 6910` on EVERY unit. Both KP DEEPGRADED pools WERE in that campaign --
they measured VIS. The pool key `(ion, holding, tier, instrument, ROUTE)` carries no band,
so `FeI_solar_kpno_molecfit_corrected_DEEPGRADED_SYNTH` names the near-UV product and the
VIS one alike. Nothing cross-attributed only because the harness matches on
(A, n_lines, holding) rather than on that key.

🔴 AND THE PAIRED DIFFERENTIAL IS NOT OPTIONAL HERE. Seven of these eight pools are
RYA-1083 collisions. The worst, Fe II molecfit 1D-LTE, reads +0.0020 by the naive
difference-of-aggregates -- near zero, and the WRONG SIGN -- against a paired median of
-0.0120. A reader subtracting the two published numbers would conclude near-UV Fe II has
no xi sensitivity at all.

⚠️ THE FIRST RUN OF THIS WAS THROWN AWAY, AND THE REASON IS THE ARTIFACT'S PROVENANCE.
It used ONE worktree for four concurrent jobs, so each job's move-out swept up whatever
was new in the shared `data/results/band_products/` -- including a sibling's output. The
kurucz2005 job's directory came back holding the molecfit job's line files. A stolen file
could have come from either xi, so the pairing was unknowable. `run_uv.py` now claims one
of four isolated worktrees per worker. Nothing in a run's own artifacts records which xi
produced it, so ISOLATION IS THE ONLY THING THAT MAKES THE PAIRING TRUE -- there is no
downstream check that would have caught it.
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

CAMPAIGN = pathlib.Path.home() / "xi_campaign" / "out_uv"
FEED = REPO / "data" / "products" / "solar" / "Fe.json"
OUT = REPO / "data" / "results" / "rya1168" / "nearuv_xi_dadxi.json"

XI_SPAN_KMS = 0.20
DELTA_XI_KMS = 0.2912
#: Below this many paired lines a slope is not publishable (RYA-907: unmeasured is not
#: a value). Every near-UV pool clears it; the guard is kept so it stays true by CHECK.
MIN_PAIRED = 3

POOLS = {
    "kur_FeI":  ("I",  "solar_kpno_kurucz2005_corrected"),
    "mol_FeI":  ("I",  "solar_kpno_molecfit_corrected"),
    "kur_FeII": ("II", "solar_kpno_kurucz2005_corrected"),
    "mol_FeII": ("II", "solar_kpno_molecfit_corrected"),
}
TREATMENTS = ("1D-LTE", "ENGINE-A")


def _lines(ion: str, tag: str, xi: str, treat: str) -> pd.DataFrame | None:
    d = CAMPAIGN / f"Fe{ion}_{tag}_nearUV_DEEPGRADED_SYNTH_xi{xi}"
    hits = [p for p in d.glob("*_lines.csv") if p.name.endswith(f"_{treat}_lines.csv")]
    if len(hits) > 1:
        raise RuntimeError(
            f"{d.name} holds {len(hits)} candidate line files for {treat}. That is the "
            f"shared-worktree contamination this run was redone to remove — refusing to "
            f"guess which one is this pool's.")
    return pd.read_csv(hits[0]) if hits else None


def main() -> int:
    feed = json.loads(FEED.read_text())
    nuv = {(x["ion"], x["holding"], x["treatment"]): x
           for x in feed["products"] if x.get("band") == "near-UV"}
    rows = []
    for tag, (ion, holding) in POOLS.items():
        for treat in TREATMENTS:
            hi, lo = _lines(ion, tag, "1.10", treat), _lines(ion, tag, "0.90", treat)
            if hi is None or lo is None:
                print(f"  {tag} {treat}: no per-line pair — SKIP")
                continue
            r = paired_differential(hi, lo)
            p = nuv[(ion, holding, treat)]
            d_a_dxi = r.median / XI_SPAN_KMS
            sigma_xi = abs(d_a_dxi) * DELTA_XI_KMS
            publishable = r.n_paired >= MIN_PAIRED
            syst_new = math.sqrt(p["sigma_syst"] ** 2 + sigma_xi ** 2)
            rows.append({
                "pool": tag, "ion": ion, "holding": holding, "band": "near-UV",
                "tier": p.get("tier"), "treatment": treat,
                "A": p["A"], "n_lines": p["n_lines"], "n_paired": r.n_paired,
                "median_delta_dex": r.median,
                "dA_dxi": round(d_a_dxi, 4),
                "sigma_xi": round(sigma_xi, 4),
                "xi_state": "MEASURED" if publishable else "UNMEASURED",
                "xi_note": (f"|dA/dxi|={abs(d_a_dxi):.4f} x delta_xi={DELTA_XI_KMS} on "
                            f"{r.n_paired} paired lines"),
                "sigma_stat": p["sigma_stat"],
                "sigma_syst_published": p["sigma_syst"],
                "sigma_syst_with_xi": round(syst_new, 4),
                "sigma_reported_published": round(
                    math.hypot(p["sigma_stat"], p["sigma_syst"]), 4),
                "sigma_reported_with_xi": round(
                    math.hypot(p["sigma_stat"], syst_new), 4),
                "difference_of_aggregates": r.difference_of_aggregates,
                "collision": r.collision,
            })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "ticket": "RYA-1168", "read_only": True, "star": "solar",
        "band": "near-UV", "window_A": [3000, 3780],
        "xi_span_kms": XI_SPAN_KMS, "delta_xi_kms": DELTA_XI_KMS,
        "min_paired": MIN_PAIRED, "n_pools": len(rows),
        "n_collisions": sum(1 for r in rows if r["collision"]),
        "pools": rows,
    }, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(REPO)}  ({len(rows)} pools, "
          f"{sum(1 for r in rows if r['collision'])} collisions)")
    for r in rows:
        print(f"  Fe {r['ion']:3s} {r['pool']:9s} {r['treatment']:9s} "
              f"n_paired={r['n_paired']:3d} dA/dxi={r['dA_dxi']:+.4f} "
              f"sigma_xi={r['sigma_xi']:.4f} {r['xi_state']}"
              f"{'  COLLISION' if r['collision'] else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
