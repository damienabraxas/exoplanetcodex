"""
scripts/redoptical_xi.py — red-optical dA/dxi
=======================================
dA/dxi for the three red-optical (6910-9199 A) GRADED Fe I pools, as a PER-LINE PAIRED
DIFFERENTIAL. Third sibling of `scripts/rya1163_ir_xi.py` and `rya1168_nearuv_xi.py`,
and the LAST band: with this, all 66 live Fe products carry a measured xi term.

WHY THIS EXISTS. The census behind RYA-1114 F2 found dA/dxi measured on 46 of 66 live
products -- ALL of them VIS. red-optical was 0 of 6 and, unlike IR and near-UV, had no
audit ticket naming it at all; it surfaced only because the census asked every band
rather than the one it was sent to look at.

🔴 RED-OPTICAL IS THE MOST xi-SENSITIVE BAND MEASURED. |dA/dxi| 0.135-0.170, against
near-UV Fe I 0.105-0.113 and IR 0.020-0.053. Its bars move +0.008 to +0.016 dex -- an
order of magnitude more than IR's +0.0002..+0.0007. The band nobody had a ticket for is
the one where the missing term mattered most.

🔴 THE ARM WAS NEVER SKIPPED; THE BAND WAS. `~/xi_campaign/run_campaign.py:132` hardcodes
`--lo 4200 --hi 6910` on EVERY unit. All three of these holdings WERE in that campaign --
they measured VIS. `--hi 6910` is exactly this band's LOWER edge, so red-optical begins
one Angstrom past where every campaign unit stopped.

⚠️ THE WINDOW IS 6910-9199, NOT 6910-10000. `band_policy.POLICIES` and
`config/synth_bands.yaml` disagree about this band's red edge (RYA-1114 F9). The
PRODUCTS are cut at 9199, so that is the window a perturbation of them must use.

🔴 AND HERE THE PAIRED DIFFERENTIAL CATCHES SOMETHING A COUNT CHECK CANNOT. The IAG
1D-LTE pool holds 69 in-aggregate lines at BOTH xi -- and they are not the same 69.
9119.97 A leaves and 8432.17 A enters, so only 68 pair. Both are marginal NON-MINIMUM
fits (frac_rise -6.5e-05 and -3.6e-05, i.e. zero to within noise) whose acceptance flips
on numerical noise. A guard comparing n_lines would read 69 == 69 and conclude the pool
was stable. RYA-1083 says line acceptance moves with xi; this is that, with the COUNT
held fixed, which is the version no aggregate check can see.

⚠️ RUN UNDER PER-WORKER WORKTREE ISOLATION FROM THE START (`run_redopt.py`, workers
rw0..rw2). The near-UV run had to be thrown away once because four concurrent jobs shared
one worktree and each job's move-out swept up a sibling's output. Nothing in a run's own
artifacts records which xi produced it, so ISOLATION IS THE ONLY THING THAT MAKES THE
PAIRING TRUE -- there is no downstream check that would have caught it.
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

CAMPAIGN = pathlib.Path.home() / "xi_campaign" / "out_redopt"
FEED = REPO / "data" / "products" / "solar" / "Fe.json"
OUT = REPO / "data" / "results" / "redopt_xi" / "redopt_xi_dadxi.json"

XI_SPAN_KMS = 0.20
DELTA_XI_KMS = 0.2912
#: Below this many paired lines a slope is not publishable (RYA-907: unmeasured is not
#: a value). Every near-UV pool clears it; the guard is kept so it stays true by CHECK.
MIN_PAIRED = 3

POOLS = {
    "iag": ("I", "solar_iag"),
    "kur": ("I", "solar_kpno_kurucz2005_corrected"),
    "mol": ("I", "solar_kpno_molecfit_corrected"),
}
TREATMENTS = ("1D-LTE", "ENGINE-A")


def _lines(ion: str, tag: str, xi: str, treat: str) -> pd.DataFrame | None:
    d = CAMPAIGN / f"Fe{ion}_{tag}_redopt_GRADED_SYNTH_xi{xi}"
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
           for x in feed["products"] if x.get("band") == "red-optical"}
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
                "pool": tag, "ion": ion, "holding": holding, "band": "red-optical",
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
        "band": "red-optical", "window_A": [3000, 3780],
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
