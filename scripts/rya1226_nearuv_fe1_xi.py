#!/usr/bin/env python3
"""
RYA-1226 part C -- near-UV Fe I dA/dxi, measured on the CURRENT molecular synthesis.

WHY A NEW RUN AND NOT AN ADOPTION. RYA-1225 could adjudicate near-UV Fe II by adopting
RYA-1213's post-opacity slope because the two tiers were provably ONE pool. Fe I is not
that case -- RYA-1209 dropped a line, so DEEPGRADED is 55/54 against REFERENCE's 56/55 --
and the ticket forbids same-pool adoption here for exactly that reason. So the slope was
measured on the products' OWN DEEPGRADED pool.

PROVENANCE. Sirius, ~/scratch/rya1213 @ 97cde328 (post-RYA-1207: `use_molecules: true` and
`write_molecular_bsyn` both present, RYA-1207 in history), venv_rya1103 / Python 3.12.13,
ISPEC_DIR=/mnt/codex-data/engines/ispec_src. The checkout carried 6 uncommitted files at
run time; all six are NIR IAG band_products artifacts (9199-11081 A) plus one untracked
audit directory, and NONE touches pipeline/, scripts/, config/ or the linelists -- so the
synthesis path is the committed one. 4 legs, xi 0.90/1.10, ~80 and ~57 minutes per pair.

🔴 THE DERIVATIVE USES pipeline.paired_differential, NOT A HAND MERGE. That function
filters both legs to `in_aggregate`; a raw merge does not, and the excluded rows still
carry abundances. My first pass merged everything and read n_paired = 57 against a 55-line
product -- which looked like a moved pool and was my own arithmetic. Same function both
prior runs used, same 0.20 span, same 0.2912 delta_xi.

🔴 AND ONE POOL REALLY DID MOVE. kurucz2005 pairs 55 of 55 -- the product's own pool
exactly. molecfit pairs 53 of 54: Fe I 3427.119 A is in the product's aggregate and in the
xi=0.90 leg, but the xi=1.10 leg drops it ("NON-MINIMUM: chi2 at a bracket end is not above
the reported minimum"). `paired_response`'s own rule is that a moved pool is HOLD at the
caller, not an automatic restriction to whichever lines survived -- so molecfit is reported
UNMEASURED with its float carried, and kurucz2005 is the only cell this run can publish.
"""

from __future__ import annotations

import json
import pathlib
import sys

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.paired_differential import paired_differential  # noqa: E402
from pipeline.xi_pairing import assert_pair                    # noqa: E402

LEGS = ROOT / "data/results/rya1226/legs"
FEED = ROOT / "data/products/solar/Fe.json"
OUT = ROOT / "data/results/rya1226/nearuv_fe1_xi_dadxi.json"

XI_SPAN_KMS = 0.20
DELTA_XI_KMS = 0.2912
POOLS = {"kur": "solar_kpno_kurucz2005_corrected",
         "kpmf": "solar_kpno_molecfit_corrected"}
TREATMENTS = ("1D-LTE", "ENGINE-A")


def main() -> int:
    feed = json.loads(FEED.read_text())
    live = {(p["ion"], p["holding"], p["tier"], p["treatment"], p["band"]): p
            for p in feed["products"]}
    pools = []
    for tag, holding in POOLS.items():
        lo_d, hi_d = LEGS / f"nuvI_{tag}_dg_xi0.9000", LEGS / f"nuvI_{tag}_dg_xi1.1000"
        try:
            assert_pair(lo_d, hi_d)
            pairing = True
        except Exception as exc:                                  # noqa: BLE001
            pairing = False
            print(f"  {tag}: PAIRING REFUSED -- {exc}")
        for treat in TREATMENTS:
            stem = (f"FeI_3000_3780_kpno_solar_atlas_{holding}_SYNTH_DEEPGRADED_"
                    f"{treat}_lines.csv")
            r = paired_differential(pd.read_csv(hi_d / stem), pd.read_csv(lo_d / stem))
            product = live[("I", holding, "DEEPGRADED", treat, "near-UV")]
            slope = r.median / XI_SPAN_KMS
            moved = r.n_paired != product["n_lines"]
            pools.append({
                "pool": f"near-UV/{holding}/deep-graded", "band": "near-UV", "ion": "I",
                "instrument": "kpno_solar_atlas", "holding": holding,
                "tier": "DEEPGRADED", "treatment": treat,
                "n_paired": r.n_paired, "product_n_lines": product["n_lines"],
                "pool_moved": moved,
                "median_delta_dex": r.median,
                "difference_of_aggregates": r.difference_of_aggregates,
                "n_nonzero": r.n_nonzero, "sd": r.sd,
                "dA_dxi": None if moved else round(slope, 4),
                "dA_dxi_unpublishable_float": round(slope, 4) if moved else None,
                "sigma_xi": None if moved else round(abs(slope) * DELTA_XI_KMS, 6),
                "xi_state": "UNMEASURED" if moved else "MEASURED",
                "xi_pairing_verified": pairing,
                "superseded_dA_dxi": product.get("dA_dxi_dex_per_kms"),
                "superseded_sigma_xi": product.get("sigma_xi"),
                "A": product["A"],
                "xi_note": (
                    f"|dA/dxi|={abs(slope):.4f} x delta_xi={DELTA_XI_KMS} on {r.n_paired} "
                    f"paired lines, measured on THIS product's own DEEPGRADED pool on the "
                    f"current molecular synthesis (RYA-1226 part C)"
                    if not moved else
                    f"UNMEASURED: the pool moved under perturbation -- {r.n_paired} of "
                    f"{product['n_lines']} lines pair. Fe I 3427.119 A is in the product's "
                    f"aggregate and in the xi=0.90 leg, but the xi=1.10 leg drops it "
                    f"(NON-MINIMUM: chi2 at a bracket end). A moved pool is HOLD, not a "
                    f"restriction to the survivors (pipeline.uncertainty_contract."
                    f"paired_response)."),
            })

    doc = {
        "ticket": "RYA-1226", "part": "C", "star": "solar", "band": "near-UV", "ion": "I",
        "read_only": True,
        "xi_span_kms": XI_SPAN_KMS, "delta_xi_kms": DELTA_XI_KMS, "min_paired": 3,
        "method": ("RYA-1120's, reused: central difference at xi 0.90/1.10, per-line paired "
                   "differential via pipeline.paired_differential (which filters both legs "
                   "to in_aggregate), pairing proven from each leg's xi_run.json stamp via "
                   "pipeline.xi_pairing.assert_pair."),
        "why_not_adopted": ("RYA-1209 dropped a line, so the DEEPGRADED and REFERENCE pools "
                            "are 55/54 against 56/55. No same-pool adoption is permitted "
                            "here -- this is a measurement on the products' own pool."),
        "provenance": {
            "host": "Sirius", "checkout": "~/scratch/rya1213", "commit": "97cde328",
            "post_molecular": ("use_molecules: true and write_molecular_bsyn both present; "
                               "RYA-1207 in history"),
            "python": "venv_rya1103 / Python 3.12.13",
            "ispec_dir": "/mnt/codex-data/engines/ispec_src",
            "uncommitted_at_run_time": ("6 files, all NIR IAG band_products artifacts "
                                        "(9199-11081 A) plus one untracked audit directory; "
                                        "none touches pipeline/, scripts/, config/ or the "
                                        "linelists, so the synthesis path is the committed one"),
            "legs": ["nuvI_kur_dg_xi0.9000", "nuvI_kur_dg_xi1.1000",
                     "nuvI_kpmf_dg_xi0.9000", "nuvI_kpmf_dg_xi1.1000"],
            "wall_clock_utc": "13:19:17 -> 15:35:27 on 2026-09-18, P=2",
        },
        "n_pools": len(pools),
        "n_measured": sum(1 for p in pools if p["xi_state"] == "MEASURED"),
        "pools": pools,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")
    for p in pools:
        print(f"  {p['holding'][:31]:31} {p['treatment']:9} n_paired={p['n_paired']:>3}/"
              f"{p['product_n_lines']:<3} dA/dxi={str(p['dA_dxi']):>8} "
              f"(was {p['superseded_dA_dxi']})  {p['xi_state']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
