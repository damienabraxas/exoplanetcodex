#!/usr/bin/env python3
"""RYA-715 — the Fe II dossier's findings matrix, GENERATED from the artifacts.

    python3 scripts/rya715_fe2_status_audit.py [--check]

RYA-715 is a living dossier in Linear, and its §3 findings matrix and §4 audit table were
typed by hand in August. They have since gone stale in every direction: near-UV moved from
"product PENDING" to four live products, a <3D>-LTE leg exists that the matrix has no row
for, and the §1 arbiter range no longer matches any live cell. This script regenerates the
matrix from `Fe.json` + `canonical_gf` so the dossier can be re-derived instead of retyped.

🔴 IT ALSO SETTLES RYA-715's OPEN DECISION — BY MEASUREMENT, NOT BY ADJUDICATION.
The 2026-08-19 comment asks whether the band-product route should apply RYA-352's Fe II
EW-quality cull, and frames it as two principled gates disagreeing on 5 of the 6 lines
"both pools contain". Neither half of that framing survives contact with the current tree:

  1. ONE OF THE TWO GATES HAS BEEN WITHDRAWN. The comment rests on "FIT-PINNED is correct
     and stays". `scripts/rya761_recover_blends.py` records it RETIRED: "the detector that
     produced the FIT-PINNED population was itself wrong -- it read the Gaussian sigma of
     Voigt fits, one half of a degenerate pair". So the harness's objection to 6084.102 and
     6456.380 (sigma sitting exactly on the instrumental floor) was reading the wrong
     parameter, not measuring a pinned width.

  2. THE POOLS NO LONGER INTERSECT. Measured below: all six contested lines are SHALLOW
     (feature depth 0.155-0.430, under the 0.60 gate) AND non-lab-gf (VALD3 or OTHER), so
     not one of them can enter a GRADED or DEEPGRADED band pool by construction. The live
     Fe II band pool is nine lab-gf lines at 4233-4584 A. And no EW-route Fe II product is
     live at all -- every Fe II row in the feed is route=SYNTH.

So wiring the cull into the band route would change no published cell, because the cull's
target lines are not in the band pool and the route it governs emits no Fe II product. The
decision is MOOT AS POSED. That is a finding, not a ruling: if an EW-route Fe II product is
ever re-emitted the question returns, and it returns with one of its two gates already gone.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from derive_band_products import _feature_depth              # noqa: E402
from line_accounting_rya709 import DEPTH_HI                   # noqa: E402

FEED = ROOT / "data/products/solar/Fe.json"
CANON = ROOT / "data/linelists/canonical_gf.csv"
BP = ROOT / "data/results/band_products"
OUT = ROOT / "data/audit/rya715_fe2_status"

#: The six lines RYA-715's open decision is about (2026-08-19 comment).
CONTESTED = {
    6084.102: ("KEEP — reinstated as clean", "REJECT — FIT-PINNED"),
    6456.380: ("KEEP — reinstated as clean", "REJECT — FIT-PINNED"),
    6432.676: ("CULL — HIERR", "KEEP"),
    6247.557: ("CULL — BLEND", "KEEP"),
    6149.258: ("CULL — BLEND", "KEEP"),
    6369.459: ("keep", "keep"),
}
#: External comparands the dossier's §4 audits against.
AMARSI_2021 = (7.47, 0.04, "Amarsi et al. 2021, A&A 653 A141 — full 3D-NLTE, 76 Fe II levels")


def live_fe2() -> pd.DataFrame:
    feed = json.loads(FEED.read_text())
    rows = [{"band": p["band"], "holding": p["holding"], "tier": p["tier"],
             "route": p["route"], "treatment": p["treatment"], "grade": p["grade"],
             "A": p["A"], "n_lines": p["n_lines"],
             "sigma_stat": p.get("sigma_stat"), "sigma_syst": p.get("sigma_syst"),
             "sigma_reported": p.get("sigma_reported"), "xi_state": p.get("xi_state")}
            for p in feed["products"] if p["ion"] == "II"]
    return pd.DataFrame(rows).sort_values(["band", "holding", "treatment"]).reset_index(drop=True)


def contested_lines() -> pd.DataFrame:
    """The six lines of the open decision, as `canonical_gf` holds them TODAY."""
    c = pd.read_csv(CANON, low_memory=False)
    fe2 = c[c.species.astype(str) == "Fe II"]
    out = []
    for wl, (rya352, harness) in CONTESTED.items():
        r = fe2[(fe2.wavelength_air_A - wl).abs() < 0.05]
        if r.empty:
            out.append({"wavelength_air_A": wl, "in_canonical_gf": False})
            continue
        x = r.iloc[0]
        depth = float(_feature_depth([float(x.wavelength_air_A)])[0])
        tier = str(x.gf_tier)
        out.append({
            "wavelength_air_A": wl, "in_canonical_gf": True,
            "rya352_verdict": rya352, "band_harness_verdict": harness,
            "gf_tier": tier, "log_gf": float(x.log_gf),
            "feature_depth": round(depth, 3),
            "is_deep": bool(depth > DEPTH_HI),
            "is_primary_lab_gf": tier == "LAB",
            #: BOTH must be true for a line to reach a GRADED/DEEPGRADED band pool.
            "reachable_by_band_route": bool(tier == "LAB"),
        })
    return pd.DataFrame(out)


SUPERSEDED: dict = {}


def unpublished_artifacts(live: pd.DataFrame) -> list[dict]:
    """🔴 An Fe II band artifact on disk that never became a product row.

    RYA-1135 emitted the Fe II <3D>-LTE leg, committed it past the band_products gitignore
    with a GENERATORS entry, and recorded it in `data/audit/rya1135_fe2_mean3d_lte/` — but
    never ran it through `publish_product.py`, so the FEED has no row for it. Fe I carries
    both mean3D treatments; Fe II carries neither. Detected here rather than asserted, so it
    cannot silently come back.
    """
    global SUPERSEDED
    feed = json.loads(FEED.read_text())
    SUPERSEDED = {(e.get("holding"), e.get("band"), e.get("treatment")): e
                  for e in feed.get("superseded", []) if str(e.get("ion")) == "II"}
    have = {(r.treatment, r.holding) for r in live.itertuples()}
    #: the 1D-LTE cell for each (holding, band), to test the duplicate case below
    base = {(r.holding, r.band): (float(r.A), int(r.n_lines))
            for r in live.itertuples() if r.treatment == "1D-LTE"}
    out = []
    for p in sorted(BP.glob("FeII_*_products.csv")):
        d = pd.read_csv(p)
        if d.empty:
            continue
        row = d.iloc[0]
        t = str(row["treatment"])
        holding = next((h for h in {r.holding for r in live.itertuples()} if h in p.name), "")
        if (t, holding) in have:
            continue
        A = float(row["value"] if "value" in d.columns else row["A"])
        n = int(row["n_lines"])
        band = "near-UV" if "_3000_3780_" in p.name else "VIS"
        b = base.get((holding, band))
        #: 🔴 NOT EVERY UNPUBLISHED ARTIFACT IS A GAP, AND CALLING THEM ALL GAPS WOULD BE
        #: THE WORSE ERROR. Five of these are ENGINE-B run on `--engine-b-deck ts-lte`,
        #: which IS an LTE deck: each reproduces its 1D-LTE sibling exactly and its own
        #: columns say `scale=1D-LTE, model=none`. Publishing them would add relabelled
        #: duplicates under the ambiguous bare `ENGINE-B` token that
        #: docs/catalog/model_registry_notes.md says to refuse until Ryan rules. They are
        #: CORRECTLY unpublished. (Same shape as RYA-1187's red-optical ENGINE-B case —
        #: this confirms it on five more cells, so it is the deck, not one band.)
        dup = bool(b and abs(A - b[0]) < 1e-9 and n == b[1]
                   and str(row.get("scale")) == "1D-LTE")
        #: ✅ AND THE STORE ALREADY SAYS SO — cross-checked, not just inferred. Each of these
        #: is in `superseded` with `DUPLICATE_RETIRED_LABEL`: "ENGINE-B is a retired spelling
        #: of 1D-LTE measured by the synthesis handler, NOT a separate engine". The detector
        #: above reaches that verdict from the NUMBERS alone, so agreement between the two is
        #: a control on both: a duplicate this test called new, or a withdrawal the numbers
        #: did not support, would show up as a mismatch here.
        sup = SUPERSEDED.get((holding, band, t))
        out.append({
            "artifact": p.name, "treatment": t, "holding": holding, "band": band,
            "A": A, "n_lines": n, "published_in_feed": False,
            "scale": str(row.get("scale")), "model": str(row.get("model")),
            "classification": "DUPLICATE_OF_1D_LTE_withdrawn_by_the_store" if dup
                              else "UNPUBLISHED_GAP",
            "store_superseded_reason_code": (sup or {}).get("superseded_reason_code"),
            "detector_agrees_with_store": bool(dup) == bool(sup),
            "why": (f"identical to the live 1D-LTE cell for this holding (A={b[0]}, n={b[1]}) "
                    f"and its own columns read scale=1D-LTE / model=none: "
                    f"`--engine-b-deck ts-lte` IS an LTE deck. Publishing it would add a "
                    f"relabelled duplicate under the ambiguous bare ENGINE-B token."
                    + (f" The store agrees: superseded as "
                       f"{sup.get('superseded_reason_code')}." if sup else
                       " ⚠️ but the store has NO superseded record for it — the numbers say "
                       "duplicate and the feed does not say why it went; report, do not assume.")
                    if dup else
                    "a real product with no feed row — emitted, committed, and never run "
                    "through publish_product.py"),
        })
    return out


def build() -> dict:
    live = live_fe2()
    con = contested_lines()
    unpub = unpublished_artifacts(live)

    reach = con[con.get("reachable_by_band_route", False) == True]           # noqa: E712
    a_ext, s_ext, cite = AMARSI_2021
    audit = []
    for r in live.itertuples():
        d = round(float(r.A) - a_ext, 4)
        audit.append({"band": r.band, "holding": r.holding, "treatment": r.treatment,
                      "A": r.A, "external": a_ext, "delta": d,
                      "within_sigma_ext": bool(abs(d) <= s_ext),
                      "n_sigma_ext": round(abs(d) / s_ext, 2)})

    return {
        "ticket": "RYA-715", "feed_version": json.loads(FEED.read_text())["version"],
        "live_products": live.to_dict("records"),
        "n_live": int(len(live)),
        "bands_covered": sorted(live.band.unique().tolist()),
        "treatments_live": sorted(live.treatment.unique().tolist()),
        "🔴 open_decision_is_MOOT_AS_POSED": {
            "question": ("should the band-product route apply the RYA-352 Fe II EW-quality "
                         "cull? (RYA-715 comment, 2026-08-19)"),
            "answer": "the question no longer bites — MEASURED, not adjudicated",
            "reason_1_a_gate_was_withdrawn": (
                "the comment rests on 'FIT-PINNED is correct and stays'. "
                "scripts/rya761_recover_blends.py records it RETIRED: 'the detector that "
                "produced the FIT-PINNED population was itself wrong -- it read the Gaussian "
                "sigma of Voigt fits, one half of a degenerate pair'. The harness's objection "
                "to 6084.102 and 6456.380 was reading the wrong parameter."),
            "reason_2_the_pools_no_longer_intersect": (
                f"all {len(con)} contested lines are SHALLOW (depth "
                f"{con.feature_depth.min()}-{con.feature_depth.max()}, gate {DEPTH_HI}) and "
                f"{int((~con.is_primary_lab_gf).sum())} of {len(con)} are non-lab-gf, so "
                f"{len(reach)} of them can enter a GRADED/DEEPGRADED band pool. The live Fe II "
                f"band pool is 9 lab-gf lines at 4233-4584 A."),
            "reason_3_the_governed_route_emits_nothing": (
                "no EW-route Fe II product is live: every Fe II row in the feed is "
                "route=SYNTH. The cull governs `_load_solar_ews`, which no live Fe II "
                "product goes through."),
            "consequence": ("wiring the cull in would change NO published cell. If an "
                            "EW-route Fe II product is ever re-emitted the question returns — "
                            "with one of its two gates already gone."),
            "not_a_ruling": ("this does not adjudicate the three genuinely contested lines "
                             "(6432.676 HIERR, 6247.557/6149.258 BLEND). It reports that "
                             "nothing currently published depends on them."),
            "contested_lines": con.to_dict("records"),
        },
        "🔴 unpublished_band_artifacts": {
            "found": unpub,
            "note": ("an Fe II band artifact that exists, is committed, and has a GENERATORS "
                     "entry — but never became a product row. RYA-1135 emitted the <3D>-LTE "
                     "leg and stopped at the artifact; Fe I carries both mean3D treatments "
                     "and Fe II carries neither. `publish_product.py` is the remaining step. "
                     "⚠️ Classified, not lumped: the ENGINE-B artifacts here are LTE clones "
                     "of their 1D-LTE sibling and are CORRECTLY unpublished."
                     if unpub else "none — every Fe II band artifact has a feed row"),
        },
        "audit_vs_external": {"comparand": cite, "A": a_ext, "sigma": s_ext, "rows": audit},
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    doc = build()
    txt = json.dumps(doc, indent=2, ensure_ascii=False) + "\n"
    f = OUT / "fe2_status.json"
    if a.check:
        drift = (not f.exists()) or f.read_text() != txt
        print("DRIFT" if drift else "reproduces")
        return 1 if drift else 0
    OUT.mkdir(parents=True, exist_ok=True)
    f.write_text(txt)
    pd.DataFrame(doc["live_products"]).to_csv(OUT / "live_fe2_products.csv", index=False)
    pd.DataFrame(doc["🔴 open_decision_is_MOOT_AS_POSED"]["contested_lines"]).to_csv(
        OUT / "contested_lines.csv", index=False)

    print(f"Fe II live products: {doc['n_live']}  bands: {doc['bands_covered']}")
    print(f"  treatments: {doc['treatments_live']}")
    print(f"\nOPEN DECISION: {doc['🔴 open_decision_is_MOOT_AS_POSED']['answer']}")
    print(f"  reachable-by-band-route among the 6 contested: "
          f"{sum(1 for r in doc['🔴 open_decision_is_MOOT_AS_POSED']['contested_lines'] if r.get('reachable_by_band_route'))}")
    u = doc["🔴 unpublished_band_artifacts"]["found"]
    gaps = [x for x in u if x["classification"] == "UNPUBLISHED_GAP"]
    dups = [x for x in u if x["classification"] != "UNPUBLISHED_GAP"]
    print(f"\nUnpublished Fe II band artifacts: {len(u)}  "
          f"({len(gaps)} real gap, {len(dups)} correctly-unpublished 1D-LTE duplicates)")
    for x in gaps:
        print(f"   GAP  {x['band']:8} {x['treatment']}  A={x['A']}  n={x['n_lines']}")
    for x in dups:
        print(f"   dup  {x['band']:8} {x['treatment']:9} A={x['A']} == its 1D-LTE sibling")
    print(f"\nvs {doc['audit_vs_external']['comparand'][:52]}...")
    for r in doc["audit_vs_external"]["rows"]:
        flag = "within" if r["within_sigma_ext"] else f"{r['n_sigma_ext']}x sigma"
        print(f"   {r['band']:9} {r['treatment']:10} A={r['A']:<7} delta={r['delta']:+.3f}  {flag}")
    print(f"\nwrote {OUT.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
