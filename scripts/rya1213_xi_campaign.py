#!/usr/bin/env python3
"""RYA-1213 — dA/dxi for every Reference Grade pool. Ryan's 2026-09-12 requirement.

    python3 scripts/rya1213_xi_campaign.py --runs <dir> [--check]

WHY THIS EXISTS. Reference Grade is a NEW TIER, and the RYA-1120 campaign keys on
(ion x holding x tier x treatment x route). So every Reference product reached the feed
with `xi_state: NOT_IN_CAMPAIGN` and `sigma_xi: None` -- the xi term absent from its bar
rather than zero. On the flagship VIS Fe I pool that is a 40% understatement: 0.0524
published against its Deep sibling's 0.0862, with the gf half of the budget HONEST
(sigma_syst 0.0504 against Deep's 0.0402) and the whole difference the missing term.

🔴 RYAN'S RULING (2026-09-12): "Reference Grade gets the SAME uncertainty treatment as
Codex and Deep. The xi campaign must run for Reference Grade pools -- measured, not
borrowed from a sibling tier." Refusing to borrow was right; flagging it was not enough.

NOTHING HERE IS A NEW METHOD. Every choice below is RYA-1120's, reused rather than
re-decided:

  * CENTRAL DIFFERENCE at xi = 0.90 / 1.10, step 0.10 km/s. The production run at the
    pinned nominal 1.0 sits BETWEEN the legs, so it is not one of them and cannot be
    reused -- which is why the campaign costs two runs per unit and not one.
  * sigma_xi = |dA/dxi| * delta_xi, with delta_xi = 0.2912 (RYA-1089's SOURCED value).
    The slope times the STAR'S UNCERTAINTY, never times the perturbation step.
  * the derivative is a PER-LINE PAIRED DIFFERENTIAL through
    `pipeline.paired_differential`, not a difference of two aggregates. Line ACCEPTANCE
    depends on xi, so the two legs need not hold the same lines.
  * `min_paired = 3`. A pool that pairs fewer is reported UNMEASURED with its float
    carried but not published -- exactly as RYA-1163 handled the 2-line CRIRES+ ENGINE-A
    NIR pool. `xi_terms` reads the VERDICT, not the presence of a number.
  * the pairing is proven from each leg's own `xi_run.json` stamp via
    `pipeline.xi_pairing.assert_pair`. Worktree isolation is not a proof of pairing.

🔴 AND THE PAIRED METHOD EARNED ITSELF ON THE FIRST CELL MEASURED. NIR CRIRES+ 1D-LTE:
the paired median gives dA/dxi = -0.0200 and differencing the aggregates gives +0.0100 --
OPPOSITE SIGN. A difference-of-aggregates derivative would have charged that product a xi
term pointing the wrong way. See xi_control.json.

OUTPUT is the BAND-KEYED schema `rya1178_emit_fe_schema.xi_band_index` already consumes
(`pools[]` keyed on ion/holding/tier/treatment/band), so wiring it in is adding the path
to XI_BAND_RUNS and nothing else. A band-keyed run answers for the product's OWN band and
wins outright over the campaign, which is the right precedence here for the same reason
RYA-1114 F2 gave: the campaign key carries no band.
"""
from __future__ import annotations

import argparse
import glob
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from rya1120_xi_campaign import dA_dxi, sigma_xi          # noqa: E402  REUSED, not copied

OUT = ROOT / "data" / "results" / "rya1213" / "reference_xi_dadxi.json"

XI_NOMINAL, XI_STEP, DELTA_XI = 1.0, 0.10, 0.2912
#: RYA-1163's floor, adopted verbatim. Below it the derivative exists as a float and is
#: NOT a measurement; publishing it would be the RYA-1031 shape on a different quantity.
MIN_PAIRED = 3

#: job-directory stem -> the product axes it measured. Derived from the job name the
#: campaign runner used, so the mapping cannot drift from the runs.
BANDS = {"nuv": "near-UV", "vis": "VIS", "ropt": "red-optical", "nir": "NIR", "h": "H"}
HOLDINGS = {"kur": "solar_kpno_kurucz2005_corrected",
            "kpmf": "solar_kpno_molecfit_corrected",
            "harps": "solar_harps_molecfit_corrected",
            "iag": "solar_iag",
            "cry": "solar_crires_plus_y_wide_rya1054",
            "crh": "solar_crires_plus_h_rya1094"}
INSTRUMENTS = {"solar_kpno_kurucz2005_corrected": "kpno_solar_atlas",
               "solar_kpno_molecfit_corrected": "kpno_solar_atlas",
               "solar_harps_molecfit_corrected": "harps",
               "solar_iag": "iag_fts_solar_atlas",
               "solar_crires_plus_y_wide_rya1054": "crires_plus",
               "solar_crires_plus_h_rya1094": "crires_plus"}

JOB = re.compile(r"^(?P<band>nuv|vis|ropt|nir|h)(?P<ion>I{1,2})_(?P<holding>[a-z0-9]+)_"
                 r"(?P<deck>g1d|gnl|m3d|m3l)_xi(?P<xi>[01]\.\d+)$")


def units(runs: Path) -> dict:
    """(band, ion, holding, deck) -> {minus: dir, plus: dir}, from the run dirs alone."""
    found: dict = {}
    for d in sorted(p for p in runs.iterdir() if p.is_dir()):
        m = JOB.match(d.name)
        if not m or not (d / "DONE").exists():
            continue
        key = (BANDS[m.group("band")], m.group("ion"),
               HOLDINGS[m.group("holding")], m.group("deck"))
        leg = "minus" if float(m.group("xi")) < XI_NOMINAL else "plus"
        found.setdefault(key, {})[leg] = d
    return found


def _per_line(d: Path, treatment: str) -> pd.DataFrame | None:
    """The per-line CSV for ONE treatment in ONE leg's own directory.

    ⚠️ MATCHED ON THE FULL SUFFIX. `_1D-LTE_lines.csv` is a substring of nothing else
    here, but `ENGINE-A` is a substring of `ENGINE-A-3DNLTE`, so a loose `in` test would
    silently pair a product with a different engine's per-line file.
    """
    hits = [f for f in glob.glob(str(d / "*_lines.csv"))
            if f.endswith(f"_{treatment}_lines.csv")]
    if len(hits) != 1:
        return None
    return pd.read_csv(hits[0])


def treatments(d: Path) -> list[str]:
    out = set()
    for f in glob.glob(str(d / "*_lines.csv")):
        m = re.search(r"_SYNTH_REFERENCE_(.+)_lines\.csv$", Path(f).name)
        if m:
            out.add(m.group(1))
        elif Path(f).name.endswith("_SYNTH_REFERENCE_1D-LTE_lines.csv"):
            out.add("1D-LTE")
    return sorted(out)


def build(runs: Path) -> dict:
    u = units(runs)
    pools, incomplete = [], []
    for (band, ion, holding, deck), legs in sorted(u.items()):
        if set(legs) != {"minus", "plus"}:
            # 🔴 NEVER HALF-DIFFERENCED. RYA-1120 §4: if one side of a pair is missing,
            # differencing what is left is not a smaller measurement, it is a different
            # quantity.
            incomplete.append({"band": band, "ion": f"Fe {ion}", "holding": holding,
                               "deck": deck, "legs_present": sorted(legs),
                               "verdict": "UNMEASURED — a pair is missing a leg"})
            continue
        lo, hi = legs["minus"], legs["plus"]
        for t in treatments(hi):
            L, H = _per_line(lo, t), _per_line(hi, t)
            if L is None or H is None:
                incomplete.append({"band": band, "ion": f"Fe {ion}", "holding": holding,
                                   "deck": deck, "treatment": t,
                                   "verdict": "UNMEASURED — a leg has no per-line file "
                                              "for this treatment"})
                continue
            d = dA_dxi(H, L, step_kms=XI_STEP, minus_dir=lo, plus_dir=hi,
                       xi_nominal=XI_NOMINAL)
            n_paired = int(d.get("n_paired") or 0)
            slope = round(float(d["dA_dxi_paired"]), 4)
            ok = n_paired >= MIN_PAIRED
            pools.append({
                "pool": f"{band}/{holding}/{deck}", "instrument": INSTRUMENTS[holding],
                "holding": holding, "band": band, "ion": ion,
                "tier": "REFERENCE", "treatment": t,
                "n_paired": n_paired,
                "median_delta_dex": round(float(d["median"]), 4) if d.get("median") is not None else None,
                "dA_dxi": slope,
                "sigma_xi": round(sigma_xi(slope, DELTA_XI), 6),
                "xi_state": "MEASURED" if ok else "UNMEASURED",
                "xi_note": (f"|dA/dxi|={abs(slope):.4f} x delta_xi={DELTA_XI} on "
                            f"{n_paired} paired lines"
                            if ok else
                            f"only {n_paired} paired line(s), below min_paired="
                            f"{MIN_PAIRED} — the derivative exists as a float and is not "
                            f"a measurement (RYA-1163's floor, RYA-1031's reason)"),
                "difference_of_aggregates": round(float(d["dA_dxi_from_aggregates"]), 4),
                "pool_moved": bool(d["pool_moved"]),
                "n_lo": d["n_lo"], "n_hi": d["n_hi"],
                "span_kms": d["span_kms"],
                "xi_pairing_verified": bool(d["xi_pairing_verified"]),
                # 🔴 CARRIED PER POOL, BECAUSE IT IS THE REASON THE PAIRED ROUTE IS USED.
                # Where these two disagree in SIGN, a difference-of-aggregates derivative
                # would charge the product a xi term pointing the wrong way.
                "sign_disagreement_with_aggregates": bool(
                    slope * float(d["dA_dxi_from_aggregates"]) < 0),
            })
    return {
        "ticket": "RYA-1213",
        "read_only": True,
        "requirement": ("Ryan, 2026-09-12 — Reference Grade gets the SAME uncertainty "
                        "treatment as Codex and Deep; measured, not borrowed."),
        "star": "solar",
        "species": "Fe",
        "method": ("RYA-1120's, reused: central difference at xi = "
                   f"{XI_NOMINAL - XI_STEP}/{XI_NOMINAL + XI_STEP}, per-line paired "
                   "differential via pipeline.paired_differential, pairing proven from "
                   "each leg's xi_run.json stamp."),
        "xi_span_kms": 2 * XI_STEP,
        "delta_xi_kms": DELTA_XI,
        "min_paired": MIN_PAIRED,
        "n_pools": len(pools),
        "n_measured": sum(1 for p in pools if p["xi_state"] == "MEASURED"),
        "n_sign_disagreements": sum(1 for p in pools
                                    if p["sign_disagreement_with_aggregates"]),
        "incomplete": incomplete,
        "pools": pools,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--runs", required=True, help="directory of per-leg run directories")
    ap.add_argument("--check", action="store_true", help="print only; write nothing")
    a = ap.parse_args()
    doc = build(Path(a.runs).expanduser())
    print(f"\n{doc['n_pools']} pool(s), {doc['n_measured']} MEASURED, "
          f"{len(doc['incomplete'])} incomplete, "
          f"{doc['n_sign_disagreements']} sign disagreement(s) with the aggregate route\n")
    print(f"{'band':12s} {'ion':4s} {'holding':32s} {'treatment':22s} "
          f"{'dA/dxi':>8s} {'sigma_xi':>9s} {'n':>4s}  state")
    for p in doc["pools"]:
        flag = "  <- SIGN" if p["sign_disagreement_with_aggregates"] else ""
        print(f"{p['band']:12s} {p['ion']:4s} {p['holding']:32s} {p['treatment']:22s} "
              f"{p['dA_dxi']:+8.4f} {p['sigma_xi']:9.5f} {p['n_paired']:4d}  "
              f"{p['xi_state']}{flag}")
    for i in doc["incomplete"]:
        print(f"  INCOMPLETE {i}")
    if not a.check:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(doc, indent=1) + "\n")
        print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
