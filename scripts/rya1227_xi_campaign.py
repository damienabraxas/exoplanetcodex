#!/usr/bin/env python3
"""RYA-1227 — dA/dxi for the Codex/Deep 1D-engine pools the band campaigns never reached.

    python3 scripts/rya1227_xi_campaign.py --runs <dir> [--check]

WHY THIS EXISTS. A feed audit found products with `xi_state: NOT_IN_CAMPAIGN` -- the xi
perturbation was never run on them, so their `sigma_reported` carries no xi term at all.
They are not the mean-3D or full-3D routes that legitimately sit out xi; they are 1D
engines where microturbulence applies, and the SAME engines are MEASURED elsewhere.

🔴 ROOT CAUSE, MEASURED FROM THE ARTIFACTS RATHER THAN ASSUMED. Every one of these product
artifacts is stamped 2026-09-08 or 2026-09-10, and every band-keyed xi run predates it:
NIR and H on 2026-08-30 (RYA-1163), near-UV and red-optical on 2026-09-03, VIS older still
(the RYA-1120 campaign). The Gerber fan-out (RYA-1208) landed AFTER the campaigns, so the
cells it created never had a xi leg to be part of. Confirmed for all of them, not inferred
from one.

NOTHING HERE IS A NEW METHOD. Every choice is RYA-1120's, reached through RYA-1213:

  * CENTRAL DIFFERENCE at xi = 0.90 / 1.10, step 0.10 km/s.
  * sigma_xi = |dA/dxi| * delta_xi, delta_xi = 0.2912 (RYA-1089's SOURCED value).
  * a PER-LINE PAIRED DIFFERENTIAL through `pipeline.paired_differential`, never a
    difference of two aggregates -- which on the first cell RYA-1213 measured gave the
    OPPOSITE SIGN.
  * `min_paired` from the runs, and the pool that pairs fewer is UNMEASURED with its float
    carried but not published.
  * pairing proven from each leg's own `xi_run.json` stamp.

`dA_dxi`, `sigma_xi`, `_per_line` and the constants are IMPORTED from the RYA-1213
campaign, which imports them from RYA-1120's. Three copies of a derivative formula is how
two of them drift.

🔴 THIS RUN EMITS ONLY EACH DECK'S OWN TREATMENT, AND THAT IS A DELIBERATE NARROWING.
A `g1d` leg also produces 1D-LTE and ENGINE-A products, and they are measured on the same
pool -- but those cells are ALREADY MEASURED in the feed from the band campaigns, and
publishing a second derivative for them would MOVE a published sigma on products this
ticket was not asked to touch. The better-keyed base-treatment derivatives exist in the run
directories and are reported, not emitted. `xi_band_index` would refuse a duplicate key
loudly in any case; this makes the refusal unnecessary rather than relying on it.

OUTPUT is the band-keyed schema `rya1178_emit_fe_schema.xi_band_index` already consumes,
so wiring it in is adding the path to XI_BAND_RUNS and nothing else.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

# REUSED, not copied -- the derivative, the sigma, the per-line reader and every constant
# come from the campaign that established them.
from rya1213_xi_campaign import (  # noqa: E402
    BANDS, DELTA_XI, HOLDINGS, INSTRUMENTS, MIN_PAIRED, XI_NOMINAL, XI_STEP,
    _per_line, dA_dxi, sigma_xi,
)

OUT = ROOT / "data" / "results" / "rya1227" / "codex_deep_xi_dadxi.json"

#: The tier segment RYA-1213's regex has no room for: its campaign was Reference-only and
#: hardcodes the tier, which is exactly why it cannot answer for these pools.
TIERS = {"g": "GRADED", "dg": "DEEPGRADED"}

#: Each deck's OWN treatment -- the one this campaign is here to measure.
OWN = {"g1d": "synth-1D-LTE-gerber", "gnl": "ENGINE-B-NLTE"}

JOB = re.compile(r"^(?P<band>nuv|vis|ropt|nir|h)(?P<ion>I{1,2})_(?P<holding>[a-z0-9]+)_"
                 r"(?P<tier>g|dg)_(?P<deck>g1d|gnl)_xi(?P<xi>[01]\.\d+)$")


def units(runs: Path) -> dict:
    """(band, ion, holding, tier, deck) -> {minus: dir, plus: dir}, from the dirs alone."""
    found: dict = {}
    for d in sorted(p for p in runs.iterdir() if p.is_dir()):
        m = JOB.match(d.name)
        if not m or not (d / "DONE").exists():
            continue
        key = (BANDS[m.group("band")], m.group("ion"), HOLDINGS[m.group("holding")],
               TIERS[m.group("tier")], m.group("deck"))
        leg = "minus" if float(m.group("xi")) < XI_NOMINAL else "plus"
        found.setdefault(key, {})[leg] = d
    return found


def build(runs: Path) -> dict:
    u = units(runs)
    pools, incomplete, not_emitted = [], [], []
    for (band, ion, holding, tier, deck), legs in sorted(u.items()):
        if set(legs) != {"minus", "plus"}:
            #: 🔴 NEVER HALF-DIFFERENCED (RYA-1120 §4). If one side of a pair is missing,
            #: differencing what is left is not a smaller measurement, it is a different
            #: quantity.
            incomplete.append({"band": band, "ion": f"Fe {ion}", "holding": holding,
                               "tier": tier, "deck": deck,
                               "legs_present": sorted(legs),
                               "verdict": "UNMEASURED — a pair is missing a leg"})
            continue
        lo, hi = legs["minus"], legs["plus"]
        t = OWN[deck]
        L, H = _per_line(lo, t), _per_line(hi, t)
        if L is None or H is None:
            incomplete.append({"band": band, "ion": f"Fe {ion}", "holding": holding,
                               "tier": tier, "deck": deck, "treatment": t,
                               "verdict": "UNMEASURED — a leg has no per-line file for "
                                          "this treatment"})
            continue
        d = dA_dxi(H, L, step_kms=XI_STEP, minus_dir=lo, plus_dir=hi,
                   xi_nominal=XI_NOMINAL)
        n_paired = int(d.get("n_paired") or 0)
        slope = round(float(d["dA_dxi_paired"]), 4)
        ok = n_paired >= MIN_PAIRED
        pools.append({
            "pool": f"{band}/{holding}/{tier}/{deck}",
            "instrument": INSTRUMENTS[holding], "holding": holding, "band": band,
            "ion": ion, "tier": tier, "treatment": t,
            "n_paired": n_paired,
            "median_delta_dex": (round(float(d["median"]), 4)
                                 if d.get("median") is not None else None),
            "dA_dxi": slope,
            "sigma_xi": round(sigma_xi(slope, DELTA_XI), 6),
            "xi_state": "MEASURED" if ok else "UNMEASURED",
            "xi_note": (f"|dA/dxi|={abs(slope):.4f} x delta_xi={DELTA_XI} on "
                        f"{n_paired} paired lines"
                        if ok else
                        f"only {n_paired} paired line(s), below min_paired={MIN_PAIRED} "
                        f"— the derivative exists as a float and is not a measurement "
                        f"(RYA-1163's floor, RYA-1031's reason)"),
            "difference_of_aggregates": round(float(d["dA_dxi_from_aggregates"]), 4),
            "pool_moved": bool(d["pool_moved"]),
            "n_lo": d["n_lo"], "n_hi": d["n_hi"], "span_kms": d["span_kms"],
            "xi_pairing_verified": bool(d["xi_pairing_verified"]),
            #: 🔴 CARRIED PER POOL, because it is the reason the paired route is used.
            #: Where these disagree in SIGN, a difference-of-aggregates derivative would
            #: charge the product a xi term pointing the wrong way.
            "sign_disagreement_with_aggregates": bool(
                slope * float(d["dA_dxi_from_aggregates"]) < 0),
        })
        #: What the same legs measured for the base treatments, recorded so the narrowing
        #: above is auditable rather than merely asserted.
        for base in ("1D-LTE", "ENGINE-A"):
            bl, bh = _per_line(lo, base), _per_line(hi, base)
            if bl is None or bh is None:
                continue
            bd = dA_dxi(bh, bl, step_kms=XI_STEP, minus_dir=lo, plus_dir=hi,
                        xi_nominal=XI_NOMINAL)
            not_emitted.append({
                "band": band, "ion": ion, "holding": holding, "tier": tier,
                "treatment": base, "n_paired": int(bd.get("n_paired") or 0),
                "dA_dxi": round(float(bd["dA_dxi_paired"]), 4),
                "why": "already MEASURED in the feed from a band campaign; emitting a "
                       "second derivative would move a published sigma out of scope",
            })

    dupes = [k for k in
             {(p["ion"], p["holding"], p["tier"], p["treatment"], p["band"]) for p in pools}
             if sum(1 for p in pools
                    if (p["ion"], p["holding"], p["tier"], p["treatment"],
                        p["band"]) == k) > 1]
    if dupes:
        raise SystemExit(f"two derivatives for one band key: {dupes} — the feed would "
                         f"refuse this, and the run set is wrong (RYA-1227)")

    return {
        "ticket": "RYA-1227",
        "read_only": True,
        "requirement": ("Every 1D-engine product carries a MEASURED xi term before "
                        "publish. NOT_IN_CAMPAIGN is not a disposition, it is an "
                        "unfinished campaign."),
        "star": "solar", "species": "Fe",
        "method": ("RYA-1120's, reused through RYA-1213: central difference at "
                   "xi = 0.9/1.1, per-line paired differential via "
                   "pipeline.paired_differential, pairing proven from each leg's "
                   "xi_run.json stamp."),
        "root_cause": ("RYA-1208's Gerber fan-out landed AFTER the band-keyed xi "
                       "campaigns: every product artifact here is stamped 2026-09-08 or "
                       "2026-09-10, every band xi run 2026-08-30 (RYA-1163) or "
                       "2026-09-03 (RYA-1168/red-optical) or earlier. The cells did not "
                       "exist when the campaigns fanned out."),
        "xi_span_kms": 0.2, "delta_xi_kms": DELTA_XI, "min_paired": MIN_PAIRED,
        "n_pools": len(pools),
        "n_measured": sum(1 for p in pools if p["xi_state"] == "MEASURED"),
        "incomplete": incomplete,
        "emits_own_treatment_only": OWN,
        "base_treatments_measured_but_not_emitted": not_emitted,
        "pools": pools,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", required=True, help="directory of per-leg run directories")
    ap.add_argument("--check", action="store_true", help="print only; write nothing")
    a = ap.parse_args()

    doc = build(Path(a.runs).expanduser())
    for p in doc["pools"]:
        flag = "  <-- SIGN DISAGREES WITH AGGREGATES" if p["sign_disagreement_with_aggregates"] else ""
        print(f"  {p['ion']:3} {p['band']:12} {p['tier']:11} {p['treatment']:20} "
              f"n_paired={p['n_paired']:<4} dA/dxi={p['dA_dxi']:+.4f} "
              f"sigma_xi={p['sigma_xi']:.6f}  {p['xi_state']}{flag}")
    for i in doc["incomplete"]:
        print(f"  INCOMPLETE {i}")
    print(f"\npools: {doc['n_pools']}  measured: {doc['n_measured']}  "
          f"incomplete: {len(doc['incomplete'])}")
    print(f"base-treatment derivatives measured but NOT emitted: "
          f"{len(doc['base_treatments_measured_but_not_emitted'])}")
    if a.check:
        print("\n--check: nothing written")
        return 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
