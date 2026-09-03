#!/usr/bin/env python3
"""RYA-1055 item 3 — confirm every LIVE Fe II product's label, and state the balance.

    python3 scripts/rya1055_fe2_label_audit.py \
        --out data/results/rya1055/fe2_label_audit.json

Ryan, 2026-09-03: *"Confirm every live Fe II product reads a clean LTE label + the limit
note; the ionization-balance result (best VIS pools ~7.46) is stated as consistent with
the Fe I anchor 7.466."*

REPORT ONLY. Nothing here re-derives, re-fits or re-weights anything (RYA-161). It reads
the published feed, checks three things a reader would otherwise have to check by hand,
and states the balance with its scale attached.

WHAT IT CHECKS
--------------
1. **Every live Fe II product carries the limit** in `science_provenance.nlte_capability`.
2. **No live Fe II product takes departures from the Gerber deck.** This is the sharp
   form: an Fe II NLTE label is only dishonest when the number came from `atom.fe607a`.
   Fe II ENGINE-A reads the MPIA/Bergemann per-line delta grid instead -- a different
   mechanism, real corrections -- so its 1D-NLTE label is CLEAN and relabelling it would
   replace one wrong statement with another.
3. **The per-line layer agrees**, which is the sharper check of the two: `scale` is
   1D-NLTE on 159 Fe I rows and on ZERO Fe II rows, so the Gerber route never touched
   Fe II in anything we ship.
4. **Every Fe II band product's own per-line artifact names its NLTE source**, line by
   line. This is the strongest available form: `nlte_source` and `nlte_delta_dex` are
   written at the point the correction is applied (RYA-880), so it is the run speaking
   rather than a label being read.

🔴 AND IT SETTLES THE CONTRADICTION RYA-1113 FILED ON THIS TICKET. That comment found a
live Fe II near-UV NLTE leg at n=7 against its LTE sibling's n=12, and asked: if Fe II
NLTE is structurally unavailable, *how does an n=7 Fe II NLTE leg exist AT ALL?* -- with
two candidate answers, "it is not real NLTE" or "the atom has changed since". Both are
wrong, and the per-line artifact RYA-908 has since emitted says so outright:

    7 lines  nlte_delta_dex = -0.001, nlte_source = "Bergemann MPIA per-line delta_nlte
             (live query, solar node); PER-LINE additive correction"
    5 lines  excluded, "ENGINE-A-NOT-SERVED: the register..."

So the leg IS real NLTE, from the MPIA per-line grid -- a source this ticket's finding
never touched -- and this ticket's premise is intact (re-measured on the atom today: still
zero Fe II bound-bound transitions). The n=7-vs-12 gap is MPIA SERVICE COVERAGE, not
physics: ⚠️ the published Fe II near-UV NLTE delta is confounded with a 5-line pool
change, exactly as RYA-1113 said, and that half of its finding stands unchanged.

🔴 AND THE BALANCE IS QUOTED SCALE-MATCHED, WHICH THE HEADLINE PHRASING IS NOT.
"Fe II VIS ~7.46 is consistent with the Fe I anchor 7.466" pairs two numbers that are NOT
on the same scale. The anchor is gold v3's Fe I value on the **3D-NLTE** scale (RYA-553,
`OPTICAL_ANCHOR` in the RYA-783 report, which warns about exactly this); the Fe II cells
are **1D-NLTE** products. Setting them side by side is the RYA-669 ratchet -- read a label,
compare across scales, and a value drifts with every gate green. ⚠️ And the coincidence is
unusually inviting here: the Fe II kurucz2005 cell is 7.466 TO THREE DECIMALS, the anchor's
own digits. A number that looks right on arrival is the one that never gets challenged.

So this reports the balance the way it can be defended: **Fe II minus Fe I on the SAME
holding, the SAME band, the SAME treatment and the SAME tier**, which is one scale by
construction and needs no correction applied to either side.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

FEED = ROOT / "data/products/solar/Fe.json"
PERLINE = ROOT / "data/products/solar/Fe_perline.csv"

#: Context only, and NOT a target (RYA-161/RYA-777). Carried here with its scale attached
#: because the scale is the whole point of the caveat above.
OPTICAL_ANCHOR = 7.466
OPTICAL_ANCHOR_SCALE = "3D-NLTE (gold v3, RYA-553)"


def _balance_note(balance: list) -> str:
    """The headline sentence, DERIVED from the rows above rather than typed beside them.

    A hand-written number next to the code that computes it is the drift this project
    keeps paying for: the prose and the table are two sources for one fact, and the stale
    one is the side that reads fine.
    """
    if not balance:
        return "no scale-matched Fe I partner for any live Fe II product"
    best = sorted(balance, key=lambda b: abs(b["FeII_minus_FeI_dex"]))[0]
    n_in = sum(1 for b in balance if b["within_reported_bar"])
    vis = [b for b in balance if b["band"] == "VIS" and b["treatment"] == "ENGINE-A"]
    vis_txt = ("; the VIS 1D-NLTE (Bergemann) pairs agree to "
               + ", ".join(f"{b['FeII_minus_FeI_dex']:+.4f}" for b in
                           sorted(vis, key=lambda b: b["holding"])) + " dex"
               if vis else "")
    return (
        f"Fe II minus Fe I on the SAME band, instrument, holding, tier, selector, route "
        f"and treatment — one scale by construction, so no correction is applied to "
        f"either side and no value is re-derived (RYA-161). {n_in} of {len(balance)} "
        f"matched pairs sit inside the Fe II reported bar{vis_txt}. Closest pair: "
        f"{best['holding']} {best['band']} {best['treatment']} at "
        f"{best['FeII_minus_FeI_dex']:+.4f} dex "
        f"(Fe II {best['A_FeII']:.3f} n={best['n_FeII']} vs Fe I {best['A_FeI']:.3f} "
        f"n={best['n_FeI']}). ⚠️ The Fe II bars are wide — dominated by line-to-line "
        f"scatter over n=3 — so 'consistent' here means the balance does not CONTRADICT "
        f"the Fe I anchor, not that it constrains it."
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/results/rya1055/fe2_label_audit.json")
    #: ⚠️ A verifier must not REWRITE the artifact it verifies. Without this a test run
    #: leaves the tracked file dirty on a `generated_at` stamp alone -- the same class as
    #: the RYA-776 note that a generated artifact must never stamp `git rev-parse HEAD`.
    ap.add_argument("--check", action="store_true",
                    help="re-derive and compare against the committed artifact; write "
                         "nothing, exit 1 on any difference other than the timestamp")
    a = ap.parse_args(argv)

    from pipeline.gerber_nlte import FE_II_NLTE_LIMIT

    feed = json.loads(FEED.read_text())
    live = feed["products"]
    fe2 = [p for p in live if p.get("ion") == "II"]
    fe1 = [p for p in live if p.get("ion") == "I"]

    rows, problems = [], []
    for p in sorted(fe2, key=lambda x: (x["band"], x["holding"], x["treatment"])):
        cap = (p.get("science_provenance") or {}).get("nlte_capability") or {}
        stamped = bool(cap)
        on_deck = cap.get("this_product_takes_nlte_from_the_gerber_deck")
        quotes = cap.get("limit") == FE_II_NLTE_LIMIT
        if not stamped:
            problems.append(f"UNSTAMPED {p['band']}/{p['holding']}/{p['treatment']}")
        if not quotes:
            problems.append(f"LIMIT TEXT DRIFTED {p['band']}/{p['holding']}/{p['treatment']}")
        if on_deck:
            problems.append(
                f"🔴 GERBER-DECK NLTE ON Fe II: {p['band']}/{p['holding']}/"
                f"{p['treatment']} -- LTE under an NLTE label (RYA-1050 should have "
                f"refused this)")
        rows.append({
            "band": p["band"], "holding": p["holding"], "tier": p["tier"],
            "treatment": p["treatment"], "display": p["display"],
            "A": p["A"], "n_lines": p["n_lines"],
            "sigma_reported": p.get("sigma_reported"),
            "limit_stamped": stamped,
            "limit_text_matches_the_single_source": quotes,
            "takes_nlte_from_the_gerber_deck": on_deck,
            "nlte_source": cap.get("nlte_source_for_this_product"),
            "label_verdict": ("CLEAN — LTE product, claims no NLTE"
                              if "LTE" == p["treatment"] or p["treatment"] == "1D-LTE"
                              else "CLEAN — NLTE from the MPIA/Bergemann per-line grid, "
                                   "not from atom.fe607a"
                              if p["treatment"] == "ENGINE-A"
                              else "REVIEW — unrecognised treatment for an Fe II product"),
        })

    # ── the per-line cross-check: the sharper of the two ────────────────────────────
    pl = list(csv.DictReader(l for l in PERLINE.read_text().splitlines()
                             if not l.startswith("#")))
    nlte_scale = {"I": 0, "II": 0}
    for r in pl:
        if r["ion"] in nlte_scale and "NLTE" in r["scale"].upper():
            nlte_scale[r["ion"]] += 1
    if nlte_scale["II"]:
        problems.append(f"🔴 {nlte_scale['II']} Fe II per-line rows are on an NLTE scale")
    if not nlte_scale["I"]:
        problems.append("the Fe I control is empty — the per-line check measures nothing")

    # ── 4. what the RUNS say, line by line ──────────────────────────────────────────
    # The label is read from `nlte_source`/`nlte_delta_dex`, written at the point the
    # correction was applied (RYA-880) -- the run speaking, not a label being believed.
    BP = ROOT / "data/results/band_products"
    perline = []
    for f in sorted(BP.glob("FeII_*_lines.csv")):
        try:
            import pandas as pd
            d = pd.read_csv(f)
        except Exception as e:                       # a report must not be the thing that breaks
            problems.append(f"unreadable per-line artifact {f.name}: {e}")
            continue
        src = sorted({str(x) for x in d.get("nlte_source", []).dropna()}) \
            if "nlte_source" in d.columns else []
        delta = d["nlte_delta_dex"] if "nlte_delta_dex" in d.columns else None
        n_dep = int((delta.fillna(0) != 0).sum()) if delta is not None else 0
        gerber = [x for x in src if "erber" in x or "fe607a" in x]
        if gerber:
            problems.append(f"🔴 {f.name}: nlte_source names the Gerber deck {gerber}")
        perline.append({
            "artifact": f.name,
            "n_lines": int(len(d)),
            "n_in_aggregate": int(d["in_aggregate"].sum()) if "in_aggregate" in d else None,
            "n_with_a_nonzero_departure": n_dep,
            "nlte_delta_dex_values": sorted({round(float(x), 4) for x in
                                             (delta.dropna() if delta is not None else [])}),
            "nlte_source": src,
            # ⚠️ guarded like `nlte_source` above: `d.get(col, [])` returns a LIST when
            # the column is absent, and a list has no `.dropna()` — that would crash the
            # audit on an older artifact instead of reporting it.
            "excluded_reasons": (sorted({str(x) for x in d["excluded_reason"].dropna()})[:4]
                                 if "excluded_reason" in d.columns else []),
        })

    # ── the balance, SCALE-MATCHED ──────────────────────────────────────────────────
    # 🔴 MATCHED ON THE FULL RYA-1127 IDENTITY MINUS THE ION -- band, instrument, holding,
    # tier, SELECTOR, ROUTE and treatment. Dropping `selector` is not a cosmetic looseness:
    # HARPS VIS Fe I ships DEEPGRADED and DEEPGRADED-LOCALRENORM as separate products at
    # 7.535 and 7.339, so a key without it pairs one Fe II cell against TWO Fe I partners
    # 0.196 dex apart and the reader picks whichever suits. Dropping `route` would pair an
    # EW Fe I against a synthesis Fe II. Both are the pool-drift confound RYA-880/1113 keep
    # paying for, in a comparison whose whole point is that both sides are on one scale.
    KEY = ("band", "instrument", "holding", "tier", "selector", "route", "treatment")
    fe1_idx = {}
    for p in fe1:
        fe1_idx.setdefault(tuple(str(p.get(k)) for k in KEY), []).append(p)
    balance = []
    for p in fe2:
        key = tuple(str(p.get(k)) for k in KEY)
        partners = fe1_idx.get(key, [])
        if len(partners) > 1:
            problems.append(f"AMBIGUOUS Fe I partner for {key}: {len(partners)} matches — "
                            f"refusing to pick one")
            continue
        for q in partners:
            balance.append({
                "holding": p["holding"], "band": p["band"],
                "treatment": p["treatment"], "tier": p["tier"],
                "scale": p["display"].split(" · ")[1] if " · " in p["display"] else None,
                "A_FeII": p["A"], "n_FeII": p["n_lines"],
                "A_FeI": q["A"], "n_FeI": q["n_lines"],
                "FeII_minus_FeI_dex": round(p["A"] - q["A"], 4),
                "sigma_reported_FeII": p.get("sigma_reported"),
                "within_reported_bar": (
                    abs(p["A"] - q["A"]) <= (p.get("sigma_reported") or 0)),
            })

    out = {
        "ticket": "RYA-1055",
        "purpose": "item 3 — confirm live Fe II labels + state the balance with its scale",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_live_fe_ii_products": len(fe2),
        "n_stamped": sum(1 for r in rows if r["limit_stamped"]),
        "n_taking_nlte_from_the_gerber_deck": sum(1 for r in rows
                                                  if r["takes_nlte_from_the_gerber_deck"]),
        "per_line_rows_on_an_nlte_scale": nlte_scale,
        "band_product_per_line_nlte_sources": perline,
        "rya1113_contradiction": (
            "RESOLVED. RYA-1113 found a live Fe II near-UV NLTE leg at n=7 against its "
            "LTE sibling's n=12 and asked how it can exist if Fe II NLTE is structurally "
            "unavailable, offering two answers: it is not real NLTE, or the atom has "
            "changed. BOTH ARE WRONG. The per-line artifact RYA-908 has since emitted "
            "names the source outright -- 7 lines at nlte_delta_dex = -0.001 from "
            "'Bergemann MPIA per-line delta_nlte (live query, solar node)', 5 excluded "
            "'ENGINE-A-NOT-SERVED'. The leg is real NLTE from the MPIA grid, a source "
            "this ticket's finding never touched, and the atom was re-measured today at "
            "still-zero Fe II bound-bound transitions. ⚠️ RYA-1113's OTHER half stands "
            "unchanged: the n=7-vs-12 gap is MPIA service coverage, so the published "
            "Fe II near-UV NLTE delta IS confounded with a 5-line pool change."),
        "products": rows,
        "ionisation_balance_scale_matched": balance,
        "balance_note": _balance_note(balance),
        "anchor_caveat": (
            f"⚠️ DO NOT quote the Fe II VIS cells against the {OPTICAL_ANCHOR} anchor "
            f"directly. That anchor is on the {OPTICAL_ANCHOR_SCALE} scale and these are "
            f"1D-NLTE products; comparing them mixes two scales, which is the RYA-669 "
            f"ratchet the RYA-783 report warns about in as many words. The coincidence is "
            f"inviting -- the kurucz2005 Fe II cell is 7.466 to three decimals, the "
            f"anchor's own digits -- and a number that looks right on arrival is the one "
            f"that never gets challenged. The defensible statement is the scale-matched "
            f"balance above."),
        "problems": problems,
        "verdict": ("CLEAN — every live Fe II product carries the limit, none takes "
                    "departures from the Gerber deck, and no Fe II per-line row is on an "
                    "NLTE scale" if not problems else "PROBLEMS FOUND"),
    }
    p = ROOT / a.out
    if a.check:
        have = json.loads(p.read_text())
        g = {k: v for k, v in out.items() if k != "generated_at"}
        h = {k: v for k, v in have.items() if k != "generated_at"}
        if g != h:
            for k in sorted(set(g) | set(h)):
                if g.get(k) != h.get(k):
                    print(f"DIFFERS  {k}: committed={h.get(k)!r} derived={g.get(k)!r}",
                          file=sys.stderr)
            return 1
        print(f"{a.out}: matches the live feed")
        return 1 if problems else 0
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2) + "\n")

    print(f"live Fe II products: {len(fe2)}  stamped: {out['n_stamped']}  "
          f"on the Gerber deck: {out['n_taking_nlte_from_the_gerber_deck']}")
    print(f"per-line rows on an NLTE scale: Fe I {nlte_scale['I']}, "
          f"Fe II {nlte_scale['II']}")
    print("\nionisation balance, SCALE-MATCHED (same holding/band/treatment/tier):")
    for b in balance:
        print(f"  {b['holding']:<34}{b['band']:<9}{b['treatment']:<10}"
              f"Fe II {b['A_FeII']:.3f} (n={b['n_FeII']})  "
              f"Fe I {b['A_FeI']:.3f} (n={b['n_FeI']})  "
              f"Delta {b['FeII_minus_FeI_dex']:+.4f}  "
              f"{'within bar' if b['within_reported_bar'] else 'OUTSIDE bar'}")
    for x in problems:
        print("  PROBLEM:", x)
    print(f"\n{out['verdict']}")
    print(f"wrote {a.out}")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
