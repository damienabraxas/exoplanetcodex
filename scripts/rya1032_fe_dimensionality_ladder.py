"""
scripts/rya1032_fe_dimensionality_ladder.py
===========================================
RYA-1032 — is the solar Fe I 1D→⟨3D⟩ climb DIMENSIONALITY, or MODEL FAMILY?

THE QUESTION. Two 1D-NLTE model families (Bergemann, Gerber) disagree on solar
A(Fe I), and the ⟨3D⟩ route sits higher still. If the ⟨3D⟩ climb were the same
size as the family disagreement, "⟨3D⟩" would be telling us nothing that choosing
a different 1D atom does not already tell us. It is not: the atmosphere step is
~2.4x the family spread, and it survives de-confounding.

WHAT THIS TICKET NO LONGER DOES. RYA-1032 was written to *solve the Fe departures
ourselves*, on the stated fact that "there is no Fe ⟨3D⟩ departure deck, anywhere."
That fact was a disk scan reported as a fact about the source, and it is wrong:
`NLTEgrid4TS_Fe_STAGGERmean3D_May-21-2021.bin` (93 MB) has been on the MPG Keeper
share since 2021. RYA-1035 found it, RYA-710 (PR #385) wired it, and model 6 has
been live since 2026-08-25. So no solver is built here — the deck is consumed.

EVERY NUMBER IS DERIVED, NOT TYPED. Values come from the live published feed
`data/products/solar/Fe.json`; the axes (scale / model family / atmosphere) come
from `data/catalog/model_registry.csv`. Nothing is hardcoded but the token names.

🔴 WHAT THIS DELIBERATELY REFUSES TO DO. It will not compute the ⟨3D⟩ NLTE effect
(model 6 − model 5) by differencing the two published cells. Those cells are
MEDIANS, and on this holding the two medians collided on exactly 7.552 — so the
subtraction returns 0.000 for a real per-line shift of +0.032. That is RYA-1099
Finding 3, and the correct statistic already exists in
`pipeline.paired_differential` (RYA-1083). Reimplementing it here would be a
second implementation of one thing. This script reports the collision and points
at the artifact instead.

Usage
-----
    python scripts/rya1032_fe_dimensionality_ladder.py
    python scripts/rya1032_fe_dimensionality_ladder.py --holding solar_kpno_kurucz2005_corrected
    python scripts/rya1032_fe_dimensionality_ladder.py --json out.json
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
FEED = REPO / "data" / "products" / "solar" / "Fe.json"
REGISTRY = REPO / "data" / "catalog" / "model_registry.csv"

DEFAULT_HOLDING = "solar_kpno_molecfit_corrected"

# The rungs, by (treatment token, route). Treatment alone is ambiguous: "1D-LTE"
# and "ENGINE-A" each appear on both the SYNTH and PROFILEFIT routes.
LADDER = [
    ("1D-LTE", "SYNTH"),
    ("synth-1D-LTE-gerber", "SYNTH"),
    ("ENGINE-A", "SYNTH"),
    ("ENGINE-B-NLTE", "SYNTH"),
    ("synth-mean3D-LTE-gerber-stagger", "SYNTH"),
    ("synth-mean3D-NLTE-gerber-stagger", "SYNTH"),
    ("ENGINE-A-3DNLTE", "EW-3D"),
]


class LadderError(RuntimeError):
    """A rung is missing, or a comparison does not hold the axis it claims."""


def load_registry():
    rows = {}
    with open(REGISTRY, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            tok = (r.get("stored_token") or "").strip()
            if tok and tok != "-":
                rows[tok] = r
    return rows


def load_rungs(holding):
    feed = json.loads(FEED.read_text(encoding="utf-8"))
    want = {k: None for k in LADDER}
    for p in feed["products"]:
        if (p.get("element") == "Fe" and p.get("ion") == "I"
                and p.get("band") == "VIS" and p.get("tier") == "GRADED"
                and p.get("holding") == holding):
            key = (p.get("treatment"), p.get("route"))
            if key in want:
                if want[key] is not None:
                    raise LadderError(f"duplicate rung {key} on holding {holding}")
                want[key] = p
    missing = [k for k, v in want.items() if v is None]
    if missing:
        raise LadderError(
            f"missing rungs on holding {holding}: {missing}. "
            "Refusing to report a partial ladder."
        )
    return feed, want


def _split_scale(scale):
    """`scale` is TWO axes wedged into one string -- split them.

    The registry writes 1D-LTE / 1D-NLTE / <3D>-LTE / <3D>-NLTE / 3D-NLTE. That is
    a dimensionality axis (1D, <3D>, 3D) crossed with an LTE/NLTE axis. Asserting
    on `scale` as a single thing is what let an early version of this script claim
    to hold "scale" fixed across 1D-LTE vs <3D>-LTE, which holds nothing.
    """
    head, _, tail = scale.rpartition("-")
    if tail not in ("LTE", "NLTE") or not head:
        raise LadderError(f"cannot split scale {scale!r} into (dim, nlte)")
    return head, tail


def _axes(reg, token):
    r = reg.get(token)
    if r is None:
        raise LadderError(f"token {token!r} is not in model_registry.csv")
    if r["status"] != "live":
        raise LadderError(f"token {token!r} is status={r['status']}, not live")
    dim, nlte = _split_scale(r["scale"])
    return {"scale": r["scale"], "dim": dim, "nlte": nlte,
            "family": r["model_family"], "atmosphere": r["atmosphere"]}


def step(reg, rungs, hi_tok, lo_tok, *, route="SYNTH",
         fixed=(), varies=None, label=""):
    """A difference between two rungs, with the axis contract ASSERTED.

    `fixed` names axes that must be identical on both sides; `varies` names the
    one axis that must differ. A comparison whose axes do not hold is not a
    weaker result, it is a different quantity -- so this raises rather than
    returning a number that would be quietly mislabelled (RYA-542).
    """
    hi, lo = rungs[(hi_tok, route)], rungs[(lo_tok, route)]
    ah, al = _axes(reg, hi_tok), _axes(reg, lo_tok)
    for ax in fixed:
        if ah[ax] != al[ax]:
            raise LadderError(
                f"{label}: {ax} must be FIXED but differs "
                f"({ah[ax]!r} vs {al[ax]!r})")
    if varies and ah[varies] == al[varies]:
        raise LadderError(
            f"{label}: {varies} must VARY but is identical ({ah[varies]!r}) "
            "-- this comparison measures nothing")
    return {
        "label": label,
        "minuend": hi_tok, "subtrahend": lo_tok,
        "A_hi": hi["A"], "A_lo": lo["A"],
        "delta": round(hi["A"] - lo["A"], 4),
        "n_hi": hi["n_lines"], "n_lo": lo["n_lines"],
        "same_pool": hi["n_lines"] == lo["n_lines"],
        "held_fixed": list(fixed), "varied": varies,
        "axes_hi": ah, "axes_lo": al,
    }


def analyse(holding=DEFAULT_HOLDING):
    reg = load_registry()
    feed, rungs = load_rungs(holding)

    # (1) DIMENSIONALITY: family and LTE-ness held, dimensionality varies 1D -> <3D>
    #     (and with it the atmosphere, which IS the dimensionality change here).
    dim = step(reg, rungs, "synth-mean3D-LTE-gerber-stagger", "synth-1D-LTE-gerber",
               fixed=("family", "nlte"), varies="dim",
               label="atmosphere step (1D-LTE -> <3D>-LTE, family fixed)")

    # (2) The 1D NLTE step: family, dimensionality and atmosphere all held.
    nlte1d = step(reg, rungs, "ENGINE-B-NLTE", "synth-1D-LTE-gerber",
                  fixed=("family", "atmosphere", "dim"), varies="nlte",
                  label="1D NLTE step (Gerber, marcs-ges)")

    # (3) MODEL FAMILY: dimensionality and LTE-ness held at 1D-NLTE, family varies.
    #     NOTE this also crosses atlas9 -> marcs-ges; step (4) measures that nuisance.
    fam = step(reg, rungs, "ENGINE-B-NLTE", "ENGINE-A",
               fixed=("dim", "nlte"), varies="family",
               label="model-family spread (Gerber - Bergemann, 1D-NLTE)")

    # (4) 🔴 THE DE-CONFOUNDING CONTROL. The family spread in (3) compares
    #     bergemann@atlas9 against gerber@marcs-ges, so it carries an atmosphere
    #     difference too. Measure that atmosphere difference on its own, in LTE,
    #     where no NLTE physics can contribute. If it is small, (3) is genuinely
    #     family. Without this control the headline ratio is not defensible.
    #     `family` is deliberately NOT held: in LTE there is no NLTE model family,
    #     so the registry's none-vs-gerber label carries no physics here -- the
    #     only thing that differs in substance is the atmosphere.
    nuisance = step(reg, rungs, "synth-1D-LTE-gerber", "1D-LTE",
                    fixed=("dim", "nlte"), varies="atmosphere",
                    label="atlas9 -> marcs-ges nuisance (LTE, no NLTE physics)")

    fam_deconf = round(fam["delta"] - nuisance["delta"], 4)

    # (5) The overshoot past the full-3D reference.
    m3 = rungs[("synth-mean3D-NLTE-gerber-stagger", "SYNTH")]
    am = rungs[("ENGINE-A-3DNLTE", "EW-3D")]
    overshoot = round(m3["A"] - am["A"], 4)

    # (6) The published <3D> NLTE effect -- reported, deliberately NOT computed.
    lte = rungs[("synth-mean3D-LTE-gerber-stagger", "SYNTH")]
    collision = (m3["A"] == lte["A"])

    return {
        "ticket": "RYA-1032",
        "feed_version": feed.get("version"),
        "feed_updated_at": feed.get("updated_at"),
        "holding": holding,
        "ladder": [
            {"token": t, "route": r, "A": rungs[(t, r)]["A"],
             "n_lines": rungs[(t, r)]["n_lines"],
             "display": rungs[(t, r)].get("display")}
            for (t, r) in LADDER
        ],
        "dimensionality_step": dim,
        "nlte_1d_step": nlte1d,
        "model_family_spread": fam,
        "atmosphere_nuisance": nuisance,
        "model_family_spread_deconfounded": fam_deconf,
        "ratio_dimensionality_over_family": round(dim["delta"] / fam["delta"], 2),
        "ratio_dimensionality_over_family_deconfounded": (
            round(dim["delta"] / fam_deconf, 2) if fam_deconf else None),
        "overshoot_vs_amarsi_full3d": overshoot,
        "ladder_monotonic_toward_amarsi": overshoot <= 0,
        "mean3d_nlte_effect": {
            "published_difference_of_medians": round(m3["A"] - lte["A"], 4),
            "median_collision": collision,
            "computed_here": False,
            "why": ("the published cells are MEDIANS; differencing them returns "
                    "0.000 for a real per-line shift (RYA-1099 Finding 3). The "
                    "correct statistic is the per-line paired differential in "
                    "pipeline.paired_differential (RYA-1083), emitted as "
                    "{stem}_{treatment}_nlte_effect.json on a re-derive."),
        },
    }


def render(a):
    out = []
    out.append(f"RYA-1032 — solar Fe I VIS, GRADED, holding={a['holding']}")
    out.append(f"feed data/products/solar/Fe.json v{a['feed_version']} "
               f"({a['feed_updated_at']})")
    out.append("")
    out.append(f"  {'A':<9}{'n':<6}model")
    for r in a["ladder"]:
        out.append(f"  {r['A']:<9}{r['n_lines']:<6}{r['display']}")
    out.append("")
    for k in ("dimensionality_step", "nlte_1d_step", "model_family_spread",
              "atmosphere_nuisance"):
        s = a[k]
        pool = "" if s["same_pool"] else f"  ⚠️ pools differ (n={s['n_hi']} vs {s['n_lo']})"
        out.append(f"  {s['delta']:+.4f}  {s['label']}{pool}")
    out.append("")
    out.append(f"  model-family spread, de-confounded : {a['model_family_spread_deconfounded']:+.4f}")
    out.append(f"  dimensionality / family            : {a['ratio_dimensionality_over_family']}x"
               f"  ({a['ratio_dimensionality_over_family_deconfounded']}x de-confounded)")
    out.append("")
    verdict = ("DIMENSIONALITY, not model family"
               if a["dimensionality_step"]["delta"] > 2 * abs(a["model_family_spread"]["delta"])
               else "NOT separable at this margin")
    out.append(f"  VERDICT: the 1D->(3D) Fe climb is {verdict}.")
    out.append("")
    o = a["overshoot_vs_amarsi_full3d"]
    out.append(f"  <3D>-NLTE minus Amarsi full-3D     : {o:+.4f}"
               f"   -> ladder is {'MONOTONIC' if a['ladder_monotonic_toward_amarsi'] else 'NON-MONOTONIC'}"
               " toward Amarsi")
    if not a["ladder_monotonic_toward_amarsi"]:
        out.append("           <3D> OVERSHOOTS full 3D-NLTE — the mean-3D approximation")
        out.append("           over-corrects; the full cube (RYA-1119) is the check.")
    out.append("")
    e = a["mean3d_nlte_effect"]
    out.append(f"  <3D> NLTE effect from the published cells: "
               f"{e['published_difference_of_medians']:+.4f}"
               + ("   🔴 MEDIAN COLLISION — not the effect" if e["median_collision"] else ""))
    out.append("           not computed here by design — see RYA-1083/1099.")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--holding", default=DEFAULT_HOLDING)
    ap.add_argument("--json", help="also write the analysis to this path")
    args = ap.parse_args()
    try:
        a = analyse(args.holding)
    except LadderError as e:
        print(f"FAIL: {e}", file=sys.stderr)
        return 2
    print(render(a))
    if args.json:
        pathlib.Path(args.json).write_text(json.dumps(a, indent=2) + "\n",
                                           encoding="utf-8")
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
