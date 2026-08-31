"""
scripts/rya1120_sigma_reported.py
=================================
RYA-282 §3 — assemble sigma_reported per product, and print the breakdown that IS the
peer-review artifact.

    sigma_reported = sqrt( sigma_SE^2 + sigma_params^2 )
    sigma_params   = sqrt( sum_p (|dA/dp| * delta_p)^2 )

Reads three things and invents none of them:
  * the live feed `data/products/solar/Fe.json` — A, sigma_stat (already a standard
    error on all 46, per RYA-1112), n_lines;
  * the RYA-1120 campaign outputs — per-pool dA/dxi, measured by perturb-and-re-derive
    and differenced PER LINE through `pipeline.paired_differential`;
  * `pipeline.uncertainty_stack` — delta_p, the STAR's own parameter uncertainties.

🔴 WHAT IT WILL NOT DO

  * It will not borrow a derivative. A product whose own pool was not measured is
    reported UNMEASURED — its bar is an honest FLOOR — never given RYA-1089's 62-line
    number. Measured so far: dA/dxi ranges 7x across pools, so one number for Fe is not
    an approximation, it is a different quantity.
  * It will not difference two aggregates. Line acceptance moves with xi (measured:
    317->318 and 257->259 on the two large DEEPGRADED pools), so the derivative is a
    per-line paired differential or it is nothing.
  * It will not charge full 3D a xi term (Ryan's ruling), and it records that exemption
    WITH the contrary measurement rather than silently.
  * It will not quietly drop Teff. delta_Teff = 1.0 K makes that term 0.000665 dex
    against a 1e-4 published rounding step; it is carried as a DECLARED BOUND with the
    arithmetic and an explicit invalidated_if. SOLAR ONLY.

Usage
-----
    python scripts/rya1120_sigma_reported.py --campaign ~/xi_campaign/out
    python scripts/rya1120_sigma_reported.py --campaign ~/xi_campaign/out --json out.json
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

FEED = REPO / "data" / "products" / "solar" / "Fe.json"
XI_SPAN_KMS = 0.20          # central difference: (xi+0.10) - (xi-0.10)
FULL_3D = {"ENGINE-A-3DNLTE"}


#: The artifact stem is `<band>_<instrument>_<holding>_<ROUTE>[_<SELECTOR>]_<treatment>`.
#: The selector is ABSENT on the default (all-lines) product, which is why a regex keyed
#: on `_GRADED_` silently missed every DEEPGRADED file and reported those products
#: UNMEASURED. That accident happened to prevent a wrong-pool derivative from attaching —
#: but a matching rule must not depend on luck, so it now anchors on the ROUTE, which is
#: always present, and treats the selector as optional.
_STEM = re.compile(r"_(?P<route>PROFILEFIT|SYNTH|EW3D)"
                   r"(?:_(?P<sel>GRADED|DEEPGRADED|CONSISTENT|FROMEW[A-Z-]*))?"
                   r"_(?P<treat>.+)_lines\.csv$")


def _treat(fname: str) -> str:
    m = _STEM.search(fname)
    return m.group("treat") if m else fname


def _route(fname: str) -> str:
    """The ROUTE the artifact was produced by — PROFILEFIT / SYNTH / EW3D.

    🔴 ROUTE IS PART OF A PRODUCT'S IDENTITY, and the feed proves it: it carries TWO live
    products under the same (ion, holding, tier, treatment) key, differing only by route —

        harps GRADED 1D-LTE  PROFILEFIT  n= 6
        harps GRADED 1D-LTE  SYNTH       n=67

    Those are different POOLS by an order of magnitude, so they cannot share a dA/dxi
    (RYA-1093), and RYA-1071 exists to separate the route systematic from the pool
    systematic precisely because the two are not interchangeable. Keying without the route
    let a PROFILEFIT derivative measured on 5-14 lines attach to a SYNTH product of 37-67.
    """
    m = _STEM.search(fname)
    return m.group("route") if m else "UNKNOWN"


def _selector(fname: str) -> str:
    """The SELECTOR the artifact was produced with — 'ALL' when unlabelled.

    Returned so the caller can refuse to attach a derivative whose selector does not
    match the product's tier. A pool is part of a product's identity (RYA-984).
    """
    m = _STEM.search(fname)
    return (m.group("sel") or "ALL") if m else "UNKNOWN"


def campaign_derivatives(root: pathlib.Path) -> dict:
    """{(ion, holding, tier, treatment): dA/dxi} from COMPLETE pairs only."""
    import pandas as pd
    from pipeline.paired_differential import paired_differential, UnpairableProducts

    legs: dict[str, dict[str, pathlib.Path]] = {}
    for d in sorted(root.iterdir()):
        if (d / "DONE").exists():
            base, xi = d.name.rsplit("_xi", 1)
            legs.setdefault(base, {})[xi] = d

    out = {}
    for base, l in legs.items():
        if {"0.90", "1.10"} - set(l):
            continue                       # half a pair is not a smaller measurement
        # 🔴 THE ROUTE IS IN THE DIRECTORY NAME because it is part of the unit. A GRADED
        # pool run through PROFILEFIT and through SYNTH are different products with pools
        # an order of magnitude apart (n=6 vs n=67), so they cannot share an output dir
        # and must not share a key.
        m = re.match(r"Fe(I{1,2})_(solar_.+)_(GRADED|DEEPGRADED)_(PROFILEFIT|SYNTH)$", base)
        if not m:
            continue
        ion, holding, tier = m.group(1), m.group(2), m.group(3)
        dir_route = m.group(4)
        for f in sorted(l["1.10"].glob("*_lines.csv")):
            g = l["0.90"] / f.name
            if not g.exists():
                continue
            hi, lo = pd.read_csv(f), pd.read_csv(g)
            try:
                pdiff = paired_differential(hi, lo)
            except UnpairableProducts:
                continue
            if _route(f.name) != dir_route:
                # the run emitted an artifact from a route the unit did not ask for --
                # refuse rather than record it under the unit's name
                print(f"  SKIP {f.name}: route {_route(f.name)} != unit route "
                      f"{dir_route} (refusing to attach)", file=sys.stderr)
                continue
            sel = _selector(f.name)
            if sel != tier:
                # 🔴 A derivative measured with one selector may not be charged to a
                # product built with another: they are different LINE SETS, and that is
                # the whole premise of the campaign. Skip loudly rather than attach.
                print(f"  SKIP {f.name}: selector {sel} != tier {tier} "
                      f"(different pool — refusing to attach)", file=sys.stderr)
                continue
            out[(ion, holding, tier, _route(f.name), _treat(f.name))] = {
                "dA_dxi": pdiff.median / XI_SPAN_KMS,
                "n_paired": pdiff.n_paired,
                "pool_moved": len(hi) != len(lo),
                "dA_dxi_from_aggregates": pdiff.difference_of_aggregates / XI_SPAN_KMS,
            }
    return out


def build(campaign_root: pathlib.Path) -> dict:
    from pipeline.stellar_params import for_product, FULL_3D_XI_EXEMPTION
    from pipeline.uncertainty_stack import params_and_deltas

    _, deltas = params_and_deltas("solar")
    delta_xi = float(deltas["vturb_kms"])
    derivs = campaign_derivatives(campaign_root)
    feed = json.loads(FEED.read_text(encoding="utf-8"))

    rows = []
    for p in feed["products"]:
        if p.get("band") != "VIS":
            continue
        key_exact = (p["ion"], p["holding"], p["tier"], p.get("route"), p["treatment"])
        # 🔴 EXACT MATCH ONLY -- NO FALLBACK TO "SOME DERIVATIVE FROM THIS POOL".
        # This previously fell back to the unit when the exact treatment key missed,
        # taking the single candidate if there was one. Two things make that unsound:
        #
        #   * dA/dxi is NOT shared across treatments. MEASURED on harps GRADED: ENGINE-B
        #     drifts -0.001 dex over the 0.20 km/s span where 1D-LTE on the SAME pool
        #     drifts -0.036 -- a factor of 36. Borrowing across treatments is not an
        #     approximation, it is a different quantity (the RYA-1093 argument, one axis
        #     over from the pool it was made about).
        #   * `derive_band_products` emits EVERY treatment for a pool in one run, and 8 of
        #     the 11 units carry a SUPERSEDED ENGINE-B in the same (ion, holding, tier)
        #     cell. So the candidate pool contains derivatives belonging to products that
        #     are not live, and the fallback could charge one to a live product.
        #
        # A product whose own treatment was not measured is UNMEASURED. That is the
        # honest answer and the bar stays a FLOOR.
        d = derivs.get(key_exact)

        if p["treatment"] in FULL_3D:
            state, sig_xi, note = "NOT_APPLICABLE", None, FULL_3D_XI_EXEMPTION
        elif d is None:
            state, sig_xi, note = "UNMEASURED", None, (
                "no perturb-and-re-derive on this pool; dA/dxi is a property of the "
                "LINE SET (RYA-1093) so no element-level number may stand in")
        else:
            state, sig_xi = "MEASURED", abs(d["dA_dxi"]) * delta_xi
            note = (f"|dA/dxi|={abs(d['dA_dxi']):.4f} x delta_xi={delta_xi:.4f} on "
                    f"{d['n_paired']} paired lines"
                    + ("  [POOL MOVED with xi]" if d["pool_moved"] else ""))

        sig_stat = float(p["sigma_stat"])
        sig_rep = math.hypot(sig_stat, sig_xi) if sig_xi is not None else None
        rows.append({
            "ion": p["ion"], "holding": p["holding"], "tier": p["tier"],
            "treatment": p["treatment"], "route": p.get("route"),
            "A": p["A"], "n_lines": p["n_lines"],
            "sigma_SE": sig_stat, "sigma_syst_published": float(p["sigma_syst"]),
            "dA_dxi": (d or {}).get("dA_dxi"),
            # 🔴 THE TRIPWIRE. The exact product key is the GATE, but three different
            # axes have now leaked a wrong-pool derivative past a key that looked right
            # (selector, treatment, route), and the thing that exposed every one of them
            # was the same: the derivative's pool did not match the product's. So the two
            # counts travel side by side on every row and no threshold is applied
            # (RYA-161: a cut must be derived, not picked). A reader seeing `7` against
            # `67` needs no rule to know something is wrong.
            "n_paired": (d or {}).get("n_paired"),
            "sigma_xi": sig_xi, "xi_state": state, "xi_note": note,
            "sigma_params": sig_xi,          # xi only; Teff is a declared bound
            "sigma_reported": sig_rep,
            "dominant": (None if sig_rep is None else
                         ("xi" if (sig_xi or 0) > sig_stat else "statistics")),
        })
    return {"delta_xi_kms": delta_xi, "n_products": len(rows), "products": rows}


def render(doc):
    L = [f"sigma_reported = sqrt(sigma_SE^2 + sigma_params^2)   delta_xi = "
         f"{doc['delta_xi_kms']:.4f} km/s   ({doc['n_products']} VIS products)", ""]
    L.append(f"{'ion':<4}{'tier':<12}{'route':<11}{'treatment':<26}{'n':>4}{'npair':>6}"
             f"{'sig_SE':>9}{'dA/dxi':>9}{'sig_xi':>9}{'sig_rep':>9}  dominant")
    for r in sorted(doc["products"], key=lambda r: (r["ion"], r["tier"], r["treatment"])):
        f = lambda v, w=9, p=4: (f"{v:>{w}.{p}f}" if isinstance(v, float) else f"{'-':>{w}}")
        np_ = r.get('n_paired')
        L.append(f"{r['ion']:<4}{r['tier']:<12}{str(r.get('route'))[:10]:<11}"
                 f"{r['treatment'][:25]:<26}"
                 f"{r['n_lines']:>4}{(str(np_) if np_ else '-'):>6}"
                 f"{f(r['sigma_SE'])}{f(r['dA_dxi'])}"
                 f"{f(r['sigma_xi'])}{f(r['sigma_reported'])}  "
                 f"{r['dominant'] or r['xi_state']}")
    n_m = sum(1 for r in doc["products"] if r["xi_state"] == "MEASURED")
    n_u = sum(1 for r in doc["products"] if r["xi_state"] == "UNMEASURED")
    n_n = sum(1 for r in doc["products"] if r["xi_state"] == "NOT_APPLICABLE")
    L += ["", f"MEASURED {n_m}   UNMEASURED {n_u} (bar is a FLOOR)   NOT APPLICABLE {n_n}"]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--json")
    a = ap.parse_args()
    doc = build(pathlib.Path(a.campaign).expanduser())
    print(render(doc))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(doc, indent=2) + "\n")
        print(f"\nwrote {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
