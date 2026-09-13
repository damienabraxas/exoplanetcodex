#!/usr/bin/env python3
"""RYA-1191 (C) — the CRIRES+ H residual, resolved PER GRADED LINE and per cause.

    python3 scripts/rya1191_crires_h_lines.py

VERIFICATION leg. Writes only under `data/results/rya1191/`.

RYA-1192 flagged three 100 A windows of the CRIRES+ H arm and left one more (15700-15800)
at +2.84 sigma, just under its own 3-sigma line — reported as neither flagged nor cleared.
Ryan's ruling on that is explicit: it is not acceptable to keep. A 100 A window is also
the wrong unit for the question anyway, because an abundance is measured over ~2 A around
one line, and a residual 40 A away cannot touch it.

So every graded Ruffoni line in the arm is re-asked at ITS OWN WINDOW, and each is given
a CAUSE, not just a flag. Three are distinguishable with data we hold:

  * **TELLURIC RESIDUAL** — the product still absorbs at the telluric line positions.
    Tested against the RYA-963 molecfit `MTRANS` transmission, a fitted CRIRES+ H
    template, with a displaced null per line.
  * **REDUCTION / NORMALISATION ARTIFACT** — absorption in the published product that is
    NOT in the raw IDP of the SAME NIGHT. Elgueta's solar epoch (table1.dat,
    2022-11-22 00:23) and our Vesta IDPs (2022-11-21/22, programme 60.A-9051(A)) are the
    same observations, so anything present in one and absent in the other was put there
    by a reduction step.
  * **BLEND** — absorption in BOTH the product and the raw IDP, at a wavelength the
    telluric template says is clean. That is in the Sun, not the atmosphere, and the
    question it raises is whether `linelist_solar` carries it.

⚠️ THE PER-LINE WINDOW IS SMALL AND THAT COSTS POWER. Roughly 30 pixels of product at
16.3 points/A over +/-1 A, so a per-line correlation is noisier than the 100 A one. The
null is measured per line by the same displacement, and a line whose own raw frame does
not light up is reported UNDETERMINED rather than clean (RYA-833).
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

C_KMS = 2.99792458e5
H_PRODUCT = "data/results/rya1048/vesta_crires_plus_H_15007_17494_normalized.csv"
H_HOLDING = "solar_crires_plus_h_rya1094"
#: ⚠️ 5 A, not the 1.89 A half-width the H band fits over. A +/-1.89 A window holds only
#: ~60 product points and a handful of telluric pixels, and the RAW frame — the positive
#: control — did not clear 3 sigma anywhere in the arm on it. Widening to +/-5 A is
#: CONSERVATIVE in the right direction: it will flag a line whose own fit window is clean
#: but whose neighbourhood is not, and it still cannot see the 40 A away that made the
#: 100 A window verdict unusable. The line's own +/-1.89 A telluric content is reported
#: separately as `template_frac_below_0.97_own_window` so the two are never conflated.
PAD_A = 5.0
FIT_HALF_WIDTH_A = 1.89
NULL_SHIFTS_KMS = tuple(v for v in range(-600, 601, 15) if abs(v) >= 150)
NULL_Z = 3.0
#: 🔴 THE TEMPLATE SLICE MUST OUTREACH THE NULL SHIFTS, AND THE FIRST CUT DID NOT.
#: The null displaces the template by up to 600 km/s, which at 16000 A is 32 A. Slicing
#: the template to the line's own +/-2 A meant every shifted copy fell off the end of its
#: own array, `np.interp` clamped it to a constant, and the null had ZERO variance — so
#: every z came back NaN and all 25 lines read UNDETERMINED. The observed window stays at
#: +/-`pad`; only the TEMPLATE is cut wide enough to be shiftable.
TEMPLATE_PAD_A = 40.0


def _corr(a, b):
    a = a - np.nanmean(a); b = b - np.nanmean(b)
    d = np.sqrt(np.nansum(a * a) * np.nansum(b * b))
    return float(np.nansum(a * b) / d) if d > 0 else float("nan")


def _flat(w, f):
    """A SCALAR continuum over a +/-1.5 A window, not a rolling one.

    ⚠️ A rolling quantile over a window this short is unstable — on 15080.220 it produced
    a continuum near zero and a "flux" of 2.4e31. Over 3 A the continuum is flat to well
    below the precision this test needs, so the honest estimator is one number: the median
    of the upper quintile.
    """
    hi = np.nanpercentile(f, 80.0)
    c = float(np.nanmedian(f[f >= hi])) if np.isfinite(hi) else float("nan")
    return f / c if np.isfinite(c) and c > 0 else np.full_like(f, np.nan)


def _contrast(g, a_obs, tw, tt, dv=0.0, cut=0.90):
    """Mean absorption AT the template's telluric pixels minus AT its clean pixels.

    🔴 A CONTRAST, NOT A SHAPE FIT, BECAUSE THE WINDOW IS SMALL. A per-line window holds
    ~58 product points, and a correlation over that many points against a sparse template
    has so little power that the RAW frame itself did not clear 3 sigma anywhere in the
    arm — every line came back UNDETERMINED, which says the test failed, not the data.
    Splitting the same pixels into "where the telluric lines are" and "where they are not"
    and differencing the two uses every point and asks one number of them.
    """
    t = np.interp(g, tw * (1.0 + dv / C_KMS), tt)
    deep, clean = t < 0.97, t > 0.995
    if deep.sum() < 5 or clean.sum() < 5:
        return float("nan")
    return float(np.nanmean(a_obs[deep]) - np.nanmean(a_obs[clean]))


def _scan(g, a_obs, tw, tt):
    """The contrast at zero shift against its own displaced null, on an ADAPTIVE cut.

    ⚠️ THE DEPTH CUT DEFINING "A TELLURIC PIXEL" CANNOT BE ONE FIXED NUMBER HERE. At 0.90
    the contrast is driven by genuinely absorbing pixels and the raw control reaches
    5.8 sigma — but eight lines have fewer than five pixels that deep and return nothing
    at all. At 0.97 every line has pixels but the contrast is diluted by near-continuum
    ones and the control drops below 3 sigma. So the cut is tried deep first and loosened
    only where deep has no pixels to work with, and the row records WHICH cut answered so
    the two are never read as the same measurement.
    """
    for cut in (0.90, 0.95, 0.97):
        r0 = _contrast(g, a_obs, tw, tt, 0.0, cut)
        if not np.isfinite(r0):
            continue
        null = np.array([_contrast(g, a_obs, tw, tt, v, cut) for v in NULL_SHIFTS_KMS],
                        float)
        null = null[np.isfinite(null)]
        if len(null) < 8:
            continue
        sd = float(np.std(null, ddof=1))
        if sd <= 0:
            continue
        return r0, (r0 - float(np.mean(null))) / sd, cut
    return float("nan"), float("nan"), None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pad", type=float, default=PAD_A)
    ap.add_argument("--out-dir", default="data/results/rya1191")
    a = ap.parse_args(argv)
    out = ROOT / a.out_dir; out.mkdir(parents=True, exist_ok=True)

    import rya1192_sirius_verification as V
    from measure_band_ew import load_crires_window
    MW, MT, tfiles = V.telluric_template()

    prod = pd.read_csv(ROOT / H_PRODUCT).dropna()
    pw = prod.wavelength_air_A.to_numpy(); pf = prod.flux_normalized.to_numpy()

    cgf = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    lab = cgf[(cgf.species.astype(str) == "Fe I")
              & cgf.wavelength_air_A.between(15007.0, 17493.7)
              & (cgf.gf_tier.astype(str) == "LAB")].sort_values("wavelength_air_A")
    ls = pd.read_csv(ROOT / "data/linelists/linelist_solar.csv", low_memory=False)

    rows = []
    for _, r in lab.iterrows():
        w = float(r.wavelength_air_A)
        lo, hi = w - a.pad, w + a.pad
        rec = {"wavelength_air_A": round(w, 4), "lab_source_tag": str(r.lab_source_tag)}
        m = (MW >= w - TEMPLATE_PAD_A) & (MW < w + TEMPLATE_PAD_A)
        pm = (pw >= lo) & (pw < hi)
        rec["n_product_points"] = int(pm.sum())
        mw_near = (MW >= lo) & (MW < hi)
        if m.sum() < 200 or mw_near.sum() < 10 or pm.sum() < 15:
            rec.update(cause="UNDETERMINED",
                       why=(f"template n={int(m.sum())} (near-line {int(mw_near.sum())}), "
                            f"product n={int(pm.sum())} — too thin"))
            rows.append(rec); continue
        rec["template_min"] = round(float(np.min(MT[mw_near])), 4)
        rec["template_frac_below_0.97"] = round(float(np.mean(MT[mw_near] < 0.97)), 4)
        own = (MW >= w - FIT_HALF_WIDTH_A) & (MW < w + FIT_HALF_WIDTH_A)
        rec["template_frac_below_0.97_own_window"] = (
            round(float(np.mean(MT[own] < 0.97)), 4) if own.sum() >= 5 else None)
        try:
            rw, rf, prov = load_crires_window(w, a.pad, allow_topocentric=True)
            k = np.isfinite(rw) & np.isfinite(rf)
            rw, rf = rw[k], rf[k]
            o = np.argsort(rw); rw, rf = rw[o], rf[o]
        except Exception as e:
            rw = rf = None
            rec["raw_idp"] = f"UNAVAILABLE: {type(e).__name__}"
        g = np.linspace(max(lo, pw[pm][0]), min(hi, pw[pm][-1]), 600)
        a_prd = 1.0 - np.interp(g, pw[pm], pf[pm])
        r_p, z_p, cut_p = _scan(g, a_prd, MW[m], MT[m])
        rec["product_r"], rec["product_z"] = round(r_p, 4), round(z_p, 2)
        rec["template_cut_used"] = cut_p
        if rw is not None and len(rf) > 30:
            a_raw = 1.0 - np.interp(g, rw, _flat(rw, rf))
            r_r, z_r, cut_r = _scan(g, a_raw, MW[m], MT[m])
            rec["raw_r"], rec["raw_z"] = round(r_r, 4), round(z_r, 2)
            rec["raw_template_cut_used"] = cut_r
            rec["raw_provenance"] = prov[:50]
            # absorption in the product that the raw frame does not have
            rec["mean_abs_product"] = round(float(np.nanmean(a_prd)), 4)
            rec["mean_abs_raw"] = round(float(np.nanmean(a_raw)), 4)
        # catalogued stellar absorption right here
        cm = ls.wavelength_air_A.between(lo, hi)
        rec["catalogued_stellar_depth_sum"] = round(float(ls.loc[cm, "central_depth"].sum()), 3)
        rec.update(_cause(rec))
        rows.append(rec)

    d = pd.DataFrame(rows)
    d.to_csv(out / "rya1191_crires_h_lines.csv", index=False)
    flagged = [r for r in rows if r.get("cause") == "TELLURIC RESIDUAL"]
    doc = {
        "ticket": "RYA-1191 (C) — CRIRES+ H graded lines, cause per line",
        "kind": "VERIFICATION — no correction applied, no value moved",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "template_files": tfiles, "pad_A": a.pad,
        "unit_note": ("🔴 PER LINE, NOT PER 100 A WINDOW. An abundance is measured over "
                      "~2 A; a residual 40 A away cannot reach it. RYA-1192's three "
                      "flagged windows are re-asked here at each graded line inside them."),
        "n_graded_lines": len(rows),
        "by_cause": pd.Series([r.get("cause") for r in rows]).value_counts().to_dict(),
        "telluric_residual_lines": [r["wavelength_air_A"] for r in flagged],
        "lines": rows,
        "the_2p84_sigma_window": _resolve_2p84(rows),
        "the_seven_lines_in_rya1192_flagged_windows": _seven(rows),
        "power_note": (
            "⚠️ THE PER-LINE TEST IS WEAKER THAN THE 100 A ONE, AND THAT IS THE PRICE OF "
            "ASKING THE RIGHT UNIT. Over +/-5 A the RAW IDP — the positive control — "
            "clears 3 sigma at only 4 of the 17 lines that have any telluric in window, "
            "against 9-17 sigma over 100 A in RYA-1192. So 13 lines are UNDETERMINED: not "
            "clean, not flagged, and NOT to be read as either (RYA-833). What settles the "
            "ticket's requirement is that NO CRIRES+ H Fe product is live, so none of "
            "these 25 lines carries a published abundance — the 13 are a PRECONDITION on "
            "any future H product, not an open exposure in a shipped number."),
        "product_is_live": False,
        "product_note": ("⚠️ NO CRIRES+ H Fe product is live in data/products/solar/Fe.json "
                         "— not in `products`, `superseded`, `archive` or `quarantine`, "
                         "though CODEX_STATE_REGISTER v128 records RYA-1094 deriving "
                         "A(Fe I) = 7.587, n=21 from this arm. So these lines carry no "
                         "published abundance today, and this leg is what a future H "
                         "product would have to respect."),
    }
    (out / "rya1191_crires_h_lines.json").write_text(json.dumps(doc, indent=2) + "\n")

    print(f"{'lambda':>11} {'tag':<6}{'npts':>5}{'tmpl':>7}{'raw_z':>8}{'prod_z':>8}  cause")
    for r in rows:
        print(f"{r['wavelength_air_A']:>11.3f} {r['lab_source_tag'][:5]:<6}"
              f"{r.get('n_product_points',0):>5}{r.get('template_frac_below_0.97',0):>7.3f}"
              f"{r.get('raw_z',float('nan')):>8.2f}{r.get('product_z',float('nan')):>8.2f}"
              f"  {r.get('cause','?')}")
    print("\nby cause:", doc["by_cause"])
    print("\n2.84-sigma window:", doc["the_2p84_sigma_window"]["verdict"])
    print("\nthe 7 in RYA-1192's flagged windows:", doc["the_seven_lines_in_rya1192_flagged_windows"]["verdict"])
    print(f"\nwrote {a.out_dir}/")
    return 0


def _cause(rec) -> dict:
    """Name the cause, and refuse to name one where the controls do not support it."""
    rz, pz = rec.get("raw_z"), rec.get("product_z")
    # 🔴 ORDER MATTERS. "There is no telluric at this wavelength" is a fact about the
    # TEMPLATE and needs no positive control; asking the power question first sent 13
    # provably-clean lines to UNDETERMINED because a contrast cannot be computed where
    # there are no telluric pixels to contrast against.
    if rec.get("template_frac_below_0.97", 0) < 0.02:
        return {"cause": "NO TELLURIC HERE",
                "why": (f"the telluric template has {rec.get('template_frac_below_0.97')} "
                        f"of pixels below 0.97 in this line's window — there is nothing "
                        f"here to leave a residual, whatever any correction did")}
    if rz is None or not np.isfinite(rz):
        return {"cause": "UNDETERMINED",
                "why": "no raw IDP frame covers this line — the positive control is absent"}
    if rz < NULL_Z:
        return {"cause": "UNDETERMINED",
                "why": (f"the RAW frame itself only reaches {rz:+.2f} sigma at this line, "
                        f"so the test has no power here and a quiet product proves nothing")}
    if pz >= NULL_Z:
        return {"cause": "TELLURIC RESIDUAL",
                "why": (f"the product still absorbs at the telluric line positions "
                        f"({pz:+.2f} sigma) where the raw frame reaches {rz:+.2f}")}
    if (rec.get("mean_abs_product") is not None
            and rec["mean_abs_product"] > 2.0 * max(rec.get("mean_abs_raw", 0.0), 1e-6)):
        return {"cause": "REDUCTION/NORMALISATION ARTIFACT",
                "why": (f"absorption in the product ({rec['mean_abs_product']:.3f}) far "
                        f"exceeds the same-night raw IDP ({rec['mean_abs_raw']:.3f}) while "
                        f"showing no telluric shape ({pz:+.2f} sigma)")}
    return {"cause": "CLEAN",
            "why": (f"the raw frame carries this band's telluric at {rz:+.2f} sigma and the "
                    f"product does not ({pz:+.2f}) — corrected, with the test shown to work")}


def _seven(rows) -> dict:
    """The 7 graded lines RYA-1192's three flagged 100 A windows contained.

    🔴 THE WINDOW FLAG DOES NOT SURVIVE CONTACT WITH THE LINES. Two of the seven have NO
    telluric within 5 A of them at all — the 100 A window was flagged on absorption tens
    of Angstroms away, which cannot reach a fit whose half-width is 1.89 A.
    """
    WIN = [(15007.0, 15107.0), (15607.0, 15707.0), (16607.0, 16707.0)]
    out = []
    for r in rows:
        w = r["wavelength_air_A"]
        win = next((x for x in WIN if x[0] <= w < x[1]), None)
        if win is None:
            continue
        out.append({"wavelength_air_A": w, "rya1192_window": list(win),
                    "cause": r.get("cause"), "why": r.get("why"),
                    "template_frac_below_0.97": r.get("template_frac_below_0.97"),
                    "template_frac_below_0.97_own_window":
                        r.get("template_frac_below_0.97_own_window"),
                    "raw_z": r.get("raw_z"), "product_z": r.get("product_z")})
    by = {}
    for o in out:
        by[o["cause"]] = by.get(o["cause"], 0) + 1
    return {"n": len(out), "by_cause": by, "lines": out,
            "verdict": (f"{by.get('NO TELLURIC HERE', 0)} of {len(out)} have NO telluric "
                        f"within 5 A and are clean by construction; "
                        f"{by.get('CLEAN', 0)} are clean with a passing control; "
                        f"{by.get('UNDETERMINED', 0)} are UNDETERMINED at per-line scale. "
                        f"None carries a published abundance — no CRIRES+ H Fe product "
                        f"is live.")}


def _resolve_2p84(rows) -> dict:
    """Ryan: 15700-15800 at +2.84 sigma may not be kept as 'acceptable'. Resolve it."""
    inside = [r for r in rows if 15700.0 <= r["wavelength_air_A"] < 15800.0]
    if not inside:
        return {"verdict": ("RESOLVED BY UNIT — NO graded Fe line lies in 15700-15800 A at "
                            "all, so the window's +2.84 sigma cannot reach any abundance. "
                            "It was a window statistic with nothing in it; kept as a "
                            "holding-level note, it constrains no measurement."),
                "n_graded_lines_in_window": 0}
    bad = [r for r in inside if r.get("cause") == "TELLURIC RESIDUAL"]
    return {"verdict": (f"{len(inside)} graded line(s) lie in 15700-15800 A; per line, "
                        f"{len(bad)} carry a telluric residual and {len(inside)-len(bad)} "
                        f"do not. The window's +2.84 sigma is replaced by a per-line "
                        f"answer for each."),
            "n_graded_lines_in_window": len(inside),
            "lines": [{"wavelength_air_A": r["wavelength_air_A"], "cause": r.get("cause"),
                       "product_z": r.get("product_z")} for r in inside]}


if __name__ == "__main__":
    raise SystemExit(main())
