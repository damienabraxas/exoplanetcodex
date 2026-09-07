#!/usr/bin/env python3
"""RYA-1191 (B) — is a "raw-like" spot a telluric RESIDUAL or an uncatalogued BLEND?

    python3 scripts/rya1191_telluric_residual.py

VERIFICATION leg. Writes only under `data/results/rya1191/`.

🔴 WHY A THIRD TEST, WHEN RYA-1192 ALREADY HAD TWO
---------------------------------------------------
RYA-1192 judged a band on one of two statistics, and BOTH are blind in a way this ticket
walked straight into:

* **DEPTH** (fraction of pixels below 0.8) only sees absorption deep enough to cross the
  cut. The O2 gamma band at 6270-6300 is shallow, so every holding looked flat there --
  Kitt Peak RAW sits 1.07x its own side windows -- and the band was written off as
  "the anchor shows no local telluric". It has telluric in it.
* **DIFFERENCE** (corrected vs raw, against a reduction floor) cannot tell a LINE REMOVAL
  from a CONTINUUM RESCALE. Measured here on the RYA-940 o2gamma product: corrected/raw
  has **min 1.0000 and NOT ONE pixel below 0.98** over the whole 6270-6300 window -- a
  featureless ~3.6% brightening. RYA-1192 read its large mean |diff| as a correction and
  returned VERIFIED-CORRECTED. Nothing was removed.

Neither error is a threshold that needed tuning; both are the wrong QUESTION. The right
one is whether the flux still absorbs AT THE WAVELENGTHS WHERE THE TELLURIC LINES ARE.

THE TEMPLATE, AND WHY THIS ONE
-------------------------------
`solar_iag_reiners2016` (uncorrected) and `solar_iag` (Baker+2020, corrected) are the SAME
IAG FTS atlas, so their RATIO is the telluric transmission that correction removed --
measured, not modelled, on a real solar spectrum at 5001-11083 A. That covers every
registered telluric band in VIS, red-optical and the near IR.

⚠️ ITS DEPTHS ARE NOT TRANSFERABLE, ONLY ITS LINE POSITIONS. Airmass and precipitable
water differ between IAG's night and Kitt Peak's, so a holding is scored by SHAPE
CORRELATION against the template, never by matching depth.

⚠️ AND THE SIGN IS EASY TO READ BACKWARDS. The template is raw/corrected, so it dips
BELOW 1 where telluric was removed. A holding that still carries telluric therefore
correlates NEGATIVELY with it. `residual_r` is reported already flipped, so that a LARGER
POSITIVE number always means MORE RESIDUAL, and the raw references are printed beside
every row as the positive control.

🔴 THE SECOND NULL, AND LEG B SHIPPED WITHOUT IT. The displaced-template null below
controls for REGISTRATION — whether a correlation survives when the template is slid off
the lines. It does NOT control for the template carrying structure that is not telluric at
all. Measured on windows with no telluric complex in them, the raw/corrected ratio of two
SEPARATELY REDUCED IAG products is not flat:

    5500-5600  clean   template min 0.985   frac<0.95 0.0000
    6000-6100  clean   template min 0.946   frac<0.95 0.0007
    6400-6500  clean   template min 0.476   frac<0.95 0.0566   <-- in a CLEAN window
    6270-6300  O2 gamma             0.172             0.1313
    6867-6884  O2 B                 0.005             0.7435
    7160-7340  H2O                  0.005             0.5092

⚠️ SO THE O2 GAMMA VERDICT WAS NOT SAFE AND IS WITHDRAWN. Its template depth is barely
above the 6400-6500 CONTROL, and molecfit's own physical model of that band — fitted on
HARPS's own night with real GDAS — puts the deepest O2 gamma pixel at **1.78% absorption**
against this template's 82.76%, a factor of 47, with only r = +0.21 shared between them.
The "+12.0 sigma, essentially uncorrected" reading was measuring HARPS against a
Baker-vs-Reiners REDUCTION difference, not against telluric.

The DEEP bands are unaffected: O2 B, H2O 7160-7340 and H2O 9280-9600 sit one to two orders
above the control level, so the template there is dominated by the real correction. A band
is now only judged if it clears the clean-window null, and o2gamma does not.

TELLURIC OR BLEND — THE DISCRIMINATOR
--------------------------------------
A stellar blend missing from `linelist_solar.csv` is IN THE SUN: it appears in every
holding at the same wavelength, corrected or not. A telluric residual appears only in the
holdings whose correction did not reach it. So the two are separated by asking the SAME
question of a holding VERIFIED corrected in that band -- if it is clean and the suspect
is not, the absorption is atmospheric, not solar.
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
#: Displaced-template null, as RYA-1192 settled: same data, same window, registration
#: destroyed. 120 km/s is 2.9 A at 7250 A -- far outside any telluric line width.
NULL_SHIFTS_KMS = tuple(v for v in range(-400, 401, 10) if abs(v) >= 120)
NULL_Z = 3.0

TEMPLATE = ("iag_fts_solar_atlas", "solar_iag_reiners2016", "solar_iag")

#: Windows with NO registered telluric complex, used to measure how much structure the
#: template carries that is not telluric. A band whose template is no deeper than these
#: cannot be judged by it.
TEMPLATE_CONTROL_WINDOWS = ((5500.0, 5600.0), (6000.0, 6100.0), (6400.0, 6500.0))
#: A band's template must be this many times deeper than the worst clean window before the
#: shape test is allowed to speak. 3x is deliberately loose: the deep bands clear it by
#: 9-17x and o2gamma fails it, so nothing hinges on the exact value.
TEMPLATE_OVER_CONTROL = 3.0

HOLDINGS = [("kpno_solar_atlas", "solar_kpno", "UNCORRECTED REFERENCE"),
            ("iag_fts_solar_atlas", "solar_iag_reiners2016", "UNCORRECTED REFERENCE"),
            ("kpno_solar_atlas", "solar_kpno_molecfit_corrected", "claims molecfit"),
            ("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected", "claims delivered"),
            ("harps", "solar_harps", "UNCORRECTED REFERENCE"),
            ("harps", "solar_harps_molecfit_corrected", "claims molecfit"),
            ("iag_fts_solar_atlas", "solar_iag", "claims delivered")]


def _prof(V, inst, hold, lo, hi):
    w, f = V._load_raw(inst, hold, 0.5 * (lo + hi), 0.5 * (hi - lo))
    o = np.argsort(w)
    return w[o], f[o]


def _flat(w, f, frac=0.25, q=0.95):
    n = len(w)
    win = max(31, int(n * frac) | 1)
    c = (pd.Series(f).rolling(win, center=True, min_periods=win // 4)
         .quantile(q).bfill().ffill().to_numpy())
    return f / np.where(c > 0, c, np.nan)


def _corr(a, b):
    a = a - np.nanmean(a)
    b = b - np.nanmean(b)
    d = np.sqrt(np.nansum(a * a) * np.nansum(b * b))
    return float(np.nansum(a * b) / d) if d > 0 else float("nan")


def template(V, lo, hi):
    """Telluric transmission the IAG correction removed, on its own wavelength grid."""
    inst, raw, corr_ = TEMPLATE
    wr, fr = _prof(V, inst, raw, lo - 2, hi + 2)
    wc, fc = _prof(V, inst, corr_, lo - 2, hi + 2)
    tel = np.clip(np.interp(wc, wr, fr) / np.clip(fc, 1e-6, None), 0.0, 2.0)
    return wc, tel


_ctrl_cache: dict = {}


def _control_level(V) -> dict:
    """How much structure the template carries where there is NO telluric.

    🔴 THE NULL LEG B DID NOT HAVE. The displaced null asks whether a correlation survives
    sliding the template off the lines; it cannot notice that the template itself is not
    telluric. `solar_iag` and `solar_iag_reiners2016` are SEPARATELY REDUCED, so their
    ratio carries reduction differences everywhere — up to 0.476 minimum transmission and
    5.7% of pixels below 0.95 in the clean 6400-6500 A window. Any band whose template is
    not clearly deeper than that is being judged against an artifact.
    """
    if "c" not in _ctrl_cache:
        fracs, mins = [], []
        for lo, hi in TEMPLATE_CONTROL_WINDOWS:
            try:
                _, tel = template(V, lo, hi)
            except Exception:
                continue
            fracs.append(float(np.mean(tel < 0.95)))
            mins.append(float(np.nanmin(tel)))
        _ctrl_cache["c"] = {
            "windows": [list(w) for w in TEMPLATE_CONTROL_WINDOWS],
            "worst_frac": (max(fracs) if fracs else 0.0),
            "worst_min": (min(mins) if mins else 1.0),
            "note": ("the template is NOT flat where there is no telluric — this is the "
                     "level a band must clear before its shape test means anything"),
        }
    return _ctrl_cache["c"]


def measure(V, inst, hold, lo, hi, tw, tel):
    """`residual_r` > 0 means the holding STILL carries the template's telluric lines."""
    try:
        w, f = _prof(V, inst, hold, lo - 2, hi + 2)
    except Exception as e:
        return {"state": "OUT-OF-SPAN-OR-REFUSED", "why": f"{type(e).__name__}: {str(e)[:70]}"}
    if len(f) < 200:
        return {"state": "TOO-THIN", "n": int(len(f))}
    g = np.linspace(lo, hi, 6000)
    a = 1.0 - np.interp(g, w, _flat(w, f))
    def r_at(dv):
        return -_corr(a, np.interp(g, tw * (1.0 + dv / C_KMS), tel) - 1.0)
    r0 = r_at(0.0)
    null = np.array([r_at(v) for v in NULL_SHIFTS_KMS], float)
    null = null[np.isfinite(null)]
    mu, sd = float(np.mean(null)), float(np.std(null, ddof=1))
    z = (r0 - mu) / sd if sd > 0 else float("nan")
    return {"state": "MEASURED", "residual_r": round(r0, 4), "residual_z": round(z, 2),
            "null_mean": round(mu, 4), "null_sd": round(sd, 4), "n": int(len(f))}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="data/results/rya1191")
    a = ap.parse_args(argv)
    out = ROOT / a.out_dir
    out.mkdir(parents=True, exist_ok=True)

    import rya1192_sirius_verification as V
    from pipeline.telluric_policy import TELLURIC_BANDS
    ls = pd.read_csv(ROOT / "data/linelists/linelist_solar.csv", low_memory=False)

    bands = []
    for lo, hi, name in TELLURIC_BANDS:
        lo, hi = float(lo), float(hi)
        try:
            tw, tel = template(V, lo, hi)
        except Exception as e:
            bands.append({"band": str(name), "lo_A": lo, "hi_A": hi,
                          "state": "NO-TEMPLATE",
                          "why": f"the IAG pair does not reach this band: {str(e)[:70]}"})
            continue
        depth_frac = float(np.mean(tel < 0.95))
        ctrl = _control_level(V)
        row_ok_by_template = depth_frac >= TEMPLATE_OVER_CONTROL * ctrl["worst_frac"]
        row = {"band": str(name), "lo_A": lo, "hi_A": hi, "state": "MEASURED",
               "template_min": round(float(np.nanmin(tel)), 4),
               "template_frac_below_0.95": round(depth_frac, 4),
               "template_frac_below_0.95_in_clean_controls": ctrl["worst_frac"],
               "template_clears_the_clean_window_null": bool(row_ok_by_template),
               "catalogued_stellar_depth_per_A": round(float(
                   ls.loc[ls.wavelength_air_A.between(lo, hi), "central_depth"].sum()
                   / (hi - lo)), 3),
               "holdings": {}}
        for inst, hold, note in HOLDINGS:
            m = measure(V, inst, hold, lo, hi, tw, tel)
            m["note"] = note
            # the DEPTH statistic beside it, because the two are blind in opposite places
            d_in, d_side = V._depth(inst, hold, lo, hi), V._side_depth(inst, hold, lo, hi)
            m["pct_below_0.8"] = (round(d_in, 3) if d_in is not None else None)
            m["side_pct_below_0.8"] = (round(d_side, 3) if d_side is not None else None)
            m["depth_over_side"] = (round(d_in / d_side, 2)
                                    if d_in is not None and d_side else None)
            row["holdings"][hold] = m
        row.update(_verdict(row))
        bands.append(row)

    doc = {
        "ticket": "RYA-1191 (B) — telluric residual vs uncatalogued blend",
        "kind": "VERIFICATION — no correction applied, no value moved",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "template": {
            "raw": TEMPLATE[1], "corrected": TEMPLATE[2],
            "what": ("their ratio IS the telluric transmission the Baker+2020 correction "
                     "removed — one atlas, one instrument, so the ratio carries no "
                     "instrumental difference"),
            "caveat": ("⚠️ LINE POSITIONS ONLY. Airmass and PWV differ between nights, so "
                       "scoring is by shape correlation and never by depth."),
        },
        "clean_window_null": _control_level(V),
        "null_estimator": (f"per band, r(0) against r(dv) for {len(NULL_SHIFTS_KMS)} "
                           f"shifts with |dv| >= 120 km/s; residual at z >= {NULL_Z}"),
        "sign_convention": ("residual_r is FLIPPED so that LARGER POSITIVE = MORE "
                            "RESIDUAL TELLURIC. The uncorrected references are in every "
                            "row as the positive control."),
        "bands": bands,
        "supersedes": _supersedes(bands),
    }
    (out / "rya1191_telluric_residual.json").write_text(json.dumps(doc, indent=2) + "\n")

    print(f"{'band':<16}{'span':<14}{'tmpl':>6}  " +
          "".join(f"{h[1].replace('solar_','')[:13]:>15}" for h in HOLDINGS))
    for b in bands:
        if b["state"] != "MEASURED":
            print(f"{b['band'][:15]:<16}{b['lo_A']:.0f}-{b['hi_A']:<8.0f} {b['state']}")
            continue
        cells = []
        for _, hold, _ in HOLDINGS:
            c = b["holdings"][hold]
            cells.append(f"{c['residual_z']:>+8.1f}s" if c["state"] == "MEASURED"
                         else f"{'--':>9}")
        print(f"{b['band'][:15]:<16}{b['lo_A']:.0f}-{b['hi_A']:<8.0f}"
              f"{b['template_frac_below_0.95']:>6.3f}  " +
              "".join(f"{c:>15}" for c in cells))
    print("\n=== verdicts ===")
    for b in bands:
        if b.get("verdict"):
            print(f"  {b['band']} {b['lo_A']:.0f}-{b['hi_A']:.0f}: {b['verdict']}")
            for g in b.get("gaps", []):
                print(f"      🔴 {g}")
    print(f"\nwrote {a.out_dir}/")
    return 0


def _verdict(row) -> dict:
    """Which holdings still carry this band's telluric, judged against the raw references.

    ⚠️ A band can only judge anything where the UNCORRECTED references themselves light
    up. Where they do not, the template has no lines here and the row says so rather than
    reporting every holding as clean.
    """
    # 🔴 FIRST: is the TEMPLATE telluric here at all? A band whose template is no deeper
    # than the clean-window controls is being judged against a reduction difference, and
    # no amount of displaced-null significance rescues that (o2gamma: template 13.1% of
    # pixels below 0.95 against a 5.7% CONTROL, while molecfit's physical model of the
    # same band absorbs at most 1.78%).
    if not row.get("template_clears_the_clean_window_null", True):
        return {"verdict": (
            f"NO TEMPLATE — this band's template is {row['template_frac_below_0.95']:.4f} "
            f"of pixels below 0.95 against {row['template_frac_below_0.95_in_clean_controls']:.4f} "
            f"in windows with NO telluric, so it is not clearly deeper than the "
            f"reduction difference between the two IAG products. ⚠️ NOT a clean verdict: "
            f"nothing here can be judged by shape (RYA-833)."),
            "judged_by": None, "gaps": []}
    refs = [c for h, c in row["holdings"].items()
            if c.get("note") == "UNCORRECTED REFERENCE" and c["state"] == "MEASURED"]
    powered = [c for c in refs if c["residual_z"] >= 3 * NULL_Z]
    if not powered:
        # 🔴 THE TEMPLATE TEST LOSES POWER EXACTLY WHERE A BAND SATURATES, and that is not
        # a failure to work around -- it is the reason the DEPTH test exists. In the O2 A
        # and B cores the raw flux is near zero over broad stretches, so "absorption" is
        # ~1 everywhere in both spectra and there is almost no shape left to correlate.
        # Depth is blind to a SHALLOW band and sharp on a saturated one; the template is
        # the reverse. Neither is the better test, and which one answers is decided here
        # by whether the uncorrected references light it up, never by preference.
        ref_depth = [c["depth_over_side"] for c in refs
                     if c.get("depth_over_side") is not None]
        if ref_depth and max(ref_depth) >= 2.0:
            gaps, clean = [], []
            for hold, c in row["holdings"].items():
                if c.get("note") == "UNCORRECTED REFERENCE" or c["state"] != "MEASURED":
                    continue
                r = c.get("depth_over_side")
                if r is None:
                    continue
                (gaps if r >= 2.0 else clean).append(
                    f"{hold} ({r:.2f}x its own side windows)")
            return {"verdict": (
                f"TEMPLATE TEST HAS NO POWER (saturated band — the references reach only "
                f"{max(x['residual_z'] for x in refs):+.1f} sigma on shape), so this band "
                f"is judged on DEPTH, where they reach {max(ref_depth):.1f}x their own "
                f"side windows. {len(clean)} corrected holding(s) clean, {len(gaps)} not."),
                "judged_by": "DEPTH (template saturated)",
                "gaps": [f"{g} — retains saturated absorption" for g in gaps]}
        return {"verdict": ("NO POWER BY EITHER TEST — no uncorrected reference stands "
                            "above its own displaced null on shape, and none is deeper "
                            "than 2x its own side windows. Nothing here can be judged, "
                            "which is not the same as nothing being here (RYA-833)"),
                "judged_by": None, "gaps": []}
    gaps = []
    for hold, c in row["holdings"].items():
        if c.get("note") == "UNCORRECTED REFERENCE" or c["state"] != "MEASURED":
            continue
        if c["residual_z"] >= NULL_Z:
            worst = max(x["residual_z"] for x in powered)
            gaps.append(f"{hold}: residual at {c['residual_z']:+.1f} sigma "
                        f"(uncorrected references reach {worst:+.1f}) — "
                        + ("PARTIAL correction" if c["residual_z"] < 0.5 * worst
                           else "essentially UNCORRECTED"))
    clean = [h for h, c in row["holdings"].items()
             if c.get("note") != "UNCORRECTED REFERENCE" and c["state"] == "MEASURED"
             and c["residual_z"] < NULL_Z]
    v = (f"references light up at {max(x['residual_z'] for x in powered):+.1f} sigma; "
         f"{len(clean)} corrected holding(s) clean, {len(gaps)} carrying a residual")
    row["judged_by"] = "TEMPLATE (line positions)"
    if gaps:
        v += (". ⚠️ NOT A STELLAR BLEND: a line missing from linelist_solar would be in "
              "the Sun and would show in EVERY holding, and the clean ones are clean.")
    return {"verdict": v, "judged_by": "TEMPLATE (line positions)", "gaps": gaps}


def _supersedes(bands) -> dict:
    """What this leg re-judges in RYA-1192's band verdicts, and why that test failed."""
    return {
        "artifact": "data/results/rya1192/rya1192_sirius_verification.json",
        "findings_rejudged": [
            ("solar_kpno_kurucz2005_corrected | H2O 7160-7340 — RYA-1192 returned "
             "VERIFIED-RAW-LIKE from a DEPTH ratio of 4.2 against a 5x bar. Wrong "
             "comparand: against the telluric TEMPLATE it is clean, and its own side "
             "windows put it level with solar_iag, which is verified corrected there."),
            ("solar_kpno_molecfit_corrected | O2 gamma 6270-6300 — RYA-1192 returned "
             "VERIFIED-CORRECTED from a mean |diff| of 0.030 against a zero floor. That "
             "difference is a FEATURELESS ~3.6% rescale (corrected/raw min 1.0000, zero "
             "pixels below 0.98); the O2 gamma lines are still there."),
        ],
        "why_the_earlier_tests_missed_it": (
            "A depth cut cannot see a shallow band and a difference cannot tell a "
            "rescale from a line removal. Both are answered by asking whether the flux "
            "absorbs AT THE TELLURIC LINE POSITIONS."),
    }


if __name__ == "__main__":
    raise SystemExit(main())
