#!/usr/bin/env python3
"""RYA-1192 — is the telluric correction ACTUALLY applied? Verify, do not trust the label.

    python3 scripts/rya1192_telluric_verification.py

VERIFICATION ONLY. Applies no correction, writes nothing outside `data/results/rya1192/`,
moves no value. The fixes are separate tickets (RYA-1191).

WHY THE LABEL CANNOT BE TRUSTED
--------------------------------
RYA-1190 found `solar_kpno_molecfit_corrected` byte-identical to the raw KP atlas over 7
of 10 sampled red-optical windows. The mechanism is in the label itself: molecfit ran on
a FIXED SET OF BANDS, and any window outside them is raw under a "corrected" name. This
generalises that check to every live Fe product.

TWO CHECKS, AND THEY ANSWER DIFFERENT QUESTIONS
------------------------------------------------
(a) DIRECT COMPARISON against the holding's own raw sibling. `max|corrected - raw|` over a
    window: exactly zero means the correction never ran there. This is decisive and needs
    no model -- but it needs the raw sibling to be READABLE, which is where most of the
    scope is lost on this machine (see below).

(b) TELLURIC-FOREST SIGNATURE, for holdings with no readable raw. Fraction of pixels below
    0.80 inside a telluric band. ⚠️ A bare threshold would be a chosen number, so it is
    CALIBRATED on a pair where the answer is already known: KP1984 raw vs KP1984 molecfit,
    same atlas, same window, one corrected and one not.

        O2 A-band 7594-7685   raw 0.555   corrected 0.048
        O2 B-band 6867-6884   raw 0.291   corrected 0.043
        control  5000-5040    raw 0.203   corrected 0.203  (byte-identical, as it must be)

    The control is the part that makes it a calibration rather than an assertion: outside
    the fitted bands the two are the SAME SPECTRUM, so the statistic is measuring the
    correction and not a difference between the products.

🔴 THE GROUND TRUTH FOR "WHICH BANDS" IS THE FIT MANIFESTS, NOT THE PROSE, AND THEY
DISAGREE. The holding's `telluric_state` says "six registered bands (RYA-940)" and
`telluric_policy.TELLURIC_BANDS` lists six. `data/audit/rya940_kp1984_telluric/` holds
SEVEN `fit_manifest.json` files -- there is an `o2gamma` band at 6270-6300 A that molecfit
fitted and the policy does not list. Both are read here and the disagreement is reported
rather than resolved by preferring one.

⚠️ WHAT THIS MACHINE CANNOT SEE. `solar_kpno_kurucz2005_corrected`, `solar_iag`,
`solar_iag_reiners2016` and `solar_vesta_crires_plus_idp` are staged on Sirius/srv and are
UNREADABLE here. They are reported UNVERIFIABLE-HERE, which is NOT the same as unverified
in principle and must not be recorded as either "corrected" or "raw" (RYA-833). Losing the
CRIRES+ raw sibling is why that arm gets check (b) rather than check (a).
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

FEED = ROOT / "data/products/solar/Fe.json"
BP = ROOT / "data/results/band_products"
FITS = ROOT / "data/audit/rya940_kp1984_telluric"

#: corrected holding -> (instrument, its RAW sibling or None)
RAW_SIBLING = {
    "solar_kpno_molecfit_corrected":   ("kpno_solar_atlas", "solar_kpno"),
    "solar_kpno_kurucz2005_corrected": ("kpno_solar_atlas", None),
    "solar_harps_molecfit_corrected":  ("harps", "solar_harps"),
    "solar_iag":                       ("iag_fts_solar_atlas", None),
    "solar_crires_plus_y_wide_rya1054": ("crires_plus", "solar_vesta_crires_plus_idp"),
}
DEEP = 0.80          # "deep absorption" cut, calibrated below on the KP pair
IDENTICAL = 1e-9     # float equality

#: 🔴 A NON-ZERO DIFFERENCE IS NOT AUTOMATICALLY A CORRECTION, AND TREATING IT AS ONE
#: OVERSTATED HARPS BY 184 LINES. KP's corrected product is the SAME ARRAY as its raw with
#: molecfit applied in-band, so an exact zero there is decisive. HARPS is different: the
#: corrected and raw products are SEPARATELY REDUCED, and they differ by ~1e-4 EVERYWHERE
#: as a reduction baseline. Measured on HARPS:
#:
#:     O2 B-band 6867-6884 (real telluric)   mean|diff| 0.264432   max 1.007
#:     control   5000-5100                   mean|diff| 0.000069
#:     control   4500-4600                   mean|diff| 0.000114
#:
#: So the correction is FOUR ORDERS larger where telluric actually bites, and a 1e-4
#: difference at a measured line is the reduction, not a telluric correction. Every
#: holding therefore gets its own baseline from control windows, and a line is called
#: corrected only if it clears that baseline by MATERIAL_MULTIPLE.
CONTROL_WINDOWS = ((5000.0, 5100.0), (4500.0, 4600.0))
MATERIAL_MULTIPLE = 10.0


def fitted_bands() -> list[dict]:
    """The bands molecfit ACTUALLY fitted, read from its own manifests."""
    out = []
    for f in sorted(FITS.glob("*/fit_manifest.json")):
        d = json.loads(f.read_text())
        out.append({"band": f.parent.name, "band_A": d["band_A"],
                    "fit_window_A": d["fit_window_A"],
                    "molecules": d.get("molecules"),
                    "n_segments": len(d.get("segments", []))})
    return out


def policy_bands() -> list[list]:
    from pipeline.telluric_policy import TELLURIC_BANDS
    return [[float(a), float(b), c] for a, b, c in TELLURIC_BANDS]


def _load(inst, hold, centre, pad):
    from measure_band_ew import load_window_ex
    w = load_window_ex(inst, centre, pad, holding=hold, allow_uncorrected=True)
    a = np.asarray(w.wave, float); f = np.asarray(w.flux, float)
    k = np.isfinite(a) & np.isfinite(f)
    return a[k], f[k]


def baseline(inst, corrected, raw) -> float | None:
    """The holding's own corrected-vs-raw difference OUTSIDE any telluric band.

    For a product pair that is one array plus an in-band correction (KP) this is exactly
    0. For a pair that was separately reduced (HARPS) it is the reduction floor, and it is
    what a per-line difference has to clear before it can be called a correction.
    """
    vals = []
    for lo, hi in CONTROL_WINDOWS:
        try:
            wc, fc = _load(inst, corrected, 0.5 * (lo + hi), 0.5 * (hi - lo))
            wr, fr = _load(inst, raw, 0.5 * (lo + hi), 0.5 * (hi - lo))
        except Exception:
            continue
        n = min(len(fc), len(fr))
        if n > 50:
            vals.append(float(np.nanmean(np.abs(fc[:n] - fr[:n]))))
    return float(np.mean(vals)) if vals else None


def compare(inst, corrected, raw, centre, pad, base=None) -> dict:
    """max|corrected - raw| over one window, judged against the holding's own baseline."""
    try:
        wc, fc = _load(inst, corrected, centre, pad)
    except Exception as e:
        return {"state": "UNVERIFIABLE-HERE", "why": f"corrected unreadable: {str(e)[:70]}"}
    if raw is None:
        return {"state": "NO-RAW-SIBLING", "n": int(len(fc))}
    try:
        wr, fr = _load(inst, raw, centre, pad)
    except Exception as e:
        return {"state": "UNVERIFIABLE-HERE", "why": f"raw unreadable: {str(e)[:70]}",
                "n": int(len(fc))}
    n = min(len(fc), len(fr))
    if n < 5:
        return {"state": "UNVERIFIABLE-HERE", "why": f"only {n} comparable points"}
    dmax = float(np.nanmax(np.abs(fc[:n] - fr[:n])))
    if dmax <= IDENTICAL:
        state = "VERIFIED-RAW"
    elif base is not None and base > 0 and dmax <= MATERIAL_MULTIPLE * base:
        # a difference no bigger than this pair's own reduction floor
        state = "NO-MATERIAL-DIFFERENCE"
    else:
        state = "VERIFIED-CORRECTED"
    return {"state": state, "max_abs_diff": round(dmax, 8),
            "baseline_abs_diff": (round(base, 8) if base is not None else None),
            "n": int(n)}


_CAT = {}


def _catalogue_predicted(w, centre, R):
    """Deep-absorption the STELLAR line list alone predicts here — the comparand for a
    window with no readable raw sibling. Same accounting as RYA-1190: exp(-SUM tau_i), so
    overlap and saturation behave. NOT radiative transfer, and it over-predicts, which is
    why only an EXCESS of observed over predicted is treated as a signal."""
    if not _CAT:
        d = pd.read_csv(ROOT / "data/linelists/linelist_solar.csv", low_memory=False)[
            ["wavelength_air_A", "central_depth"]].dropna()
        _CAT["w"] = d.wavelength_air_A.to_numpy()
        _CAT["t"] = -np.log(1.0 - np.clip(d.central_depth.to_numpy(), 0.0, 0.999))
    sig = float(np.hypot(centre * 1.6 / 2.998e5, (centre / R) / 2.355))
    m = (_CAT["w"] > w[0] - 2.0) & (_CAT["w"] < w[-1] + 2.0)
    tau = np.zeros_like(w)
    for a, t in zip(_CAT["w"][m], _CAT["t"][m]):
        z = (w - a) / sig
        k = np.abs(z) < 6
        if k.any():
            tau[k] += t * np.exp(-0.5 * z[k] ** 2)
    return np.exp(-tau)


def forest(inst, hold, lo, hi, R=None) -> dict:
    """Fraction of pixels below DEEP inside a window — check (b).

    With `R` (a resolving power) it also computes what the STELLAR catalogue predicts, so
    an EXCESS of observed over predicted can be separated from a band that is simply full
    of stellar lines. That excess is the only telluric-shaped signal available where no
    raw sibling can be read.
    """
    try:
        w, f = _load(inst, hold, 0.5 * (lo + hi), 0.5 * (hi - lo))
    except Exception as e:
        return {"state": "UNVERIFIABLE-HERE", "why": str(e)[:80]}
    if len(f) < 20:
        return {"state": "UNVERIFIABLE-HERE", "why": f"only {len(f)} points"}
    out = {"n": int(len(f)), "frac_below_0.80": round(float((f < DEEP).mean()), 4),
           "frac_below_0.95": round(float((f < 0.95).mean()), 4),
           "mean_flux": round(float(f.mean()), 5)}
    if R:
        p = _catalogue_predicted(w, 0.5 * (lo + hi), R)
        out["catalogue_frac_below_0.80"] = round(float((p < DEEP).mean()), 4)
        out["excess_over_catalogue"] = round(out["frac_below_0.80"]
                                             - out["catalogue_frac_below_0.80"], 4)
    return out


def measured_lines() -> pd.DataFrame:
    """Every live Fe product's measured lines, from the per-line artifacts on disk.

    ⚠️ The feed's `provenance.copied_to` resolves for only 18 of 70 live products, so the
    stems are globbed out of band_products instead and the shortfall is REPORTED. No
    per-line artifact exists on disk for ANY red-optical, NIR or CRIRES+ Fe product, which
    is why those are verified at WINDOW level and said to be.
    """
    feed = json.loads(FEED.read_text())
    want = {(p["holding"], p["band"], p["instrument"]) for p in feed["products"]}
    rows = []
    for f in sorted(BP.glob("Fe*_lines.csv")):
        d = pd.read_csv(f)
        if "wavelength_air_A" not in d.columns:
            continue
        for (hold, band, inst) in want:
            if hold in f.name:
                for w in d.wavelength_air_A.dropna():
                    rows.append((hold, inst, round(float(w), 4)))
                break
    return pd.DataFrame(rows, columns=["holding", "instrument",
                                       "wavelength_air_A"]).drop_duplicates()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pad", type=float, default=0.30,
                    help="half-window around each measured line (A)")
    ap.add_argument("--out-dir", default="data/results/rya1192")
    a = ap.parse_args(argv)
    out_dir = ROOT / a.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    fb, pb = fitted_bands(), policy_bands()
    fitted_names = {b["band"] for b in fb}
    policy_spans = {(round(x[0], 1), round(x[1], 1)) for x in pb}
    fitted_spans = {(round(b["band_A"][0], 1), round(b["band_A"][1], 1)) for b in fb}

    # ── per measured line ──────────────────────────────────────────────────────────
    lines = measured_lines()
    bases = {}
    for hold, (inst, raw) in RAW_SIBLING.items():
        bases[hold] = baseline(inst, hold, raw) if raw else None
    rows = []
    for _, r in lines.iterrows():
        inst, raw = RAW_SIBLING.get(r.holding, (r.instrument, None))
        res = compare(inst, r.holding, raw, float(r.wavelength_air_A), a.pad,
                      base=bases.get(r.holding))
        # 🔴 `band_A` IS THE CORRECTED SPAN; `fit_window_A` IS 25 A WIDER EACH SIDE AND IS
        # ONLY THE DATA THE FIT WAS CONSTRAINED ON. Using the fit window as the membership
        # test called 4 KP lines "inside a fitted band" while the comparison found them
        # byte-identical to raw -- and the measurement is right, the test was wrong: all
        # four (6315.811, 6843.655, 6855.161, 6858.148) sit in the PAD, never in a band
        # core. Recorded as both, because "the pad is not corrected" is itself a fact a
        # reader of the manifests would otherwise assume the other way.
        in_core = any(b["band_A"][0] <= r.wavelength_air_A <= b["band_A"][1] for b in fb)
        in_pad = any(b["fit_window_A"][0] <= r.wavelength_air_A <= b["fit_window_A"][1]
                     for b in fb)
        rows.append(dict(holding=r.holding, instrument=inst,
                         wavelength_air_A=r.wavelength_air_A,
                         inside_corrected_band=in_core,
                         inside_fit_window_only=bool(in_pad and not in_core), **res))
    per_line = pd.DataFrame(rows)
    per_line.to_csv(out_dir / "rya1192_per_line.csv", index=False)

    # ── the calibration + the forest statistic ────────────────────────────────────
    cal_bands = [(b["band_A"][0], b["band_A"][1], b["band"]) for b in fb]
    cal_bands.append((5000.0, 5040.0, "CONTROL (no telluric band)"))
    calibration = []
    for lo, hi, name in cal_bands:
        calibration.append({
            "band": name, "lo_A": lo, "hi_A": hi,
            "kp_raw": forest("kpno_solar_atlas", "solar_kpno", lo, hi),
            "kp_molecfit": forest("kpno_solar_atlas", "solar_kpno_molecfit_corrected", lo, hi),
        })

    # ── CRIRES+, the priority: EVERY arm across its FULL span ──────────────────────
    #
    # 🔴 THE FIRST CUT OF THIS AUDIT STOPPED AT 10796 A AND THAT WAS WRONG. It scanned one
    # of the three CRIRES+ holdings -- the Y-wide arm -- and called that "the full measured
    # span". CRIRES+ reaches 53000 A per instrument_catalog.csv, the H arm
    # (solar_crires_plus_h_rya1094) was already loading in RYA-1189/1190, and canonical_gf
    # carries 29 GRADED Fe lines beyond 10796 of which 27 are Ruffoni-2013 in the H window.
    # Auditing "is the correction applied" over a window that excludes every graded line
    # past 1 micron answers a question nobody asked.
    CRIRES_ARMS = (("solar_crires_plus_y_wide_rya1054", 9800, 10796),
                   ("solar_crires_plus_h_rya1094", 15000, 17500))
    crires = []
    for hold, a0, a1 in CRIRES_ARMS:
        for lo in range(a0, a1, 100):
            hi = min(lo + 100, a1)
            row = {"holding": hold, "lo_A": lo, "hi_A": hi,
                   "crires": forest("crires_plus", hold, lo, hi, R=70000.0)}
            # KP1984 stops at ~13000 A, so it is a comparand for the Y arm only. Where it
            # cannot reach, say so rather than omitting the column.
            row["kp_raw_same_window"] = (forest("kpno_solar_atlas", "solar_kpno", lo, hi)
                                         if hi <= 13000 else
                                         {"state": "NO-COMPARAND", "why": "KP1984 ends ~13000 A"})
            crires.append(row)

    # ── beyond 10796: the graded lines the first cut never reached ─────────────────
    cgf = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    lab = cgf[(cgf.species.astype(str).str.startswith("Fe"))
              & (cgf.wavelength_air_A > 10796)
              & (cgf.gf_tier.astype(str) == "LAB")]
    from pipeline.telluric_policy import in_telluric_band
    beyond = []
    for _, r in lab.sort_values("wavelength_air_A").iterrows():
        w = float(r.wavelength_air_A)
        servedby, npts = None, 0
        for hold, _, _ in CRIRES_ARMS:
            try:
                win = _load("crires_plus", hold, w, 0.5)
                if len(win[1]) >= 20:
                    servedby, npts = hold, len(win[1]); break
            except Exception:
                pass
        rec = {"wavelength_air_A": round(w, 4), "species": str(r.species),
               "lab_source_tag": str(r.lab_source_tag),
               "policy_says_telluric_band": bool(in_telluric_band(w)),
               "served_by": servedby, "n_points": npts}
        if servedby is None:
            # not on CRIRES+ — is it on KP, and is KP corrected there?
            rec.update(compare("kpno_solar_atlas", "solar_kpno_molecfit_corrected",
                               "solar_kpno", w, 0.5,
                               base=bases.get("solar_kpno_molecfit_corrected")))
        beyond.append(rec)

    # the H arm's ACTUAL coverage, probed rather than read off the registry
    h_lo = h_hi = None
    for c in range(14900, 17600, 25):
        try:
            _load("crires_plus", "solar_crires_plus_h_rya1094", float(c), 1.0)
            h_lo = c if h_lo is None else h_lo
            h_hi = c
        except Exception:
            pass

    # ── coverage: does the corrected span reach each product's band at all? ────────
    #
    # 🔴 SPEC 1(b), AND IT IS THE ONLY CHECK AVAILABLE FOR THE BANDS WITH NO PER-LINE
    # ARTIFACT. No _lines.csv exists on disk for ANY red-optical, NIR or CRIRES+ Fe
    # product, so those cannot be verified line by line here. What CAN be measured is the
    # fraction of each product's declared band that molecfit corrected at all -- which for
    # a band-limited correction is the number that decides whether ANY of its lines could
    # have been touched.
    feed = json.loads(FEED.read_text())
    cov = []
    seen = set()
    for pr in feed["products"]:
        rng = pr.get("wavelength_range_A")
        if not rng:
            continue
        key = (pr["holding"], pr["band"], tuple(rng))
        if key in seen:
            continue
        seen.add(key)
        lo, hi = float(rng[0]), float(rng[1])
        span = hi - lo
        if pr["holding"] in ("solar_kpno_molecfit_corrected",):
            covered = 0.0
            for b in fb:
                a0, a1 = max(lo, b["band_A"][0]), min(hi, b["band_A"][1])
                covered += max(0.0, a1 - a0)
            frac = covered / span if span else 0.0
            basis = "molecfit band_A spans (RYA-940 manifests)"
        elif pr["holding"] == "solar_harps_molecfit_corrected":
            frac, basis = None, "global per-exposure correction — not band-limited"
        else:
            frac, basis = None, "no band registry available for this holding"
        cov.append({"holding": pr["holding"], "band": pr["band"],
                    "window_A": [lo, hi], "span_A": round(span, 1),
                    "fraction_of_band_corrected": (round(frac, 4) if frac is not None else None),
                    "basis": basis,
                    "n_lines_in_product": pr.get("n_lines")})

    # ── availability, stated ───────────────────────────────────────────────────────
    availability = {}
    for hold, (inst, raw) in RAW_SIBLING.items():
        c = compare(inst, hold, raw, 5500.0 if inst != "crires_plus" else 10300.0, 1.0)
        availability[hold] = {"instrument": inst, "raw_sibling": raw,
                              "probe": c.get("state"), "why": c.get("why")}

    doc = {
        "ticket": "RYA-1192",
        "kind": "VERIFICATION ONLY — no correction applied, no value moved",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "bands_molecfit_ACTUALLY_fitted": fb,
        "bands_the_policy_LISTS": pb,
        "band_registry_disagreement": {
            "fitted_not_in_policy": sorted(fitted_spans - policy_spans),
            "policy_not_fitted": sorted(policy_spans - fitted_spans),
            "note": ("The holding's telluric_state says SIX registered bands and "
                     "telluric_policy.TELLURIC_BANDS lists six; RYA-940 shipped SEVEN "
                     "fit manifests. Reported, not resolved."),
        },
        "forest_calibration": calibration,
        "crires_full_span": crires,
        "crires_h_arm_coverage": {
            "declared": [15007, 17494],
            "probed_actual": [h_lo, h_hi],
            "note": ("probed at 25 A steps. The declared span comes from "
                     "prenormalised_guard's comment; the arm serves a narrower range, and "
                     "an echelle also has INTER-ORDER GAPS inside it -- see how many "
                     "graded lines below come back unserved despite sitting in range."),
        },
        "beyond_10796_graded_lines": beyond,
        "telluric_policy_reach": {
            "TELLURIC_BANDS_max_A": 11560.0,
            "crires_plus_instrument_reach_A": [9500.0, 53000.0],
            "graded_Fe_lines_beyond_10796": int(len(lab)),
            "of_those_in_a_declared_telluric_band": int(sum(
                1 for b in beyond if b["policy_says_telluric_band"])),
            "finding": ("🔴 telluric_policy.TELLURIC_BANDS stops at 11560 A while CRIRES+ "
                        "reaches 53000 A, so in_telluric_band() answers FALSE for every "
                        "graded Fe line past 1 micron -- and FALSE reads as 'clean', not "
                        "as 'no band declared here'. The H window it covers is bracketed "
                        "by the 1.38 and 1.9 micron H2O bands with CO2 and CH4 inside it; "
                        "the registry simply does not describe that region."),
        },
        "band_coverage_by_product": cov,
        "crires_excess_windows": None,   # filled below
        "crires_verdict_legacy": (
            "NO EVIDENCE OF UNCORRECTED TELLURIC across the full 9800-10796 A measured "
            "span. Window by window the deep-absorption fraction tracks the KP1984 RAW "
            "atlas to within ~0.003 -- and KP raw is uncorrected, so the deep features "
            "both show are STELLAR. For contrast the same statistic separates raw from "
            "corrected by an order of magnitude where telluric actually bites (O2 A-band: "
            "raw 0.555 vs corrected 0.048). ⚠️ THIS IS CHECK (b), NOT PROOF THE "
            "CORRECTION RAN: solar_vesta_crires_plus_idp is staged on Sirius and "
            "unreadable here, so a direct corrected-vs-raw comparison is OWED there. In a "
            "window this clean the two would look alike either way -- what is established "
            "is that no telluric residual REMAINS, which is the science-relevant half."),
        "reduction_baselines": {k: (round(v, 8) if v is not None else None)
                                for k, v in bases.items()},
        "material_multiple": MATERIAL_MULTIPLE,
        "harps_materially_corrected_caveat": (
            "⚠️ The handful of HARPS lines that clear the reduction floor do so by 11-19x "
            "at ~0.001 in flux -- still THREE ORDERS below the real O2 B-band correction "
            "(mean |diff| 0.264). Two of them (6855.161, 6858.148) sit just blueward of "
            "the O2 B-band edge at 6867 A, which is physically sensible; the rest are at "
            "the 4203-4245 A blue edge where the two reductions diverge most. Read them "
            "as 'above this pair's noise floor', NOT as 'a telluric correction of "
            "consequence' -- nothing in the HARPS Fe pool sits in a band where molecfit "
            "removed anything material."),
        "availability_on_this_machine": availability,
        "crires_verdict_by_arm": None,   # filled below
        "pad_is_not_corrected": {
            "what": ("fit_window_A is 25 A wider each side than band_A. Every measured "
                     "line falling in the pad but not the core came back byte-identical "
                     "to raw, so the pad constrains the fit and is not itself corrected."),
            "lines_in_pad_only": sorted(
                float(x) for x in per_line[per_line.inside_fit_window_only]
                .wavelength_air_A.unique()),
        },
        "per_line_summary": (per_line.groupby(["holding", "state"]).size()
                             .rename("n").reset_index().to_dict("records")
                             if len(per_line) else []),
    }
    # ── verdict per ARM, computed from the excess ──────────────────────────────────
    EXCESS = 0.02
    by_arm = {}
    exc_all = []
    for hold, _, _ in CRIRES_ARMS:
        rows = [r for r in crires if r["holding"] == hold
                and "excess_over_catalogue" in r["crires"]]
        hot = [r for r in rows if r["crires"]["excess_over_catalogue"] > EXCESS]
        exc_all.extend({"holding": hold, "lo_A": r["lo_A"], "hi_A": r["hi_A"],
                        "excess": r["crires"]["excess_over_catalogue"]} for r in hot)
        if not rows:
            by_arm[hold] = "NO DATA"
        elif not hot:
            by_arm[hold] = (
                f"NO EVIDENCE OF UNCORRECTED TELLURIC over {rows[0]['lo_A']}-"
                f"{rows[-1]['hi_A']} A: no window's observed deep-absorption fraction "
                f"exceeds what the stellar line list alone predicts. ⚠️ Check (b) only — "
                f"the raw sibling is on Sirius, so this is 'no residual REMAINS', not "
                f"proof the correction ran.")
        else:
            by_arm[hold] = (
                f"🔴 LOCALISED CANDIDATE RESIDUAL in {len(hot)} of {len(rows)} windows: "
                + ", ".join(f"{r['lo_A']}-{r['hi_A']} (+{r['crires']['excess_over_catalogue']:.4f})"
                            for r in hot)
                + ". Observed deep absorption EXCEEDS the stellar catalogue there while "
                  "every other window sits below it, and those windows are where telluric "
                  "CO2 (15700-16100) and the 1.9 micron H2O wing are expected. ⚠️ NOT "
                  "PROOF: the accounting is not radiative transfer and over-predicts in "
                  "general, so only the EXCESS and its localisation carry weight. A raw "
                  "comparison on Sirius would settle it.")
    doc["crires_verdict_by_arm"] = by_arm
    doc["crires_excess_windows"] = exc_all
    # which graded lines fall in an excess window
    for b in doc["beyond_10796_graded_lines"]:
        b["in_a_candidate_residual_window"] = any(
            e["lo_A"] <= b["wavelength_air_A"] <= e["hi_A"] for e in exc_all)

    (out_dir / "rya1192_verification.json").write_text(json.dumps(doc, indent=2) + "\n")

    print("=== bands molecfit ACTUALLY fitted (from its own manifests) ===")
    for b in fb:
        print(f"  {b['band']:<10} {b['band_A'][0]:>8.0f}-{b['band_A'][1]:<8.0f} "
              f"fit window {b['fit_window_A'][0]:>8.0f}-{b['fit_window_A'][1]:<8.0f} "
              f"{b['molecules']}  {b['n_segments']} segments")
    print(f"\n  fitted but NOT in telluric_policy: {doc['band_registry_disagreement']['fitted_not_in_policy']}")
    print(f"  in policy but NOT fitted        : {doc['band_registry_disagreement']['policy_not_fitted']}")

    print("\n=== per measured line ===")
    if len(per_line):
        print(per_line.groupby(["holding", "state"]).size().to_string())
    print("\n=== availability on this machine ===")
    for h, v in availability.items():
        print(f"  {h:<36} raw={str(v['raw_sibling']):<30} {v['probe']}")
    print(f"\nwrote {a.out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
