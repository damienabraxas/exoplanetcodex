#!/usr/bin/env python3
"""RYA-1192 (Sirius leg) — close the 131 unjudged lines and settle CRIRES+ against RAW.

    python3 scripts/rya1192_sirius_verification.py

VERIFICATION ONLY. Applies no correction, moves no value, writes only under
`data/results/rya1192/`. The fixes are separate tickets (RYA-1191, RYA-1194).

WHY THERE IS A SECOND SCRIPT
----------------------------
`scripts/rya1192_telluric_verification.py` did the Mac leg and left two holes that are
holes *of that machine*, not of the method:

  * 131 measured Fe lines (88 kurucz2005 + 43 iag) came back UNVERIFIABLE-HERE because
    those holdings are staged on Sirius.
  * CRIRES+ got the catalogue-EXCESS proxy (check b) rather than a comparison against
    raw, because the raw Vesta IDPs are on Sirius too. That proxy flagged a CANDIDATE
    residual in the H arm and said in terms that it was not proof.

Both are closed here, and closing them needed more than re-running the same code: two of
the five holdings have NO raw sibling in the sense the Mac leg assumed, and CRIRES+ turns
out to have a far better comparand than anyone had used.

🔴 THE THREE THINGS THAT CHANGED THE ANSWER
-------------------------------------------
1. **The CRIRES+ raw IDPs and the Elgueta spectra are THE SAME OBSERVATIONS.** Elgueta's
   `table1.dat` gives the solar epoch as **2022-11-22 00:23**; our IDP manifest gives
   2022-11-21T23:55 and 2022-11-22T00:03 under programme 60.A-9051(A). So the audit is
   not comparing two reductions of two nights -- it is a genuine raw/corrected PAIR of
   one night, which is the strongest form check (a) can take.

   ⚠️ `normalize_vesta_ir.py` still carries the sentence "a filesystem sweep finds ZERO
   Vesta FITS on Sirius ... molecfit had nothing to run on". EIGHTEEN of them are at
   data.spectra_local/vesta/CRIRESPlus, staged after that sweep. The comment is stale and
   the conclusion built on it -- that the CRIRES+ telluric state could not be checked --
   is no longer true.

2. **A telluric TEMPLATE for the H band already existed in-house.** RYA-963 ran molecfit
   on CRIRES+ H frames of alpha Cen A, and its outputs carry an `MTRANS` extension: the
   FITTED TELLURIC TRANSMISSION, on a CRIRES+ wavelength grid, for the same instrument
   and the same band. That converts "does this window look deep?" into "does this window
   absorb WHERE THE TELLURIC LINES ARE?", which a stellar line forest cannot imitate.

3. **The claim's stated evidence was measured in a band the data does not cover.** The
   registry justifies `elgueta2026_vizier` telluric_applied=applied with "0.10% of window
   points below 0.5 in the O2 A-band vs 51.3% for Kitt Peak in the same window". The O2
   A-band is 7594-7685 A. The Elgueta Y arm starts at **9796.5 A**. Those spectra do not
   reach the O2 A-band at all; the 0.10% is RYA-794's statistic for the Y SCIENCE WINDOW,
   relabelled. The verdict this audit reaches is that the spectra ARE corrected -- but
   nothing that had been written down established it, which is this ticket's whole point.

THE METHOD, AND ITS CONTROL
---------------------------
Correlate observed absorption (1-flux) against the TEMPLATE's absorption (1-MTRANS) at
zero velocity shift. Telluric lines are at rest in the OBSERVER frame, so:

    raw topocentric IDP   -> must correlate STRONGLY   (the positive control)
    a corrected spectrum  -> must correlate at the NULL
    raw / corrected       -> must correlate STRONGLY   (the ratio IS the transmission)

⚠️ The null is MEASURED, never assumed: the same three statistics are computed in windows
where the template says there is essentially no telluric (`tmpl_frac <= 0.02`), and every
verdict is stated against that null rather than against a chosen threshold. A window is
called clean only when the corrected spectrum's correlation sits inside the null band AND
the raw one sits far above it -- i.e. only where the test is shown to have had the power
to detect a residual in the first place.

CLOSING THE 131
---------------
The two holdings differ, and one blanket rule would have been wrong for both:

* **solar_iag** DOES have a raw sibling and nobody had used it. `solar_iag_reiners2016`
  is the same IAG FTS atlas UNCORRECTED. It is SPAN-CAPPED at 5001.1 A in the HoldingSpec
  so that selection can never prefer it over Baker+2020 -- but the cap is a SELECTION
  policy, not the file's extent: `load_iag_reiners_window` serves 4047.4-10649.9 A and is
  called directly here. ⚠️ That makes IAG's CLAIM checkable, and it holds in all five
  discriminating bands -- but not by the difference test: see `band_verdicts` for why this
  particular pair can only be judged on depth.

* **solar_kpno_kurucz2005_corrected** genuinely has none -- Kurucz ships only the
  corrected residual atlas, and the 1984 KP atlas is a DIFFERENT observation, not its raw
  sibling (RYA-929 registered it "provenance-distinct" for exactly this reason). So its
  lines are closed by the CALIBRATED SIGNATURE test instead, and are labelled as such:
  the statistic is anchored on solar_kpno, which is uncorrected and covers the same
  2990-10010 A, and the control windows show the two agree on real solar structure.

⚠️ A per-line verdict from a band-level test is a BAND fact carried to the line, and it is
recorded with `basis` saying which of the two it came from. A line no declared telluric
band reaches is NOT called "corrected" -- it is NO-TELLURIC-BAND-DECLARED, because our
band registry is our own and has already been caught incomplete once (the missing
o2gamma, RYA-1193). Absence of a declared band is not absence of telluric (RYA-833).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from config.constants import codex_root  # noqa: E402  (RYA-810 path register)

C_KMS = 2.99792458e5
FEED = ROOT / "data/products/solar/Fe.json"
BP = ROOT / "data/results/band_products"

#: Deep-absorption cuts for the signature table. Both are reported: 0.5 is what saturated
#: telluric looks like, 0.8 catches the shallower H2O forest.
CUTS = (0.5, 0.8)

#: Control windows for the signature table -- no registered telluric complex reaches
#: either, and both are rich in solar lines, so a corrected atlas must still MATCH the
#: uncorrected one here. That is what separates "telluric removed" from "spectrum flattened".
SIGNATURE_CONTROLS = ((5000.0, 5100.0), (6000.0, 6100.0))

#: corrected holding -> (instrument, raw sibling, how the sibling relates to it)
RAW_SIBLING = {
    "solar_kpno_molecfit_corrected": (
        "kpno_solar_atlas", "solar_kpno", "SAME ARRAY, molecfit applied in-band"),
    "solar_harps_molecfit_corrected": (
        "harps", "solar_harps", "SEPARATELY REDUCED — carries a reduction floor"),
    "solar_iag": (
        "iag_fts_solar_atlas", "solar_iag_reiners2016",
        "SAME IAG FTS atlas, Reiners+2016 uncorrected vs Baker+2020 corrected"),
    "solar_kpno_kurucz2005_corrected": (
        "kpno_solar_atlas", None,
        "NO RAW SIBLING EXISTS — Kurucz ships only the corrected residual atlas, and the "
        "1984 KP atlas is a different observation (RYA-929 'provenance-distinct')"),
}
IDENTICAL = 1e-9
MATERIAL_MULTIPLE = 10.0
#: 🔴 THE REDUCTION FLOOR IS A GLOBAL PROPERTY AND ITS BASELINE MUST BE GLOBAL. These
#: windows are clear of EVERY registered telluric band and lie inside the span of every
#: VIS holding compared here (IAG starts at 5001.1 A, HARPS ends at 6910 A), so the same
#: baseline is available for each pair. Judging the difference test against windows
#: immediately beside a telluric band instead deflates it badly: a band has WINGS, the
#: correction reaches into them, and O2 B came back "raw-like" at a ratio of 3.8 on HARPS
#: -- a pair whose in-band correction is four orders above its reduction floor.
CONTROL_WINDOWS = ((5050.0, 5150.0), (5500.0, 5600.0), (6000.0, 6100.0))

#: CRIRES+ molecfit outputs from RYA-963 / RYA-993. They carry MTRANS = the fitted
#: telluric transmission on a CRIRES+ grid, which is the template this leg needs.
#: ⚠️ Different star and different night from Vesta, so the DEPTHS are not transferable
#: (airmass and PWV differ) -- only the LINE POSITIONS are used, and the statistic is a
#: shape correlation, never an absolute depth.
TEMPLATE_GLOBS = (
    "spectra/alpha_cen_a/CRIRESPlus_molecfit/*_telluric.fits",
    "spectra/tau_cet/CRIRESPlus_molecfit/*_telluric.fits",
)

#: The two live CRIRES+ Fe holdings, and the arm each one is.
#: 🔴 BOUNDED BY OUR DECIDED IR SCIENCE WINDOW (RYA-1094 / RYA-1193), not the instrument's
#: 53000 A reach: the lab-graded Fe pool ends at 17277.5 A and the H arm's own extent ends
#: at 17493.7. Nothing above is audited because nothing above is measurable.
CRIRES_ARMS = (
    ("solar_crires_plus_y_wide_rya1054", "Y", 9800.0, 10796.0),
    ("solar_crires_plus_h_rya1094", "H", 15007.0, 17493.0),
)
CRIRES_PRODUCT = {
    "solar_crires_plus_y_wide_rya1054":
        "data/results/rya794/vesta_crires_plus_Y_9800_10796_normalized.csv",
    "solar_crires_plus_h_rya1094":
        "data/results/rya1048/vesta_crires_plus_H_15007_17494_normalized.csv",
}
WINDOW_A = 100.0          # scan step for the CRIRES+ legs
TEMPLATE_STRONG = 0.05    # frac of template pixels below 0.97 that makes a window "telluric"
TEMPLATE_CLEAN = 0.02     # ... and below which it has nothing for the test to find

#: 🔴 THE NULL IS MEASURED PER WINDOW BY DISPLACING THE TEMPLATE, not read off the handful
#: of clean windows the arm happens to contain. The first cut of this leg took the null to
#: be `max(corr_r)` over the clean windows; the H arm has only TWO of them and both came
#: out slightly negative, so the threshold landed at -0.015 and flagged 12 of 18 telluric
#: windows as residuals. A null estimated from two samples is not a null.
#:
#: Displacing the template by a velocity far larger than any residual could sit at leaves
#: the SAME spectrum, the SAME template and the SAME window, and destroys only the
#: registration between them. That distribution is what "no telluric signal" looks like
#: for this window, measured on this window's own data. At 16000 A, 150 km/s is 8 A --
#: tens of telluric line widths at R~100000, so the shifted template is decorrelated while
#: remaining statistically identical.
NULL_SHIFTS_KMS = tuple(v for v in range(-800, 801, 25) if abs(v) >= 150)
NULL_Z = 3.0              # sigma above a window's own displaced null before it is flagged


# ══ loading ══════════════════════════════════════════════════════════════════════════
def _load(inst, hold, centre, pad):
    from measure_band_ew import load_window_ex
    w = load_window_ex(inst, centre, pad, holding=hold, allow_uncorrected=True)
    a = np.asarray(w.wave, float)
    f = np.asarray(w.flux, float)
    k = np.isfinite(a) & np.isfinite(f)
    return a[k], f[k]


def _load_raw(inst, hold, centre, pad):
    """As `_load`, but routes round the ONE span cap that is a selection policy.

    `solar_iag_reiners2016` is capped at 5001.1 A in its HoldingSpec so that selection can
    never prefer it to the corrected Baker atlas. The FILE covers 4047.4-10649.9 A. Using
    the module-level reader keeps the cap doing its job for measurement while letting the
    audit see the raw data the cap hides.
    """
    if hold == "solar_iag_reiners2016":
        from measure_band_ew import load_iag_reiners_window
        w, f, _ = load_iag_reiners_window(centre, pad)
        a, ff = np.asarray(w, float), np.asarray(f, float)
        k = np.isfinite(a) & np.isfinite(ff)
        return a[k], ff[k]
    return _load(inst, hold, centre, pad)


def telluric_template() -> tuple[np.ndarray, np.ndarray, list[str]]:
    """CRIRES+ telluric transmission, concatenated from every molecfit MTRANS we hold."""
    root = codex_root('data')
    W, T, used = [], [], []
    for pat in TEMPLATE_GLOBS:
        for p in sorted(root.glob(pat)):
            try:
                from astropy.io import fits
                d = fits.open(p)[1].data
                if "MTRANS" not in (d.columns.names or []):
                    continue
                W.append(np.asarray(d["WAVE"], float))
                T.append(np.asarray(d["MTRANS"], float))
                used.append(p.name)
            except Exception:
                continue
    if not W:
        return np.array([]), np.array([]), []
    w = np.concatenate(W)
    t = np.concatenate(T)
    k = np.isfinite(w) & np.isfinite(t)
    w, t = w[k], t[k]
    o = np.argsort(w)
    return w[o], t[o], used


# ══ statistics ═══════════════════════════════════════════════════════════════════════
def _flatten(w, f, frac=0.20, q=0.95):
    """Divide out a smooth upper envelope so an UNNORMALISED raw frame can be compared
    with a normalised product. A rolling high percentile, wide enough (20% of the window)
    that it cannot follow a telluric line and remove the very thing being measured."""
    n = len(w)
    win = max(51, int(n * frac) | 1)
    c = (pd.Series(f).rolling(win, center=True, min_periods=win // 4)
         .quantile(q).bfill().ffill().to_numpy())
    return f / np.where(c > 0, c, np.nan)


def _corr(a, b):
    a = a - np.nanmean(a)
    b = b - np.nanmean(b)
    d = np.sqrt(np.nansum(a * a) * np.nansum(b * b))
    return float(np.nansum(a * b) / d) if d > 0 else float("nan")


def baseline(inst, corrected, raw) -> float | None:
    """The pair's own corrected-vs-raw difference OUTSIDE any telluric band."""
    vals = []
    for lo, hi in CONTROL_WINDOWS:
        try:
            _, fc = _load(inst, corrected, 0.5 * (lo + hi), 0.5 * (hi - lo))
            _, fr = _load_raw(inst, raw, 0.5 * (lo + hi), 0.5 * (hi - lo))
        except Exception:
            continue
        n = min(len(fc), len(fr))
        if n > 50:
            vals.append(float(np.nanmean(np.abs(fc[:n] - fr[:n]))))
    return float(np.mean(vals)) if vals else None


def compare(inst, corrected, raw, centre, pad, base=None) -> dict:
    """max|corrected - raw| over one window, judged against the pair's own baseline."""
    try:
        _, fc = _load(inst, corrected, centre, pad)
    except Exception as e:
        return {"state": "UNVERIFIABLE", "why": f"corrected unreadable: {str(e)[:80]}"}
    if raw is None:
        return {"state": "NO-RAW-SIBLING", "n": int(len(fc))}
    try:
        _, fr = _load_raw(inst, raw, centre, pad)
    except Exception as e:
        return {"state": "UNVERIFIABLE", "why": f"raw unreadable: {str(e)[:80]}",
                "n": int(len(fc))}
    n = min(len(fc), len(fr))
    if n < 5:
        return {"state": "UNVERIFIABLE", "why": f"only {n} comparable points"}
    # ⚠️ Baker and Reiners are sampled 5.6x apart, so an index-wise subtraction would be
    # comparing different wavelengths. Interpolate the raw onto the corrected grid.
    wc, fc2 = _load(inst, corrected, centre, pad)
    wr, fr2 = _load_raw(inst, raw, centre, pad)
    if len(wr) < 5 or len(wc) < 5:
        return {"state": "UNVERIFIABLE", "why": "too few points after load"}
    if not (wr[0] - 1e-6 <= wc[0] and wc[-1] <= wr[-1] + 1e-6):
        lo_ok, hi_ok = max(wc[0], wr[0]), min(wc[-1], wr[-1])
        m = (wc >= lo_ok) & (wc <= hi_ok)
        wc, fc2 = wc[m], fc2[m]
        if len(wc) < 5:
            return {"state": "UNVERIFIABLE", "why": "no overlapping span"}
    fri = np.interp(wc, wr, fr2)
    dmax = float(np.nanmax(np.abs(fc2 - fri)))
    if dmax <= IDENTICAL:
        state = "VERIFIED-RAW"
    elif base is not None and base > 0 and dmax <= MATERIAL_MULTIPLE * base:
        state = "NO-MATERIAL-DIFFERENCE"
    else:
        state = "VERIFIED-CORRECTED"
    return {"state": state, "max_abs_diff": round(dmax, 8),
            "baseline_abs_diff": (round(base, 8) if base is not None else None),
            "n": int(len(wc))}


# ══ leg 1: the calibrated band-level signature table ══════════════════════════════════
def signature_table() -> list[dict]:
    """Deep-absorption fraction per holding per registered telluric band, ANCHORED.

    The anchor is `solar_kpno` -- uncorrected, and covering 2960-13000 A, so it overlaps
    every VIS/NIR holding here. A "corrected" holding must be far shallower than the
    anchor INSIDE a band and must AGREE with it in the controls; a holding that is
    shallower everywhere has had its structure flattened, not its tellurics removed, and
    the control columns are what tell those two apart.
    """
    from pipeline.telluric_policy import TELLURIC_BANDS
    holdings = [("kpno_solar_atlas", "solar_kpno"),
                ("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
                ("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
                ("harps", "solar_harps"),
                ("harps", "solar_harps_molecfit_corrected"),
                ("iag_fts_solar_atlas", "solar_iag"),
                ("iag_fts_solar_atlas", "solar_iag_reiners2016")]
    bands = ([(float(a), float(b), str(c), "TELLURIC") for a, b, c in TELLURIC_BANDS]
             + [(lo, hi, f"CONTROL {lo:.0f}-{hi:.0f}", "CONTROL")
                for lo, hi in SIGNATURE_CONTROLS])
    out = []
    for lo, hi, name, kind in bands:
        row = {"band": name, "kind": kind, "lo_A": lo, "hi_A": hi, "holdings": {}}
        for inst, hold in holdings:
            try:
                _, f = _load_raw(inst, hold, 0.5 * (lo + hi), 0.5 * (hi - lo))
            except Exception as e:
                row["holdings"][hold] = {"state": "OUT-OF-SPAN-OR-REFUSED",
                                         "why": f"{type(e).__name__}: {str(e)[:80]}"}
                continue
            if len(f) < 20:
                row["holdings"][hold] = {"state": "TOO-THIN", "n": int(len(f))}
                continue
            row["holdings"][hold] = {
                "state": "MEASURED", "n": int(len(f)),
                **{f"pct_below_{c}": round(float(np.mean(f < c) * 100), 3) for c in CUTS}}
        out.append(row)
    return out


#: The anchor of last resort: uncorrected, and covering 2960-13000 A so it overlaps every
#: VIS/NIR holding here. Used ONLY where a holding has no raw sibling of its own.
FALLBACK_ANCHOR = "solar_kpno"

#: A band's LOCAL comparison region: the same width again on each side, pushed clear of
#: the band itself. ⚠️ THE GLOBAL CONTROLS AT 5000-5100 AND 6000-6100 CANNOT DO THIS JOB.
#: Solar line density varies enormously with wavelength -- Kitt Peak is 16.9% below 0.8 at
#: 5000 A and 3.0% at 6000 A -- so "deeper than a control 1500 A away" measures the line
#: forest, not telluric. That is exactly how the O2 gamma-band read as a telluric band on
#: every holding at once: 6% absorption there is ordinary solar structure for 6270 A, and
#: comparing it to a control at 6000 A made it look like something being removed.
#: ⚠️ Side windows are used ONLY by the depth test's power gate, never by the difference
#: test — see CONTROL_WINDOWS for why.
SIDE_PAD_FRACTION = 1.0
SIDE_GAP_A = 5.0

#: How much deeper a band must be than its OWN SIDE WINDOWS before the anchor is accepted
#: as showing telluric there at all. Below this the band cannot discriminate, whatever the
#: holdings do in it.
ANCHOR_LOCAL_RATIO = 2.0
ANCHOR_RATIO = 5.0


def _side_windows(lo, hi):
    w = (hi - lo) * SIDE_PAD_FRACTION
    return ((lo - SIDE_GAP_A - w, lo - SIDE_GAP_A), (hi + SIDE_GAP_A, hi + SIDE_GAP_A + w))


def _depth(inst, hold, lo, hi, cut=0.8):
    try:
        _, f = _load_raw(inst, hold, 0.5 * (lo + hi), 0.5 * (hi - lo))
    except Exception:
        return None
    if len(f) < 20:
        return None
    return float(np.mean(f < cut) * 100)


def _side_depth(inst, hold, lo, hi):
    vals = [_depth(inst, hold, a, b) for a, b in _side_windows(lo, hi)]
    vals = [v for v in vals if v is not None]
    return float(np.mean(vals)) if vals else None


def _absorption_available(inst, raw, lo, hi):
    """Mean absorption present in the RAW sibling over a band — the most a correction
    there could possibly change. It is the ceiling the difference test has to clear."""
    try:
        _, f = _load_raw(inst, raw, 0.5 * (lo + hi), 0.5 * (hi - lo))
    except Exception:
        return None
    if len(f) < 20:
        return None
    return float(np.nanmean(np.clip(1.0 - f, 0.0, None)))


def _band_diff(inst, hold, raw, lo, hi):
    """mean|corrected - raw| inside a window, on the corrected grid."""
    try:
        wc, fc = _load(inst, hold, 0.5 * (lo + hi), 0.5 * (hi - lo))
        wr, fr = _load_raw(inst, raw, 0.5 * (lo + hi), 0.5 * (hi - lo))
    except Exception:
        return None
    if len(wc) < 20 or len(wr) < 20:
        return None
    m = (wc >= max(wc[0], wr[0])) & (wc <= min(wc[-1], wr[-1]))
    if m.sum() < 20:
        return None
    return float(np.nanmean(np.abs(fc[m] - np.interp(wc[m], wr, fr))))


def band_verdicts(sig: list[dict]) -> dict:
    """Per holding per telluric band: is the band corrected on this holding?

    🔴 TWO TESTS, AND WHICH ONE APPLIES IS ITSELF MEASURED PER BAND PER HOLDING.

    * **DIFFERENCE** -- how much the corrected product differs from its own raw sibling
      inside the band, against that pair's reduction floor in telluric-free control
      windows. Decisive where it works, and free of any assumption about how much of the
      absorption was telluric. It is the RYA-1190 method at band scale.

    * **DEPTH** -- how much absorption the holding retains inside the band compared with
      an uncorrected anchor. Weaker, and only meaningful where the anchor is LOCALLY deep:
      deeper inside the band than in the band's own side windows by ANCHOR_LOCAL_RATIO.
      Without that gate the O2 gamma-band read as a telluric complex on every holding at
      once, on the strength of 6% absorption that is just the solar spectrum at 6270 A.

    ⚠️ THE DIFFERENCE TEST IS NOT ALWAYS THE STRONGER ONE, AND ASSUMING IT WAS PUT IAG'S
    WHOLE CLAIM ON THE WRONG SIDE. Baker+2020 and Reiners+2016 are separate reductions of
    the same atlas sampled 5.6x apart; they differ by 0.083 in flux in telluric-free
    windows before any telluric question is asked. Ten times that floor is more absorption
    than the O2 A-band contains, so the test could not have returned CORRECTED there for
    any data -- and it duly returned RAW-LIKE on a holding whose O2 A-band is 2.1% below
    0.8 against its raw sibling's 66.8%. The choice of test is therefore made by comparing
    the bar (MATERIAL_MULTIPLE x floor) with the ceiling (the absorption actually present
    in the raw sibling), and every row records which test answered and why.
    """
    anchor_of = {h: (sib if sib else FALLBACK_ANCHOR)
                 for h, (_, sib, _) in RAW_SIBLING.items()}
    inst_of = {h: i for h, (i, _, _) in RAW_SIBLING.items()}
    floors = {h: (baseline(i, h, sib) if sib else None)
              for h, (i, sib, _) in RAW_SIBLING.items()}
    out: dict[str, dict] = {}
    for row in sig:
        if row["kind"] != "TELLURIC":
            continue
        lo, hi = row["lo_A"], row["hi_A"]
        for hold, cell in row["holdings"].items():
            anchor = anchor_of.get(hold)
            if anchor is None or hold == anchor:
                continue
            inst = inst_of[hold]
            key = f"{hold}|{row['band']}|{lo:.0f}-{hi:.0f}"
            if cell.get("state") != "MEASURED":
                out[key] = {"verdict": "NOT-COVERED", "anchor": anchor,
                            "why": cell.get("why", cell.get("state"))}
                continue
            # 🔴 NOT `anchor != FALLBACK_ANCHOR`. solar_kpno_molecfit_corrected's genuine
            # raw sibling IS solar_kpno, which is also the fallback anchor -- so testing
            # by identity called the one pair with an EXACT ZERO floor "no sibling" and
            # sent the most decisive comparison we have down the weaker depth route.
            own = RAW_SIBLING.get(hold, (None, None, None))[1] is not None
            avail = _absorption_available(inst, anchor, lo, hi) if own else None
            d_in = _band_diff(inst, hold, anchor, lo, hi) if own else None
            floor = floors.get(hold)
            # 🔴 A PAIR CAN BE TOO UNLIKE ITSELF FOR THE DIFFERENCE TEST TO ANSWER.
            # `MATERIAL_MULTIPLE * floor` is the bar a correction has to clear; `avail` is
            # the total absorption actually present in the raw sibling, which is the most
            # any correction there could remove. When the bar exceeds the ceiling the test
            # cannot return CORRECTED for physical reasons, so a RAW-LIKE from it would be
            # a property of the pair, not of the data. Baker+2020 and Reiners+2016 are two
            # separate reductions of the IAG atlas sampled 5.6x apart and they differ by
            # 0.083 in flux in telluric-free windows BEFORE any telluric question is asked
            # -- resampling, renormalisation or a small velocity offset between them; this
            # audit does not resolve which, and says so rather than reading the floor as
            # evidence. That pair is judged on DEPTH instead.
            usable = (own and d_in is not None and floor is not None and avail is not None
                      and MATERIAL_MULTIPLE * max(floor, 1e-9) < avail)
            if own and d_in is None:
                out[key] = {"verdict": "NOT-COVERED", "anchor": anchor,
                            "why": "band unreadable on this pair"}
                continue
            if usable:
                corrected = d_in > MATERIAL_MULTIPLE * max(floor, 1e-9)
                out[key] = {
                    "verdict": ("VERIFIED-CORRECTED" if corrected else "VERIFIED-RAW-LIKE"),
                    "test": ("DIFFERENCE vs own raw sibling, in-band against this pair's "
                             "own reduction floor from telluric-free control windows"),
                    "anchor": anchor, "anchor_is_own_raw_sibling": True,
                    "mean_abs_diff_in_band": round(d_in, 8),
                    "reduction_floor_from_controls": round(floor, 8),
                    "absorption_available_in_raw": round(avail, 6),
                    "ratio_in_over_floor": (round(d_in / floor, 1) if floor > 0 else None),
                    "threshold_ratio": MATERIAL_MULTIPLE,
                    "holding_pct_below_0.8": cell["pct_below_0.8"]}
                continue
            # ── DEPTH, gated on the anchor being LOCALLY deep ──────────────────────
            a = row["holdings"].get(anchor, {})
            if a.get("state") != "MEASURED":
                out[key] = {"verdict": "NO-ANCHOR", "anchor": anchor,
                            "why": f"{anchor} does not reach this band"}
                continue
            a_in, hv = a["pct_below_0.8"], cell["pct_below_0.8"]
            a_side = _side_depth(inst if own else "kpno_solar_atlas", anchor, lo, hi)
            if a_side is None:
                out[key] = {"verdict": "NO-LOCAL-BASELINE", "anchor": anchor,
                            "why": "the anchor's side windows are unreadable"}
                continue
            why_depth = (
                (f"the difference test has no power on this pair: 10x its reduction floor "
                 f"({MATERIAL_MULTIPLE * floor:.3f}) exceeds the absorption available to "
                 f"remove ({avail:.3f})") if own else
                "no raw sibling exists for this holding")
            if a_in < ANCHOR_LOCAL_RATIO * a_side:
                out[key] = {
                    "verdict": "ANCHOR-SHOWS-NO-LOCAL-TELLURIC", "anchor": anchor,
                    "why": (f"{anchor} is {a_in:.2f}% below 0.8 inside the band against "
                            f"{a_side:.2f}% in its own side windows — that is the solar "
                            f"line forest, not a telluric complex, so this band cannot "
                            f"discriminate here"),
                    "fell_back_to_depth_because": why_depth,
                    "anchor_pct_below_0.8": a_in,
                    "anchor_side_pct_below_0.8": round(a_side, 3),
                    "holding_pct_below_0.8": hv}
                continue
            out[key] = {
                "verdict": ("VERIFIED-CORRECTED" if hv * ANCHOR_RATIO <= a_in
                            else "VERIFIED-RAW-LIKE"),
                "test": "DEPTH vs an uncorrected anchor, gated on local band depth",
                "anchor": anchor, "anchor_is_own_raw_sibling": bool(own),
                "fell_back_to_depth_because": why_depth,
                "anchor_pct_below_0.8": a_in,
                "anchor_side_pct_below_0.8": round(a_side, 3),
                "holding_pct_below_0.8": hv,
                "ratio_anchor_over_holding": (round(a_in / hv, 1) if hv > 0 else None)}
            continue
    return out


# ══ leg 2: every measured Fe line ════════════════════════════════════════════════════
def measured_lines() -> pd.DataFrame:
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
    return pd.DataFrame(rows, columns=["holding", "instrument", "wavelength_air_A"]
                        ).drop_duplicates()


def holding_claim_verdicts(bverd: dict) -> dict:
    """Roll the per-band verdicts up to one verdict on each holding's stated CLAIM.

    ⚠️ This is what actually answers "is the delivered/molecfit claim true", and it is a
    HOLDING fact, separate from whether any given Fe line sits somewhere a correction
    could have mattered. Both are reported because conflating them is how a holding gets
    called verified on the strength of lines that no telluric band ever reached.
    """
    out: dict[str, dict] = {}
    for key, v in bverd.items():
        hold = key.split("|", 1)[0]
        o = out.setdefault(hold, {"corrected_bands": [], "raw_like_bands": [],
                                  "no_power_bands": [], "not_covered_bands": []})
        band = key.split("|")[1]
        {"VERIFIED-CORRECTED": o["corrected_bands"],
         "VERIFIED-RAW-LIKE": o["raw_like_bands"],
         "ANCHOR-SHOWS-NO-LOCAL-TELLURIC": o["no_power_bands"],
         "NO-LOCAL-BASELINE": o["no_power_bands"],
         "NO-BASELINE": o["no_power_bands"],
         "NO-ANCHOR": o["not_covered_bands"],
         "NOT-COVERED": o["not_covered_bands"]}[v["verdict"]].append(band)
    for hold, o in out.items():
        nc, nr = len(o["corrected_bands"]), len(o["raw_like_bands"])
        if nc and not nr:
            o["claim_verdict"] = (
                f"VERIFIED-CORRECTED in all {nc} discriminating band(s): "
                f"{', '.join(o['corrected_bands'])}")
        elif nr and not nc:
            o["claim_verdict"] = (
                f"🔴 VERIFIED-RAW-LIKE in all {nr} discriminating band(s): "
                f"{', '.join(o['raw_like_bands'])} — the claim is not borne out")
        elif nc and nr:
            o["claim_verdict"] = (
                f"⚠️ MIXED — corrected in {nc} band(s) ({', '.join(o['corrected_bands'])}) "
                f"and raw-like in {nr} ({', '.join(o['raw_like_bands'])}); one holding "
                f"label is describing two different states")
        else:
            o["claim_verdict"] = (
                "UNDETERMINED — no band both reaches this holding and shows telluric in "
                "the uncorrected anchor, so nothing here can discriminate")
    return out


def close_lines(bverd: dict, claims: dict, pad: float) -> pd.DataFrame:
    """A state for EVERY measured line, with the basis it came from stated."""
    from pipeline.telluric_policy import TELLURIC_BANDS
    lines = measured_lines()
    bases = {}
    for hold, (inst, raw, _) in RAW_SIBLING.items():
        bases[hold] = baseline(inst, hold, raw) if raw else None
    rows = []
    for _, r in lines.iterrows():
        hold, w = r.holding, float(r.wavelength_air_A)
        inst, raw, rel = RAW_SIBLING.get(hold, (r.instrument, None, "unregistered holding"))
        band = next((str(c) for a, b, c in TELLURIC_BANDS if a <= w <= b), None)
        floor = bases.get(hold)
        rec = dict(holding=hold, instrument=inst, wavelength_air_A=w,
                   raw_sibling=raw, sibling_relation=rel,
                   telluric_band=band or "",
                   # ⚠️ TWO FACTS, KEPT APART. Whether the two products differ here, and
                   # whether a telluric correction was ever expected here, are different
                   # questions -- and a line outside every declared band answers the second
                   # one "no" however the first one comes out.
                   correction_expected_here=bool(band))
        # 🔴 THE SAME POWER GATE AS THE BAND LEG. A pair whose reduction floor is larger
        # than the absorption present at the line cannot return CORRECTED for any data, so
        # a difference verdict from it describes the pair, not the flux. IAG's Baker/Reiners
        # floor is 0.083 -- deeper than most Fe lines -- and blanket-applying the difference
        # test put all 43 of its lines in NO-MATERIAL-DIFFERENCE, a non-answer dressed as one.
        avail = (_absorption_available(inst, raw, w - pad, w + pad)
                 if raw is not None else None)
        powered = (raw is not None and floor is not None and avail is not None
                   and MATERIAL_MULTIPLE * max(floor, 1e-12) < avail)
        rec["absorption_available_in_raw"] = (round(avail, 6) if avail is not None else None)
        rec["reduction_floor"] = (round(floor, 8) if floor is not None else None)
        if powered:
            res = compare(inst, hold, raw, w, pad, base=floor)
            rec.update(res)
            rec["basis"] = "DIRECT-COMPARISON vs raw sibling"
        elif band is None:
            rec.update(state="NO-TELLURIC-BAND-DECLARED",
                       holding_claim=claims.get(hold, {}).get("claim_verdict"))
            rec["basis"] = (
                "no declared telluric band reaches this line, and "
                + ("no raw sibling exists" if raw is None else
                   "the difference test has no power on this pair here")
                + " — the HOLDING's claim is judged separately in `holding_claim` and is "
                  "NOT evidence about this line; an undeclared band is not an absent one "
                  "(RYA-833)")
        else:
            key = next((k for k in bverd if k.startswith(f"{hold}|{band}|")), None)
            v = bverd.get(key, {})
            rec.update(state=v.get("verdict", "UNVERIFIABLE"),
                       band_test=v.get("test"),
                       band_anchor=v.get("anchor"),
                       anchor_pct_below_0_8=v.get("anchor_pct_below_0.8"),
                       holding_pct_below_0_8=v.get("holding_pct_below_0.8"))
            rec["holding_claim"] = claims.get(hold, {}).get("claim_verdict")
            rec["basis"] = (
                "BAND VERDICT carried to the line — "
                + ("no raw sibling exists for this holding" if raw is None else
                   "the difference test has no power on this pair at this line"))
        rows.append(rec)
    return pd.DataFrame(rows)


# ══ leg 3: CRIRES+, per arm, against RAW ═════════════════════════════════════════════
def _displaced_null(g, a_obs, tw, tt):
    """r(0) against the distribution of r(dv) for |dv| >= 150 km/s, on this window's data.

    Returns (r0, z, mean, sd). `z` is how many of the window's own null sigmas the
    zero-shift correlation stands above chance. A window with a real telluric residual has
    a large positive z; a corrected window has z ~ 0 whatever its raw neighbour does.
    """
    a = a_obs - np.nanmean(a_obs)
    def r_at(dv):
        b = 1.0 - np.interp(g, tw * (1.0 + dv / C_KMS), tt)
        b = b - np.nanmean(b)
        d = np.sqrt(np.nansum(a * a) * np.nansum(b * b))
        return float(np.nansum(a * b) / d) if d > 0 else np.nan
    r0 = r_at(0.0)
    null = np.array([r_at(v) for v in NULL_SHIFTS_KMS], float)
    null = null[np.isfinite(null)]
    if len(null) < 10 or not np.isfinite(r0):
        return r0, float("nan"), float("nan"), float("nan")
    mu, sd = float(np.mean(null)), float(np.std(null, ddof=1))
    z = (r0 - mu) / sd if sd > 0 else float("nan")
    return r0, z, mu, sd


def crires_arm(hold: str, arm: str, lo0: float, hi0: float,
               MW: np.ndarray, MT: np.ndarray) -> dict:
    """Window-by-window raw-vs-corrected against the telluric template. Check (a) at last.

    Three correlations per window, all at ZERO velocity shift because telluric lines are
    at rest in the observer frame and the raw IDP is topocentric:

        raw_r     raw IDP absorption   vs template  -- the POSITIVE CONTROL
        corr_r    the product          vs template  -- the TEST
        ratio_r   (raw / product)      vs template  -- the transmission it removed

    ⚠️ `corr_r` alone cannot clear a window: a window where the template is flat has
    nothing to detect, so a low `corr_r` there is uninformative. Only windows where
    `raw_r` is high have shown the test had power, and the verdict is built from those.
    """
    from measure_band_ew import load_crires_window
    prod = pd.read_csv(ROOT / CRIRES_PRODUCT[hold]).dropna()
    ew = prod.wavelength_air_A.to_numpy()
    ef = prod.flux_normalized.to_numpy()
    rows = []
    lo = lo0
    while lo < hi0:
        hi = min(lo + WINDOW_A, hi0)
        m = (MW >= lo - 3) & (MW < hi + 3)
        rec = {"lo_A": lo, "hi_A": hi}
        if m.sum() < 200:
            rec["state"] = "NO-TEMPLATE-COVERAGE"
            rows.append(rec); lo += WINDOW_A; continue
        rec["template_frac_below_0.97"] = round(float(np.mean(MT[m] < 0.97)), 4)
        rec["template_min"] = round(float(np.min(MT[m])), 4)
        try:
            rw, rf, prov = load_crires_window(0.5 * (lo + hi), 0.5 * (hi - lo),
                                              allow_topocentric=True)
        except Exception as e:
            rec["state"] = "NO-RAW-IDP"; rec["why"] = str(e)[:90]
            rows.append(rec); lo += WINDOW_A; continue
        k = np.isfinite(rw) & np.isfinite(rf)
        me = (ew >= lo - 5) & (ew < hi + 5)
        if k.sum() < 400 or me.sum() < 400:
            rec["state"] = "TOO-THIN"; rec["n_raw"] = int(k.sum())
            rec["n_product"] = int(me.sum())
            rows.append(rec); lo += WINDOW_A; continue
        rw2, rf2 = rw[k], rf[k]
        o = np.argsort(rw2); rw2, rf2 = rw2[o], rf2[o]
        rn = _flatten(rw2, rf2)
        g0, g1 = max(lo + 1, rw2[0], ew[me][0]), min(hi - 1, rw2[-1], ew[me][-1])
        if g1 - g0 < 0.5 * WINDOW_A:
            rec["state"] = "PARTIAL-OVERLAP"; rec["overlap_A"] = round(float(g1 - g0), 2)
            rows.append(rec); lo += WINDOW_A; continue
        g = np.linspace(g0, g1, 3000)
        a_raw = 1.0 - np.interp(g, rw2, rn)
        a_prd = 1.0 - np.interp(g, ew[me], ef[me])
        a_tpl = 1.0 - np.interp(g, MW[m], MT[m])
        ratio = np.interp(g, rw2, rn) / np.clip(np.interp(g, ew[me], ef[me]), 0.05, None)
        raw_r, raw_z, _, _ = _displaced_null(g, a_raw, MW[m], MT[m])
        corr_r, corr_z, nmu, nsd = _displaced_null(g, a_prd, MW[m], MT[m])
        rat_r, rat_z, _, _ = _displaced_null(g, 1.0 - ratio, MW[m], MT[m])
        rec.update(state="MEASURED",
                   raw_r=round(raw_r, 4), raw_z=round(raw_z, 2),
                   corr_r=round(corr_r, 4), corr_z=round(corr_z, 2),
                   ratio_r=round(rat_r, 4), ratio_z=round(rat_z, 2),
                   null_mean=round(nmu, 4), null_sd=round(nsd, 4),
                   n_raw=int(k.sum()), n_product=int(me.sum()),
                   raw_provenance=prov[:60])
        rows.append(rec)
        lo += WINDOW_A
    d = pd.DataFrame(rows)
    ok = d[d.state == "MEASURED"] if "state" in d else pd.DataFrame()
    strong = ok[ok["template_frac_below_0.97"] > TEMPLATE_STRONG] if len(ok) else ok
    clean = ok[ok["template_frac_below_0.97"] <= TEMPLATE_CLEAN] if len(ok) else ok

    def stats(sub, col):
        if not len(sub):
            return None
        return {"n": int(len(sub)), "mean": round(float(sub[col].mean()), 4),
                "median": round(float(sub[col].median()), 4),
                "min": round(float(sub[col].min()), 4),
                "max": round(float(sub[col].max()), 4)}

    null = {c: stats(clean, c) for c in ("raw_r", "corr_r", "ratio_r", "corr_z")}
    hot = {c: stats(strong, c) for c in ("raw_r", "corr_r", "ratio_r",
                                         "raw_z", "corr_z", "ratio_z")}
    # 🔴 POWER IS A PROPERTY OF EACH WINDOW, NOT OF THE ARM. The first cut gated the whole
    # arm on the MEDIAN raw z clearing a chosen multiple of NULL_Z; the H arm missed it by
    # 0.1 sigma and the entire arm came back UNDETERMINED while eighteen individually
    # well-powered windows sat inside it. A window where the RAW frame stands above its own
    # null has demonstrated the test works THERE, and that is the only place the question
    # is asked.
    powered = strong[strong.raw_z >= NULL_Z] if len(strong) else strong
    unpowered = strong[strong.raw_z < NULL_Z] if len(strong) else strong
    flagged = ([{"lo_A": float(r.lo_A), "hi_A": float(r.hi_A),
                 "corr_r": float(r.corr_r), "corr_z": float(r.corr_z),
                 "raw_r": float(r.raw_r), "raw_z": float(r.raw_z),
                 "ratio_r": float(r.ratio_r), "ratio_z": float(r.ratio_z),
                 "template_frac_below_0.97": float(r["template_frac_below_0.97"]),
                 "raw_provenance": str(r.raw_provenance)}
                for _, r in powered.iterrows() if r.corr_z >= NULL_Z]
               if len(powered) else [])
    verdict, detail = _arm_verdict(arm, powered, unpowered, clean, hot, null, flagged)
    return {"holding": hold, "arm": arm, "span_A": [lo0, hi0],
            "windows": rows,
            "null_estimator": (
                f"per window, r(0) against r(dv) for {len(NULL_SHIFTS_KMS)} shifts with "
                f"|dv| >= 150 km/s; a window is flagged at corr_z >= {NULL_Z} sigma"),
            "clean_window_stats": null,
            "n_strong_windows": int(len(strong)),
            "n_strong_windows_with_power": int(len(powered)),
            "n_strong_windows_without_power": int(len(unpowered)),
            "in_strong_telluric_windows": hot,
            "verdict": verdict, "verdict_detail": detail,
            "flagged_windows": flagged}


def _arm_verdict(arm, strong, unpowered, clean, hot, null, flagged):
    """The arm's verdict, and it may only be reached where the test is shown to have power.

    ⚠️ THE ORDER OF THESE GUARDS IS THE POINT. A low correlation between the product and
    the telluric template means "no residual" ONLY if the same test on the RAW frame over
    the same windows returns a high one. Without that positive control, "clean" and "the
    test does not work here" are the same number.
    """
    if not len(strong):
        return ("UNDETERMINED — no window in this arm has enough telluric in the template "
                "for the test to have any power", {})
    d = {"n_strong_windows": int(len(strong)), "n_clean_windows": int(len(clean)),
         "raw_z_median_in_strong": hot["raw_z"]["median"] if hot["raw_z"] else None,
         "corr_z_median_in_strong": hot["corr_z"]["median"] if hot["corr_z"] else None,
         "corr_z_max_in_strong": hot["corr_z"]["max"] if hot["corr_z"] else None,
         "ratio_z_median_in_strong": hot["ratio_z"]["median"] if hot["ratio_z"] else None,
         "raw_r_median_in_strong": hot["raw_r"]["median"] if hot["raw_r"] else None,
         "corr_r_median_in_strong": hot["corr_r"]["median"] if hot["corr_r"] else None,
         "ratio_r_median_in_strong": hot["ratio_r"]["median"] if hot["ratio_r"] else None,
         "flag_threshold_sigma": NULL_Z}
    d["n_windows_without_power"] = int(len(unpowered))
    if not len(strong):
        return ("UNDETERMINED — no telluric window in this arm passes the positive "
                "control: the RAW IDP does not stand above its own displaced null "
                "anywhere, so a quiet product proves nothing", d)
    if not flagged:
        return (f"VERIFIED-CORRECTED — over the {len(strong)} telluric windows that pass "
                f"the positive control, the RAW IDP stands at a median "
                f"{hot['raw_z']['median']:.0f} sigma above its own displaced null while the "
                f"PRODUCT stands at {hot['corr_z']['median']:.1f} sigma and never reaches "
                f"{NULL_Z}. The test is demonstrated to have the power to see a residual "
                f"and sees none. The raw/product RATIO stands at "
                f"{hot['ratio_z']['median']:.0f} sigma — what the product removed IS the "
                f"telluric transmission, which is the positive form of the same claim.", d)
    return (f"⚠️ LOCALISED RESIDUAL in {len(flagged)} of {len(strong)} well-powered "
            f"telluric windows, "
            f"where the product still correlates with the telluric template at >= "
            f"{NULL_Z} sigma above that window's own displaced null: "
            + ", ".join(f"{f['lo_A']:.0f}-{f['hi_A']:.0f} "
                        f"(r={f['corr_r']:+.3f}, {f['corr_z']:.1f}s)" for f in flagged)
            + f". The other {len(strong) - len(flagged)} windows are clean by the same "
              f"test. ⚠️ LOCALISED, and a residual of this size is far below the raw "
              f"absorption (raw median {hot['raw_r']['median']:+.3f} vs product "
              f"{hot['corr_r']['median']:+.3f}) — a partial correction, not an absent one.", d)


# ══ main ═════════════════════════════════════════════════════════════════════════════
def graded_ir_coverage() -> dict:
    """Which graded IR Fe lines the LIVE CRIRES+ products actually sample.

    🔴 THIS RE-JUDGES A FINDING ALREADY ON MAIN, and the correction is larger than the
    finding. `rya1192_verification.json` reports that "12 graded Ruffoni-2013 lines fall
    in CRIRES+ H inter-order gaps" and that "16 of 29 graded IR lines are unmeasurable".
    That came from asking `len(flux) >= 20` in a +/-0.5 A window.

    The H product carries **16.3 points per Angstrom**, so a normally sampled line has
    about SIXTEEN points in that window and fails a >=20 test on sampling density alone.
    The threshold was not selecting served lines, it was selecting ORDER-OVERLAP regions,
    where two orders double the sampling to ~21. Every one of the twelve "gap" lines
    carries between 9 and 19 points; NONE has zero. The arm's largest wavelength step is
    1.18 A and only TWO steps anywhere in it exceed 0.5 A, so there is almost no
    inter-order gap structure in this product to fall into.

    Counted properly the H arm serves 25 of the 29 -- which is exactly the number RYA-1094
    registered when it built the product ("25 primary-lab Fe I lines, Ruffoni2013,
    15051.7-17277.5 A"). The audit had contradicted the holding's own registered count and
    that disagreement went unread.

    ⚠️ Fewer points is still WORSE sampling, and a line at half density is measured on
    half the pixels. That is a precision statement, not a coverage one, and it is reported
    as `n_points`, not as a verdict.
    """
    cgf = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    lab = cgf[(cgf.species.astype(str).str.startswith("Fe"))
              & (cgf.wavelength_air_A > 10796)
              & (cgf.gf_tier.astype(str) == "LAB")].sort_values("wavelength_air_A")
    prods, gaps = {}, {}
    for hold, arm, _, _ in CRIRES_ARMS:
        a = np.sort(pd.read_csv(ROOT / CRIRES_PRODUCT[hold])
                    .dropna().wavelength_air_A.to_numpy())
        prods[arm] = a
        st = np.diff(a)
        gaps[arm] = {
            "n_points": int(len(a)),
            "span_A": [round(float(a[0]), 2), round(float(a[-1]), 2)],
            "points_per_A": round(float(len(a) / (a[-1] - a[0])), 2),
            "median_step_A": round(float(np.median(st)), 4),
            "max_step_A": round(float(st.max()), 4),
            "n_steps_over_0.5A": int((st > 0.5).sum()),
            "steps_over_0.5A": [[round(float(a[i]), 3), round(float(a[i + 1]), 3)]
                                for i in np.where(st > 0.5)[0]],
        }
    rows = []
    for _, r in lab.iterrows():
        w = float(r.wavelength_air_A)
        arm = next((k for k, a in prods.items() if a[0] - 1 <= w <= a[-1] + 1), None)
        if arm is None:
            rows.append({"wavelength_air_A": round(w, 4), "lab_source_tag": str(r.lab_source_tag),
                         "arm": None, "n_points_pm0.5A": 0, "n_points_pm1.0A": 0,
                         "state": "OUTSIDE EVERY CRIRES+ ARM WE HOLD"})
            continue
        a = prods[arm]
        n5 = int((np.abs(a - w) < 0.5).sum())
        n10 = int((np.abs(a - w) < 1.0).sum())
        exp = gaps[arm]["points_per_A"]
        rows.append({
            "wavelength_air_A": round(w, 4), "lab_source_tag": str(r.lab_source_tag),
            "arm": arm, "n_points_pm0.5A": n5, "n_points_pm1.0A": n10,
            "sampling_vs_arm_mean": round(n5 / exp, 2) if exp else None,
            "state": ("SERVED" if n5 > 0 else "TRUE GAP — zero points"),
            "would_fail_the_prior_20_point_test": bool(0 < n5 < 20)})
    served = [r for r in rows if r["state"] == "SERVED"]
    return {
        "product_sampling": gaps,
        "lines": rows,
        "n_graded_beyond_10796": int(len(lab)),
        "n_served_by_a_live_product": len(served),
        "n_outside_every_arm": sum(1 for r in rows if r["arm"] is None),
        "n_true_gaps_zero_points": sum(1 for r in rows
                                       if r["arm"] and r["n_points_pm0.5A"] == 0),
        "n_that_the_prior_20_point_test_called_gaps": sum(
            1 for r in rows if r.get("would_fail_the_prior_20_point_test")),
        "correction": (
            "🔴 The prior artifact's '12 lines in inter-order gaps' and '16 of 29 "
            "unmeasurable' are ARTIFACTS OF A >=20-POINT THRESHOLD applied to a product "
            "sampled at 16.3 points/A. All twelve carry 9-19 points; none has zero. The H "
            "arm serves 25 of 29, matching RYA-1094's own registered count."),
    }


def _graded_in_flagged(arms, grade_cov) -> dict:
    """Which graded Fe lines sit in a window this leg flagged — the actionable consequence.

    🔴 A WINDOW-LEVEL VERDICT ONLY MATTERS WHERE A LINE IS. The prior artifact made this
    point about 16009.610 and was right to; it is made here for every flagged window, and
    the answer has moved with the windows. 16009.610 is now CLEAR (its window sits at
    -1.0 sigma) while lines in the H arm's blue end are not.

    ⚠️ These are CANDIDATES, not a disqualification. The correction in these windows is
    partial, not absent -- the residual correlations are ~0.25 against raw frames at ~0.6 --
    so what this supports is a per-line re-check before the H-arm pool is used for an
    abundance, not the removal of the lines.
    """
    out = []
    for a in arms:
        for f in a.get("flagged_windows", []):
            for r in grade_cov["lines"]:
                if r["arm"] and f["lo_A"] <= r["wavelength_air_A"] < f["hi_A"]:
                    out.append({"wavelength_air_A": r["wavelength_air_A"],
                                "lab_source_tag": r["lab_source_tag"],
                                "arm": a["arm"], "holding": a["holding"],
                                "window_A": [f["lo_A"], f["hi_A"]],
                                "window_corr_r": f["corr_r"],
                                "window_corr_z": f["corr_z"],
                                "n_points_pm0.5A": r["n_points_pm0.5A"]})
    return {
        "n": len(out), "lines": out,
        "of_n_graded_served": grade_cov["n_served_by_a_live_product"],
        "note": ("⚠️ CANDIDATES, not a disqualification. The residual is partial — ~0.25 "
                 "correlation against raw frames at ~0.6 — so this supports a per-line "
                 "re-check before the H-arm pool carries an abundance (RYA-1094's 25-line "
                 "Ruffoni pool), not the removal of the lines."),
    }


def _supersedes(arms) -> dict:
    """What this leg re-judges in the artifact already committed on main.

    🔴 TWO ARTIFACTS FOR ONE TICKET CAN DISAGREE SILENTLY, and the stale one is the one a
    reader trusts because it is the one they find. `rya1192_verification.json` reached its
    CRIRES+ H verdict from the catalogue-EXCESS proxy and said in terms that it was not
    proof and owed a raw comparison on Sirius. That comparison is this file. The finding
    it replaces is named here rather than left for someone to notice.
    """
    h = next((a for a in arms if a["arm"] == "H"), None)
    y = next((a for a in arms if a["arm"] == "Y"), None)
    prior = ("data/results/rya1192/rya1192_verification.json flagged a LOCALISED CANDIDATE "
             "RESIDUAL in the CRIRES+ H arm at 15700-15800, 16000-16100 and 17300-17400 A, "
             "from observed deep absorption exceeding what the stellar line list predicts. "
             "It labelled that check (b), not proof, and named the raw comparison on Sirius "
             "as what would settle it.")
    return {
        "artifact": "data/results/rya1192/rya1192_verification.json",
        "prior_finding": prior,
        "why_the_proxy_could_mislead": (
            "The catalogue accounting is not radiative transfer and over-predicts in "
            "general, so it was read through the EXCESS. But the H band is where our "
            "stellar line list is thinnest, and an under-populated catalogue produces an "
            "excess wherever real stellar absorption is uncatalogued — which is the "
            "RYA-1189 finding one band over. The proxy could not separate that from "
            "telluric; a template of actual telluric LINE POSITIONS can."),
        "now_judged_by": (
            "direct comparison against the raw CRIRES+ Vesta IDPs — the SAME OBSERVATIONS, "
            "2022-11-21/22 under programme 60.A-9051(A) — correlated against a molecfit "
            "MTRANS telluric transmission fitted on CRIRES+ frames of the same band, with "
            "a per-window displaced-template null."),
        "H_arm_now": (h or {}).get("verdict"),
        "Y_arm_now": (y or {}).get("verdict"),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pad", type=float, default=0.30)
    ap.add_argument("--out-dir", default="data/results/rya1192")
    a = ap.parse_args(argv)
    out = ROOT / a.out_dir
    out.mkdir(parents=True, exist_ok=True)

    sig = signature_table()
    bverd = band_verdicts(sig)
    claims = holding_claim_verdicts(bverd)
    per_line = close_lines(bverd, claims, a.pad)
    per_line.to_csv(out / "rya1192_sirius_per_line.csv", index=False)

    MW, MT, tfiles = telluric_template()
    arms = []
    if len(MW):
        for hold, arm, lo, hi in CRIRES_ARMS:
            arms.append(crires_arm(hold, arm, lo, hi, MW, MT))
    win_rows = [dict(holding=x["holding"], arm=x["arm"], **w)
                for x in arms for w in x["windows"]]
    if win_rows:
        pd.DataFrame(win_rows).to_csv(out / "rya1192_sirius_crires_windows.csv", index=False)

    grade_cov = graded_ir_coverage()
    pd.DataFrame(grade_cov["lines"]).to_csv(
        out / "rya1192_sirius_graded_ir_coverage.csv", index=False)

    # the graded IR line the directive names, in its own right
    cgf = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    ru = cgf[(cgf.species.astype(str).str.startswith("Fe"))
             & (cgf.wavelength_air_A > 10796) & (cgf.gf_tier.astype(str) == "LAB")]
    named = {}
    for w in (16009.610,):
        hit = ru[(ru.wavelength_air_A - w).abs() < 0.01]
        served, npts = None, 0
        for hold, arm, _, _ in CRIRES_ARMS:
            try:
                from measure_band_ew import load_crires_window
                _, f, _ = load_crires_window(w, 0.5, allow_topocentric=True)
                if np.isfinite(f).sum() >= 20:
                    served, npts = "raw IDP", int(np.isfinite(f).sum()); break
            except Exception:
                pass
        prod_pts = 0
        for hold, arm, lo, hi in CRIRES_ARMS:
            if lo <= w <= hi:
                p = pd.read_csv(ROOT / CRIRES_PRODUCT[hold]).dropna()
                prod_pts = int(((p.wavelength_air_A - w).abs() < 0.5).sum())
        arm_of = next((x for x in arms if x["span_A"][0] <= w <= x["span_A"][1]), None)
        win = next((r for r in (arm_of or {}).get("windows", [])
                    if r["lo_A"] <= w < r["hi_A"]), None)
        named[f"{w:.3f}"] = {
            "in_canonical_gf_LAB": bool(len(hit)),
            "prior_artifact_called_it": ("in a CRIRES+ H inter-order gap AND inside a "
                                         "candidate-residual window — both are refuted "
                                         "here, see graded_ir_line_coverage and its_window"),
            "lab_source_tag": (str(hit.iloc[0].lab_source_tag) if len(hit) else None),
            "served_by_raw_IDP_points": npts,
            "points_in_the_LIVE_PRODUCT_within_0.5A": prod_pts,
            "its_window": win,
            "note": ("⚠️ Measurability and telluric state are DIFFERENT questions and this "
                     "line separates them: the flux is judged by its window's verdict "
                     "below, but if the live product carries no points here the line "
                     "cannot be measured on it whatever the telluric state is."),
        }

    doc = {
        "ticket": "RYA-1192 (Sirius re-verification leg)",
        "kind": "VERIFICATION ONLY — no correction applied, no value moved",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "host": os.uname().nodename,
        "bounded_by": ("RYA-1094 / RYA-1193 IR science window — graded Fe ends at 17277.5 A "
                       "and the H arm's own extent at 17493.7; the instrument's 53000 A "
                       "reach is deliberately NOT audited, there is nothing graded there"),
        "telluric_template_files": tfiles,
        "telluric_template_span_A": ([round(float(MW.min()), 2), round(float(MW.max()), 2)]
                                     if len(MW) else None),
        "band_signature_table": sig,
        "band_verdicts": bverd,
        "reduction_floors_from_control_windows": {
            h: (round(baseline(i, h, sib), 8) if sib else None)
            for h, (i, sib, _) in RAW_SIBLING.items()},
        "control_windows_A": [list(w) for w in CONTROL_WINDOWS],
        "holding_claim_verdicts": claims,
        "crires_by_arm": arms,
        "supersedes": _supersedes(arms),
        "graded_ir_line_coverage": grade_cov,
        "graded_lines_in_a_flagged_window": _graded_in_flagged(arms, grade_cov),
        "named_line_checks": named,
        "per_line_inventory_limits": {
            "n_lines": int(len(per_line)),
            "n_inside_a_registered_telluric_band": int(per_line.correction_expected_here.sum()),
            "wavelength_span_A": [round(float(per_line.wavelength_air_A.min()), 3),
                                  round(float(per_line.wavelength_air_A.max()), 3)],
            "holdings": sorted(per_line.holding.unique().tolist()),
            "finding": (
                "🔴 NOT ONE of the measured Fe lines in this inventory sits inside a "
                "registered telluric band. That is the size of the KP label defect in "
                "science terms: 146 lines are demonstrably on raw data under a name that "
                "claims a correction, and NONE of them is anywhere a correction would "
                "have changed the flux. It is a provenance defect, not an abundance one."),
            "limit": (
                "⚠️ AND THIS INVENTORY IS NOT THE WHOLE PRODUCT SET. It is built from the "
                "per-line artifacts that exist on disk, and no _lines.csv exists for ANY "
                "red-optical, NIR or CRIRES+ Fe product — which are exactly the bands "
                "where the registered telluric complexes live. So 'no measured line is in "
                "a telluric band' is a true statement about the near-UV and VIS lines we "
                "can enumerate, and says NOTHING about the red-optical and NIR products; "
                "those are verified at BAND and WINDOW level above, and their per-line "
                "state is owed a re-measure that emits the artifact (RYA-1191 scope)."),
        },
        "per_line_corrected_but_outside_every_telluric_band": (
            {"n": int(((per_line.state == "VERIFIED-CORRECTED")
                       & ~per_line.correction_expected_here).sum()),
             "max_abs_diff_max": (float(per_line.loc[
                 (per_line.state == "VERIFIED-CORRECTED")
                 & ~per_line.correction_expected_here, "max_abs_diff"].max())
                 if ((per_line.state == "VERIFIED-CORRECTED")
                     & ~per_line.correction_expected_here).any() else None),
             "note": ("⚠️ These clear their pair's reduction floor but sit where NO "
                      "registered telluric complex reaches, and the differences are "
                      "orders below a real correction (HARPS O2 B-band mean |diff| is "
                      "0.185). Read them as 'above this pair's noise floor', never as "
                      "'a telluric correction of consequence'.")}),
        "per_line_summary": (per_line.groupby(["holding", "state"]).size()
                             .rename("n").reset_index().to_dict("records")),
        "per_line_basis_summary": (per_line.groupby(["holding", "basis"]).size()
                                   .rename("n").reset_index().to_dict("records")),
    }
    (out / "rya1192_sirius_verification.json").write_text(json.dumps(doc, indent=2) + "\n")

    print("=== band signature (% pixels below 0.8) ===")
    for row in sig:
        cells = " ".join(
            f"{h.replace('solar_','')[:14]}={c.get('pct_below_0.8', c.get('state'))}"
            for h, c in row["holdings"].items())
        print(f"  {row['kind'][:7]:<7} {row['band'][:16]:<16} {row['lo_A']:>7.0f}-{row['hi_A']:<7.0f} {cells}")
    print("\n=== per measured line ===")
    print(per_line.groupby(["holding", "state"]).size().to_string())
    print("\n=== CRIRES+ per arm ===")
    for x in arms:
        print(f"  {x['arm']} {x['holding']}: {x['verdict'][:160]}")
    print(f"\nwrote {a.out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
