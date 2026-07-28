#!/usr/bin/env python3
"""
AUDIT (recon only): Solar NIRPS reflected-solar YJH reference set.

Discovery-mode: dumps the real NIRPS HDU/column/header structure on the first
file, then best-effort generic extraction with LOUD fallbacks. RECON ONLY --
no normalization, EW, abundance, telluric correction, or frame conditioning.
Headline check: telluric verification (permanent rule -- no IR abundance
without confirmed telluric correction; UNVERIFIED == CRITICAL).
Empirical measure-the-line frame check is DEFERRED to conditioning (project
line list), NOT hardcoded here. nm-native is EXPECTED for NIRPS (info, not a
GES-trap). The 2.3 um CO overtone is outside NIRPS range -- not assessed here.
Headers are PROVEN fallible (RYA-481): report distinct OBJECTs, flag MULTI-BODY.
"""

import argparse, hashlib, sys
from pathlib import Path
import numpy as np
from astropy.io import fits

# --- Named recon boundaries (scratch-local; cited, not pipeline constants) ---
SNR_FLOOR = 200            # PIPELINE['snr_min_science']
# Approx atmospheric band edges (air A) for TAGGING ONLY -- not science values.
Y_BAND   = (9800.0, 11100.0)
J_BAND   = (11500.0, 13500.0)
H_BAND   = (14500.0, 18000.0)
DEEP_H2O = [(13500.0, 14500.0), (18000.0, 19500.0)]  # exclude/telluric-verify
TELL_NAME_HINTS = ("TELL", "RECON", "TRANS", "ATMO", "CORR")


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def hdr_get(h, *keys, default=None):
    for k in keys:
        if k in h:
            return h[k]
    return default


def discover(path, out_lines):
    """Dump the real HDU/column/WCS/telluric-header map for the first file."""
    out_lines.append(f"### Structure discovery -- {path.name}\n")
    with fits.open(path) as hdul:
        for i, hdu in enumerate(hdul):
            kind = type(hdu).__name__
            out_lines.append(f"- HDU[{i}] name={hdu.name!r} type={kind}")
            if getattr(hdu, "columns", None) is not None:
                out_lines.append(f"    columns: {list(hdu.columns.names)}")
            hh = hdu.header
            wcs = {k: hh[k] for k in ("CRVAL1","CDELT1","CD1_1","CRPIX1","CUNIT1","CTYPE1") if k in hh}
            if wcs:
                out_lines.append(f"    wcs: {wcs}")
            tell = [c for c in hh.cards if any(t in str(c.keyword).upper() for t in TELL_NAME_HINTS)]
            for c in tell:
                out_lines.append(f"    tell-hdr: {c.keyword} = {c.value}")
    out_lines.append("")


def read_spectrum(hdul):
    """Return (wave_A, flux, snr, native_nm, flux_cols). Handles BINTABLE or image+WCS."""
    flux_cols = []
    # Try a BINTABLE with array columns first.
    for hdu in hdul[1:]:
        if getattr(hdu, "columns", None) is None:
            continue
        names = {c.upper(): c for c in hdu.columns.names}
        wkey = next((names[k] for k in ("WAVE","WAVELENGTH","LAMBDA") if k in names), None)
        if wkey is None:
            continue
        data = hdu.data
        wave = np.asarray(data[wkey][0], dtype=float)
        fkey = next((names[k] for k in ("FLUX_CORR","FLUX_TELL","FLUX_EL","FLUX") if k in names), None)
        flux = np.asarray(data[fkey][0], dtype=float) if fkey else None
        flux_cols = [names[k] for k in names if "FLUX" in k]
        skey = next((names[k] for k in ("SNR","SNR_FLUX") if k in names), None)
        snr = np.asarray(data[skey][0], dtype=float) if skey else None
        cunit = str(hdu.header.get("TUNIT" + str(list(hdu.columns.names).index(wkey)+1), "")).lower()
        native_nm = ("nm" in cunit) or (np.nanmax(wave) < 3000.0)
        return (wave*10.0 if native_nm else wave), flux, snr, native_nm, flux_cols
    # Fall back to a 1D image + WCS.
    for hdu in hdul:
        if getattr(hdu, "data", None) is not None and np.ndim(hdu.data) == 1:
            h = hdu.header
            n = hdu.data.size
            crval = float(hdr_get(h, "CRVAL1", default=np.nan))
            cdelt = float(hdr_get(h, "CDELT1", "CD1_1", default=np.nan))
            crpix = float(hdr_get(h, "CRPIX1", default=1.0))
            if not np.isfinite(crval) or not np.isfinite(cdelt):
                continue
            wave = crval + (np.arange(n) - (crpix - 1)) * cdelt
            cunit = str(h.get("CUNIT1", "")).lower()
            native_nm = ("nm" in cunit) or (np.nanmax(wave) < 3000.0)
            return (wave*10.0 if native_nm else wave), np.asarray(hdu.data, float), None, native_nm, ["IMAGE"]
    return None, None, None, False, []


def telluric_status(wave_A, flux, flux_cols):
    """CORRECTED / UNCORRECTED / UNKNOWN -- heuristic, reported LOUD on ambiguity."""
    name_evidence = any(any(t in c.upper() for t in TELL_NAME_HINTS) for c in flux_cols)
    if flux is None:
        return "UNKNOWN(no-flux)", np.nan
    # Deep-water recovery heuristic: corrected spectra recover flux in the bands.
    ratio = np.nan
    for lo, hi in DEEP_H2O:
        win = (wave_A > lo) & (wave_A < hi)
        cont = (wave_A > lo - 400) & (wave_A < lo - 50)
        if win.sum() >= 10 and cont.sum() >= 10:
            dm = np.nanmedian(flux[win]); cm = np.nanmedian(flux[cont])
            if np.isfinite(dm) and np.isfinite(cm) and cm > 0:
                ratio = dm / cm
                break
    if name_evidence:
        return "CORRECTED(named-col)", ratio
    if np.isfinite(ratio):
        return ("UNCORRECTED(band-zeroed)" if ratio < 0.05 else "MAYBE-CORRECTED(band-recovered)"), ratio
    return "UNKNOWN", ratio


def coverage(wave_A):
    tags = []
    wmin, wmax = np.nanmin(wave_A), np.nanmax(wave_A)
    for nm, (lo, hi) in (("Y",Y_BAND),("J",J_BAND),("H",H_BAND)):
        if wmin < hi and wmax > lo:
            tags.append(nm)
    if any(wmin < hi and wmax > lo for lo, hi in DEEP_H2O):
        tags.append("DEEP-H2O")
    return ",".join(tags) if tags else "none"


def audit_file(path, first):
    rec = {"file": path.name, "md5": md5(path), "bytes": path.stat().st_size}
    try:
        with fits.open(path) as hdul:
            h = hdul[0].header
            rec["OBJECT"]   = str(hdr_get(h, "OBJECT", default="?"))
            rec["INSTRUME"] = str(hdr_get(h, "INSTRUME", default="?"))
            rec["DATE_OBS"] = str(hdr_get(h, "DATE-OBS", default="?"))
            rec["EXPTIME"]  = float(hdr_get(h, "EXPTIME", default=np.nan))
            rec["AIRMASS"]  = float(hdr_get(h, "AIRMASS", "HIERARCH ESO TEL AIRM START", default=np.nan))
            rec["PRODCATG"] = str(hdr_get(h, "PRODCATG", "HIERARCH ESO PRO CATG", default="?"))
            rec["PROTYPE"]  = str(hdr_get(h, "PRO TYPE", "HIERARCH ESO PRO TYPE", default="?"))
            rec["SPECSYS"]  = str(hdr_get(h, "SPECSYS", default="?")).upper()
            rec["BERV"]     = float(hdr_get(h, "HIERARCH ESO QC BERV", "HIERARCH ESO DRS BERV", default=np.nan))
            rec["DRS_RV"]   = float(hdr_get(h, "HIERARCH ESO QC CCF RV", "HIERARCH ESO DRS CCF RVC", default=np.nan))
            w, flux, snr, native_nm, fcols = read_spectrum(hdul)
            if w is None:
                rec["READ_ERROR"] = "no recognizable wave array (table or WCS)"
                print(f"[LOUD] {path.name}: {rec['READ_ERROR']}", file=sys.stderr)
                return rec
            rec["WMIN_A"] = float(np.nanmin(w)); rec["WMAX_A"] = float(np.nanmax(w))
            rec["NATIVE_NM"] = bool(native_nm)
            rec["SNR_MED"] = float(np.nanmedian(snr)) if snr is not None else np.nan
            rec["FLUX_COLS"] = ",".join(fcols)
            rec["TELLURIC"], rec["TELL_RATIO"] = telluric_status(w, flux, fcols)
            rec["COVERAGE"] = coverage(w)
            if first:
                rec["_DISCOVER"] = path
    except Exception as e:
        rec["READ_ERROR"] = repr(e)
        print(f"[LOUD] failed to audit {path.name}: {e!r}", file=sys.stderr)
    return rec


def flags_for(rec, expect_object):
    f = []
    if "READ_ERROR" in rec:
        return ["READ-ERROR (CRITICAL)"]
    if "NIRPS" not in rec["INSTRUME"].upper():
        f.append(f"NOT-NIRPS (INSTRUME={rec['INSTRUME']})")
    if expect_object and expect_object.upper() not in rec["OBJECT"].upper():
        f.append(f"WRONG-OBJECT ({rec['OBJECT']})")
    t = rec.get("TELLURIC", "UNKNOWN")
    if not t.startswith("CORRECTED"):
        f.append(f"TELLURIC-UNVERIFIED (CRITICAL: {t})")
    snr = rec.get("SNR_MED", np.nan)
    if np.isfinite(snr) and snr < SNR_FLOOR:
        f.append(f"SNR<{SNR_FLOOR} ({snr:.0f})")
    if rec.get("NATIVE_NM"):
        f.append("nm-native (INFO: expected for NIRPS)")
    if "DEEP-H2O" in rec.get("COVERAGE", ""):
        f.append("DEEP-H2O (info: exclude/telluric-verify)")
    return f


def md_table(recs):
    cols = ["file","OBJECT","DATE_OBS","WMIN_A","WMAX_A","SNR_MED","SPECSYS",
            "TELLURIC","COVERAGE","FLAGS"]
    head = "| " + " | ".join(cols) + " |\n| " + " | ".join("---" for _ in cols) + " |"
    rows = []
    for r in recs:
        cells = []
        for c in cols:
            v = r.get(c, "")
            if c == "FLAGS":
                v = "; ".join(v) if v else "ok"
            elif isinstance(v, float):
                v = "" if not np.isfinite(v) else (f"{v:.0f}" if c.startswith(("W","SNR")) else f"{v}")
            cells.append(str(v))
        rows.append("| " + " | ".join(cells) + " |")
    return head + "\n" + "\n".join(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nirps-dir", required=True)
    ap.add_argument("--expect-object", default="", help="blank = report-only (recommended); detects MULTI-BODY")
    ap.add_argument("--out", default="scratch/nirps_solar_audit.md")
    a = ap.parse_args()

    files = sorted(Path(a.nirps_dir).glob("*.fits"))
    recs, disc, first = [], [], True
    for p in files:
        r = audit_file(p, first)
        if "_DISCOVER" in r:
            discover(r.pop("_DISCOVER"), disc); first = False
        r["FLAGS"] = flags_for(r, a.expect_object)
        recs.append(r)

    n = len(recs)
    crit = sum(1 for r in recs if any("CRITICAL" in x for x in r.get("FLAGS", [])))
    tell_ok = sum(1 for r in recs if r.get("TELLURIC","").startswith("CORRECTED"))
    lowsnr = sum(1 for r in recs if any(x.startswith("SNR<") for x in r.get("FLAGS", [])))
    foot = sum(r.get("bytes",0) for r in recs)/1e6
    dates = sorted(r["DATE_OBS"] for r in recs if r.get("DATE_OBS","?") != "?")
    drange = f"{dates[0]} .. {dates[-1]}" if dates else "?"

    # Header-fallibility check (RYA-481): distinct OBJECTs + per-body counts.
    bodies = {}
    for r in recs:
        o = r.get("OBJECT", "?")
        bodies[o] = bodies.get(o, 0) + 1
    body_str = ", ".join(f"{k}={v}" for k, v in sorted(bodies.items()))
    multi = "MULTI-BODY (FLAG: folder mixes bodies)" if len(bodies) > 1 else "single-body"

    lines = ["# Solar NIRPS audit (recon only) -- YJH reflected-solar\n",
             "_RECON ONLY. Verdict GO / GO WITH CAVEATS / NO-GO from flags. Telluric"
             " UNVERIFIED == CRITICAL (permanent rule). K-band CO is OUTSIDE NIRPS"
             " range and is NOT assessed here._\n",
             "## Summary",
             f"N={n}, dates {drange}, {foot:.1f} MB; telluric-CORRECTED={tell_ok}/{n},"
             f" below-SNR-{SNR_FLOOR}={lowsnr}, CRITICAL-flagged={crit}.",
             f"Bodies present ({multi}): {body_str}\n",
             "## Structure discovery", *disc,
             "## Per-frame audit", md_table(recs)]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text("\n".join(lines))
    print(f"Wrote {a.out}: N={n}, telluric-CORRECTED={tell_ok}, CRITICAL={crit}, bodies=[{body_str}]")


if __name__ == "__main__":
    main()
