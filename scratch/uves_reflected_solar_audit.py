#!/usr/bin/env python3
"""
AUDIT (recon only): reflected-solar UVES reference sets.

Part A  Vesta UVES  -- RE-VERIFICATION of a re-downloaded set against the
                       RYA-370/372 cleared baseline ("don't trust the re-pull").
Part B  Ceres UVES  -- FRESH audit (new body, outside RYA-370 scope).

This is RECON ONLY. No normalization, EW, abundance, or velocity *conditioning*
(that is RYA-372's job downstream). The audit measures and reports; it does not
fix. Reuses the codex-data-audit UVES procedure + RYA-271 GES-trap rule.
"""

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

# --- Named recon boundaries (scratch-local; cited, not pipeline constants) ---
NUV_BLUE_A   = 3780.0   # below HARPS blue floor -> UVES346 NUV gain (the reason UVES is here)
OVERLAP_LO_A = 4800.0   # HARPS/UVES overlap low edge -> proxy validation window
OVERLAP_HI_A = 6800.0   # overlap high edge; also red-arm per-epoch telluric-watch threshold
SNR_FLOOR    = 200      # PIPELINE['snr_min_science']
C_KMS        = 299792.458
# Coarse frame-sanity lines (air A). First one present in a frame's coverage is used.
FRAME_LINES_A = [6562.8, 6173.34, 5250.21, 4383.55]  # Halpha, Fe I 6173, Fe I 5250, Fe I 4383
REFLEX_SANE_KMS = 60.0


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hdr_get(h, *keys, default=None):
    for k in keys:
        if k in h:
            return h[k]
    return default


def read_spectrum(hdul):
    """UVES Phase 3 IDP: BINTABLE in ext 1, single-row arrays WAVE/FLUX/ERR/SNR."""
    data = hdul[1].data
    names = {c.upper(): c for c in data.columns.names}
    wave = np.asarray(data[names["WAVE"]][0], dtype=float)
    flux = np.asarray(data[names["FLUX"]][0], dtype=float) if "FLUX" in names else None
    snr = None
    for cand in ("SNR", "SNR_FLUX"):
        if cand in names:
            snr = np.asarray(data[names[cand]][0], dtype=float)
            break
    cunit = str(hdul[1].header.get("TUNIT1", "")).lower()
    native_nm = ("nm" in cunit) or (np.nanmax(wave) < 2000.0)
    wave_A = wave * 10.0 if native_nm else wave
    return wave_A, flux, snr, native_nm


def frame_velocity_kms(wave_A, flux):
    """Coarse reflex-velocity sanity: parabolic min of the first in-range clean line."""
    if flux is None or wave_A.size < 16:
        return np.nan, None
    for lam0 in FRAME_LINES_A:
        if not (wave_A[0] + 2 < lam0 < wave_A[-1] - 2):
            continue
        win = (wave_A > lam0 - 1.5) & (wave_A < lam0 + 1.5)
        if win.sum() < 7:
            continue
        w, fl = wave_A[win], flux[win]
        if not np.all(np.isfinite(fl)):
            continue
        i = int(np.argmin(fl))
        if i == 0 or i == len(fl) - 1:
            continue
        x0, x1, x2 = w[i - 1], w[i], w[i + 1]
        y0, y1, y2 = fl[i - 1], fl[i], fl[i + 1]
        denom = (y0 - 2 * y1 + y2)
        lam_obs = x1 if denom == 0 else x1 - 0.5 * (x2 - x0) * (y2 - y0) / (2 * denom)
        return C_KMS * (lam_obs - lam0) / lam0, lam0
    return np.nan, None


def audit_file(path):
    rec = {"file": path.name, "md5": md5(path), "bytes": path.stat().st_size}
    try:
        with fits.open(path) as hdul:
            h = hdul[0].header
            rec["OBJECT"]   = str(hdr_get(h, "OBJECT", default="?"))
            rec["INSTRUME"] = str(hdr_get(h, "INSTRUME", default="?"))
            rec["DATE_OBS"] = str(hdr_get(h, "DATE-OBS", default="?"))
            rec["EXPTIME"]  = float(hdr_get(h, "EXPTIME", default=np.nan))
            rec["PRODCATG"] = str(hdr_get(h, "PRODCATG", "HIERARCH ESO PRO CATG", default="?"))
            rec["PROTYPE"]  = str(hdr_get(h, "PRO TYPE", "HIERARCH ESO PRO TYPE", default="?"))
            rec["SPECSYS"]  = str(hdr_get(h, "SPECSYS", default="?")).upper()
            rec["OBSTECH"]  = str(hdr_get(h, "OBSTECH", "HIERARCH ESO PRO TECH", default="?")).upper()
            rec["FLUXCAL"]  = str(hdr_get(h, "FLUXCAL", default="?"))
            rec["WLEN"]     = hdr_get(h, "HIERARCH ESO INS GRAT1 WLEN",
                                      "HIERARCH ESO INS GRAT2 WLEN", default="?")
            rec["BERV"]     = float(hdr_get(h, "HIERARCH ESO QC BERV",
                                            "HIERARCH ESO DRS BERV", default=np.nan))
            rec["TUNIT1"]   = str(hdul[1].header.get("TUNIT1", "?")) if len(hdul) > 1 else "?"
            wave_A, flux, snr, native_nm = read_spectrum(hdul)
            rec["WMIN_A"]  = float(np.nanmin(wave_A))
            rec["WMAX_A"]  = float(np.nanmax(wave_A))
            rec["SNR_MED"] = float(np.nanmedian(snr)) if snr is not None else np.nan
            rec["NATIVE_NM"] = bool(native_nm)
            v, lam = frame_velocity_kms(wave_A, flux)
            rec["VREFLEX_KMS"] = v
            rec["VLINE_A"] = lam
    except Exception as e:  # LOUD -- never silently skip
        rec["READ_ERROR"] = repr(e)
        print(f"[LOUD] failed to audit {path.name}: {e!r}", file=sys.stderr)
    return rec


def flags_for(rec):
    f = []
    if "READ_ERROR" in rec:
        return ["READ-ERROR (CRITICAL)"]
    if rec["SPECSYS"] == "HELIOCEN" or rec["OBSTECH"] == "MOS" or rec.get("NATIVE_NM"):
        f.append("GES-TRAP (CRITICAL: heliocen/MOS/nm)")
    if rec["PROTYPE"].upper() != "REDUCED" or "SCIENCE.SPECTRUM" not in rec["PRODCATG"].upper():
        f.append("PRODUCT-TYPE (not REDUCED SCIENCE.SPECTRUM)")
    snr = rec.get("SNR_MED", np.nan)
    if np.isfinite(snr) and snr < SNR_FLOOR:
        f.append(f"SNR<{SNR_FLOOR} ({snr:.0f})")
    v = rec.get("VREFLEX_KMS", np.nan)
    if np.isfinite(v) and abs(v) > REFLEX_SANE_KMS:
        f.append(f"FRAME-SUSPECT (v={v:.1f} km/s)")
    cov = []
    if rec.get("WMIN_A", np.inf) < NUV_BLUE_A:
        cov.append("NUV")
    if rec.get("WMIN_A", np.inf) < OVERLAP_HI_A and rec.get("WMAX_A", 0) > OVERLAP_LO_A:
        cov.append("OVERLAP")
    if rec.get("WMAX_A", 0) > OVERLAP_HI_A:
        cov.append("RED>6800(telluric-watch)")
    rec["COVERAGE"] = ",".join(cov) if cov else "none"
    return f


def audit_dir(d, expect_object):
    recs = []
    for p in sorted(Path(d).glob("*.fits")):
        r = audit_file(p)
        r["FLAGS"] = flags_for(r)
        if "OBJECT" in r and expect_object.upper() not in r["OBJECT"].upper():
            r["FLAGS"].append(f"WRONG-FOLDER? OBJECT={r['OBJECT']} in {expect_object} dir")
        recs.append(r)
    return recs


def md_table(recs):
    cols = ["file", "OBJECT", "DATE_OBS", "WLEN", "WMIN_A", "WMAX_A",
            "SNR_MED", "SPECSYS", "VREFLEX_KMS", "COVERAGE", "FLAGS"]
    head = "| " + " | ".join(cols) + " |\n| " + " | ".join("---" for _ in cols) + " |\n"
    rows = []
    for r in recs:
        cells = []
        for c in cols:
            v = r.get(c, "")
            if c == "FLAGS":
                v = "; ".join(v) if v else "ok"
            elif isinstance(v, float):
                v = "" if not np.isfinite(v) else (f"{v:.1f}" if c == "VREFLEX_KMS"
                     else f"{v:.0f}" if (c.startswith("W") or c.startswith("SNR")) else f"{v}")
            cells.append(str(v))
        rows.append("| " + " | ".join(cells) + " |")
    return head + "\n".join(rows)


def summary(name, recs):
    n = len(recs)
    crit = sum(1 for r in recs if any("CRITICAL" in x for x in r.get("FLAGS", [])))
    nuv = sum(1 for r in recs if "NUV" in r.get("COVERAGE", ""))
    ovl = sum(1 for r in recs if "OVERLAP" in r.get("COVERAGE", ""))
    lowsnr = sum(1 for r in recs if any(x.startswith("SNR<") for x in r.get("FLAGS", [])))
    foot = sum(r.get("bytes", 0) for r in recs) / 1e6
    dates = sorted(r["DATE_OBS"] for r in recs if r.get("DATE_OBS", "?") != "?")
    drange = f"{dates[0]} .. {dates[-1]}" if dates else "?"
    return (f"**{name}** -- N={n}, dates {drange}, {foot:.1f} MB; "
            f"NUV-capable={nuv}, in-overlap={ovl}, below-SNR-{SNR_FLOOR}={lowsnr}, "
            f"CRITICAL-flagged={crit}.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vesta-dir", required=True)
    ap.add_argument("--ceres-dir", required=True)
    ap.add_argument("--out", default="scratch/uves_reflected_solar_audit.md")
    a = ap.parse_args()

    vesta = audit_dir(a.vesta_dir, "VESTA")
    ceres = audit_dir(a.ceres_dir, "CERES")

    seen, dups = {}, []
    for r in vesta + ceres:
        m = r.get("md5")
        if m in seen:
            dups.append((seen[m], r["file"]))
        else:
            seen[m] = r["file"]

    lines = ["# UVES reflected-solar audit -- Vesta (re-verify) + Ceres (fresh)\n",
             "_RECON ONLY. Verdicts (GO / GO WITH CAVEATS / NO-GO) filled from flags below._\n",
             "## Summary", summary("Vesta UVES (re-verify vs RYA-370/372)", vesta),
             summary("Ceres UVES (fresh)", ceres),
             f"\nByte-duplicate file pairs: {len(dups)}"]
    if dups:
        lines += [f"- {x} == {y}" for x, y in dups]
    lines += ["\n## Part A -- Vesta UVES (re-verification)", md_table(vesta),
              "\n## Part B -- Ceres UVES (fresh audit)", md_table(ceres)]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text("\n".join(lines))
    print(f"Wrote {a.out}: Vesta N={len(vesta)}, Ceres N={len(ceres)}, byte-dups={len(dups)}")


if __name__ == "__main__":
    main()
