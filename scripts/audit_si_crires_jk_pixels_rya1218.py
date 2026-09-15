"""Validate exact Si NIR line windows against corrected CRIRES+ J products.

This is a coverage/conditioning audit only. It deliberately does not measure EWs
or abundances. A line is REACHED only when the complete +/- window is present on
one detector with finite corrected flux and a non-unity MTRANS transmission.
"""
from __future__ import annotations
import argparse, csv, glob, json
from pathlib import Path
import numpy as np
from astropy.io import fits

LINES = [
    (11991.57, 4.92, -0.109, "Si I 11991.57"),
    (11984.20, 4.93, 0.239, "Si I 11984.20"),
    (12103.54, 4.93, -0.351, "Si I 12103.54"),
    (12031.50, 4.95, 0.477, "Si I 12031.50"),
]

def audit(path: Path, half_window: float) -> list[dict]:
    with fits.open(path, memmap=False) as hdul:
        d = hdul["SPECTRUM"].data
        wave = np.asarray(d["WAVE"], float)
        flux = np.asarray(d["FLUX"], float)
        trans = np.asarray(d["MTRANS"], float)
        det = np.asarray(d["DETEC"], int)
        applied = bool(hdul[0].header.get("TELLAPP", False)) and "MTRANS" in [h.name for h in hdul]
    rows=[]
    for wavelength, ep, loggf, label in LINES:
        lo, hi = wavelength-half_window, wavelength+half_window
        m = np.isfinite(wave) & (wave >= lo) & (wave <= hi)
        n = int(m.sum())
        dets = sorted(set(det[m].tolist()))
        finite_flux = bool(np.isfinite(flux[m]).all()) if n else False
        nonunity = bool(n and np.any(np.abs(trans[m]-1.0) > 1e-8))
        step = float(np.median(np.diff(np.sort(wave[m])))) if n > 2 else 0.0
        edge_tol = max(0.08, 1.5 * step)
        reached = bool(n and finite_flux and applied and nonunity and len(dets)==1 and wave[m].min() <= lo + edge_tol and wave[m].max() >= hi - edge_tol)
        if not applied: status, reason = "HOLD_TELLURIC_EVIDENCE", "corrected product lacks TELLAPP/MTRANS"
        elif n == 0: status, reason = "HOLD_NO_PIXELS", "no valid pixels in +/-1 A window"
        elif len(dets) != 1: status, reason = "HOLD_DETECTOR_SPLIT", "window spans detector segments"
        elif wave[m].min() > lo + edge_tol or wave[m].max() < hi - edge_tol: status, reason = "HOLD_TRUNCATED_WINDOW", "available pixels do not cover complete +/-1 A window"
        elif not finite_flux: status, reason = "HOLD_NONFINITE_FLUX", "non-finite corrected flux in window"
        elif not nonunity: status, reason = "HOLD_UNITY_MTRANS", "MTRANS is unity in window"
        else: status, reason = "REACHED_EXACT_WINDOW", "complete corrected pixel window on one detector"
        rows.append(dict(product=path.as_posix(), setting=path.name.split("_") [2] if "_" in path.name else "", line=label, wavelength_air_A=wavelength, ep_eV=ep, loggf=loggf, window_lo_A=lo, window_hi_A=hi, n_pixels=n, detectors=";".join(map(str,dets)), observed_lo_A=float(wave[m].min()) if n else "", observed_hi_A=float(wave[m].max()) if n else "", telluric_applied=applied, nonunity_mtrans=nonunity, status=status, reason=reason))
    return rows

def main() -> None:
    ap=argparse.ArgumentParser(); ap.add_argument("--products", default="data/results/rya1219_crires_products/J/*.fits"); ap.add_argument("--out", default="data/audit/rya1218_si_protocol/si_jk_exact_pixel_validation.csv"); ap.add_argument("--summary", default="data/audit/rya1218_si_protocol/si_jk_exact_pixel_validation.json"); ap.add_argument("--half-window", type=float, default=1.0); args=ap.parse_args()
    rows=[]
    for name in sorted(glob.glob(args.products)): rows.extend(audit(Path(name), args.half_window))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    summary={"ticket":"RYA-1218","scope":"corrected CRIRES+ J products","window_half_width_A":args.half_window,"n_products":len(set(r['product'] for r in rows)),"n_lines":len(LINES),"counts":{},"measurement_status":"HOLD_NO_ABUNDANCE"}
    for r in rows: summary["counts"][r["status"]]=summary["counts"].get(r["status"],0)+1
    Path(args.summary).write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))
if __name__ == "__main__": main()
