#!/usr/bin/env python3
"""RYA-1232 -- write the full-arm molecfit-corrected CRIRES+ H frames in RYA-1219's product
format, so `rya1214_condition_crires_jk.py --arm H` rest-frame conditions them exactly as
it did J and K.

    python3 scripts/rya1232_crires_h_products.py --work <rya1191 --full-arm work dir> \
        --out data/results/rya1232_crires_products --manifest data/audit/rya1232_crires_h/corrected_products_manifest.csv

WHY. The H holding every product used, `solar_crires_plus_h_rya1094`, is Elgueta+2026's own
reduced spectrum (sp/Sun_H_rv.dat). Its telluric removal left the CH4 2nu3 Q-branch at
16656 A 36.5% deep (raw Vesta frames: 56%), which drove the OH oxygen fit to 9.47; and our
own RYA-1191 correction covered only segments holding a graded Fe line (16645-16810 A
never corrected). `rya1191_crires_h_correct.py --full-arm` now corrects all 116 segments
of the 7 H frames; this script assembles them.

FORMAT (RYA-1219, read by the conditioner): one FITS per frame, primary header with
TELLAPP=True, SPECSYS=TOPOCENT, MJD-OBS, WLEN, BAND; a SPECTRUM table with WAVE (vacuum A,
topocentric), FLUX (corrected, continuum-normalised PER SEGMENT with
`pipeline.crires_telluric.continuum_normalize` -- RYA-1219's "normalized per segment after
correction"), FLUX_RAW, ERR (raw units: the conditioner weights by (FLUX_RAW/ERR)^2),
MTRANS, ORDER, DETEC; and an MTRANS image. Saturated cores are NaN in FLUX already
(molecfit masks below its floor) and the conditioner drops MTRANS < 0.5.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline.crires_telluric import continuum_normalize  # noqa: E402


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def frame_tag(f) -> str:
    """The tag rya1191_crires_h_correct --full-arm wrote, without the segment suffix."""
    return f"{f.wlen_id}_{f.path.stem.split('T')[-1].replace(':', '').replace('.', '')}"


def main() -> int:
    from rya1191_crires_h_correct import h_frames
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path)
    a = ap.parse_args()
    (a.out / "H").mkdir(parents=True, exist_ok=True)
    a.manifest.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in h_frames():
        tag = frame_tag(f)
        cols = {k: [] for k in ("WAVE", "FLUX", "FLUX_RAW", "ERR", "MTRANS", "ORDER", "DETEC")}
        n_seg, missing = 0, []
        for s in f.segments:
            z = a.work / f"{tag}_o{int(s.order)}d{int(s.detector)}_corrected.npz"
            if not z.exists():
                w = np.asarray(s.wave_A, float)
                if np.isfinite(w).sum() >= 50:
                    missing.append(f"o{int(s.order)}d{int(s.detector)}")
                continue
            d = np.load(z)
            w = np.asarray(d["wave_A"], float)
            corr = np.asarray(d["flux_corr"], float)
            raw = np.asarray(d["flux_raw"], float)
            mt = np.asarray(d["mtrans"], float)
            sw, se = np.asarray(s.wave_A, float), np.asarray(s.err, float)
            ok = np.isfinite(sw) & np.isfinite(se)
            err = np.interp(w, sw[ok], se[ok], left=np.nan, right=np.nan)
            norm = continuum_normalize(w, corr)
            for k, v in (("WAVE", w), ("FLUX", norm), ("FLUX_RAW", raw), ("ERR", err),
                         ("MTRANS", mt), ("ORDER", np.full(w.size, int(s.order), float)),
                         ("DETEC", np.full(w.size, int(s.detector), float))):
                cols[k].append(v)
            n_seg += 1
        if missing:
            raise SystemExit(f"{f.path.name}: {len(missing)} segment(s) have no corrected "
                             f"file {missing} -- refusing to write a partially corrected frame")
        arr = {k: np.concatenate(v) for k, v in cols.items()}
        ph = fits.PrimaryHDU()
        h = ph.header
        h["RYA"] = "RYA-1232"
        h["ORIGIN"] = "exoplanetcodex"
        h["BASEFILE"] = f.path.name
        h["RAWOBJ"] = "Vesta"
        h["WLEN"] = f.wlen_id
        h["BAND"] = "H"
        h["MJD-OBS"] = float(f.mjd)
        h["SPECSYS"] = "TOPOCENT"
        h["TELLAPP"] = (True, "molecfit over EVERY H segment (rya1191 --full-arm)")
        h["TELLENG"] = "molecfit"
        h["MOLEC"] = "H2O,CO2,CH4"
        h["NSEGMENT"] = n_seg
        tab = fits.BinTableHDU.from_columns(
            [fits.Column(name=k, format="D", array=arr[k]) for k in
             ("WAVE", "FLUX", "FLUX_RAW", "ERR", "MTRANS", "ORDER", "DETEC")], name="SPECTRUM")
        img = fits.ImageHDU(data=arr["MTRANS"], name="MTRANS")
        out = a.out / "H" / f"vesta_crires_{f.wlen_id}_{f.path.stem}_telluric.fits"
        fits.HDUList([ph, tab, img]).writeto(out, overwrite=True)
        rows.append({"arm": "H", "product": str(out), "source": f.path.name, "setting": "",
                     "wlen": f.wlen_id, "telluric_applied": True, "mtrans": True,
                     "telluric_engine": "molecfit", "n_segments": n_seg,
                     "sha256": _sha(out)})
        print(f"  {f.path.name} {f.wlen_id}: {n_seg} segments -> {out.name}")
    pd.DataFrame(rows).to_csv(a.manifest, index=False)
    print(f"wrote {a.manifest} ({len(rows)} frames)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
