#!/usr/bin/env python3
"""RYA-794 / RYA-377 — inventory the 18 CRIRES+ Vesta IDPs and what they actually reach.

    python3 scripts/inventory_vesta_crires.py

WHY THIS EXISTS. RYA-794 recorded that our own Vesta IDPs "do not exist on Sirius" — true
at the time, and the reason its Step 1 was blocked. They were never lost from ESO, only
from us: all 18 re-pull anonymously from the ESO archive in 37 MB. This script writes the
inventory so the next person does not have to rediscover where they came from.

THE COVERAGE TEST IS AGAINST PIXELS, NOT HEADERS. RYA-377's warning is the point: CRIRES+
settings are narrow chunks with detector gaps, so "Y/J/H/K are tiled" does not mean a
given line landed on a detector. Every line here is tested against the real WAVE array
with QUAL==0 and finite non-zero flux, never against WAVELMIN/WAVELMAX.

⚠️ WAVE is in NANOMETRES in these IDPs (TUNIT 'nm'), which is the same trap RYA-794 hit
from the other direction with Elgueta's sp/ files. Converted explicitly, once, here.

⚠️ SPECSYS is TOPOCENT. These are NOT barycentric; anything measuring a wavelength off
them owes a barycentric correction first.
"""
from __future__ import annotations

import csv
import glob
from pathlib import Path

import numpy as np
from astropy.io import fits
# Standalone-script bootstrap (RYA-313): repo root on sys.path BEFORE importing
# config/pipeline, so this runs from any cwd. Derived from __file__, never cwd.
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from config.constants import codex_path, codex_root  # RYA-810 path register

ROOT = Path(__file__).resolve().parents[1]
STAGE = Path(str(codex_path('data.spectra_vesta_crires')))
VIZ = Path(str(codex_root('work') / 'rya789' / 'data' / 'reference' / 'elgueta2026_vizier'))
GD_ROB = 295                     # G dwarf robust flag; the Sun's type. See RYA-794.
TOL_A = 0.3
OUT = ROOT / "data" / "audit" / "vesta_crires_plus"


def load_idps():
    out = []
    for f in sorted(glob.glob(str(STAGE / "*.fits"))):
        if open(f, "rb").read(6) != b"SIMPLE":
            raise SystemExit(f"{f} is not FITS — an HTML error page can masquerade "
                             f"as a download (RYA-791)")
        with fits.open(f) as hd:
            h, d = hd[0].header, hd[1].data[0]
            cols = hd[1].columns.names
            w = np.asarray(d["WAVE"], float) * 10.0      # nm -> Angstrom
            fl = np.asarray(d["FLUX"], float)
            q = np.asarray(d["QUAL"]) if "QUAL" in cols else np.zeros(len(w))
            ok = np.isfinite(w) & np.isfinite(fl) & (fl != 0) & (q == 0)
            out.append(dict(
                file=Path(f).name, setting=h.get("HIERARCH ESO INS WLEN ID"),
                date_obs=h.get("DATE-OBS"), exptime=h.get("EXPTIME"),
                snr=h.get("SNR"), pipeline=h.get("HIERARCH ESO PRO REC1 PIPE ID"),
                prog_id=h.get("HIERARCH ESO OBS PROG ID"), specsys=h.get("SPECSYS"),
                instrume=h.get("INSTRUME"),
                lam_min=float(w[ok].min()), lam_max=float(w[ok].max()),
                n_valid=int(ok.sum()), n_total=int(len(w)), _w=np.sort(w[ok])))
    return out


def certified_lines():
    rows = []
    for b in ("y", "j", "h"):
        p = VIZ / f"atomic{b}.dat"
        if not p.exists():
            continue
        for l in p.read_text(errors="replace").splitlines():
            if len(l) != 805 or l[GD_ROB] != "Y":
                continue
            try:
                rows.append(dict(band=b.upper(), species=l[13:17].strip(),
                                 wave_A=float(l[0:12]), ep_eV=float(l[18:27]),
                                 loggf=float(l[28:37])))
            except ValueError:
                continue
    return rows


def main() -> None:
    idps = load_idps()
    allw = np.sort(np.concatenate([d["_w"] for d in idps]))
    print(f"{len(idps)} CRIRES+ Vesta IDPs at {STAGE}")
    print(f"  {len(allw)} valid pixels, {allw.min():.1f}-{allw.max():.1f} A")
    print(f"  settings : {sorted({d['setting'] for d in idps})}")
    print(f"  pipeline : {sorted({str(d['pipeline']) for d in idps})}  "
          f"(cr2res -> CRIRES+, not the 2012-13 CRIRES)")
    print(f"  dates    : {min(d['date_obs'] for d in idps)[:10]} .. "
          f"{max(d['date_obs'] for d in idps)[:10]}")
    print(f"  prog id  : {sorted({str(d['prog_id']) for d in idps})}")
    print(f"  specsys  : {sorted({str(d['specsys']) for d in idps})}  "
          f"-> barycentric correction still OWED")

    def reached(lam):
        i = np.searchsorted(allw, lam)
        return any(0 <= j < len(allw) and abs(allw[j] - lam) <= TOL_A
                   for j in (i - 1, i))

    lines = certified_lines()
    for r in lines:
        r["reached"] = reached(r["wave_A"])

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "vesta_crires_plus_idp_manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[k for k in idps[0] if not k.startswith("_")])
        w.writeheader()
        w.writerows([{k: v for k, v in d.items() if not k.startswith("_")} for d in idps])
    with open(OUT / "vesta_crires_plus_certified_line_reach.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(lines[0]))
        w.writeheader()
        w.writerows(lines)

    import collections
    hit = collections.Counter(r["species"] for r in lines if r["reached"])
    tot = collections.Counter(r["species"] for r in lines)
    print(f"\nElgueta-certified SOLAR (G-dwarf) lines vs our actual pixels:")
    print(f"  {'species':8s} {'certified':>9s} {'reached':>8s}")
    for sp, n in tot.most_common():
        print(f"  {sp:8s} {n:9d} {hit[sp]:8d}")
    print(f"  {'TOTAL':8s} {sum(tot.values()):9d} {sum(hit.values()):8d}")
    print(f"\nwrote {OUT.relative_to(ROOT)}/ (2 files)")


if __name__ == "__main__":
    main()
