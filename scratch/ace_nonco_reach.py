#!/usr/bin/env python3
"""
SCOPING (recon only): non-CO molecular reach of the ACE-FTS solar atlas.

Extends the proven CO infrastructure (RYA-440/441/442/443) to the rest of C/N/O.
Reuses the RYA-392 verified atlas loader. RECON ONLY -- reports which C/N/O
molecular band systems (OH->O, NH->N, CH->C, CO 1-0 fundamental->C) fall in the
700-4430 cm^-1 ACE coverage, at what SNR and line density. No abundances: final
values inherit the CO-leg disk-center mu~1 intensity geometry + 1D->3D framework
(gated on RYA-444), not measured here.

Band windows are APPROXIMATE scan windows to localize the search; authoritative
band-heads + line positions come from HITRAN/ExoMol and the RYA-236 lists --
document the source per band. They are NOT science constants.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

# Single source: reuse the RYA-392 atlas loader -- do NOT duplicate the reader.
# Add the repo root to sys.path so the import resolves when run as `python scratch/...`.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pipeline.audits.audit_reference_atlases import read_any, sniff_xaxis, to_cm  # noqa: E402

# species, system, lo_cm, hi_cm  (approximate scan windows; VERIFY vs HITRAN/ExoMol)
BAND_WINDOWS = [
    ("C", "CO 1-0 fundamental",       1900.0, 2250.0),
    ("C", "CO 2-0 overtone (anchor)", 4150.0, 4360.0),
    ("O", "OH 1-0 fundamental",       2600.0, 3600.0),
    ("C", "CH 1-0 fundamental",       2650.0, 3100.0),
    ("N", "NH 1-0 fundamental",       3000.0, 3500.0),
]
# NOTE: ~3000-3500 cm^-1 carries OH+NH+CH simultaneously -> O/N/C are BLENDED there.


def band_report(xcm, y, lo, hi):
    m = (xcm >= lo) & (xcm <= hi)
    if int(m.sum()) < 10:
        return dict(npts=int(m.sum()), in_cov=False)
    yy = y[m]
    cont = np.nanpercentile(yy, 95)
    depth = float(1.0 - np.nanmin(yy) / cont) if cont > 0 else float("nan")
    hi_pts = yy[yy > np.nanpercentile(yy, 80)]
    noise = np.nanstd(hi_pts) if hi_pts.size else np.nan
    thr = cont - 3.0 * noise if np.isfinite(noise) else cont
    nlines = int(np.sum(yy < thr))
    snr = float(cont / noise) if (np.isfinite(noise) and noise > 0) else float("nan")
    return dict(npts=int(m.sum()), in_cov=True, depth=depth, nlines=nlines, snr=snr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True, help="ACE-FTS atlas file (RYA-390 store)")
    ap.add_argument("--out", default="scratch/ace_nonco_reach.md")
    a = ap.parse_args()

    # RYA-392 read_any returns (cols, header_text) with cols['x'], cols['y'] = the first
    # two numeric columns -- adapt to that signature (do NOT re-implement the reader).
    cols, hdr = read_any(a.atlas)
    x, y = cols["x"], cols["y"]
    kind, (lo, hi) = sniff_xaxis(x)
    if kind == "UNKNOWN":
        print(f"[LOUD] ACE atlas axis UNKNOWN (range {lo}-{hi}); abort", file=sys.stderr)
        sys.exit(2)
    xcm = to_cm(x, kind)
    cmin, cmax = float(np.nanmin(xcm)), float(np.nanmax(xcm))

    out = [f"# ACE-FTS non-CO molecular reach (recon only)\n",
           f"Atlas axis: {kind}; coverage {cmin:.1f}-{cmax:.1f} cm^-1 "
           f"({1e4/cmax:.2f}-{1e4/cmin:.2f} um).\n",
           "_Windows are approximate; confirm band-heads vs HITRAN/ExoMol. Abundances inherit "
           "the CO-leg disk-center mu~1 + 1D->3D framework (RYA-444); NOT measured here._\n",
           "| species | system | window cm^-1 | in_cov | npts | depth | n_lines | ~SNR |",
           "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for sp, name, blo, bhi in BAND_WINDOWS:
        r = band_report(xcm, y, blo, bhi)
        if not r.get("in_cov"):
            out.append(f"| {sp} | {name} | {blo:.0f}-{bhi:.0f} | NO | {r['npts']} |  |  |  |")
        else:
            out.append(f"| {sp} | {name} | {blo:.0f}-{bhi:.0f} | yes | {r['npts']} | "
                       f"{r['depth']:.2f} | {r['nlines']} | {r['snr']:.0f} |")
    out.append("\n**Blend flag:** OH + NH + CH all populate ~3000-3500 cm^-1 -> O/N/C "
               "disentangling there is a simultaneous-synthesis problem (Turbospectrum), never isolated EW.")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text("\n".join(out))
    print(f"Wrote {a.out}: coverage {cmin:.0f}-{cmax:.0f} cm^-1, {len(BAND_WINDOWS)} windows scanned")


if __name__ == "__main__":
    main()
