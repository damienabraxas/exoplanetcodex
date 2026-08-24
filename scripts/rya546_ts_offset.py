#!/usr/bin/env python3
"""
RYA-546 Part B — does the RYA-545 production-TS-synth-EW absolute offset (Ti ran ~+0.3 dex high vs
Asplund) RIDE ALONG on the vintage-audit elements? Run the SAME production path (RYA-285
`_bisect_synth_abundance`, iSpec+Turbospectrum, MARCS.GES) on a few clean weak solar lines of
Mn I / Cr I / Ni I and compare the median absolute A(X) to Asplund-2021. A flat per-element offset
cancels in [X/H] (star-sun through the same code); a varying one would bias it.

Runs on Sirius: ISPEC_DIR=/srv/codex/engines/ispec_src, venv312 (iSpec TS compiled RYA-545).
"""
from __future__ import annotations
import csv, os, sys
from pathlib import Path
REPO = str(Path(__file__).resolve().parents[1])
if REPO not in sys.path:
    sys.path.insert(0, REPO)
import numpy as np
import pipeline.abundances_derive as ad
from pipeline import line_match
import ispec

STAR = dict(teff=5772.0, logg=4.44, feh=0.0, vturb=1.0)
HALFWIN_A = 0.3
# element -> (ion 'I', Z, Asplund-2021 solar A(X))
ELEMS = {'Mn': ('I', 25, 5.42), 'Cr': ('I', 24, 5.62), 'Ni': ('I', 28, 6.20), 'Ti': ('I', 22, 4.97)}


def _nonK10_wls(Z, ion_code):
    """Wavelengths whose canonical gf is NOT the ungraded K10 synth-gf (isolate lab/literature gf,
    so the TS absolute offset isn't confounded by the K10 gf-scale problem — RYA-521/545)."""
    keep = []
    with open(os.path.join(REPO, 'data/linelists/canonical_gf.csv')) as f:
        for r in csv.DictReader(f):
            if r['key_z'] == str(Z) and r['ion'] == str(ion_code) and r['loggf_reference'] != 'K10':
                keep.append(float(r['wavelength_air_A']))
    return np.sort(np.asarray(keep, dtype=float))


def _pool(el, ion, Z, amax_ew=55):
    # 🔴 RYA-1033: `wl in lab` used to compare 2-dp ROUNDED wavelengths across two files.
    # canonical_gf and the EW pool store the same line to different precision, so lines that
    # ARE lab-gf were read as K10 and silently left out of the offset pool. Tolerance match.
    lab = _nonK10_wls(Z, 1 if ion == 'I' else 2)
    cand = []
    with open(os.path.join(REPO, 'data/measured/sol_ew_results_v1.csv')) as f:
        for r in csv.DictReader(f):
            if r['element'] == el and r['ion'] == ion and r['blend_flag'] == 'False':
                e, er = float(r['ew_mA']), float(r['ew_err_mA'])
                if 8 <= e <= amax_ew and er / max(e, 1e-6) < 0.4:
                    cand.append(float(r['wavelength_air_A']))
    if not cand:
        return []
    res = line_match.match(np.asarray(cand, dtype=float), lab)
    return sorted({cand[i] for i in range(len(cand)) if res.index[i] >= 0})


def main():
    os.makedirs('/tmp/ispec_codex_synth', exist_ok=True)
    atm = ad._load_atmosphere(STAR['teff'], STAR['logg'], STAR['feh'], STAR['vturb'], model_grid='MARCS.GES')
    ll, iso, ch = ad._load_synth_resources()
    sa = ispec.read_solar_abundances(ad._ISPEC_SOLAR_ABUND_FILE)
    print("=== RYA-546 Part B — production TS synth-EW absolute offset vs Asplund-2021 (MARCS.GES) ===")
    ewmap = {}
    with open(os.path.join(REPO, 'data/measured/sol_ew_results_v1.csv')) as f:
        for r in csv.DictReader(f):
            if r['blend_flag'] == 'False':
                # Raw wavelength: `_pool` returns values read from THIS file, so the key
                # is exact by construction and needs no rounding (RYA-1033).
                ewmap[(r['element'], r['ion'], float(r['wavelength_air_A']))] = float(r['ew_mA'])
    for el, (ion, Z, asp) in ELEMS.items():
        pool = _pool(el, ion, Z)[:12]     # a handful of clean weak lab-gf lines (non-K10)
        A = []
        for wl in pool:
            ew = ewmap.get((el, ion, wl))
            if ew is None:
                continue
            w = np.linspace((wl - HALFWIN_A) / 10, (wl + HALFWIN_A) / 10, int(2 * HALFWIN_A * 240))
            try:
                a, conv, _ = ad._bisect_synth_abundance(w, ew, atm, STAR['teff'], STAR['logg'],
                                                        STAR['feh'], STAR['vturb'], ll, iso, sa, el, Z)
            except Exception:
                continue
            if conv and np.isfinite(a) and asp - 1.2 <= a <= asp + 1.2:   # physical bracket
                A.append(a)
        A = np.array(A)
        if len(A):
            med = float(np.median(A)); sem = float(np.std(A, ddof=1)/np.sqrt(len(A))) if len(A) > 1 else float('nan')
            print(f"  {el} {ion}: n={len(A)}/{len(pool)}  A_TS median={med:.3f}  vs Asplund {asp:.2f}  "
                  f"OFFSET={med-asp:+.3f}  (SEM {sem:.3f})")
        else:
            print(f"  {el} {ion}: no lines resolved")


if __name__ == '__main__':
    main()
