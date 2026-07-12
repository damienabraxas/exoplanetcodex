#!/usr/bin/env python3
"""
RYA-545 Route 2 — Ti I/II ionization balance via the PRODUCTION synthesis-EW bisection
(RYA-285, iSpec + Turbospectrum), atmosphere-matched to MARCS.GES. Independent cross-check of
Route 1 (PySME full-window blended synthesis): same blend-aware differential-EW idea, but the real
pipeline machinery — GES linelist (canonical single-source gf), iSpec `generate_spectrum(code=
'turbospectrum')`, `_bisect_synth_abundance` (ew_target = ew_floor + ew_obs, blend floor at a_lo).

Pre-declared (RYA-545): Ti I = Lawler-2013/lab-gf lines only, ALL EW; Ti II = untouched validator.
Firewall: fix the instrument, never widen/de-blend until it balances.
Runs on Sirius: ISPEC_DIR=/srv/codex/engines/ispec_src, venv312 (iSpec TS compiled RYA-545). NO MERGE.
"""
from __future__ import annotations
import csv
import os
import sys
from pathlib import Path

REPO = str(Path(__file__).resolve().parents[1])
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import numpy as np
import pipeline.abundances_derive as ad
from config.constants import FE_IONISATION_GATE

STAR = dict(teff=5772.0, logg=4.44, feh=0.0, vturb=1.0)
A_TI, ATOM_TI = 4.97, 22
DELTA_TI1, DELTA_TI2 = 0.0506, -0.002
HALFWIN_A = 0.3    # narrow per-line window (matches the pipeline's ~±0.2 A GES regions); the
                   # bisection uses ABSOLUTE window EW (ew_floor+ew_obs), so a wide window lets
                   # blends/neighbour-Ti dominate and runs A away (RYA-545 driver fix).
LAB_REFS = {'LGWSC', '2013ApJS..', '1982MNRAS.', '1983MNRAS.', '1986MNRAS.', 'NWL', 'MFW', 'BLNP', 'SK'}


def _lab_wls():
    keep = set()
    with open(os.path.join(REPO, 'data/linelists/canonical_gf.csv')) as f:
        for r in csv.DictReader(f):
            if r['key_z'] == '22' and r['ion'] == '1' and (r['loggf_reference'] in LAB_REFS or r['nist_grade'].strip()):
                keep.add(round(float(r['wavelength_air_A']), 2))
    return keep


def _pool(ion, graded):
    lab = _lab_wls() if graded else None
    out = []
    with open(os.path.join(REPO, 'data/measured/sol_ew_results_v1.csv')) as f:
        for r in csv.DictReader(f):
            if r['element'] == 'Ti' and r['ion'] == ion and r['blend_flag'] == 'False':
                e, er = float(r['ew_mA']), float(r['ew_err_mA'])
                wl = round(float(r['wavelength_air_A']), 2)
                if 5 <= e <= 130 and er / max(e, 1e-6) < 0.5 and (lab is None or wl in lab):
                    out.append((wl, e))
    return sorted(set(out))


def main():
    print("=== RYA-545 Route 2 — production Turbospectrum synth-EW bisection, MARCS.GES ===")
    os.makedirs('/tmp/ispec_codex_synth', exist_ok=True)
    atm = ad._load_atmosphere(STAR['teff'], STAR['logg'], STAR['feh'], STAR['vturb'],
                              model_grid='MARCS.GES')
    import ispec
    linelist, isotopes, chem = ad._load_synth_resources()
    solar_abund = ispec.read_solar_abundances(ad._ISPEC_SOLAR_ABUND_FILE)

    out = {}
    for ion, graded in (('I', True), ('II', False)):
        pool = _pool(ion, graded)
        wa = []
        for wl, ew in pool:
            wobs_nm = np.linspace((wl - HALFWIN_A) / 10.0, (wl + HALFWIN_A) / 10.0,
                                  int(2 * HALFWIN_A * 240))
            try:
                A, conv, nit = ad._bisect_synth_abundance(
                    wobs_nm, ew, atm, STAR['teff'], STAR['logg'], STAR['feh'], STAR['vturb'],
                    linelist, isotopes, solar_abund, 'Ti', ATOM_TI)
            except Exception as e:
                print(f"    Ti {ion} {wl}: fail ({e})"); continue
            # Same physical-plausibility bracket Route 1 uses (A(Ti)=4.97 ± ~0.9): the pipeline
            # bisection returns conv=True even when a saturated line runs A away to 7-8 (unlike
            # Route 1's rail), so drop those non-physical inversions for an apples-to-apples pool.
            if conv and np.isfinite(A) and 4.07 <= A <= 5.87:
                wa.append((wl, float(A)))
            else:
                why = "not converged" if not (conv and np.isfinite(A)) else f"A={A:.2f} out of physical bracket"
                print(f"       {wl:9.2f}  {why} — dropped")
        a = np.array([x[1] for x in wa])
        med = float(np.median(a)) if len(a) else float('nan')
        mad = float(np.median(np.abs(a - med))) if len(a) else float('nan')
        sem = float(np.std(a, ddof=1) / np.sqrt(len(a))) if len(a) > 1 else float('nan')
        out[ion] = (med, sem, len(a))
        print(f"  Ti {ion}: n={len(a)}/{len(pool)}  A_LTE={med:.3f} (MAD {mad:.3f}, SEM {sem:.3f})")
        for wl, ab in sorted(wa):
            print(f"       {wl:9.2f}  A={ab:.3f}")
    (a1, s1, n1), (a2, s2, n2) = out['I'], out['II']
    a1n, a2n = a1 + DELTA_TI1, a2 + DELTA_TI2
    sem_bal = float(np.hypot(s1, s2))
    print(f"\n  A(Ti I)_LTE={a1:.3f} (SEM {s1:.3f})  A(Ti II)_LTE={a2:.3f} (SEM {s2:.3f})  LTE bal={a1-a2:+.3f}")
    print(f"  A(Ti I)_NLTE={a1n:.3f}  A(Ti II)_NLTE={a2n:.3f}  NLTE bal={a1n-a2n:+.3f} ± {sem_bal:.3f}")
    precise = sem_bal < 0.5 * FE_IONISATION_GATE
    ok = abs(a1n - a2n) < FE_IONISATION_GATE
    v = ("INCONCLUSIVE (SEM not << gate)" if not precise else
         "ACCEPT — balances within gate" if ok else "STOP — physics-limited, no wire")
    print(f"  gate |.|<{FE_IONISATION_GATE}, SEM<<gate  =>  VERDICT: {v}")


if __name__ == '__main__':
    main()
