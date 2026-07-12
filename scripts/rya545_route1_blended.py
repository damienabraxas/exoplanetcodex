#!/usr/bin/env python3
"""
RYA-545 Route 1 — Ti I/II ionization balance via PySME FULL-WINDOW BLENDED synthesis (MARCS).

Blend-aware (not the single-line COG stand-in): each target line's EW is measured as the DIFFERENTIAL
absorption it adds on top of everything else in the window —
    EW_target(A) = INT_window [ flux(window WITHOUT target line) - flux(window WITH target at A) ] dl
so the blend that sits under/next to the line is properly accounted (it depresses the local
continuum the target rides on). All blends come from linelist_full (atomic species, proper per-line
log-gamma vdW = damping_vdW); molecules are skipped (noted). Invert measured EW -> A(Ti) per line.

Pre-declared (RYA-545): Ti I pool = Lawler-2013/lab-gf lines only, ALL EW (widened by quality).
Ti II = untouched validator. Firewall: fix the instrument, never widen/de-blend until it balances.
Runs on Sirius (venv_pysme, MARCS marcs2012). Route 2 (production Turbospectrum synth-EW) cross-checks.
"""
from __future__ import annotations
import csv
import os
import sys
from pathlib import Path

REPO = str(Path(__file__).resolve().parents[1])
if REPO not in sys.path:
    sys.path.insert(0, REPO)
from pipeline import _runtime as _rt   # RYA-514 force-fork + single-thread BLAS (before numpy)
import numpy as np

from config.constants import FE_IONISATION_GATE

STAR = dict(teff=5772, logg=4.44, feh=0.0, vmic=1.0)
A_TI = 4.97
DELTA_TI1, DELTA_TI2 = 0.0506, -0.002
HALFWIN = 1.5          # window half-width (A) for the blended synthesis
EW_INT_HW = 1.0        # integration half-width (A) for the differential EW

LAB_REFS = {'LGWSC', '2013ApJS..', '1982MNRAS.', '1983MNRAS.', '1986MNRAS.', 'NWL', 'MFW', 'BLNP', 'SK'}

# element symbol -> Z (atomic blends only; molecules skipped)
_Z = {'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,'Na':11,'Mg':12,'Al':13,
      'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,
      'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'Sr':38,'Y':39,'Zr':40,'Nb':41,
      'Mo':42,'Ba':56,'La':57,'Ce':58,'Pr':59,'Nd':60,'Sm':62,'Eu':63,'Gd':64,'Dy':66}
_ROMAN = {'I':1,'II':2,'III':3}


def _lab_wls():
    keep = set()
    with open(os.path.join(REPO, 'data/linelists/canonical_gf.csv')) as f:
        for r in csv.DictReader(f):
            if r['key_z'] == '22' and r['ion'] == '1' and (r['loggf_reference'] in LAB_REFS or r['nist_grade'].strip()):
                keep.add(round(float(r['wavelength_air_A']), 2))
    return keep


def _pool(ion, graded):
    """(wl, elo, loggf, ew) for Ti of `ion` from sol_ew + canonical gf; Ti I optionally lab-graded."""
    lab = _lab_wls() if graded else None
    ew = {}
    with open(os.path.join(REPO, 'data/measured/sol_ew_results_v1.csv')) as f:
        for r in csv.DictReader(f):
            if r['element'] == 'Ti' and r['ion'] == ion and r['blend_flag'] == 'False':
                e, er = float(r['ew_mA']), float(r['ew_err_mA'])
                if 5 <= e <= 130 and er / max(e, 1e-6) < 0.5:
                    ew[round(float(r['wavelength_air_A']), 2)] = e
    gf = {}
    with open(os.path.join(REPO, 'data/linelists/canonical_gf.csv')) as f:
        for r in csv.DictReader(f):
            if r['key_z'] == '22' and r['ion'] == ('1' if ion == 'I' else '2'):
                try: gf[round(float(r['wavelength_air_A']), 2)] = (float(r['excitation_potential_eV']), float(r['log_gf']))
                except (ValueError, KeyError): pass
    out = []
    for wl, e in sorted(ew.items()):
        if lab is not None and wl not in lab:
            continue
        if wl in gf:
            out.append((wl, gf[wl][0], gf[wl][1], e))
    return out


def _window_rows(wl0):
    """All ATOMIC lines in [wl0-HALFWIN, wl0+HALFWIN] from linelist_full, as PySME short rows
    (proper per-line vdW). Molecules (species not in _Z) are skipped. Returns (rows, target_index)."""
    import pandas as pd
    rows = []
    with open(os.path.join(REPO, 'data/linelists/linelist_full.csv')) as f:
        for r in csv.DictReader(f):
            el, ion = r['element'], r.get('ion', 'I')
            if el not in _Z:
                continue                      # molecule / unsupported -> skip
            try:
                w = float(r['wavelength_air_A'])
            except ValueError:
                continue
            if abs(w - wl0) > HALFWIN:
                continue
            ionz = _ROMAN.get(ion, 1) if ion in _ROMAN else (int(ion) if str(ion).isdigit() else 1)
            try:
                vdw = float(r['damping_vdW'])
            except (ValueError, KeyError):
                vdw = 0.0
            rows.append(dict(species=f'{el} {ionz}', wlcent=w, excit=float(r['excitation_potential_eV']),
                             gflog=float(r['log_gf']), gamrad=7.8, gamqst=0.0, gamvw=vdw,
                             atom_number=_Z[el], ionization=ionz, lande=0.0, depth=0.5,
                             reference='RYA-545', error=0.0))
    return rows


def _synth(rows, offset):
    import pandas as pd
    from pysme.sme import SME_Structure
    from pysme.abund import Abund
    from pysme.linelist.linelist import LineList
    from pysme.synthesize import synthesize_spectrum
    if not rows:
        return None, None
    sme = SME_Structure()
    sme.teff, sme.logg, sme.monh = STAR['teff'], STAR['logg'], STAR['feh']
    sme.vmic, sme.vmac, sme.vsini = STAR['vmic'], 0.0, 0.0
    ab = Abund.solar(); ab['Ti'] = A_TI + offset; sme.abund = ab
    df = pd.DataFrame(rows).sort_values('wlcent').reset_index(drop=True)
    sme.linelist = LineList(df, lineformat='short')
    wc = float(np.median([r['wlcent'] for r in rows]))
    lo, hi = min(r['wlcent'] for r in rows) - 1.0, max(r['wlcent'] for r in rows) + 1.0
    sme.wave = np.linspace(lo, hi, int((hi - lo) * 240))
    sme.atmo.source = 'marcs2012.sav'; sme.atmo.method = 'grid'; sme.atmo.geom = 'PP'
    sme = synthesize_spectrum(sme)
    return np.asarray(sme.wave[0]), np.asarray(sme.synth[0])


def _blended_abund(wl, ew_meas, brackets=(-0.6, -0.3, 0.0, 0.3, 0.6, 0.9)):
    """Differential blended EW inversion for the Ti line at `wl`."""
    rows = _window_rows(wl)
    ti_idx = [i for i, r in enumerate(rows) if r['atom_number'] == 22 and abs(r['wlcent'] - wl) < 0.03]
    if not ti_idx:
        return float('nan')
    rows_noT = [r for i, r in enumerate(rows) if i not in ti_idx]     # drop target Ti line(s)
    ews, As = [], []
    for off in brackets:
        w1, f1 = _synth(rows, off)          # window WITH target Ti at A
        w2, f2 = _synth(rows_noT, off)      # window WITHOUT target Ti
        if w1 is None or w2 is None:
            continue
        f2i = np.interp(w1, w2, f2)
        m = (w1 > wl - EW_INT_HW) & (w1 < wl + EW_INT_HW)
        ews.append(float(np.trapz(f2i[m] - f1[m], w1[m]) * 1000.0))   # differential = target EW
        As.append(A_TI + off)
    ews, As = np.array(ews), np.array(As)
    o = np.argsort(ews); ews, As = ews[o], As[o]
    if len(ews) < 2 or ew_meas < ews[0] or ew_meas > ews[-1]:
        return float('nan')
    return float(np.interp(ew_meas, ews, As))


def main():
    print("=== RYA-545 Route 1 — PySME full-window BLENDED synthesis, MARCS, reference-blind ===")
    out = {}
    for ion, graded in (('I', True), ('II', False)):
        pool = _pool(ion, graded)
        wa = []
        for wl, elo, lgf, ew in pool:
            try:
                a = _blended_abund(wl, ew)
            except Exception as e:
                print(f"    Ti {ion} {wl}: fail ({e})"); continue
            if np.isfinite(a):
                wa.append((wl, a))
            else:
                print(f"       {wl:9.2f}  railed/no-target — dropped")
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
