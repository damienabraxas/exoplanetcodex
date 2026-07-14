#!/usr/bin/env python3
"""
RYA-524 — LTE ionization-balance instrument for the master 27x2 audit (Cr gate + others).

Ports the RYA-544/545 Ti balance machinery, element-general and LTE-only. For the folded-in
RYA-546 Cr drop-to-LTE gate: measure A(Cr I)_LTE vs A(Cr II)_LTE on the SAME MARCS solar model
(Cr II = the LTE-robust majority-ion validator, UNTOUCHED). Classify:

  - |balance| < gate AND SEM << gate  -> LTE balances; the SH=0 scaled-Drawin +0.05..0.10 NLTE
    correction would BREAK balance by over-correcting -> RETIRE-TO-LTE (GO).
  - Cr I sits ~the SH=0 delta BELOW Cr II (NLTE would fix)                  -> STOP-needs-NLTE.
  - SEM >> gate (thin pools)                                               -> STOP-precision-limited.

FIREWALL: read + classify only. No gf edits, no gold changes, no wiring. Cr II untouched.
Weak/moderate lines only (linear COG, well-conditioned EW->A); railed inversions dropped, not clamped.
Runs on Sirius (PySME/MARCS). venv_pysme + ISPEC_DIR/PYTHONPATH.

  python scripts/rya524_ion_balance.py --element Cr --a-sun 5.62 --z 24
"""
from __future__ import annotations
import argparse, csv, os, sys
from pathlib import Path
REPO = str(Path(__file__).resolve().parents[1])
if REPO not in sys.path:
    sys.path.insert(0, REPO)
import numpy as np
from config.constants import FE_IONISATION_GATE

STAR = {'teff': 5772, 'logg': 4.44, 'feh': 0.0, 'vmic': 1.0}


def _load_pool(element, z, ion, ew_min, ew_max, relerr_max=0.5):
    """Join OUR measured solar EWs with canonical gf for one ion. Weak/moderate cut (linear COG)."""
    ew = {}
    with open(os.path.join(REPO, 'data/measured/sol_ew_results_v1.csv')) as f:
        for r in csv.DictReader(f):
            if r['element'] == element and r['ion'] == ion and r['blend_flag'] == 'False':
                e, er = float(r['ew_mA']), float(r['ew_err_mA'])
                if ew_min <= e <= ew_max and er / max(e, 1e-6) < relerr_max:
                    ew[round(float(r['wavelength_air_A']), 2)] = e
    ion_code = '1' if ion == 'I' else '2'
    gf = {}
    with open(os.path.join(REPO, 'data/linelists/canonical_gf.csv')) as f:
        for r in csv.reader(f):
            if len(r) > 7 and r[1] == str(z) and r[2] == ion_code:
                try:
                    gf[round(float(r[4]), 2)] = (float(r[5]), float(r[7]))
                except (ValueError, IndexError):
                    pass
    pool = []
    for wl, e in sorted(ew.items()):
        if wl in gf:
            elo, lgf = gf[wl]
            pool.append((wl, elo, lgf, e))
    return pool


def _lte_abund(element, z, ion_code, wl, elo, loggf, ew_meas, a_sun,
               brackets=(-0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2)):
    """LTE abundance reproducing the measured EW for one line, MARCS solar, multi-point COG.
    NaN (not clamped) if the EW is outside the synth COG range (railed = visible)."""
    import pandas as pd
    from pysme.sme import SME_Structure
    from pysme.abund import Abund
    from pysme.linelist.linelist import LineList
    from pysme.synthesize import synthesize_spectrum
    ews, As = [], []
    for off in brackets:
        sme = SME_Structure()
        sme.teff, sme.logg, sme.monh = STAR['teff'], STAR['logg'], STAR['feh']
        sme.vmic, sme.vmac, sme.vsini = STAR['vmic'], 0.0, 0.0
        ab = Abund.solar(); ab[element] = a_sun + off; sme.abund = ab
        row = dict(species=f'{element} {ion_code}', wlcent=wl, excit=elo, gflog=loggf,
                   gamrad=7.8, gamqst=0.0, gamvw=0.0, atom_number=z, ionization=ion_code,
                   lande=0.0, depth=0.6, reference='RYA-524', error=0.0)
        sme.linelist = LineList(pd.DataFrame([row]), lineformat='short')
        sme.wave = np.linspace(wl - 2.0, wl + 2.0, int(4 * 220))
        sme.atmo.source = 'marcs2012.sav'; sme.atmo.method = 'grid'; sme.atmo.geom = 'PP'
        sme = synthesize_spectrum(sme)
        w = np.asarray(sme.wave[0]); fl = np.asarray(sme.synth[0])
        m = (w > wl - 0.8) & (w < wl + 0.8)
        ews.append(float(np.trapz(1 - fl[m], w[m]) * 1000.0)); As.append(a_sun + off)
    ews, As = np.array(ews), np.array(As)
    order = np.argsort(ews); ews, As = ews[order], As[order]
    if ew_meas < ews[0] or ew_meas > ews[-1]:
        return float('nan')
    return float(np.interp(ew_meas, ews, As))


def balance(element, z, a_sun, ew1=(5.0, 60.0), ew2=(5.0, 95.0)):
    print(f"=== RYA-524 LTE ionization balance — {element} I vs {element} II — MARCS solar, reference-blind ===")
    print(f"  A({element})_sun = {a_sun}  gate = FE_IONISATION_GATE = {FE_IONISATION_GATE} dex\n")
    out = {}
    for ion, code, ewc in (('I', 1, ew1), ('II', 2, ew2)):
        pool = _load_pool(element, z, ion, ewc[0], ewc[1])
        vals = []
        for wl, elo, lgf, ew in pool:
            try:
                a = _lte_abund(element, z, code, wl, elo, lgf, ew, a_sun)
            except Exception as e:
                print(f"    {element} {ion} {wl}: inversion failed ({e})"); continue
            if not np.isfinite(a):
                print(f"    {element} {ion} {wl:9.2f}  railed (outside COG) — dropped"); continue
            vals.append((wl, a))
        A = np.array([a for _, a in vals])
        med = float(np.median(A)) if len(A) else float('nan')
        mad = float(np.median(np.abs(A - med))) if len(A) else float('nan')
        sem = float(np.std(A, ddof=1) / np.sqrt(len(A))) if len(A) > 1 else float('nan')
        out[ion] = (med, mad, sem, len(A))
        print(f"  {element} {ion}: n={len(A)}/{len(pool)}  A_LTE = {med:.3f}  MAD {mad:.3f}  SEM {sem:.3f}")
        for wl, a in vals:
            print(f"       {wl:9.2f}  A_LTE = {a:.3f}")
    a1, mad1, sem1, n1 = out['I']; a2, mad2, sem2, n2 = out['II']
    bal = a1 - a2
    print(f"\n  A({element} I)_LTE = {a1:.3f} (n{n1}, SEM {sem1:.3f})   "
          f"A({element} II)_LTE = {a2:.3f} (n{n2}, SEM {sem2:.3f})")
    print(f"  LTE ionization balance A({element} I) - A({element} II) = {bal:+.3f}  (gate |.| < {FE_IONISATION_GATE})")
    sem_max = np.nanmax([sem1, sem2])
    print(f"\n  --- classification ---")
    if not np.isfinite(sem_max) or sem_max >= FE_IONISATION_GATE:
        print(f"  STOP-PRECISION-LIMITED: max SEM {sem_max:.3f} NOT << {FE_IONISATION_GATE} "
              f"(thin pool: {element} I n{n1} / {element} II n{n2}). Cannot certify a retire on an "
              f"uncertain balance (RYA-546 Cr gate discipline). Flag for Ryan.")
    elif abs(bal) < FE_IONISATION_GATE:
        print(f"  RETIRE-TO-LTE (GO): LTE balances (|{bal:+.3f}| < {FE_IONISATION_GATE}); the SH=0 "
              f"scaled-Drawin +0.05..0.10 NLTE correction would BREAK balance by over-correcting.")
    else:
        print(f"  STOP-NEEDS-NLTE / OTHER: |{bal:+.3f}| > {FE_IONISATION_GATE}. If {element} I sits "
              f"~the SH=0 delta below {element} II, NLTE is genuinely needed (don't drop); else a "
              f"gf-scale/pool issue. Report, do not force.")
    return {'a1': a1, 'a2': a2, 'balance': bal, 'sem1': sem1, 'sem2': sem2, 'n1': n1, 'n2': n2}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--element', default='Cr')
    ap.add_argument('--z', type=int, default=24)
    ap.add_argument('--a-sun', type=float, default=5.62)
    a = ap.parse_args()
    balance(a.element, a.z, a.a_sun)
