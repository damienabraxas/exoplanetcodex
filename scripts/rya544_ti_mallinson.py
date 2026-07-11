#!/usr/bin/env python3
"""
RYA-544 — wire the Mallinson-2024 Ti departure grid + ionization-balance acceptance gate.

RYA-542's 3rd-reference adjudication established that the correct solar Ti I NLTE correction
is ~+0.04 (Mallinson 2022/2024 + Sitnova 2020, ab-initio Grumer-Barklem-2020 H collisions,
ionization-balance-validated), NOT our Engine-A MAFAGS-OS +0.108 or Engine-B MARCS +0.20 (both
on the outdated Bergemann-2011 scaled-Drawin atom `atom.ti503b`). This script wires the *grid*
(Mallinson 2024, Zenodo 10753497, PySME-compatible) and lets it PRODUCE the correction, then
gates on OUR ionization balance.

FIREWALL (RYA-544): wire the grid and let it produce balance. Do NOT tune to +0.04. The winner
is decided by whether OUR Ti I(new NLTE) balances OUR Ti II in the SAME MARCS atmosphere within
FE_IONISATION_GATE — reference-blind. Ti II is the VALIDATOR, not a third candidate.

Runs on Sirius only (grids Sirius-only, RYA-526); PySME + MARCS via venv_tsfitpy/venv_pysme.

  --derive : Part 1. Register Ti on the Mallinson grid + derive OUR per-line Ti I NLTE delta.
  --gate   : Part 2. LTE EW-invert OUR solar Ti I + Ti II pools on MARCS, apply the delta,
             test |A(Ti I)_NLTE - A(Ti II)| < FE_IONISATION_GATE. ACCEPT / STOP.
"""
from __future__ import annotations
import argparse
import csv
import os
import sys
from pathlib import Path

# repo root (this file: <repo>/scripts/rya544_ti_mallinson.py)
REPO = str(Path(__file__).resolve().parents[1])
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import numpy as np
from pipeline import pysme_nlte as H          # the PySME/MARCS NLTE-derivation harness
from config.constants import FE_IONISATION_GATE

# --- Mallinson-2024 solar node + reference abundance --------------------------
STAR = {'teff': 5772, 'logg': 4.44, 'feh': 0.0, 'vmic': 1.0}   # solar (matches harness default)
A_TI = 4.97                                                    # our solar A(Ti) (constants.py)
TIII_NLTE = -0.002                                             # Mallinson-2024 solar Ti II corr (~0)

# Ti I diagnostic lines for the DERIVATION — the RYA-542 gate lines (clean, well-measured,
# GES-identified). (wl_air, loggf, Elow_eV, gamvw_ABO from GES). Level labels/J are resolved
# from the grid BY ENERGY (auto_labels) — the Mallinson README warns VALD labels differ from
# the grid's NIST labels, so energy-matching is the robust path.
TI_DIAG = [
    (5689.460, -0.360, 2.297, 794.242),
    (5648.565, -0.161, 2.495, 926.229),
    (5662.150,  0.010, 2.318, 808.240),
]
_HC = 12398.419843   # eV*Angstrom (for Eup = Elow + hc/lambda)


def _register_ti():
    """Make Ti a first-class PySME element on the Mallinson grid (idempotent monkeypatch;
    the permanent harness edit mirrors this)."""
    H._A_SUN.setdefault('Ti', A_TI)
    H._GRID_FILENAME.setdefault('Ti', 'nlte_Ti_pysme.grd')
    _orig_Z = H._Z
    H._Z = lambda el, _o=_orig_Z: 22 if el == 'Ti' else _o(el)


# Mallinson's authoritative NIST level table (ships with the grid). Our Amarsi-format .grd
# reader mis-parses the Mallinson grid's internal level array, so we resolve level labels from
# this ASCII table (energy-matched) — the grid was "mapped to labels given in NIST" (README),
# and PySME's native set_nlte reads the .grd's departures directly for the synthesis.
LABEL_FILE = '/srv/codex/grids/nlte/amarsi_galah/label_Ti.txt'


def _label_index(ion_code: int = 1):
    """[(energy_eV, 'config term', J)] for Ti of the given ion from label_Ti.txt."""
    rows = []
    with open(LABEL_FILE) as f:
        for ln in f:
            if ln.startswith('#'):
                continue
            t = ln.split()
            if len(t) < 7 or t[1] != 'Ti' or t[2] != str(ion_code):
                continue
            rows.append((float(t[6]), f"{t[3]} {t[4]}", float(t[5])))
    return rows


def _resolve(idx, e_eV, tol=0.03):
    best = min(idx, key=lambda r: abs(r[0] - e_eV))
    if abs(best[0] - e_eV) > tol:
        raise ValueError(f"no Ti label within {tol} eV of {e_eV:.3f} (nearest {best[0]:.3f})")
    return best[1], best[2]     # 'config term', J


def _build_ti_lines():
    """Resolve NLTE_LINES['Ti'] from label_Ti.txt by energy (grid-native NIST labels)."""
    idx = _label_index(1)
    lines = []
    for wl, gf, elo, vw in TI_DIAG:
        eup = elo + _HC / wl
        tl, jl = _resolve(idx, elo)
        tu, ju = _resolve(idx, eup)
        lines.append((wl, gf, elo, jl, eup, ju, tl, tu, vw))
    H.NLTE_LINES['Ti'] = lines
    return lines


def derive():
    _register_ti()
    lines = _build_ti_lines()
    print("=== RYA-544 Part 1 — OUR Ti I NLTE correction from the Mallinson-2024 grid ===")
    print(f"  grid = {H._GRID_FILENAME['Ti']}  atmosphere = MARCS (marcs2012, PP)  node = {STAR}")
    for ln in lines:
        print(f"  line {ln[0]:.3f}: Elow={ln[2]:.3f} Eup={ln[4]:.3f}  levels [{ln[6]} -> {ln[7]}] J {ln[3]}->{ln[5]}")
    res = H.nlte_delta('Ti', star=STAR)
    print("\n  per-line delta = A_NLTE - A_LTE:")
    for wl, d in res['per_line'].items():
        print(f"    Ti I {wl:.3f}: {d:+.4f}")
    print(f"  MEDIAN Ti I NLTE delta (Mallinson-2024, ab-initio) = {res['delta_median']:+.4f}")
    print(f"  [context: Engine-A MAFAGS-OS +0.108, Engine-B MARCS +0.20, lit ~+0.04 — do NOT tune]")
    return res


# ---------------------------------------------------------------------------
# Part 2 — ionization-balance gate: LTE EW-invert OUR Ti I + Ti II pools on MARCS.
# ---------------------------------------------------------------------------
def _load_pool(ion: str, ew_min=5.0, ew_max=60.0, relerr_max=0.5):
    """Join OUR measured solar Ti EWs with canonical gf. Quality cuts: unblended, EW in
    [ew_min,ew_max] mA (weak/moderate = LINEAR COG, well-conditioned inversion, vdW-insensitive;
    saturated >60 mA lines have a flat COG so EW->A is ill-posed), relative EW error < relerr_max.
    Returns [(wl, elow_eV, loggf, ew_mA), ...]."""
    ew = {}
    with open(os.path.join(REPO, 'data/measured/sol_ew_results_v1.csv')) as f:
        for r in csv.DictReader(f):
            if r['element'] == 'Ti' and r['ion'] == ion and r['blend_flag'] == 'False':
                e, er = float(r['ew_mA']), float(r['ew_err_mA'])
                if ew_min <= e <= ew_max and er / max(e, 1e-6) < relerr_max:
                    ew[round(float(r['wavelength_air_A']), 2)] = e
    ion_code = '1' if ion == 'I' else '2'
    gf = {}
    with open(os.path.join(REPO, 'data/linelists/canonical_gf.csv')) as f:
        for r in csv.reader(f):
            if len(r) > 7 and r[1] == '22' and r[2] == ion_code:
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


def _lte_abund(ion_code: int, wl: float, elo: float, loggf: float, ew_meas: float,
               brackets=(-0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2)) -> float:
    """LTE abundance that reproduces the measured EW for one Ti line, on the MARCS solar
    model (atmosphere-matched to the grid). Multi-point COG interpolation. PySME, no NLTE.
    Returns NaN (not a silently-clamped endpoint) if the measured EW is outside the synth
    COG range — so a railed inversion is visible, not mistaken for a real abundance."""
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
        ab = Abund.solar(); ab['Ti'] = A_TI + off; sme.abund = ab
        row = dict(species=f'Ti {ion_code}', wlcent=wl, excit=elo, gflog=loggf,
                   gamrad=7.8, gamqst=0.0, gamvw=0.0, atom_number=22, ionization=ion_code,
                   lande=0.0, depth=0.6, reference='RYA-544', error=0.0)
        sme.linelist = LineList(pd.DataFrame([row]), lineformat='short')
        sme.wave = np.linspace(wl - 2.0, wl + 2.0, int(4 * 220))
        sme.atmo.source = 'marcs2012.sav'; sme.atmo.method = 'grid'; sme.atmo.geom = 'PP'
        sme = synthesize_spectrum(sme)
        w = np.asarray(sme.wave[0]); fl = np.asarray(sme.synth[0])
        m = (w > wl - 0.8) & (w < wl + 0.8)
        ews.append(float(np.trapz(1 - fl[m], w[m]) * 1000.0))
        As.append(A_TI + off)
    ews, As = np.array(ews), np.array(As)
    order = np.argsort(ews)
    ews, As = ews[order], As[order]
    if ew_meas < ews[0] or ew_meas > ews[-1]:
        return float('nan')                       # railed — outside COG range, do NOT clamp
    return float(np.interp(ew_meas, ews, As))


def gate(delta_ti1: float):
    print("\n=== RYA-544 Part 2 — ionization-balance acceptance gate (reference-blind) ===")
    print(f"  OUR Ti I(NLTE) = OUR Ti II on MARCS; gate = FE_IONISATION_GATE = {FE_IONISATION_GATE} dex")
    diag_wl = {ln[0] for ln in TI_DIAG}   # the 3 clean RYA-542 lines (good gf)
    out = {}
    for ion, code in (('I', 1), ('II', 2)):
        pool = _load_pool(ion)
        abunds, diag = [], []
        for wl, elo, lgf, ew in pool:
            try:
                a = _lte_abund(code, wl, elo, lgf, ew)
            except Exception as e:
                print(f"    Ti {ion} {wl}: inversion failed ({e})"); continue
            if not np.isfinite(a):
                continue                                   # railed (outside COG) — dropped
            abunds.append(a)
            if round(wl, 2) in {round(x, 2) for x in diag_wl}:
                diag.append(a)
        abunds = np.array(abunds)
        med = float(np.median(abunds)) if len(abunds) else float('nan')
        mad = float(np.median(np.abs(abunds - med))) if len(abunds) else float('nan')
        out[ion] = (med, mad, len(abunds))
        n_pool = len(pool)
        print(f"  Ti {ion}: n={len(abunds)}/{n_pool} lines resolved  A_LTE = {med:.3f} (MAD {mad:.3f})")
        if ion == 'I' and diag:
            print(f"         clean RYA-542 diag lines (good gf): A_LTE median = {np.median(diag):.3f} (n={len(diag)})")
    a1_lte, _, n1 = out['I']; a2_lte, _, n2 = out['II']
    a1_nlte = a1_lte + delta_ti1
    a2_nlte = a2_lte + TIII_NLTE
    bal_lte = a1_lte - a2_lte
    bal_nlte = a1_nlte - a2_nlte
    print(f"\n  A(Ti I)_LTE  = {a1_lte:.3f}   A(Ti II)_LTE = {a2_lte:.3f}   LTE  balance = {bal_lte:+.3f}")
    print(f"  A(Ti I)_NLTE = {a1_nlte:.3f} (LTE {a1_lte:.3f} + delta {delta_ti1:+.3f})   "
          f"A(Ti II)_NLTE = {a2_nlte:.3f}")
    print(f"  NLTE ionization balance A(Ti I)_NLTE - A(Ti II) = {bal_nlte:+.3f}  "
          f"(gate |.| < {FE_IONISATION_GATE})")
    accept = abs(bal_nlte) < FE_IONISATION_GATE
    if accept:
        print(f"\n  VERDICT: ACCEPT — the Mallinson-2024 grid brings OUR Ti I into balance with OUR Ti II.")
    else:
        print(f"\n  VERDICT: STOP — OUR Ti I(NLTE) does NOT balance OUR Ti II (|{bal_nlte:+.3f}| > "
              f"{FE_IONISATION_GATE}). Surface: our Ti I/Ti II gf-scale or pool is the problem, "
              f"not the NLTE grid. Do NOT wire until resolved.")
    return {'accept': accept, 'balance_nlte': bal_nlte, 'balance_lte': bal_lte,
            'a1_lte': a1_lte, 'a2_lte': a2_lte, 'delta_ti1': delta_ti1, 'n1': n1, 'n2': n2}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--derive', action='store_true')
    ap.add_argument('--gate', action='store_true')
    ap.add_argument('--delta', type=float, help='Ti I NLTE delta for the gate (else derive it)')
    a = ap.parse_args()
    d = None
    if a.derive or (a.gate and a.delta is None):
        d = derive()
    if a.gate:
        _register_ti()
        gate(a.delta if a.delta is not None else d['delta_median'])
