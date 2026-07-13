#!/usr/bin/env python3
"""
RYA-546 Mn re-derivation — Gate 0 + triplet-exact ab-initio δ.

Reuses RYA-473's HFS machinery + RYA-476's Amarsi GALAH grid wiring (no new synthesis path):
derives Δ(Mn) = A(Mn)_NLTE − A(Mn)_LTE on the Den Hartog e6S→z6P triplet (6013/16/21) from the
LIVE `nlte_Mn_scatt_pysme.grd`, HFS-resolved IDENTICALLY on both legs (the harness `_synth_ew`
expands the HFS components for the NLTE EW and the LTE COG alike — Addition A). The harness
`nlte_delta` now runs the RYA-546 in-hull guard (Addition B) before interpolating.

GATE 0: the triplet must return a NON-NaN δ (these lines NaN'd on the MPIA high-EP grid, RYA-476).
Expect order ~+0.02–0.04 (≈ the Engine-B TS-Gerber +0.043), NOT +0.108. If NaN / grid offline →
STOP, no silent fallback to +0.108. Solar node from STAR_PARAMS (`get_star_params`, no hardcode).
Runs on Sirius (venv_pysme, grid Sirius-only).
"""
from __future__ import annotations
import sys
from pathlib import Path
REPO = str(Path(__file__).resolve().parents[1])
if REPO not in sys.path:
    sys.path.insert(0, REPO)
from pipeline import _runtime as _rt   # RYA-514 force-fork + single-thread BLAS (before numpy)
import numpy as np

from config.constants import SOLAR_ASPLUND2021, get_star_params
from pipeline import pysme_nlte as P
from scripts.resolve_camn_stops_rya411 import _build_lines

TRIPLET = [6013.51, 6016.67, 6021.82]     # Den Hartog e6S→z6P (RYA-473); same lower level (EP ~3.07)
MPIA_HIEP = [6304.79, 6306.29, 6867.90]   # the high-EP MPIA-node lines (RYA-411) — for the EP trend
MN_ENGINE_A = 0.108                        # registered MPIA/Bergemann scaled-Drawin δ (the suspect)
MN_ENGINE_B = 0.043                        # RYA-534 TS-Gerber ab-initio corroboration


def _solar_star():
    p = get_star_params('solar')
    teff = float(p.get('teff', p.get('teff_K', 5772)))
    logg = float(p.get('logg', 4.44))
    feh = float(p.get('feh', p.get('feh_ref', 0.0)))
    vmic = float(p.get('xi', p.get('vmic', p.get('vturb_kms', 1.0))))
    return {'teff': teff, 'logg': logg, 'feh': feh, 'vmic': vmic}


def _derive(waves, star):
    P._A_SUN.setdefault('Mn', SOLAR_ASPLUND2021['Mn'])
    lines, _ = _build_lines('Mn', waves, hfs=True)   # HFS-resolved, both legs (Addition A)
    P.NLTE_LINES['Mn'] = lines
    return P.nlte_delta('Mn', star=star)


def main():
    star = _solar_star()
    print("=== RYA-546 Mn Gate 0 — triplet-exact ab-initio δ (live Amarsi GALAH grid, HFS both legs) ===")
    print(f"  solar node (get_star_params, no hardcode): {star}")
    print(f"  context: Engine-A MPIA scaled-Drawin δ={MN_ENGINE_A:+.3f} (suspect); "
          f"Engine-B TS-Gerber ab-initio δ={MN_ENGINE_B:+.3f} (RYA-534)")

    res = _derive(TRIPLET, star)
    print("\n  TRIPLET (Den Hartog e6S→z6P, low-EP ~3.07 eV):")
    for wl, d in sorted(res['per_line'].items()):
        print(f"    Mn I {wl:.3f}: δ = {d:+.4f}")
    med = res['delta_median']
    print(f"  triplet MEDIAN δ = {med:+.4f}")

    nan = any(not np.isfinite(v) for v in res['per_line'].values())
    print(f"\n  GATE 0: triplet δ {'contains NaN — STOP (grid offline)' if nan else 'NON-NaN — PASS'}")
    if not nan:
        print(f"  vs Engine-A +0.108 (scaled-Drawin): {'MUCH SMALLER (ab-initio confirmed)' if med < 0.07 else 'NOT smaller — CHECK'}; "
              f"vs Engine-B +0.043: Δ={med-MN_ENGINE_B:+.3f}")

    # EP trend (Addition C): the high-EP MPIA-node lines should give a LARGER δ than the low-EP
    # triplet — a FLAT δ across EP would be the red flag that grid/HFS handling is wrong.
    try:
        hi = _derive(MPIA_HIEP, star)
        print("\n  HIGH-EP MPIA-node lines (EP trend check):")
        for wl, d in sorted(hi['per_line'].items()):
            print(f"    Mn I {wl:.3f}: δ = {d:+.4f}")
        print(f"  high-EP MEDIAN δ = {hi['delta_median']:+.4f}  "
              f"(EP trend low→high: {med:+.3f} → {hi['delta_median']:+.3f}; "
              f"{'trend present (good)' if abs(hi['delta_median']-med) > 0.01 else 'FLAT — red flag'})")
    except Exception as e:
        print(f"\n  high-EP EP-trend check skipped: {type(e).__name__}: {e}")


if __name__ == '__main__':
    main()
