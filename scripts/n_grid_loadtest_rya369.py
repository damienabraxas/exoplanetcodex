"""
scripts/n_grid_loadtest_rya369.py
=================================
RYA-369 — load-test the N I NLTE grid (Amarsi/GALAH PySME `nlte_N_scatt_pysme.grd`,
Zenodo 3982506, staged to Sirius under RYA-477 and to the Mac transiently for this test).
The strategy gate: N's hard blocker (the N I NLTE grid) clears IFF this passes.

Three checks (Ryan's greenlit list):
  1. FORMAT/WIRING — the pipeline reads the PySME `.grd` (already true for 9 elements
     incl. Mn, RYA-476). Register N in pysme_nlte and confirm read_amarsi_grid('N') loads.
  2. RESOLVE a real N I line — N I 7468 (+ 8216/8683) at solar params (5772/4.44/0.0)
     returns a SANE, NEGATIVE, sane-magnitude NLTE delta (sign-discipline, RYA-339). High-
     excitation N I is expected to carry a real (negative) NLTE correction.
  3. COVERAGE — the grid axes span Sun -> Procyon (6554) -> aCenA/B (+0.20/+0.30 dex) ->
     55 Cnc (+0.32); report min/max, NO silent clamp.

validate-don't-tune: this only verifies the grid loads + returns a sane delta; it does not
measure or tune A(N).

    python -m scripts.n_grid_loadtest_rya369
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import pysme_nlte as P                              # noqa: E402
from pipeline.nlte_bfactor_synth import read_amarsi_grid          # noqa: E402
from scripts.resolve_camn_stops_rya411 import _build_lines        # noqa: E402

RESULTS = ROOT / 'data' / 'results'
N_SUN = 7.83                                  # Asplund 2021 A(N), the COG zero point (ref only)
N_LINES = [7468.31, 8216.34, 8683.40]         # primary N I atomic red-optical (RYA-369)
# stars to bracket coverage: (teff, logg, feh)
COVER = {'Sun': (5772, 4.44, 0.0), 'Procyon': (6554, 4.0, 0.03),
         'aCenA': (5790, 4.31, 0.20), 'aCenB': (5230, 4.53, 0.23), '55Cnc': (5196, 4.45, 0.32)}


def _rule(c='='):
    print(c * 78)


def main():
    rep = {'ticket': 'RYA-369 N grid load-test'}
    _rule(); print("  RYA-369 — N I NLTE grid load-test (Amarsi GALAH PySME)"); _rule()

    # N is now registered in production pysme_nlte (RYA-369): _GRID_FILENAME['N'],
    # _A_SUN['N']=7.83, NLTE_LINES['N'] (7468/8216/8683), and the _Z atomic-number map
    # (the one gap the first load-test found: _Z lacked N=7). Confirm the wiring is live.
    assert P._GRID_FILENAME.get('N') == 'nlte_N_scatt_pysme.grd', 'N grid filename not registered'
    assert P._Z('N') == 7, 'N missing from the _Z atomic-number map'
    assert 'N' in P.NLTE_LINES and 'N' in P._A_SUN, 'N not registered in NLTE_LINES/_A_SUN'
    rep['production_wiring'] = "registered: _GRID_FILENAME/_A_SUN/NLTE_LINES['N'] + _Z 'N':7 (RYA-369)"

    # ── 1. FORMAT / WIRING ─────────────────────────────────────────────────────
    print("\n  [1] format/wiring — read the PySME .grd")
    try:
        g = read_amarsi_grid('N')
        teff = np.asarray(g.get('teff'), float); grav = np.asarray(g.get('grav'), float)
        feh = np.asarray(g.get('feh'), float); nlev = g.nlevels() if hasattr(g, 'nlevels') else None
        print(f"      LOADED ✓  PySME DirectAccessFile; {nlev} model-atom levels; "
              f"axes teff[{teff.min():.0f}-{teff.max():.0f}] logg[{grav.min():.2f}-{grav.max():.2f}] "
              f"feh[{feh.min():+.2f}-{feh.max():+.2f}]")
        rep['format'] = dict(loaded=True, nlevels=int(nlev) if nlev else None,
                             teff=[float(teff.min()), float(teff.max())],
                             logg=[float(grav.min()), float(grav.max())],
                             feh=[float(feh.min()), float(feh.max())])
    except Exception as e:
        print(f"      FAILED: {type(e).__name__}: {e}")
        rep['format'] = dict(loaded=False, error=f"{type(e).__name__}: {e}")
        (RESULTS / 'n_grid_loadtest_rya369.json').write_text(json.dumps(rep, indent=2, default=str))
        return rep

    # ── 3. COVERAGE (cheap, do before the slow synthesis) ──────────────────────
    print("\n  [3] coverage — benchmark set inside the grid box?")
    tmin, tmax = rep['format']['teff']; gmin, gmax = rep['format']['logg']; fmin, fmax = rep['format']['feh']
    cov = {}
    for name, (te, lg, fe) in COVER.items():
        inside = (tmin <= te <= tmax) and (gmin <= lg <= gmax) and (fmin <= fe <= fmax)
        cov[name] = inside
        print(f"      {name:8} ({te},{lg},{fe:+.2f}) -> {'inside ✓' if inside else 'OUTSIDE ⚠'}")
    rep['coverage'] = cov
    rep['coverage_all_inside'] = all(cov.values())

    # ── 2. RESOLVE a real N I line -> sane negative delta ──────────────────────
    print("\n  [2] resolve N I lines -> NLTE delta at solar (sane + negative?)")
    try:
        lines, _ = _build_lines('N', N_LINES, hfs=False)         # auto_labels reads the grid levels
        P.NLTE_LINES['N'] = lines
        star = {'teff': 5772, 'logg': 4.44, 'feh': 0.0, 'vmic': 1.0}
        res = P.nlte_delta('N', star=star)
        per = {float(k): round(float(v), 4) for k, v in res['per_line'].items()}
        med = round(float(res['delta_median']), 4)
        for wl, d in per.items():
            print(f"      N I {wl}: δ {d:+.4f}")
        print(f"      median δ = {med:+.4f}")
        sane = np.isfinite(med) and -0.6 < med < 0.05        # negative, sane magnitude
        print(f"      -> {'SANE + negative ✓ (sign-discipline RYA-339)' if sane else 'OUT OF RANGE ⚠ — investigate'}")
        rep['resolve'] = dict(per_line=per, median=med, sane_negative=bool(sane))
    except Exception as e:
        print(f"      FAILED: {type(e).__name__}: {e}")
        rep['resolve'] = dict(error=f"{type(e).__name__}: {e}")

    # ── verdict ────────────────────────────────────────────────────────────────
    _rule()
    passed = (rep['format']['loaded'] and rep.get('coverage_all_inside')
              and rep.get('resolve', {}).get('sane_negative'))
    rep['loadtest_passed'] = bool(passed)
    print(f"  LOAD-TEST: {'PASS — N hard blocker (grid) CLEARED; primary N I atomic becomes executable' if passed else 'INCOMPLETE — see failing check'}")
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'n_grid_loadtest_rya369.json').write_text(json.dumps(rep, indent=2, default=str))
    print(f"  [out] {RESULTS / 'n_grid_loadtest_rya369.json'}")
    return rep


if __name__ == '__main__':
    main()
