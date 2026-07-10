"""
scripts/wire_li_cu_rya540.py
============================
RYA-540 finalize — derive + (conditionally) vendor the Li and Cu NLTE departure grids
onto the validated PySME b-factor path (RYA-402), same recipe as N/RYA-526 and the
RYA-409 Family-A re-source.

DISCIPLINE (validate-don't-tune): derive the SOLAR delta via PySME and REPRODUCE the
element's published solar anchor within tol BEFORE writing any production grid. If the
anchor is not reproduced, STOP and surface it — never tune, never register a grid that
does not validate.

Li: grid staged on Sirius (nlte_Li_scatt_pysme.grd, Amarsi-2020 GALAH). Li I 6707 is a
    ground-state resonance doublet; for the SUN it is weak (A(Li)~1.05) so NLTE is small.
Cu: grid nlte_Cu_caliskan_Oct2024_pysme.grd is NOT staged on Sirius and has no pinned
    source URL (ledger: Zenodo/request from Amarsi) -> derivation is BLOCKED on acquiring
    the grid; reported, not faked.

    python -m scripts.wire_li_cu_rya540 --elements Li --crosscheck-only
    python -m scripts.wire_li_cu_rya540 --elements Li,Cu --derive
"""
from __future__ import annotations
import argparse, sys, json, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
from pipeline import pysme_nlte as P   # noqa: E402

# Diagnostic-line atomic data for the new elements (well-established; the grid LEVEL
# labels are resolved at runtime by nearest energy via P.auto_labels — no hand-mapping).
#   Li I 6707.8 resonance doublet 2s 2S -> 2p 2P* (Elow=0). The two fine-structure
#   components (2P*_3/2 6707.76 loggf +0.174, 2P*_1/2 6707.91 loggf -0.127) split the
#   upper level by ~0.34 cm^-1 -> identical departures, so they are one NLTE feature;
#   gf-summed loggf = log10(10^0.174 + 10^-0.127) = +0.350. vdW: Unsold (weak solar line).
_LI_ELOW, _LI_EUP = 0.000, 1.8477
_NEW_LINES = {
    'Li': dict(wl=6707.81, loggf=0.350, elow=_LI_ELOW, eup=_LI_EUP, vw=0.0),
}
# Reuse the RYA-409 node set (solar + cool + metal-rich to +0.6, covers 55 Cnc +0.31).
_NODES = [(5772, 4.44, 0.0), (5100, 4.44, 0.0), (6200, 4.44, 0.0), (5772, 4.0, 0.0),
          (5772, 4.7, 0.0), (5772, 4.44, -0.3), (5172, 4.43, 0.31), (5100, 4.5, 0.5),
          (6000, 4.3, -0.2), (5772, 4.44, 0.6), (5400, 4.5, 0.5)]
_CSV_NAME = {'Li': 'Li_Amarsi2020_PySME.csv', 'Cu': 'Cu_Caliskan2024_PySME.csv'}


def _setup_lines(element):
    """Populate P.NLTE_LINES[element] from the diagnostic atomic data, resolving the
    grid level labels by nearest energy (raises if a level is off-grid)."""
    d = _NEW_LINES[element]
    tl, tu, jl, ju = P.auto_labels(element, d['elow'], d['eup'])
    P.NLTE_LINES[element] = [(d['wl'], d['loggf'], d['elow'], jl, d['eup'], ju, tl, tu, d['vw'])]
    print(f"  {element} {d['wl']:.2f}: loggf={d['loggf']:+.3f} EP={d['elow']:.3f} "
          f"-> {tl} -> {tu} (J{jl:.0f}->{ju:.0f})")


def _grid_present(element):
    return (P._GRID_DIR / P._GRID_FILENAME[element]).exists()


def process(element, derive=False):
    print(f"\n══ {element} ══")
    if not _grid_present(element):
        print(f"  BLOCKED: grid {P._GRID_FILENAME[element]} not staged at {P._GRID_DIR} — "
              f"cannot derive. Acquire the grid first (no pinned source URL in repo).")
        return {'element': element, 'passed': False, 'blocked': True}
    if element in _NEW_LINES:
        _setup_lines(element)
    # validate-don't-tune: reproduce the published solar anchor
    res = P.validate(element)
    print(f"  SOLAR delta (PySME) = {res['delta_median']:+.4f}  vs anchor "
          f"{res['anchor']:+.3f} +/- {res['anchor_tol']}  -> {'PASS' if res['passed'] else 'STOP'}")
    print(f"  anchor ref: {res['ref']}")
    if not res['passed']:
        print(f"  !! {element}: PySME did NOT reproduce the published solar anchor — STOPPING, "
              f"NOT registering (surface for adjudication; do not tune).")
        return {'element': element, 'passed': False, 'delta': res['delta_median']}
    if not derive:
        return {'element': element, 'passed': True, 'delta': res['delta_median']}
    # PASS -> derive the multi-node production grid
    rows = []
    for (te, lg, fe) in _NODES:
        r = P.nlte_delta(element, star={'teff': te, 'logg': lg, 'feh': fe, 'vmic': 1.0})
        for wl, dd in r['per_line'].items():
            rows.append(dict(element=element, ion=1, wave_A=round(wl, 3),
                             teff_K=te, logg=lg, feh=fe, delta_nlte=round(float(dd), 4)))
    out = _REPO / 'data' / 'nlte_grids' / _CSV_NAME[element]
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"  wrote {out.name} ({len(rows)} rows, [Fe/H] up to "
          f"{max(n[2] for n in _NODES):+.1f})  solar delta {res['delta_median']:+.4f}")
    return {'element': element, 'passed': True, 'delta': res['delta_median'],
            'grid': out.name, 'solar_delta': res['delta_median']}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--elements', default='Li')
    ap.add_argument('--derive', action='store_true', help='write the production CSV on PASS')
    ap.add_argument('--crosscheck-only', action='store_true')
    a = ap.parse_args(argv)
    results = [process(el, derive=a.derive and not a.crosscheck_only)
               for el in a.elements.split(',')]
    print("\n── SUMMARY ──")
    for r in results:
        tag = 'BLOCKED' if r.get('blocked') else ('PASS' if r['passed'] else 'STOP')
        print(f"  {r['element']:3} {tag}  delta={r.get('delta','?')}  {r.get('grid','')}")
    print(json.dumps(results, indent=2, default=float))


if __name__ == '__main__':
    main()
