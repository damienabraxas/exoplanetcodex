"""
scripts/resource_clamped_grids_rya409.py
========================================
RYA-409 Part B — re-source the 55-Cnc-clamped Family-A grids (Na/Mg/Si/Ca/Mn) onto
the validated PySME b-factor path (RYA-402), over the in-repo Amarsi-2020 grids which
cover [Fe/H] -> +1. Closes the metal-rich clamp for 5 of the 7 (Cr/Ti not in Amarsi).

For each element:
  1. resolve the diagnostic line atomic data (loggf/EP/Eup/ABO vdW) from the GES
     synthesis linelist (HFS summed), and the grid level labels by nearest energy
     (pysme_nlte.auto_labels) -- no hand-mapping, no hardcoded line data;
  2. derive the SOLAR delta via PySME (pipeline.pysme_nlte.nlte_delta);
  3. CROSS-CHECK (validate-don't-tune): the PySME solar delta MUST reproduce the
     currently-registered MPIA/INSPECT solar value within tol -- else STOP and surface
     it as an MPIA-vs-Amarsi finding (never silently swap engines);
  4. on PASS: generate the full grid over a node set spanning solar + 55 Cnc + metal
     -rich ([Fe/H] to +0.6), write {El}_Amarsi2020_PySME.csv + provenance.

Registration (constants + Beast YAML) is done separately after review of the cross-checks.

    python -m scripts.resource_clamped_grids_rya409 --elements Na,Mg,Si,Ca,Mn
    python -m scripts.resource_clamped_grids_rya409 --elements Na --crosscheck-only
"""
from __future__ import annotations
import argparse, sys, warnings, json
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
from config.constants import ISPEC_DIR, NLTE_CORRECTION_ELEMENTS  # noqa: E402
from pipeline import pysme_nlte as P                              # noqa: E402

# Diagnostic wavelengths (air, A) per element — standard NLTE optical lines.
DIAG_WAVES = {
    'Na': [5682.633, 5688.205],
    'Mg': [5711.088, 4730.029],
    'Si': [5772.146, 6125.021],
    'Ca': [6166.439, 6169.563, 6455.598],
    'Mn': [6013.51, 6016.67, 6021.80],     # the Mn I triplet (HFS — summed per feature)
}
CROSSCHECK_TOL = 0.04   # PySME-Amarsi solar delta must reproduce registered MPIA within this

_atomic_cache = {}


def _ges():
    if 'll' not in _atomic_cache:
        sys.path.insert(0, str(ISPEC_DIR)); import ispec
        import pipeline.abundances_derive as ad
        ll = ispec.read_atomic_linelist(ad._SYNTH_LINELIST_FILE)
        _atomic_cache['ll'] = ll
        _atomic_cache['notes'] = np.array([str(x) for x in ll['element']])
        _atomic_cache['w'] = np.asarray(ll['wave_A'], float)
    return _atomic_cache['ll'], _atomic_cache['notes'], _atomic_cache['w']


def _resolve_atomic(element, wl, ion='1', dw=0.25):
    """(loggf_total, EP, Eup, gamvw) for the feature at `wl` from GES (HFS summed)."""
    ll, notes, w = _ges()
    m = np.where((notes == f'{element} {ion}') & (np.abs(w - wl) < dw))[0]
    if len(m) == 0:
        raise ValueError(f"{element} {wl}: not in GES linelist (±{dw} A).")
    gf = np.array([10 ** float(ll['loggf'][i]) for i in m])
    ep = float(np.median([float(ll['lower_state_eV'][i]) for i in m]))
    vw = float(np.median([float(ll['waals'][i]) for i in m]))
    loggf = float(np.log10(gf.sum()))
    eup = ep + 12398.42 / wl
    return loggf, ep, eup, (vw if vw > 0 else 0.0)


def _registered_solar_delta(element):
    """Median solar delta of the CURRENTLY registered grid (the value to reproduce)."""
    from pipeline.nlte_corrections import _load_mpia_element_grid, _mpia_element_delta
    c = _load_mpia_element_grid(element)
    ds = [_mpia_element_delta(element, float(w), 5772, 4.44, 0.0) for w in c['waves']]
    ds = [d for d in ds if np.isfinite(d)]
    return float(np.median(ds)) if ds else np.nan


_NODES = [(5772, 4.44, 0.0), (5100, 4.44, 0.0), (6200, 4.44, 0.0), (5772, 4.0, 0.0),
          (5772, 4.7, 0.0), (5772, 4.44, -0.3), (5172, 4.43, 0.31), (5100, 4.5, 0.5),
          (6000, 4.3, -0.2), (5772, 4.44, 0.6), (5400, 4.5, 0.5)]


def process(element, crosscheck_only=False):
    print(f"\n══ {element} ══")
    lines = []
    for wl in DIAG_WAVES[element]:
        loggf, ep, eup, vw = _resolve_atomic(element, wl)
        tl, tu, jl, ju = P.auto_labels(element, ep, eup)
        lines.append((wl, loggf, ep, jl, eup, ju, tl, tu, vw))
        print(f"  {wl:.3f}: loggf={loggf:+.3f} EP={ep:.3f} -> {tl} -> {tu} (J{jl:.0f}->{ju:.0f}) vdW={vw:.1f}")
    P.NLTE_LINES[element] = lines
    from config.constants import SOLAR_ASPLUND2021
    P._A_SUN.setdefault(element, SOLAR_ASPLUND2021[element])

    res = P.nlte_delta(element, star={'teff': 5772, 'logg': 4.44, 'feh': 0.0, 'vmic': 1.0})
    d_pysme = res['delta_median']
    d_reg = _registered_solar_delta(element)
    diff = d_pysme - d_reg
    ok = abs(diff) <= CROSSCHECK_TOL
    print(f"  CROSS-CHECK: PySME-Amarsi {d_pysme:+.3f}  vs registered MPIA {d_reg:+.3f}  "
          f"(diff {diff:+.3f}, tol {CROSSCHECK_TOL}) -> {'PASS' if ok else 'STOP — MPIA-vs-Amarsi DISCREPANCY'}")
    if not ok:
        print(f"  !! {element}: PySME-Amarsi does NOT reproduce the registered value — STOPPING, "
              f"surface for adjudication (do NOT swap engines silently).")
        return {'element': element, 'd_pysme': d_pysme, 'd_reg': d_reg, 'passed': False}
    if crosscheck_only:
        return {'element': element, 'd_pysme': d_pysme, 'd_reg': d_reg, 'passed': True}

    rows = []
    for (te, lg, fe) in _NODES:
        r = P.nlte_delta(element, star={'teff': te, 'logg': lg, 'feh': fe, 'vmic': 1.0})
        for wl, dd in r['per_line'].items():
            rows.append(dict(element=element, ion=1, wave_A=round(wl, 3),
                             teff_K=te, logg=lg, feh=fe, delta_nlte=round(float(dd), 4)))
    out = _REPO / 'data' / 'nlte_grids' / f'{element}_Amarsi2020_PySME.csv'
    pd.DataFrame(rows).to_csv(out, index=False)
    feh_max = max(n[2] for n in _NODES)
    print(f"  wrote {out.name} ({len(rows)} rows, [Fe/H] up to {feh_max:+.1f} — covers 55 Cnc +0.31)")
    return {'element': element, 'd_pysme': d_pysme, 'd_reg': d_reg, 'passed': True, 'grid': out.name}


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('--elements', default='Na,Mg,Si,Ca,Mn')
    p.add_argument('--crosscheck-only', action='store_true')
    a = p.parse_args(argv)
    results = [process(el, a.crosscheck_only) for el in a.elements.split(',')]
    print("\n── SUMMARY ──")
    for r in results:
        print(f"  {r['element']:3} PySME {r['d_pysme']:+.3f} vs MPIA {r['d_reg']:+.3f}  "
              f"{'PASS' if r['passed'] else 'STOP'}  {r.get('grid','')}")
    if any(not r['passed'] for r in results):
        raise SystemExit("RYA-409 Part B: cross-check STOP — adjudicate MPIA-vs-Amarsi before swapping.")


if __name__ == '__main__':
    main()
