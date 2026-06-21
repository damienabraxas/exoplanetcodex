"""
scripts/fetch_inspect_nlte.py
=============================
Build NLTE-correction grids from the INSPECT database (www.inspect-stars.com) by
driving its query endpoints programmatically (no GUI). Na (Lind 2011), Mg (Osorio
2015); Li (Lind 2009) is deferred (Tier-3 + RYA-103 blend).

Two input modes (RYA-396 — Ryan: support both, the scraper is reused for the
multi-star phase; pick one per element and apply it IDENTICALLY Sun↔star so the
model-transfer systematic cancels in the differential [X/H], RYA-282):
  --mode ew : endpoint A_from_e         (input EW  → A_LTE, A_NLTE, Δ)
  --mode a  : endpoint nonlte_from_lte  (input A_LTE → A_NLTE, Δ)

KEY CORRECTION (RYA-396): INSPECT **does** take a metallicity input — the form field
`f` is [Fe/H], a REAL axis (verified live: Na D2 Δ moves −0.117→−0.060 dex across
[Fe/H] −1.0→+0.3). The RYA-241b/242 premise "INSPECT has no [Fe/H]" was a mis-read of
the form (the scraper never sent `f`), which is why 242 *faked* the axis with
`δ×(1−0.04·feh)`. That hack is unnecessary and forbidden: we query the real (Teff,
logg, [Fe/H]) grid directly — matching the RYA-235 interpolator natively.

Endpoint params: element_name, e|A_lte, t(Teff), g(logg), f([Fe/H]), x(vmic), wi(line).
Result row (last 5 floats): EW, A_LTE, A_NLTE, Δ, [Fe/H]_NLTE.

Line source (RYA-396): committed data/measured/sol_ew_results_v1.csv — the element's
measured lines (resolved wavelengths + EWs), matched to INSPECT's wi dropdown.

Output: data/nlte_grids/{El}_{Ref}_INSPECT.csv (+ .prov.json), schema
  element, ion, wave_A, teff_K, logg, feh, delta_nlte.

Caveat (logged): Δ is on INSPECT's model atmosphere, applied to our ATLAS9/Castelli
A_LTE — standard differential-transfer, ~0.02–0.05 dex systematic (same as MPIA, RYA-245).

Usage:
  python scripts/fetch_inspect_nlte.py --element Na --mode ew
  python scripts/fetch_inspect_nlte.py --element Na --mode a
"""
import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator

ROOT = Path(__file__).resolve().parents[1]
MEASURED = ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv'
OUT_DIR = ROOT / 'data' / 'nlte_grids'
BASE = 'https://www.inspect-stars.com'

ENDPOINT = {'ew': 'A_from_e', 'a': 'nonlte_from_lte'}
STRENGTH_PARAM = {'ew': 'e', 'a': 'A_lte'}

ELEMENTS = {
    'Na': {'ion': 'I', 'ref': 'Lind et al. 2011',   'outfile': 'Na_Lind2011_INSPECT.csv'},
    'Mg': {'ion': 'I', 'ref': 'Osorio et al. 2015', 'outfile': 'Mg_Osorio2015_INSPECT.csv'},
}

# Real axes — f is [Fe/H] (the RYA-396 correction). Same nodes as the MPIA grids.
TEFF_GRID = [5000, 5500, 5750, 5772, 6000, 6500]
LOGG_GRID = [3.5, 4.0, 4.4, 4.5]
FEH_GRID  = [-0.5, 0.0, 0.3]
VMIC = 1.0
MATCH_TOL_A = 0.15
SOLAR = (5772.0, 4.438, 0.0)
SLEEP_S = 0.1


def _get(path):
    req = urllib.request.Request(BASE + path, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=25) as r:
        return r.read().decode('utf-8', 'replace')


def line_options(element):
    """wi → wavelength (Å) from the element's INSPECT dropdown."""
    h = _get(f'/A_from_e?element_name={element}')
    return {int(wi): float(w) for wi, w in
            re.findall(r'<option value="(\d+)"[^>]*>\s*([0-9.]+)', h)}


def query(element, mode, t, g, f, wi, strength, x=VMIC):
    """One INSPECT query → (ew, a_lte, a_nlte, delta) or None on miss."""
    params = {'element_name': element, STRENGTH_PARAM[mode]: strength,
              't': t, 'g': g, 'f': f, 'x': x, 'wi': wi}
    try:
        h = _get(f'/{ENDPOINT[mode]}?{urllib.parse.urlencode(params)}')
    except Exception as e:
        print(f'  HTTP {element} t{t} g{g} f{f} wi{wi}: {e}')
        return None
    nums = re.findall(r'-?\d+\.\d+', re.sub(r'<[^>]+>', ' ', h))
    if len(nums) < 5:
        return None
    ew, a_lte, a_nlte, delta = (float(nums[-5]), float(nums[-4]),
                                float(nums[-3]), float(nums[-2]))
    return ew, a_lte, a_nlte, delta


def matched_lines(element):
    """Element's measured lines (wavelength, EW) matched to INSPECT wi options."""
    if not MEASURED.exists():
        sys.exit(f'CRITICAL: measured-line source missing: {MEASURED}')
    ew = pd.read_csv(MEASURED)
    ion = ELEMENTS[element]['ion']
    ours = ew[(ew['element'] == element) & (ew['ion'] == ion)]
    opts = line_options(element)
    out = []
    for _, r in ours.iterrows():
        w = float(r['wavelength_air_A'])
        wi, best = min(opts.items(), key=lambda kv: abs(kv[1] - w))
        if abs(best - w) <= MATCH_TOL_A:
            out.append({'wi': wi, 'wave_A': best, 'our_wave_A': round(w, 3),
                        'ew_mA': float(r['ew_mA'])})
    return out


def build(element, mode):
    spec = ELEMENTS[element]
    print(f'=== {element} {spec["ion"]} NLTE grid — INSPECT {ENDPOINT[mode]} '
          f'({spec["ref"]}) ===')
    lines = matched_lines(element)
    if not lines:
        sys.exit(f'CRITICAL: no measured {element} lines matched INSPECT options')
    print(f'matched {len(lines)} measured {element} lines: '
          f'{[l["wave_A"] for l in lines]}')

    # representative strength per line: measured solar EW (ew-mode) or the
    # solar A_LTE INSPECT reports for that EW at solar params (a-mode anchor).
    rows = []
    for ln in lines:
        if mode == 'ew':
            strength = ln['ew_mA']
        else:
            anchor = query(element, 'ew', SOLAR[0], SOLAR[1], 0.0, ln['wi'], ln['ew_mA'])
            if anchor is None:
                print(f'  {ln["wave_A"]}: could not anchor A_LTE — skipped'); continue
            strength = anchor[1]                      # solar A_LTE for this line
        for t in TEFF_GRID:
            for g in LOGG_GRID:
                for f in FEH_GRID:
                    res = query(element, mode, t, g, f, ln['wi'], strength)
                    time.sleep(SLEEP_S)
                    if res is None:
                        continue
                    rows.append({'element': element, 'ion': spec['ion'],
                                 'wave_A': ln['our_wave_A'], 'teff_K': t, 'logg': g,
                                 'feh': f, 'delta_nlte': round(res[3], 4),
                                 'strength_mode': mode, 'strength': round(strength, 3)})
    grid = pd.DataFrame(rows)
    if grid.empty:
        sys.exit(f'CRITICAL: {element} INSPECT returned no rows')

    out = OUT_DIR / spec['outfile']
    grid[['element', 'ion', 'wave_A', 'teff_K', 'logg', 'feh', 'delta_nlte']].to_csv(out, index=False)
    print(f'saved {len(grid)} rows → {out}')

    # solar anchor (mean over lines at the solar grid point)
    g0 = grid.groupby(['teff_K', 'logg', 'feh'])['delta_nlte'].mean().reset_index()
    f_int = LinearNDInterpolator(g0[['teff_K', 'logg', 'feh']].values, g0['delta_nlte'].values)
    solar_delta = float(f_int([SOLAR])[0])
    print(f'solar anchor Δ = {solar_delta:+.4f} dex (Teff={SOLAR[0]}, logg={SOLAR[1]}, feh={SOLAR[2]})')

    prov = {
        'element': element, 'ion': spec['ion'], 'reference': spec['ref'],
        'source': 'INSPECT database (inspect-stars.com)',
        'endpoint': f'{BASE}/{ENDPOINT[mode]}', 'input_mode': mode,
        'convention': 'delta_nlte = A_NLTE - A_LTE (dex)',
        'feh_axis': 'REAL — INSPECT form field `f` is [Fe/H] (RYA-396 correction; not faked)',
        'strength': ('measured solar EW per line (mA)' if mode == 'ew'
                     else 'solar A_LTE per line (anchored via ew-mode)'),
        'line_source': str(MEASURED.relative_to(ROOT)),
        'lines': [{'our_wave_A': l['our_wave_A'], 'inspect_wave_A': l['wave_A'],
                   'wi': l['wi'], 'ew_mA': l['ew_mA']} for l in lines],
        'grid_axes': {'teff_K': TEFF_GRID, 'logg': LOGG_GRID, 'feh': FEH_GRID},
        'vmic_kms': VMIC, 'extraction_date': date.today().isoformat(),
        'n_rows': int(len(grid)),
        'caveat': 'Δ on INSPECT model atmosphere applied to our A_LTE — ~0.02-0.05 dex '
                  'differential-transfer systematic (RYA-245), cancels Sun↔star if used consistently.',
        'solar_anchor': {'teff_K': SOLAR[0], 'logg': SOLAR[1], 'feh': SOLAR[2],
                         'delta_nlte': round(solar_delta, 4)},
    }
    out.with_suffix('.prov.json').write_text(json.dumps(prov, indent=2))
    print(f'saved provenance → {out.with_suffix(".prov.json")}')
    return solar_delta


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--element', required=True, choices=sorted(ELEMENTS))
    ap.add_argument('--mode', default='ew', choices=['ew', 'a'])
    a = ap.parse_args()
    build(a.element, a.mode)
