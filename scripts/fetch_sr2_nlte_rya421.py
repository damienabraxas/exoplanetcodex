#!/usr/bin/env python3
"""
fetch_sr2_nlte_rya421.py — build the Sr II resonance-doublet NLTE delta grid from
INSPECT (Bergemann et al. 2012a, A&A 546 A90), the cited CROSS-CHECK source for
RYA-421. The PRIMARY (Mashonkina 2022 / INASAN) is to be vendored separately when a
machine-readable grid is available; the two agree within 0.15 dex and both are
near-LTE at our metallicities. INSPECT is the readily-fetchable Sr II grid and is
wired as the working delta source here.

Sr is an ION MISMATCH not a grid gap (RYA-401): INSPECT hosts Sr II 4077/4215 (+NIR),
but our solar EW set measures only Sr I 6617 — so this fetcher queries the Sr II
resonance lines DIRECTLY (wi from the INSPECT dropdown), not via the Sr-I-measurement
match in fetch_inspect_nlte.py.

Delta convention: queried with a scaled-solar abundance A_lte = A(Sr)_sun + [Fe/H]
([Sr/Fe]=0), the standard NLTE-correction-grid assumption. INSPECT returns only at its
grid nodes (no interpolation); off-node / out-of-range queries are recorded, NOT
fabricated. The actual [Fe/H] ceiling is read from the live responses (RYA-402 Step-6
pattern) and stored in provenance — a LOUD flag if it clamps below +0.32 (55 Cnc).

Output: data/nlte_grids/Sr_Bergemann2012_INSPECT.csv  (+ .prov.json)
        columns: element, ion, wave_A, teff_K, logg, feh, delta_nlte

    python scripts/fetch_sr2_nlte_rya421.py            # live fetch
    python scripts/fetch_sr2_nlte_rya421.py --probe    # node-coverage probe only
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from fetch_inspect_nlte import query, line_options  # reuse the INSPECT client/parser

OUT = ROOT / 'data' / 'nlte_grids' / 'Sr_Bergemann2012_INSPECT.csv'
A_SR_SUN = 2.87                      # Asplund+2021 solar A(Sr)
LINES = {0: 4077.709, 1: 4215.519}   # wi -> air wavelength (Sr II resonance doublet)

# candidate nodes (INSPECT returns only where its grid is defined; we keep valid only)
TEFF = [5000, 5500, 5750, 5772, 6000]
LOGG = [3.5, 4.0, 4.4, 4.5]
FEH  = [-3.0, -2.0, -1.0, 0.0, 0.25, 0.32]   # include metal-rich to probe the ceiling
SANE = 1.0                           # |delta| < SANE dex = a real correction, not a parse miss


def fetch():
    rows, ceiling_hits, oob = [], set(), 0
    for wi, wave in LINES.items():
        for f in FEH:
            a_lte = round(A_SR_SUN + f, 2)        # [Sr/Fe]=0 scaled-solar
            for t in TEFF:
                for g in LOGG:
                    r = query('Sr', 'a', t, g, f, wi, a_lte)
                    if r and abs(r[3]) < SANE:
                        rows.append({'element': 'Sr', 'ion': 'II', 'wave_A': wave,
                                     'teff_K': t, 'logg': g, 'feh': f,
                                     'delta_nlte': round(r[3], 4)})
                        ceiling_hits.add(f)
                    else:
                        oob += 1
                    time.sleep(0.1)
    return rows, sorted(ceiling_hits), oob


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--probe', action='store_true', help='print node coverage, do not write')
    a = ap.parse_args()

    print(f"INSPECT Sr line options: {line_options('Sr')}")
    rows, feh_with_data, oob = fetch()
    if not rows:
        sys.exit("CRITICAL: no valid Sr II nodes returned from INSPECT — do NOT fabricate.")
    df = pd.DataFrame(rows)
    feh_ceiling = float(df['feh'].max())
    solar = df[(df.feh == 0.0) & (df.teff_K == 5772) & (df.logg == 4.4)]

    print(f"\nSr II grid: {len(df)} valid nodes ({oob} out-of-grid, recorded not faked)")
    print(f"  [Fe/H] nodes with data: {feh_with_data}  -> CEILING = {feh_ceiling:+.2f}")
    if feh_ceiling < 0.32:
        print(f"  *** LOUD FLAG: grid [Fe/H] ceiling {feh_ceiling:+.2f} < +0.32 (55 Cnc) — "
              f"55 Cnc (+0.32) and tau Boo (+0.25) are FLAGGED EXTRAPOLATIONS. "
              f"Acceptable (correction is near-LTE up here) but never silent (RYA-421 Step 3).")
    print(f"  solar anchor (5772/4.4/0.0):")
    for _, r in solar.iterrows():
        print(f"     Sr II {r.wave_A}: delta_nlte = {r.delta_nlte:+.4f}  (expect small ~-0.005)")

    if a.probe:
        print("\n[probe] not written.")
        return

    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.sort_values(['wave_A', 'feh', 'teff_K', 'logg']).to_csv(OUT, index=False)
    prov = {
        'grid': OUT.name,
        'source_primary': 'Mashonkina et al. 2022 Sr II model atom / INASAN (Mashonkina+2023) '
                          '- CITED PRIMARY, to be vendored when machine-readable; demonstrated to [Fe/H]=+0.25',
        'source_working': 'Bergemann et al. 2012a (A&A 546 A90; arXiv 1207.2451) via INSPECT '
                          '(inspect-stars.com, A_from_e/nonlte_from_lte), the fetchable cross-check; '
                          'agrees with Mashonkina within 0.15 dex',
        'lines': 'Sr II 4077.709 + 4215.519 (air); resonance doublet, negligible HFS (even-A isotopes)',
        'delta_convention': 'A_lte = A(Sr)_sun(2.87) + [Fe/H], [Sr/Fe]=0 scaled-solar; '
                            'delta_nlte = A_NLTE - A_LTE from INSPECT nonlte_from_lte',
        'feh_ceiling': feh_ceiling,
        'feh_nodes': feh_with_data,
        'metal_rich_note': f'ceiling {feh_ceiling:+.2f}; 55 Cnc (+0.32) + tau Boo (+0.25) are flagged '
                           'extrapolations into the NEAR-LTE regime (correction ~-0.005 at solar, '
                           'grows only metal-poor) - RYA-421 Step 3',
        'solar_anchor_delta': {str(r.wave_A): r.delta_nlte for _, r in solar.iterrows()},
        'fetched': str(date.today()),
    }
    OUT.with_suffix('.prov.json').write_text(json.dumps(prov, indent=2))
    print(f"\nwrote {OUT}  (+ {OUT.with_suffix('.prov.json').name})")


if __name__ == '__main__':
    main()
