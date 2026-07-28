#!/usr/bin/env python3
"""
RYA-545 — derive the production Ti_Mallinson2024_PySME.csv multi-node NLTE correction grid.

The RYA-544 decision (corroboration-accept, Ryan 2026-07-13) authorized wiring the Mallinson-2024
Ti departure grid (Zenodo 10753497, ab-initio Grumer-Barklem-2020 H collisions) as the production
Ti Engine-A source, replacing the outdated Bergemann-2011 scaled-Drawin atom. This script derives
the per-line NLTE delta = A(NLTE) - A(LTE) across the benchmark node box via the EXISTING PySME
harness (pipeline.pysme_nlte.nlte_delta), the same path as Na/Mg/Si/N (RYA-409/410/526).

FIREWALL: the grid PRODUCES the correction; nothing is tuned. The solar node MUST reproduce the
RYA-544-derived +0.0506 (ab-initio, corroborated by two blend-aware instruments in RYA-545) or the
script aborts — a guard that the permanent NLTE_LINES['Ti'] entry + grid wiring are self-consistent
with the validated derivation, not silently drifted.

Ti is a first-class PySME element now (permanent NLTE_LINES['Ti'] / _A_SUN / _GRID_FILENAME / _Z in
pipeline/pysme_nlte.py); no monkeypatch. Runs on Sirius only (11.4 GB grid, PySME/MARCS).

  python scripts/rya545_derive_ti_grid.py            # derive + write the CSV + provenance
  python scripts/rya545_derive_ti_grid.py --solar-only   # just re-derive the solar node (firewall check)
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import numpy as np
import pandas as pd

from pipeline import pysme_nlte as P

# Benchmark node box (RYA-409 _NODES): solar + Procyon (~6200/4.0), alpha Cen A/B, 55 Cnc
# (5172/4.43/+0.31), tau Boo, metal-rich to +0.6. All well inside the Mallinson Ti grid hull
# (Teff 2500-8000, logg -0.5..5.5, [Fe/H] -5..+1), so the in-hull guard passes at every node.
_NODES = [(5772, 4.44, 0.0), (5100, 4.44, 0.0), (6200, 4.44, 0.0), (5772, 4.0, 0.0),
          (5772, 4.7, 0.0), (5772, 4.44, -0.3), (5172, 4.43, 0.31), (5100, 4.5, 0.5),
          (6000, 4.3, -0.2), (5772, 4.44, 0.6), (5400, 4.5, 0.5)]

_SOLAR_DELTA_EXPECT = 0.0506       # RYA-544 derived, RYA-545 corroborated (ab-initio)
_SOLAR_TOL = 0.003                 # numerical reproducibility band

_GRID_SRC = Path('/srv/codex/grids/nlte/amarsi_galah/nlte_Ti_pysme.grd')
_OUT = REPO / 'data' / 'nlte_grids' / 'Ti_Mallinson2024_PySME.csv'


def _solar_node():
    r = P.nlte_delta('Ti', star={'teff': 5772, 'logg': 4.44, 'feh': 0.0, 'vmic': 1.0})
    return r


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--solar-only', action='store_true')
    a = ap.parse_args(argv)

    print("=== RYA-545 — derive Ti_Mallinson2024_PySME.csv (ab-initio, firewall: no tuning) ===")
    print(f"  grid = {P._GRID_FILENAME['Ti']}  atmosphere = MARCS (marcs2012, PP)")
    print(f"  Ti lines = {[l[0] for l in P.NLTE_LINES['Ti']]}  A(Ti)_sun = {P._A_SUN['Ti']}")

    # Firewall: the solar node must reproduce the validated +0.0506.
    sol = _solar_node()
    print("\n  solar node per-line delta:")
    for wl, d in sol['per_line'].items():
        print(f"    Ti I {wl:.3f}: {d:+.4f}")
    d_sol = sol['delta_median']
    print(f"  solar MEDIAN delta = {d_sol:+.4f}  (expect {_SOLAR_DELTA_EXPECT:+.4f} +/- {_SOLAR_TOL})")
    if abs(d_sol - _SOLAR_DELTA_EXPECT) > _SOLAR_TOL:
        raise SystemExit(f"FIREWALL ABORT: solar Ti delta {d_sol:+.4f} != validated "
                         f"{_SOLAR_DELTA_EXPECT:+.4f} (tol {_SOLAR_TOL}). Do NOT write the grid.")
    print("  [firewall OK] solar node reproduces the RYA-544/545 ab-initio value.")
    if a.solar_only:
        return

    # Multi-node derivation.
    rows = []
    for (te, lg, fe) in _NODES:
        r = P.nlte_delta('Ti', star={'teff': te, 'logg': lg, 'feh': fe, 'vmic': 1.0})
        for wl, dd in r['per_line'].items():
            rows.append(dict(element='Ti', ion=1, wave_A=round(wl, 3),
                             teff_K=te, logg=lg, feh=fe, delta_nlte=round(float(dd), 4)))
        print(f"  node Teff={te} logg={lg} [Fe/H]={fe:+.2f}: median {r['delta_median']:+.4f} "
              f"({len(r['per_line'])} lines)")

    df = pd.DataFrame(rows)
    df.to_csv(_OUT, index=False)
    feh_max = max(n[2] for n in _NODES)
    print(f"\n  wrote {_OUT.name} ({len(rows)} rows, {len(_NODES)} nodes x {len(P.NLTE_LINES['Ti'])} "
          f"lines; [Fe/H] up to {feh_max:+.1f} — covers 55 Cnc +0.31)")

    # Provenance sidecar (md5-pinned source grid + derivation params).
    md5 = hashlib.md5(_GRID_SRC.read_bytes()).hexdigest() if _GRID_SRC.exists() else 'GRID_ABSENT'
    prov = {
        'grid_csv': _OUT.name,
        'source_grid': str(_GRID_SRC),
        'source_grid_md5': md5,
        'source_doi': '10.5281/zenodo.10753497',
        'reference': 'Mallinson et al. 2024 (A&A 687, A5), ab-initio Grumer & Barklem 2020 H collisions',
        'derivation': 'pipeline.pysme_nlte.nlte_delta (PySME COG, MARCS marcs2012 PP), RYA-545',
        'lines': [l[0] for l in P.NLTE_LINES['Ti']],
        'solar_delta_median': round(d_sol, 4),
        'nodes': _NODES,
        'supersedes': 'Ti_Bergemann2011_MPIA.csv (Bergemann-2011 scaled-Drawin, RYA-546 vintage audit)',
        'firewall': 'derived not tuned; solar node reproduces RYA-544 +0.0506 (RYA-545 corroboration-accept)',
        'ticket': 'RYA-545',
    }
    prov_path = _OUT.with_suffix('.csv.prov.json')
    prov_path.write_text(json.dumps(prov, indent=2))
    print(f"  wrote {prov_path.name} (source md5 {md5[:8]}...)")


if __name__ == '__main__':
    main()
