#!/usr/bin/env python3
"""
scripts/audit_fe_leg_placeholders_rya417.py
===========================================
RYA-417 Part A — per-line LIVE verification of the placeholder-zero lines in the
*separate* MPIA Bergemann Fe NLTE leg (Fe_Bergemann_MPIA.csv) + the stale Si 6253
(Si_Bergemann_MPIA.csv), the path the RYA-413 registry guard does NOT cover.

A line that is identically 0 across ALL its (teff,logg,feh) nodes in the committed
grid is a *candidate*. The load-bearing step (RYA-413's RCA method) is to RE-QUERY
MPIA live and classify — never drop on the stored zero alone:

  * TRUE_PLACEHOLDER  — MPIA offers the line but returns literal 0.0 at every queried
                        node (served-not-modelled, like Ca 6166). → safe to DROP.
  * GENUINE_NEAR_ZERO — MPIA returns small-but-real nonzero corrections (the line is
                        genuinely NLTE-insensitive, or the grid 0 is a build artifact).
                        → KEEP. Dropping it would delete real physics.
  * NOT_SERVED        — MPIA returns only error codes (10/20/30 → NaN) at every node;
                        it cannot model the line at all (distinct from a served 0).

Live re-query reuses the committed scraper (`_submit_batch`, the same POST the grids
were built with). No tuning, no drop without a live TRUE_PLACEHOLDER confirmation.

    python scripts/audit_fe_leg_placeholders_rya417.py            # live (hits MPIA)
    python scripts/audit_fe_leg_placeholders_rya417.py --offline  # scan only, no query
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / 'scripts'))

from pipeline.nlte_corrections import detect_placeholder_zero_lines  # RYA-413 detector

FE_GRID = _REPO / 'data' / 'nlte_grids' / 'Fe_Bergemann_MPIA.csv'
SI_GRID = _REPO / 'data' / 'nlte_grids' / 'Si_Bergemann_MPIA.csv'
OUT_CSV = _REPO / 'data' / 'audit' / 'fe_leg_placeholders_rya417.csv'

# SIU Z.01 species codes (Fe I 26.01, Fe II 26.02, Si I 14.01).
SPECIES = {('Fe', 'I'): '26.01', ('Fe', 'II'): '26.02', ('Si', 'I'): '14.01'}
EPS = 1e-9
MATCH_TOL = 0.15   # Å — returned-vs-queried wavelength match

# Representative live-query nodes spanning the committed grid box
# (teff 5000-6600 / logg 4.0-4.6 / feh -0.2..+0.35), incl. a near-solar node.
QUERY_NODES = [
    {'name': 'n5800g44f00', 'teff_K': 5800, 'logg': 4.40, 'feh': 0.00},
    {'name': 'n5400g43f00', 'teff_K': 5400, 'logg': 4.30, 'feh': 0.00},
    {'name': 'n6200g45f02', 'teff_K': 6200, 'logg': 4.50, 'feh': 0.20},
    {'name': 'n5000g40f02m', 'teff_K': 5000, 'logg': 4.00, 'feh': -0.20},
    {'name': 'n6600g46f035', 'teff_K': 6600, 'logg': 4.60, 'feh': 0.35},
    {'name': 'n5800g43f00', 'teff_K': 5800, 'logg': 4.30, 'feh': 0.00},
]


def _scan(path: Path, ions: list[str]) -> list[tuple]:
    """Return [(element, ion, wave)] that are 0 across all nodes (per ion if present)."""
    df = pd.read_csv(path)
    out = []
    if 'ion' in df.columns:
        for ion in ions:
            sub = df[df['ion'] == ion]
            for w in detect_placeholder_zero_lines(sub):
                out.append((str(df['element'].iloc[0]), ion, w))
    else:
        for w in detect_placeholder_zero_lines(df):
            out.append((str(df['element'].iloc[0]), ions[0], w))
    return out


def _live_query(lines: list[float], species_code: str) -> dict:
    """{queried_wave: [delta per QUERY_NODE]} via the committed MPIA scraper."""
    from build_nlte_grids_mpia import _submit_batch
    res = _submit_batch(lines, species_code, QUERY_NODES)   # {star: {mpia_wave: delta}}
    mpia_waves = sorted({w for d in res.values() for w in d})
    out = {}
    for qw in lines:
        if not mpia_waves:
            out[qw] = [np.nan] * len(QUERY_NODES)
            continue
        mw = min(mpia_waves, key=lambda m: abs(m - qw))
        if abs(mw - qw) > MATCH_TOL:
            out[qw] = [np.nan] * len(QUERY_NODES)        # MPIA did not offer this wave
            continue
        out[qw] = [res.get(n['name'], {}).get(mw, np.nan) for n in QUERY_NODES]
    return out


def _classify(deltas: list[float]) -> str:
    arr = np.array(deltas, dtype=float)
    finite = arr[np.isfinite(arr)]
    if len(finite) == 0:
        return 'NOT_SERVED'                  # only error codes → MPIA cannot model it
    if (np.abs(finite) < EPS).all():
        return 'TRUE_PLACEHOLDER'            # served, but literal 0 at every node
    return 'GENUINE_NEAR_ZERO'               # small-but-real correction exists → KEEP


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--offline', action='store_true', help='scan only, skip the live MPIA query')
    args = ap.parse_args()

    cand = _scan(FE_GRID, ['I', 'II']) + _scan(SI_GRID, ['I'])
    print(f"RYA-417 Part A — placeholder-zero candidates in the Fe leg + Si grid: {len(cand)}")
    n_fe1 = sum(1 for e, i, _ in cand if e == 'Fe' and i == 'I')
    n_fe2 = sum(1 for e, i, _ in cand if e == 'Fe' and i == 'II')
    n_si = sum(1 for e, _, _ in cand if e == 'Si')
    print(f"  Fe I: {n_fe1}   Fe II: {n_fe2}   Si I: {n_si}")

    if args.offline:
        for e, i, w in cand:
            print(f"  CANDIDATE {e} {i} {w}")
        return 0

    # live re-query, grouped by (element, ion) species code → one POST each
    verdicts = {}
    served = {}
    for (el, ion), code in SPECIES.items():
        waves = [w for e, i, w in cand if e == el and i == ion]
        if not waves:
            continue
        print(f"\n  live MPIA query: {el} {ion} ({code}) — {len(waves)} lines × {len(QUERY_NODES)} nodes …")
        q = _live_query(waves, code)
        for w in waves:
            served[(el, ion, w)] = q[w]
            verdicts[(el, ion, w)] = _classify(q[w])
        time.sleep(1.0)

    # table
    print(f"\n{'='*92}")
    print(f"  {'line':18s} {'served Δ across nodes (live MPIA)':50s} {'verdict'}")
    print(f"  {'-'*88}")
    rows = []
    for (el, ion, w) in sorted(served, key=lambda k: (k[0], k[1], k[2])):
        d = served[(el, ion, w)]
        v = verdicts[(el, ion, w)]
        shown = ' '.join('NaN' if not np.isfinite(x) else f'{x:+.3f}' for x in d)
        print(f"  {el+' '+ion+' '+str(w):18s} {shown:50s} {v}")
        rows.append({'element': el, 'ion': ion, 'wave_A': w, 'verdict': v,
                     'served_deltas': shown,
                     'query_nodes': ';'.join(n['name'] for n in QUERY_NODES)})

    tp = [r for r in rows if r['verdict'] == 'TRUE_PLACEHOLDER']
    gn = [r for r in rows if r['verdict'] == 'GENUINE_NEAR_ZERO']
    ns = [r for r in rows if r['verdict'] == 'NOT_SERVED']
    print(f"\n  SUMMARY: TRUE_PLACEHOLDER={len(tp)} (DROP)  GENUINE_NEAR_ZERO={len(gn)} (KEEP)  "
          f"NOT_SERVED={len(ns)}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    print(f"  wrote provenance → {OUT_CSV.relative_to(_REPO)}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
