#!/usr/bin/env python3
"""
scripts/audit_physics_regime_rya400.py
======================================
RYA-400 "The Beast" — make the per-element physics-regime map a LIVING, code-verified
artifact instead of a static doc. Cross-checks config/physics_regime_rya400.yaml against
the actual repository wiring so a regime CALL can never silently drift from the code:

  1. COVERAGE   — every TARGET_ELEMENTS symbol has a regime row (loud-fail on a gap).
  2. LOCKED     — every element verdict-LOCKED is ACTUALLY wired: in
                  NLTE_CORRECTION_ELEMENTS with its grid file present, OR the Fe leg
                  (apply_fe_nlte_corrections + Fe grid), OR the C/N/O nlte_cno leg.
                  A LOCKED row that is not wired → loud-fail (the whole point: a green
                  matrix means the physics it claims is really applied).
  3. NOT-WIRED  — a non-LOCKED verdict (GET-GRID/GET-DATA/HARD/LTE-OK) must NOT be
                  registered as an applied NLTE element (consistency; no false "done").
  4. MEAS-LINES — each row's `measured_ew_lines` matches the live count in
                  data/measured/sol_ew_results_v1.csv (catches data drift; 0-line
                  elements like Co/K/Sc/N are flagged for the GET-DATA routing).

Exit 0 = matrix consistent with code. Exit 1 = a LOCKED claim is unwired, a target is
unmapped, or a measured-line count drifted. Validate-don't-tune: this NEVER changes an
abundance — it only asserts each element got the treatment the map claims.

Usage:  python scripts/audit_physics_regime_rya400.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import TARGET_ELEMENTS, NLTE_CORRECTION_ELEMENTS  # noqa: E402

REGIME_YAML = ROOT / 'config' / 'physics_regime_rya400.yaml'
GRID_DIR = ROOT / 'data' / 'nlte_grids'
MEASURED = ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv'

LOCKED_VERDICTS = {'LOCKED'}
# elements whose NLTE leg is NOT the generic registry (handled by dedicated modules)
_FE_LEG = {'Fe'}
_CNO_LEG = {'C', 'N', 'O'}


def _load_map() -> dict:
    return yaml.safe_load(REGIME_YAML.read_text())['elements']


def _measured_counts() -> dict:
    import pandas as pd
    if not MEASURED.exists():
        return {}
    df = pd.read_csv(MEASURED)
    return df.groupby('element').size().to_dict()


def _grid_present(spec_grid) -> bool:
    return bool(spec_grid) and (GRID_DIR / str(spec_grid)).exists()


def _is_wired(el: str) -> "tuple[bool, str]":
    """Is `el` actually NLTE-wired in the live code? Returns (wired, how)."""
    if el in NLTE_CORRECTION_ELEMENTS:
        g = NLTE_CORRECTION_ELEMENTS[el].get('grid')
        return (_grid_present(g), f"registry[{el}] grid={g}")
    if el in _FE_LEG:
        return (_grid_present('Fe_Bergemann_MPIA.csv'), "Fe leg (apply_fe_nlte_corrections)")
    if el in _CNO_LEG:
        return ((GRID_DIR / 'amarsi2019_cno').is_dir(), "nlte_cno (amarsi2019_cno)")
    return (False, "not wired")


def run() -> int:
    emap = _load_map()
    meas = _measured_counts()
    errors, warns, rows = [], [], []

    # 1. coverage
    for el in TARGET_ELEMENTS:
        if el not in emap:
            errors.append(f"COVERAGE: target element {el!r} has no regime row")

    for el, spec in emap.items():
        verdict = str(spec.get('verdict', '?'))
        wired, how = _is_wired(el)
        claim = spec.get('verified_in_repo', {}) or {}

        # 2. LOCKED must be wired
        if verdict in LOCKED_VERDICTS and not wired:
            errors.append(f"LOCKED-NOT-WIRED: {el} verdict=LOCKED but not NLTE-wired ({how})")
        # 3. non-LOCKED must not be registered as applied
        if verdict not in LOCKED_VERDICTS and el in NLTE_CORRECTION_ELEMENTS:
            errors.append(f"FALSE-DONE: {el} verdict={verdict} but is in NLTE_CORRECTION_ELEMENTS")

        # 4. measured-line count consistency
        live = int(meas.get(el, 0))
        claimed = claim.get('measured_ew_lines')
        if claimed is not None and int(claimed) != live:
            errors.append(f"MEAS-DRIFT: {el} map says {claimed} measured lines, live CSV has {live}")
        if live == 0 and verdict in LOCKED_VERDICTS and el not in _CNO_LEG:
            warns.append(f"{el}: LOCKED but 0 measured EW lines (synthesis/other path?)")

        rows.append((el, spec.get('ion', '?'), verdict, 'wired' if wired else '-', live,
                     str(spec.get('routing', '') or '')))

    # report
    print("=" * 96)
    print("  RYA-400 — per-element physics-regime map ↔ live-code audit")
    print("=" * 96)
    print(f"  {'el':<4}{'ion':<7}{'verdict':<20}{'nlte':<7}{'meas':>5}  routing")
    print("  " + "-" * 92)
    order = {'LOCKED': 0, 'GET-GRID': 1, 'GET-DATA': 2, 'GET-3D': 3,
             'HARD-carry-forward': 4, 'LTE-OK': 5}
    for el, ion, verdict, w, live, routing in sorted(rows, key=lambda r: (order.get(r[2], 9), r[0])):
        print(f"  {el:<4}{str(ion):<7}{verdict:<20}{w:<7}{live:>5}  {routing[:46]}")

    # routing tally
    from collections import Counter
    tally = Counter(spec.get('verdict', '?') for spec in emap.values())
    print("  " + "-" * 92)
    print("  verdict tally: " + "  ".join(f"{k}={v}" for k, v in sorted(tally.items())))
    print(f"  elements mapped: {len(emap)} (TARGET_ELEMENTS: {len(TARGET_ELEMENTS)}; Fe I/II → 27 species)")

    if warns:
        print("\n  NOTES:")
        for w in warns:
            print(f"    • {w}")
    print("\n" + "-" * 96)
    if errors:
        print(f"  RESULT: FAIL — {len(errors)} inconsistency(ies) between the map and the code:")
        for e in errors:
            print(f"    ✗ {e}")
        print("=" * 96)
        return 1
    print("  RESULT: PASS — every LOCKED regime is wired; coverage complete; measured-line "
          "counts match. The map is consistent with the code.")
    print("=" * 96)
    return 0


if __name__ == '__main__':
    sys.exit(run())
