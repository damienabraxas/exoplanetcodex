#!/usr/bin/env python3
"""
scripts/audit_physics_regime_rya400.py
======================================
RYA-400 "The Beast" — make the per-element physics-regime map a LIVING, code-verified
artifact instead of a static doc. Cross-checks config/physics_regime_rya400.yaml against
the actual repository wiring so a regime CALL can never silently drift from the code:

  1. COVERAGE   — every TARGET_ELEMENTS symbol has a regime row (loud-fail on a gap).
  2. HANDLED    — every element verdict-LOCKED is ACTUALLY handled. RYA-428: "handled"
                  requires BOTH (i) a wired grid (NLTE_CORRECTION_ELEMENTS + grid file,
                  OR the Fe leg, OR the C/N/O nlte_cno leg) AND (ii) at least one MEASURED
                  line of the prescribed ion in the star's pool. Grid registration alone
                  is NOT handled (the paper-done gap: a Sr II grid can be registered while
                  only Sr I is measured). LOCKED-without-a-measured-prescribed-ion →
                  loud-fail, naming the element. (Fe/CNO use their own non-EW legs → exempt.)
  3. PENDING-OK — a registered element with a non-LOCKED verdict is the LEGITIMATE
                  grid-registered/measurement-owed state (RYA-428 relaxes the old rule that
                  equated registration with "done"). Recorded as a note, never a false done.
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


def _measured_df():
    import pandas as pd
    if not MEASURED.exists():
        return None
    return pd.read_csv(MEASURED)


def _measured_counts(df=None) -> dict:
    if df is None:
        df = _measured_df()
    if df is None:
        return {}
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


_ION_INT_TO_STR = {1: 'I', 2: 'II', 3: 'III'}


def _prescribed_ion(el: str) -> "str | None":
    """The ion the registered NLTE grid corrects — the ion that must be MEASURED for
    `el` to count as handled (RYA-428). For a registry element, read it from the live
    registration (single source of truth); else fall back to the map's ion field."""
    if el in NLTE_CORRECTION_ELEMENTS:
        return _ION_INT_TO_STR.get(int(NLTE_CORRECTION_ELEMENTS[el]['ion']))
    return None


def _measured_lines_of_ion(el: str, ion: "str | None", df) -> int:
    """Count measured EW lines of `el` in ionisation stage `ion` (ew_mA > 0) in the
    star's canonical measured pool. The pool's single source of truth is sol_ew_results."""
    if df is None or ion is None or 'ion' not in df.columns:
        return 0
    sub = df[(df['element'] == el) & (df['ion'].astype(str) == ion) & (df['ew_mA'] > 0)]
    return int(len(sub))


def run() -> int:
    emap = _load_map()
    meas_df = _measured_df()
    meas = _measured_counts(meas_df)
    errors, warns, rows, false_handled = [], [], [], []

    # 1. coverage
    for el in TARGET_ELEMENTS:
        if el not in emap:
            errors.append(f"COVERAGE: target element {el!r} has no regime row")

    for el, spec in emap.items():
        verdict = str(spec.get('verdict', '?'))
        wired, how = _is_wired(el)
        claim = spec.get('verified_in_repo', {}) or {}

        # 2. LOCKED ("handled") requires BOTH a wired grid AND a MEASURED line of the
        # prescribed ion in the pool — registration alone is not "handled" (RYA-428). This
        # closes the paper-done gap: a grid can be registered without the abundance existing
        # (e.g. Sr II grid registered while only Sr I 6617 is measured). The generic-registry
        # elements are bound here; Fe (gate) and C/N/O (nlte_cno synthesis) use their own legs
        # and do not draw the prescribed ion from the EW pool, so they are exempt.
        if verdict in LOCKED_VERDICTS and not wired:
            errors.append(f"LOCKED-NOT-WIRED: {el} verdict=LOCKED but not NLTE-wired ({how})")
        if verdict in LOCKED_VERDICTS and wired and el in NLTE_CORRECTION_ELEMENTS:
            ion = _prescribed_ion(el)
            n_ion = _measured_lines_of_ion(el, ion, meas_df)
            if n_ion == 0:
                msg = (f"LOCKED-NOT-MEASURED: {el} verdict=LOCKED (grid {NLTE_CORRECTION_ELEMENTS[el].get('grid')} "
                       f"registered) but ZERO measured {el} {ion} lines in the pool — grid registration "
                       f"is not measurement; revert to GET-DATA-pending until {el} {ion} is measured (RYA-428)")
                errors.append(msg)
                false_handled.append((el, ion))
        # 3. A registered element with a non-LOCKED verdict is the LEGITIMATE
        # grid-registered/measurement-owed state (RYA-428 relaxes the old FALSE-DONE rule,
        # which equated registration with "done"). Recorded as a note, not an error.
        if verdict not in LOCKED_VERDICTS and el in NLTE_CORRECTION_ELEMENTS:
            ion = _prescribed_ion(el)
            n_ion = _measured_lines_of_ion(el, ion, meas_df)
            if n_ion == 0:
                warns.append(f"{el}: grid registered ({NLTE_CORRECTION_ELEMENTS[el].get('grid')}) but "
                             f"verdict={verdict} -- measurement of {el} {ion} OWED (correct pending state)")
            else:
                warns.append(f"{el}: registered + {n_ion} measured {el} {ion} line(s) but verdict={verdict} "
                             f"(under-claim? could be promotable to LOCKED) -- review")

        # 4. measured-line count consistency
        live = int(meas.get(el, 0))
        claimed = claim.get('measured_ew_lines')
        if claimed is not None and int(claimed) != live:
            errors.append(f"MEAS-DRIFT: {el} map says {claimed} measured lines, live CSV has {live}")
        if live == 0 and verdict in LOCKED_VERDICTS and el not in _CNO_LEG:
            warns.append(f"{el}: LOCKED but 0 measured EW lines (synthesis/other path?)")

        # 5. GRID-DRIFT (RYA-418): the map's recorded grid_file MUST equal the live registry
        # grid. The live registration (NLTE_CORRECTION_ELEMENTS) is the single source of truth;
        # the map conforms to it, never the reverse. A silent divergence (e.g. RYA-410 re-sourced
        # Na/Mg/Si onto Amarsi PySME while the map still named Bergemann/INSPECT) is the defect.
        if el in NLTE_CORRECTION_ELEMENTS:
            reg_grid = NLTE_CORRECTION_ELEMENTS[el].get('grid')
            map_grid = claim.get('grid_file')
            if map_grid is not None and map_grid != reg_grid:
                errors.append(f"GRID-DRIFT: {el} map grid_file={map_grid!r} != live registry "
                              f"grid={reg_grid!r} (RYA-418: registry is source of truth; conform the map)")

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

    # RYA-428 RCA sweep: registry elements marked LOCKED with no measured prescribed-ion line
    print("\n  RYA-428 handled-precondition sweep (LOCKED requires a measured prescribed-ion line):")
    if false_handled:
        for el, ion in false_handled:
            print(f"    ✗ FALSE-HANDLED: {el} (needs measured {el} {ion}) — see errors")
    else:
        print("    ✓ no registry element is marked LOCKED without a measured prescribed-ion line")

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
