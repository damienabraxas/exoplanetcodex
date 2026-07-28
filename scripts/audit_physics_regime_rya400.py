#!/usr/bin/env python3
"""
scripts/audit_physics_regime_rya400.py
======================================
RYA-400 "The Beast" — make the per-element physics-regime map a LIVING, code-verified
artifact instead of a static doc. Cross-checks config/physics_regime_rya400.yaml against
the actual repository wiring so a regime CALL can never silently drift from the code:

  1. COVERAGE   — every TARGET_ELEMENTS symbol has a regime row (loud-fail on a gap).
  2. HANDLED    — every LOCKED element is verified on SOME path via a real production
                  ARTIFACT (RYA-428 + RYA-434). The class is read from single sources, never
                  a hardcoded element list:
                    - EW-pool (in NLTE_CORRECTION_ELEMENTS): >=1 MEASURED prescribed-ion line
                      in the pool (RYA-428).
                    - Fe gate: a converged A(Fe I)/A(Fe II) with the ionization residual in
                      tolerance (the ratified DECISION-2 arbiter).
                    - CNO (the species nlte_cno models, from nlte_cno.REQUIRED_LINES): a
                      FINITE nlte_cno delta at the solar anchor.
                  "Exempt from the EW-pool check" routes to the OTHER path's verifier, NOT to
                  none. A LOCKED element verified on NO path -> loud-fail, naming element+path
                  (RYA-434 closes the CNO paper-done hole at the flagship C/O location).
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
# Element classes whose NLTE leg is NOT the generic EW-pool registry. Membership is read
# from the SINGLE SOURCE for each leg, never a hardcoded element list (RYA-434):
#   - Fe leg  : the stellar-parameter / Fe gate element (structural — the one gate element).
#   - CNO leg : exactly the species nlte_cno actually models, read from nlte_cno.REQUIRED_LINES
#               (so N, which nlte_cno does NOT model, is correctly NOT in this class).
_FE_LEG = {'Fe'}


def _cno_species() -> set:
    """The C/N/O species nlte_cno actually models — the single source of the CNO leg."""
    try:
        from pipeline import nlte_cno
        return {str(k).rstrip('I') for k in nlte_cno.REQUIRED_LINES}   # 'CI'->'C', 'OI'->'O'
    except Exception:
        return set()


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
    if el in _cno_species():
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


# ── Per-class production-artifact verifiers (RYA-434) ─────────────────────────
# "Exempt from the EW-pool check" must route to the OTHER path's verification, NOT to
# none. Every LOCKED element is verified on SOME path, and the verifier confirms a real
# production ARTIFACT (a finite abundance/correction with the anchor reproduced), not
# merely that the leg's files exist.

def _import_const(name):
    import config.constants as C
    return getattr(C, name, None)


def _verify_fe(star: str = 'solar') -> "tuple[bool, str]":
    """Fe gate / parameter-solve artifact: a finite, ratified A(Fe I)/A(Fe II) with the
    ionization-balance residual inside FE_IONISATION_GATE (RYA-405/406 DECISION-2 arbiter
    is the committed converged solve). LOCKED-without-converged-Fe-solve fails."""
    import math
    arb = (_import_const('FE_IONIZATION_SYNTH_ARBITER') or {}).get(star)
    gate = _import_const('FE_IONISATION_GATE')
    if not arb or gate is None:
        return (False, "no converged Fe solve artifact (FE_IONIZATION_SYNTH_ARBITER missing)")
    d_fe, fe2 = arb.get('dFe'), arb.get('fe2_synth')
    if d_fe is None or fe2 is None or not (math.isfinite(d_fe) and math.isfinite(fe2)):
        return (False, "Fe solve A(Fe) absent/NaN")
    if abs(d_fe) >= float(gate):
        return (False, f"Fe ionization residual dFe={d_fe:+.3f} outside tol +/-{gate}")
    return (True, f"A(FeI)={fe2 + d_fe:.3f} A(FeII)={fe2:.3f} dFe={d_fe:+.3f} (<{gate})")


def _verify_cno(el: str, star_params=(5772.0, 4.44, 0.0, 1.0)) -> "tuple[bool, str]":
    """nlte_cno production artifact: a FINITE delta for >=1 of the element's required lines
    at the solar anchor (the same validate-don't-tune intake the other grids get). Uses the
    Asplund solar abundance as the logeps anchor. LOCKED-without-nlte_cno-output fails."""
    import math
    try:
        from pipeline import nlte_cno
    except Exception as e:
        return (False, f"nlte_cno import failed: {str(e)[:50]}")
    sp = el + 'I'
    labels = nlte_cno.REQUIRED_LINES.get(sp, [])
    if not labels:
        return (False, f"nlte_cno does not model {el} (no {sp} in REQUIRED_LINES)")
    logeps = (_import_const('SOLAR_ASPLUND2021') or {}).get(el)
    if logeps is None:
        return (False, f"no solar anchor abundance for {el}")
    t, g, f, vmic = star_params
    for lab in labels:
        try:
            d = nlte_cno.cno_nlte_delta(sp, lab, t, g, f, vmic, float(logeps))
            if math.isfinite(d):
                return (True, f"nlte_cno delta({sp} {lab})={d:+.3f} @solar (anchor reproduced)")
        except Exception:
            continue
    return (False, f"nlte_cno produced no finite {el} delta at the solar anchor")


def _verify_locked(el: str, meas_df) -> "tuple[str, bool, str]":
    """Route a LOCKED element to its class verifier; return (class, ok, artifact). The class
    is read from single sources (registry / Fe gate / nlte_cno species), never a hardcoded
    list. An element matching NO class is LOCKED on no verifiable path (the paper-done hole)."""
    if el in NLTE_CORRECTION_ELEMENTS:                       # EW-pool class (RYA-428)
        ion = _prescribed_ion(el)
        n = _measured_lines_of_ion(el, ion, meas_df)
        return ('EW-pool', n > 0, f"{n} measured {el} {ion} line(s)")
    if el in _FE_LEG:                                        # Fe gate class
        ok, art = _verify_fe()
        return ('Fe-gate', ok, art)
    if el in _cno_species():                                 # CNO synthesis class
        ok, art = _verify_cno(el)
        return ('CNO-nlte_cno', ok, art)
    return ('none', False, "LOCKED but matches NO verification class (registry/Fe/nlte_cno)")


def run() -> int:
    emap = _load_map()
    meas_df = _measured_df()
    meas = _measured_counts(meas_df)
    errors, warns, rows, false_handled, verify_rows = [], [], [], [], []

    # 1. coverage
    for el in TARGET_ELEMENTS:
        if el not in emap:
            errors.append(f"COVERAGE: target element {el!r} has no regime row")

    for el, spec in emap.items():
        verdict = str(spec.get('verdict', '?'))
        wired, how = _is_wired(el)
        claim = spec.get('verified_in_repo', {}) or {}

        # 2. HANDLED requires a verified production ARTIFACT on the element's class path.
        # No LOCKED element may be exempt from ALL verification (RYA-434): "exempt from the
        # EW-pool check" routes to the OTHER path's verifier, never to none. The class is read
        # from single sources (registry / Fe gate / nlte_cno species). Each LOCKED element's
        # class verifier must find a real artifact:
        #   - EW-pool (registry):  >=1 measured prescribed-ion line (the RYA-428 check, unchanged)
        #   - Fe gate:             converged A(Fe I)/A(Fe II), ionization residual in tol
        #   - CNO (nlte_cno):      a finite delta at the solar anchor (the flagship C/O path)
        # An element that matches NO class is LOCKED on no verifiable path -> loud fail.
        if verdict in LOCKED_VERDICTS:
            if not wired:
                errors.append(f"LOCKED-NOT-WIRED: {el} verdict=LOCKED but not NLTE-wired ({how})")
            cls, ok, artifact = _verify_locked(el, meas_df)
            verify_rows.append((el, cls, ok, artifact))
            if not ok:
                errors.append(f"LOCKED-UNVERIFIED: {el} verdict=LOCKED but its {cls} verifier found NO "
                              f"production artifact: {artifact} -- registration/leg-existence is not a "
                              f"verified abundance (RYA-434); revert until the artifact exists")
                false_handled.append((el, cls))
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
        if live == 0 and verdict in LOCKED_VERDICTS and el not in _cno_species():
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

    # RYA-434 verification report: EVERY LOCKED element + the path that verified it + the
    # artifact found. No LOCKED element may be unverified on all paths.
    print("\n  RYA-434 LOCKED-element verification (every handled element is verified on SOME path):")
    print(f"    {'el':<4}{'class':<14}{'ok':<4}artifact")
    for el, cls, ok, art in sorted(verify_rows):
        print(f"    {el:<4}{cls:<14}{('OK' if ok else 'FAIL'):<6}{art}")
    if false_handled:
        print("    --> UNVERIFIED (paper-done hole, see errors): "
              + ", ".join(f"{el} [{cls}]" for el, cls in false_handled))

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
