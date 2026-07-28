"""
tests/test_grid_drift_guard_rya418.py
=====================================
RYA-418 — the anti-drift guard. The physics-regime map records each registered element's
NLTE grid in verified_in_repo.grid_file; the LIVE registry (NLTE_CORRECTION_ELEMENTS) is the
single source of truth. The audit's GRID-DRIFT invariant fails loud if the map's grid_file
diverges from the registry (the exact defect RYA-410 left: Na/Mg/Si re-sourced onto Amarsi
PySME while the map still named Bergemann/INSPECT).
"""
import copy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import audit_physics_regime_rya400 as A            # noqa: E402
from config.constants import NLTE_CORRECTION_ELEMENTS  # noqa: E402


def test_map_grid_files_match_registry():
    # after the RYA-418 reconcile, every registered element's map grid_file == registry grid.
    MAP = A._load_map()
    for el, spec in MAP.items():
        if el not in NLTE_CORRECTION_ELEMENTS:
            continue
        gf = (spec.get('verified_in_repo') or {}).get('grid_file')
        if gf is not None:
            assert gf == NLTE_CORRECTION_ELEMENTS[el]['grid'], f"{el} map grid_file drifted"


def test_audit_passes_clean():
    assert A.run() == 0


def test_guard_fires_on_injected_drift(monkeypatch, capsys):
    # diverge a registry grid to ANOTHER EXISTING grid file (so _is_wired still passes and
    # GRID-DRIFT is isolated from LOCKED-NOT-WIRED), then the audit must FAIL loud on it.
    reg = copy.deepcopy(A.NLTE_CORRECTION_ELEMENTS)
    reg['Ca'] = dict(reg['Ca']); reg['Ca']['grid'] = 'Ti_Bergemann2011_MPIA.csv'
    monkeypatch.setattr(A, 'NLTE_CORRECTION_ELEMENTS', reg)
    rc = A.run()
    out = capsys.readouterr().out
    assert rc == 1
    assert 'GRID-DRIFT' in out
    assert 'Ca' in out
