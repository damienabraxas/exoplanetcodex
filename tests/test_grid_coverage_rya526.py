"""
tests/test_grid_coverage_rya526.py
==================================
RYA-526 — grid coverage program.

Part A (N wiring): N is registered in NLTE_CORRECTION_ELEMENTS; the applied solar delta is
the VENDORED PySME-derived grid value read through the production interpolation subsystem
(validate-don't-tune — it reproduces the committed RYA-369 load-test, not fitted to Asplund);
the excluded cool-metal-rich weak-line rails fall out-of-hull (loud-flagged, never a spurious
value). NOTE: applying this registered N I NLTE delta in the phase_c KP channel (N NLTE-OWED ->
CURATION-OWED, off the NLTE-OWED bucket — NOT PASS, the +0.36 gf/data residual stays owed) is
RYA-556's scope; here the grid delta is only verified to resolve, verdict wiring untouched.

Part C (stewardship): every grid referenced by NLTE_CORRECTION_ELEMENTS is PRESENT on disk —
the RYA-525 loud-fail guard's "something real to check" (Mn's PASS depends on the COMMITTED
Mn_Bergemann_MPIA.csv, not a missing 8.4 GB file — the ticket's GB was a KB transcription).

Coverage: every canonical TARGET_ELEMENT has a row in the living coverage ledger.

Runs on the COMMITTED grid CSVs — no slow synthesis, no gitignored .grd needed.
"""
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import NLTE_CORRECTION_ELEMENTS, TARGET_ELEMENTS  # noqa: E402
from pipeline import nlte_corrections as N                              # noqa: E402

GRID_DIR = ROOT / 'data' / 'nlte_grids'
LEDGER = ROOT / 'data' / 'curation' / 'nlte_grid_availability.csv'


# ── Part A: N wiring ────────────────────────────────────────────────────────────
def test_N_registered_in_nlte_elements():
    assert 'N' in NLTE_CORRECTION_ELEMENTS
    spec = NLTE_CORRECTION_ELEMENTS['N']
    assert spec['grid'] == 'N_Amarsi2020_PySME.csv'
    assert int(spec['ion']) == 1
    assert spec['flag'] == 'NLTE_Amarsi2020_PySME_1D'
    assert (GRID_DIR / spec['grid']).exists()


def test_N_solar_delta_reproduces_rya369_loadtest():
    """validate-don't-tune: the wired solar deltas ARE the committed RYA-369 load-test values."""
    N._mpia_element_cache.clear()
    expect = {7468.31: -0.0115, 8216.34: -0.0145, 8683.40: -0.0154}
    for w, e in expect.items():
        d = N._mpia_element_delta('N', w, 5772.0, 4.44, 0.0)
        assert np.isfinite(d), f'N {w} solar delta not finite'
        assert abs(d - e) < 1e-3, f'N {w}: {d} != load-test {e}'
        assert d < 0.0                       # N I red is near-LTE, small NEGATIVE (RYA-339 sign)


def test_N_cool_metalrich_rail_is_out_of_hull_not_spurious():
    """The excluded weak-line rails must NOT resolve to a value — they loud-flag (RYA-409)."""
    N._mpia_element_cache.clear()
    d = N._mpia_element_delta('N', 7468.31, 5100.0, 4.5, 0.5)   # an excluded rail node
    assert not np.isfinite(d)                # out-of-hull -> nan -> flagged, never +0.2


def test_N_grid_has_no_railed_values():
    """No |delta| may sit at the +-0.2 bracket edge (a rail sentinel) in the committed grid."""
    N._mpia_element_cache.clear()
    c = N._load_mpia_element_grid('N')
    with (GRID_DIR / 'N_Amarsi2020_PySME.csv').open() as fh:
        deltas = [float(r['delta_nlte']) for r in csv.DictReader(fh)]
    assert deltas
    assert max(abs(d) for d in deltas) < 0.199   # every retained node is a physical value


# ── Part C: stewardship — every wired grid is present on disk ───────────────────
def test_every_wired_grid_is_present_on_disk():
    """The RYA-525 loud-fail guard's real check: a wired PASS depends on a file that EXISTS."""
    missing = []
    for el, spec in NLTE_CORRECTION_ELEMENTS.items():
        grid = spec.get('grid')
        if not grid or grid.endswith('/'):    # CNO-dir style handled elsewhere
            continue
        if not (GRID_DIR / grid).exists():
            missing.append(f'{el} -> {grid}')
    assert not missing, f'wired NLTE grids missing from disk: {missing}'


def test_mn_pass_depends_on_committed_csv_not_missing_giant():
    """Confirm-don't-assume: Mn's grid is the COMMITTED small CSV (ticket's 8.4 GB was KB)."""
    mn = GRID_DIR / 'Mn_Bergemann_MPIA.csv'
    assert mn.exists()
    assert mn.stat().st_size < 1_000_000     # KB-scale committed file, not a GB grid


# ── Coverage: every canonical element is accounted for in the ledger ────────────
def test_every_target_element_in_coverage_ledger():
    with LEDGER.open() as fh:
        rows = list(csv.DictReader(fh))
    elements_in_ledger = {r['element'] for r in rows}
    missing = [el for el in TARGET_ELEMENTS if el not in elements_in_ledger]
    assert not missing, f'canonical elements absent from the coverage ledger: {missing}'
