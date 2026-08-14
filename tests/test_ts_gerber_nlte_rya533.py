"""
tests/test_ts_gerber_nlte_rya533.py
===================================
RYA-533 — THE TEST TELL for the TS-native Gerber NLTE synthesis deck (Engine-B NLTE
for the two-engine floor RYA-525). This is the reproduction test whose ABSENCE was the
decisive Step-0 forensic signal that the deck had never been built; it now exists so the
next capability sweep (RYA-530 class) finds it.

What the deck is: Turbospectrum bsyn NLTE with the Gerber 2023 (A&A 669, A43) departure
grids + TS-native model atom, driven via NLTEINFOFILE / MODELATOMFILE / DEPARTUREFILE,
with 4D (Teff/logg/[Fe/H]/A(X)) departure interpolation. This is TS-native in-synthesis
NLTE — distinct from the PySME / Amarsi-.grd path (RYA-402/529) and the MPIA/Korotin/
INSPECT/nlte_cno 1D-delta paths (Engine A).

THE GATE (validate-don't-tune): reproduce the solar Na correction through the bsyn NLTE
departure path. Anchor -0.107 (INSPECT/Lind 2011 / Amarsi); a ~0.02-0.05 dex model-atom
difference vs the Gerber atom is path-correctness, not a tune.

The live deck runs only on Sirius (needs the compiled Turbospectrum_NLTE binaries + the
15.9 GB Gerber Na departure grid). In CI the live test skips; the recorded-result test
always runs and guards the validated number.
"""
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
# Standalone-script bootstrap (RYA-313): repo root on sys.path BEFORE importing
# config/pipeline, so this runs from any cwd. Derived from __file__, never cwd.
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from config.constants import codex_path  # RYA-810 path register

_REPO = Path(__file__).resolve().parents[1]
_RESULT = _REPO / 'data' / 'audit' / 'rya533_ts_gerber_build' / 'na_gate_result.txt'
_DECK = _REPO / 'scripts' / 'ts_gerber_na_gate_rya533.py'

# Sirius provisioning (the live gate only runs where these exist).
_TS_NLTE_BSYN = Path(str(codex_path('engines.ts_bsyn')))
_GERBER_NA_GRID = codex_path('grids.gerber_ts') / 'NLTEgrid4TS_Na_MARCS_Jul-14-2023.bin'

_ANCHOR = -0.107
_MODEL_ATOM_TOL = 0.05


def _parse_median(text: str) -> float:
    m = re.search(r'median delta\s*=\s*([+-]?\d+\.\d+)', text)
    if not m:
        raise AssertionError(f"no 'median delta' line in gate output:\n{text[-500:]}")
    return float(m.group(1))


def test_ts_gerber_na_gate_recorded():
    """The recorded TS-Gerber Na gate reproduces the anchor within model-atom tolerance.

    Guards the banked result (data/audit/rya533_ts_gerber_build/na_gate_result.txt) so a
    regression in the deck or a bad re-bank is caught, and so this reproduction test is
    grep-discoverable (gerber / ts-native / bsyn nlte / DEPARTUREFILE) by the next sweep."""
    assert _RESULT.exists(), f"missing recorded gate result: {_RESULT}"
    txt = _RESULT.read_text()
    assert 'NLTE departure engaged: True' in txt, "departures must be engaged (no silent LTE)"
    assert 'VERDICT: PASS' in txt
    median = _parse_median(txt)
    assert abs(median - _ANCHOR) <= _MODEL_ATOM_TOL, (
        f"TS-Gerber Na median {median:+.4f} not within {_MODEL_ATOM_TOL} of anchor {_ANCHOR}")


@pytest.mark.skipif(
    not (_TS_NLTE_BSYN.exists() and _GERBER_NA_GRID.exists()),
    reason="live TS-Gerber deck runs only on Sirius (Turbospectrum_NLTE build + 15.9 GB "
           "Gerber Na departure grid); see data/nlte_grids/gerber_ts/Na_gerber2023.prov.json")
def test_ts_gerber_na_gate_live():
    """Run the actual bsyn NLTE deck end-to-end and assert the Na anchor reproduces.

    RAISES inside the deck if departures do not engage (silent-LTE guard)."""
    venv_py = str(codex_path('engines.venv_pysme_python'))
    py = venv_py if os.path.exists(venv_py) else sys.executable
    r = subprocess.run([py, str(_DECK)], capture_output=True, text=True, timeout=1800)
    assert 'NLTE departure engaged: True' in r.stdout, (
        f"departures not engaged:\n{r.stdout[-1500:]}\n{r.stderr[-500:]}")
    median = _parse_median(r.stdout)
    assert abs(median - _ANCHOR) <= _MODEL_ATOM_TOL, (
        f"live TS-Gerber Na median {median:+.4f} not within {_MODEL_ATOM_TOL} of {_ANCHOR}")
