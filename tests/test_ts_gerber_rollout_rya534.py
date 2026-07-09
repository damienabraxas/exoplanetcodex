"""
tests/test_ts_gerber_rollout_rya534.py
======================================
RYA-534 — Family-A TS-Gerber NLTE rollout guard. Auto-discovers every landed
element's recorded gate result (data/audit/rya534_family_a_rollout/<El>_gate_result.txt)
and asserts it reproduced its anchor with departures engaged (no silent LTE). Extends
the RYA-533 Na "test tell" to the rolled-out Family-A elements; grows automatically as
each element lands. The live per-element deck runs on Sirius (scripts/ts_gerber_gate.py).
"""
import re
from pathlib import Path

import pytest

_ROLLOUT = Path(__file__).resolve().parents[1] / 'data' / 'audit' / 'rya534_family_a_rollout'
_RESULTS = sorted(_ROLLOUT.glob('*_gate_result.txt')) if _ROLLOUT.exists() else []


# Elements the deck provisions + validates-runs (departures engaged, consistent lines) but
# whose TS-Gerber value diverges from our registered Engine-A anchor by a documented
# model-atom difference — a real RYA-525 cross-engine diagnostic, NOT a deck failure, so
# they are allowed to be CHECK rather than PASS. Ti: +0.221 vs Bergemann-2011 +0.107 (~2x).
_MODEL_ATOM_DIFF = {'Ti'}


@pytest.mark.parametrize('result_file', _RESULTS, ids=[p.stem for p in _RESULTS])
def test_family_a_gate_recorded(result_file):
    """Each landed Family-A element: the deck ran (departures engaged) and either reproduced
    its anchor (PASS) or is a documented model-atom-difference CHECK (never a silent LTE)."""
    el = result_file.stem.split('_')[0]
    txt = result_file.read_text()
    # the deck must have engaged departures for every element (no silent LTE)
    assert 'NLTE departure engaged: True' in txt, f"{result_file.name}: departures not engaged"
    med = re.search(r'median delta\s*=\s*([+-]?\d+\.\d+)', txt)
    anc = re.search(r'anchor\s*([+-]?\d+\.\d+)\s*\+/-\s*([\d.]+)', txt)
    assert med and anc, f"{result_file.name}: could not parse median/anchor"
    median, anchor, tol = float(med.group(1)), float(anc.group(1)), float(anc.group(2))
    if el in _MODEL_ATOM_DIFF:
        # documented cross-engine model-atom difference: must still have RUN + be recorded CHECK
        assert 'VERDICT: CHECK' in txt, f"{el}: expected documented CHECK (model-atom diff)"
        assert abs(median - anchor) > tol, f"{el}: recorded as model-atom-diff but is within tol"
    else:
        assert 'VERDICT: PASS' in txt, f"{result_file.name}: gate not PASS"
        assert abs(median - anchor) <= tol, (
            f"{result_file.name}: median {median:+.4f} not within {tol} of anchor {anchor:+.4f}")


def test_rollout_has_landed_elements():
    """Sanity: at least O has landed (guards against an empty/relocated rollout dir)."""
    assert (_ROLLOUT / 'O_gate_result.txt').exists(), "O gate result missing"
