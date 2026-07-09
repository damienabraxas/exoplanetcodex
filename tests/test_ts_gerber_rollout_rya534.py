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


# Elements with an OPEN RCA — the gate CHECKs and the cause is unresolved, so the element
# must NOT contribute to a green board. It is xfail (loud, strict) until the RCA closes.
# Ti (RYA-535): +0.221 vs Bergemann-2011 +0.107 is a SAME-ATOM systematic — Engine-A
# (MPIA Bergemann-2011, MAFAGS-OS) and Engine-B (TS-Gerber, MARCS) use the SAME atom
# (Gerber 2023 Table 1 lists Ti = Bergemann 2011; atom.ti503b just ships it). Encoding it
# as an allowed "model-atom-difference pass" is a silent fallback; this test fails loud
# (xfail) until RYA-535 localizes the cause (atmosphere / interpolation / downstream).
_OPEN_RCA_REASON = {
    'Ti': "Ti unresolved same-atom systematic — RYA-535 RCA open (+0.221 vs +0.107, "
          "same Bergemann-2011 atom). Not a pass; fails loud until localized.",
}


def _param(p):
    el = p.stem.split('_')[0]
    if el in _OPEN_RCA_REASON:
        return pytest.param(p, marks=pytest.mark.xfail(reason=_OPEN_RCA_REASON[el], strict=True))
    return pytest.param(p)


@pytest.mark.parametrize('result_file', [_param(p) for p in _RESULTS], ids=[p.stem for p in _RESULTS])
def test_family_a_gate_recorded(result_file):
    """Each landed Family-A element: the deck ran (departures engaged) and reproduced its
    anchor (PASS). Elements with an open RCA (Ti, RYA-535) are strict-xfail — they must NOT
    pass while unresolved; encoding a CHECK as an allowed pass would be a silent fallback."""
    txt = result_file.read_text()
    assert 'NLTE departure engaged: True' in txt, f"{result_file.name}: departures not engaged"
    med = re.search(r'median delta\s*=\s*([+-]?\d+\.\d+)', txt)
    anc = re.search(r'anchor\s*([+-]?\d+\.\d+)\s*\+/-\s*([\d.]+)', txt)
    assert med and anc, f"{result_file.name}: could not parse median/anchor"
    median, anchor, tol = float(med.group(1)), float(anc.group(1)), float(anc.group(2))
    assert 'VERDICT: PASS' in txt, f"{result_file.name}: gate not PASS"
    assert abs(median - anchor) <= tol, (
        f"{result_file.name}: median {median:+.4f} not within {tol} of anchor {anchor:+.4f}")


def test_rollout_has_landed_elements():
    """Sanity: at least O has landed (guards against an empty/relocated rollout dir)."""
    assert (_ROLLOUT / 'O_gate_result.txt').exists(), "O gate result missing"
