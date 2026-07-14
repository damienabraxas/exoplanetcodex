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


# Elements that must NOT contribute to a green board here, kept as loud strict-xfail.
# Ti: the RCA is now CLOSED (RYA-542/544/546) — the +0.22 this Engine-B TS-Gerber deck
# produces is the outdated Bergemann-2011 SCALED-DRAWIN atom (`atom.ti503b`), inflated ~2x
# over the correct ab-initio ~+0.05 (Grumer-Barklem-2020). Ti's PRODUCTION NLTE is now
# registered on the Engine-A Mallinson-2024 grid (+0.0506, RYA-545, corroboration-accept).
# But THIS gate exercises the Engine-B TS-Gerber deck, which STILL ships the outdated atom,
# so it honestly CHECKs — asserting a PASS would falsely certify the TS-Gerber Ti deck.
# Stays strict-xfail (Engine-B Ti atom is superseded, not validated) until a follow-on swaps
# the TS-Gerber Ti atom to an ab-initio one — a firewall-clean flag, not a silent fallback.
_OPEN_RCA_REASON = {
    'Ti': "Ti Engine-B TS-Gerber leg still on the outdated atom.ti503b (Bergemann-2011 "
          "scaled-Drawin) -> +0.22, KNOWN-inflated ~2x (RYA-546). Production Ti NLTE is "
          "registered on Engine-A Mallinson-2024 ab-initio (+0.0506, RYA-545); this TS-Gerber "
          "deck is superseded, not validated. Not a pass until the Engine-B Ti atom is swapped (RYA-548).",
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
