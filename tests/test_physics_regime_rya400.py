"""
tests/test_physics_regime_rya400.py
===================================
RYA-400 — the per-element physics-regime map is a LIVING artifact, not a static doc.
Guards: full coverage of TARGET_ELEMENTS; every regime call is cited (no regime from
memory); the live-code audit PASSES (every LOCKED element is actually NLTE-wired); and
the audit's central guard BITES — a LOCKED row that is not wired loud-fails (so a green
map can never overstate the physics).
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import audit_physics_regime_rya400 as A            # noqa: E402
from config.constants import TARGET_ELEMENTS, NLTE_CORRECTION_ELEMENTS  # noqa: E402

MAP = A._load_map()


def test_full_coverage_of_target_elements():
    missing = [el for el in TARGET_ELEMENTS if el not in MAP]
    assert not missing, f"unmapped target elements: {missing}"


def test_every_regime_call_is_cited():
    # "no regime asserted from memory" — dimensionality + departure each carry a cite
    for el, spec in MAP.items():
        assert spec['dimensionality'].get('cite'), f"{el} dimensionality missing cite"
        assert spec['departure'].get('cite'), f"{el} departure missing cite"
        assert spec.get('verdict'), f"{el} missing verdict"


def test_registry_locked_requires_measured_or_pending():
    # RYA-428: registration alone is NOT "handled". A registered element is LOCKED only if
    # a measured line of its prescribed ion exists; otherwise it is the legitimate
    # grid-registered/measurement-owed state (verdict GET-DATA), NOT LOCKED.
    import audit_physics_regime_rya400 as AUD
    df = AUD._measured_df()
    for el in NLTE_CORRECTION_ELEMENTS:
        ion = AUD._prescribed_ion(el)
        n = AUD._measured_lines_of_ion(el, ion, df)
        v = MAP[el]['verdict']
        if n > 0:
            assert v == 'LOCKED', f"{el} has {n} measured {ion} line(s) but verdict={v}"
        else:
            assert v != 'LOCKED', \
                f"{el} verdict=LOCKED but ZERO measured {ion} lines — registration != handled (RYA-428)"
    for el in ('Fe', 'C', 'N', 'O'):
        assert MAP[el]['verdict'] == 'LOCKED'


def test_handled_precondition_guard_bites(monkeypatch):
    # RYA-428 central contract: if a registered element is marked LOCKED but no measured line
    # of its prescribed ion exists, the audit FAILS loudly (the paper-done gap Sr exposed).
    import audit_physics_regime_rya400 as AUD
    bad = {el: dict(s) for el, s in AUD._load_map().items()}
    bad['Sr'] = dict(bad['Sr']); bad['Sr']['verdict'] = 'LOCKED'   # registered, 0 Sr II measured
    monkeypatch.setattr(AUD, '_load_map', lambda: bad)
    assert AUD.run() == 1                                          # loud-fail


def test_sr_is_grid_registered_but_pending_not_locked():
    # Sr: grid registered (RYA-421) but measurement owed -> GET-DATA, NOT LOCKED (RYA-428).
    assert 'Sr' in NLTE_CORRECTION_ELEMENTS                        # registration preserved
    assert MAP['Sr']['verdict'] == 'GET-DATA'
    assert MAP['Sr']['verified_in_repo']['grid_file'] == 'Sr_Bergemann2012_INSPECT.csv'


def test_known_lurkers_and_data_gaps():
    # spotlight lurkers flagged; the 0-measured-line elements routed GET-DATA / synthesis
    for el in ('S', 'K', 'Co'):
        assert MAP[el].get('spotlight') is True, f"{el} should be spotlight"
    assert MAP['Co']['verdict'] == 'GET-DATA' and MAP['Sc']['verdict'] == 'GET-DATA'
    assert MAP['Ba']['ion'] == 'II' and MAP['Sr']['ion'] == 'II'   # ion-II elements
    # HARD-carry-forward trio tied to their open bugs
    assert MAP['Li']['verdict'] == 'HARD-carry-forward'
    assert MAP['Eu']['verdict'] == 'HARD-carry-forward'
    assert MAP['Zr']['verdict'] == 'HARD-carry-forward'


def test_live_audit_passes():
    assert A.run() == 0           # every LOCKED is wired, coverage + counts consistent


def test_is_wired_reflects_real_code():
    # RYA-412: Al and S NLTE-wired (RYA-402 PySME). RYA-421: Sr II grid registered + wired.
    for el in ('Fe', 'C', 'O', 'Mg', 'Ca', 'Ti', 'Cr', 'Na', 'Ba', 'Mn', 'Si', 'Al', 'S', 'Sr'):
        assert A._is_wired(el)[0], f"{el} should be wired"
    for el in ('K', 'Co', 'Cu', 'Ni', 'Li', 'Eu', 'Zr', 'V'):
        assert not A._is_wired(el)[0], f"{el} should NOT be wired"


def test_locked_not_wired_guard_bites(monkeypatch):
    # the central contract: if a NOT-wired element is mislabelled LOCKED, the audit FAILS.
    # RYA-412: Al is now genuinely wired, so use Cu (HAVE grid but unregistered/not-wired).
    bad = dict(MAP)
    cu = dict(bad['Cu']); cu['verdict'] = 'LOCKED'        # Cu is not NLTE-wired
    bad['Cu'] = cu
    monkeypatch.setattr(A, '_load_map', lambda: bad)
    assert A.run() == 1


def test_measured_count_drift_guard_bites(monkeypatch):
    bad = dict(MAP)
    fe = dict(bad['Fe'])
    fe['verified_in_repo'] = dict(fe['verified_in_repo'])
    fe['verified_in_repo']['measured_ew_lines'] = 99999    # wrong on purpose
    bad['Fe'] = fe
    monkeypatch.setattr(A, '_load_map', lambda: bad)
    assert A.run() == 1
