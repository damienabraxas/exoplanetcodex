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
    # RYA-540 refined this: "has measured lines" does NOT imply LOCKED. Cu has 2 measured
    # Cu I lines AND a registered grid, yet stays GET-DATA because those lines fail
    # QUALITY (5 lines at red_chi2 8-127, RYA-395/466). Measured-but-quality-blocked is a
    # third legitimate state, distinct from both LOCKED and no-data.
    #
    # The strong direction is kept intact and un-weakened:
    #     LOCKED  =>  n > 0        (registration alone is never "handled", RYA-428)
    # The converse now carries an ADJUDICATED exception list. Anything that acquires
    # measured lines without being LOCKED and without an entry here still FAILS — the
    # exception must be argued per element, not assumed.
    QUALITY_BLOCKED = {
        'Cu': 'measured lines fail quality (red_chi2 8-127); GET-DATA per RYA-395/466/540',
    }
    import audit_physics_regime_rya400 as AUD
    df = AUD._measured_df()
    for el in NLTE_CORRECTION_ELEMENTS:
        ion = AUD._prescribed_ion(el)
        n = AUD._measured_lines_of_ion(el, ion, df)
        v = MAP[el]['verdict']
        if n > 0:
            assert v == 'LOCKED' or el in QUALITY_BLOCKED, \
                f"{el} has {n} measured {ion} line(s) but verdict={v} with no adjudicated reason"
        else:
            assert v != 'LOCKED', \
                f"{el} verdict=LOCKED but ZERO measured {ion} lines — registration != handled (RYA-428)"
    # the exception list may not rot into a blanket excuse: every entry must still be
    # non-LOCKED and still have the measured lines that made it interesting.
    for el, why in QUALITY_BLOCKED.items():
        assert MAP[el]['verdict'] != 'LOCKED', f"{el} is LOCKED now — drop it from QUALITY_BLOCKED ({why})"
        assert AUD._measured_lines_of_ion(el, AUD._prescribed_ion(el), df) > 0, \
            f"{el} no longer has measured lines — drop it from QUALITY_BLOCKED ({why})"
    # Fe + the CNO species nlte_cno actually models (C, O) are LOCKED — verified on their legs.
    for el in ('Fe', 'C', 'O'):
        assert MAP[el]['verdict'] == 'LOCKED'
    # RYA-434: N is NOT an nlte_cno species (C/O only) and has no other NLTE path -> demoted
    # from the premature LOCKED to GET-DATA (the flagship CNO paper-done hole this ticket closed).
    assert MAP['N']['verdict'] == 'GET-DATA'


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
    # RYA-412: Al/S NLTE-wired (RYA-402). RYA-421: Sr registered. RYA-462: K registered
    # (grid now applied; the EW-pool measurement is still owed -> GET-DATA-pending, not LOCKED).
    # RYA-434: the CNO leg is read from nlte_cno.REQUIRED_LINES (C/O only).
    #
    # SUPERSEDED since: N and Cu both moved onto the wired side.
    #   N  — RYA-526 registered N_Amarsi2020_PySME.csv in NLTE_CORRECTION_ELEMENTS and
    #        RYA-556 wired it into the phase_c Kitt Peak channel, clearing N off NLTE-OWED.
    #        (N is wired via the registry, not via nlte_cno.REQUIRED_LINES — the CNO
    #        single-source contract below is unchanged.)
    #   Cu — RYA-540 derived + registered Cu_Caliskan2024_PySME.csv.
    # Being wired is NOT being finished: both remain CURATION-OWED / GET-DATA on
    # data-channel grounds, which is what the verdict tests above assert.
    for el in ('Fe', 'C', 'O', 'Mg', 'Ca', 'Ti', 'Cr', 'Na', 'Ba', 'Mn', 'Si', 'Al', 'S',
               'Sr', 'K', 'N', 'Cu'):
        assert A._is_wired(el)[0], f"{el} should be wired"
    for el in ('Co', 'Ni', 'Li', 'Eu', 'Zr', 'V'):
        assert not A._is_wired(el)[0], f"{el} should NOT be wired"


def test_cno_leg_read_from_nlte_cno_single_source():
    # the CNO class is the species nlte_cno models, read from its own REQUIRED_LINES -- not a
    # hardcoded list. C and O are in it; N is not (Amarsi 2019 is a C/O grid).
    assert A._cno_species() == {'C', 'O'}
    assert 'N' not in A._cno_species()


def test_cno_guard_bites_no_nlte_cno_output(monkeypatch):
    # RYA-434 (a): if nlte_cno emits no finite delta, the CNO-class verifier finds no artifact
    # and the audit FAILS -- O (the flagship) LOCKED-without-nlte_cno-output -> exit 1.
    import audit_physics_regime_rya400 as AUD
    from pipeline import nlte_cno
    monkeypatch.setattr(nlte_cno, 'cno_nlte_delta', lambda *a, **k: float('nan'))
    assert AUD.run() == 1


def test_fe_guard_bites_no_converged_solve(monkeypatch):
    # RYA-434 (b): if the Fe parameter-solve artifact (the ratified ionization arbiter) is
    # absent, the Fe-class verifier finds no A(Fe) and the audit FAILS -> exit 1.
    import audit_physics_regime_rya400 as AUD
    import config.constants as C
    monkeypatch.setattr(C, 'FE_IONIZATION_SYNTH_ARBITER', {})   # no converged Fe solve
    assert AUD.run() == 1


def test_locked_not_wired_guard_bites(monkeypatch):
    # the central contract: if a NOT-wired element is mislabelled LOCKED, the audit FAILS.
    # RYA-412: Al wired; RYA-530: Cu wired (Caliskan2024) — the old Cu example went stale
    # (it only "passed" because the pre-fix Ti GRID-DRIFT was independently failing the audit,
    # RYA-552). Use V, the genuine NLTE_VOID (no grid, never wired), so this guard cannot
    # silently go stale again.
    bad = dict(MAP)
    v = dict(bad['V']); v['verdict'] = 'LOCKED'           # V is a genuine NLTE_VOID, not NLTE-wired
    bad['V'] = v
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
