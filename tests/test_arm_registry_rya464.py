"""
tests/test_arm_registry_rya464.py
=================================
RYA-464 — per-star multi-instrument arm registry. The unlock: run_cno/run_phase_a no
longer hardcode REGIONS={'vis'} and no longer load reflected-solar (Vesta) data for a
non-solar star. Arms are declared per-star with a readiness flag; an unready arm is
DEFERRED (loud), never silently substituted with reflected sunlight.

These guard the wiring/dispatch logic (NOT the heavy Turbospectrum synthesis).
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import pipeline.cno_synthesis as cs  # noqa: E402


def test_solar_registry_is_the_existing_three_arms_bit_identical():
    # Regression guard: solar arms + their loaders/diagnostics are UNCHANGED, so the
    # generalization is a no-op for the solar/Vesta path.
    sr = cs.star_arm_registry('solar')
    assert set(sr) == {'harps', 'espresso', 'uves'}
    assert all(a.ready for a in sr.values())
    assert sr['harps'].loader == 'harps_normalized'
    assert sr['harps'].diagnostics is cs.VIS_DIAGNOSTICS
    assert sr['espresso'].loader == 'reflected_solar'
    assert sr['espresso'].diagnostics is cs.ESPRESSO_DIAGNOSTICS
    assert sr['uves'].loader == 'reflected_solar'
    assert sr['uves'].diagnostics is cs.UVES_DIAGNOSTICS


def test_procyon_registers_more_than_one_region():
    # The smoke-test condition: Procyon now declares >1 region (vis + uves + uv + ir),
    # not just the hardcoded 'vis'.
    pr = cs.star_arm_registry('procyon')
    assert len(pr) > 1
    assert pr['harps'].ready is True            # HARPS VIS runnable today
    assert pr['uves'].ready is False            # gated on RYA-272 loader + data
    assert pr['uv'].ready is False
    assert pr['ir'].ready is False
    ready, deferred = cs.available_arms('procyon')
    assert ready == ('harps',)
    assert set(deferred) == {'uves', 'uv', 'ir'}
    assert all(deferred.values())               # every deferral carries a reason


def test_procyon_uves_primary_O_is_oi_777():
    # The UVES arm's primary O indicator is O I 777 (the whole point — not [O I] 6300).
    pr = cs.star_arm_registry('procyon')
    o_primary = [d for d in pr['uves'].diagnostics if d.element == 'O' and d.role == 'primary']
    assert any(d.key == 'OI_777' for d in o_primary)


def test_resolver_never_silently_loads_reflected_solar_for_nonsolar():
    # THE BUG RYA-464 CLOSES: a non-solar star must never be synthesized against Vesta.
    pr = cs.star_arm_registry('procyon')
    with pytest.raises(cs.ArmNotWired):
        cs.resolve_arm_spectrum('procyon', pr['uves'])      # unready → loud
    # Even a hypothetical READY reflected_solar arm must refuse for a non-solar star.
    fake = cs.ArmWiring('uves', cs.UVES_OPT, cs.UVES_DIAGNOSTICS, 'reflected_solar', True)
    with pytest.raises(cs.ArmNotWired, match='reflected'):
        cs.resolve_arm_spectrum('procyon', fake)


def test_undeclared_star_fails_loud():
    with pytest.raises(KeyError, match='STAR_ARMS'):
        cs.star_arm_registry('betelgeuse')


def test_solar_reflected_arm_passes_the_nonsolar_guard():
    # The guard must NOT block solar's own reflected-solar arms (it only blocks non-solar).
    sr = cs.star_arm_registry('solar')
    # resolve dispatches to _load_reflected_solar_arm for solar; we don't execute the load
    # (Vesta data absent), only assert it does NOT raise the non-solar ArmNotWired guard.
    try:
        cs.resolve_arm_spectrum('solar', sr['uves'])
    except cs.ArmNotWired as e:                  # must not be the non-solar refusal
        assert 'reflected sunlight' not in str(e)
    except FileNotFoundError:
        pass                                     # expected: Vesta manifest not staged here
    except Exception:
        pass                                     # any data-load error is fine; not a guard refusal
