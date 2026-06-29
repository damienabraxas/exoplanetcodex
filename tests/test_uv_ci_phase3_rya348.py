"""
tests/test_uv_ci_phase3_rya348.py
=================================
RYA-348 Phase 3 — UV arm C I discriminator. The Phase-3 run found the FUV C I 1657
synthesis is NOT science-grade (railed under both continuum choices) and the solar UV
reference cannot resolve UV C I, so the arm stays deferred. These tests pin that honest
state in code (no external data needed). The real-spectrum fit lives in the harness
scripts/procyon_uv_ci_phase3_rya348.py (skipped here when the STIS frame is absent).
"""
from pathlib import Path

import numpy as np
import pytest

import pipeline.cno_synthesis as cs
from pipeline import uv_conditioning as uvc

ROOT = Path(__file__).resolve().parents[1]


def test_uv_arm_stays_deferred_after_phase3():
    # Running the synthesis end-to-end produced garbage (railed) -> honest: NOT ready.
    uv = cs.star_arm_registry('procyon')['uv']
    assert uv.ready is False
    ready, deferred = cs.available_arms('procyon')
    assert 'uv' not in ready
    assert 'uv' in deferred


def test_uv_defer_reason_names_the_two_phase3_blockers():
    reason = cs.star_arm_registry('procyon')['uv'].defer_reason.lower()
    # blocker 1: the FUV pseudo-continuum / railed fit
    assert 'pseudo-continuum' in reason or 'line forest' in reason
    assert 'rail' in reason
    # blocker 2: the solar UV reference cannot resolve UV C I
    assert 'calspec' in reason or 'solar uv reference' in reason


def test_uv_arm_resolver_never_silently_loads_vesta():
    uv = cs.star_arm_registry('procyon')['uv']
    forced_solar = cs.ArmWiring('uv', uv.region, uv.diagnostics, 'hst_stis', True)
    with pytest.raises(cs.ArmNotWired):
        cs.resolve_arm_spectrum('solar', forced_solar)


def test_fuv_ci_carries_a_grid_owed_flag_not_an_applied_correction():
    # C I 1657 is a diagnostic but its NLTE is owed, never silently applied as a scalar.
    diags = {d.key: d for d in cs.star_arm_registry('procyon')['uv'].diagnostics}
    assert 'CI_1657_fuv' in diags
    assert 'owed' in diags['CI_1657_fuv'].nlte_flag.lower()
    # and the FUV key is unknown to the cited-correction grid (no fabricated delta)
    assert 'CI_1657_fuv' not in cs._ATOMIC_GRID_KEYS


def test_vac_to_air_named_line_mg_ii_k():
    air = float(uvc.to_pipeline_frame(np.array([2796.3543]))[0])
    assert abs(air - 2795.528) < 0.01            # known air wavelength
    # FUV line below 2000 A stays vacuum (unchanged)
    assert abs(float(uvc.to_pipeline_frame(np.array([1657.380]))[0]) - 1657.380) < 1e-6
