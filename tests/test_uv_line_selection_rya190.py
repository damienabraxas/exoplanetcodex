"""
tests/test_uv_line_selection_rya190.py
======================================
RYA-190 — the UV line selection + NLTE policy (single source of truth). Pins the cited
verdicts, the FUV synthesis discipline, the grids-not-scalars NLTE policy, the optical
cross-check pairing the RYA-471 leg validation needs, and the RYA-426 window-map
integration. Implements policy only — asserts NO production science change.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import uv_line_selection as sel  # noqa: E402


# ── cited verdicts (the policy) ──────────────────────────────────────────────
def test_carbon_uv_1657_usable_synthesis_nlte_cited():
    d = sel.lookup('C I', 1657.38)
    assert d['verdict'] == 'USE' and d['method'] == 'synthesis'
    assert d['nlte_status'] == 'GRID_OWED' and abs(d['nlte_expected_dex'] - 0.10) < 1e-9
    assert 'Amarsi' in d['nlte_ref']


def test_oxygen_1355_use_but_1302_resonance_is_a_trap():
    assert sel.is_usable('O I', 1355.60)
    o1302 = sel.lookup('O I', 1302.17)
    assert o1302['verdict'] == 'DO_NOT_USE' and o1302['nlte_status'] == 'DO_NOT_USE'


def test_nitrogen_uv_1200_is_a_trap():
    n = sel.lookup('N I', 1199.55)
    assert n['verdict'] == 'DO_NOT_USE'
    assert n in sel.traps()


def test_sulphur_uv_usable():
    assert sel.is_usable('S I', 1473.99)


# ── FUV synthesis discipline (mirrors RYA-426 gate 5) ────────────────────────
def test_every_usable_fuv_line_is_synthesis_not_ew():
    for d in sel.usable_diagnostics(regime='FUV'):
        assert d['method'] == 'synthesis', d


# ── grids-not-scalars NLTE policy ────────────────────────────────────────────
def test_uv_nlte_is_grid_owed_not_applied():
    # every usable atomic UV line is GRID_OWED (no UV grid yet) — magnitude is informational
    for d in sel.usable_diagnostics():
        if d['regime'] in ('FUV', 'NUV') and d['nlte_status'] not in ('NA_MOLECULAR',):
            assert d['nlte_status'] == 'GRID_OWED', d
    pol = sel.nlte_policy('C I', 1657.38)
    assert pol['nlte_status'] == 'GRID_OWED' and 'grid owed' in pol['note']


def test_module_applies_no_correction_to_production_abundances():
    # discipline: this module is policy-only — it must NOT wire a scalar into abundances_derive
    import inspect
    import pipeline.uv_line_selection as m
    src = inspect.getsource(m)
    assert 'abundances_derive' not in src           # no production-path coupling
    assert 'import abundances_derive' not in src


# ── optical cross-check (the RYA-471 leg-validation pair) ────────────────────
def test_optical_cross_check_for_carbon():
    waves = sorted(round(d['wavelength_A'], 1) for d in sel.optical_cross_check('C'))
    assert waves == [5052.2, 5380.3]
    for d in sel.optical_cross_check('C'):
        assert d['regime'] == 'optical' and d['verdict'] == 'USE'


# ── RYA-426 window-map integration ───────────────────────────────────────────
def test_anchors_are_vacuum_and_feed_window_map_shape():
    anchors = sel.anchors_for_window_map(usable_only=True)
    # optical lines excluded; each anchor has the {species, element, lambda_A} shape
    assert anchors and all({'species', 'element', 'lambda_A'} <= set(a) for a in anchors)
    assert all(a['lambda_A'] < 3500 for a in anchors)          # UV only, no optical 5052
    # NH (NUV, air 3360) converted to vacuum
    nh = [a for a in anchors if a['species'] == 'NH'][0]
    assert nh['lambda_A'] > 3360.0                              # air->vac shifts up


def test_air_to_vac_known_line():
    # Mg II k air 2795.528 -> vacuum 2796.354 (inverse check on the BD1994 helper)
    assert abs(sel._air_to_vac(2795.528) - 2796.354) < 0.01
    assert sel._air_to_vac(1657.38) == 1657.38                  # FUV identity


# ── controlled vocabulary / fail-loud ────────────────────────────────────────
def test_validation_rejects_unknown_vocab():
    bad = dict(sel.UV_DIAGNOSTICS[0])
    bad['verdict'] = 'MAYBE'
    with pytest.raises(AssertionError):
        # re-run the same invariant the module asserts at import
        assert bad['verdict'] in sel.VERDICTS


def test_lookup_miss_returns_none_and_unknown_policy():
    assert sel.lookup('Fe I', 1500.0) is None
    assert sel.nlte_policy('Fe I', 1500.0)['nlte_status'] == 'UNKNOWN'
