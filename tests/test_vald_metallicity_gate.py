# =============================================================================
# THE EXOPLANET CODEX — tests
# -----------------------------------------------------------------------------
# RYA-321 — VALD intake metallicity gate.
#
# VALD3 "Extract Stellar" sets composition via the free-text Chemical
# composition field; if it is left blank the extraction runs at the solar
# default and a metal-rich star's weak-line pool is under-selected. The applied
# composition is recorded only in the footer model-atmosphere grid filename
# (e.g. 'castelli_ap00k2_T05750G45.krz' -> [M/H] +0.0), NOT the header. These
# tests lock the parse + gate so a wrong-composition delivery is REJECTed like a
# truncated or HFS-OFF file.
# =============================================================================

import pytest

from data.linelists.vald_parse import (
    parse_vald_model_metallicity, verify_metallicity, _nearest_castelli_node)


def _write_delivery(tmp_path, grid_name):
    """Minimal VALD-shaped file: a couple of transition records then the footer
    model-atmosphere block carrying the Castelli grid filename."""
    p = tmp_path / 'delivery.txt'
    p.write_text(
        " 3780.00000, 6910.00000, 100, 17107375, 1.0 Wavelength region, ...\n"
        "Spec Ion       WL_air(A)  log gf* E_low(eV) ...\n"
        "'Fe 1',       3780.25500, -2.876,  3.9286, ...\n"
        f"'{grid_name}',\n"
        "'H :  0.92','He: -1.11',\n"
        "'Fe: -4.54','END'\n"
    )
    return str(p)


# --- parse_vald_model_metallicity --------------------------------------------

@pytest.mark.parametrize('grid, mh', [
    ('castelli_ap00k2_T05750G45.krz', 0.0),
    ('castelli_ap02k2_T05750G45.krz', 0.2),
    ('castelli_ap05k2_T05250G45.krz', 0.5),
    ('castelli_am05k2_T05250G45.krz', -0.5),
    ('castelli_am25k2_T05250G45.krz', -2.5),
])
def test_parse_grid_token(tmp_path, grid, mh):
    assert parse_vald_model_metallicity(_write_delivery(tmp_path, grid)) == pytest.approx(mh)


def test_parse_returns_none_when_absent(tmp_path):
    p = tmp_path / 'no_model.txt'
    p.write_text("'Fe 1', 3780.0, -2.0, ...\n'H :  0.92','END'\n")
    assert parse_vald_model_metallicity(str(p)) is None


# --- nearest Castelli node ----------------------------------------------------

@pytest.mark.parametrize('feh, node', [
    (0.00, 0.0), (0.03, 0.0),       # Procyon ~solar snaps to the 0.0 node
    (0.20, 0.2), (0.31, 0.2),       # 55 Cnc +0.31 -> nearest node +0.2
    (0.45, 0.5), (-0.4, -0.5),
])
def test_nearest_node(feh, node):
    assert _nearest_castelli_node(feh) == pytest.approx(node)


# --- verify_metallicity gate --------------------------------------------------

STAR_PARAMS = {
    'solar':       {'feh_ref': 0.00},
    'procyon':     {'feh_ref': 0.03},
    'alpha_cen_a': {'feh_ref': 0.20},
    '55cnc_a':     {'feh_ref': 0.31},
    'no_feh':      {'teff': 5000},          # record without feh_ref
}


def test_solar_default_metal_rich_rejects(tmp_path):
    """The core RYA-321 catch: a metal-rich star extracted at solar -> REJECT."""
    f = _write_delivery(tmp_path, 'castelli_ap00k2_T05750G45.krz')
    verdict, msg = verify_metallicity(f, 'alpha_cen_a', STAR_PARAMS)
    assert verdict == 'REJECT'
    assert 'METALLICITY MISMATCH' in msg


def test_solar_anchor_accepts(tmp_path):
    f = _write_delivery(tmp_path, 'castelli_ap00k2_T05750G45.krz')
    assert verify_metallicity(f, 'solar', STAR_PARAMS)[0] == 'ACCEPT'


def test_procyon_near_solar_accepts(tmp_path):
    """+0.03 snaps to the 0.0 node, so a solar-grid Procyon delivery is fine."""
    f = _write_delivery(tmp_path, 'castelli_ap00k2_T06500G40.krz')
    assert verify_metallicity(f, 'procyon', STAR_PARAMS)[0] == 'ACCEPT'


def test_correct_metal_rich_accepts(tmp_path):
    """A faithful +0.31 re-extraction is served as the +0.2 node -> ACCEPT
    (the gate compares to the nearest node, not the raw catalog value)."""
    f = _write_delivery(tmp_path, 'castelli_ap02k2_T05250G45.krz')
    assert verify_metallicity(f, '55cnc_a', STAR_PARAMS)[0] == 'ACCEPT'


def test_missing_grid_record_rejects(tmp_path):
    p = tmp_path / 'no_model.txt'
    p.write_text("'Fe 1', 3780.0, -2.0, ...\n'END'\n")
    verdict, msg = verify_metallicity(str(p), 'alpha_cen_a', STAR_PARAMS)
    assert verdict == 'REJECT' and 'no model-atmosphere metallicity' in msg


def test_unwired_star_rejects(tmp_path):
    """No STAR_PARAMS feh_ref -> cannot verify -> REJECT (no silent fallback)."""
    f = _write_delivery(tmp_path, 'castelli_ap00k2_T05750G45.krz')
    for sid in ('not_a_star', 'no_feh'):
        verdict, msg = verify_metallicity(f, sid, STAR_PARAMS)
        assert verdict == 'REJECT' and 'feh_ref' in msg
