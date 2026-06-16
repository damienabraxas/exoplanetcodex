# =============================================================================
# THE EXOPLANET CODEX — tests
# -----------------------------------------------------------------------------
# RYA-321 — VALD intake metallicity gate.
#
# VALD3 "Extract Stellar" sets composition via the free-text Chemical
# composition field; if it is left blank the extraction runs at the solar
# default and a metal-rich star's weak-line pool is under-selected. The applied
# composition is recorded in the footer ABUNDANCE BLOCK (all metals offset from
# VALD's solar scale by the requested [M/H]) — NOT the Castelli model filename,
# which VALD substitutes to 'castelli_ap00...' (solar structure) when it lacks
# the exact metal-rich model even though M/H was applied (job.019562, 2026-06-15).
# These tests lock the abundance-block read so (a) a genuine solar-default
# delivery REJECTs and (b) a correct metal-rich delivery ACCEPTs even when its
# structure-model filename reads ap00.
# =============================================================================

import pytest

from data.linelists.vald_parse import (
    parse_vald_applied_metallicity, parse_vald_structure_grid,
    verify_metallicity, _VALD_SOLAR_ABUND)


def _abundance_block(mh):
    """The footer abundance lines for an extraction at applied [M/H]=mh: every
    metal in VALD's solar scale offset by mh, H/He fixed."""
    lines = ["'H :  0.92','He: -1.11',"]
    for el, sol in _VALD_SOLAR_ABUND.items():
        lines.append(f"'{el}: {sol + mh:.2f}',")
    lines.append("'END'")
    return "\n".join(lines)


def _write_delivery(tmp_path, mh, structure_grid='castelli_ap00k2_T05750G45.krz'):
    """Minimal VALD-shaped file: a transition record, the substituted structure
    model filename, then the abundance block at applied [M/H]=mh."""
    p = tmp_path / 'delivery.txt'
    p.write_text(
        " 3780.0, 6910.0, 100, 1, 1.0 Wavelength region, ...\n"
        "Spec Ion       WL_air(A)  log gf* ...\n"
        "'Fe 1',       3780.25500, -2.876, ...\n"
        f"'{structure_grid}',\n"
        f"{_abundance_block(mh)}\n"
    )
    return str(p)


# --- parse_vald_applied_metallicity (abundance block) ------------------------

@pytest.mark.parametrize('mh', [0.0, 0.2, 0.35, -0.5, 0.31])
def test_parse_applied_metallicity(tmp_path, mh):
    got = parse_vald_applied_metallicity(_write_delivery(tmp_path, mh))
    assert got == pytest.approx(mh, abs=1e-6)


def test_applied_metallicity_ignores_structure_filename(tmp_path):
    """The whole RYA-321 correction: +0.35 composition with an ap00 structure
    model (the substitution case) reads as +0.35, not 0.0."""
    f = _write_delivery(tmp_path, 0.35, structure_grid='castelli_ap00k2_T05250G45.krz')
    assert parse_vald_applied_metallicity(f) == pytest.approx(0.35)
    assert parse_vald_structure_grid(f) == pytest.approx(0.0)  # informational only


def test_parse_returns_none_when_no_block(tmp_path):
    p = tmp_path / 'no_block.txt'
    p.write_text("'Fe 1', 3780.0, -2.0, ...\n'castelli_ap00k2_T05750G45.krz',\n")
    assert parse_vald_applied_metallicity(str(p)) is None


def test_parse_returns_none_on_nonuniform_scaling(tmp_path):
    """A half-corrupted block where metals disagree must fail loud, not average."""
    p = tmp_path / 'bad.txt'
    p.write_text("'Fe: -4.34','Mg: -4.46','Si: -4.49','Ca: -5.68',"
                 "'Na: -5.71','Ti: -7.02','Cr: -6.37','Ni: -5.79','END'\n")
    assert parse_vald_applied_metallicity(str(p)) is None


# --- verify_metallicity gate --------------------------------------------------

STAR_PARAMS = {
    'solar':       {'feh_ref': 0.00},
    'procyon':     {'feh_ref': 0.03},
    'alpha_cen_a': {'feh_ref': 0.20},
    '55cnc_a':     {'feh_ref': 0.31},
    'no_feh':      {'teff': 5000},
}


def test_solar_default_metal_rich_rejects(tmp_path):
    """The core catch: a metal-rich star whose composition ran solar -> REJECT."""
    f = _write_delivery(tmp_path, 0.0)
    verdict, msg = verify_metallicity(f, 'alpha_cen_a', STAR_PARAMS)
    assert verdict == 'REJECT' and 'METALLICITY MISMATCH' in msg


def test_correct_metal_rich_accepts_despite_ap00_structure(tmp_path):
    """Regression for the false-reject bug: α Cen extracted at +0.20 with the
    substituted ap00 structure model must ACCEPT."""
    f = _write_delivery(tmp_path, 0.20, structure_grid='castelli_ap00k2_T05750G45.krz')
    assert verify_metallicity(f, 'alpha_cen_a', STAR_PARAMS)[0] == 'ACCEPT'


def test_solar_anchor_accepts(tmp_path):
    f = _write_delivery(tmp_path, 0.0)
    assert verify_metallicity(f, 'solar', STAR_PARAMS)[0] == 'ACCEPT'


def test_reference_difference_within_tol_accepts(tmp_path):
    """55 Cnc applied +0.35 vs catalog +0.31 (Δ0.04 < 0.05 tol) -> ACCEPT."""
    f = _write_delivery(tmp_path, 0.35)
    assert verify_metallicity(f, '55cnc_a', STAR_PARAMS)[0] == 'ACCEPT'


def test_missing_block_rejects(tmp_path):
    p = tmp_path / 'no_block.txt'
    p.write_text("'Fe 1', 3780.0, -2.0, ...\n'END'\n")
    verdict, msg = verify_metallicity(str(p), 'alpha_cen_a', STAR_PARAMS)
    assert verdict == 'REJECT' and 'abundance block' in msg


def test_unwired_star_rejects(tmp_path):
    """No STAR_PARAMS feh_ref -> cannot verify -> REJECT (no silent fallback)."""
    f = _write_delivery(tmp_path, 0.0)
    for sid in ('not_a_star', 'no_feh'):
        verdict, msg = verify_metallicity(f, sid, STAR_PARAMS)
        assert verdict == 'REJECT' and 'feh_ref' in msg
