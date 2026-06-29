"""
tests/test_uv_conditioning_rya426.py
====================================
RYA-426 — the UV (FUV/NUV) conditioning stage. Pins the DISCIPLINE, not science numbers:
the vacuum->air convention (NUV converts, FUV stays vacuum), the FUV synthesis-not-EW
refusal, chromospheric-core masking, scattered-light verification, ISM + NLTE-coverage
flags, the window-to-diagnostic map, the analysis_ready composite, and the loud-failure
rules (excluded/unknown grating, binary identity). Synthetic arrays only — runs in CI.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import uv_conditioning as uv  # noqa: E402


# ── gate 2: vacuum -> air (Birch & Downs 1994) ───────────────────────────────
def test_vac_to_air_known_mg_ii_lines():
    air = uv.vac_to_air(np.array([2796.3543, 2803.5324]))   # Mg II k, h (vacuum)
    assert abs(air[0] - 2795.528) < 0.01                    # known air (Morton/NIST)
    assert abs(air[1] - 2802.705) < 0.01


def test_fuv_stays_vacuum_below_2000():
    w = np.array([1355.6, 1657.38, 1999.9])
    assert np.allclose(uv.vac_to_air(w), w)                 # identity below the boundary
    # at/above 2000 the wavelength shifts (air < vacuum)
    assert uv.vac_to_air(np.array([2500.0]))[0] < 2500.0


def test_air_vac_roundtrip_nuv():
    w = np.array([2500.0, 3000.0])
    assert np.allclose(uv.air_to_vac(uv.vac_to_air(w)), w, atol=1e-3)


# ── regime + gate 5: FUV synthesis-not-EW ────────────────────────────────────
def test_classify_regime_boundary():
    assert uv.classify_regime(1657.0) == 'FUV'
    assert uv.classify_regime(2000.0) == 'NUV'
    assert uv.classify_regime(2796.0) == 'NUV'


def test_refuse_ew_in_fuv_raises_for_fuv_only():
    with pytest.raises(uv.FUVSynthesisRequired):
        uv.refuse_ew_in_fuv(np.array([1657.38]))
    with pytest.raises(uv.FUVSynthesisRequired):
        uv.refuse_ew_in_fuv(np.array([2500.0, 1900.0]))     # any FUV pixel refuses
    uv.refuse_ew_in_fuv(np.array([2796.0, 3000.0]))         # pure NUV — no raise


# ── gate 4: chromospheric-core masking ───────────────────────────────────────
def test_chromospheric_mask_catches_mg_ii_in_pipeline_frame():
    # grid in pipeline (air) frame around Mg II k air 2795.53
    w = np.linspace(2794.0, 2797.0, 301)
    mask = uv.chromospheric_core_mask(w)
    assert mask.any()
    assert mask[np.argmin(np.abs(w - 2795.53))]            # core pixel masked


def test_chromospheric_mask_catches_fuv_core_in_vacuum():
    w = np.linspace(1334.0, 1336.0, 201)                   # C II 1334/1335 (vacuum, FUV)
    assert uv.chromospheric_core_mask(w).any()


# ── gate 3: scattered light ──────────────────────────────────────────────────
def test_scattered_light_flags_high_negative_fraction():
    clean = uv.scattered_light_check(np.array([1.0, 1.1, 0.9, 1.0]))
    assert clean['passed'] and clean['neg_flux_fraction'] == 0.0
    over = uv.scattered_light_check(np.array([-1.0, -1.0, -1.0, 1.0]))   # 75% negative
    assert not over['passed']


# ── gate 6: ISM ──────────────────────────────────────────────────────────────
def test_ism_flag_by_distance():
    assert uv.ism_review_flag(3.51)['review_required'] is False    # Procyon, nearby
    assert uv.ism_review_flag(120.0)['review_required'] is True


# ── gate 7: NLTE coverage (no UV grids exist yet -> LTE loud) ─────────────────
def test_nlte_coverage_all_lte_flagged_loud():
    flags = uv.nlte_coverage_flags()
    assert flags and all(f['lte_flag'] == 'LTE_ASSUMED_LOUD' for f in flags)
    assert all(f['uv_nlte_grid'] == 'NONE' for f in flags)


# ── gate 8: window-to-diagnostic map ─────────────────────────────────────────
def test_window_map_flags_coverage_and_gaps():
    windows = [{'grating': 'E230H', 'lo_A': 2620.0, 'hi_A': 3158.0},   # NUV, no FUV anchor
               {'grating': 'E140M', 'lo_A': 1140.0, 'hi_A': 1730.0}]   # FUV C I 1657 etc.
    m = uv.window_to_diagnostic_map(windows)
    covered = {c['species'] + str(c['lambda_A']) for c in m['covered']}
    assert 'C I1657.38' in covered                          # covered by E140M
    # an NUV-only window covers no FUV anchors -> appears as spectrum-without-lines if isolated
    only_nuv = uv.window_to_diagnostic_map([{'grating': 'E230H', 'lo_A': 2620.0, 'hi_A': 3158.0}])
    assert 'E230H' in only_nuv['spectra_without_lines']
    assert only_nuv['lines_without_spectra']                # FUV anchors not covered


# ── full chain: condition_uv_spectrum ────────────────────────────────────────
def _fuv():
    w = np.linspace(1163.0, 1357.0, 2000)                   # E140H FUV (vacuum)
    return w, np.ones_like(w)


def _nuv():
    w = np.linspace(2620.0, 3158.0, 2000)                   # E230H NUV (vacuum)
    return w, np.ones_like(w)


def test_condition_fuv_refuses_ew_and_is_ready(tmp_path, monkeypatch):
    monkeypatch.setattr(uv, 'CONDITIONING_DIR', tmp_path)
    w, f = _fuv()
    res = uv.condition_uv_spectrum(w, f, star='procyon', grating='E140H', distance_pc=3.51)
    assert res.regime == 'FUV'
    assert res.gates['synthesis_gate']['ew_allowed'] is False    # FUV -> synthesis only
    assert res.telluric_applied is False                         # space UV, always
    assert res.analysis_ready is True


def test_condition_nuv_allows_ew_and_masks_mg_ii(tmp_path, monkeypatch):
    monkeypatch.setattr(uv, 'CONDITIONING_DIR', tmp_path)
    w, f = _nuv()
    res = uv.condition_uv_spectrum(w, f, star='procyon', grating='E230H', distance_pc=3.51)
    assert res.regime == 'NUV'
    assert res.gates['synthesis_gate']['ew_allowed'] is True
    names = [l['name'] for l in res.gates['chromospheric_mask']['lines_masked_in_range']]
    assert 'Mg II k' in names and 'Mg II h' in names
    assert res.analysis_ready is True


def test_condition_writes_manifest(tmp_path, monkeypatch):
    monkeypatch.setattr(uv, 'CONDITIONING_DIR', tmp_path)
    w, f = _nuv()
    uv.condition_uv_spectrum(w, f, star='procyon', grating='E230H')
    out = tmp_path / 'procyon' / 'procyon_E230H_uv_conditioning.json'
    assert out.exists()
    import json
    m = json.loads(out.read_text())
    assert m['schema'] == uv._MANIFEST_VERSION and m['analysis_ready'] is True


# ── loud failure over silent fallback ────────────────────────────────────────
def test_excluded_grating_raises():
    with pytest.raises(uv.UVConditioningError):
        uv.condition_uv_spectrum(np.array([2500.0]), np.array([1.0]), star='x', grating='WFC3')


def test_unknown_grating_raises():
    with pytest.raises(uv.UVConditioningError):
        uv.condition_uv_spectrum(np.array([2500.0]), np.array([1.0]), star='x', grating='NOPE')


def test_binary_indeterminate_identity_raises():
    w, f = _nuv()
    with pytest.raises(uv.UVConditioningError):
        uv.condition_uv_spectrum(w, f, star='alpha_cen_a', grating='E230H',
                                 is_binary=True, identity_confirmed=False, write=False)


def test_single_star_skips_identity_gate():
    w, f = _nuv()
    res = uv.condition_uv_spectrum(w, f, star='procyon', grating='E230H', write=False)
    assert res.gates['identity']['note'].startswith('single star')
