"""
tests/test_solar_ba_deblend_rya581.py
=====================================
RYA-581 — Ba II 5853 in-window DEBLEND, superseding the RYA-559 EW->COG value.

RYA-559 landed A(Ba)_NLTE 2.410 by inverting the OBSERVED pool EW (74.62 mA) through an
HFS-resolved curve of growth. That EW carries blend_flag=True — ~10 mA over the clean
solar Ba II 5853 (~64 mA) — and an EW inversion CANNOT deblend: a single scalar has no
way to separate barium from the rest of the absorption inside the integration window, so
the neighbours were charged to Ba. RYA-559 recorded that and routed the debt here.

RYA-581 fits the PROFILE instead (the RYA-551 Sr II in-window pattern): the full VALD3
in-window block is synthesised alongside the Ba II HFS/isotope components and A(Ba) is
fitted by chi2, with the Engine-A Korotin2015 delta read at the solar node as before.

CI-safe (no synthesis): the measurement runs on Sirius via
scripts/rya581_ba2_deblend_sirius.py; this tests the committed result JSON and the
phase_c fold. Mirrors test_solar_ba_synthesis_rya559.py.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import phase_c_verdict_rya371 as V                                  # noqa: E402
from config.constants import SOLAR_ASPLUND2021                      # noqa: E402

JSON = ROOT / 'data' / 'results' / 'solar_ba_deblend_rya581.json'
JSON_559 = ROOT / 'data' / 'results' / 'solar_ba_synthesis_rya559.json'


def _data():
    return json.loads(JSON.read_text())


# ── the measurement product exists and is well-formed ─────────────────────────

def test_measurement_json_present_and_schema():
    d = _data()
    assert d['ticket'] == 'RYA-581'
    assert d['element'] == 'Ba' and d['ion'] == 'II' and d['line_A'] == 5853.668
    assert d['A_nlte'] is not None and d['A_lte'] is not None
    # Same HFS treatment and same canonical gf scale as RYA-559 — the deblend is the
    # ONLY thing that changed, so the comparison against 2.410 is apples-to-apples.
    assert d['hfs']['n_components'] >= 20
    assert abs(d['hfs']['canonical_loggf'] - (-1.0)) < 0.02        # NIST / Davidson 1992
    assert d['hfs']['canonical_line_id'] == 'gf_001470'            # single source of truth
    # Engine-A production Korotin delta, read at the solar node (validate-don't-tune).
    assert d['engineA_in_bounds'] is True
    assert -0.05 < d['engineA_korotin_delta'] < 0.0    # solar Ba II 5853 small-negative
    assert abs(d['A_nlte'] - (d['A_lte'] + d['engineA_korotin_delta'])) < 0.005


# ── the blend was actually MODELLED — the RYA-581 stop condition ───────────────

def test_in_window_blend_is_present_not_culled():
    # RYA-581's whole point: "if the fit still returns ~2.41 the blend was not modelled".
    # The blend has to be IN the window, in quantity.
    d = _data()
    bm = d['blend_model']
    assert bm['n_rows'] > 100 and bm['n_species'] > 10
    assert bm['ba_ii_rows_removed'] > 0        # Ba comes from the HFS list, never doubled
    assert bm['blend_core_EW_mA'] > 1.0        # real absorption under the core
    # Fe I is the known contaminant of the solar Ba II 5853 window; it must dominate.
    per_species = bm['per_species_core_EW_mA']
    assert max(per_species, key=per_species.get) == 'Fe I'


def test_deblend_is_demonstrated_by_fit_quality():
    # The evidence that the blend is modelled rather than culled is that modelling it
    # EXPLAINS observed structure — i.e. the profile fit gets materially better. The
    # shift in A alone would be a weak proxy: this blend sits largely in the wings.
    ev = _data()['deblend_evidence']
    assert ev['red_chi2_blends_modelled'] < ev['red_chi2_ba_alone']
    assert ev['chi2_improvement_factor'] >= 1.5


def test_result_moved_off_the_rya559_ew_cog_value():
    # The RYA-581 stop condition itself: a fit that still returns ~2.41 means the blend
    # was culled again and must NOT be shipped.
    d = _data()
    assert d['rya559_ew_cog_A_nlte'] == 2.410
    assert abs(d['A_nlte'] - 2.410) > 0.03


def test_correction_budget_is_honest_about_where_the_shift_came_from():
    # Most of 2.410 -> 2.237 is abandoning the EW inversion, NOT the in-window blend
    # model. Recording that split is what keeps the result honest.
    cb = _data()['correction_budget_dex']
    assert cb['rya559_ew_cog'] == 2.410
    assert abs(cb['total'] - (cb['from_dropping_the_EW_inversion']
                              + cb['from_modelling_the_in_window_blend'])) < 0.002
    # The independent check that the profile fit found the CLEAN line: its synthetic core
    # EW must land on the literature / RYA-559-calibration clean EW (~64-66 mA), NOT on
    # the blend-inflated pool EW (74.62 mA) that produced 2.410.
    assert 62.0 < cb['fitted_core_EW_mA'] < 69.0


# ── validate-don't-tune ────────────────────────────────────────────────────────

def test_not_tuned_toward_asplund():
    d = _data()
    # The answer is whatever chi2 returned; it is NOT pinned to Asplund. If it had been
    # tuned, delta_vs_asplund would be ~0 and the two arms would not disagree.
    assert d['delta_vs_asplund'] != 0.0
    assert 'never fitted toward Asplund' in d['validate_dont_tune'] \
        or 'Never fitted toward Asplund' in d['validate_dont_tune']
    # Two independent solar arms (HARPS + IAG) fitted separately, and the fit-window
    # choice swept — a value that only exists at one window would not be trustworthy.
    harps = d['per_arm']['harps']['deblended']
    assert len(harps) >= 3
    vals = [v['A_NLTE'] for v in harps.values()]
    assert max(vals) - min(vals) < 0.05        # window-choice insensitive


def test_fit_is_reliable_not_railed():
    d = _data()
    prim = d['per_arm']['harps']['deblended'][f"{d['fit']['half_window_A']:.1f}"]
    assert prim['railed'] is False
    assert prim['gsig_railed'] is False        # RYA-643: gsig railing was a real defect
    assert prim['dEW_dA_mA_dex'] >= d['fit']['reliable_dEW_dA_floor']
    assert prim['red_chi2'] <= d['fit']['reliable_red_chi2_ceiling']
    assert d['reliable'] is True
    # RYA-643: the residual velocity must be MEASURED per arm, never a silent zero.
    for arm in d['arms'].values():
        assert arm['n_lines'] >= 4 and arm['v_kms'] != 0.0


# ── the phase_c fold consumes the DEBLEND, not the superseded 559 value ────────

def test_loader_prefers_the_deblend_over_rya559():
    assert JSON_559.exists()                   # the superseded record is still on disk
    assert V._load_ba_synthesis()['ticket'] == 'RYA-581'


def test_ba_reclassify_lands_the_deblended_value():
    d = _data()
    ov = V._ba_reclassify(d)
    assert set(ov) == {'Ba'}
    assert abs(ov['Ba']['A_measured'] - d['A_nlte']) < 1e-6
    assert abs(ov['Ba']['A_measured'] - 2.410) > 0.03      # NOT the superseded value
    assert 'in-window' in ov['Ba']['channel'] and '5853' in ov['Ba']['channel']
    assert 'Korotin' in ov['Ba']['channel']
    assert 'RYA-581' in ov['Ba']['provenance']


def test_verdict_follows_the_measurement_and_names_its_caveats():
    d = _data()
    ov = V._ba_reclassify(d)['Ba']
    reconciled = abs(d['A_nlte'] - SOLAR_ASPLUND2021['Ba']) <= V.TOL_PASS
    # RYA-581 says explicitly: if it reconciles, report PASS honestly — do not
    # force-hold it owed. The verdict must track the measurement, not a preference.
    assert ov['verdict'] == ('PASS' if (reconciled and d['reliable']) else 'CURATION-OWED')
    owed = ov['owed']
    assert 'SUPERSEDES' in owed and '2.410' in owed
    assert 'do NOT tune' in owed
    if ov['verdict'] == 'PASS':
        # A single-line PASS must carry its caveats, including the RYA-669 finding that
        # gate 3 still cannot see this value (Ba has no two-engine record).
        assert 'ONE line' in owed
        assert 'RYA-527' in owed               # flagged for the freeze review
        assert 'RYA-669' in owed               # gate 3 UNEVALUABLE


def test_rya559_record_still_folds_on_its_own_terms():
    # The 559 branch must keep working for a checkout that predates the deblend —
    # superseding the value must not break the fallback.
    ov = V._ba_reclassify(json.loads(JSON_559.read_text()))
    assert ov['Ba']['verdict'] == 'CURATION-OWED'
    assert abs(ov['Ba']['A_measured'] - 2.410) < 1e-6
