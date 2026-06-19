"""
tests/test_nlte_cno.py
======================
Tests for pipeline/nlte_cno.py — the C I / O I non-LTE correction grids
(RYA-359, Amarsi 2019 A&A 630 A104, CDS J/A+A/630/A104).

Covers: vendoring + provenance, the Phase-1 coverage STOP-GATE, the nm→Å line
resolver (incl. O I 777 multiplet averaging), 3D/1D leg selection, correct
(negative) sign + sign guard, the table7 self-validation, and the public
apply_cno_nlte_corrections wiring.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipeline import nlte_cno as cno

DATA_DIR = Path(__file__).parent.parent / 'data' / 'nlte_grids' / 'amarsi2019_cno'


# ── Vendoring + provenance ──────────────────────────────────────────────────────

class TestVendoring:
    def test_all_tables_vendored(self):
        for f in ['ReadMe', 'table1.dat', 'table2.dat.gz', 'table3.dat.gz',
                  'table5.dat.gz', 'table6.dat.gz', 'table7.dat', 'provenance.json']:
            assert (DATA_DIR / f).exists(), f"missing vendored file {f}"

    def test_provenance_manifest(self):
        prov = json.loads((DATA_DIR / 'provenance.json').read_text())
        assert 'J/A+A/630/A104' in prov['source']
        assert 'Amarsi' in prov['citation']
        assert prov['download_url'].startswith('https://cdsarc')
        # every vendored data file carries an md5
        for fname, meta in prov['files'].items():
            assert (DATA_DIR / fname).exists()
            assert len(meta['md5']) == 32

    def test_fe2_table_not_vendored(self):
        # table4 (Fe II) is intentionally excluded — Fe NLTE is the MPIA grid.
        assert not (DATA_DIR / 'table4.dat').exists()
        assert not (DATA_DIR / 'table4.dat.gz').exists()


# ── Coverage / STOP-GATE ────────────────────────────────────────────────────────

class TestCoverage:
    def test_required_lines_and_metallicity(self):
        ok, problems = cno.coverage_ok()
        assert ok, f"coverage STOP-GATE failed: {problems}"

    def test_3d_ceiling_and_1d_reach(self):
        cov3 = cno.grid_coverage('OI', '3D')
        cov1 = cno.grid_coverage('OI', '1D')
        assert cov3['teff'][1] <= 6600       # 3D STAGGER ceiling ~6513
        assert cov1['teff'][1] >= 7500       # 1D leg reaches the F-dwarf regime
        assert cov3['feh'][1] >= 0.32        # covers 55 Cnc (+0.32)
        assert cov1['feh'][1] >= 0.32

    @pytest.mark.parametrize('star,teff,feh', [
        ('solar', 5772, 0.0), ('alpha_cen_a', 5792, 0.20), ('55cnc_a', 5196, 0.31),
    ])
    def test_benchmark_dwarfs_in_3d_box(self, star, teff, feh):
        v = cno.star_in_grid('OI', teff, 4.4, feh, 1.0)
        assert v['leg'] == '3D' and v['in_box']
        assert not v['required_missing']

    def test_procyon_routes_to_1d_leg(self):
        v = cno.star_in_grid('CI', 6554, 4.0, 0.03, 1.8)
        assert v['leg'] == '1D' and v['in_box']


# ── Leg selection ───────────────────────────────────────────────────────────────

def test_select_leg():
    assert cno.select_leg(5772) == '3D'
    assert cno.select_leg(6500) == '3D'
    assert cno.select_leg(6554) == '1D'   # Procyon above the 3D ceiling
    assert cno.select_leg(7000) == '1D'


# ── Line resolver (nm→Å + multiplet averaging) ──────────────────────────────────

class TestLineResolver:
    def test_optical_c_lines(self):
        assert cno.resolve_line('CI', 5052.16) == '505.2nm'
        assert cno.resolve_line('CI', 5380.34) == '538.0nm'

    def test_ir_c_lines(self):
        assert cno.resolve_line('CI', 9094.82) == '909.5nm'
        assert cno.resolve_line('CI', 9111.80) == '911.2nm'
        assert cno.resolve_line('CI', 9405.72) == '940.6nm'

    def test_oi_777_triplet_routes_to_multiplet_average(self):
        # all three 777 components → the single multiplet-averaged label
        for w in (7771.94, 7774.16, 7775.38):
            assert cno.resolve_line('OI', w) == '777nm'

    def test_oi_844_926_multiplets(self):
        assert cno.resolve_line('OI', 8446.4) == '844nm'
        assert cno.resolve_line('OI', 9262.6) == '926nm'

    def test_no_match_returns_none(self):
        assert cno.resolve_line('CI', 4000.0) is None

    def test_nm_to_angstrom_landed_right(self):
        # label '505.2nm' must resolve near 5052 Å, not 505 Å (unit catch)
        anchors = cno._load_line_table('CI')
        assert abs(anchors['505.2nm'] - 5052.16) < 1.0


# ── Species mapping ──────────────────────────────────────────────────────────────

def test_species_of():
    assert cno.species_of('C', 'I') == 'CI'
    assert cno.species_of('O', 'I') == 'OI'
    assert cno.species_of('C', 'II') is None     # only neutral tabulated
    assert cno.species_of('Fe', 'I') is None


# ── Correction values: sign + magnitude ─────────────────────────────────────────

class TestCorrectionSign:
    def test_solar_o777_large_negative(self):
        d = cno.cno_nlte_delta('OI', '777nm', 5772, 4.44, 0.0, 1.0, 8.69)
        assert np.isfinite(d) and -0.30 < d < -0.05   # large negative

    def test_solar_optical_c_small_negative(self):
        d = cno.cno_nlte_delta('CI', '505.2nm', 5772, 4.44, 0.0, 1.0, 8.46)
        assert np.isfinite(d) and -0.10 < d < 0.02    # near-LTE, small

    def test_procyon_o777_1d_leg_large_negative(self):
        d = cno.cno_nlte_delta('OI', '777nm', 6554, 4.0, 0.03, 1.8, 8.70)
        assert np.isfinite(d) and d < -0.2            # F dwarf: very negative

    def test_out_of_hull_returns_nan(self):
        # absurd logg far outside [3,5] (3D leg) → outside the convex hull
        d = cno.cno_nlte_delta('OI', '777nm', 5772, 0.0, 0.0, 1.0, 8.69, leg='3D')
        assert np.isnan(d)


class TestSignGuard:
    def test_guard_passes_on_negative(self):
        cno.assert_cno_sign('OI', '777nm', -0.18)        # no raise
        cno.assert_cno_sign('CI', '940.6nm', -0.06)

    def test_guard_raises_on_flipped_large_line(self):
        with pytest.raises(ValueError, match='sign flip|POSITIVE'):
            cno.assert_cno_sign('OI', '777nm', +0.18)
        with pytest.raises(ValueError, match='sign flip|POSITIVE'):
            cno.assert_cno_sign('CI', '909.5nm', +0.20)

    def test_guard_tolerates_tiny_positive_on_optical_c(self):
        cno.assert_cno_sign('CI', '505.2nm', +0.02)      # near-LTE, no raise


# ── table7 self-validation ───────────────────────────────────────────────────────

class TestTable7Validation:
    def test_sun_anchor_and_science_range(self):
        res = cno.validate_against_table7()
        # Sun row present with finite residuals
        sun = [r for r in res['rows'] if r['star'].lower() == 'sun'][0]
        assert np.isfinite(sun['o_resid']) and abs(sun['o_resid']) < 0.05
        assert np.isfinite(sun['c_resid']) and abs(sun['c_resid']) < 0.05
        # science-range pass
        assert res['passed'], (
            f"self-validation failed: C max(sci)={res['c_max_sci']}, "
            f"O max(sci)={res['o_max_sci']}")
        assert res['o_max_sci'] <= res['tol']
        assert res['c_max_sci'] <= res['tol']


# ── Public apply API ─────────────────────────────────────────────────────────────

class TestApply:
    def _solar_inputs(self):
        ab = pd.DataFrame([
            {'element': 'O', 'ion': 'I', 'A_X': 8.69},
            {'element': 'C', 'ion': 'I', 'A_X': 8.46},
            {'element': 'Fe', 'ion': 'I', 'A_X': 7.46},   # passthrough
        ])
        per_line = pd.DataFrame([
            {'element': 'O', 'ion': 'I', 'wavelength_air_A': 7771.94, 'a_1dlte': 8.69},
            {'element': 'O', 'ion': 'I', 'wavelength_air_A': 7774.16, 'a_1dlte': 8.70},
            {'element': 'C', 'ion': 'I', 'wavelength_air_A': 5052.16, 'a_1dlte': 8.46},
            {'element': 'C', 'ion': 'I', 'wavelength_air_A': 5380.34, 'a_1dlte': 8.45},
        ])
        params = {'teff_K': 5772, 'logg': 4.44, 'feh': 0.0, 'vturb_kms': 1.0}
        return ab, per_line, params

    def test_apply_solar(self):
        ab, per_line, params = self._solar_inputs()
        out = cno.apply_cno_nlte_corrections(ab, params, per_line_df=per_line)
        o = out[out['element'] == 'O'].iloc[0]
        c = out[out['element'] == 'C'].iloc[0]
        # O corrected, negative, flagged 3D
        assert o['n_nlte_lines'] == 2
        assert o['A_X_nlte'] < o['A_X']                 # negative correction
        assert o['nlte_flag'] == 'NLTE_Amarsi2019_3D'
        assert o['delta_nlte_mean'] < 0
        # C corrected
        assert c['n_nlte_lines'] == 2
        assert c['nlte_flag'] == 'NLTE_Amarsi2019_3D'

    def test_fe_row_passthrough(self):
        ab, per_line, params = self._solar_inputs()
        out = cno.apply_cno_nlte_corrections(ab, params, per_line_df=per_line)
        fe = out[out['element'] == 'Fe'].iloc[0]
        assert fe['A_X_nlte'] == fe['A_X']              # untouched
        assert fe['nlte_flag'] == '1D_LTE'

    def test_procyon_uses_1d_leg_flag(self):
        ab = pd.DataFrame([{'element': 'O', 'ion': 'I', 'A_X': 8.70}])
        per_line = pd.DataFrame([
            {'element': 'O', 'ion': 'I', 'wavelength_air_A': 7771.94, 'a_1dlte': 8.70},
        ])
        params = {'teff_K': 6554, 'logg': 4.0, 'feh': 0.03, 'vturb_kms': 1.8}
        out = cno.apply_cno_nlte_corrections(ab, params, per_line_df=per_line)
        assert out.iloc[0]['nlte_flag'] == 'NLTE_Amarsi2019_1D'

    def test_no_per_line_flags_unavailable(self):
        ab, _, params = self._solar_inputs()
        out = cno.apply_cno_nlte_corrections(ab, params, per_line_df=None)
        o = out[out['element'] == 'O'].iloc[0]
        assert o['nlte_flag'] == 'NLTE_unavailable'      # flagged, not silent NLTE
        assert o['A_X_nlte'] == o['A_X']
