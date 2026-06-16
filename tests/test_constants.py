"""
tests/test_constants.py
=======================
Tests for config/constants.py.
"""

import pytest
from pathlib import Path

from config.constants import (
    PHYSICS, ASTRO, SOLAR_ASPLUND2021, STAR_55CNC,
    PIPELINE, PATHS, MODEL, validate_constants,
    STAR_PARAMS, get_star_params,
)


class TestPhysics:
    def test_speed_of_light_kms(self):
        assert PHYSICS['c_kms'] == pytest.approx(299792.458, rel=1e-6)

    def test_speed_of_light_ms(self):
        assert PHYSICS['c_ms'] == pytest.approx(299792458.0, rel=1e-6)

    def test_planck_constant_eV(self):
        assert PHYSICS['h_eV'] == pytest.approx(4.135667696e-15, rel=1e-6)

    def test_required_keys_present(self):
        for key in ('c_kms', 'c_ms', 'h_eV'):
            assert key in PHYSICS


class TestAstro:
    def test_solar_mass(self):
        assert ASTRO['Msun_g'] == pytest.approx(1.989e33, rel=1e-3)

    def test_solar_teff(self):
        assert ASTRO['Teff_sun'] == pytest.approx(5778.0, rel=1e-3)

    def test_solar_logg(self):
        assert ASTRO['logg_sun'] == pytest.approx(4.438, rel=1e-3)


class TestSolarAbundances:
    def test_iron_asplund2021(self):
        # Fe = 7.46 (Asplund 2021). Note: older Asplund 2009 value is 7.50.
        assert SOLAR_ASPLUND2021['Fe'] == pytest.approx(7.46, rel=1e-3)

    def test_magnesium(self):
        assert SOLAR_ASPLUND2021['Mg'] == pytest.approx(7.55, rel=1e-3)

    def test_oxygen(self):
        assert SOLAR_ASPLUND2021['O'] == pytest.approx(8.69, rel=1e-3)

    def test_all_expected_elements_present(self):
        expected = {'Fe', 'Mg', 'Si', 'Ca', 'O', 'Ni', 'Na', 'Al', 'Ti', 'Cr'}
        missing = expected - set(SOLAR_ASPLUND2021.keys())
        assert not missing, f"Missing elements: {missing}"

    def test_abundance_values_in_plausible_range(self):
        # A(X) = log N(X) + 12 with A(H) ≡ 12.00 by definition (the normalization
        # anchor), so 12 is the inclusive ceiling — nothing should exceed it, but
        # H sits exactly on it. Upper bound is therefore <= 12, not < 12 (RYA-295).
        for element, abundance in SOLAR_ASPLUND2021.items():
            assert 0 < abundance <= 12.0, (
                f"{element} abundance {abundance} outside plausible range (0 < A(X) <= 12 dex)"
            )


class TestAbundanceScaleGuard:
    """RYA-334: the range-sanity tripwire for the relative/absolute double-add."""

    def test_passes_normal_absolute_values(self):
        import math
        from config.constants import assert_abundance_on_scale
        for v in (7.46, 7.516, 12.0, 0.52, -1.0):   # Fe, MPIA solar, H, Eu, faint
            assert assert_abundance_on_scale(v, 'test') == v

    def test_nan_and_none_pass_through(self):
        import math
        from config.constants import assert_abundance_on_scale
        assert math.isnan(assert_abundance_on_scale(float('nan'), 'test'))
        assert assert_abundance_on_scale(None, 'test') is None

    def test_raises_on_double_add(self):
        # A 7.46 solar offset re-added to an already-absolute ~7.5 lands ~15 — the
        # exact RYA-267/320/330 signature. Must fail loud, not pass.
        from config.constants import assert_abundance_on_scale
        with pytest.raises(ValueError, match="off the A.H.=12 scale"):
            assert_abundance_on_scale(7.495 + 7.46, 'Fe I gate (simulated double-add)')

    def test_raises_below_floor(self):
        from config.constants import assert_abundance_on_scale
        with pytest.raises(ValueError):
            assert_abundance_on_scale(-5.0, 'test')


class TestScaleAwareFeGate:
    """RYA-336: absolute Fe diagnostic centred on 3D-true + published 1D-3D offset."""

    def test_diagnostic_centre_and_window(self):
        from config.constants import (
            SOLAR_ASPLUND2021, FE_1D3D_SOLAR_OFFSET, FE_ABS_DIAG_HALFWIDTH)
        centre = SOLAR_ASPLUND2021['Fe'] + FE_1D3D_SOLAR_OFFSET
        assert centre == pytest.approx(7.51, abs=1e-9)          # 7.46 + 0.05
        lo, hi = centre - FE_ABS_DIAG_HALFWIDTH, centre + FE_ABS_DIAG_HALFWIDTH
        assert (lo, hi) == pytest.approx((7.44, 7.58), abs=1e-9)
        # both grids straddle the centre and fall inside the window
        for grid_solar in (7.495, 7.516):                       # Amarsi, MPIA
            assert lo <= grid_solar <= hi
        # a gross zero-point error (e.g. loggf slip ~+0.3) must fall outside
        assert not (lo <= 7.46 + 0.30 <= hi)

    def test_offset_is_independent_of_our_output(self):
        # The offset is a fixed published quantity, NOT read from any run result.
        from config.constants import FE_1D3D_SOLAR_OFFSET
        assert FE_1D3D_SOLAR_OFFSET == pytest.approx(0.05, abs=1e-9)


class TestStar55Cnc:
    def test_teff(self):
        assert STAR_55CNC['teff_K'] == pytest.approx(5196.0, rel=1e-3)

    def test_logg(self):
        assert STAR_55CNC['logg'] == pytest.approx(4.41, rel=1e-3)

    def test_feh(self):
        assert STAR_55CNC['feh'] == pytest.approx(0.32, rel=1e-2)

    def test_rv_kms(self):
        # Positive RV = receding. Important sign convention for Doppler correction.
        assert STAR_55CNC['rv_kms'] > 0
        assert STAR_55CNC['rv_kms'] == pytest.approx(27.58, rel=1e-2)

    def test_hd_identifier(self):
        assert STAR_55CNC['hd'] == 'HD 75732'

    def test_hip_identifier(self):
        assert STAR_55CNC['hip'] == 'HIP 43587'


class TestPipeline:
    def test_wavelength_range(self):
        assert PIPELINE['wav_min_A'] < PIPELINE['wav_max_A']
        assert PIPELINE['wav_min_A'] == pytest.approx(3780.0)
        assert PIPELINE['wav_max_A'] == pytest.approx(6910.0)

    def test_ew_range(self):
        assert PIPELINE['ew_min_mA'] < PIPELINE['ew_max_mA']

    def test_snr_min(self):
        assert PIPELINE['snr_min_science'] >= 100

    def test_min_nist_grade(self):
        assert PIPELINE['min_nist_grade'] in ('A+', 'A', 'B', 'C', 'D', 'E')


class TestPaths:
    def test_paths_are_path_objects(self):
        for key, val in PATHS.items():
            assert isinstance(val, Path), f"PATHS['{key}'] is not a Path object"

    def test_key_paths_present(self):
        for key in ('raw_spectra', 'plots', 'linelist_master'):
            assert key in PATHS


class TestModel:
    def test_model_type(self):
        assert MODEL['type'] == 'ATLAS9'

    def test_geometry(self):
        assert MODEL['geometry'] == 'plane-parallel'

    def test_radiative_transfer(self):
        assert MODEL['radiative_transfer'] == 'LTE'


class TestValidateConstants:
    def test_validate_passes(self):
        """validate_constants() should not raise for a well-formed constants module."""
        validate_constants()  # raises if any checks fail

    def test_speed_of_light_consistency(self):
        # c_ms should equal c_kms * 1000 to within floating-point tolerance
        assert PHYSICS['c_ms'] == pytest.approx(PHYSICS['c_kms'] * 1000, rel=1e-9)


class TestStarParamsPolicy:
    """Per-star fundamental-param record + pin/solve policy (RYA-292 / RYA-293)."""

    def test_55cnc_a_pins_teff_and_logg(self):
        # RYA-293: 55 Cnc A graduates to a pinned target (von Braun et al. 2011).
        rec = get_star_params('55cnc_a')
        assert rec['teff'] == pytest.approx(5196.0, rel=1e-3)
        assert rec['logg'] == pytest.approx(4.45, rel=1e-3)
        assert rec['pin'] == ['teff', 'logg']
        assert 'feh' in rec['solve'] and 'xi' in rec['solve']
        # [Fe/H] must stay solved (spectroscopic), not pinned — it is a sanity target.
        assert 'feh' not in rec['pin']

    def test_55cnc_a_provenance_is_honest_about_isochrone_mass(self):
        # The mass is model-dependent — the basis string must not imply a dynamical mass.
        rec = get_star_params('55cnc_a')
        assert 'von Braun' in rec['source']
        assert 'isochrone' in rec['logg_basis'].lower()
        assert 'interferometric' in rec['logg_basis'].lower()

    def test_solve_path_still_covered(self):
        # RYA-293 Step 2: after 55 Cnc A became pinned, a no-fundamental-logg star
        # must still exercise the spectroscopic-solve path (logg solved, pin empty).
        rec = get_star_params('synthetic_no_logg')
        assert rec['pin'] == []
        assert 'logg' in rec['solve']
        assert 'teff' in rec['solve']

    def test_benchmarks_pin_teff_and_logg(self):
        for star in ('solar', 'procyon', 'alpha_cen_a', 'alpha_cen_b'):
            rec = get_star_params(star)
            assert 'teff' in rec['pin'] and 'logg' in rec['pin'], star

    def test_missing_record_raises(self):
        # No silent default — fail loud for an unknown star (guard intact).
        with pytest.raises(KeyError):
            get_star_params('kepler-10')

    def test_every_record_has_pin_and_solve(self):
        for key, rec in STAR_PARAMS.items():
            assert 'pin' in rec and 'solve' in rec, key
            assert isinstance(rec['pin'], list) and isinstance(rec['solve'], list), key
