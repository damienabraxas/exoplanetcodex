"""
tests/test_constants.py
=======================
Tests for config/constants.py.
"""

import pytest
from pathlib import Path

from config.constants import (
    PHYSICS, ASTRO, SOLAR_ASPLUND2021, STAR_55CNC,
    PIPELINE, PATHS, MODEL, validate_constants
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
