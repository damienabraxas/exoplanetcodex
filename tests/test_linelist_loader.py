"""
tests/test_linelist_loader.py
=============================
Tests for data/linelists/loader.py

Run with: pytest tests/test_linelist_loader.py -v
"""

import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import patch

from data.linelists.loader import load_linelist, get_element_lines, summarize_linelist


# ── Shared test data ──────────────────────────────────────────────────────────

COLUMNS = [
    'element', 'ion', 'wavelength_air_A', 'excitation_potential_eV',
    'log_gf', 'loggf_source', 'nist_grade', 'blend_flag', 'priority', 'notes',
]

SAMPLE_ROWS = [
    ('Fe', 'I',  5270.36, 1.608, -1.337, 'NIST', 'B',  False, 1, ''),
    ('Fe', 'I',  5328.04, 0.915, -1.466, 'NIST', 'A',  False, 1, ''),
    ('Fe', 'II', 5018.44, 2.891, -1.220, 'NIST', 'B',  False, 2, ''),
    ('Ca', 'I',  6122.22, 1.886, -0.316, 'NIST', 'A',  False, 1, ''),
    ('Ca', 'I',  6162.17, 1.899, -0.167, 'NIST', 'B',  False, 1, ''),
    ('Mg', 'I',  5183.60, 2.717, -0.239, 'NIST', 'A+', False, 1, ''),
    ('Na', 'I',  5889.95, 0.000,  0.108, 'NIST', 'A',  True,  1, 'IS blend'),
    ('O',  'I',  6300.30, 0.000, -9.717, 'NIST', 'C',  False, 2, 'Forbidden'),
]

PATCH_PIPELINE = {
    'wav_min_A': 3780.0,
    'wav_max_A': 6910.0,
    'min_nist_grade': 'B',
}


@pytest.fixture
def linelist_csv(tmp_path) -> Path:
    """Write the sample rows to a temporary CSV and return the path."""
    df = pd.DataFrame(SAMPLE_ROWS, columns=COLUMNS)
    path = tmp_path / 'linelist_master.csv'
    df.to_csv(path, index=False)
    return path


def _patches(linelist_csv):
    """Return a context manager that patches both PATHS and PIPELINE in the loader."""
    import contextlib

    @contextlib.contextmanager
    def _ctx():
        with patch('data.linelists.loader.PATHS', {'linelist_master': linelist_csv}), \
             patch('data.linelists.loader.PIPELINE', PATCH_PIPELINE):
            yield

    return _ctx()


# ── load_linelist ─────────────────────────────────────────────────────────────

class TestLoadLinelist:

    def test_returns_dataframe(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist()
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_has_expected_columns(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist(min_nist_grade='C', exclude_blends=False)
        for col in COLUMNS:
            assert col in result.columns, f"Missing column: {col}"

    def test_element_filtering_single(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist(elements=['Fe'])
        assert set(result['element'].unique()) == {'Fe'}
        assert len(result) > 0

    def test_element_filtering_multiple(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist(elements=['Fe', 'Ca'])
        assert set(result['element'].unique()) == {'Fe', 'Ca'}

    def test_element_filtering_excludes_others(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist(elements=['Mg'])
        assert 'Fe' not in result['element'].values
        assert 'Ca' not in result['element'].values

    def test_ion_filtering_neutral(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist(ion='I')
        assert set(result['ion'].unique()) == {'I'}

    def test_ion_filtering_ionized(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist(ion='II')
        assert set(result['ion'].unique()) == {'II'}

    def test_grade_filter_excludes_lower_grades(self, linelist_csv):
        """min_nist_grade='A' should exclude B, C, D, E."""
        with _patches(linelist_csv):
            result = load_linelist(min_nist_grade='A')
        assert 'B' not in result['nist_grade'].values
        assert 'C' not in result['nist_grade'].values

    def test_grade_filter_includes_better_grades(self, linelist_csv):
        """min_nist_grade='B' should include A+ and A lines."""
        with _patches(linelist_csv):
            result = load_linelist(min_nist_grade='B')
        assert 'A+' in result['nist_grade'].values
        assert 'A' in result['nist_grade'].values

    def test_grade_filter_c_includes_all(self, linelist_csv):
        """min_nist_grade='C' should include A+, A, B, and C."""
        with _patches(linelist_csv):
            result = load_linelist(min_nist_grade='C', exclude_blends=False)
        assert set(result['nist_grade'].unique()) >= {'A+', 'A', 'B', 'C'}

    def test_exclude_blends_true(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist(exclude_blends=True)
        assert result['blend_flag'].sum() == 0

    def test_exclude_blends_false_includes_blended(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist(exclude_blends=False, min_nist_grade='C')
        assert result['blend_flag'].sum() > 0

    def test_priority_filter(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist(priority=1, min_nist_grade='C', exclude_blends=False)
        assert (result['priority'] <= 1).all()

    def test_wavelength_filter(self, linelist_csv):
        with _patches(linelist_csv):
            result = load_linelist(wav_min=5000.0, wav_max=5500.0)
        assert result['wavelength_air_A'].between(5000.0, 5500.0).all()

    def test_wavelength_filter_excludes_out_of_range(self, linelist_csv):
        """Ca lines at 6122/6162 should not appear when wav_max=5500."""
        with _patches(linelist_csv):
            result = load_linelist(wav_max=5500.0)
        assert (result['wavelength_air_A'] <= 5500.0).all()

    def test_returns_reset_index(self, linelist_csv):
        """Index should be 0-based after filtering."""
        with _patches(linelist_csv):
            result = load_linelist(elements=['Fe'])
        assert list(result.index) == list(range(len(result)))

    def test_file_not_found(self, tmp_path):
        missing = tmp_path / 'does_not_exist.csv'
        with patch('data.linelists.loader.PATHS', {'linelist_master': missing}), \
             patch('data.linelists.loader.PIPELINE', PATCH_PIPELINE):
            with pytest.raises(FileNotFoundError, match='linelist_master'):
                load_linelist()

    def test_no_elements_returns_all(self, linelist_csv):
        """elements=None should return lines for all elements (subject to grade filter)."""
        with _patches(linelist_csv):
            result = load_linelist(min_nist_grade='C', exclude_blends=False)
        assert len(result['element'].unique()) > 1


# ── get_element_lines ─────────────────────────────────────────────────────────

class TestGetElementLines:

    def test_matches_load_linelist(self, linelist_csv):
        """get_element_lines should return same result as load_linelist with matching args."""
        with _patches(linelist_csv):
            via_load = load_linelist(elements=['Ca'], ion='I')
            via_func = get_element_lines('Ca', ion='I')
        pd.testing.assert_frame_equal(
            via_load.reset_index(drop=True),
            via_func.reset_index(drop=True),
        )

    def test_defaults_to_neutral(self, linelist_csv):
        """Default ion should be 'I' (neutral)."""
        with _patches(linelist_csv):
            result = get_element_lines('Fe')
        assert set(result['ion'].unique()) == {'I'}

    def test_kwargs_pass_through(self, linelist_csv):
        """Extra kwargs like min_nist_grade should be forwarded."""
        with _patches(linelist_csv):
            result = get_element_lines('Fe', ion='I', min_nist_grade='A')
        assert 'B' not in result['nist_grade'].values


# ── summarize_linelist ────────────────────────────────────────────────────────

class TestSummarizeLinelist:

    def test_prints_summary(self, linelist_csv, capsys):
        with _patches(linelist_csv):
            df = load_linelist(min_nist_grade='C', exclude_blends=False)
            summarize_linelist(df)
        out = capsys.readouterr().out
        assert 'Total lines' in out
        assert 'Elements' in out
        assert 'NIST grades' in out
        assert 'Priorities' in out

    def test_does_not_raise(self, linelist_csv):
        with _patches(linelist_csv):
            df = load_linelist(min_nist_grade='C', exclude_blends=False)
        summarize_linelist(df)  # should not raise
