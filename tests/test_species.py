"""
tests/test_species.py
======================
Tests for pipeline/species.py — the canonical species/ion encoding normalizer
(RYA-345). Every inventoried encoding must round-trip to the same canonical key,
and unresolvable / inconsistent ions must RAISE (no silent mis-keying).
"""
import pytest

from pipeline.species import (
    species_key, species_note, element_z, z_symbol, parse_ion, MOLECULE,
)

FE_II = (26, 2)
FE_I = (26, 1)


class TestCanonicalKeyAllEncodings:
    """The RYA-345 acceptance set: every Fe II encoding → (26, 2)."""

    def test_ges_combined_string(self):
        assert species_key('Fe 2') == FE_II            # GES element col

    def test_combined_string_plus_redundant_int(self):
        assert species_key('Fe 2', 2) == FE_II         # ('Fe 2', 2)

    def test_symbol_plus_roman(self):
        assert species_key('Fe', 'II') == FE_II         # linelist_solar / solar_ew

    def test_symbol_plus_int(self):
        assert species_key('Fe', 2) == FE_II            # GES regions int32 ion

    def test_combined_roman_string(self):
        assert species_key('Fe II') == FE_II            # MPIA-style roman embedded

    def test_moog_code_one_digit(self):
        assert species_key('26.1') == FE_II             # MOOG species code

    def test_moog_code_two_digit(self):
        assert species_key('26.01') == FE_II            # 2-digit variant

    def test_numeric_z_plus_ion(self):
        assert species_key(26, 2) == FE_II              # (Z, ion)

    def test_ispec_region_note(self):
        assert species_key('Fe 1') == FE_I              # iSpec 'note' field

    def test_all_fe2_encodings_identical(self):
        keys = {
            species_key('Fe 2'),
            species_key('Fe 2', 2),
            species_key('Fe', 'II'),
            species_key('Fe', 2),
            species_key('Fe II'),
            species_key('26.1'),
            species_key('26.01'),
            species_key(26, 2),
        }
        assert keys == {FE_II}


class TestNeutralAndOther:
    def test_neutral_codes(self):
        assert species_key('Fe 1') == FE_I
        assert species_key('Fe', 'I') == FE_I
        assert species_key('26.0') == FE_I
        assert species_key('26.00') == FE_I
        assert species_key(26, 1) == FE_I

    def test_other_species(self):
        assert species_key('Ca', 'II') == (20, 2)
        assert species_key('Ti 1') == (22, 1)
        assert species_key('22.0') == (22, 1)
        assert species_key('8.0') == (8, 1)       # O I


class TestMolecules:
    def test_molecule_flag(self):
        assert species_key('CH', molecule='T') == (MOLECULE, 'CH')

    def test_molecule_flag_false_is_atomic(self):
        # 'F' flag → not molecular; still an atom if a valid symbol
        assert species_key('Fe 2', molecule='F') == FE_II

    def test_unknown_symbol_is_molecule(self):
        assert species_key('MgH 1') == (MOLECULE, 'MgH')
        assert species_key('CN') == (MOLECULE, 'CN')

    def test_molecular_numeric_code(self):
        assert species_key('822.0')[0] == MOLECULE     # TiO-style code, Z>99


class TestLoudFailures:
    """No silent fallback — bad input raises, never mis-keys to a default."""

    def test_unresolvable_ion_raises(self):
        with pytest.raises(ValueError):
            species_key('Fe')              # bare symbol, no ion anywhere

    def test_numeric_code_without_ion_raises(self):
        with pytest.raises(ValueError):
            species_key('26')              # bare Z, no ionization

    def test_inconsistent_embedded_vs_explicit_raises(self):
        with pytest.raises(ValueError):
            species_key('Fe 2', 1)         # embedded II vs explicit I

    def test_inconsistent_code_vs_explicit_raises(self):
        with pytest.raises(ValueError):
            species_key('26.1', 1)         # code II vs explicit I

    def test_unknown_symbol_with_ion_resolves_molecule_not_raise(self):
        # 'Xx' is not an atom; with no recognizable element it falls to molecule
        assert species_key('Xx')[0] == MOLECULE

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            species_key('')

    def test_bad_ion_token_raises(self):
        with pytest.raises(ValueError):
            parse_ion('banana')

    def test_zero_ion_raises(self):
        with pytest.raises(ValueError):
            parse_ion(0)


class TestParseIon:
    def test_roman(self):
        assert parse_ion('I') == 1
        assert parse_ion('II') == 2
        assert parse_ion('III') == 3

    def test_int_and_digit_string(self):
        assert parse_ion(2) == 2
        assert parse_ion('2') == 2

    def test_lowercase_roman(self):
        assert parse_ion('ii') == 2


class TestElementTable:
    def test_symbol_to_z(self):
        assert element_z('Fe') == 26
        assert element_z('H') == 1
        assert element_z('O') == 8

    def test_z_to_symbol(self):
        assert z_symbol(26) == 'Fe'
        assert z_symbol(1) == 'H'

    def test_roundtrip(self):
        for sym in ('Fe', 'Ca', 'Ti', 'Ni', 'O', 'Mg'):
            assert z_symbol(element_z(sym)) == sym


class TestSpeciesNote:
    """species_note renders the iSpec region note and fixes the old ion>II bug."""

    def test_renders_note(self):
        assert species_note('Fe', 'II') == 'Fe 2'
        assert species_note('Fe 2') == 'Fe 2'
        assert species_note('26.0') == 'Fe 1'

    def test_higher_ion_not_collapsed_to_2(self):
        # old _ours_to_ispec_note returned 2 for anything != 'I'; here III → 3
        assert species_note('Fe', 'III') == 'Fe 3'


class TestZeroMatchGuard:
    """RYA-345 Step 4: a species with measured lines but 0 matched regions must
    fail LOUD (warn), not silently look like 'no such line'."""

    def _fake_regions(self):
        import numpy as np
        dt = np.dtype([('note', 'U10'), ('wave_A', 'f8'),
                       ('ew', 'f8'), ('ew_err', 'f8')])
        # only an Fe I region exists; no Fe II region anywhere
        return np.array([('Fe 1', 5000.0, 0.0, 0.0)], dtype=dt)

    def test_zero_match_species_warns(self, monkeypatch):
        import warnings
        import numpy as np
        import pandas as pd
        import pipeline.abundances_derive as ad

        monkeypatch.setattr(ad.ispec, 'read_line_regions',
                            lambda _p: self._fake_regions())
        ew_df = pd.DataFrame({
            'element': ['Fe', 'Fe'],
            'ion':     ['I', 'II'],         # Fe I matches; Fe II has no region
            'wavelength_air_A': [5000.0, 6000.0],
            'ew_mA':   [50.0, 40.0],
        })
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            lm = ad._build_ispec_line_regions(ew_df, line_regions_path='x')
        msgs = [str(x.message) for x in w if issubclass(x.category, RuntimeWarning)]
        assert any('0-match guard' in m and '(26, 2)' in m for m in msgs), msgs
        assert len(lm) == 1            # the Fe I line still matched (no regression)
