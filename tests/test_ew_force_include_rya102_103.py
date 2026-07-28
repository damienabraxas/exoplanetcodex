"""
tests/test_ew_force_include_rya102_103.py
=========================================
RYA-371 Phase 1 — guard the RYA-102/103 fold-in: the abundance EW pre-filter must
FORCE-INCLUDE Eu II 6645 (RYA-102, weak HFS-total) and Li I 6707 (RYA-103, CN-blend
upper-limit), which are each their element's only solar indicator and are blend-
flagged (Li also below the 5 mA floor) — without that exemption they vanish into a
DATA-GAP. The exemption must be additive: every non-force-include line is filtered
exactly as before (Fe/metal pools unchanged).
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline.abundances_derive import _ew_prefilter, EW_FORCE_INCLUDE  # noqa: E402

EW_MIN, EW_MAX = 5.0, 300.0


def _df():
    return pd.DataFrame([
        # element, ion, wave, ew_mA, blend_flag, notes
        ('Fe', 'I',  5576.089, 60.0, False, ''),                 # clean -> kept
        ('Fe', 'I',  6065.482, 40.0, True,  'blend'),            # blended non-special -> DROPPED
        ('Si', 'I',  5772.146,  2.0, False, ''),                 # below ew_min, non-special -> DROPPED
        ('Eu', 'II', 6645.127,  6.80, True, 'blend_flag; hfs_total_ew'),   # RYA-102 -> force-kept
        ('Li', 'I',  6707.840,  1.77, True, 'upper_limit; CN_blend_possible'),  # RYA-103 (blend AND <5) -> force-kept
    ], columns=['element', 'ion', 'wavelength_air_A', 'ew_mA', 'blend_flag', 'notes'])


def test_eu_li_force_included_with_flags():
    clean, n_forced = _ew_prefilter(_df(), EW_MIN, EW_MAX)
    el = set(clean['element'])
    assert 'Eu' in el and 'Li' in el          # both kept despite blend_flag (Li also <5)
    assert n_forced == 2
    li = clean[clean['element'] == 'Li'].iloc[0]
    assert 'upper_limit' in li['notes']       # flag preserved for the downstream verdict


def test_additive_non_special_lines_filtered_as_before():
    clean, _ = _ew_prefilter(_df(), EW_MIN, EW_MAX)
    el = set(clean['element'])
    assert 'Fe' in el                          # clean Fe kept
    # the blended Fe (non-special) and the weak Si (non-special) are still dropped
    assert not ((clean['element'] == 'Fe') & (clean['wavelength_air_A'] == 6065.482)).any()
    assert 'Si' not in el


def test_force_include_set_is_rya_102_103():
    assert ('Eu', 'II', 6645.127) in EW_FORCE_INCLUDE   # RYA-102
    assert ('Li', 'I', 6707.840) in EW_FORCE_INCLUDE    # RYA-103


def test_empty_df_safe():
    clean, n = _ew_prefilter(pd.DataFrame(columns=['element', 'ion', 'wavelength_air_A',
                                                   'ew_mA', 'blend_flag']), EW_MIN, EW_MAX)
    assert len(clean) == 0 and n == 0
