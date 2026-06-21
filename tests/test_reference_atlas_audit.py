"""
tests/test_reference_atlas_audit.py
===================================
RYA-392 — guard the read-only IR reference-atlas audit: axis classification (no
silent UNKNOWN), cm^-1 <-> wavelength conversion, comma/whitespace parsing, the
ACE telluric-free-vs-occultation gate (incl. the "wallACE" substring trap), the
endpoint coverage tolerance, and the real-store integration verdict (GO).
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline.audits import audit_reference_atlases as a  # noqa: E402

STORE = ROOT / 'data' / 'solar_reference' / 'ir_atlases'


def test_sniff_xaxis_wavenumber_vs_wavelength():
    assert a.sniff_xaxis(np.array([4255.0, 4367.0]))[0] == 'wavenumber_cm'
    assert a.sniff_xaxis(np.array([2.29, 2.35]))[0] == 'wavelength_um'
    assert a.sniff_xaxis(np.array([22893.0, 23495.0]))[0] == 'wavelength_AA'


def test_sniff_xaxis_unknown_is_loud():
    # a value range that matches nothing must NOT be silently guessed
    assert a.sniff_xaxis(np.array([40.0, 60.0]))[0] == 'UNKNOWN'


def test_to_cm_roundtrip():
    assert a.to_cm(np.array([4360.0]), 'wavenumber_cm')[0] == pytest.approx(4360.0)
    # 4360 cm^-1 -> 2.2936 um -> back to 4360 cm^-1
    um = 1.0e4 / 4360.0
    assert a.to_cm(np.array([um]), 'wavelength_um')[0] == pytest.approx(4360.0)
    aa = 1.0e8 / 4360.0
    assert a.to_cm(np.array([aa]), 'wavelength_AA')[0] == pytest.approx(4360.0)


def test_tokens_handles_comma_and_whitespace():
    assert a._tokens('4255.0,23501.7,0.994') == ['4255.0', '23501.7', '0.994']
    assert a._tokens('4255.0  23501.7  0.994') == ['4255.0', '23501.7', '0.994']


def test_read_any_parses_csv(tmp_path):
    p = tmp_path / 'mini.csv'
    p.write_text('wavenumber_cm-1,wavelength_vac_A,intensity\n'
                 '4255.0,23501.76,0.99\n4256.0,23496.24,0.98\n')
    cols, _ = a.read_any(str(p))
    assert list(cols['x']) == [4255.0, 4256.0]          # not collapsed by comma-strip
    assert 'wavenumber_cm-1' in cols and 'intensity' in cols


def test_ace_truth_source_accepts_derived_rejects_occultation():
    derived = dict(name='ACE-FTS solar atlas',
                   citation='Hase+2010 — A complete solar spectrum based on ACE data',
                   source_url='https://databace.scisat.ca/solarspectrum/',
                   raw_file='ace-solar-spectrum.txt',
                   telluric_status='telluric-FREE (space occultation)')
    ok, verdict, _ = a._ace_truth_source(derived)
    assert ok and 'DERIVED' in verdict

    occ = dict(derived, raw_file='ace_occultation_transmission.txt',
               telluric_status='raw occultation product')
    ok2, verdict2, _ = a._ace_truth_source(occ)
    assert not ok2 and 'OCCULTATION' in verdict2


def test_wallace_does_not_trip_the_ace_gate():
    # "wallace" contains the substring "ace" — the gate must key off the atlas name,
    # not a naive substring, or Wallace would be flagged RAW-OCCULTATION SUSPECT.
    r = a.audit_one(str(STORE / 'wallace_telluric_co_ratio.csv'), _prov())
    assert r['ok']
    assert not any('OCCULTATION' in n for n in r['notes'])


def test_coverage_tolerance_absorbs_sampling_grid():
    # photatl samples 4255.007..4366.998 — inside the FTS step of the exact band edges
    r = a.audit_one(str(STORE / 'nso_photatl_co_4255_4367.csv'), _prov())
    assert r['covers_full'] and r['covers_ok']


def _prov():
    by_file, _ = a.load_provenance(str(STORE))
    return by_file


def test_real_store_verdict_is_go():
    assert a.main(['--store', str(STORE)]) == 0


def test_missing_store_is_nogo(tmp_path):
    assert a.main(['--store', str(tmp_path)]) == 2
