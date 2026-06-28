"""
tests/test_solar_reference_rya459.py
====================================
RYA-459 (under RYA-162) — guard the solar reference library: registry/provenance
wiring, the cited-vs-measured discipline (the cardinal rule), the coverage of the
solar-N channels, and that the committed extracted segments are real solar spectra.

These run on the COMMITTED extracted CSVs + constants — they do NOT need the raw
external atlas (43 MB, outside the repo).
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from config.constants import SOLAR_REFERENCE_SPECTRA  # noqa: E402
from pipeline import audit_solar_reference as A        # noqa: E402

KP = ROOT / 'data' / 'solar_reference' / 'kpno_flux_atlas'
UV = ROOT / 'data' / 'solar_reference' / 'uv_composite'


def test_registry_has_tier1_sources():
    assert 'kpno_flux_atlas' in SOLAR_REFERENCE_SPECTRA
    assert 'uv_composite' in SOLAR_REFERENCE_SPECTRA
    kp = SOLAR_REFERENCE_SPECTRA['kpno_flux_atlas']
    assert kp['provenance'] == 'measured'
    assert kp['wavelength_coverage_nm'] == [296.0, 1300.0]


def test_uv_is_cited_composite_never_measured():
    # The cardinal rule: Hubble can't see the Sun → UV is cited, not measured.
    uv = SOLAR_REFERENCE_SPECTRA['uv_composite']
    assert uv['provenance'] == 'cited-composite'
    assert uv['provenance'] != 'measured'


def test_provenance_gate_passes():
    assert A.assert_provenance() is True


def test_provenance_gate_fails_if_uv_tagged_measured(monkeypatch):
    # Flip the UV flag to 'measured' and prove the gate catches it.
    bad = {k: dict(v) for k, v in SOLAR_REFERENCE_SPECTRA.items()}
    bad['uv_composite']['provenance'] = 'measured'
    monkeypatch.setattr(A, 'SOLAR_REFERENCE_SPECTRA', bad)
    assert A.assert_provenance() is False


def test_every_uv_row_carries_cited_flag():
    df = pd.read_csv(UV / 'sun_calspec_composite.csv')
    assert (df['provenance'] == 'cited-composite').all()
    assert (df['flux_erg_s_cm2_A'] > 0).all()


def test_solar_N_channels_are_measured_segments():
    # The RYA-369 unblock: all five solar-N channels staged from Kitt Peak (measured).
    for f in ('kpno_NH_3360.csv', 'kpno_CN_violet_3883.csv', 'kpno_NI_7442_7468.csv',
              'kpno_NI_8216_8223.csv', 'kpno_NI_8680_8718.csv'):
        df = pd.read_csv(KP / f)
        assert len(df) > 100
        assert {'wavelength_air_A', 'residual_flux'} <= set(df.columns)
        # real solar absorption present (line cores below continuum)
        assert df['residual_flux'].min() < 0.97


def test_oi_777_triplet_is_real():
    # Three O I lines near 7772/7774/7775 with cores well below continuum.
    df = pd.read_csv(KP / 'kpno_OI_777_triplet.csv')
    for line in (7771.94, 7774.17, 7775.39):
        core = df[(df['wavelength_air_A'] > line - 0.3) & (df['wavelength_air_A'] < line + 0.3)]
        assert core['residual_flux'].min() < 0.85, line


def test_coverage_matrix_all_diagnostics_measured():
    # Build the matrix from the committed provenance + assert full KP coverage.
    import json
    kp_prov = json.loads((KP / 'kpno_provenance_rya459.json').read_text())
    uv_prov = json.loads((UV / 'uv_provenance_rya459.json').read_text())
    uv_range_nm = [uv_prov['wl_vac_A_range'][0] / 10.0, uv_prov['wl_vac_A_range'][1] / 10.0]
    rows = A.coverage_matrix(kp_prov['extracted_diagnostics'], uv_range_nm)
    assert all(r['kitt_peak_measured'] for r in rows)        # 11/11 measured
    n_chan = [r for r in rows if r['element'] == 'N' and r['kitt_peak_measured']]
    assert len(n_chan) == 5                                   # the solar-N unblock


def test_irtf_tier2_is_deferred_not_faked():
    irtf = SOLAR_REFERENCE_SPECTRA['ir_atlas_irtf']
    assert irtf.get('status') == 'deferred'
    # placeholder documented, no fabricated spectrum
    readme = ROOT / 'data' / 'solar_reference' / 'ir_atlas' / 'README_SOURCE.md'
    assert readme.exists()
    assert not list((ROOT / 'data' / 'solar_reference' / 'ir_atlas').glob('*.csv'))
