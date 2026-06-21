"""
tests/test_ir_atlases_rya390.py
===============================
RYA-390 — guard the IR reference-atlas intake: vac→air conversion + the three
extracted CO segments (columns, CO-band coverage, physical ranges).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))
import intake_ir_atlases_rya390 as intake  # noqa: E402

SEG = ROOT / 'data' / 'solar_reference' / 'ir_atlases'


def test_vac_to_air_birch_downs():
    # air < vac, and the index of refraction n−1 ≈ 2.7e-4 in the near-IR (2.3 µm)
    vac = np.array([23501.76])              # Å (≈ 4255 cm⁻¹, vacuum)
    air = intake.vac_to_air_A(vac)
    assert air[0] < vac[0]
    n_minus_1 = vac[0] / air[0] - 1.0
    assert 2.5e-4 < n_minus_1 < 3.0e-4
    # H-alpha anchor (well-known): vac 6564.61 Å → air 6562.81 Å (Birch & Downs)
    ha = intake.vac_to_air_A(np.array([6564.61]))[0]
    assert abs(ha - 6562.81) < 0.05


@pytest.mark.skipif(not (SEG / 'ace_fts_solar_co_4255_4367.csv').exists(),
                    reason="CO segments not intaked")
@pytest.mark.parametrize('fname,cols,role', [
    ('ace_fts_solar_co_4255_4367.csv',
     {'wavenumber_cm-1', 'wavelength_vac_A', 'wavelength_air_A', 'intensity'}, 'solar'),
    ('nso_photatl_co_4255_4367.csv',
     {'wavenumber_cm-1', 'wavelength_vac_A', 'wavelength_air_A',
      'solar', 'atmospheric', 'total'}, 'solar'),
    ('wallace_telluric_co_ratio.csv',
     {'wavenumber_cm-1', 'wavelength_vac_A', 'wavelength_air_A', 'telluric_ratio'}, 'telluric'),
])
def test_segment_schema_and_band(fname, cols, role):
    df = pd.read_csv(SEG / fname)
    assert set(df.columns) == cols
    assert len(df) > 1000
    wn = df['wavenumber_cm-1']
    # within the CO band, monotonic, air < vac
    assert wn.min() >= intake.CO_LO_CM - 1 and wn.max() <= intake.CO_HI_CM + 1
    assert (df['wavelength_air_A'] < df['wavelength_vac_A']).all()
    # physical intensity/transmission range
    val = df['intensity'] if role == 'solar' and 'intensity' in df else (
        df['telluric_ratio'] if role == 'telluric' else df['solar'])
    assert val.min() >= -0.05 and val.max() <= 1.1


def test_provenance_manifest_complete():
    import json
    prov = json.loads((SEG / 'ir_atlases_provenance_rya390.json').read_text())
    assert len(prov['sources']) == 3
    statuses = ' '.join(s['telluric_status'] for s in prov['sources']).lower()
    assert 'telluric-free' in statuses and 'pure telluric' in statuses
    for s in prov['sources']:
        assert s['citation'] and s['source_url'] and s['raw_sha256']
