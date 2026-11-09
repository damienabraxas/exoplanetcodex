"""Nonblind source accounting must not masquerade as laboratory validation."""
import json
from pathlib import Path
import sys
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_reference_matrix import build, triplet_reference
from rya1354_native_summary import main as native_summary


def test_native_receipts_replay_and_keep_failed_ladder(tmp_path):
    output=tmp_path/'native.json'
    native_summary(output)
    assert output.read_bytes()==(ROOT/'data/audit/rya1354/nonblind_v1/native_summary.json').read_bytes()
    d=json.loads(output.read_text())
    assert d['receipts']['runs/valence_sd_v5']['status']=='FAILED_RETAINED'
    assert 'n5' not in d['ladder']
    assert len(d['coverage']['IR_3d_E1'])==9
    assert d['forbidden']['status']=='NATIVE_M1_AND_E2_COMPUTED_NOT_ADOPTED'


def test_matrix_replays_and_preserves_missing_and_ambiguous_rows(tmp_path):
    output=tmp_path/'matrix.json'
    build(destination=output)
    held=ROOT/'data/audit/rya1354/nonblind_v1/reference_matrix.json'
    assert output.read_bytes()==held.read_bytes()
    d=json.loads(output.read_text())
    assert len(d['rows'])==26
    assert sum(r['prediction_status']=='HOLD_AMBIGUOUS_GRID_COMPONENT_IDENTITY' for r in d['rows'])==9
    assert sum('exposed_grid_loggf_residual_dex' in r for r in d['rows'])==15
    assert all(r['experimental_gf_accuracy_score'] is None for r in d['rows'])
    forbidden=next(r for r in d['rows'] if r['source_grid']['wavelength_air_A']=='6300.300')
    assert {r['Type'] for r in forbidden['nist_candidates']}=={'M1','E2'}
    assert forbidden['uncertainty'] is None


def test_wrong_term_and_duplicate_identity_refused_even_at_same_wavelength():
    r={'Lower level':'2s2.2p3.(4S*).3s | 5S* | 2',
       'Upper level':'2s2.2p3.(4S*).3p | 5P | 1','Ritz':'7777.527'}
    assert triplet_reference([r]) is r
    bad=dict(r,**{'Upper level':'2s2.2p3.(4S*).3p | 3P | 1'})
    with pytest.raises(ValueError):
        triplet_reference([bad])
    with pytest.raises(ValueError):
        triplet_reference([r,r.copy()])
