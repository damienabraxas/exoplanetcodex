"""Native frozen-layer replay and adversarial source admission checks."""
import hashlib
import json
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_frozen_layer import validate,PARENT


@pytest.fixture(scope='module')
def native():
    base=ROOT/'data/audit/rya1354'
    return load_parts(base/'from_n3_terminal_v12')[0],load_parts(base/'lifetime_aux_v4b')[0]


def test_all_six_native_states_are_restricted_not_full_relaxation(native):
    report=validate(*native)
    assert len(report['states'])==6
    assert report['parent_status']=='FAILED_RETAINED'
    assert not report['fully_released_qualified'] and not report['adoption_allowed']
    assert report['accepted_lifetime'] is report['uncertainty'] is report['calibrated_sigma'] is None
    assert all(s['exact_frozen_n_le_3_records']==9 for s in report['states'].values())


@pytest.mark.parametrize('case',['outer_n4_quintet','outer_n4_triplet_verify'])
def test_changed_odd_state_is_refused_even_with_updated_receipt(native,case):
    files,target=native;files=dict(files)
    key=case+'/rmcdhf.stdin'
    files[PARENT+key]=files[PARENT+key].replace(b'4*\n',b'2*,3*,4*\n')
    receipt=json.loads(files[PARENT+'receipt.json'])
    receipt['files'][key]=dict(sha256=hashlib.sha256(files[PARENT+key]).hexdigest(),bytes=len(files[PARENT+key]))
    files[PARENT+'receipt.json']=json.dumps(receipt).encode()
    with pytest.raises(ValueError,match='control differs'):
        validate(files,target)


def test_parent_full_release_cannot_be_relabelled_success(native):
    files,target=native;files=dict(files)
    receipt=json.loads(files[PARENT+'receipt.json']);receipt['status']='COMPLETED'
    files[PARENT+'receipt.json']=json.dumps(receipt).encode()
    with pytest.raises(ValueError,match='scope'):
        validate(files,target)


def test_complete_native_decay_replays_without_accuracy_admission():
    from rya1354_frozen_summary import summarize
    base=ROOT/'data/audit/rya1354'
    report=summarize(base)
    assert report==json.loads((base/'source_curation_v2/frozen_lifetime_v12.json').read_text())
    assert len(report['decay']['channels'])==13
    assert report['decay']['model']['B']['model_one_photon_lifetime_ns']==pytest.approx(2.40239334941)
    assert report['decay']['model']['C']['model_one_photon_lifetime_ns']==pytest.approx(1.68034761772)
    assert report['accepted_lifetime'] is None and not report['fully_released_qualified']


def test_rate_replay_refuses_replaced_wavefunction(monkeypatch):
    import rya1354_frozen_summary as audit
    original=audit.load_parts
    def changed(path):
        files,manifest=original(path)
        if path.name=='frozen_complete_v12':
            files[audit.RUN+'to_ground/upper.w']=files[audit.RUN+'ground/state.w']
        return files,manifest
    monkeypatch.setattr(audit,'load_parts',changed)
    with pytest.raises(ValueError,match='transition input differs'):
        audit.summarize(ROOT/'data/audit/rya1354')
