"""Source provenance and declared n5 radial scope, using retained native bytes."""
import json
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_frozen_extend import qualify_source, controls, PARENT, TARGET
from rya1354_frozen_layer import STATUS


@pytest.fixture
def sources(tmp_path):
    base=ROOT/'data/audit/rya1354'
    source=tmp_path/'source';parent=tmp_path/'parent';core=tmp_path/'core'
    for bundle,prefix,dest in [('frozen_complete_v12','runs/lifetime_frozen_prior_n4_v1/',source),
                               ('from_n3_terminal_v12',PARENT,parent),('lifetime_aux_v4b',TARGET,core)]:
        files,_=load_parts(base/bundle)
        for name,value in files.items():
            if name.startswith(prefix):
                p=dest/name.removeprefix(prefix);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(value)
    return source,parent,core


def test_qualified_n4_model_is_valid_extension_source(sources):
    result=qualify_source(*sources)
    assert result['source_audit']['parent_status']=='FAILED_RETAINED'
    assert not result['source_audit']['fully_released_qualified']
    for state in ('ground','quintet','triplet'):
        for verify in (False,True):
            text=controls(state,verify)
            assert '\n5*\n' in text and '2*,3*,4*,5*' not in text
            assert '\n439\ny\n5.9604644775390625e-08\n' in text
            assert '\ny\n24\nn\n1\n' in text
            assert '\n\n2s,2p-,2p' in controls(state,True,True)


def test_failed_or_changed_n4_source_is_refused(sources):
    source,parent,core=sources
    receipt=json.loads((source/'receipt.json').read_text())
    receipt['status']='FAILED_RETAINED'
    (source/'receipt.json').write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match='Unqualified'):
        qualify_source(source,parent,core)
    receipt['status']=STATUS;(source/'receipt.json').write_text(json.dumps(receipt))
    p=source/'triplet/state.w';p.write_bytes(p.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='artifact differs'):
        qualify_source(source,parent,core)


def test_native_layer_change_identity_and_corruption_refusal():
    from rya1354_layer_change import report,decompose
    from rya1354_lifetime_summary import decay
    base=ROOT/'data/audit/rya1354'
    result=report(base)
    assert result==json.loads((base/'source_curation_v2/layer_change_v13.json').read_text())
    assert all(.956<r['native_energy_cubed_ratio']<.957 for r in result['components'])
    files,_=load_parts(base/'lifetime_relax_v4')
    model=decay(files,'runs/lifetime_relax_n3_v3/')
    changed=json.loads(json.dumps(model));changed['channels'][0]['native_S_B_au']*=2
    with pytest.raises(ValueError,match='identity differs'):
        decompose(model,changed)


def test_collector_archives_frozen_terminal_without_admitting_science(tmp_path):
    from rya1354_collect_pending import observe
    run=tmp_path/'runs/lifetime_frozen_prior_n5_v1';run.mkdir(parents=True)
    (run/'receipt.json').write_text(json.dumps(dict(status=STATUS,files={})))
    result=observe(tmp_path,tmp_path/'evidence', ['lifetime_frozen_prior_n5_v1'])
    record=result['campaigns']['lifetime_frozen_prior_n5_v1']
    assert record['collection']=='ARCHIVED_TERMINAL_OUTPUTS'
    assert not record['adoption_allowed'] and record['scientific_admission']=='HOLD'
