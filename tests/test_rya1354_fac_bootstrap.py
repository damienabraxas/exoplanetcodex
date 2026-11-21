"""Pinned upstream FAC inventory parsing; no oxygen/accuracy admission."""
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_fac_bootstrap import level_inventory,PIN


@pytest.fixture(scope='module')
def native_reference():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_intake_v17')
    return files['evidence/fac_bootstrap_v17/source_documents/demo/structure/ref/ne.lev'].decode()


def test_upstream_reference_inventory_is_explicit(native_reference):
    rows=level_inventory(native_reference)
    assert len(rows)==37 and rows[0]['index']==0
    assert rows[0]['parity']==rows[0]['two_J']==0
    assert rows[0]['labels'][1]=='2p6'
    assert PIN=='9ab6defb6af62e1bc8251a35ff995809abef717c'


def test_inventory_rejects_nonfinite_and_duplicate_native_levels(native_reference):
    with pytest.raises(ValueError,match='quantum numbers'):
        level_inventory(native_reference.replace('0.00000000E+00','NaN',1))
    line=next(line for line in native_reference.splitlines() if line.split()[:2]==['0','-1'])
    with pytest.raises(ValueError,match='duplicate'):
        level_inventory(native_reference+'\n'+line+'\n')


def test_completed_native_fac_inventory_does_not_admit_rates_or_convergence():
    import json
    from rya1354_fac_summary import summarize
    base=ROOT/'data/audit/rya1354'
    result=summarize(base)
    assert result==json.loads((base/'source_curation_v2/fac_feasibility_v17.json').read_text())
    assert result['upstream_example_levels']==37 and len(result['oxygen_levels'])==15
    assert not result['numerical_convergence_admitted'] and not result['physical_LS_terms_admitted']
    assert not result['adoption_allowed'] and result['accepted_lifetime'] is None


@pytest.mark.parametrize('change',['physics','electron_count'])
def test_fac_feasibility_refuses_changed_model_or_species(monkeypatch,change):
    import rya1354_fac_summary as audit
    original=audit.load_parts
    def changed(path):
        files,manifest=original(path)
        if path.name=='fac_oxygen_complete_v17':
            if change=='physics':
                key=audit.OXYGEN+'input.sf';files[key]=files[key].replace(b'SetBreit(0)',b'SetBreit(5)')
            else:
                key=audit.OXYGEN+'o.lev';files[key]=files[key].replace(b'NELE\t= 8',b'NELE\t= 7')
        return files,manifest
    monkeypatch.setattr(audit,'load_parts',changed)
    with pytest.raises(ValueError,match='physics/configuration input differs|electron/parity/J inventory differs'):
        audit.summarize(ROOT/'data/audit/rya1354')
