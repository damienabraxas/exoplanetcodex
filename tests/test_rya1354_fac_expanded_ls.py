"""Metadata/native record replay only on Mac; numerical transforms Sirius-only."""
import json
from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_fac_expanded_ls_audit import ls_blocks,label_rows,review
from rya1354_fac_recoupling_audit import ls_basis
from rya1354_fac_expanded_rates import rate_input


def evidence():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_expanded_angular_complete_v24')
    report=json.loads(files['evidence/fac_expanded_replay_v24/report.json'])
    return files,report


def test_new_label_formatter_reproduces_actual_minimal_control():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_recoupling_complete_v19')
    text=files['runs/fac_recoupling_v19/coherent10/state.lsj.c'].decode()
    assert [[r['name'] for r in b] for b in ls_blocks(text,1)]==ls_basis(text)


def test_all_expanded_native_packed_names_and_original_ids():
    files,report=evidence()
    for name,r in report['cases'].items():
        rank=1 if '_r1_' in name else 2;p='runs/fac_expanded_recoupling_v24/'+name+'/'
        blocks=ls_blocks(files[p+'state.lsj.c'].decode(),rank)
        names={v['name'] for b in blocks for v in b};labels=label_rows(files[p+'state.lsj.lbl'].decode(),184 if rank==1 else 788)
        assert r['targets']=={'3S':[6],'5S':[5]}
        assert r['full_LS_coefficients_reconstructed']==(3896 if rank==1 else 63316)
        for i,row in labels.items():
            assert row['components'][0][2]==r['levels'][str(i)]['dominant']
            assert all(component[2] in names for component in row['components'])


def test_expanded_rate_requests_derived_from_all_states_and_parity():
    _,report=evidence()
    for name,mapping in report['cases'].items():
        rank=1 if '_r1_' in name else 2;mult,scale=(.03,100) if name.endswith('100') else (.3,10)
        for gauge in (1,2):
            text,requests,lower=rate_input(rank,mult,scale,gauge,mapping)
            assert lower==list(range(6)) and len(requests)==13
            assert {r['operator'] for r in requests}=={'E1','M2','E3','M1','E2','M3'}
            assert f'SetTransitionOptions({gauge},0)' in text and text.count('PrintTable(\'rates.tr.b\'')==1
            assert len({(r['upper'],r['lower'],r['multipole']) for r in requests})==13
        altered=json.loads(json.dumps(mapping));altered['levels']['787']=dict(two_J=8,parity='+',energy_au=mapping['levels']['6']['energy_au']-.001)
        _,requests,lower=rate_input(rank,mult,scale,2,altered)
        assert 787 in lower and {r['operator'] for r in requests if r['lower']==787}=={'E3','M4','E5'}


def test_full_angular_replay_refuses_mac_execution():
    import platform
    if platform.node().lower().split('.')[0]=='sirius':pytest.skip('Mac host guard')
    with pytest.raises(RuntimeError,match='Sirius'):review({}, {},1)
