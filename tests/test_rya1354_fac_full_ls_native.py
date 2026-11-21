"""Metadata validation of actual full LS outputs; transforms remain Sirius-only."""
from pathlib import Path
import json
import struct
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_fac_recoupling_audit import records
from rya1354_fac_full_ls_audit import ls_blocks,label_rows,stdout_ids
from rya1354_fac_full_rates import rate_input,AUDIT_SHA


def test_all_full_native_seniority_names_and_overflowed_ID_authorities():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_full_angular_complete_v29');report=json.loads(files['evidence/fac_full_replay_v29b/report.json'])
    import hashlib
    assert hashlib.sha256(files['evidence/fac_full_replay_v29b/report.json']).hexdigest()==AUDIT_SHA
    for name,r in report['cases'].items():
        p='runs/fac_full_recoupling_v29/'+name+'/';cm=records(files[p+'state.cm']);_,_,_,nlev,_,nb=struct.unpack('<6i',cm[1]);expected=[]
        for b in range(nb):
            _,n,nv,jp,parity=struct.unpack('<5i',cm[2+4*b]);expected.extend((i,jp-1,'+' if parity==1 else '-') for i in struct.unpack('<'+'i'*nv,cm[3+4*b]))
        lab=label_rows(files[p+'state.lsj.lbl'].decode(),nlev,expected,stdout_ids(files[p+'jj2lsj.stdout'].decode()))
        names={row['name'] for block in ls_blocks(files[p+'state.lsj.c'].decode()) for row in block}
        assert nlev==1858 and r['label_ID_overflows']==858 and r['full_LS_coefficients_reconstructed']==329198
        for i,row in lab.items():
            m=r['levels'][str(i)]
            assert row['display_overflow']==m['label_display_id_overflow']
            assert row['components'][0][2]==m['dominant'] and len(row['components'])==m['native_component_count']
            assert all(name in names for _,_,name in row['components'])


def test_full_native_positive_and_four_negative_replay_controls():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_full_angular_complete_v29');r=json.loads(files['evidence/fac_full_canaries_v29/report.json'])
    assert r['status']=='POSITIVE_AND_FOUR_NEGATIVE_CONTROLS_PASS'
    assert set(r['results'])=={'d4_seniority','vector_sign','matrix_sign','explicit_stdout_ID'}
    assert not files['evidence/fac_full_canaries_v29/stderr.log']
    assert b'Incomplete LS IDs' in files['evidence/fac_full_replay_v29/replay.stderr']


def test_actual_full_mapping_rate_builder_reenumerates_both_gauges():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_full_angular_complete_v29');r=json.loads(files['evidence/fac_full_replay_v29b/report.json'])
    for name,m in r['cases'].items():
        policy=name.split('_')[0];mult,scale=(.03,100) if name.endswith('100') else (.3,10)
        for gauge in (1,2):
            text,requests,lower=rate_input(policy,mult,scale,gauge,m)
            assert text.count('Config(')==31 and len(requests)==13 and lower==list(range(6))
            assert len({(v['upper'],v['lower'],v['multipole']) for v in requests})==len(requests)
            assert f'SetTransitionOptions({gauge},0)' in text
