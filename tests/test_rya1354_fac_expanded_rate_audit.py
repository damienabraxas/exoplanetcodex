"""Expanded native raw/formatted rates and guarded full-LS canary receipts."""
from pathlib import Path
import json
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_fac_expanded_rate_audit import summarize,review
from rya1354_fac_expanded_rates import rate_input
from rya1354_fac_rate_audit import constants


def test_complete_eight_case_rate_replay_and_no_adoption():
    report=summarize(ROOT/'data/audit/rya1354')
    assert len(report['cases'])==8 and report['accepted_lifetime'] is None
    for row in report['cases'].values():
        assert row['channels']==13 and row['printed_rows_verified']==12
        assert len(row['writer_suppressed'])==1 and row['writer_suppressed'][0]['multipole']==3
        assert all(v>0 for v in row['E3_A_s1'].values())


@pytest.mark.parametrize('mutation',['missing','status','gauge','energy','table','E3_row'])
def test_corruption_and_refusal_never_become_complete_total(mutation):
    base=ROOT/'data/audit/rya1354';f,_=load_parts(base/'fac_expanded_rates_complete_v25');a,_=load_parts(base/'fac_expanded_angular_complete_v24');p,_=load_parts(base/'fac_rates_complete_v21')
    name='n3_r1_coherent10';prefix='runs/fac_expanded_rates_v25/'+name+'_B/'
    mapping=json.loads(a['evidence/fac_expanded_replay_v24/report.json'])['cases'][name]
    _,requests,_=rate_input(1,.3,10,2,mapping);raw=f[prefix+'raw.tsv'].decode();table=f[prefix+'rates.tr'].decode()
    c=constants(p['runs/fac_model_rates_v21/source_documents/faclib/consts.h'].decode())
    lines=raw.splitlines();i=next(i for i,line in enumerate(lines) if line.strip() and not line.startswith('#'));fields=lines[i].split()
    if mutation=='missing':lines.pop(i)
    elif mutation=='status':fields[8]='-1';lines[i]=' '.join(fields)
    elif mutation=='gauge':fields[6]='1';lines[i]=' '.join(fields)
    elif mutation=='energy':fields[11]=str(float(fields[11])+.01);lines[i]=' '.join(fields)
    elif mutation=='table':table=table.replace('MODE\t= 0','MODE\t= 1',1)
    elif mutation=='E3_row':
        before,rest=table.split('MULTIP\t= -3',1);r=rest.splitlines();j=next(j for j,line in enumerate(r) if len(line.split())==8);r.pop(j);table=before+'MULTIP\t= -3'+'\n'.join(r)
    with pytest.raises(ValueError):review('\n'.join(lines),table,requests,2,mapping,c)


def test_sirius_full_replay_negative_cases_rejected():
    f,_=load_parts(ROOT/'data/audit/rya1354/fac_expanded_rates_complete_v25')
    r=json.loads(f['evidence/fac_expanded_canaries_v24/report.json'])
    assert r['status']=='ALL_FOUR_MUTATIONS_REJECTED'
    assert set(r['results'])=={'vector_sign','matrix_sign','original_ID','occupation'}
    assert not f['evidence/fac_expanded_canaries_v24/stderr.log']
