"""Balanced space enumeration and independent state-count guard."""
from collections import Counter
from math import comb,prod
from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_space import configurations,expected_inventory,input_text,REFERENCES,SHELLS,inspect
from rya1354_native_parts import load_parts


def test_determinant_count_reproduces_actual_minimal_native_inventory():
    from rya1354_fac_bootstrap import level_inventory
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_coherent_complete_v18')
    levels=level_inventory(files['runs/fac_coherent_controls_v18/coherent100/o.lev'].decode())
    assert expected_inventory(REFERENCES)==Counter((r['parity'],r['two_J']) for r in levels)


@pytest.mark.parametrize('rank,count,levels',[(1,7,184),(2,15,788)])
def test_balanced_spaces_and_complete_magnetic_degeneracy(rank,count,levels):
    configs=configurations(rank);inventory=expected_inventory(configs)
    assert len(configs)==count and sum(inventory.values())==levels
    assert set(REFERENCES)<=set(configs)
    assert sum(n*(j+1) for (p,j),n in inventory.items())==sum(prod(comb(cap,q) for (_,_,cap),q in zip(SHELLS,c)) for c in configs)
    if rank==2:assert set(configurations(1))<set(configs)
    for mult,scale in ((.3,10),(.03,100)):
        text=input_text(rank,mult,scale)
        assert text.count('Config(')==count
        assert "OptimizeRadial(['ground'])" in text and 'TransitionTable' not in text
        assert 'SetBreit(0)' in text and 'SetSE(0)' in text and 'SetVP(0)' in text


def test_old_minimal_inventory_cannot_pass_as_expanded_space():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_coherent_complete_v18')
    with pytest.raises(ValueError,match='determinant count'):
        inspect(files['runs/fac_coherent_controls_v18/coherent100/o.lev'].decode(),1)


def test_actual_expanded_native_transport_all_four_cases():
    from rya1354_fac_space_audit import summarize
    r=summarize(ROOT/'data/audit/rya1354')
    assert len(r['cases'])==4 and r['accepted_lifetime'] is None
    for name,case in r['cases'].items():
        tr=case['transport'];rank=1 if '_r1_' in name else 2
        assert tr['levels']==(184 if rank==1 else 788)
        assert tr['declared_configurations_present']==len(configurations(rank))
        assert not tr['LS_mapping_admitted'] and tr['signed_printed_coefficients_verified']>100


@pytest.mark.parametrize('mutation',['sign','occupation','truncated'])
def test_expanded_transport_rejects_corruption(mutation):
    from rya1354_fac_space_audit import transport
    import re
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_space_complete_v23');p='runs/fac_space_v23/n3_r1_coherent10/'
    basis=files[p+'o.basis'].decode();c=files[p+'o_export_000.c'].decode();cm=files[p+'o_export_000.cm']
    if mutation=='sign':
        before,after=basis.split('Mixing Coefficients',1)
        lines=after.splitlines()
        for i,line in enumerate(lines):
            fields=line.split()
            if len(fields)==11:
                value=fields[9];new=value[1:] if value.startswith('-') else '-'+value
                lines[i]=line.replace(value,new,1);break
        after='\n'.join(lines)
        basis=before+'Mixing Coefficients'+after
    if mutation=='occupation':c=c.replace('2p ( 2)','2p ( 1)',1)
    if mutation=='truncated':cm=cm[:-1]
    with pytest.raises(ValueError):transport(basis,c,cm,1)
