import numpy as np
import pytest
from scripts.probe_bride_cn_native_binding import reporting_sources,parse_binding,summarize,EQ_ANCHOR,ATM_ANCHOR


def test_reporting_insertion_is_unique_and_fixed_form():
    eq='        implicit none\n'+EQ_ANCHOR+'end\n';bsyn=ATM_ANCHOR+'end\n';a,b,block=reporting_sources(eq,bsyn)
    assert 'integer bride_hcount'in a and 'BRIDE_BIND_ATMOS'in b
    assert all(len(line)<=72 for line in block.splitlines()if not line.startswith('*'))
    with pytest.raises(ValueError):reporting_sources(eq+EQ_ANCHOR,bsyn)


def test_hydrogen_stoichiometry_and_incomplete_inventory_refusal():
    rho=(20.+.001)*float(np.float32(1.2123e-8))/5000
    text=f'BRIDE_BIND_EQ 1 1 1 5000 1 16 16 {rho} 1 20 1 .001\nBRIDE_BIND_ATMOS 1 5000 1 16 {rho} 100000 1.38066e-16 1.6735e-24 {rho} 1\nBRIDE_BIND_COMP 1 1 1 1\n'
    text+=''.join(f'BRIDE_BIND_SPEC 1 {i} {1 if i<5 else 2} {i} {1 if i<5 else 2}\n'for i in range(1,6))
    eq,atm,table,_=parse_binding(text,ndep=1);r=summarize(eq,atm,table)[0]
    assert r['H_nuclei_sum_dyn_cm2']==20
    assert r['H_in_molecules_fraction']==.5 and r['arithmetic_passed']
    for bad in [text+'BRIDE_BIND_SPEC 1 1 1 1 1\n',text.replace('BRIDE_BIND_SPEC 1 5 2 5 2\n','')]:
        with pytest.raises(ValueError):parse_binding(bad,ndep=1)
