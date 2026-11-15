"""Complete valence closure count, exact input addition and transport controls."""
from itertools import product
from math import comb
from pathlib import Path
import ast
import re
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_full_space import configurations,input_text
from rya1354_fac_space import configurations as truncated,SHELLS,expected_inventory
from rya1354_fac_no_correction import input_text as sd_input
from rya1354_native_parts import load_parts
from rya1354_fac_full_space_audit import transport
from rya1354_fac_space_audit import transport as sd_transport


def test_full_space_is_exhaustive_and_counts_all_determinants():
    q=configurations();assert len(q)==31 and set(truncated(2))<set(q)
    independent={v for v in product(range(7),range(3),range(7),range(11)) if sum(v)==4}
    assert set(q)==independent
    counts=expected_inventory(q)
    assert sum(counts.values())==1858 and max(counts.values())==253
    assert sum(n*(j+1) for (p,j),n in counts.items())==comb(24,4)


def test_inputs_only_add_sixteen_configuration_groups():
    for policy in ('ground','equal'):
        for mult,scale in ((.3,10),(.03,100)):
            old=sd_input(policy,mult,scale);new=input_text(policy,mult,scale)
            assert new.count('Config(')==31 and new.count('ConfigDiagnostics(')==2
            lines=[l for l in new.splitlines() if not l.startswith("Config('closure")]
            old_structure=next(l for l in old.splitlines() if l.startswith('Structure('))
            lines=[old_structure if l.startswith('Structure(') else l for l in lines]
            assert '\n'.join(lines)+'\n'==old
            structure=next(l for l in new.splitlines() if l.startswith('Structure('))
            groups=ast.literal_eval(re.search(r'(\[.*\])',structure)[0]);assert len(groups)==len(set(groups))==31


def test_generic_transport_preserves_previous_actual_sd_control():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_space_complete_v23');p='runs/fac_space_v23/n3_r2_coherent100/'
    args=(files[p+'o.basis'].decode(),files[p+'o_export_000.c'].decode(),files[p+'o_export_000.cm'])
    assert transport(*args,declared=truncated(2))==sd_transport(*args,2)


def test_nk_d3_d4_seniority_labels_are_distinct_and_source_consistent():
    from rya1354_fac_full_ls_audit import packed,nk_d_term
    assert nk_d_term('2D1')=='2D1' and nk_d_term('2D3')=='2D2'
    assert nk_d_term('1S0')=='1S1' and nk_d_term('1S4')=='1S2'
    occ=[(1,'s',2),(2,'s',2),(3,'d',4)]
    a=packed(occ,['1S0','1S0','3P2','1S','3P'])
    b=packed(occ,['1S0','1S0','3P4','1S','3P'])
    assert a=='1s(2).2s(2).3d(4)3P1_3P' and b=='1s(2).2s(2).3d(4)3P2_3P'


def test_full_ls_formatter_preserves_actual_previous_sd_names():
    from rya1354_fac_full_ls_audit import ls_blocks as full_ls
    from rya1354_fac_expanded_ls_audit import ls_blocks as sd_ls
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_expanded_angular_complete_v24')
    text=files['runs/fac_expanded_recoupling_v24/n3_r2_coherent100/state.lsj.c'].decode()
    assert [[r['name'] for r in b] for b in full_ls(text)]==[[r['name'] for r in b] for b in sd_ls(text,2)]


def test_overflowed_display_ids_require_explicit_stdout_and_binary_agreement():
    import pytest
    from rya1354_fac_full_ls_audit import label_rows
    expected=[(i,0,'+') for i in range(1001)]
    text='\n'.join(f'{i if i<=999 else "***"} 0 + -74.000000000 100.000%\n 1.00000000 1.00000000 1s(2).2s(2).3d(4)1S1_1S' for i in range(1001))
    rows=label_rows(text,1001,expected,expected)
    assert rows[1000]['display_overflow'] and not rows[999]['display_overflow']
    bad=list(expected);bad[-1]=(999,0,'+')
    with pytest.raises(ValueError,match='stdout'):label_rows(text,1001,expected,bad)
    with pytest.raises(ValueError,match='Unqualified'):label_rows(text,1001)


def test_actual_full_space_all_four_cases_preserve_potentials_and_corrections():
    from rya1354_fac_full_space_audit import summarize
    r=summarize(ROOT/'data/audit/rya1354')
    assert len(r['cases'])==4 and r['accepted_lifetime'] is None
    for case in r['cases'].values():
        assert case['paired_saved_potential_byte_identical']
        assert case['transport']['levels']==1858 and case['transport']['declared_configurations_present']==31
        for row in case['configuration_state'].values():
            assert row['groups']==31 and row['configurations']==189 and row['max_abs_delta']==row['max_abs_shift']==0
