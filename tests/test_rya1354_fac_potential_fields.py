from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_fac_potential_fields import layout,compare
from rya1354_fac_radial_preflight import review


def actual():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_radial_probe_failed_v31');root='runs/fac_radial_probe_v31/'
    spec=layout(files[root+'source_documents/faclib/radial.c'].decode(),files[root+'source_documents/faclib/orbital.c'].decode(),files[root+'source_documents/faclib/consts.h'].decode())
    p=root+'ground_coherent10/';data={k.removeprefix(p):v for k,v in files.items() if k.startswith(p)}
    return spec,data


def test_source_mapped_changes_are_only_named_energy_kappa_work_arrays():
    spec,d=actual();assert spec['header']==348 and spec['arrays']==69
    assert spec['work_arrays']=={'W':22,'dW':23,'dW2':24}
    assert compare(d['source.pot'],d['restored.pot'],spec)['serialized_byte_identical']
    for fn in ('old_solved.pot','extended.pot'):
        r=compare(d['source.pot'],d[fn],spec)
        assert r['protected_header_grid_potential_byte_identical'] and not r['serialized_byte_identical']
        assert set(r['derived_energy_kappa_work_arrays_changed'])=={'W','dW','dW2'}


@pytest.mark.parametrize('offset',[0,348+3*1500*8,348+19*1500*8,348+60*1500*8])
def test_any_header_grid_central_or_QED_array_change_is_rejected(offset):
    spec,d=actual();changed=bytearray(d['old_solved.pot']);changed[offset]^=1
    with pytest.raises(ValueError,match='Protected'):compare(d['source.pot'],bytes(changed),spec)


def test_failed_v31_remains_failed_but_physical_orbital_diagnostics_replay():
    spec,d=actual();r=review(d,10,spec)
    import json
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_radial_probe_failed_v31')
    assert json.loads(files['runs/fac_radial_probe_v31/receipt.json'])['status']=='FAILED_RETAINED'
    assert len(r['old_orbitals'])==9 and len(r['new_orbitals'])==3
    assert all(v['energy_eV']<0 for v in r['new_orbitals'].values())
    assert sum(v['after']['pairs'] for v in r['overlaps'].values())==36


def test_all_four_fresh_preflights_replay_with_complete_overlap_coverage():
    from rya1354_fac_radial_preflight_audit import summarize
    report=summarize(ROOT/'data/audit/rya1354')
    assert report['adoption_allowed'] is False and report['accepted_lifetime'] is None
    assert len(report['cases'])==4
    for case in report['cases'].values():
        assert len(case['old_orbitals'])==9 and len(case['new_orbitals'])==3
        assert all(v['energy_eV']<0 for v in case['new_orbitals'].values())
        assert sum(v['before']['pairs'] for v in case['overlaps'].values())==19
        assert sum(v['after']['pairs'] for v in case['overlaps'].values())==36
        assert all(v['protected_header_grid_potential_byte_identical'] for v in case['potential_comparisons'].values())
