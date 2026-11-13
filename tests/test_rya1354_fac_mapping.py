"""Signed native/export FAC coefficient transport, not physical LS labels."""
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_fac_mapping import verify_transport


@pytest.fixture(scope='module')
def native():
    return load_parts(ROOT/'data/audit/rya1354/fac_coherent_complete_v18')[0]


@pytest.mark.parametrize('case',['coherent10','coherent100'])
def test_real_native_signed_coefficients_survive_basis_reordering(native,case):
    p='runs/fac_coherent_controls_v18/'+case+'/'
    r=verify_transport(native[p+'o.basis'].decode(),native[p+'o_export_000.c'].decode(),native[p+'o_export_000.cm'])
    assert r['signed_coefficients_checked']==43 and r['exported_CSFs']==15
    assert r['level_indices_preserved']==15


def test_native_coefficient_sign_substitution_is_refused(native):
    p='runs/fac_coherent_controls_v18/coherent100/'
    text=native[p+'o.basis'].decode()
    before,after=text.split('============Mixing Coefficients===================')
    lines=after.splitlines()
    i=next(i for i,line in enumerate(lines) if line.strip() and not line.startswith('#'))
    fields=lines[i].split();fields[9]=format(-float(fields[9]),'.8E');lines[i]=' '.join(fields)
    changed=before+'============Mixing Coefficients==================='+'\n'.join(lines)
    with pytest.raises(ValueError,match='coefficient/sign/order differs'):
        verify_transport(changed,native[p+'o_export_000.c'].decode(),native[p+'o_export_000.cm'])


def test_complete_native_angular_replay_maps_targets_without_rates():
    import json
    from rya1354_fac_recoupling_audit import summarize
    base=ROOT/'data/audit/rya1354'
    result=summarize(base)
    assert result==json.loads((base/'source_curation_v2/fac_recoupling_v19.json').read_text())
    for case in result['cases'].values():
        angular=case['recoupling']
        assert angular['targets']=={'3S':[6],'5S':[5]}
        assert angular['full_LS_coefficients_reconstructed']==43 and angular['unique_term_controls']==3
    assert not result['rates_admitted'] and not result['independent_angular_implementation']
    assert result['accepted_lifetime'] is result['uncertainty'] is None


@pytest.mark.parametrize('change',['LS_sign','level_ID','input'])
def test_recoupling_audit_refuses_sign_id_and_input_substitutions(monkeypatch,change):
    import rya1354_fac_recoupling_audit as audit
    original=audit.load_parts
    def changed(path):
        files,manifest=original(path)
        if path.name=='fac_recoupling_complete_v19':
            prefix=audit.RUN+'coherent100/'
            if change=='LS_sign':
                key=prefix+'state.lsj.j';value=files[key]
                files[key]=value.replace(b'-0.00564776',b' 0.00564776',1)
                assert files[key]!=value
            elif change=='level_ID':
                key=prefix+'state.lsj.lbl';value=files[key]
                files[key]=value.replace(b'  0    2',b'999    2',1)
                assert files[key]!=value
            else:
                files[prefix+'state.cm']=files[audit.RUN+'coherent10/state.cm']
        return files,manifest
    monkeypatch.setattr(audit,'load_parts',changed)
    with pytest.raises(ValueError,match='coefficient/sign/order reconstruction differs|level IDs incomplete|exact input/controls differ'):
        audit.summarize(ROOT/'data/audit/rya1354')
