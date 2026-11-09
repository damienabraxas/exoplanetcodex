"""Reference adapter input identity, resource and ownership refusals."""
import copy,json
from pathlib import Path
import pytest
from pipeline.reference_engines import prepare_run,validate_plan,fingerprint,inside,quote_nml
ROOT=Path(__file__).resolve().parents[1]


def fixture_plan(tmp_path):
    source=tmp_path/'source';source.mkdir();(source/'atmos.mod').write_text('test atmosphere receipt only\n');atom=source/'atom';atom.write_text(' O\n 8.70 15.999\n 41 256 22 0\n')
    p=json.loads((ROOT/'config/experiments/rya1371_openmulti3d_oxygen_1d.json').read_text());p['atom']=dict(resource='atom',receipt=fingerprint(atom));p['atmosphere']=dict(resource='atmos',filename='atmos.mod',receipt=fingerprint(source/'atmos.mod'))
    app=tmp_path/'app';(app/'input').mkdir(parents=True)
    return p,app,lambda key:atom if key=='atom'else source


def test_wrong_element_or_changed_atom_fails_before_allocating_output(tmp_path):
    p,app,resolve=fixture_plan(tmp_path);p['expected_element']='Fe';out=tmp_path/'run'
    with pytest.raises(ValueError,match='element'):prepare_run(p,out,app,resolver=resolve)
    assert not out.exists();p['expected_element']='O';resolve('atom').write_text(' O\n 8.71 15.999\n 41 256 22 0\n')
    with pytest.raises(ValueError,match='Changed'):prepare_run(p,out,app,resolver=resolve)
    assert not out.exists()


def test_preparation_is_fresh_and_preserves_source_and_explicit_input_identity(tmp_path):
    p,app,resolve=fixture_plan(tmp_path);source=resolve('atom').read_bytes();out=tmp_path/'run';r=prepare_run(p,out,app,resolver=resolve)
    assert r['element']=='O'and r['atom_header']['levels']==41
    assert resolve('atom').read_bytes()==source
    text=(out/'input.nml').read_text();assert 'include_absmet=F'in text and 'include_linelist=F'in text and str(resolve('atom'))in text
    with pytest.raises(FileExistsError):prepare_run(p,out,app,resolver=resolve)
    with pytest.raises(ValueError,match='outside'):prepare_run(p,app/'bad',app,resolver=resolve)


def test_unknown_dimensions_and_over_budget_runs_cannot_execute(tmp_path):
    p,_,_=fixture_plan(tmp_path)
    for key,value in [('mode','3D-NLTE'),('engine','unknown')]:
        bad=copy.deepcopy(p);bad[key]=value
        with pytest.raises(ValueError):validate_plan(bad)
    for key,value in [('threads',2),('seconds',601),('rtol',1e-3),('maxiter',1501)]:
        bad=copy.deepcopy(p);bad['controls'][key]=value
        with pytest.raises(ValueError):validate_plan(bad)
    for path in ['../escape','/absolute']:
        with pytest.raises(ValueError):inside(tmp_path,path)
    with pytest.raises(ValueError):quote_nml("unsafe'input")
