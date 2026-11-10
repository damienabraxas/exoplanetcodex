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
    for key,value in [('mode','3D-NLTE'),('engine','unknown'),('profile_variant','unreviewed')]:
        bad=copy.deepcopy(p);bad[key]=value
        with pytest.raises(ValueError):validate_plan(bad)
    for key,value in [('threads',2),('seconds',601),('rtol',1e-3),('maxiter',1501),('ntheta',7),('molecular_equilibrium','disabled-control'),('output_solver','unsupported'),('disable_stark_control',1)]:
        bad=copy.deepcopy(p);bad['controls'][key]=value
        with pytest.raises(ValueError):validate_plan(bad)
    for path in ['../escape','/absolute']:
        with pytest.raises(ValueError):inside(tmp_path,path)
    with pytest.raises(ValueError):quote_nml("unsafe'input")


def test_explicit_output_grid_is_canonical_vacuum_and_gap_preserving(tmp_path):
    import numpy as np
    from pipeline.wavelength_util import vac_to_air,air_to_vac
    p,app,resolve=fixture_plan(tmp_path)
    p['output_grid']=dict(frame='air_A',step_A=.005,windows=[[6300.,6300.02],[6363.,6363.02]])
    out=tmp_path/'dense';prepare_run(p,out,app,resolver=resolve)
    vac=np.loadtxt(out/'output_wavelengths_vac_A.txt');air=vac_to_air(vac)
    np.testing.assert_allclose(vac[[0,-1]],air_to_vac(np.array([6300.,6363.02])),rtol=0,atol=1e-12)
    assert np.max(np.diff(air))>60
    assert "in_air=F" in (out/'input.nml').read_text()
    p['output_grid']['windows']=[[6363.,6364.],[6300.,6301.]]
    with pytest.raises(ValueError):validate_plan(p)


def test_engine_pin_and_native_build_mutations_are_refused(tmp_path,monkeypatch):
    import yaml
    from pipeline.reference_engines import load_engine
    checkout=tmp_path/'checkout';app=checkout/'app';app.mkdir(parents=True)
    exe=app/'dispatch.x';exe.write_bytes(b'native executable')
    lib=app/'HOP/standalone/build/libhop.so';lib.parent.mkdir(parents=True);lib.write_bytes(b'native library')
    build=tmp_path/'build';build.mkdir();receipt=build/'receipt.json'
    receipt.write_text(json.dumps(dict(source_commit='pin',tracked_changes='',executable=fingerprint(exe),libraries=dict(hop=fingerprint(lib)))))
    registry=tmp_path/'registry.yaml';registry.write_text(yaml.safe_dump(dict(schema='codex.reference_engines/1',engines=dict(openmulti3d=dict(status='experimental-reference',production_admission='HOLD',checkout_resource='checkout',build_resource='build',application_rel='app',executable_rel='dispatch.x',source_commit='pin')))))
    state=dict(commit='pin',dirty='')
    monkeypatch.setattr('pipeline.reference_engines.subprocess.check_output',lambda argv,**kw:state['commit']if argv[1]=='rev-parse'else state['dirty'])
    resolve=lambda key:checkout if key=='checkout'else build
    load_engine('openmulti3d',registry=registry,resolver=resolve)
    state['commit']='other'
    with pytest.raises(ValueError,match='source'):load_engine('openmulti3d',registry=registry,resolver=resolve)
    state['commit']='pin';state['dirty']=' M tracked-source'
    with pytest.raises(ValueError,match='source'):load_engine('openmulti3d',registry=registry,resolver=resolve)
    state['dirty']='';exe.write_bytes(b'mutated')
    with pytest.raises(ValueError,match='executable'):load_engine('openmulti3d',registry=registry,resolver=resolve)
    exe.write_bytes(b'native executable');lib.write_bytes(b'mutated')
    with pytest.raises(ValueError,match='HOP'):load_engine('openmulti3d',registry=registry,resolver=resolve)


def test_profile_control_is_explicit_and_original_dispatch_stays_identical(monkeypatch):
    import pipeline.reference_engines as module
    original=(object(),object(),object(),{'original':'sentinel'})
    monkeypatch.setattr(module,'load_engine',lambda *args,**kwargs:original)
    assert module.load_plan_engine({'engine':'openmulti3d'})==original
    with pytest.raises(ValueError,match='profile variant'):
        module.load_plan_engine({'engine':'openmulti3d','profile_variant':'unreviewed'})


def test_derived_profile_refuses_any_unexpected_source_diff_including_staged_changes(tmp_path,monkeypatch):
    import pipeline.reference_engines as module
    pin='32085d26c46e9919a986397b5514a9728abe2074'
    spec=dict(source_commit=pin,application_rel='application',executable_rel='dispatch.x')
    monkeypatch.setattr(module,'load_engine',lambda *args,**kwargs:(spec,None,None,{'original':'sentinel'}))
    source=tmp_path/'source';source.mkdir();build=tmp_path/'build';build.mkdir();(build/'receipt.json').write_text('{}')
    calls=[]
    def output(args,**kwargs):
        calls.append(args)
        return pin+'\n'if args[1]=='rev-parse'else b'unexpected source diff'
    monkeypatch.setattr(module.subprocess,'check_output',output)
    resolve=lambda key:source if key.endswith('profile_reference')else build
    with pytest.raises(ValueError,match='derivative source'):
        module.load_plan_engine(dict(engine='openmulti3d',profile_variant='accurate-kernel-control'),resolver=resolve)
    assert ['git','diff','--binary','HEAD']in calls
