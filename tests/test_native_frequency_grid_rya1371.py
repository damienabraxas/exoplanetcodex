import copy,json
from pathlib import Path
import numpy as np
import pytest
from pipeline.reference_engines import ROOT,fingerprint
from pipeline.native_reference_controls import validate_native_grid_plan
from pipeline.native_frequency_grid import load_native_frequency_grid


def fixture(tmp_path):
    qualified=ROOT/'config/experiments/rya1371_openmulti3d_oxygen_shared16.json'
    plan=json.loads(qualified.read_text());source=tmp_path/'source';(source/'native/save').mkdir(parents=True)
    wave=np.arange(39841,dtype=float)*2+100
    binary=source/'native/save/all_lam.bin';wave.astype('<f8').tofile(binary)
    report=source/'report.json';report.write_text(json.dumps(dict(computational_completed=True,plan=fingerprint(qualified),output_files=[dict(path='native/save/all_lam.bin',receipt=fingerprint(binary))])))
    directory=tmp_path/'grid';directory.mkdir();np.savetxt(directory/'frequency_grid_vac_A.txt',wave,fmt='%.17e')
    producer=directory/'report.json';producer.write_text(json.dumps(dict(computational_completed=True,all_original_nodes_exact=True,proposed_points=len(wave),grid=fingerprint(directory/'frequency_grid_vac_A.txt'),inputs={str(report):fingerprint(report),str(binary):fingerprint(binary)})))
    plan['native_frequency_grid_control']=dict(kind='baseline',resource='grid',producer=fingerprint(producer),grid=fingerprint(directory/'frequency_grid_vac_A.txt'),source_resource='source',source_report=fingerprint(report),source_wave=fingerprint(binary))
    return plan,lambda key:{'grid':directory,'source':source}[key]


def update_grid(plan,resolve,wave):
    directory=resolve('grid');grid=directory/'frequency_grid_vac_A.txt';np.savetxt(grid,wave,fmt='%.17e')
    p=directory/'report.json';d=json.loads(p.read_text());d.update(proposed_points=len(wave),grid=fingerprint(grid));p.write_text(json.dumps(d))
    plan['native_frequency_grid_control'].update(grid=fingerprint(grid),producer=fingerprint(p))


def test_exact_original_grid_and_valid_refinement(tmp_path):
    plan,resolve=fixture(tmp_path);path,audit=load_native_frequency_grid(plan,resolver=resolve)
    assert audit['points']==39841 and audit['all_original_nodes_exact']
    plan['native_frequency_grid_control']['kind']='continuum-refinement'
    wave=np.loadtxt(path);update_grid(plan,resolve,np.sort(np.r_[wave,101.]))
    assert load_native_frequency_grid(plan,resolver=resolve)[1]['points']==39842


@pytest.mark.parametrize('change',['atom','rtol','seconds','unknown','kind'])
def test_other_axes_and_incomplete_control_refuse(tmp_path,change):
    plan,resolve=fixture(tmp_path);qualified=json.loads((ROOT/'config/experiments/rya1371_openmulti3d_oxygen_shared16.json').read_text())
    if change=='atom':plan['abundance']=8.8
    elif change=='rtol':plan['controls']['rtol']=1e-5
    elif change=='seconds':plan['controls']['seconds']=601
    elif change=='unknown':plan['native_frequency_grid_control']['extra']=True
    else:plan['native_frequency_grid_control']['kind']='other'
    with pytest.raises(ValueError):validate_native_grid_plan(plan,qualified)


@pytest.mark.parametrize('change',['changed_bytes','removed','endpoint','close','baseline_enriched','missing_file'])
def test_grid_identity_and_retention_refusals(tmp_path,change):
    plan,resolve=fixture(tmp_path);path=resolve('grid')/'frequency_grid_vac_A.txt';wave=np.loadtxt(path)
    if change=='changed_bytes':path.write_text('100\n102\n')
    elif change=='missing_file':path.rename(path.with_suffix('.missing'))
    else:
        if change=='removed':wave=np.sort(np.r_[np.delete(wave,100),101.])
        elif change=='endpoint':wave=np.r_[wave,80000.]
        elif change=='close':wave=np.sort(np.r_[wave,100.0001]);plan['native_frequency_grid_control']['kind']='continuum-refinement'
        else:wave=np.sort(np.r_[wave,101.])
        update_grid(plan,resolve,wave)
    with pytest.raises((ValueError,FileNotFoundError)):load_native_frequency_grid(plan,resolver=resolve)


def test_wrong_source_producer_refuses_even_with_fresh_receipt(tmp_path):
    plan,resolve=fixture(tmp_path);p=resolve('grid')/'report.json';d=json.loads(p.read_text());d['inputs']={};p.write_text(json.dumps(d));plan['native_frequency_grid_control']['producer']=fingerprint(p)
    with pytest.raises(ValueError,match='another source'):load_native_frequency_grid(plan,resolver=resolve)


def second_fixture(tmp_path):
    plan,resolve=fixture(tmp_path);original=np.loadtxt(resolve('grid')/'frequency_grid_vac_A.txt')
    parent_wave=np.sort(np.r_[original,101.+2*np.arange(19162)])
    parent=tmp_path/'parent';(parent/'native/save').mkdir(parents=True)
    binary=parent/'native/save/all_lam.bin';parent_wave.astype('<f8').tofile(binary)
    report=parent/'report.json';report.write_text(json.dumps(dict(computational_completed=True,plan=fingerprint(ROOT/'config/experiments/rya1371_openmulti3d_oxygen_grid_refined.json'),output_files=[dict(path='native/save/all_lam.bin',receipt=fingerprint(binary))])))
    plan['native_frequency_grid_control']['kind']='continuum-refinement2'
    update_grid(plan,resolve,np.sort(np.r_[parent_wave,40101.]))
    p=resolve('grid')/'report.json';doc=json.loads(p.read_text());doc['predecessor']=dict(resource='experiments.openmulti3d_oxygen_grid_refined_run',report=fingerprint(report),wave=fingerprint(binary));doc['inputs'].update({str(report):fingerprint(report),str(binary):fingerprint(binary)});p.write_text(json.dumps(doc));plan['native_frequency_grid_control']['producer']=fingerprint(p)
    return plan,lambda key:parent if key=='experiments.openmulti3d_oxygen_grid_refined_run'else resolve(key)


def test_second_refinement_requires_every_parent_node_and_separate_scope(tmp_path):
    plan,resolve=second_fixture(tmp_path)
    assert load_native_frequency_grid(plan,resolver=resolve)[1]['points']==59004
    first=copy.deepcopy(plan);first['native_frequency_grid_control']['kind']='continuum-refinement'
    with pytest.raises(ValueError,match='scoped'):load_native_frequency_grid(first,resolver=resolve)
    wave=np.loadtxt(resolve('grid')/'frequency_grid_vac_A.txt');wave=np.sort(np.r_[wave[wave!=101.],40103.]);update_grid(plan,resolve,wave)
    with pytest.raises(ValueError,match='predecessor nodes'):load_native_frequency_grid(plan,resolver=resolve)


def test_second_refinement_refuses_a_changed_parent_receipt(tmp_path):
    plan,resolve=second_fixture(tmp_path);p=resolve('experiments.openmulti3d_oxygen_grid_refined_run')/'report.json';p.write_text(p.read_text()+'\n')
    with pytest.raises(ValueError,match='predecessor'):load_native_frequency_grid(plan,resolver=resolve)
