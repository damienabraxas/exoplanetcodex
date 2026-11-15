import time
from types import SimpleNamespace
import numpy as np
import pytest
from pipeline.native_continuum_control import NativeContinuumControl


def fixture(tmp_path,monkeypatch):
    value=2.**-60;tables={(3,40):np.array([[943.253,value],[1064.179,value]])}
    exe=tmp_path/'oracle';exe.write_text('fixture executable identity only\n');calls=[]
    def run(args,*,input,**kwargs):
        records=input.splitlines();n,m=map(int,records[0].split());calls.append(records[n+1:]);return SimpleNamespace(returncode=0,stdout=('%.17e\n'%value)*m,stderr='')
    monkeypatch.setenv('OMP_NUM_THREADS','1');monkeypatch.setenv('OPENBLAS_NUM_THREADS','1')
    monkeypatch.setattr('pipeline.native_continuum_control.subprocess.run',run)
    return NativeContinuumControl(tables,exe,{},tmp_path/'queries',time.perf_counter()+60),calls,value


def test_canonical_source_knots_units_support_and_owned_cache(tmp_path,monkeypatch):
    control,calls,value=fixture(tmp_path,monkeypatch);edge=control.tables[3,40][:,0]*.1
    wave=np.array([edge[0]-.01,edge[0],edge.mean(),edge[1],edge[1]+.01])
    result=control.evaluate((3,40),wave);np.testing.assert_array_equal(result,np.array([0,value,value,value,0])*1e-4)
    assert float(calls[0][0])==943.253 and float(calls[0][-1])==1064.179
    result[1]=-1;assert control.evaluate((3,40),wave)[1]>0 and len(calls)==1
    control.verify_inputs()


def test_single_point_duplicate_query_and_empty_support(tmp_path,monkeypatch):
    control,calls,value=fixture(tmp_path,monkeypatch)
    assert control.evaluate((3,40),np.array([100.]))[0]==value*1e-4
    assert len(calls[0])==2 and calls[0][0]==calls[0][1]
    np.testing.assert_array_equal(control.evaluate((3,40),np.array([80.,90.])),[0.,0.]);assert len(calls)==1


@pytest.mark.parametrize('wave',[np.array([np.nan]),np.array([-1.]),np.ones(20001),np.ones((2,2))])
def test_invalid_queries_refuse_without_fallback(tmp_path,monkeypatch,wave):
    control,calls,_=fixture(tmp_path,monkeypatch)
    with pytest.raises(ValueError):control.evaluate((3,40),wave)
    assert not calls


@pytest.mark.parametrize('mode',['negative','nan','shape','process'])
def test_oracle_errors_are_preserved_and_refuse(tmp_path,monkeypatch,mode):
    control,_,_=fixture(tmp_path,monkeypatch)
    def bad(*args,**kwargs):return SimpleNamespace(returncode=1 if mode=='process'else 0,stdout='-1\n-1\n'if mode=='negative'else'nan\nnan\n'if mode=='nan'else'1\n',stderr='preserved fixture error')
    monkeypatch.setattr('pipeline.native_continuum_control.subprocess.run',bad)
    with pytest.raises((ValueError,RuntimeError)):control.evaluate((3,40),np.array([99.,100.]))
    assert list((tmp_path/'queries').glob('*.err')) and not control.cache


def test_changed_runtime_evidence_refuses(tmp_path,monkeypatch):
    control,_,_=fixture(tmp_path,monkeypatch);control.evaluate((3,40),np.array([99.,100.]))
    next((tmp_path/'queries').glob('*.out')).write_text('altered')
    with pytest.raises(ValueError,match='evidence changed'):control.verify_inputs()
