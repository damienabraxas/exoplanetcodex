"""Actual Git-pinned continuum policy controls; no duplicate estimator."""
import json
from pathlib import Path
import numpy as np
from scripts.analyze_bride_orch_continuum_parity import load_git_module
ROOT=Path(__file__).resolve().parents[1]


def test_corroborated_sag_is_applied_but_local_dip_and_missing_neighbours_are_refused():
    p=json.loads((ROOT/'config/experiments/rya1281_orch_continuum_parity.json').read_text())
    new,_=load_git_module(p['orch_commits']['current'],'test_bride_current_continuum')
    old,_=load_git_module(p['orch_commits']['previous'],'test_bride_previous_continuum')
    w=np.linspace(-60,60,12001);ref=np.ones_like(w);cut=abs(w)<=5
    for smooth,expected in [(True,True),(False,False)]:
        f=np.full_like(w,.9)if smooth else np.where(cut,.9,1.)
        before=f.copy();callback=new.neighbourhood_check(w,f,0.,w,ref,exclude_half_width_A=.5)
        previous,pc=old.fit_model_guided(w[cut],f[cut],0.,w,ref,exclude_half_width_A=.5)
        current,cc=new.fit_model_guided(w[cut],f[cut],0.,w,ref,exclude_half_width_A=.5,corroborate=callback)
        assert previous.applied and current.applied==expected
        np.testing.assert_array_equal(f,before)
        if expected:np.testing.assert_array_equal(pc,cc)
        else:assert cc is None and 'CONTINUUM_LOCAL'in current.reason
    r,c=new.fit_model_guided(w[cut],np.full(cut.sum(),.9),0.,w,ref,exclude_half_width_A=.5)
    assert not r.applied and c is None and 'no neighbourhood measured'in r.reason


def test_sub_three_percent_fit_is_identical_and_never_requests_neighbourhood():
    p=json.loads((ROOT/'config/experiments/rya1281_orch_continuum_parity.json').read_text())
    new,_=load_git_module(p['orch_commits']['current'],'test_bride_small_current')
    old,_=load_git_module(p['orch_commits']['previous'],'test_bride_small_previous')
    w=np.linspace(-5,5,1001);f=.99+.0001*w;ref=np.ones_like(w)
    def forbidden(level):raise AssertionError('Unnecessary wide-spectrum request')
    a,ac=old.fit_model_guided(w,f,0.,w,ref,exclude_half_width_A=.5)
    b,bc=new.fit_model_guided(w,f,0.,w,ref,exclude_half_width_A=.5,corroborate=forbidden)
    assert a.__dict__==b.__dict__;np.testing.assert_array_equal(ac,bc)
