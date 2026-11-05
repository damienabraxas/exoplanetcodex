import numpy as np
import pytest


def test_trace_matches_native_estimator_and_keeps_initial_inputs_immutable():
    from pipeline import local_continuum as lc
    from pipeline.bride_continuum_audit import selection,trace_fit
    w=np.arange(4994.,5006.,.01);centre=5000.
    model=1-.006*(1+np.sin(w*2))
    obs=model*(.997+.0002*(w-centre));old=obs.copy()
    mask,_=selection(w,obs,centre,model)
    trace,cont=trace_fit(w,obs,centre,model,mask)
    rec,native=lc.fit_model_guided(w,obs,centre,w,model,exclude_half_width_A=.5)
    assert trace['applied']==rec.applied and trace['reported_n_bins']==rec.n_bins
    assert trace['level']==rec.level_at_centre and trace['slope_per_A']==rec.slope_per_A
    np.testing.assert_array_equal(cont,native);np.testing.assert_array_equal(obs,old)


def test_refusal_and_fixed_support_extrapolation_are_explicit():
    from pipeline.bride_continuum_audit import trace_fit,support_geometry
    w=np.linspace(-5,5,1001);ref=np.ones_like(w)
    mask=(w<-4)
    trace,cont=trace_fit(w,.92*ref,0.,ref,mask)
    assert not trace['applied'] and cont is None
    g=support_geometry(w,0.,np.flatnonzero(mask))
    assert g['blue_count']==100 and g['red_count']==0 and g['centre_outside_hull_A']>4
    assert g['centre_weight_l1']>1
    small=np.zeros(len(w),bool);small[:3]=True
    trace,cont=trace_fit(w,ref,0.,ref,small)
    assert trace['reason']=='INSUFFICIENT_ORIGINAL_ANCHORS' and cont is None


def test_factorial_accounting_is_exact_and_includes_interaction():
    from pipeline.bride_continuum_audit import two_factor_changes
    r=two_factor_changes({k:dict(level=v) for k,v in zip(('synthesis_synthesis','synthesis_iag','iag_synthesis','iag_iag'),(1.,3.,2.,7.))},'level')
    assert r['selector']==2.5 and r['denominator_plus_induced_clipping']==3.5 and r['total']==6 and r['closure']==0
    assert r['interaction']==3


def test_archived_ratio_exposes_cap_support_without_adopting_refit():
    """Derived ratios reproduce the discovered cap defect; no raw atlas needed."""
    import csv
    from pathlib import Path
    from pipeline.bride_continuum_audit import trace_fit
    path=Path(__file__).resolve().parents[1]/'data/results/rya1281/continuum_anchor_v1/pixel_diagnostics.csv'
    with path.open() as stream:
        rows=[r for r in csv.DictReader(stream) if r['case']=='canonical_refined'
              and r['centre_air_A']=='7774.166' and r['reference']=='synthesis']
    w=np.array([float(r['offset_air_A']) for r in rows])
    ratio=np.array([float(r['relative_ratio']) for r in rows])
    mask=np.array([r['initially_selected']=='True' for r in rows])
    trace,_=trace_fit(w,ratio,0.,np.ones_like(w),mask)
    assert len(trace['passes'])==5 and trace['postfit_survivors_differ']
    assert len(trace['fit_indices'])==73 and len(trace['surviving_indices'])==70
    assert trace['coefficient_support_count']==73 and trace['postfit_survivor_count']==70
    assert trace['clipping_candidate_count']==70 and trace['clipping_pass_count']==5
    assert trace['clipping_limit_reached'] and trace['clipping_stop_reason']=='ITERATION_LIMIT_AFTER_CLIP'
    actual=np.polyfit(w[trace['fit_indices']],ratio[trace['fit_indices']],1)
    np.testing.assert_array_equal(actual,trace['passes'][-1]['coefficients'])
    frozen=np.zeros(len(w),bool);frozen[trace['surviving_indices']]=True
    diagnostic,_=trace_fit(w,ratio,0.,np.ones_like(w),frozen,clip=False)
    assert diagnostic['level']-trace['level']==pytest.approx(-.000990121732,abs=1e-12)
    assert not diagnostic['postfit_survivors_differ']
    assert not diagnostic['clipping_limit_reached']
    assert diagnostic['clipping_stop_reason']=='CLIPPING_DISABLED'


def test_clipping_rejection_keeps_fit_support_and_candidate_distinct():
    from pipeline.bride_continuum_audit import trace_fit
    w=np.linspace(-5,5,20);ratio=np.ones(20);ratio[-1]=1.02
    trace,_=trace_fit(w,ratio,0.,np.ones(20),np.ones(20,bool))
    assert trace['clipping_stop_reason']=='CLIP_WOULD_UNDERRUN_MINIMUM'
    assert trace['coefficient_support_count']==trace['postfit_survivor_count']==20
    assert trace['clipping_candidate_count']<20
    assert not trace['clipping_limit_reached'] and not trace['postfit_survivors_differ']


def test_retained_and_insufficient_support_do_not_claim_clipping_limit():
    from pipeline.bride_continuum_audit import trace_fit
    w=np.linspace(-5,5,100);one=np.ones(100)
    trace,_=trace_fit(w,one,0.,one,np.ones(100,bool))
    assert trace['clipping_stop_reason']=='ALL_RETAINED'
    assert trace['coefficient_support_count']==trace['postfit_survivor_count']==100
    assert trace['clipping_candidate_count']==100 and not trace['clipping_limit_reached']
    trace,cont=trace_fit(w,one,0.,one,w<-4)
    assert trace['clipping_stop_reason']=='INSUFFICIENT_ORIGINAL_ANCHORS'
    assert trace['coefficient_support_count']==0 and trace['clipping_pass_count']==0
    assert trace['clipping_candidate_count'] is None and cont is None


def test_explicit_support_archive_preserves_all_prior_scientific_fields():
    import json
    from collections import Counter
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]/'data/results/rya1281'
    old=json.loads((root/'continuum_anchor_v1/report.json').read_text())
    new=json.loads((root/'continuum_support_v1/report.json').read_text())
    added={'coefficient_support_count','postfit_survivor_count','clipping_candidate_count',
           'clipping_pass_count','clipping_stop_reason','clipping_limit_reached'}
    stops=Counter()
    def compare(a,b,path=''):
        if path=='/adapter':return  # This audit adapter changed; native code did not.
        if isinstance(a,dict):
            assert set(a)<=set(b) and set(b)-set(a)<=added
            if 'clipping_stop_reason' in b:
                stops[b['clipping_stop_reason']]+=1
                if b['clipping_limit_reached']:
                    assert b['coefficient_support_count']==73 and b['postfit_survivor_count']==70
            for key,value in a.items():compare(value,b[key],path+'/'+key)
        elif isinstance(a,list):
            assert len(a)==len(b)
            for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
        else:assert a==b,path
    compare(old,new)
    assert stops=={'ALL_RETAINED':38,'CLIPPING_DISABLED':26,
                   'INSUFFICIENT_ORIGINAL_ANCHORS':6,'ITERATION_LIMIT_AFTER_CLIP':4}
