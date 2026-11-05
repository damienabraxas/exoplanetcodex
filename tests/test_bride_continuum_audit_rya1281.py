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
    actual=np.polyfit(w[trace['fit_indices']],ratio[trace['fit_indices']],1)
    np.testing.assert_array_equal(actual,trace['passes'][-1]['coefficients'])
    frozen=np.zeros(len(w),bool);frozen[trace['surviving_indices']]=True
    diagnostic,_=trace_fit(w,ratio,0.,np.ones_like(w),frozen,clip=False)
    assert diagnostic['level']-trace['level']==pytest.approx(-.000990121732,abs=1e-12)
    assert not diagnostic['postfit_survivors_differ']
