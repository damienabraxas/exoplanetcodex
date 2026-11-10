import numpy as np
import pytest
from pipeline.bride_crossed_continuum import crossed_fits,finite_accounting


def test_finite_paths_interaction_closure_and_refusal():
    corners={k:dict(area=v) for k,v in zip(('before_before','before_after','after_before','after_after'),(1.,3.,2.,7.))}
    r=finite_accounting(corners,'area')
    assert r['total']==6 and r['support_path_average']==2.5 and r['denominator_path_average']==3.5
    assert r['interaction']==3 and r['closure']==0
    assert r['support_first']+r['denominator_after_support']==r['total']
    assert r['denominator_first']+r['support_after_denominator']==r['total']
    corners['before_after']['area']=None
    assert finite_accounting(corners,'area')==dict(available=False,reason='REFUSED_CORNER')
    corners['before_after']['area']=np.nan
    with pytest.raises(ValueError):finite_accounting(corners,'area')


def test_fixed_support_denominator_effect_and_identity_keep_arrays_immutable():
    w=np.linspace(-5,5,1001);obs=.997+.0002*w
    refs={'before':np.ones_like(w),'after':np.full_like(w,1.001)}
    masks={'before':w<-1,'after':w>1};old=[a.copy()for a in (w,obs,*refs.values(),*masks.values())]
    results=crossed_fits(w,obs,0.,refs,masks,clip=False)
    for key,(trace,cont) in results.items():
        factor=1.001 if key.endswith('_after')else 1.
        assert trace['level']==pytest.approx(.997/factor,abs=1e-14)
        np.testing.assert_allclose(cont,obs/factor,rtol=0,atol=1e-14)
        assert trace['clipping_stop_reason']=='CLIPPING_DISABLED'
    for before,after in zip(old,(w,obs,*refs.values(),*masks.values())):np.testing.assert_array_equal(before,after)
    masks['before']=masks['before'].astype(int)
    with pytest.raises(ValueError):crossed_fits(w,obs,0.,refs,masks,clip=False)


def test_selection_clipping_differs_from_frozen_coefficient_support():
    w=np.linspace(-5,5,1001);obs=np.full_like(w,.997);obs[::97]+=.01
    refs={'before':np.ones_like(w),'after':np.ones_like(w)}
    masks={'before':abs(w)>1,'after':abs(w)>1.5}
    dynamic=crossed_fits(w,obs,0.,refs,masks,clip=True)
    fixed=crossed_fits(w,obs,0.,refs,masks,clip=False)
    assert dynamic['before_before'][0]['coefficient_support_count']<fixed['before_before'][0]['coefficient_support_count']
    assert dynamic['before_before'][0]['level']==pytest.approx(.997)
    assert fixed['before_before'][0]['level']>.997
    for corner in dynamic:
        support=corner.split('_')[0]
        assert set(dynamic[corner][0]['fit_indices'])<=set(np.flatnonzero(masks[support]))
