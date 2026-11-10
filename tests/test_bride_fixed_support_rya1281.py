import numpy as np
from pipeline.bride_fixed_support import shared_masks,fixed_fit


def test_bright_global_intersection_can_fail_while_per_side_support_brackets():
 w=np.linspace(-5,5,1001);s=1+.001*w;i=1+.002*w
 masks=shared_masks(w,0.,s,i,inner=.5,outer=5.,percentile=90.)
 args=dict(maximum_shift=.15,minimum_pixels=8,minimum_each_side=8,require_bracket=True)
 r,c=fixed_fit(w,i,0.,i,masks['shared_top_decile_intersection'],**args)
 assert not r['applied']and r['reason']=='REFUSED_BRACKETING'and c is None
 before=[x.copy()for x in [w,i,masks['balanced_iag_per_side']]]
 r,c=fixed_fit(w,i,0.,i,masks['balanced_iag_per_side'],**args)
 assert r['applied']and r['blue_count']>=8 and r['red_count']>=8 and r['geometry']['centre_outside_hull_A']==0
 np.testing.assert_allclose(c,1,atol=1e-14)
 for actual,old in zip([w,i,masks['balanced_iag_per_side']],before):np.testing.assert_array_equal(actual,old)


def test_fixed_coefficients_recover_known_slope_and_keep_bound_refusal():
 w=np.linspace(-5,5,1001);ref=1-.001*np.cos(w);mask=(abs(w)>1)&(abs(w)<4);args=dict(maximum_shift=.15,minimum_pixels=8,minimum_each_side=8,require_bracket=True)
 r,c=fixed_fit(w,ref*(1.02+.001*w),0.,ref,mask,**args)
 assert r['applied']and r['indices']==np.flatnonzero(mask).tolist()
 np.testing.assert_allclose([r['level'],r['slope_per_A']],[1.02,.001],rtol=0,atol=1e-14)
 r,c=fixed_fit(w,ref*1.2,0.,ref,mask,**args)
 assert not r['applied']and r['reason']=='REFUSED_PINNED_LEVEL_BOUND'and c is None
