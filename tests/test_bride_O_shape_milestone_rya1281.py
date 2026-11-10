import numpy as np
from scripts.analyze_bride_O_shape_milestone import shape_metrics


def test_shape_statistics_cancel_common_scale_but_retain_core_changes():
 w=np.linspace(-1,1,2001);f=1-.2*np.exp(-(w/.12)**2);kw=dict(core_half_A=.1,feature_half_A=.5,shoulder_inner_A=.35,shoulder_outer_A=.55)
 a=shape_metrics(w,f,0.,**kw);b=shape_metrics(w,f*.92,0.,**kw)
 np.testing.assert_allclose([a['contrast'],a['relative_feature_area_mA']],[b['contrast'],b['relative_feature_area_mA']],rtol=0,atol=1e-10)
 deeper=shape_metrics(w,1-.25*np.exp(-(w/.12)**2),0.,**kw)
 assert deeper['contrast']>a['contrast']and deeper['relative_feature_area_mA']>a['relative_feature_area_mA']
