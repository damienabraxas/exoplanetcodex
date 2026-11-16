import numpy as np
import pytest
from pipeline.saved_field_quadrature import gaussian_panel_batches,interpolate_saved_field


def test_positive_gaussian_weights_integrate_known_moments_on_uneven_panels():
    panels=np.array([1.,1.1,1.8,3.])
    for order in [4,8,16]:
        batches=list(gaussian_panel_batches(panels,order,batch_panels=2))
        assert all(len(q)<=2*order and np.all(w>0)for q,w in batches)
        for degree in range(2*order):
            # Rescale to unit domain to avoid large powers obscuring errors.
            actual=sum(np.sum(w*((q-1)/2)**degree)for q,w in batches)
            np.testing.assert_allclose(actual,2/(degree+1),rtol=2e-14,atol=1e-15)


def test_interpolated_field_preserves_nodes_and_depth_axes_without_extrapolation():
    x=np.array([1.,2.,3.]);J=np.array([[1.,10.],[3.,5.],[2.,8.]])
    np.testing.assert_array_equal(interpolate_saved_field(x,J,x),J)
    np.testing.assert_array_equal(interpolate_saved_field(x,J,[1.5,2.5]),[[2.,7.5],[2.5,6.5]])
    with pytest.raises(ValueError):interpolate_saved_field(x,J,[.9])
    with pytest.raises(ValueError):interpolate_saved_field(x,-J,[2.])


@pytest.mark.parametrize('panels,order,batch',[([1.,1.],4,2),([1.,2.],True,2),([1.,2.],16,2000),([1.,2.],8,0)])
def test_invalid_gaussian_contract_refuses(panels,order,batch):
    with pytest.raises(ValueError):list(gaussian_panel_batches(panels,order,batch_panels=batch))
