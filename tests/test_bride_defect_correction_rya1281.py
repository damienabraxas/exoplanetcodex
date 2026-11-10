import numpy as np
from pipeline.bride_defect_correction import corrected_populations

def test_preconditioner_cannot_move_a_physical_fixed_point():
    # Deliberately different physical and approximate stationary distributions.
    r=np.array([[-2.,3.],[2.,-3.]])[:,:,None]
    p=np.array([[-7.,1.],[7.,-1.]])[:,:,None]
    n=np.array([[3.],[2.]])
    assert np.max(abs(np.einsum('ijd,jd->id',p,n)))>1
    np.testing.assert_array_equal(corrected_populations(p,r,n,[5.]),n)

def test_physical_operator_step_recovers_solution_and_conserves_each_depth():
    r=np.repeat(np.array([[-2.,3.],[2.,-3.]])[:,:,None],2,axis=2)
    n=np.array([[1.,8.],[4.,2.]])
    before=n.copy();out=corrected_populations(r,r,n,[5.,10.])
    np.testing.assert_allclose(out,[[3.,6.],[2.,4.]])
    np.testing.assert_allclose(np.einsum('ijd,jd->id',r,out),0,atol=1e-14)
    np.testing.assert_array_equal(n,before)
