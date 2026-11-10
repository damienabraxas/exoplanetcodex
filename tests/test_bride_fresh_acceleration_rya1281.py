"""Fresh nonlinear rates discriminate a stale exact solve from a physical root."""
import numpy as np
from pipeline.bride_fresh_acceleration import select_fresh_population


def evaluator():
    states=[]
    def evaluate(n):
        states.append(n.copy());x=float(n[1,0])
        # Independent two-level example: radiation feedback changes the up rate.
        incoming=(.1+.4*x)*(1-x);outgoing=x
        return dict(residual=abs(incoming-outgoing)/(incoming+outgoing))
    return states,evaluate


def test_stale_zero_residual_is_rejected_and_actual_state_is_restored():
    plain=np.array([[.85],[.15]])
    # At the old state x=.5, up=.3; this is exact for that stale rate matrix.
    x=.3/1.3;proposal=np.array([[1-x],[x]])
    assert abs(.3*proposal[0,0]-proposal[1,0])<1e-15
    states,evaluate=evaluator();chosen,result,info=select_fresh_population(plain,proposal,evaluate)
    assert not info['accepted']and info['restored_plain']
    assert info['proposal_fresh_residual']>info['plain_fresh_residual']
    np.testing.assert_array_equal(chosen,plain)
    np.testing.assert_array_equal(states[-1],chosen)
    assert len(states)==3


def test_independent_nonlinear_root_is_accepted_with_current_radiation_state():
    x=(-.7+np.sqrt(.65))/.8;root=np.array([[1-x],[x]])
    plain=np.array([[.85],[.15]]);states,evaluate=evaluator()
    chosen,result,info=select_fresh_population(plain,root,evaluate)
    assert info['accepted']and not info['restored_plain']
    assert result['residual']<1e-14
    np.testing.assert_array_equal(states[-1],root)
    np.testing.assert_array_equal(chosen,root)
