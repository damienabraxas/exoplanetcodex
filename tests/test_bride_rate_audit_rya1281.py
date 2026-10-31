import numpy as np
import pytest
from pipeline.bride_rate_audit import process_flow


def test_destination_source_flow_balance_and_wrong_direction_canary():
    p=np.array([[10.,20.],[2.,4.]])
    r=np.zeros((2,2,2));r[1,0]=2.;r[0,1]=10.
    net,gross=process_flow(r,p)
    np.testing.assert_array_equal(net,0.)
    assert np.all(gross>0)
    bad=r.swapaxes(0,1)
    wrong,_=process_flow(bad,p)
    assert np.max(abs(wrong))>0
    np.testing.assert_array_equal(wrong.sum(axis=0),0.)
    r[0,0]=-99
    np.testing.assert_array_equal(process_flow(r,p)[0],net)
    r[0,1]=-1
    with pytest.raises(ValueError):process_flow(r,p)
