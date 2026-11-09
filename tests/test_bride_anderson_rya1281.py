import numpy as np
from pipeline.bride_anderson import LogPopulationAnderson


def test_accelerates_slow_positive_conservative_map_and_keeps_its_root():
    root=np.array([[.7,.2],[.3,.8]]);total=np.array([4.,9.])
    n=np.array([[.6,.3],[.4,.7]])*total
    plain=n.copy();a=LogPopulationAnderson(memory=3,delay=2,rcond=1e-10,max_log_step=.5)
    for _ in range(12):
        plain=.98*plain+.02*root*total
        direct=.98*n+.02*root*total
        n,info=a.update(n,direct,total)
        assert np.all(n>0)
        np.testing.assert_allclose(n.sum(axis=0),total,rtol=1e-14)
    assert np.max(abs(n/total-root))<1e-7
    assert np.max(abs(plain/total-root))>.05
    held=root*total;out,_=a.update(held,held,total)
    np.testing.assert_allclose(out,held,rtol=1e-12,atol=1e-12)


def test_oversized_proposal_returns_exact_direct_map_without_clipping():
    a=LogPopulationAnderson(memory=2,delay=0,rcond=1e-10,max_log_step=1e-12)
    n=np.array([[.6],[.4]]);g=.99*n+.01*np.array([[.7],[.3]])
    a.update(n,g,[1.]);h=.99*g+.01*np.array([[.7],[.3]])
    out,info=a.update(g,h,[1.])
    assert info['attempted']and not info['accepted']
    np.testing.assert_array_equal(out,h)
