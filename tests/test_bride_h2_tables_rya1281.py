import numpy as np
import pytest


def table():
    from config.constants import require_codex_path
    from pipeline.bride_h2_tables import H2Tables
    return H2Tables(require_codex_path('experiments.bride_molecular_tables'))


def test_published_knots_units_and_strict_domain():
    from scipy.constants import k
    t=table()
    np.testing.assert_allclose(t.phi(t.temperature),k*t.temperature/10**t.log10_Kp_Pa,rtol=5e-14)
    for value in (np.nan,0.,np.nextafter(1000.,-np.inf),np.nextafter(10000.,np.inf)):
        with pytest.raises(ValueError):t.phi(value)
    assert t.phi(5000.)==pytest.approx(1.742421763773e-26,rel=1e-11)


def test_network_independent_root_and_molecular_limit():
    from scipy.optimize import brentq
    t=table()
    for T in (1000.,4062.,5772.,9930.,10000.):
        for N in (1e12,1e20,1e28):
            for f,minus in ((1.,0.),(.95,1e-8),(.01,.001)):
                got=t.hydrogen_network(T,N,f,minus);phi=t.phi(T)
                root=brentq(lambda x:(1/f+minus)*x+2*phi*N*x*x-1,0.,1.,xtol=1e-15)
                assert got['HI_m3']/N==pytest.approx(root,abs=2e-15)
                assert abs(got['nuclei_residual'])<1e-14
                assert got['H2_m3']/N<=.5
    for N,f,m in ((0.,1.,0.),(1.,0.,0.),(1.,1.01,0.),(1.,1.,-1.)):
        with pytest.raises(ValueError):t.hydrogen_network(5000.,N,f,m)


def test_source_tampering_is_refused(tmp_path):
    import shutil
    from config.constants import require_codex_path
    from pipeline.bride_h2_tables import H2Tables
    source=require_codex_path('experiments.bride_molecular_tables')
    for file in source.iterdir():shutil.copyfile(file,tmp_path/file.name)
    with (tmp_path/'table7_vNov2022.dat').open('a') as f:f.write('\n')
    with pytest.raises(ValueError,match='Changed chemistry source'):H2Tables(tmp_path)
