import numpy as np
import pytest
from pipeline.bride_observables import flux_quadrature,integrate_flux,window_integral


@pytest.mark.parametrize('rays',[3,5,10])
def test_flux_matches_independent_isotropic_and_linear_limb_laws(rays):
    mu,w=flux_quadrature(rays)
    assert integrate_flux(np.full((2,rays),4.),mu,w)==pytest.approx([4*np.pi]*2)
    u=.6
    assert integrate_flux(4*(1-u+u*mu),mu,w)==pytest.approx(4*np.pi*(1-u/3))


def test_flux_is_not_disk_centre_intensity_and_rejects_bad_shapes():
    mu,w=flux_quadrature(5)
    assert integrate_flux(mu,mu,w)==pytest.approx(2*np.pi/3)
    with pytest.raises(ValueError):integrate_flux(np.ones(4),mu,w)


def test_window_area_and_uncovered_input_fail_closed():
    assert window_integral(np.arange(4.),np.full(4,.8),.5,2.5)==pytest.approx(.4)
    with pytest.raises(ValueError):window_integral(np.arange(4.),np.ones(4),-1,2)


def test_corrected_kitt_reader_keeps_published_continuum_and_converts_vacuum(tmp_path):
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.wavelength_util import air_to_vac,vac_to_air
    p=tmp_path/'residual.dat'
    air=np.array([7771.,7772.,7773.]);vac_nm=air_to_vac(air)/10
    p.write_text('KURUCZ RESIDUAL IRRADIANCE ATLAS 2005\n'+
                 '\n'.join(f'{x:.12f} {y}'for x,y in zip(vac_nm,[.99,.7,1.01])))
    w,f=read_kurucz2005_residual(p,7770,7774)
    # Check the loader's delegated conversion/units, not exact invertibility of
    # the canonical approximate air_to_vac formula (round-trip differs ~1e-5 A).
    np.testing.assert_allclose(w,vac_to_air(vac_nm*10),rtol=0,atol=1e-9)
    np.testing.assert_array_equal(f,[.99,.7,1.01])
