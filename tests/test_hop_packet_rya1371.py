"""Native source-density units and wrong-state refusals, independent of a live HOP run."""
from types import SimpleNamespace
import numpy as np
import pytest
from pipeline.hop_packet import native_planck,compute_packet


def test_native_planck_matches_wavelength_form_in_si_and_safe_wien_tail(tmp_path):
    source=tmp_path/'constants.f90';source.write_text('hh = 6.626176d-27\nbk = 1.380662d-16\ncc = 2.99792458d10\n')
    wave=np.array([200.,500.,777.,15000.]);T=np.array([4000.,6000.]);actual,c=native_planck(wave,T,source)
    h,k,speed=c['hh']*1e-7,c['bk']*1e-7,c['cc']*.01;lam=wave[:,None]*1e-9
    wavelength_form=2*h*speed**2/lam**5/np.expm1(h*speed/(lam*k*T[None,:]));expected=wavelength_form*lam**2/speed
    np.testing.assert_allclose(actual,expected,rtol=3e-14)
    with np.errstate(over='raise',invalid='raise'):
        tail,_=native_planck(np.array([.01]),np.array([4000.]),source)
    np.testing.assert_array_equal(tail,[[0.]])


def test_invalid_thermodynamics_and_wavelengths_refuse_before_native_access(monkeypatch):
    monkeypatch.setenv('OMP_NUM_THREADS','1')
    def forbidden():raise AssertionError('Native API must not be loaded')
    monkeypatch.setattr('pipeline.hop_packet.checked_runtime',forbidden)
    args=[np.array([5000.,6000.]),np.array([1e16,2e16]),np.array([1e-6,2e-6]),np.array([500.,600.]),np.full(92,1.),SimpleNamespace(symbol='O')]
    bad=list(args);bad[0]=np.array([5000.,np.nan])
    with pytest.raises(ValueError,match='thermodynamic'):compute_packet(*bad,vmic_m_s=1000.)
    bad=list(args);bad[3]=np.array([600.,500.])
    with pytest.raises(ValueError,match='vacuum-nm'):compute_packet(*bad,vmic_m_s=1000.)
