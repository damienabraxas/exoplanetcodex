"""Non-square independent binary fixtures discriminate axes, sentinels and units."""
from pathlib import Path
import numpy as np
import pytest
from pipeline.openmulti3d_output import read_rate_pair,read_patch_set,interpolate_positive
from pipeline.bride_rate_audit import process_flow


def test_native_wavelength_density_conversion_matches_independent_planck_laws():
    from pipeline.openmulti3d_output import per_angstrom_to_per_hz
    wave=np.array([6158.,7774.,8446.]);meters=wave*1e-10;nu=299792458./meters
    h=6.62607015e-34;k=1.380649e-23;T=5772.
    frequency_SI=2*h*nu**3/299792458.**2/np.expm1(h*nu/(k*T))
    wavelength_SI_per_m=2*h*299792458.**2/meters**5/np.expm1(h*299792458./(meters*k*T))
    native_cgs_per_A=wavelength_SI_per_m*1e-7
    np.testing.assert_allclose(per_angstrom_to_per_hz(wave,native_cgs_per_A)*1e-3,frequency_SI,rtol=2e-15)


def write_fixture(tmp_path,*,sentinel_rate=0.):
    n=np.array([[1.,2.],[3.,6.],[7.,14.]]);up=np.array([[3.,6.],[7.,14.],[sentinel_rate,sentinel_rate]]);down=np.array([[1.,2.],[1.,2.],[0.,0.]])
    r=tmp_path/'rates.bin';c=tmp_path/'coll.bin'
    with r.open('wb')as f:
        np.array([3,3,1,1,2,5],'<i4').tofile(f);np.array([1,1,3],'<i4').tofile(f);np.array([2,3,0],'<i4').tofile(f)
        for a in [n,n]:a.T.astype('<f8').ravel(order='F').tofile(f)
        for a in [up,down]:a.T.astype('<f4').ravel(order='F').tofile(f)
    with c.open('wb')as f:np.array([3,1,1,2],'<i4').tofile(f);np.array([[.3,.7,2.1],[.6,1.4,4.2]],'<f4').ravel(order='F').tofile(f)
    return r,c,n


def test_rate_layout_and_detailed_balance_with_explicit_inactive_top_ion(tmp_path):
    r,c,n=write_fixture(tmp_path);d=read_rate_pair(r,c);np.testing.assert_array_equal(d['populations_cm3'],n);assert list(d['inactive_sentinel_slots'])==[2]
    net,gross=process_flow(d['collision_rates_s_inv']+d['radiative_rates_s_inv'],n)
    np.testing.assert_allclose(net,0,atol=1e-13);assert np.all(gross>0)
    assert d['radiative_rates_s_inv'][1,0,1]==6 and d['radiative_rates_s_inv'][0,1,1]==2


def test_sentinel_with_nonzero_rate_and_truncated_data_are_refused(tmp_path):
    r,c,_=write_fixture(tmp_path,sentinel_rate=1.)
    with pytest.raises(ValueError,match='invalid/sentinel'):read_rate_pair(r,c)
    r.write_bytes(r.read_bytes()[:-1])
    with pytest.raises(ValueError,match='Truncated'):read_rate_pair(r,c)


def test_patch_axes_and_byte_count_are_explicit(tmp_path):
    a=np.arange(12.).reshape(2,3,2);p=tmp_path/'atmos_001000.bin';a.astype('<f4').ravel(order='F').tofile(p);(tmp_path/'atmos_meta.txt').write_text('* metadata\natmos_001000.bin 32 2 3 2\n');d=read_patch_set(tmp_path,'atmos');np.testing.assert_array_equal(d[1,0],a)
    p.write_bytes(p.read_bytes()[:-4])
    with pytest.raises(ValueError,match='byte-count'):read_patch_set(tmp_path,'atmos')


def test_reference_mapping_forbids_physical_extrapolation_and_retains_positivity():
    x=np.array([0.,2.,4.]);y=np.array([[1.,4.,16.],[2.,8.,32.]]);before=y.copy();r,receipt=interpolate_positive(x,y,[1.,3.]);np.testing.assert_allclose(r,[[2.,8.],[4.,16.]]);np.testing.assert_array_equal(y,before)
    assert receipt['max_endpoint_rounding_clamp_m']==0
    with pytest.raises(ValueError,match='extrapolation'):interpolate_positive(x,y,[-1.])
