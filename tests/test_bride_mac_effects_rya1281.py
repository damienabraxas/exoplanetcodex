"""Exercise the same installed iSpec post-effects used by synthesis; no fitting."""
import sys
import numpy as np
import pytest


@pytest.fixture(scope='module')
def effects():
    from config.constants import require_codex_path
    sys.path.insert(0,str(require_codex_path('engines.ispec')))
    import ispec
    return ispec.apply_post_fundamental_effects


def test_ispec_rv_sign_and_broadening_do_not_manufacture_absorption(effects):
    from config.constants import PHYSICS
    from pipeline.bride_observables import window_integral
    wave=np.arange(7768.,7780.,.0005);step=wave[1]-wave[0]
    f=1-.3*np.exp(-.5*((wave-7774.)/.05)**2)
    segments=np.array([(wave[0]/10,wave[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
    shifted=effects(wave/10,f.copy(),segments,macroturbulence=0,vsini=0,R=0,vrad=(1.,))
    # Exclude the API's zero-filled outer edges when locating the diagnostic line.
    inside=(wave>=7773)&(wave<=7775)
    centre=wave[inside][shifted[inside].argmin()]
    assert centre==pytest.approx(7774*(1+1/PHYSICS['c_kms']),abs=2*step)
    broad=effects(wave/10,f.copy(),segments,macroturbulence=3.8,vsini=1.8,R=500000,vrad=(0.,))
    assert broad[inside].min()>f[inside].min()
    # A two-endpoint grid-area bound, not a physical tolerance fit to the Sun.
    assert abs(window_integral(wave,broad,7773,7775)-window_integral(wave,f,7773,7775))<2*step*.3


def test_orchestrator_keeps_the_corrected_pre_normalized_kitt_product():
    from pipeline.telluric_policy import gate_holding,holding_basis,applied_state
    from pipeline.band_policy import check_intake
    holding='solar_kpno_kurucz2005_corrected';instrument='kpno_solar_atlas'
    assert gate_holding(holding,instrument)[0]
    assert holding_basis(holding,instrument)=='corrected'
    assert applied_state(holding)=='applied'
    policy=check_intake(7774.166,'synthesis',instrument=instrument)
    assert 'atlas continuum trusted' in policy.continuum_treatment


def test_air_comparison_is_not_vacuum_or_double_converted(tmp_path):
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.wavelength_util import air_to_vac,vac_to_air
    air=np.array([7771.944,7774.166,7775.388])
    vacuum=air_to_vac(air)
    path=tmp_path/'irradrelwl.dat'
    path.write_text('KURUCZ RESIDUAL IRRADIANCE ATLAS 2005\n'+
        '\n'.join(f'{v/10:.12f} 0.9'for v in vacuum))
    actual,_=read_kurucz2005_residual(path,7770,7777)
    # Agreement at the published three-decimal air-line precision; the canonical
    # inverse approximation itself has an ~1e-5 A round-trip remainder.
    np.testing.assert_allclose(actual,air,rtol=0,atol=.001)
    assert np.all(vacuum-actual>2)
    assert np.all(actual-vac_to_air(actual)>2)
