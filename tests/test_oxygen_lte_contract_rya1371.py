import numpy as np
import pytest
from pipeline.oxygen_lte_contract import normalized_weights,native_lte_weights,held_lw_weights
C=dict(hh=6.626176e-27,bk=1.380662e-16,ee=1.602189e-12,em=9.109534e-28,ec=4.80325e-10,pi=3.14159265359)


def test_two_stage_saha_boltzmann_and_normalization_independent_formula():
    T=np.array([4000.,6000.]);ne=np.array([1e10,1e13]);E=np.array([0.,1.,13.6]);g=np.array([2.,4.,1.]);stage=[1,1,2]
    log,aux=native_lte_weights(E,g,stage,T,ne,C,round_native=False,shield=False)
    actual=normalized_weights(log,np.array([1e8,2e8]))
    expected=2*(2*np.pi*C['em']*C['bk']*T/C['hh']**2)**1.5/ne*.5*np.exp(-13.6*C['ee']/(C['bk']*T))
    np.testing.assert_allclose(actual[2]/actual[0],expected,rtol=2e-13)
    np.testing.assert_allclose(actual[1]/actual[0],2*np.exp(-C['ee']/(C['bk']*T)),rtol=1e-14)
    np.testing.assert_allclose(actual.sum(axis=0),[1e8,2e8],rtol=2e-16)
    assert np.all(aux['dxi_eV']==0)


def test_debye_and_precision_axes_are_discriminating():
    args=([0.,13.6],[2.,1.],[1,2],[5000.,6000.],[1e13,1e16],C)
    log,aux=native_lte_weights(*args);bare,_=native_lte_weights(*args,shield=False);full,_=native_lte_weights(*args,round_native=False)
    assert np.all(log[1]>bare[1])
    assert np.max(abs(log-full))>1e-8
    for value in aux.values():np.testing.assert_array_equal(value,value.astype(np.float32).astype(float))


def test_invalid_population_and_stage_states_refuse():
    with pytest.raises(ValueError):normalized_weights([[0.,np.nan]],[1.,1.])
    with pytest.raises(ValueError):normalized_weights([[0.]],[0.])
    with pytest.raises(ValueError):native_lte_weights([0.,10.],[1.,1.],[1,3],[5000.],[1e10],C)


def test_independent_held_lw_formula_matches_actual_library_including_debye():
    lw=pytest.importorskip('lightweaver');from lightweaver import rh_atoms,atomic_set,constants as c
    atom=rh_atoms.H_6_atom();T=np.array([4500.,6000.]);ne=np.array([1e16,1e19]);total=np.array([1e20,2e20]);E=[l.E_SI for l in atom.levels];g=[l.g for l in atom.levels];stage=np.array([l.stage+1 for l in atom.levels])
    constants=dict(h=c.HPlanck,k=c.KBoltzmann,m=c.MElectron,q=c.QElectron,epsilon=c.Epsilon0)
    for shield in [False,True]:
        actual=normalized_weights(held_lw_weights(E,g,stage,T,ne,constants,shield=shield),total)
        expected=atomic_set.lte_pops(atom,T,ne,total,debye=shield)
        np.testing.assert_allclose(actual,expected,rtol=5e-14,atol=0)
