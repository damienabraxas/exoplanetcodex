"""Behavioral canaries for native collider roles, packet axes and SE refusal."""
import numpy as np
import pytest
from types import SimpleNamespace as NS
pytest.importorskip('lightweaver')
from lightweaver.atomic_table import PeriodicTable
from pipeline.bride_hop_transfer import Perturbers,NativeHopBackground


def test_native_perturber_roles_distinguish_neutral_and_proton_without_fake_levels():
    p={'H_I_m3':np.array([10.,20.]),'H_II_m3':np.array([1.,2.]),'H_-_m3':np.array([.01,.02]),'He_I_m3':np.array([3.,8.])};v=Perturbers(p)
    np.testing.assert_array_equal(v['H'][0,:],[10.,20.]);np.testing.assert_array_equal(v['H'][-1,:],[1.,2.]);np.testing.assert_array_equal(v['H-'],[.01,.02])
    np.testing.assert_allclose(v.abundance[PeriodicTable[2]],[.3,.4])
    with pytest.raises(ValueError,match='excited-state'):v['H'][1,:]
    with pytest.raises(ValueError,match='Unqualified'):v['C']


def test_background_purpose_refuses_coupled_se_before_any_native_call(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('Native packet must not be requested')
    monkeypatch.setattr('pipeline.bride_hop_transfer.compute_packet',forbidden)
    with pytest.raises(ValueError,match='no coupled SE'):
        NativeHopBackground(None,NS(activeAtoms=[object()],detailedAtoms=[]),np.array([500.]),state={},atom=None,abundances=[],records=[])


def test_nonsquare_packet_copy_and_changed_thermodynamic_state_refusal(monkeypatch):
    state=dict(temperature_K=np.array([5000.,6000.]),ne_m3=np.array([1e16,2e16]),rho_kg_m3=np.array([1e-6,2e-6]),nHTot_m3=np.array([1e20,2e20]),vturb_m_s=np.array([1000.,1000.]),height_m=np.array([100.,0.]))
    atm=NS(**{name:state[key].copy()for name,key in [('temperature','temperature_K'),('ne','ne_m3'),('nHTot','nHTot_m3'),('vturb','vturb_m_s'),('height','height_m')]})
    values=np.arange(1.,7.).reshape(3,2);packet=dict(nHTot_m3=state['nHTot_m3'].copy(),chi_m_inv=values,eta_SI=values*3,scat_m_inv=values*.1)
    monkeypatch.setattr('pipeline.bride_hop_transfer.compute_packet',lambda *args,**kw:(packet,{'purpose':'synthetic axis fixture only'}))
    bg=NativeHopBackground(NS(atmosphere=atm),NS(activeAtoms=[],detailedAtoms=[NS(element=PeriodicTable['O'])]),np.array([500.,600.,700.]),state=state,atom=None,abundances=[],records=[])
    outputs=[np.zeros((3,2))for _ in range(3)];bg.compute_background(NS(pyAtmos=atm),*outputs)
    for actual,key in zip(outputs,['chi_m_inv','eta_SI','scat_m_inv']):np.testing.assert_array_equal(actual,packet[key])
    atm.ne[0]*=2
    with pytest.raises(ValueError,match='Changed bound HOP state ne'):bg.compute_background(NS(pyAtmos=atm),*outputs)


def test_miswired_helium_control_refuses_actual_he_on_reads():
    from pipeline.bride_hop_transfer import verify_helium_reads
    p={'H_I_m3':np.array([10.,20.]),'He_I_m3':np.array([3.,8.])}
    on=Perturbers(p);on.abundance[PeriodicTable[2]]
    with pytest.raises(ValueError,match='did not reach'):verify_helium_reads(on.reads,False)
    off=Perturbers(p,helium=False)
    np.testing.assert_array_equal(off.abundance[PeriodicTable[2]],[0.,0.])
    verify_helium_reads(off.reads,False)
    np.testing.assert_array_equal(p['He_I_m3'],[3.,8.])
    with pytest.raises(ValueError,match='did not reach'):verify_helium_reads(off.reads,True)
