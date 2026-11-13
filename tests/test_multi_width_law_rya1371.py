"""Source-law density/T/micro roles and fail-closed source-domain/prescription guards."""
from types import SimpleNamespace as NS
import numpy as np
import pytest
pytest.importorskip('lightweaver')
from pipeline.multi_width_law import OriginalMultiVdwLaw,WidthRecord

def record():return WidthRecord(0,1,2,1.,1e7,2e-15,.3,'legacy-Unsold')

def test_source_law_uses_total_H_with_density_T_response_and_no_level_or_He_lookup():
    class ForbiddenLevels:
        def __getitem__(self,key):raise AssertionError('Width law must not infer H levels/extra He')
    r=record();density=np.array([1e20,2e20]);w=OriginalMultiVdwLaw(r,density,(4000.,10000.));w.setup(NS(i=0,j=1,atom=NS(element=NS(Z=8))))
    atm=NS(temperature=np.array([5000.,6000.]),vturb=np.array([250.,250.]));g=w.broaden(atm,ForbiddenLevels())
    np.testing.assert_allclose(g,2e-15*density*atm.temperature**.3)
    half=OriginalMultiVdwLaw(r,density/2,(4000.,10000.));np.testing.assert_array_equal(half.broaden(atm,None),g/2)
    atm.vturb[:]=4000;np.testing.assert_array_equal(w.broaden(atm,None),g)
    atm.temperature+=10;np.testing.assert_allclose(w.broaden(atm,None)/g,(atm.temperature/(atm.temperature-10))**.3)
    assert w.audit['extra_helium_factor']is False
    with pytest.raises(ValueError):w.neutral_H_m3[0]=5


def test_source_law_refuses_wrong_species_transition_negative_density_and_numeric_domain():
    w=OriginalMultiVdwLaw(record(),[1e20,2e20],(4000.,10000.))
    for line in [NS(i=1,j=2,atom=NS(element=NS(Z=8))),NS(i=0,j=1,atom=NS(element=NS(Z=26)))]:
        with pytest.raises(ValueError,match='model-line'):w.setup(line)
    for t in [[3999.,5000.],[5000.,float('nan')],[5000.]]:
        with pytest.raises(ValueError,match='numeric T/shape'):w.broaden(NS(temperature=np.array(t)),None)
    with pytest.raises(ValueError,match='neutral-H'):OriginalMultiVdwLaw(record(),[-1.,2.],(4000.,10000.))


def test_context_refuses_double_He_or_width_omission_before_any_native_packet(monkeypatch):
    from pipeline.bride_hop_transfer import make_context
    def forbidden(*a,**kw):raise AssertionError('Must refuse before native packet')
    monkeypatch.setattr('pipeline.bride_hop_transfer.compute_packet',forbidden)
    for h,v in [(True,True),(False,False)]:
        with pytest.raises(ValueError,match='owns full VDW/He'):
            make_context(None,None,None,None,None,extra_vac_nm=None,helium_broadening=h,vdw_broadening=v,original_width_law=object())
