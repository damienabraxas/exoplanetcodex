"""Source-determined transfer overrides remain separate from original level energies."""
import numpy as np
import pytest
pytest.importorskip('lightweaver')
from lightweaver.atomic_model import AtomicModel,AtomicLevel,VoigtLine,LineType,LinearCoreExpWings
from lightweaver.atomic_table import PeriodicTable
from lightweaver.broadening import LineBroadening,RadiativeBroadening
from pipeline.bride_native_table_control import apply_table_control
from pipeline.wavelength_util import air_to_vac


def fixture_model(tmp_path):
    lines=[VoigtLine(j=1,i=0,f=.1,type=LineType.CRD,quadrature=LinearCoreExpWings(qCore=3,qWing=10,Nlambda=9),broadening=LineBroadening(natural=[RadiativeBroadening(1e7)],elastic=[]))for _ in range(10)]
    model=AtomicModel(element=PeriodicTable['O'],levels=[AtomicLevel(E=0,g=5,label='lower',stage=0),AtomicLevel(E=12860,g=3,label='upper',stage=0)],lines=lines,continua=[],collisions=[])
    path=tmp_path/'O.txt';path.write_text('i j lam_air evi loggf grad gvdw gstk nq qmax q0 iso\n'+'\n'.join(f'{i} 0 1 7774.166 0.0 0.000 1e7 0 0 9 10 3 0'for i in range(10))+'\n');return model,path


def test_centre_only_changes_transfer_frequency_and_einstein_a_without_changing_levels_or_f(tmp_path):
    model,path=fixture_model(tmp_path);prior=model.lines[0];energy=[x.E for x in model.levels];old_A=prior.Aji;old_lam=prior.lambda0
    audit=apply_table_control(model,path,'centres');line=model.lines[0]
    assert len(audit['changes'])==10 and [x.E for x in model.levels]==energy
    assert line.f==prior.f and prior.lambda0==old_lam
    np.testing.assert_allclose(line.lambda0,float(air_to_vac(7774.166))/10,rtol=1e-15)
    np.testing.assert_allclose(line.Aji/old_A,(old_lam/line.lambda0)**2,rtol=1e-14)


def test_strength_only_keeps_transfer_frequency_and_original_line_object(tmp_path):
    model,path=fixture_model(tmp_path);prior=model.lines[0];old_lam=prior.lambda0
    apply_table_control(model,path,'strengths');line=model.lines[0]
    assert line.lambda0==old_lam and prior.f==.1
    assert line.f==pytest.approx(.2)
    np.testing.assert_allclose(line.Aji/prior.Aji,2.)
    path.write_text(path.read_text().replace('0 0 1','0 1 0',1))
    with pytest.raises(ValueError,match='identity'):apply_table_control(model,path,'centres')
