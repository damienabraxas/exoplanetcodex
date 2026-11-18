import numpy as np
import pytest
from types import SimpleNamespace
from lightweaver.atomic_model import LinearCoreExpWings
from pipeline.native_line_grid import NativeFractionGrid,NativeLineGridControl


class Line:
    def __init__(self):
        self.i=0;self.j=1;self.f=.1;self.lambda0=500.
        self.quadrature=LinearCoreExpWings(qCore=1.25,qWing=10.,Nlambda=9);self.quadrature.setup(self)
    def wavelength(self):return self.quadrature.wavelength(self)


@pytest.mark.parametrize('mode',['native-spacing-and-support','native-support-only'])
def test_numerical_control_preserves_centres_strengths_and_matches_source_endpoints(mode):
    model=SimpleNamespace(lines=[Line()for _ in range(256)]);offsets=np.linspace(-2e-4,2e-4,17)
    packet=dict(scope='synthetic numerical control',grids=[dict(index=k,i=0,j=1,f=.1,offsets_fraction=offsets.tolist())for k in range(256)])
    audit=NativeLineGridControl(packet,{},mode=mode).apply(model)
    for line in model.lines:
        assert line.lambda0==500. and line.f==.1
        wave=line.wavelength();assert len(wave)==17
        np.testing.assert_allclose(wave[[0,-1]],500*(1+offsets[[0,-1]]),rtol=0,atol=1e-13)
        assert wave[len(wave)//2]==500.
    assert len(audit['changes'])==256


def test_invalid_or_partial_grid_controls_refuse():
    with pytest.raises(ValueError):NativeFractionGrid([-.1,0,.1])
    with pytest.raises(ValueError):NativeFractionGrid([-.1,-.05,.01,.05,.1])
    with pytest.raises(ValueError):NativeLineGridControl(dict(grids=[],scope='invalid'),{}).apply(SimpleNamespace(lines=[]))
    with pytest.raises(ValueError):NativeLineGridControl({}, {},mode='fit-to-native')


def test_source_strength_mismatch_is_refused_before_any_override():
    lines=[Line()for _ in range(256)];offsets=np.linspace(-2e-4,2e-4,17)
    rows=[dict(index=k,i=0,j=1,f=.1,offsets_fraction=offsets.tolist())for k in range(256)]
    rows[0]['f']=.11
    with pytest.raises(ValueError,match='identity mismatch'):NativeLineGridControl(dict(grids=rows,scope='invalid'),{}).apply(SimpleNamespace(lines=lines))
    assert isinstance(lines[0].quadrature,LinearCoreExpWings)
